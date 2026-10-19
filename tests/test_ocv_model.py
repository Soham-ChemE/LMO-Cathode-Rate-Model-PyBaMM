"""Pseudo-OCV and electrode geometry sanity checks (no PyBaMM solve, fast)."""

import numpy as np
import pytest

from lmo.model import disc_geometry, ocv_to_sto
from lmo.ocv import build_pseudo_ocv, dqdv, dqdv_peaks, ocv_table_for_pybamm, sto_window, theoretical_capacity_mAh
from lmo.io import get_step


def test_pseudo_ocv_is_monotonic_and_within_window():
    ocv = build_pseudo_ocv()
    assert np.all(np.diff(ocv["ocv_V"]) > 0)
    assert 3.55 < ocv["ocv_V"].iloc[0] < 3.75
    assert 4.35 < ocv["ocv_V"].iloc[-1] < 4.5
    assert np.all(np.diff(ocv["stoichiometry"]) < 0)


def test_ocp_table_extends_past_both_cutoffs():
    table = ocv_table_for_pybamm(build_pseudo_ocv())
    assert np.all(np.diff(table[:, 0]) > 0)
    assert table[0, 1] > 4.8 and table[-1, 1] < 3.0


def test_stoichiometry_window_matches_theoretical_capacity():
    q = get_step("C/10", 2, "discharge").total_capacity_mAh
    lo, hi = sto_window(q)
    assert hi == 1.0
    assert (hi - lo) * theoretical_capacity_mAh() == pytest.approx(q)
    assert 0.2 < lo < 0.35


def test_dqdv_shows_two_spinel_plateaus():
    peaks = dqdv_peaks(dqdv(get_step("C/10", 2, "discharge")))
    assert len(peaks) == 2
    assert 3.95 < peaks[0]["voltage_V"] < 4.05
    assert 4.08 < peaks[1]["voltage_V"] < 4.2


def test_ocv_inversion_round_trip_and_clipping():
    assert ocv_to_sto(3.0) == pytest.approx(0.995)
    assert ocv_to_sto(5.0) < 0.26  # past the measured window, on the linear extension
    assert 0.5 < ocv_to_sto(4.05) < 0.7


@pytest.mark.parametrize("disc", [12, 14, 16])
def test_disc_geometry_is_self_consistent(disc):
    g = disc_geometry(disc)
    assert g["active_volume_fraction"] + g["inactive_volume_fraction"] + g["porosity"] == pytest.approx(1.0)
    assert 0 < g["porosity"] < 1
    # Active mass recovered from the volume fraction equals the mass implied by the capacity.
    from lmo.ocv import LMO_DENSITY_KG_M3, active_mass_g
    from lmo.model import COATING_THICKNESS_M

    mass = g["active_volume_fraction"] * g["area_m2"] * COATING_THICKNESS_M * LMO_DENSITY_KG_M3 * 1e3
    assert mass == pytest.approx(active_mass_g())
