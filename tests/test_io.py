"""The parser must reproduce the summary values quoted in data/README.md within 1%."""

import numpy as np
import pytest

from lmo import CVS, RATES
from lmo.io import cv_peaks, gcd_summary, get_step, load_gcd

REPORT = {
    # rate: (cycle, discharge capacity mAh, mean voltage V, energy mWh)
    "C/10": (1, 1.430, 4.04, 5.77),
    "C/5": (1, 1.422, 4.02, 5.71),
    "C/2": (1, 1.361, 3.96, 5.39),
    "2C": (2, 0.823, 3.70, 3.04),
}

REPORT_CV = {
    "0.1 mV/s": (4.36, 1.62, 3.81, -1.28),
    "0.2 mV/s": (4.42, 2.42, 3.72, -2.06),
    "0.3 mV/s": (4.49, 3.12, 3.64, -2.73),
}


@pytest.mark.parametrize("rate", list(RATES))
def test_every_gcd_file_has_two_full_cycles(rate):
    df = load_gcd(rate)
    assert sorted(df["cycle"].unique()) == [1, 2]
    for cycle in (1, 2):
        for kind in ("charge", "discharge"):
            step = get_step(rate, cycle, kind, df)
            assert step.duration_s > 60
            assert step.total_capacity_mAh > 0.5


@pytest.mark.parametrize("rate", list(RATES))
def test_discharge_summary_matches_report(rate):
    cycle, cap, vavg, energy = REPORT[rate]
    step = get_step(rate, cycle, "discharge")
    assert step.total_capacity_mAh == pytest.approx(cap, rel=0.01)
    assert step.mean_voltage_V == pytest.approx(vavg, rel=0.01)
    assert step.energy_mWh == pytest.approx(energy, rel=0.01)


def test_measured_current_close_to_setpoint():
    for rate, meta in RATES.items():
        step = get_step(rate, 1, "discharge")
        # C/10 was on the 10 mA range and reads 2% low on discharge (data/README.md).
        assert step.mean_current_mA == pytest.approx(meta["current_mA"], rel=0.03)


def test_c10_current_range_quirk_is_visible():
    chg = get_step("C/10", 1, "charge")
    dis = get_step("C/10", 1, "discharge")
    assert chg.mean_current_mA > 0.145
    assert dis.mean_current_mA < 0.142


def test_2c_first_charge_is_larger_than_second():
    assert get_step("2C", 1, "charge").total_capacity_mAh > get_step("2C", 2, "charge").total_capacity_mAh
    assert get_step("2C", 1, "charge").total_capacity_mAh == pytest.approx(1.01, abs=0.01)


def test_summary_table_shape_and_ce():
    tab = gcd_summary()
    assert len(tab) == 16
    ce = tab.set_index(["rate", "cycle", "step"])["coulombic_efficiency"]
    assert ce[("C/10", 1, "discharge")] == pytest.approx(0.91, abs=0.01)
    assert (tab["coulombic_efficiency"] < 1.0).all()


@pytest.mark.parametrize("scan", list(CVS))
def test_cv_peaks_match_report(scan):
    vox, iox, vred, ired = REPORT_CV[scan]
    p = cv_peaks(scan, cycle=2)
    assert p["oxidation_peak_V"] == pytest.approx(vox, abs=0.01)
    assert p["oxidation_peak_mA"] == pytest.approx(iox, rel=0.01)
    assert p["reduction_peak_V"] == pytest.approx(vred, abs=0.01)
    assert p["reduction_peak_mA"] == pytest.approx(ired, rel=0.01)


def test_capacity_is_cumulative_and_monotonic():
    step = get_step("C/5", 1, "discharge")
    assert step.capacity_mAh[0] == 0.0
    assert np.all(np.diff(step.capacity_mAh) >= 0)
