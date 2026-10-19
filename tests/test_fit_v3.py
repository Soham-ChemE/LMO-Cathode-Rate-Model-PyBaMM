"""v3: the oxidation side reaction has the right sign and potential dependence, the residual
carries the excess-charge constraint, and 2C still never enters the fit."""

import numpy as np
import pybamm
import pytest

from lmo import oxidation
from lmo.fit_v3 import BOUNDS, EXCESS_TRAINING, EXTRA, MODEL_OPTIONS, V3Residual, leakage_current_table, rest_rmse_split
from lmo.io import get_step
from lmo.model import build_parameters, make_model, make_simulation, ocv_to_sto


def test_v3_options_add_the_side_reaction_to_the_double_layer_model():
    assert MODEL_OPTIONS["surface form"] == "differential"
    assert MODEL_OPTIONS[oxidation.OPTION] == "Tafel"
    assert EXTRA == ("C_dl", "j_ref", "alpha_ox")
    assert all(n in BOUNDS for n in ("Ds", "k", "R_contact") + EXTRA)


def test_make_model_attaches_the_oxidation_submodel_only_when_asked():
    plain = make_model("SPMe", **{"surface form": "differential"})
    assert isinstance(plain.submodels["positive primary sei"], pybamm.sei.NoSEI)
    with_ox = make_model("SPMe", **MODEL_OPTIONS)
    assert isinstance(with_ox.submodels["positive primary sei"], oxidation.OxidationCurrent)
    assert isinstance(with_ox.submodels["negative primary sei"], pybamm.sei.NoSEI)  # nothing on the lithium metal
    assert with_ox.options["SEI film resistance"] == "none"
    with pytest.raises(ValueError):
        make_model("SPMe", **{oxidation.OPTION: "linear"})


def _rest_at(voltage, j_ref, alpha, seconds=120.0):
    params = build_parameters(initial_sto=ocv_to_sto(voltage), Ds=2.5e-14, k=1.1e-11, R_contact=89.0)
    params.update(
        {"Positive electrode double-layer capacity [F.m-2]": 3.0, **oxidation.oxidation_parameters(j_ref, alpha)},
        check_already_exists=False,
    )
    sim = make_simulation(pybamm.Experiment([f"Rest for {seconds} seconds"], period="1 second"), params, "SPMe", make_model("SPMe", **MODEL_OPTIONS))
    sol = sim.solve()
    return sol["Voltage [V]"].entries, sol[oxidation.PARASITIC_CHARGE_VARIABLE].entries, sol["X-averaged positive electrode SEI interfacial current density [A.m-2]"].entries


def test_oxidation_is_anodic_and_self_discharges_the_cell_at_open_circuit():
    v, q, j = _rest_at(4.45, j_ref=1.0, alpha=0.5)
    assert np.all(j > 0)  # anodic (oxidation) current
    assert v[-1] < v[0] - 0.01  # re-intercalation pulls the voltage down
    assert q[-1] > q[0] and q[0] == pytest.approx(0.0, abs=1e-12)  # the counter accumulates the parasitic charge
    v0, q0, _ = _rest_at(4.45, j_ref=1e-30, alpha=0.5)
    assert abs(v0[-1] - v0[0]) < 1e-3 and q0[-1] == pytest.approx(0.0, abs=1e-12)


def test_oxidation_current_grows_with_potential_at_the_tafel_slope():
    _, _, j_hi = _rest_at(4.45, j_ref=1e-3, alpha=0.5, seconds=1.0)
    _, _, j_lo = _rest_at(4.20, j_ref=1e-3, alpha=0.5, seconds=1.0)
    ratio = j_hi[0] / j_lo[0]
    # 0.25 V at 118 mV/decade is about two decades; the surface potential differs slightly from the OCV
    assert 50 < ratio < 300


def test_residual_blocks_and_the_excess_constraint():
    res = V3Residual("SPMe")
    assert len(res.blocks) == 12  # 4 curve steps + 6 rests + 2 excess-charge constraints
    assert [b for b in res.blocks if b[0] == "excess"] == [("excess", "C/10", "cycle 1"), ("excess", "C/5", "cycle 1")]
    assert all(rate != "2C" for rate, _ in EXCESS_TRAINING)
    assert all(b[1] != "2C" for b in res.blocks)
    for rate, cyc in EXCESS_TRAINING:
        chg, dis = get_step(rate, cyc, "charge"), get_step(rate, cyc, "discharge")
        assert res.excess_meas[rate] == pytest.approx(chg.total_capacity_mAh - dis.total_capacity_mAh)
        assert res.excess_meas[rate] > 0
    # Smaller measured excess -> larger weight per mAh (same relative tolerance).
    assert res.w_excess["C/5"] > res.w_excess["C/10"]
    assert res.param_names == ["Ds", "k", "R_contact", "C_dl", "j_ref", "alpha_ox"]
    assert res.block_sizes()["C/10 excess charge"] == 1


def test_rest_rmse_split_pools_the_blocks_by_what_they_follow():
    per = {"C/10 rest after charge": 10.0, "C/5 rest after charge": 10.0, "C/2 rest after charge": 10.0, "C/10 rest after discharge": 2.0, "C/5 rest after discharge": 2.0, "C/2 rest after discharge": 2.0}
    n = {k: 5 for k in per}
    out = rest_rmse_split(per, n)
    assert out["rests_after_charge_mV"] == pytest.approx(10.0)
    assert out["rests_after_discharge_mV"] == pytest.approx(2.0)


def test_leakage_table_tafel_slope_and_reference_current():
    tab = leakage_current_table({"j_ref": 1e-2, "alpha_ox": 0.5})
    assert tab["tafel_slope_mV_per_decade"] == pytest.approx(118.3, abs=0.5)
    assert tab["current_at_reference_potential_mA"] == pytest.approx(1000 * 1e-2 * tab["particle_surface_area_m2"])
    i = tab["current_mA_vs_surface_potential_V"]
    assert i["4.6"] / i["4.5"] == pytest.approx(10 ** (100 / 118.3), rel=0.02)
