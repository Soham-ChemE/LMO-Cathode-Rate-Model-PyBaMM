"""v2.1 variants: option mapping, the lognormal size distribution, and the residual's parameter list."""

import numpy as np
import pybamm
import pytest

from lmo.fit_v21 import EXTRA_BOUNDS_LOG10, R_MAX_REL, R_MIN_REL, VARIANTS, lognormal_area_weighted, param_hook, sim_kwargs
from lmo.model import PARTICLE_RADIUS_M, build_parameters


def test_variants_map_to_pybamm_options():
    assert VARIANTS["dl"]["options"] == {"surface form": "differential"}
    assert VARIANTS["psd"]["options"] == {"particle size": "distribution"}
    assert set(VARIANTS["psd+dl"]["extra"]) == {"C_dl", "sd_p"}
    for v in VARIANTS.values():
        assert all(n in EXTRA_BOUNDS_LOG10 for n in v["extra"])


@pytest.mark.parametrize("sd", [0.1, 0.3, 0.7, 1.5])
def test_lognormal_is_normalised_on_the_size_axis(sd):
    R_av = PARTICLE_RADIUS_M
    R = np.linspace(R_MIN_REL * R_av, R_MAX_REL * R_av, 20001)
    f = lognormal_area_weighted(pybamm.Array(R), R_av, sd).evaluate().flatten()
    assert np.trapezoid(f, R) == pytest.approx(1.0, abs=2e-3)
    assert np.all(f >= 0)


def test_lognormal_mean_is_the_v1_radius_when_narrow():
    R_av = PARTICLE_RADIUS_M
    R = np.linspace(R_MIN_REL * R_av, R_MAX_REL * R_av, 20001)
    f = lognormal_area_weighted(pybamm.Array(R), R_av, 0.1).evaluate().flatten()
    assert np.trapezoid(f * R, R) == pytest.approx(R_av, rel=0.02)


def test_param_hook_adds_the_variant_parameters():
    p = build_parameters()
    param_hook("psd+dl")(p)
    assert "Positive area-weighted particle-size distribution [m-1]" in p
    assert p["Positive minimum particle radius [m]"] == pytest.approx(R_MIN_REL * PARTICLE_RADIUS_M)
    assert isinstance(p["Positive electrode double-layer capacity [F.m-2]"], pybamm.InputParameter)
    kw = sim_kwargs("psd+dl")
    assert kw["var_pts_extra"]["R_p"] > 5 and kw["extra_inputs"] == ("C_dl", "sd_p")
    assert "var_pts_extra" not in sim_kwargs("dl")
