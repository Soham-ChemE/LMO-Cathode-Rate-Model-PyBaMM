"""Structure of the v2 residual: which data enter it, how the rests are weighted, and that 2C never does."""

import numpy as np
import pytest

from lmo.fit import TRAINING
from lmo.fit_v2 import HEADLINE, JUMP_VARIANT, REST_KINDS, REST_TRAINING, WEIGHTINGS, V2Residual
from lmo.io import get_rest
from lmo.model import gcd_experiment


def test_2c_is_never_in_the_training_set():
    assert all(rate != "2C" for rate, _ in TRAINING)
    assert all(rate != "2C" for rate, _ in REST_TRAINING)
    assert set(r for r, _ in REST_TRAINING) == {"C/10", "C/5", "C/2"}


def test_rest_step_gets_its_own_output_period():
    exp = gcd_experiment(0.287, 1, rest_period_s=0.05)
    rest_steps = [s for s in exp.steps if "Rest" in str(s)]
    assert len(rest_steps) == 2
    assert all(s.period == pytest.approx(0.05) for s in rest_steps)
    # Without the argument the rests share the default period (v1 behaviour).
    exp = gcd_experiment(0.287, 1)
    assert all(s.period is None or s.period > 1.0 for s in exp.steps if "Rest" in str(s))


@pytest.fixture(scope="module")
def residual():
    return V2Residual("SPMe", lam=1.0, mode="trace")


def test_residual_blocks_and_sizes(residual):
    sizes = residual.block_sizes()
    assert len(residual.blocks) == 10  # 4 curve steps + 6 rests
    assert residual.n_curve == 4 * 300
    n_rest = sum(len(get_rest(r, c, a).time_s) for r, c in REST_TRAINING for a in REST_KINDS)
    assert residual.n_rest == n_rest
    assert sizes["C/5 charge"] == 300
    assert sizes["C/10 rest after charge"] == len(get_rest("C/10", 1, "charge").time_s)


def test_rest_weight_rule(residual):
    assert residual.w_rest == pytest.approx(np.sqrt(1.0 * residual.n_curve / residual.n_rest))
    assert residual.w_rest > 3  # a rest point counts for more than three curve points at lambda = 1
    unweighted = V2Residual.__new__(V2Residual)
    unweighted.__dict__.update(residual.__dict__)
    unweighted.lam = None
    assert WEIGHTINGS["unweighted"] == (None, "trace")


def test_jump_only_variant_uses_one_point_per_rest():
    lam, mode = WEIGHTINGS[JUMP_VARIANT]
    assert mode == "jump"
    res = V2Residual("SPMe", lam=lam, mode=mode)
    assert res.n_rest == 6
    assert all(len(idx) == 1 and idx[0] == 0 for idx in res.rest_idx.values())
    assert res.w_rest == pytest.approx(np.sqrt(lam * 1200 / 6))


def test_headline_is_full_traces_at_lambda_one():
    assert WEIGHTINGS[HEADLINE] == (1.0, "trace")
