"""The rest parser must find every 30 s post-step rest and reproduce the raw-file jumps."""

import numpy as np
import pytest

from lmo import RATES
from lmo.io import all_rests, get_rest, get_step, load_gcd, rest_summary


@pytest.mark.parametrize("rate", list(RATES))
def test_every_cycle_has_a_rest_after_each_step(rate):
    df = load_gcd(rate)
    for cycle in (1, 2):
        for after in ("charge", "discharge"):
            r = get_rest(rate, cycle, after, df)
            assert r.rate == rate and r.cycle == cycle and r.after == after
            assert 29.0 < r.duration_s < 31.0
            assert len(r.time_s) >= 10
            assert np.all(np.diff(r.time_s) > 0)


def test_sixteen_rests_in_total():
    rests = all_rests()
    assert len(rests) == 16
    assert len(rest_summary()) == 16


@pytest.mark.parametrize("rate", list(RATES))
def test_rest_starts_right_after_the_step_it_follows(rate):
    df = load_gcd(rate)
    for after in ("charge", "discharge"):
        step = get_step(rate, 1, after, df)
        rest = get_rest(rate, 1, after, df)
        assert rest.v_before_V == pytest.approx(step.voltage_V[-1])
        assert rest.current_before_mA == pytest.approx(step.current_mA[-1])
        # The step ended at its cut-off voltage.
        assert rest.v_before_V == pytest.approx(4.5 if after == "charge" else 3.5, abs=2e-3)


def test_first_sample_delay_is_tens_of_milliseconds():
    for r in all_rests():
        assert 0.01 < r.first_sample_delay_s < 0.5
        assert r.time_s[0] == pytest.approx(r.first_sample_delay_s)


def test_jump_sign_and_direction_of_relaxation():
    for r in all_rests():
        sign = -1.0 if r.after == "charge" else 1.0
        assert sign * r.jump_V > 0  # voltage drops after charge, rises after discharge
        assert sign * r.total_relaxation_V > sign * r.jump_V  # and keeps relaxing the same way
        assert np.all(sign * np.diff(r.voltage_V) >= -1e-3)


def test_instantaneous_resistance_after_charge_matches_raw_files():
    # Values read directly off the raw exports: first rest sample vs last sample under current.
    expected = {"C/10": 70.7, "C/5": 67.3, "C/2": 70.0, "2C": 74.7}
    for rate, ohm in expected.items():
        r = get_rest(rate, 1, "charge")
        assert r.instantaneous_resistance_ohm == pytest.approx(ohm, abs=0.5)
        assert r.instantaneous_resistance_ohm < 96.5  # below the v1 fitted contact resistance


def test_c5_cycle2_discharge_rest_shows_the_contact_disturbance():
    # data/README.md: the second C/5 discharge settled ~22 mV lower after a contact disturbance;
    # that extra drop comes back at the switch-off, so its rest jump is larger than cycle 1's.
    r1 = get_rest("C/5", 1, "discharge")
    r2 = get_rest("C/5", 2, "discharge")
    assert r2.jump_V - r1.jump_V == pytest.approx(0.018, abs=0.005)


def test_missing_rest_raises():
    with pytest.raises(KeyError):
        get_rest("C/5", 3, "charge")
    with pytest.raises(ValueError):
        get_rest("C/5", 1, "rest")
