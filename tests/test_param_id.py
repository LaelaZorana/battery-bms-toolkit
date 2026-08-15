import numpy as np
import pytest
from bms.ecm import ECM, ECMParams
from bms.param_id import hppc_pulse, identify_1rc, initial_guess_from_pulse

def _pulse_data():
    true = ECMParams(capacity_ah=2.5, r0=0.025, r1=0.012, c1=1500.0, n_rc=1)
    i, dt = hppc_pulse(5.0, 10.0, 60.0, 0.1)
    v, _ = ECM(true).simulate(i, dt, 0.8)
    return true, i, v, dt

def test_identify_recovers_params():
    true, i, v, dt = _pulse_data()
    x0 = initial_guess_from_pulse(i, v, dt)
    r0, r1, c1 = identify_1rc(i, v, dt, 0.8, 2.5, x0=x0)
    assert abs(r0 - true.r0) / true.r0 < 0.05
    assert abs(r1 - true.r1) / true.r1 < 0.05
    assert abs(c1 - true.c1) / true.c1 < 0.05

def test_seed_is_close_under_short_pulse():
    # The default pulse is 10 s against an 18 s tau, so the uncorrected slow
    # drop identity would be about 2x off. The corrected seed must stay close.
    true, i, v, dt = _pulse_data()
    r0, r1, c1 = initial_guess_from_pulse(i, v, dt)
    assert abs(r0 - true.r0) / true.r0 < 0.10
    assert abs(r1 - true.r1) / true.r1 < 0.30
    assert abs(r1 * c1 - true.r1 * true.c1) / (true.r1 * true.c1) < 0.30

def test_seed_rejects_bad_inputs():
    _, i, v, dt = _pulse_data()
    with pytest.raises(ValueError):
        initial_guess_from_pulse(np.zeros_like(i), v, dt)
    with pytest.raises(ValueError):
        initial_guess_from_pulse(np.full(100, 5.0), v[:100], dt)  # pulse never ends
    with pytest.raises(ValueError):
        initial_guess_from_pulse(np.concatenate([[5.0], i[i == 0]]),
                                 np.full(1 + int(np.sum(i == 0)), 3.9), dt)  # starts at sample 0

def test_identify_warns_on_wrong_soc0():
    _, i, v, dt = _pulse_data()
    with pytest.warns(UserWarning):
        identify_1rc(i, v, dt, 0.7, 2.5)  # soc0 wrong by 0.1
