from bms.ecm import ECM, ECMParams
from bms.param_id import hppc_pulse, identify_1rc, initial_guess_from_pulse

def test_identify_recovers_params():
    true = ECMParams(capacity_ah=2.5, r0=0.025, r1=0.012, c1=1500.0, n_rc=1)
    i, dt = hppc_pulse(5.0, 10.0, 60.0, 0.1)
    v, _ = ECM(true).simulate(i, dt, 0.8)
    x0 = initial_guess_from_pulse(i, v, dt)
    r0, r1, c1 = identify_1rc(i, v, dt, 0.8, 2.5, x0=x0)
    assert abs(r0 - true.r0) / true.r0 < 0.05
    assert abs(r1 - true.r1) / true.r1 < 0.05
    assert abs(c1 - true.c1) / true.c1 < 0.05
