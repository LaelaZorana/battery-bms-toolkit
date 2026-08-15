import numpy as np
from bms.ecm import ECM, ECMParams, ocv

def test_1rc_step_response_matches_analytic():
    p = ECMParams(capacity_ah=1e6, r0=0.02, r1=0.015, c1=2000.0, n_rc=1)  # huge Q freezes SOC
    m = ECM(p); m.reset(0.5)
    dt, i, n = 0.1, 2.0, 3000
    v = np.array([m.step(i, dt) for _ in range(n)])
    t = dt * np.arange(1, n + 1)
    v_an = ocv(0.5) - i * p.r0 - i * p.r1 * (1 - np.exp(-t / (p.r1 * p.c1)))
    assert np.max(np.abs(v - v_an)) < 1e-6

def test_ocv_monotonic():
    s = np.linspace(0, 1, 101)
    assert np.all(np.diff(ocv(s)) > 0)
