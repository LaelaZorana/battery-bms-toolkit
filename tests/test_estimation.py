import numpy as np
from bms.ecm import ECM, ECMParams, synthetic_drive_cycle
from bms.estimation import SOCEKF, SOCUKF, CapacityRLS

def _run(filter_cls):
    p = ECMParams(); dt = 1.0
    i = synthetic_drive_cycle(1800, dt, i_max=5.0, seed=1)
    v, s = ECM(p).simulate(i, dt, soc0=0.9)
    rng = np.random.default_rng(0)
    vn = v + 0.005 * rng.standard_normal(len(v))
    f = filter_cls(p, dt, soc0=0.7)
    est = np.array([f.update(float(ik), float(vk)) for ik, vk in zip(i, vn)])
    return est, s

def test_ekf_converges_from_20pct_error():
    est, s = _run(SOCEKF)
    assert abs(est[-1] - s[-1]) < 0.02
    assert np.mean(np.abs(est[-300:] - s[-300:])) < 0.02

def test_ukf_converges_from_20pct_error():
    est, s = _run(SOCUKF)
    assert abs(est[-1] - s[-1]) < 0.02

def test_rls_tracks_capacity():
    r = CapacityRLS(2.5)
    q_true = 2.2
    rng = np.random.default_rng(0)
    for _ in range(50):
        d = rng.uniform(0.2, 0.6)
        r.update(d, q_true * d + 0.002 * rng.standard_normal())
    assert abs(r.q - q_true) / q_true < 0.02
