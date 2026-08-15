import numpy as np
from bms.ecm import ECM, ECMParams, synthetic_drive_cycle
from bms.estimation import SOCEKF, SOCUKF, CapacityRLS

def _truth(seed=1):
    p = ECMParams(); dt = 1.0
    i = synthetic_drive_cycle(1800, dt, i_max=5.0, seed=seed)
    v, s = ECM(p).simulate(i, dt, soc0=0.9)
    rng = np.random.default_rng(0)
    vn = v + 0.005 * rng.standard_normal(len(v))
    return p, dt, i, vn, s

def _run(filter_cls, filter_params=None, i_bias=0.0):
    p, dt, i, vn, s = _truth()
    f = filter_cls(filter_params or p, dt, soc0=0.7)
    est = np.array([f.update(float(ik) + i_bias, float(vk)) for ik, vk in zip(i, vn)])
    return est, s

def test_ekf_converges_from_20pct_error():
    est, s = _run(SOCEKF)
    assert abs(est[-1] - s[-1]) < 0.02
    assert np.mean(np.abs(est[-300:] - s[-300:])) < 0.02

def test_ukf_converges_from_20pct_error():
    est, s = _run(SOCUKF)
    assert abs(est[-1] - s[-1]) < 0.02

def test_ekf_readme_scenario_with_current_bias():
    # The README validation row: 20 percent initial error, 5 mV noise,
    # 0.15 A current bias fed to the filter. Enforces the published bound.
    est, s = _run(SOCEKF, i_bias=0.15)
    assert abs(est[-1] - s[-1]) < 0.02

def test_coulomb_counting_drifts_under_bias():
    from bms.ecm import coulomb_count
    p, dt, i, vn, s = _truth()
    cc = coulomb_count(i + 0.15, dt, p.capacity_ah, 0.7)
    assert abs(cc[-1] - s[-1]) > 0.20  # the README's 23 percent end error

def test_ekf_survives_model_mismatch():
    # Filter parameters perturbed 15 percent plus a 0.15 A current bias, so the
    # validation is no longer against the filter's own generating model.
    p = ECMParams()
    wrong = ECMParams(capacity_ah=p.capacity_ah * 1.15, r0=p.r0 * 0.85,
                      r1=p.r1 * 1.15, c1=p.c1 * 0.85)
    est, s = _run(SOCEKF, filter_params=wrong, i_bias=0.15)
    assert abs(est[-1] - s[-1]) < 0.05

def test_filters_skip_nan_samples():
    p = ECMParams()
    for cls in (SOCEKF, SOCUKF):
        f = cls(p, 1.0, soc0=0.8)
        f.update(1.0, 3.9)
        x_before = f.x.copy()
        f.update(float("nan"), 3.9)
        f.update(1.0, float("nan"))
        assert np.all(np.isfinite(f.x))
        assert np.allclose(f.x, x_before)

def test_rls_tracks_capacity_with_physical_signs():
    # Discharging: soc falls (delta_soc negative), ah_moved positive.
    r = CapacityRLS(2.5)
    q_true = 2.2
    rng = np.random.default_rng(0)
    for _ in range(50):
        d = rng.uniform(0.2, 0.6)
        r.update(-d, q_true * d + 0.002 * rng.standard_normal())
    assert abs(r.q - q_true) / q_true < 0.02
    assert abs(r.capacity_ah - r.q) == 0.0
    assert abs(r.soh(2.5) - q_true / 2.5) < 0.02
