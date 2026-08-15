import math

from bms.thermal import LumpedThermal, ThermalParams


def test_entropic_heating_on_discharge_hand_value():
    # dU/dT < 0 and discharge current must generate reversible heat:
    # q_ent = -I T dU/dT = -5 * 298.15 * (-1e-4) = +0.1491 W at 25 C.
    # Near adiabatic isolates the generation terms.
    p = ThermalParams(dudt_v_per_k=-1e-4, h_w_per_k=1e-9)
    th = LumpedThermal(p)
    dt = 1.0
    t0 = th.temp_c
    th.step(5.0, 0.0, dt)  # r = 0 removes joule heat, only the entropic term remains
    q_ent = -5.0 * (t0 + 273.15) * p.dudt_v_per_k
    expected = t0 + q_ent * dt / (p.mass_kg * p.cp_j_per_kgk)
    assert q_ent > 0
    assert abs(th.temp_c - expected) < 1e-6


def test_adiabatic_joule_heating_hand_value():
    p = ThermalParams(h_w_per_k=1e-12)
    th = LumpedThermal(p)
    for _ in range(60):
        th.step(10.0, 0.02, 1.0)  # 2 W for 60 s
    expected = 25.0 + 2.0 * 60.0 / (p.mass_kg * p.cp_j_per_kgk)
    assert abs(th.temp_c - expected) < 1e-6


def test_large_dt_is_stable():
    th = LumpedThermal()
    for _ in range(20):
        th.step(5.0, 0.035, 1000.0)  # far beyond the explicit Euler limit
    t_inf = 25.0 + (25.0 * 0.035) / 0.15
    assert 25.0 < th.temp_c <= t_inf + 1e-9


def test_default_params_not_shared():
    a, b = LumpedThermal(), LumpedThermal()
    assert a.p is not b.p


def test_rejects_bad_dt():
    import pytest
    with pytest.raises(ValueError):
        LumpedThermal().step(1.0, 0.02, 0.0)
