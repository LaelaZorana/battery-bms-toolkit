"""Identify R0, R1, C1 from an HPPC style pulse by least squares."""
from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares
from .ecm import ECM, ECMParams


def hppc_pulse(i_pulse: float = 5.0, t_pulse: float = 10.0, t_rest: float = 60.0, dt: float = 0.1):
    n_p, n_r = int(t_pulse / dt), int(t_rest / dt)
    return np.concatenate([np.zeros(n_r // 6), np.full(n_p, i_pulse), np.zeros(n_r)]), dt


def _model_voltage(theta, current, dt, soc0, capacity_ah):
    r0, r1, c1 = theta
    m = ECM(ECMParams(capacity_ah=capacity_ah, r0=r0, r1=r1, c1=c1, n_rc=1))
    v, _ = m.simulate(current, dt, soc0)
    return v


def identify_1rc(current, voltage, dt, soc0, capacity_ah, x0=(0.01, 0.01, 1000.0)):
    """Return (r0, r1, c1) fitted to the measured pulse response."""
    def resid(theta):
        return _model_voltage(theta, current, dt, soc0, capacity_ah) - voltage
    scale = np.array([0.01, 0.01, 1000.0])
    res = least_squares(lambda z: resid(z * scale), np.array(x0) / scale,
                        bounds=(np.array([1e-4, 1e-4, 1.0]) / scale, np.array([1.0, 1.0, 1e6]) / scale))
    return tuple(res.x * scale)


def initial_guess_from_pulse(current, voltage, dt):
    """Closed form seed: R0 from the instantaneous drop, R1 from the slow drop, C1 from a 63% time."""
    k_on = int(np.argmax(np.abs(current) > 0))
    k_off = k_on + int(np.argmax(np.abs(current[k_on:]) == 0))
    i = float(current[k_on])
    r0 = (voltage[k_on - 1] - voltage[k_on]) / i
    r1 = (voltage[k_on] - voltage[k_off - 1]) / i
    v_rest = voltage[k_off - 1:]
    target = v_rest[0] + 0.632 * (v_rest[-1] - v_rest[0])
    tau = dt * int(np.argmax(v_rest >= target))
    return max(r0, 1e-4), max(r1, 1e-4), max(tau / max(r1, 1e-4), 1.0)
