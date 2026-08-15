"""Identify R0, R1, C1 from an HPPC style pulse by least squares."""
from __future__ import annotations

import warnings

import numpy as np
from scipy.optimize import least_squares
from .ecm import ECM, ECMParams


def hppc_pulse(i_pulse: float = 5.0, t_pulse: float = 10.0, t_rest: float = 60.0, dt: float = 0.1):
    """Zero pre-rest, current pulse, then rest. The short pre-rest gives the seed a
    clean pre-pulse voltage sample. Returns (current, dt) so profile and step size
    travel together."""
    n_p, n_r = int(t_pulse / dt), int(t_rest / dt)
    return np.concatenate([np.zeros(n_r // 6), np.full(n_p, i_pulse), np.zeros(n_r)]), dt


def _model_voltage(theta, current, dt, soc0, capacity_ah, ocv_soc=None, ocv_v=None):
    r0, r1, c1 = theta
    m = ECM(ECMParams(capacity_ah=capacity_ah, r0=r0, r1=r1, c1=c1, n_rc=1,
                      ocv_soc=ocv_soc, ocv_v=ocv_v))
    v, _ = m.simulate(current, dt, soc0)
    return v


def identify_1rc(current, voltage, dt, soc0, capacity_ah, x0=(0.01, 0.01, 1000.0),
                 ocv_soc=None, ocv_v=None):
    """Return (r0, r1, c1) fitted to the measured pulse response.

    Pass the cell's OCV table through ocv_soc and ocv_v when the data did not
    come from the default NMC table. Warns when the residual stays large or a
    parameter lands on a bound, both of which mean the fit is not trustworthy,
    usually because soc0 or the OCV table is wrong.
    """
    def resid(theta):
        return _model_voltage(theta, current, dt, soc0, capacity_ah, ocv_soc, ocv_v) - voltage
    scale = np.array([0.01, 0.01, 1000.0])
    lo = np.array([1e-4, 1e-4, 1.0])
    hi = np.array([1.0, 1.0, 1e6])
    res = least_squares(lambda z: resid(z * scale), np.array(x0) / scale,
                        bounds=(lo / scale, hi / scale))
    fitted = res.x * scale
    rms = float(np.sqrt(np.mean(res.fun ** 2)))
    if not res.success or rms > 0.005:
        warnings.warn(f"identify_1rc: poor fit, residual rms {rms:.4g} V. Check soc0 and the data.")
    on_bound = (np.isclose(fitted, lo, rtol=1e-3) | np.isclose(fitted, hi, rtol=1e-3))
    if np.any(on_bound):
        warnings.warn("identify_1rc: a fitted parameter sits on a bound, the estimate is not trustworthy.")
    return tuple(fitted)


def initial_guess_from_pulse(current, voltage, dt):
    """Closed form seed from a single pulse.

    R0 comes from the instantaneous drop at pulse start. The slow drop during
    the pulse equals i R1 (1 - exp(-t_p / tau)), not i R1, so the seed solves
    that pair self consistently: estimate tau from the 63 percent relaxation
    time after the pulse, then divide the slow drop by the settling factor.
    Without the correction the R1 seed is biased low whenever the pulse is not
    long against tau, by about 2x for the default 10 s pulse on an 18 s tau.
    """
    current = np.asarray(current, dtype=float)
    voltage = np.asarray(voltage, dtype=float)
    if not np.any(np.abs(current) > 0):
        raise ValueError("initial_guess_from_pulse: current is all zero, no pulse to fit")
    k_on = int(np.argmax(np.abs(current) > 0))
    if k_on == 0:
        raise ValueError("initial_guess_from_pulse: the pulse must start after at least one rest sample")
    if not np.any(np.abs(current[k_on:]) == 0):
        raise ValueError("initial_guess_from_pulse: the pulse never ends, a rest tail is required")
    k_off = k_on + int(np.argmax(np.abs(current[k_on:]) == 0))
    i = float(current[k_on])
    t_pulse = (k_off - k_on) * dt
    r0 = (voltage[k_on - 1] - voltage[k_on]) / i
    slow_drop = (voltage[k_on] - voltage[k_off - 1]) / i
    # Measure tau on the RC relaxation only, starting after the instantaneous
    # R0 recovery jump at k_off, otherwise the 63 percent target lands inside
    # the ohmic step and tau collapses to one sample.
    v_rest = voltage[k_off:]
    target = v_rest[0] + 0.632 * (v_rest[-1] - v_rest[0])
    tau = max(dt * int(np.argmax(v_rest >= target)), dt)
    settle = 1.0 - np.exp(-t_pulse / tau)
    r1 = slow_drop / max(settle, 1e-3)
    return max(r0, 1e-4), max(r1, 1e-4), max(tau / max(r1, 1e-4), 1.0)
