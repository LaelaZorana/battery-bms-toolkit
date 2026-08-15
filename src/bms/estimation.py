"""SOC estimation (EKF, UKF) on the 1RC model and SOH capacity tracking (RLS).

Both filters use the OCV table carried by their ECMParams, so a flat chemistry
like LFP degrades observability exactly as it should. The measurement models
assume 25 C. A NaN current or voltage sample is skipped rather than folded into
the state.
"""
from __future__ import annotations

import math

import numpy as np
from .ecm import ECMParams


class SOCEKF:
    """State x = [soc, v1]. Measurement is terminal voltage."""

    def __init__(self, params: ECMParams, dt: float, soc0: float = 0.5,
                 p0=(0.05, 1e-3), q=(1e-8, 1e-6), r=1e-3):
        self.p, self.dt = params, dt
        self.x = np.array([soc0, 0.0])
        self.P = np.diag(p0)
        self.Q = np.diag(q)
        self.R = r

    def _f(self, x, i):
        p = self.p
        a = np.exp(-self.dt / (p.r1 * p.c1))
        return np.array([x[0] - i * self.dt / (p.capacity_ah * 3600.0), a * x[1] + p.r1 * (1 - a) * i]), a

    def _h(self, x, i):
        return float(self.p.ocv(x[0]) - x[1] - i * self.p.r0)

    def update(self, i: float, v_meas: float) -> float:
        if not (math.isfinite(i) and math.isfinite(v_meas)):
            return float(self.x[0])
        xp, a = self._f(self.x, i)
        F = np.array([[1.0, 0.0], [0.0, a]])
        P = F @ self.P @ F.T + self.Q
        H = np.array([[self.p.docv_dsoc(xp[0]), -1.0]])
        S = (H @ P @ H.T).item() + self.R
        K = (P @ H.T / S).ravel()
        self.x = xp + K * (v_meas - self._h(xp, i))
        self.x[0] = float(np.clip(self.x[0], 0.0, 1.0))
        self.P = (np.eye(2) - np.outer(K, H)) @ P
        return float(self.x[0])


class SOCUKF:
    """Sigma-point (unscented) filter on the same 1RC state."""

    def __init__(self, params: ECMParams, dt: float, soc0: float = 0.5,
                 p0=(0.05, 1e-3), q=(1e-8, 1e-6), r=1e-3, alpha=1e-1, beta=2.0, kappa=0.0):
        self.p, self.dt = params, dt
        self.x = np.array([soc0, 0.0])
        self.P = np.diag(p0)
        self.Q = np.diag(q)
        self.R = r
        n = 2
        self.lam = alpha ** 2 * (n + kappa) - n
        self.wm = np.full(2 * n + 1, 1.0 / (2 * (n + self.lam)))
        self.wc = self.wm.copy()
        self.wm[0] = self.lam / (n + self.lam)
        self.wc[0] = self.wm[0] + (1 - alpha ** 2 + beta)

    def _sigmas(self):
        n = 2
        L = np.linalg.cholesky((n + self.lam) * self.P)
        return np.vstack([self.x, self.x + L.T, self.x - L.T])

    def update(self, i: float, v_meas: float) -> float:
        if not (math.isfinite(i) and math.isfinite(v_meas)):
            return float(self.x[0])
        p = self.p
        a = np.exp(-self.dt / (p.r1 * p.c1))
        X = self._sigmas()
        Xp = np.column_stack([X[:, 0] - i * self.dt / (p.capacity_ah * 3600.0), a * X[:, 1] + p.r1 * (1 - a) * i])
        xm = self.wm @ Xp
        Pm = (self.wc[:, None] * (Xp - xm)).T @ (Xp - xm) + self.Q
        Y = p.ocv(Xp[:, 0]) - Xp[:, 1] - i * p.r0
        ym = self.wm @ Y
        Pyy = float(self.wc @ (Y - ym) ** 2) + self.R
        Pxy = (self.wc[:, None] * (Xp - xm)).T @ (Y - ym)
        K = Pxy / Pyy
        self.x = xm + K * (v_meas - ym)
        self.x[0] = float(np.clip(self.x[0], 0.0, 1.0))
        self.P = Pm - np.outer(K, K) * Pyy
        self.P = 0.5 * (self.P + self.P.T) + 1e-12 * np.eye(2)
        return float(self.x[0])


class CapacityRLS:
    """Recursive least squares on capacity from ah_moved = -Q * delta_soc.

    delta_soc is the physically signed change soc_end - soc_start, and ah_moved
    is the discharged charge in Ah, positive on discharge, matching the repo's
    current sign convention. Discharging by 0.3 SOC on a 2.2 Ah cell means
    update(-0.3, 0.66). The regressor is phi = -delta_soc so that y = Q * phi.
    """

    def __init__(self, q0_ah: float, lam: float = 0.98, p0: float = 10.0):
        self.q = q0_ah
        self.P = p0
        self.lam = lam

    def update(self, delta_soc: float, ah_moved: float) -> float:
        phi = -delta_soc
        if abs(phi) < 1e-6:
            return self.q
        k = self.P * phi / (self.lam + phi * self.P * phi)
        self.q += k * (ah_moved - phi * self.q)
        self.P = (self.P - k * phi * self.P) / self.lam
        return self.q

    @property
    def capacity_ah(self) -> float:
        """Current capacity estimate in Ah."""
        return self.q

    def soh(self, q_nominal_ah: float) -> float:
        """State of health as a fraction of the given nominal capacity."""
        return self.q / q_nominal_ah
