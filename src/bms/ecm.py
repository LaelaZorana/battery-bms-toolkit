"""Equivalent circuit models (Thevenin 1RC and 2RC) for a Li-ion cell.

Sign convention: current is positive on discharge.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import numpy as np

# Representative NMC OCV table (SOC 0..1, volts). Used when no per-cell table is given.
OCV_SOC = np.array([0.0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0])
OCV_V = np.array([3.00, 3.35, 3.48, 3.58, 3.64, 3.68, 3.72, 3.79, 3.87, 3.95, 4.05, 4.12, 4.20])


def ocv(soc, temp_c: float = 25.0, dvdt: float = 0.0, soc_pts=None, v_pts=None) -> np.ndarray:
    """OCV(SOC) with optional linear temperature scaling hook (dV/dT in V/K).

    Pass soc_pts and v_pts to interpolate a per-cell table instead of the
    default NMC one.
    """
    sp = OCV_SOC if soc_pts is None else np.asarray(soc_pts, dtype=float)
    vp = OCV_V if v_pts is None else np.asarray(v_pts, dtype=float)
    s = np.clip(soc, 0.0, 1.0)
    return np.interp(s, sp, vp) + dvdt * (temp_c - 25.0)


def docv_dsoc(soc, eps: float = 1e-3, soc_pts=None, v_pts=None):
    """Numerical dOCV/dSOC by central difference. Accepts scalar or array SOC."""
    s = np.clip(np.asarray(soc, dtype=float), 0.0, 1.0)
    hi = np.minimum(s + eps, 1.0)
    lo = np.maximum(s - eps, 0.0)
    slope = (ocv(hi, soc_pts=soc_pts, v_pts=v_pts) - ocv(lo, soc_pts=soc_pts, v_pts=v_pts)) / (hi - lo)
    return float(slope) if np.isscalar(soc) or np.ndim(soc) == 0 else slope


def temp_scale(value: float, temp_c: float, ea_ratio: float = 3000.0) -> float:
    """Arrhenius style scaling hook for resistances relative to 25 C.

    ea_ratio is Ea/R in kelvin. The default 3000 K, about 25 kJ/mol, is
    illustrative and sits in the range commonly reported for Li-ion charge
    transfer resistance. Fit it to data before trusting it quantitatively.
    """
    t = temp_c + 273.15
    return value * np.exp(ea_ratio * (1.0 / t - 1.0 / 298.15))


@dataclass
class ECMParams:
    capacity_ah: float = 2.5
    r0: float = 0.020
    r1: float = 0.015
    c1: float = 2000.0
    r2: float = 0.010
    c2: float = 20000.0
    dvdt: float = 0.0
    n_rc: int = 1
    ocv_soc: np.ndarray | None = None   # per-cell OCV table, defaults to the NMC table
    ocv_v: np.ndarray | None = None

    def __post_init__(self):
        if self.n_rc not in (1, 2):
            raise ValueError("n_rc must be 1 or 2")

    def ocv(self, soc, temp_c: float = 25.0):
        return ocv(soc, temp_c, self.dvdt, self.ocv_soc, self.ocv_v)

    def docv_dsoc(self, soc):
        return docv_dsoc(soc, soc_pts=self.ocv_soc, v_pts=self.ocv_v)


@dataclass
class ECM:
    params: ECMParams = field(default_factory=ECMParams)
    soc: float = 1.0
    v_rc: np.ndarray = field(default_factory=lambda: np.zeros(2))
    temp_c: float = 25.0

    def reset(self, soc: float = 1.0):
        self.soc = soc
        self.v_rc = np.zeros(2)

    def step(self, current_a: float, dt: float) -> float:
        if dt <= 0.0:
            raise ValueError("dt must be positive")
        p = self.params
        q = p.capacity_ah * 3600.0
        self.soc = float(np.clip(self.soc - current_a * dt / q, 0.0, 1.0))
        r0 = temp_scale(p.r0, self.temp_c)
        rc = [(temp_scale(p.r1, self.temp_c), p.c1), (temp_scale(p.r2, self.temp_c), p.c2)][: p.n_rc]
        for i, (r, c) in enumerate(rc):
            a = np.exp(-dt / (r * c))
            self.v_rc[i] = a * self.v_rc[i] + r * (1.0 - a) * current_a
        return float(p.ocv(self.soc, self.temp_c) - current_a * r0 - self.v_rc[: p.n_rc].sum())

    def simulate(self, current: np.ndarray, dt: float, soc0: float = 1.0):
        self.reset(soc0)
        v = np.empty(len(current))
        s = np.empty(len(current))
        for k, i in enumerate(current):
            v[k] = self.step(float(i), dt)
            s[k] = self.soc
        return v, s


def coulomb_count(current: np.ndarray, dt: float, capacity_ah: float, soc0: float) -> np.ndarray:
    return np.clip(soc0 - np.cumsum(current) * dt / (capacity_ah * 3600.0), 0.0, 1.0)


def synthetic_drive_cycle(duration_s: float = 1800.0, dt: float = 1.0, i_max: float = 5.0, seed: int = 0) -> np.ndarray:
    """Piecewise accel/cruise/regen profile with mild noise. Positive is discharge."""
    if dt <= 0.0:
        raise ValueError("dt must be positive")
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, duration_s, dt)
    i = np.zeros_like(t)
    k = 0
    while k < len(t):
        seg = max(1, int(rng.integers(20, 90) / dt))
        kind = rng.choice(["accel", "cruise", "regen", "idle"], p=[0.35, 0.35, 0.2, 0.1])
        level = {"accel": i_max, "cruise": 0.4 * i_max, "regen": -0.5 * i_max, "idle": 0.0}[kind]
        i[k : k + seg] = level
        k += seg
    i += 0.05 * i_max * rng.standard_normal(len(t))
    return i
