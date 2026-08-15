"""Passive bleed-resistor balancing of a series pack."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from .ecm import ocv


@dataclass
class BalanceResult:
    t: np.ndarray
    soc: np.ndarray          # shape (steps, n_cells)
    time_to_balance_s: float
    energy_wasted_wh: float
    converged: bool


def simulate_passive_balancing(soc0, capacity_ah=2.5, r_bleed=33.0, tol=0.005,
                               dt=10.0, t_max=6 * 3600.0) -> BalanceResult:
    """Bleed every cell above (min_soc + tol) through r_bleed until the spread is within tol."""
    soc = np.array(soc0, dtype=float)
    q = capacity_ah * 3600.0
    hist, ts = [soc.copy()], [0.0]
    energy_j, t, done = 0.0, 0.0, False
    while t < t_max:
        spread = soc.max() - soc.min()
        if spread <= tol:
            done = True
            break
        v = ocv(soc)
        active = soc > soc.min() + tol
        i_bleed = np.where(active, v / r_bleed, 0.0)
        soc = soc - i_bleed * dt / q
        energy_j += float(np.sum(v * i_bleed) * dt)
        t += dt
        hist.append(soc.copy())
        ts.append(t)
    return BalanceResult(np.array(ts), np.array(hist), t, energy_j / 3600.0, done)
