"""Lumped cell thermal model with I^2R heating, optional entropic term, and a shutdown flag.

Sign convention matches ecm.py: current is positive on discharge. With the
charge-positive Bernardi form q = I_c (V - U) + I_c T dU/dT and I_c = -I_dis,
the reversible term becomes q_ent = -I_dis T dU/dT, so a cell with dU/dT < 0
generates reversible heat on discharge.
"""
from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass
class ThermalParams:
    mass_kg: float = 0.048
    cp_j_per_kgk: float = 900.0
    h_w_per_k: float = 0.15          # convective conductance to ambient
    t_ambient_c: float = 25.0
    t_shutdown_c: float = 80.0       # protective shutdown limit, well below runaway onset
    dudt_v_per_k: float = 0.0        # entropic coefficient, optional


class LumpedThermal:
    def __init__(self, params: ThermalParams | None = None, t0_c: float | None = None):
        self.p = params if params is not None else ThermalParams()
        self.temp_c = self.p.t_ambient_c if t0_c is None else t0_c
        self.shutdown = False

    def step(self, current_a: float, r_total_ohm: float, dt: float) -> float:
        """Advance one step with the exact exponential solution of the linear ODE.

        m cp dT/dt = q_gen - h (T - T_amb) has the closed form
        T_new = T_amb + q/h + (T - T_amb - q/h) exp(-h dt / (m cp)),
        which is unconditionally stable for any dt > 0.
        """
        if dt <= 0.0:
            raise ValueError("dt must be positive")
        p = self.p
        t_k = self.temp_c + 273.15
        q_joule = current_a ** 2 * r_total_ohm
        q_ent = -current_a * t_k * p.dudt_v_per_k
        q_gen = q_joule + q_ent
        x = p.h_w_per_k * dt / (p.mass_kg * p.cp_j_per_kgk)
        if x < 1e-6:
            # Near adiabatic limit of the exact solution, avoids cancellation.
            self.temp_c += (q_gen - p.h_w_per_k * (self.temp_c - p.t_ambient_c)) * dt / (p.mass_kg * p.cp_j_per_kgk)
        else:
            t_inf = p.t_ambient_c + q_gen / p.h_w_per_k
            self.temp_c = t_inf + (self.temp_c - t_inf) * math.exp(-x)
        if self.temp_c >= p.t_shutdown_c:
            self.shutdown = True
        return self.temp_c
