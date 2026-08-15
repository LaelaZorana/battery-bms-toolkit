"""Lumped cell thermal model with I^2R heating, optional entropic term, and runaway flag."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class ThermalParams:
    mass_kg: float = 0.048
    cp_j_per_kgk: float = 900.0
    h_w_per_k: float = 0.15          # convective conductance to ambient
    t_ambient_c: float = 25.0
    t_runaway_c: float = 80.0
    dudt_v_per_k: float = 0.0        # entropic coefficient, optional


class LumpedThermal:
    def __init__(self, params: ThermalParams = ThermalParams(), t0_c: float | None = None):
        self.p = params
        self.temp_c = params.t_ambient_c if t0_c is None else t0_c
        self.runaway = False

    def step(self, current_a: float, r_total_ohm: float, dt: float, temp_k_for_entropic: float | None = None) -> float:
        p = self.p
        q_joule = current_a ** 2 * r_total_ohm
        t_k = (self.temp_c + 273.15) if temp_k_for_entropic is None else temp_k_for_entropic
        q_ent = current_a * t_k * p.dudt_v_per_k
        q_out = p.h_w_per_k * (self.temp_c - p.t_ambient_c)
        self.temp_c += (q_joule + q_ent - q_out) * dt / (p.mass_kg * p.cp_j_per_kgk)
        if self.temp_c >= p.t_runaway_c:
            self.runaway = True
        return self.temp_c
