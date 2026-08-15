"""Protection limits with hysteresis and a small BMS state machine.

States: idle, active, warning, fault. A fault latches until reset() is called and
all limits are back inside their release thresholds.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum


class State(str, Enum):
    IDLE = "idle"
    ACTIVE = "active"
    WARNING = "warning"
    FAULT = "fault"


@dataclass
class Limits:
    v_max: float = 4.20
    v_min: float = 3.00
    v_hyst: float = 0.05
    i_max_a: float = 10.0
    i_hyst_a: float = 1.0
    t_max_c: float = 60.0
    t_hyst_c: float = 5.0
    warn_frac: float = 0.95   # fraction of a limit at which warning is raised


@dataclass
class SafetyMonitor:
    limits: Limits = field(default_factory=Limits)
    state: State = State.IDLE
    flags: dict = field(default_factory=lambda: {"ov": False, "uv": False, "oc": False, "ot": False})

    def _hyst(self, name: str, tripped: bool, released: bool):
        if tripped:
            self.flags[name] = True
        elif released:
            self.flags[name] = False

    def check(self, v_cell_max: float, v_cell_min: float, current_a: float, temp_c: float) -> State:
        L = self.limits
        self._hyst("ov", v_cell_max > L.v_max, v_cell_max < L.v_max - L.v_hyst)
        self._hyst("uv", v_cell_min < L.v_min, v_cell_min > L.v_min + L.v_hyst)
        self._hyst("oc", abs(current_a) > L.i_max_a, abs(current_a) < L.i_max_a - L.i_hyst_a)
        self._hyst("ot", temp_c > L.t_max_c, temp_c < L.t_max_c - L.t_hyst_c)
        any_fault = any(self.flags.values())
        near = (v_cell_max > L.v_max * L.warn_frac + (1 - L.warn_frac) * L.v_min * 0
                or v_cell_min < L.v_min * (2 - L.warn_frac)
                or abs(current_a) > L.i_max_a * L.warn_frac
                or temp_c > L.t_max_c * L.warn_frac)
        if self.state == State.FAULT:
            return self.state
        if any_fault:
            self.state = State.FAULT
        elif near:
            self.state = State.WARNING
        elif abs(current_a) > 1e-3:
            self.state = State.ACTIVE
        else:
            self.state = State.IDLE
        return self.state

    def reset(self) -> bool:
        """Clear a latched fault if all flags have released. Returns True on success."""
        if self.state == State.FAULT and not any(self.flags.values()):
            self.state = State.IDLE
            return True
        return self.state != State.FAULT
