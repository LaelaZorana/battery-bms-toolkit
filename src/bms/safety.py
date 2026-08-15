"""Protection limits with hysteresis and a small BMS state machine.

States: idle, active, warning, fault. A fault latches until reset() is called and
all limits are back inside their release thresholds. Sign convention matches
ecm.py: current is positive on discharge, negative on charge, and the two
directions carry separate limits because real cells accept far less charge
current than they deliver. Any NaN input trips the corresponding flag, since a
safety monitor must fail toward fault, never toward silence.
"""
from __future__ import annotations

import math
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
    i_dis_max_a: float = 10.0        # discharge limit, positive current
    i_chg_max_a: float = 4.0         # charge limit, negative current, magnitude
    i_hyst_a: float = 1.0
    t_max_c: float = 60.0
    t_hyst_c: float = 5.0
    # Warning thresholds are absolute offsets from each limit. Fractions of an
    # absolute voltage are meaningless because the zero point is arbitrary.
    v_warn_margin: float = 0.05
    i_warn_margin_a: float = 0.5
    t_warn_margin_c: float = 3.0


def _bad(x: float) -> bool:
    return not math.isfinite(x)


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
        i_dis = max(current_a, 0.0)
        i_chg = max(-current_a, 0.0)
        self._hyst("ov", _bad(v_cell_max) or v_cell_max > L.v_max,
                   not _bad(v_cell_max) and v_cell_max < L.v_max - L.v_hyst)
        self._hyst("uv", _bad(v_cell_min) or v_cell_min < L.v_min,
                   not _bad(v_cell_min) and v_cell_min > L.v_min + L.v_hyst)
        oc = _bad(current_a) or i_dis > L.i_dis_max_a or i_chg > L.i_chg_max_a
        oc_rel = (not _bad(current_a) and i_dis < L.i_dis_max_a - L.i_hyst_a
                  and i_chg < L.i_chg_max_a - L.i_hyst_a)
        self._hyst("oc", oc, oc_rel)
        self._hyst("ot", _bad(temp_c) or temp_c > L.t_max_c,
                   not _bad(temp_c) and temp_c < L.t_max_c - L.t_hyst_c)
        any_fault = any(self.flags.values())
        near = (v_cell_max > L.v_max - L.v_warn_margin
                or v_cell_min < L.v_min + L.v_warn_margin
                or i_dis > L.i_dis_max_a - L.i_warn_margin_a
                or i_chg > L.i_chg_max_a - L.i_warn_margin_a
                or temp_c > L.t_max_c - L.t_warn_margin_c)
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
