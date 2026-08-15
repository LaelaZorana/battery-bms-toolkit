import math

from bms.safety import Limits, SafetyMonitor, State


def test_state_machine_transitions():
    m = SafetyMonitor()
    assert m.check(3.7, 3.6, 0.0, 25) == State.IDLE
    assert m.check(3.7, 3.6, 3.0, 25) == State.ACTIVE
    assert m.check(3.7, 3.6, 9.7, 25) == State.WARNING
    assert m.check(4.25, 3.6, 3.0, 25) == State.FAULT
    # inside hysteresis band: flag stays, fault latched
    assert m.check(4.18, 3.6, 0.0, 25) == State.FAULT
    assert m.reset() is False
    m.check(4.10, 3.6, 0.0, 25)
    assert m.reset() is True
    assert m.state == State.IDLE
    assert m.check(3.7, 3.6, 0.0, 70) == State.FAULT


def test_healthy_full_cell_is_not_warning():
    m = SafetyMonitor()
    assert m.check(4.00, 3.6, 0.0, 25) == State.IDLE
    assert m.check(4.14, 3.6, 0.0, 25) == State.IDLE
    assert m.check(4.16, 3.6, 0.0, 25) == State.WARNING  # inside the 0.05 V margin


def test_charge_limit_is_tighter_than_discharge():
    m = SafetyMonitor()
    assert m.check(3.7, 3.6, 8.0, 25) == State.ACTIVE      # 8 A discharge fine
    assert m.check(3.7, 3.6, -8.0, 25) == State.FAULT      # 8 A charge exceeds 4 A limit


def test_oc_trip_and_release():
    m = SafetyMonitor()
    assert m.check(3.7, 3.6, 12.0, 25) == State.FAULT
    m.check(3.7, 3.6, 9.5, 25)                             # inside hysteresis, flag holds
    assert m.flags["oc"] is True
    assert m.reset() is False
    m.check(3.7, 3.6, 2.0, 25)
    assert m.flags["oc"] is False
    assert m.reset() is True


def test_uv_trip_and_release():
    m = SafetyMonitor()
    assert m.check(3.7, 2.9, 0.0, 25) == State.FAULT
    m.check(3.7, 3.02, 0.0, 25)                            # inside hysteresis
    assert m.flags["uv"] is True
    m.check(3.7, 3.2, 0.0, 25)
    assert m.flags["uv"] is False
    assert m.reset() is True


def test_reset_while_unfaulted_is_true():
    m = SafetyMonitor()
    m.check(3.7, 3.6, 0.0, 25)
    assert m.reset() is True
    assert m.state == State.IDLE


def test_nan_inputs_fault():
    nan = float("nan")
    for args in [(nan, 3.6, 0.0, 25), (3.7, nan, 0.0, 25), (3.7, 3.6, nan, 25), (3.7, 3.6, 0.0, nan)]:
        m = SafetyMonitor()
        assert m.check(*args) == State.FAULT
