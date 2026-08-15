from bms.safety import SafetyMonitor, State

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
