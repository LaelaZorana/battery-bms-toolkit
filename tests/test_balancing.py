import numpy as np
from bms.balancing import simulate_passive_balancing

def test_balancing_terminates_within_tol():
    r = simulate_passive_balancing([0.80, 0.85, 0.90, 0.95], tol=0.005)
    assert r.converged
    assert r.soc[-1].max() - r.soc[-1].min() <= 0.005
    assert r.energy_wasted_wh > 0
