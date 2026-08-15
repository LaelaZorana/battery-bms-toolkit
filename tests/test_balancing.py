import math

import numpy as np
import pytest
from bms.balancing import simulate_passive_balancing

def test_balancing_matches_readme_numbers():
    r = simulate_passive_balancing([0.80, 0.85, 0.90, 0.95], tol=0.005)
    assert r.converged
    assert r.soc[-1].max() - r.soc[-1].min() <= 0.005
    assert abs(r.time_to_balance_s / 3600.0 - 2.97) < 0.05
    assert abs(r.energy_wasted_wh - 2.86) < 0.05

def test_unconverged_time_is_nan():
    r = simulate_passive_balancing([0.1, 0.9], t_max=60.0)
    assert not r.converged
    assert math.isnan(r.time_to_balance_s)

def test_empty_pack_raises():
    with pytest.raises(ValueError):
        simulate_passive_balancing([])
