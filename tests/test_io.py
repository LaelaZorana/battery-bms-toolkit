import pathlib

import numpy as np
import pytest
from bms.ecm import ECM, synthetic_drive_cycle
from bms.estimation import SOCEKF
from bms.io import load_cell, load_cells, load_fade_log, load_profile

DATA = pathlib.Path(__file__).resolve().parents[1] / "data"


def test_load_all_cells():
    cells = load_cells(DATA / "cells")
    assert set(cells) == {"nmc_21700_5ah", "lfp_prismatic_100ah", "nmc_pouch_power"}
    lfp = cells["lfp_prismatic_100ah"]
    assert lfp.capacity_ah == 100.0
    p = lfp.ecm_params(0.5)
    assert 0.0005 < p.r0 < 0.001
    assert p.ocv_soc is not None
    # OCV tables must be nondecreasing and within the cell's voltage window
    for c in cells.values():
        assert np.all(np.diff(c.ocv_v) >= 0)
        assert c.ocv_v[0] >= c.limits.v_min - 1e-9
        assert c.ocv_v[-1] <= c.limits.v_max + 1e-9


def test_ecm_params_interpolates_with_soc():
    nmc = load_cell(DATA / "cells" / "nmc-21700-5ah.yaml")
    assert nmc.ecm_params(0.05).r0 > nmc.ecm_params(0.6).r0


def test_load_profiles_scaled_per_cell():
    for name in ["drive_cycle_1hz.csv", "hppc_test.csv", "fast_charge_cccv.csv"]:
        t, i_lfp = load_profile(DATA / "profiles" / name, "lfp_prismatic_100ah")
        _, i_nmc = load_profile(DATA / "profiles" / name, "nmc_21700_5ah")
        assert len(t) == len(i_lfp) == len(i_nmc)
        assert np.max(np.abs(i_lfp)) > np.max(np.abs(i_nmc))  # scaled by capacity
    t, i = load_profile(DATA / "profiles" / "drive_cycle_1hz.csv", "nmc_21700_5ah")
    assert np.allclose(np.diff(t), 1.0)  # 1 Hz
    t, i = load_profile(DATA / "profiles" / "fast_charge_cccv.csv", "nmc_21700_5ah")
    assert np.all(i < 0)  # charge is negative
    with pytest.raises(KeyError):
        load_profile(DATA / "profiles" / "drive_cycle_1hz.csv", "nope")


def test_load_fade_log():
    cyc, ah, temp = load_fade_log(DATA / "field" / "capacity_fade_nmc21700.csv")
    assert cyc[0] == 0 and cyc[-1] == 600
    assert ah[0] > ah[-1] > 4.0  # fades but stays plausible
    assert np.all((temp > 10) & (temp < 45))


def test_lfp_flat_ocv_degrades_observability_vs_nmc():
    cells = load_cells(DATA / "cells")
    nmc = cells["nmc_21700_5ah"].ecm_params(0.5)
    lfp = cells["lfp_prismatic_100ah"].ecm_params(0.5)
    # The mid SOC OCV slope is the observability knob, and LFP's is far smaller.
    assert nmc.docv_dsoc(0.5) > 5 * lfp.docv_dsoc(0.5)
    # Same experiment on both chemistries: 20 percent initial error, 5 mV noise,
    # scaled 1C class cycles. The LFP estimate must end further from truth.
    errs = {}
    for key, p in [("nmc", nmc), ("lfp", lfp)]:
        dt = 1.0
        i = synthetic_drive_cycle(1800, dt, i_max=1.5 * p.capacity_ah, seed=3)
        v, s = ECM(p).simulate(i, dt, soc0=0.9)
        rng = np.random.default_rng(0)
        vn = v + 0.005 * rng.standard_normal(len(v))
        f = SOCEKF(p, dt, soc0=0.7)
        est = np.array([f.update(float(a), float(b)) for a, b in zip(i, vn)])
        errs[key] = float(np.mean(np.abs(est[-300:] - s[-300:])))
    assert errs["lfp"] > 2 * errs["nmc"]
    assert errs["nmc"] < 0.02
