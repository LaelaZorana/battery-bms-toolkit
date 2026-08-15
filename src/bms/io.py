"""Loaders for the cell parameter YAMLs, current profiles, and field logs in data/."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import yaml

from .ecm import ECMParams
from .safety import Limits
from .thermal import ThermalParams


@dataclass
class CellData:
    key: str
    name: str
    chemistry: str
    capacity_ah: float
    soc_pts: np.ndarray
    r0_ohm: np.ndarray
    r1_ohm: np.ndarray
    c1_f: np.ndarray
    ocv_soc: np.ndarray
    ocv_v: np.ndarray
    thermal: ThermalParams
    limits: Limits

    def ecm_params(self, soc: float = 0.5) -> ECMParams:
        """1RC parameters interpolated at the given SOC, carrying the cell's OCV table."""
        return ECMParams(
            capacity_ah=self.capacity_ah,
            r0=float(np.interp(soc, self.soc_pts, self.r0_ohm)),
            r1=float(np.interp(soc, self.soc_pts, self.r1_ohm)),
            c1=float(np.interp(soc, self.soc_pts, self.c1_f)),
            n_rc=1,
            ocv_soc=self.ocv_soc,
            ocv_v=self.ocv_v,
        )


def load_cell(path) -> CellData:
    path = Path(path)
    with open(path) as f:
        d = yaml.safe_load(f)
    rc = d["rc_vs_soc"]
    ocv = d["ocv_soc_table"]
    th = d["thermal"]
    lim = d["limits"]
    soc_pts = np.asarray(rc["soc"], dtype=float)
    ocv_soc = np.asarray(ocv["soc"], dtype=float)
    ocv_v = np.asarray(ocv["voltage_v"], dtype=float)
    if not (np.all(np.diff(soc_pts) > 0) and np.all(np.diff(ocv_soc) > 0)):
        raise ValueError(f"{path.name}: SOC axes must be strictly increasing")
    if not np.all(np.diff(ocv_v) >= 0):
        raise ValueError(f"{path.name}: OCV table must be nondecreasing in SOC")
    return CellData(
        key=d["key"],
        name=d["name"],
        chemistry=d["chemistry"],
        capacity_ah=float(d["capacity_ah"]),
        soc_pts=soc_pts,
        r0_ohm=np.asarray(rc["r0_ohm"], dtype=float),
        r1_ohm=np.asarray(rc["r1_ohm"], dtype=float),
        c1_f=np.asarray(rc["c1_f"], dtype=float),
        ocv_soc=ocv_soc,
        ocv_v=ocv_v,
        thermal=ThermalParams(
            mass_kg=float(th["mass_kg"]),
            cp_j_per_kgk=float(th["cp_j_per_kgk"]),
            h_w_per_k=float(th["h_w_per_k"]),
            t_ambient_c=float(th.get("t_ambient_c", 25.0)),
            t_shutdown_c=float(th.get("t_shutdown_c", 80.0)),
            dudt_v_per_k=float(th.get("dudt_v_per_k", 0.0)),
        ),
        limits=Limits(
            v_max=float(lim["v_max"]),
            v_min=float(lim["v_min"]),
            i_dis_max_a=float(lim["i_dis_max_a"]),
            i_chg_max_a=float(lim["i_chg_max_a"]),
            t_max_c=float(lim["t_max_c"]),
        ),
    )


def load_cells(dir_path) -> dict[str, CellData]:
    cells = [load_cell(p) for p in sorted(Path(dir_path).glob("*.yaml"))]
    return {c.key: c for c in cells}


def load_profile(path, cell_key: str):
    """Read a profile CSV with a time_s column and one current column per cell key.

    Returns (t, current_a) as float arrays. Positive current is discharge.
    """
    path = Path(path)
    with open(path) as f:
        header = f.readline().strip().split(",")
    if cell_key not in header:
        raise KeyError(f"{path.name} has no column {cell_key}, columns are {header}")
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    t = data[:, header.index("time_s")]
    i = data[:, header.index(cell_key)]
    return t, i


def load_fade_log(path):
    """Read a capacity fade log CSV with columns cycle, measured_ah, temp_c.

    Returns (cycle, measured_ah, temp_c) arrays.
    """
    data = np.genfromtxt(path, delimiter=",", skip_header=1)
    return data[:, 0].astype(int), data[:, 1], data[:, 2]
