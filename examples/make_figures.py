"""Generate all figures into figures/ from the parameter sets in data/. Run from repo root."""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bms import (ECM, coulomb_count, SOCEKF, identify_1rc, simulate_passive_balancing,
                 LumpedThermal, CapacityRLS, load_cells, load_profile, load_fade_log)
from bms.param_id import initial_guess_from_pulse

FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True)
DATA = ROOT / "data"
cells = load_cells(DATA / "cells")
nmc = cells["nmc_21700_5ah"]
lfp = cells["lfp_prismatic_100ah"]

# 1 OCV curves, NMC vs LFP
s = np.linspace(0, 1, 400)
plt.figure(figsize=(6.5, 4))
plt.plot(s * 100, nmc.ecm_params().ocv(s), label=nmc.name)
plt.plot(s * 100, lfp.ecm_params().ocv(s), label=lfp.name)
plt.xlabel("SOC (%)"); plt.ylabel("OCV (V)"); plt.legend(fontsize=8)
plt.title("OCV vs SOC, steep NMC against flat LFP"); plt.grid(alpha=.3)
plt.tight_layout(); plt.savefig(FIG / "ocv_curve.png", dpi=130); plt.close()

# 2 EKF on the drive cycle, NMC and LFP side by side. Same experiment on both:
# 20 percent initial SOC error, 5 mV voltage noise plus a 10 mV sensor offset,
# and 0.15 percent of capacity as current bias fed to the filter. The flat LFP
# curve starves the EKF of observability: a voltage offset maps to SOC error
# through the inverse OCV slope, and the LFP slope is about ten times smaller.
dt = 1.0
results = {}
for cell in (nmc, lfp):
    p = cell.ecm_params(0.6)
    t, i = load_profile(DATA / "profiles" / "drive_cycle_1hz.csv", cell.key)
    v, s_true = ECM(p).simulate(i, dt, soc0=0.9)
    rng = np.random.default_rng(0)
    vn = v + 0.005 * rng.standard_normal(len(v)) + 0.010
    bias = 0.0015 * p.capacity_ah
    f = SOCEKF(p, dt, soc0=0.7)
    ekf = np.array([f.update(float(a + bias), float(b)) for a, b in zip(i, vn)])
    cc = coulomb_count(i + bias, dt, p.capacity_ah, 0.7)
    results[cell.key] = (t, s_true, ekf, cc)
fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
for ax, cell in zip(axes, (nmc, lfp)):
    t, s_true, ekf, cc = results[cell.key]
    ax.plot(t, s_true * 100, "k", label="truth")
    ax.plot(t, cc * 100, "--", label="coulomb counting, biased sensor")
    ax.plot(t, ekf * 100, label="EKF")
    err = abs(ekf[-1] - s_true[-1]) * 100
    ax.set_title(f"{cell.chemistry}, final EKF error {err:.2f} %", fontsize=10)
    ax.set_xlabel("time (s)"); ax.grid(alpha=.3)
axes[0].set_ylabel("SOC (%)"); axes[0].legend(fontsize=8)
fig.suptitle("Same estimator, same noise and bias, two chemistries", fontsize=11)
fig.tight_layout(); fig.savefig(FIG / "ekf_soc.png", dpi=130); plt.close(fig)
nmc_err = abs(results[nmc.key][2][-1] - results[nmc.key][1][-1]) * 100
lfp_err = abs(results[lfp.key][2][-1] - results[lfp.key][1][-1]) * 100

# 3 HPPC fit on the NMC cell from the data profile
p = nmc.ecm_params(0.8)
tp, ip = load_profile(DATA / "profiles" / "hppc_test.csv", nmc.key)
dtp = float(tp[1] - tp[0])
vp, _ = ECM(p).simulate(ip, dtp, 0.8)
rng = np.random.default_rng(2)
vpn = vp + 0.001 * rng.standard_normal(len(vp))
r0, r1, c1 = identify_1rc(ip, vpn, dtp, 0.8, p.capacity_ah, x0=initial_guess_from_pulse(ip, vpn, dtp),
                          ocv_soc=nmc.ocv_soc, ocv_v=nmc.ocv_v)
pfit = nmc.ecm_params(0.8); pfit.r0, pfit.r1, pfit.c1 = r0, r1, c1
vfit, _ = ECM(pfit).simulate(ip, dtp, 0.8)
plt.figure(figsize=(7, 4))
plt.plot(tp, vpn, ".", ms=2, label="measured, noisy")
plt.plot(tp, vfit, label=f"fit R0={r0*1e3:.1f} mOhm R1={r1*1e3:.1f} mOhm C1={c1:.0f} F")
plt.xlabel("time (s)"); plt.ylabel("V"); plt.legend(fontsize=8); plt.grid(alpha=.3)
plt.title("HPPC pulse fit, NMC 21700"); plt.tight_layout()
plt.savefig(FIG / "hppc_fit.png", dpi=130); plt.close()

# 4 passive balancing of a four cell NMC pack
b = simulate_passive_balancing([0.80, 0.85, 0.90, 0.95], capacity_ah=nmc.capacity_ah)
plt.figure(figsize=(7, 4)); plt.plot(b.t / 3600, b.soc * 100)
plt.xlabel("time (h)"); plt.ylabel("cell SOC (%)")
plt.title(f"Passive balancing, {b.time_to_balance_s/3600:.2f} h, {b.energy_wasted_wh:.2f} Wh wasted")
plt.grid(alpha=.3); plt.tight_layout(); plt.savefig(FIG / "balancing.png", dpi=130); plt.close()

# 5 cell temperature under the drive cycle, NMC 21700
p = nmc.ecm_params(0.6)
t, i = load_profile(DATA / "profiles" / "drive_cycle_1hz.csv", nmc.key)
th = LumpedThermal(nmc.thermal)
m = ECM(p); m.reset(0.9); temps = []
for ik in i:
    m.step(float(ik), dt); temps.append(th.step(float(ik), p.r0 + p.r1, dt))
plt.figure(figsize=(7, 4)); plt.plot(t, temps)
plt.xlabel("time (s)"); plt.ylabel("cell temperature (C)")
plt.title("Lumped thermal response under drive cycle, NMC 21700"); plt.grid(alpha=.3)
plt.tight_layout(); plt.savefig(FIG / "pack_temperature.png", dpi=130); plt.close()

# 6 SOH tracking on the field fade log
cyc, ah, temp = load_fade_log(DATA / "field" / "capacity_fade_nmc21700.csv")
rls = CapacityRLS(nmc.capacity_ah, lam=0.95)
est = []
rng = np.random.default_rng(5)
for q_meas in ah:
    d = rng.uniform(0.3, 0.6)                       # partial cycle depth
    rls.update(-d, q_meas * d)                      # discharge: soc falls, Ah positive
    est.append(rls.capacity_ah)
est = np.array(est)
plt.figure(figsize=(7, 4))
plt.plot(cyc, ah, ".", ms=3, alpha=.6, label="measured capacity")
plt.plot(cyc, est, label="RLS estimate")
plt.xlabel("cycle"); plt.ylabel("capacity (Ah)"); plt.legend(fontsize=8); plt.grid(alpha=.3)
plt.title(f"SOH tracking on the fade log, final SOH {rls.soh(nmc.capacity_ah)*100:.1f} %")
plt.tight_layout(); plt.savefig(FIG / "soh_fade.png", dpi=130); plt.close()

print(f"EKF final error NMC {nmc_err:.2f} % SOC, LFP {lfp_err:.2f} % SOC")
print(f"HPPC fit: R0 {r0*1e3:.2f} mOhm (true {p.r0*1e3:.2f}), R1 {r1*1e3:.2f} mOhm, C1 {c1:.0f} F")
print(f"Balancing: {b.time_to_balance_s/3600:.2f} h, {b.energy_wasted_wh:.3f} Wh")
print(f"Peak temp {max(temps):.1f} C, shutdown {th.shutdown}")
print(f"RLS final capacity {rls.capacity_ah:.3f} Ah vs measured {ah[-1]:.3f} Ah, SOH {rls.soh(nmc.capacity_ah)*100:.1f} %")
