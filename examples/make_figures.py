"""Generate all figures into figures/. Run from repo root."""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bms import (ECM, ECMParams, ocv, coulomb_count, synthetic_drive_cycle,
                 SOCEKF, hppc_pulse, identify_1rc, simulate_passive_balancing,
                 LumpedThermal, ThermalParams)
from bms.param_id import initial_guess_from_pulse
FIG = ROOT / "figures"; FIG.mkdir(exist_ok=True)

# 1 OCV curve
s = np.linspace(0, 1, 200)
plt.figure(figsize=(6, 4)); plt.plot(s * 100, ocv(s)); plt.xlabel("SOC (%)"); plt.ylabel("OCV (V)")
plt.title("Representative NMC OCV vs SOC"); plt.grid(alpha=.3); plt.tight_layout(); plt.savefig(FIG / "ocv_curve.png", dpi=130); plt.close()

# 2 EKF vs truth vs coulomb counting
p = ECMParams(); dt = 1.0
i = synthetic_drive_cycle(1800, dt, i_max=5.0, seed=1)
v, s_true = ECM(p).simulate(i, dt, soc0=0.9)
rng = np.random.default_rng(0)
vn = v + 0.005 * rng.standard_normal(len(v))
i_biased = i + 0.15   # sensor bias to show coulomb-count drift
cc = coulomb_count(i_biased, dt, p.capacity_ah, 0.7)
f = SOCEKF(p, dt, soc0=0.7)
ekf = np.array([f.update(float(a + 0.15), float(b)) for a, b in zip(i, vn)])
t = np.arange(len(i)) * dt
plt.figure(figsize=(7, 4)); plt.plot(t, s_true * 100, "k", label="truth")
plt.plot(t, cc * 100, "--", label="coulomb counting (biased sensor)"); plt.plot(t, ekf * 100, label="EKF")
plt.xlabel("time (s)"); plt.ylabel("SOC (%)"); plt.legend(); plt.grid(alpha=.3); plt.title("SOC estimation, 20% initial error")
plt.tight_layout(); plt.savefig(FIG / "ekf_soc.png", dpi=130); plt.close()
ekf_final_err = abs(ekf[-1] - s_true[-1]) * 100; cc_final_err = abs(cc[-1] - s_true[-1]) * 100

# 3 HPPC fit
true = ECMParams(capacity_ah=2.5, r0=0.025, r1=0.012, c1=1500.0)
ip, dtp = hppc_pulse(5.0, 10.0, 60.0, 0.1)
vp, _ = ECM(true).simulate(ip, dtp, 0.8)
vpn = vp + 0.001 * rng.standard_normal(len(vp))
r0, r1, c1 = identify_1rc(ip, vpn, dtp, 0.8, 2.5, x0=initial_guess_from_pulse(ip, vpn, dtp))
vfit, _ = ECM(ECMParams(capacity_ah=2.5, r0=r0, r1=r1, c1=c1)).simulate(ip, dtp, 0.8)
tp = np.arange(len(ip)) * dtp
plt.figure(figsize=(7, 4)); plt.plot(tp, vpn, ".", ms=2, label="measured (noisy)"); plt.plot(tp, vfit, label=f"fit R0={r0*1e3:.1f} mOhm R1={r1*1e3:.1f} mOhm C1={c1:.0f} F")
plt.xlabel("time (s)"); plt.ylabel("V"); plt.legend(); plt.grid(alpha=.3); plt.title("HPPC pulse fit"); plt.tight_layout(); plt.savefig(FIG / "hppc_fit.png", dpi=130); plt.close()

# 4 balancing
b = simulate_passive_balancing([0.80, 0.85, 0.90, 0.95])
plt.figure(figsize=(7, 4)); plt.plot(b.t / 3600, b.soc * 100); plt.xlabel("time (h)"); plt.ylabel("cell SOC (%)")
plt.title(f"Passive balancing, {b.time_to_balance_s/3600:.2f} h, {b.energy_wasted_wh:.2f} Wh wasted"); plt.grid(alpha=.3)
plt.tight_layout(); plt.savefig(FIG / "balancing.png", dpi=130); plt.close()

# 5 pack temperature under drive cycle
th = LumpedThermal(ThermalParams(dudt_v_per_k=-1e-4))
m = ECM(p); m.reset(0.9); temps = []
for ik in i:
    m.step(float(ik), dt); temps.append(th.step(float(ik), p.r0 + p.r1, dt))
plt.figure(figsize=(7, 4)); plt.plot(t, temps); plt.xlabel("time (s)"); plt.ylabel("cell temperature (C)")
plt.title("Lumped thermal response under drive cycle"); plt.grid(alpha=.3); plt.tight_layout(); plt.savefig(FIG / "pack_temperature.png", dpi=130); plt.close()

print(f"EKF final error {ekf_final_err:.2f} % SOC; coulomb count final error {cc_final_err:.2f} %")
print(f"HPPC fit: R0 {r0*1e3:.2f} mOhm (25), R1 {r1*1e3:.2f} mOhm (12), C1 {c1:.0f} F (1500)")
print(f"Balancing: {b.time_to_balance_s/3600:.2f} h, {b.energy_wasted_wh:.3f} Wh")
print(f"Peak temp {max(temps):.1f} C, runaway {th.runaway}")
