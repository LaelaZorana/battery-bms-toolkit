# battery-bms-toolkit

You'd think battery management math needs a vendor toolbox, but everything a BMS actually computes fits in plain numpy and scipy, and this repo proves it. Every algorithm here is small enough to read in one sitting, and every one is covered by a test that checks it against a known answer. No black boxes.

## Theory summary

The cell is an equivalent circuit: an OCV source in series with an ohmic resistance R0 and one or two RC pairs, which gives you the Thevenin 1RC and 2RC models. Terminal voltage is V = OCV of SOC minus I R0 minus the sum of RC voltages, with current positive on discharge, and each RC voltage follows a first order recursion with time constant R C. The OCV curve is an interpolated table representative of an NMC cell running 3.0 to 4.2 V. Resistances pass through an Arrhenius style temperature scaling hook, and OCV takes an optional dV/dT term when you care about entropic effects.

SOC estimation is where the interesting failures live. Coulomb counting just integrates current, so any sensor bias makes it drift without bound. The Extended Kalman Filter uses the 1RC model as the process and the terminal voltage as the measurement, linearising OCV against SOC with a numerical slope, while the sigma point unscented variant propagates 2n+1 sigma points through the nonlinear OCV and never needs a Jacobian. Both start with a 20 percent SOC error on purpose and pull to truth from noisy voltage.

SOH tracking rides on a simple identity: between two rest points the charge moved equals Q times the SOC change, so a recursive least squares with a forgetting factor tracks the usable capacity Q as it fades.

Parameter identification fits an HPPC style current pulse followed by rest with nonlinear least squares over R0, R1 and C1. The seed is closed form. The instantaneous drop gives R0, the slow drop gives R1, and a 63 percent relaxation time gives the time constant.

Balancing is passive: every cell more than one tolerance above the lowest cell bleeds through a resistor, and the simulator reports how long it takes and how much energy turns into heat.

Thermal is a single lumped mass with I squared R heating, an optional entropic term I T dU/dT, convection to ambient, and a runaway threshold flag.

Safety checks cover over voltage, under voltage, over current and over temperature with hysteresis, and they feed a four state machine running idle, active, warning and fault. Faults latch until the flags release and reset() is called. Which is how a real BMS behaves.

## API

| Module | Main entry points |
|---|---|
| `bms.ecm` | `ECMParams`, `ECM.step(i, dt)`, `ECM.simulate(i, dt, soc0)`, `ocv(soc, temp_c, dvdt)`, `coulomb_count`, `synthetic_drive_cycle` |
| `bms.estimation` | `SOCEKF.update(i, v)`, `SOCUKF.update(i, v)`, `CapacityRLS.update(delta_soc, ah)` |
| `bms.param_id` | `hppc_pulse`, `initial_guess_from_pulse`, `identify_1rc(i, v, dt, soc0, capacity_ah)` |
| `bms.balancing` | `simulate_passive_balancing(soc0, capacity_ah, r_bleed, tol)` |
| `bms.thermal` | `ThermalParams`, `LumpedThermal.step(i, r_total, dt)` |
| `bms.safety` | `Limits`, `SafetyMonitor.check(v_max, v_min, i, temp)`, `SafetyMonitor.reset()` |

## Validation

| Check | Result |
|---|---|
| 1RC step response vs analytic exponential | max error 7e-8 V (limit 1e-6) |
| EKF from 20 percent initial SOC error, 5 mV noise, 30 min drive cycle | final error 0.26 percent (limit 2) |
| UKF, same scenario | within 2 percent |
| Coulomb counting with 0.15 A sensor bias, same run | 23 percent error at end |
| HPPC fit of R0 25 mOhm, R1 12 mOhm, C1 1500 F with 1 mV noise | 25.04 mOhm, 12.03 mOhm, 1505 F (all within 5 percent) |
| RLS capacity, true 2.2 Ah from 2.5 Ah start | within 2 percent after 50 cycles |
| Passive balancing of 4 cells at 80/85/90/95 percent, 33 ohm | balanced in 2.97 h, 2.86 Wh wasted, spread within 0.5 percent |
| Safety state machine | idle, active, warning, fault, latch, hysteresis, reset all exercised |

Figures produced by `examples/make_figures.py` live in `figures/`: `ocv_curve.png`, `ekf_soc.png`, `hppc_fit.png`, `balancing.png`, `pack_temperature.png`.

## Run

```
python -m pytest -q
python examples/make_figures.py
```

Requires numpy, scipy, matplotlib, pandas, pytest.
