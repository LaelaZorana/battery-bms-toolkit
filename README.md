# battery-bms-toolkit

Li-ion cell modeling and battery management system algorithms in plain numpy and scipy. Built as a portfolio piece: every algorithm is small, readable, and covered by a test that checks it against a known answer.

## Theory summary

Equivalent circuit model. The cell is an OCV source in series with an ohmic resistance R0 and one or two RC pairs (Thevenin 1RC and 2RC). Terminal voltage is V = OCV(SOC) minus I R0 minus the sum of RC voltages, with current positive on discharge. Each RC voltage follows a first order recursion with time constant R C. The OCV curve is an interpolated table representative of an NMC cell (3.0 to 4.2 V). Resistances pass through an Arrhenius style temperature scaling hook and OCV takes an optional dV/dT term.

SOC estimation. Coulomb counting integrates current and drifts with any sensor bias. The Extended Kalman Filter uses the 1RC model as the process and the terminal voltage as the measurement, linearising OCV(SOC) with a numerical slope. The sigma-point (unscented) variant propagates 2n+1 sigma points through the nonlinear OCV without a Jacobian. Both are initialised with a 20 percent SOC error and pull to truth from noisy voltage.

SOH tracking. Between two rest points the charge moved equals Q times the SOC change, so a recursive least squares with a forgetting factor tracks the usable capacity Q as it fades.

Parameter identification. An HPPC style current pulse followed by rest is fitted by nonlinear least squares over R0, R1, C1. A closed form seed uses the instantaneous drop for R0, the slow drop for R1, and a 63 percent relaxation time for the time constant.

Balancing. A passive scheme bleeds every cell that is more than one tolerance above the lowest cell through a resistor. The simulator reports time to balance and energy dissipated as heat.

Thermal. A single lumped mass with I squared R heating, an optional entropic term I T dU/dT, convection to ambient, and a runaway threshold flag.

Safety. Over voltage, under voltage, over current and over temperature checks with hysteresis feed a four state machine (idle, active, warning, fault). Faults latch until the flags release and reset() is called.

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
