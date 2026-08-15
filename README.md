# battery-bms-toolkit

Put the same 30 minute drive cycle through two SOC estimators with a 0.15 A current sensor bias sitting in the measurements, and coulomb counting ends 23 percent off while the Extended Kalman Filter, started 20 percent wrong on purpose, ends under 2 percent. One integrates the bias forever and the other keeps checking itself against terminal voltage, and closing that gap is most of what a battery management system is for. So this covers the equivalent circuit model, both SOC filters, parameter identification, balancing, thermal and safety logic, with every number below asserted by a test that checks it against a known answer. It all runs from parameter files for three representative cells rather than from constants buried in the code. Swap the cell file and everything downstream follows.

## Theory summary

The cell is an equivalent circuit: an OCV source in series with an ohmic resistance R0 and one or two RC pairs, which gives you the Thevenin 1RC and 2RC models. Terminal voltage is V equal to OCV of SOC minus I R0 minus the sum of RC voltages, with current positive on discharge, and each RC voltage follows a first order recursion with time constant R C. Each cell carries its own interpolated OCV table, and the default is a representative NMC curve running 3.0 to 4.2 V. Resistances pass through an Arrhenius style temperature scaling hook, but its default activation energy is illustrative and the estimators themselves assume 25 C, so don't read the temperature path as calibrated.

SOC estimation is where the interesting failures live. Coulomb counting looks like the honest option because it only integrates current, and that's exactly what breaks it, since any sensor bias makes it drift without bound. The Extended Kalman Filter uses the 1RC model as the process and the terminal voltage as the measurement, linearising OCV against SOC with a numerical slope, while the sigma point unscented variant propagates 2n+1 sigma points through the nonlinear OCV and never needs a Jacobian. Both start with a 20 percent SOC error on purpose and pull to truth from noisy voltage. Observability lives entirely in the OCV slope, which is why the same filter that nails a steep NMC curve struggles on a flat LFP one, and the examples show that difference on real parameter sets instead of asserting it.

SOH tracking rides on a simple identity: between two rest points the charge moved equals capacity times the SOC change. The signs are physical, because discharge moves positive amp hours while SOC falls, and a recursive least squares with a forgetting factor tracks the usable capacity as it fades.

Parameter identification fits an HPPC style current pulse followed by rest with nonlinear least squares over R0, R1 and C1. The closed form seed takes R0 from the instantaneous drop and the time constant from a 63 percent relaxation measured after the ohmic recovery jump, then divides the slow drop by its settling factor 1 minus exp of minus pulse time over tau. That factor is easy to forget. Leave it out and the R1 seed comes back biased low whenever the pulse is short against tau, by about two times for the default pulse. The fitter warns when the residual stays large or a parameter sits on a bound, which is what a wrong starting SOC looks like from the inside.

Balancing is passive: every cell more than one tolerance above the lowest cell bleeds through a resistor, and the simulator reports how long it takes and how much energy turns into heat. When the run doesn't converge the reported time is NaN rather than a number that looks valid.

Thermal is a single lumped mass with I squared R heating, an entropic term, convection to ambient, and a shutdown threshold flag. With discharge positive current the reversible term is minus I times T times dU/dT, so a cell with negative dU/dT heats on discharge, and a test pins that sign to a hand computed value. The update uses the exact exponential solution of the linear ODE, so it's stable for any step size.

Safety checks cover over voltage, under voltage, over current and over temperature with hysteresis, feeding a four state machine running idle, active, warning and fault, and faults latch until the flags release and a reset is requested, which is how a real BMS behaves. Charge and discharge carry separate current limits because real cells accept far less charge than they deliver, warning thresholds are absolute offsets from each limit rather than fractions of an arbitrary zero point, and any NaN input trips a fault. A safety monitor has to fail toward fault.

## Data

The `data/` folder holds parameter YAMLs for three representative cells, an NMC 21700 energy cell near 5 Ah, an LFP prismatic storage cell at 100 Ah with its characteristically flat OCV table, and a high power NMC pouch. Alongside those sit the current profiles, a synthetic 1 Hz EV drive cycle, an HPPC test and a CC CV fast charge with a derating step, each one scaled per cell, plus a multi month capacity fade log for the SOH tracker. It's all synthetic and representative, it's documented in `data/README.md`, and it loads through `bms.io`.

## API

| Module | Main entry points |
|---|---|
| `bms.ecm` | `ECMParams`, `ECM.step(i, dt)`, `ECM.simulate(i, dt, soc0)`, `ocv`, `coulomb_count`, `synthetic_drive_cycle` |
| `bms.estimation` | `SOCEKF.update(i, v)`, `SOCUKF.update(i, v)`, `CapacityRLS.update(delta_soc, ah)` |
| `bms.param_id` | `hppc_pulse`, `initial_guess_from_pulse`, `identify_1rc` |
| `bms.balancing` | `simulate_passive_balancing(soc0, capacity_ah, r_bleed, tol)` |
| `bms.thermal` | `ThermalParams`, `LumpedThermal.step(i, r_total, dt)` |
| `bms.safety` | `Limits`, `SafetyMonitor.check(v_max, v_min, i, temp)`, `SafetyMonitor.reset()` |
| `bms.io` | `load_cell`, `load_cells`, `load_profile`, `load_fade_log` |

## Validation

Every row below is asserted by a test, and the EKF rows state their full conditions including the deliberate current bias, because a convergence number that doesn't say what noise and bias it survived isn't telling you much.

| Check | Conditions | Result |
|---|---|---|
| 1RC step response vs analytic exponential | 300 s at 2 A | max error 7e-8 V, limit 1e-6 |
| EKF convergence | 20 percent initial SOC error, 5 mV noise, 30 min drive cycle | final error under 2 percent |
| EKF with sensor bias | same run plus 0.15 A current bias fed to the filter | final error under 2 percent |
| EKF under model mismatch | filter R0, R1, C1 and capacity each off by 15 percent, plus the 0.15 A bias | final error under 5 percent |
| UKF | same base scenario | final error under 2 percent |
| LFP vs NMC observability | identical estimator, noise and relative bias on both data cells | steady LFP error more than twice the NMC error |
| Coulomb counting with 0.15 A sensor bias | same run | 23 percent error at end |
| HPPC fit of R0 25 mOhm, R1 12 mOhm, C1 1500 F | 1 mV noise | all three within 5 percent |
| HPPC closed form seed | default 10 s pulse against an 18 s tau | R0 within 10 percent, R1 within 30 percent |
| RLS capacity with physical signs | true 2.2 Ah from a 2.5 Ah start, 50 partial cycles | within 2 percent |
| Passive balancing of 4 cells at 80, 85, 90, 95 percent | 33 ohm bleed, 2.5 Ah | 2.97 h, 2.86 Wh wasted |
| Entropic heat sign | 5 A discharge, dU/dT of minus 1e-4 V per K | heats by the hand computed 0.149 W |
| Safety state machine | trips, releases, latch, NaN inputs, asymmetric charge limit | all exercised |

The example run on the data cells adds a 10 mV voltage sensor offset to the EKF scenario, and the offset maps to SOC error through the inverse OCV slope: the NMC 21700 ends at 1.17 percent error while the LFP cell, whose mid range slope is about ten times smaller, ends at 2.07 percent and converges visibly more slowly. The filter is identical in both runs, so what you're watching there is the OCV curve rather than the code. That comparison is in `figures/ekf_soc.png`.

Figures produced by `examples/make_figures.py` live in `figures/`: `ocv_curve.png`, `ekf_soc.png`, `hppc_fit.png`, `balancing.png`, `pack_temperature.png`, `soh_fade.png`.

## Run

```
pip install -e ".[test]"
python3 -m pytest -q
python3 examples/make_figures.py
```

You'll need numpy, scipy, matplotlib and pyyaml, with pytest as a test extra.

## License

MIT, see LICENSE.
