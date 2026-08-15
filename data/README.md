# Data

These read like datasheet values, and they aren't. Everything in this folder is synthetic and representative, with parameter values sitting in the ranges published for each cell class, so the algorithms see numbers that behave the way real ones do. But no file here is a measurement of a real named product, so treat them as realistic inputs rather than as a datasheet.

## cells/

Three representative parameter sets, one YAML per cell, and each one carries capacity, R0, R1 and C1 versus SOC, its own OCV against SOC table, thermal mass and convection, and protection limits with separate charge and discharge current ceilings.

| File | Cell | Why it is here |
|---|---|---|
| `nmc-21700-5ah.yaml` | NMC811 21700 cylindrical, 4.85 Ah | The common EV energy cell, a steep OCV curve that makes voltage based SOC estimation easy |
| `lfp-prismatic-100ah.yaml` | LFP prismatic, 100 Ah | The storage workhorse, an OCV curve that is nearly flat from 20 to 90 percent SOC, which starves any voltage based estimator of observability |
| `nmc-pouch-power.yaml` | NMC622 pouch, 8 Ah | A high power cell with milliohm resistances and 10C discharge capability |

## profiles/

Current profiles as CSV, with one `time_s` column plus one current column per cell key, scaled to each cell's capacity. Positive current is discharge, and that convention carries through everything downstream, so it's the first thing to check when a result comes back with the wrong sign.

| File | Content |
|---|---|
| `drive_cycle_1hz.csv` | Synthetic 1 Hz EV drive cycle, 30 minutes, peaks near 1.5C per cell |
| `hppc_test.csv` | HPPC test, 10 s discharge pulse at 2C then 60 s rest, 10 Hz sampling |
| `fast_charge_cccv.csv` | Fast charge, constant current at the cell's charge limit, a derating step to 60 percent at 900 s emulating a thermal derate, then a CV style exponential taper |

## field/

| File | Content |
|---|---|
| `capacity_fade_nmc21700.csv` | Synthetic multi month cycling log for the NMC 21700 cell, columns cycle, measured_ah, temp_c, 600 cycles of square root plus linear fade with measurement noise |

Load all of it through `bms.io`: `load_cell`, `load_cells`, `load_profile`, `load_fade_log`.
