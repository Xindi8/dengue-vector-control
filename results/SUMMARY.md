# Key results

Written by `run_all.py`. Change settings in `dengue/config.py` and rerun.

- Headline s0 = 0.4; sensitivity band 0.2–0.6; η = 0.5, Wolbachia blocking q = 0.9.
- Calibration season: 52 weeks from 30 Apr 2010, 3,864 DENV-1-attributed cases.

## Calibration (E1)

- R_eff on 30 Apr 2010 = 1.44 at every s0 (range 1.437–1.459): m × s0 ≈ 3.94 and φ × s0 ≈ 0.0019.

| s0 | m | φ | R0 | RSS |
|---|---|---|---|---|
| 0.2 | 19.54 | 0.00942 | 3.22 | 35.53 |
| 0.39 | 10.06 | 0.00494 | 2.31 | 36.41 |
| 0.4 | 9.81 | 0.00482 | 2.28 | 36.43 |
| 0.5 | 7.87 | 0.00388 | 2.04 | 36.55 |
| 0.6 | 6.57 | 0.00324 | 1.86 | 36.59 |
| 0.64 | 6.17 | 0.00304 | 1.81 | 36.60 |

## Single levers needed to bring R below 1 (F)

| Lever | Ceiling | R_eff (stop 2010/11) | R0, s0 = 0.2 | R0, s0 = 0.4 | R0, s0 = 0.6 |
|---|---|---|---|---|---|
| Personal protection, kappa | 0.6 | 0.614 (above ceiling) | 1.378 (impossible) | 1.122 (impossible) | 0.927 (above ceiling) |
| Source reduction, delta | 0.85 | 0.520 | 0.903 (above ceiling) | 0.807 | 0.712 |
| Wolbachia coverage, p_w | 1 | 0.578 | 1.004 (impossible) | 0.897 | 0.791 |

R_eff stops the 2010/11 season as it started (the share susceptible on 30 Apr 2010). R0 prevents outbreaks in a fully susceptible population, so it depends on the s0 behind the fit.

## Packages (F), 52 weeks from 30 Apr 2010, s0 = 0.4

| ID | Package | κ | δ | p_w | R_eff | R0 at s0 = 0.2 / 0.4 / 0.6 | Cases | Averted | Peak (week) |
|---|---|---|---|---|---|---|---|---|---|
| S0 | Baseline | 0 | 0 | 0 | 1.44 | 3.22 / 2.28 / 1.86 | 3,821 | 0% | 203 (16) |
| S1 | Realistic protection | 0.35 | 0 | 0 | 1.19 | 2.65 / 1.88 / 1.54 | 2,580 | 32% | 79 (19) |
| S2 | Moderate spraying | 0 | 0.4 | 0 | 1.12 | 2.49 / 1.76 / 1.44 | 1,881 | 51% | 48 (20) |
| S3 | Proposal's combination | 0.35 | 0.4 | 0 | 0.92 | 2.05 / 1.46 / 1.19 | 526 | 86% | 30 (1) |
| S4 | Max feasible, no Wolbachia | 0.6 | 0.85 | 0 | 0.39 | 0.87 / 0.62 / 0.51 | 56 | 99% | 23 (1) |
| S5 | Balanced package | 0.35 | 0.75 | 0 | 0.59 | 1.33 / 0.94 / 0.77 | 100 | 97% | 25 (1) |
| S6 | Everything at its ceiling | 0.6 | 0.85 | 0.8 | 0.21 | 0.46 / 0.33 / 0.27 | 40 | 99% | 22 (1) |

## Equal cuts, biting rate vs mosquito density (RQ2), s0 = 0.4

| Cut | Lever | R_eff | Cases averted | Peak reduction |
|---|---|---|---|---|
| 10% | biting rate a | 1.30 | 14% | 36% |
| 10% | mosquito density m | 1.37 | 6% | 19% |
| 20% | biting rate a | 1.15 | 40% | 68% |
| 20% | mosquito density m | 1.29 | 16% | 39% |
| 30% | biting rate a | 1.01 | 72% | 83% |
| 30% | mosquito density m | 1.21 | 30% | 59% |
