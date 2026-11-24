# Dengue vector control in San Juan, Puerto Rico

MATH 574 (Mathematical Modeling of Infectious Diseases), University of Alberta, Fall 2026. Team 2 term project.

We fit a host–vector model (SEIR for people, SEI for *Aedes aegypti*) to the 2010/11 DENV-1 season in San Juan. Then we ask two questions. Which mix of personal protection, mosquito reduction and *Wolbachia* would have stopped that season, and which mix keeps R0 below 1? And does a cut in the biting rate do more than the same cut in mosquito numbers?

Everything in this repository runs from one command, and every number below comes from `results/SUMMARY.md`, which that command writes.

## Results at a glance

Headline s0 = 0.4 (share of people susceptible to DENV-1 on 30 Apr 2010), sensitivity band 0.2–0.6, protection efficacy η = 0.5.

**What the data identify.** One season pins down R_eff on 30 Apr 2010 = 1.44, not R0. Every s0 from 0.1 to 1 fits the 2010/11 curve equally well, because only m × s0 ≈ 3.9 and φ × s0 ≈ 0.0019 are identified. R0 therefore depends on s0: 3.22 at s0 = 0.2, 2.28 at 0.4, 1.86 at 0.6. We report two targets. R_eff < 1 means "this control would have stopped the 2010/11 season". R0 < 1 means "this control prevents outbreaks even in a fully susceptible population".

**RQ1: which controls bring R below 1?** No single lever meets both targets within its ceiling.

| Lever alone | Stop 2010/11 (R_eff < 1) | Keep R0 < 1, s0 = 0.4 | Keep R0 < 1, s0 = 0.2–0.6 | Ceiling |
|---|---|---|---|---|
| Personal protection, κ | 0.61 | impossible | 0.93 at s0 = 0.6; impossible for s0 ≤ 0.5 | 0.6 |
| Source reduction, δ | 0.52 | 0.81 | 0.71–0.90 | 0.85 |
| *Wolbachia* coverage, p_w | 0.58 | 0.90 | 0.79 at s0 = 0.6; impossible at s0 = 0.2 | 1 |

| Package (κ, δ, p_w) | R_eff | R0 at s0 = 0.2 / 0.4 / 0.6 | 2010/11 cases averted |
|---|---|---|---|
| S1 Realistic protection (0.35, 0, 0) | 1.19 | 2.65 / 1.88 / 1.54 | 32% |
| S2 Moderate spraying (0, 0.4, 0) | 1.12 | 2.49 / 1.76 / 1.44 | 51% |
| S3 Proposal's combination (0.35, 0.4, 0) | 0.92 | 2.05 / 1.46 / 1.19 | 86% |
| S4 Maximum feasible, no *Wolbachia* (0.6, 0.85, 0) | 0.39 | 0.87 / 0.62 / 0.51 | 99% |
| S5 Balanced (0.35, 0.75, 0) | 0.59 | 1.33 / 0.94 / 0.77 | 97% |
| S6 Everything at its ceiling (0.6, 0.85, 0.8) | 0.21 | 0.46 / 0.33 / 0.27 | 99% |

The proposal's combination would have stopped the 2010/11 season but leaves R0 above 1 for every s0 in the band. Without *Wolbachia*, only the maximum feasible package (S4) keeps R0 below 1 across the whole band. The balanced package (S5) does so at s0 = 0.4 and 0.6 but not at 0.2.

**RQ2: biting rate or mosquito numbers?** Percent for percent, cutting the biting rate does twice as much to R, because R scales with a but with √m. Over 52 weeks, 10%, 20% and 30% cuts to the biting rate avert 14%, 40% and 72% of cases; the same cuts to mosquito density avert 6%, 16% and 30%. The advantage at the peak is smaller: 1.9 times at a 10% cut and 1.4 times at 30%.

**How far to trust it (E2's validation).**
- Forecasts made after the peak is in the data are good. Fitted to week 16, the forecast for the rest of the season is 23% high; fitted to week 26, 8.5% low. Fitted to week 12, before the peak, it is 237% too high.
- The model does not carry over between seasons. Run forward from the 2010/11 fit, it predicts R_eff ≈ 0.64–0.85 in April 2012 and so no 2012/13 epidemic, but 4,091 DENV-1-attributed cases followed. The model uses up susceptibles too fast because it has no seasonality, importation or heterogeneity.
- So the results above compare controls within one season. They are not multi-year projections.

## How to run

Python 3.10 or later.

```bash
python -m venv .venv        # optional: a separate environment for the project
source .venv/bin/activate   # on Windows: .venv\Scripts\activate
pip install -r requirements.txt
python run_all.py           # about 45 s: data checks, model checks, calibration, validation, interventions, figures
python -m pytest            # 44 tests, about 6 s
```

`python run_all.py --skip-validation` leaves out E2's seven validation checks and takes about 10 s. `python run_all.py --rebuild-data` also rebuilds the data package from its raw files.

Tested on Python 3.10 (numpy 2.2, pandas 2.3, scipy 1.15) and Python 3.13 (numpy 2.5, pandas 3.1, scipy 1.18): the tests pass on both and `results/SUMMARY.md` comes out identical.

Every shared setting is in `dengue/config.py`. To reproduce F's documents exactly (s0 = 0.5, band 0.39–0.64), set `S0_HEADLINE = 0.5` and `S0_BAND = (0.39, 0.64)` and rerun.

## What is where

```
run_all.py                 the whole pipeline; writes results/ and figures/
dengue/                    the project as one Python package
    config.py              every shared setting (s0, bounds, η, q, ceilings, packages)
    model.py               D, C   7-state model, disease-free equilibrium, R0 (next-generation and closed form)
    data.py                AB     loaders for data_and_literature/
    observation.py         E2     weekly expected cases, starting state, residuals, holdout tools
    diagnostics.py         E2     fit diagnostics (RSS, AIC, standard errors, correlations, residual tests)
    calibration.py         E1     fits of m and φ for each fixed s0
    interventions.py       F      controls, packages, single-lever thresholds, equal-cut comparison
    figures.py                    the four figures
validation/                E2     seven checks (dispersion, ridge, synthetic recovery, holdout,
                                  cross-season and 2024, residuals, diagnostics); see validation/README.md
tests/                     44 tests, including regression tests against E1's and F's own results
data_and_literature/       AB     data, parameter table, literature evidence, checks, build script
results/                   tables and SUMMARY.md (written by run_all.py)
figures/                   fig1 calibration, fig2a–c R = 1 contours, fig3 equal cuts, fig4 package curves
docs/                      each role's documents, and INTEGRATION.md (how the parts were combined)
archive/                   superseded files, kept for the record
```

## The model

State y = (S_h, E_h, I_h, R_h, S_v, E_v, I_v), time in days, N_v = m N_h, biting rate a = a0 (1 − κη).

```math
\begin{aligned}
S_h' &= d_h N_h - a b \tfrac{I_v}{N_h} S_h - d_h S_h, &
S_v' &= d_v N_v - a c \tfrac{I_h}{N_h} S_v - d_v S_v,\\
E_h' &= a b \tfrac{I_v}{N_h} S_h - (\varepsilon_h + d_h) E_h, &
E_v' &= a c \tfrac{I_h}{N_h} S_v - (\varepsilon_v + d_v) E_v,\\
I_h' &= \varepsilon_h E_h - (\gamma_h + d_h) I_h, &
I_v' &= \varepsilon_v E_v - d_v I_v,\\
R_h' &= \gamma_h I_h - d_h R_h.
\end{aligned}
```

```math
R_0 = \sqrt{\frac{a^2 b c\, m\, \varepsilon_h \varepsilon_v}{(d_h+\varepsilon_h)(d_h+\gamma_h)(d_v+\varepsilon_v)\, d_v}}, \qquad R_\text{eff} = R_0 \sqrt{s_0}.
```

Reported cases in week w are φ times the new symptomatic infections (ε_h E_h) in that week. The fit starts from E2's growth-mode state on 30 Apr 2010, scaled to weeks 1–3, and uses square-root residuals.

The three controls each change one parameter: personal protection a → a0 (1 − κη); source reduction (fewer adult mosquitoes per person) m → m (1 − δ); *Wolbachia* b, c → each × √(1 − p_w q). Together, R → R (1 − κη) √(1 − δ) √(1 − p_w q).

Baseline parameter values, ranges and sources are in `data_and_literature/parameters/parameters.csv`; the model reads them from `parameters.json`, so the table is the single source.

## Team

| Role | Member | Work | Now in |
|---|---|---|---|
| AB, data and literature | Chaya, Fahad | Chaya: first-pass parameter table and cleaned DENV-1 series. Fahad: data and literature package | `data_and_literature/`, `dengue/data.py` |
| C, theory | Aileen | Equations, disease-free equilibrium, next-generation R0 | `docs/C_theory_derivation_v0.1.pdf`, `dengue/model.py` |
| D, model code | Xindi | Shared ODE and R0 code, parameter interface | `dengue/model.py` |
| E1, calibration | Zaid | Fits of m and φ at fixed s0 | `dengue/calibration.py` |
| E2, observation and validation | Fahad | Observation model, diagnostics, validation checks | `dengue/observation.py`, `dengue/diagnostics.py`, `validation/` |
| F, interventions | Deckar | Controls, packages, contours, answers to RQ1 and RQ2 | `dengue/interventions.py`, `dengue/figures.py` |
| Integration | Fahad | One package, pipeline, tests | `run_all.py`, `tests/`, `docs/INTEGRATION.md` |

`docs/INTEGRATION.md` lists what each part contributed, what changed when the parts were joined, and the decisions still open for the team.

## Assumptions and limitations

- s0 is not identified by the data. The headline 0.4 follows the AB evidence (about 0.29–0.44); the band 0.2–0.6 is the range E1 fitted. Answers about R_eff do not depend on s0; answers about R0 do.
- η = 0.5 is a middle value; the evidence runs from near 0 (field surveys) to about 0.9 (DEET in the lab). At η = 0.3 the proposal's combination leaves R_eff at 1.00.
- *Wolbachia* blocking q = 0.9 and coverage up to 0.8 are theoretical. The San Juan data contain no *Wolbachia* releases.
- Controls start on 30 Apr 2010, last all season and act at once. Source reduction removes adult mosquitoes directly, with no larval delay.
- Protection is averaged over everyone; users and non-users are not split.
- No seasonality, importation, age structure or second serotype. This is why the model cannot carry one season into the next.
- No cost data, so the "cheapest" package cannot be identified. S5 follows F's rule: keep κ at 0.35 and raise δ in steps of 0.05 until R0 < 1 at the headline s0.
- Week 35 of the season (24–31 Dec 2010) is 8 days long because the data's weeks restart on 1 January, so weekly curves show a small step there.

## Data

San Juan weekly dengue cases by serotype, 1990/91–2012/13, from the Puerto Rico Department of Health and CDC (the 2015 dengue forecasting challenge data); island-wide weekly cases 2014–2024 from OpenDengue v1.3; population from the US Census Bureau and PAHO. Sources, cleaning steps and checks are in `data_and_literature/README.md` and `data_and_literature/checks/`.
