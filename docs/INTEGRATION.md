# How the parts were combined

Integration of the team repository (github.com/Xindi8/dengue-vector-control, state of 20 Nov 2026) into one project. Fahad, 21 Nov 2026.

## Summary

Before integration, each role worked in its own folder with its own copies of the shared files. As committed, two of the analyses did not run from a fresh clone: E1's `e1_fit.py` stopped with an `ImportError`, and F's scripts read from `../inputs/` and wrote to `../outputs/` and `../figures/`, which are not in the repository. The roles also used two different s0 settings.

Now the code is one Python package (`dengue/`) with one settings file (`dengue/config.py`) and one command (`python run_all.py`) that runs the data checks, model checks, calibration, validation, interventions and figures in about 45 seconds. 44 tests pass. The package reproduces E1's fits to the precision of the files they committed, and F's scenario results to within one part in a million (`tests/test_regression.py`).

## Review of each part

### AB, data and literature (Chaya, Fahad)

**Submitted.** Chaya: a parameter table and a cleaned DENV-1 series (4 Nov). Fahad: the `data_and_literature/` package (9 Nov): calibration and secondary seasons, the full 1990–2013 series, OpenDengue 2014–2024, the 2024 validation year, demography, temperature, the parameter table, the literature evidence, checks and a build script.

**Assessment.** Chaya's extraction is correct as far as it goes, and her parameter table was a reasonable first pass. Three things kept it from being used as the fitting input. The series counts only typed DENV-1 cases, so it leaves out the DENV-1 share of untyped positives and of the estimated untested cases. It covers calendar years 2010–2011 rather than the season that starts on 30 Apr 2010. And several values in the table were later replaced by sourced ones: c = 0.5–1.0 (now 0.33), a 70-year human lifespan (now 78.7 years, PAHO 2010), and a reporting fraction of 0.035–0.09, which E2 showed forces s0 below 0.06. The AB package is what D's parameters and E1's data already came from.

**In the project.** The package is unchanged apart from its build script, where two pandas function names were updated (pandas 4 removes the old ones) and one figure setting was fixed. `dengue/data.py` loads it, `run_all.py` checks the raw files against their checksums, and the model reads its baseline values from `parameters/parameters.json`. Chaya's two files are in `archive/`.

### C, theory (Aileen)

**Submitted.** Theory Derivation v0.1 (31 Oct): the equations, the effective biting rate, the disease-free equilibrium, and R0 from the next-generation matrix.

**Assessment.** Correct, and D's code follows it section by section. The next-generation matrix and the closed form agree to 1e-10 on random parameter sets (`tests/test_model.py`). Two things are not in the repository yet: the R_eff = R0 √s0 addition the team asked for, and the controlled-R formula that F cites as C's.

**In the project.** `docs/C_theory_derivation_v0.1.pdf`; the section numbers are cited in the comments of `dengue/model.py`.

### D, model code (Xindi)

**Submitted.** `dengue_model.py` and its README (last updated 9 Nov).

**Assessment.** Correct and clean. It matches C's equations, passes E2's contract test, and its two R0 routes agree. Two gaps. There is no `r0(params)` helper, but E1's script imports one, which is why `e1_fit.py` fails. And the parameter values are typed in by hand from the AB table; they match it today but would drift if the table changed. The demo in `__main__` starts from s0 = 1 and m = 2 (R0 ≈ 1.03) and shows almost no outbreak, as D's README says.

**In the project.** `dengue/model.py` keeps D's equations line for line. It reads the parameters from `parameters.json` and adds `r0()` (closed form) and `r_effective()`. `basic_reproduction_number()` still returns both R0 values.

### E1, calibration (Zaid)

**Submitted.** `dengue_calibration_E1.py` (2 Nov, a scaffold tested on synthetic data) and `e1_fit.py` (11 Nov, Stage 2) with its outputs.

**Assessment.** The Stage 2 fit follows E2's plan exactly. It fixes s0, fits m and φ on log scales to the square-root residuals, multi-starts along the m × s0 ridge, and reports E2's diagnostics. The results reproduce. As committed, it does not run: it imports `r0` from `dengue_model`, which D's file does not define, and it looks for the data under the Section 4 file name `02_calibration_san_juan_2010_11.csv`. The 2 Nov scaffold is superseded.

**In the project.** `dengue/calibration.py` is the Stage 2 method as functions (`fit_s0`, `calibrate`). The bounds, step size and starting points are in `config.py`. One run fits a grid of ten s0 values (0.1 to 1), so every table and figure uses the same fits; before, E1 fitted three values and F refitted three others in a separate script.

### E2, observation model and validation (Fahad)

**Submitted.** `e2_observation 2/`: the observation and diagnostics modules, a provisional model, seven checks and 18 tests.

**Assessment.** The tests pass and the checks run. They were written before D's model and the AB package were final, so they used `model_provisional.py` and the Section 4 file names.

**In the project.** The two modules are in the package unchanged apart from their imports and the observation module's docstring, which now describes the final interfaces. The checks in `validation/` read the AB package through `dengue.data` and run against D's model; `run_all.py` runs them and keeps their logs in `results/validation/`. The provisional model is in `archive/`.

### F, interventions (Deckar)

**Submitted.** The F folder: five scripts, private copies of D's model (identical) and E2's observation module (one comment line differs), fit and scenario tables, and three documents (intervention setups; decisions and assumptions; answers to RQ1 and RQ2). Earlier, `r0_contour_plots.py` (3 Nov prototype).

**Assessment.** The strongest part of the repository. The code is short and clear, the scenarios answer the research questions directly, and the results reproduce exactly once the file paths are fixed. Points to fix:

1. The scripts run only inside a folder layout (`../inputs`, `../outputs`, `../figures`) that is not in the repository.
2. The week boundaries are typed by hand. They are right (they equal E2's `week_edges`), but they would not follow a change in the data.
3. The headline s0 is 0.5 with a 0.39–0.64 band. That band is the range E2 derived early on under the m ≤ 10 bound, which AB later raised to 20. The other roles use 0.4 with 0.2–0.6, and the seroprevalence evidence puts s0 at 0.29–0.44, so 0.5 lies outside it.
4. S5's δ = 0.65 is set by hand for s0 = 0.5, so it is not "balanced" at any other s0.
5. The setups document writes the controlled reproduction number as R0 (1 − κη)(1 − δ)(1 − pφ). Because R0 scales with √m and √(bc), the mosquito and *Wolbachia* terms need square roots: R0 (1 − κη) √(1 − δ) √(1 − p q). The code and every number in the results already use the square roots; only the document is wrong.
6. The same document uses φ for the *Wolbachia* blocking, but φ is the reporting fraction everywhere else (the code uses q), and it names the packages E0–E6, which clash with the role names E1 and E2 (the results use S0–S6).
7. δ is called "larval spraying", "larval source reduction", "spraying" and "adult mosquito reduction" in different places. In the model it is a cut in adult female mosquitoes per person.

**In the project.** `dengue/interventions.py` keeps F's formulas, packages and table layouts. Paths come from the package, fits from the calibration, and week boundaries from the data. S5's δ comes from F's own rule (hold κ at 0.35 and raise δ in steps of 0.05 until R0 < 1 at the headline s0), which gives F's 0.65 at s0 = 0.5. `dengue/figures.py` follows F's designs; the R0 band on the contour plots is now drawn from the closed form, so its edges are smooth.

## Old file → new file

| Repository file | Now |
|---|---|
| `dengue_model.py` (D) | `dengue/model.py` |
| `README.md` (D's model interface) | `README.md` (project) and the docstring of `dengue/model.py` |
| `Theory _Derivation_v0.1.pdf` (C) | `docs/C_theory_derivation_v0.1.pdf` |
| `data_and_literature/` (AB) | `data_and_literature/` (build script updated, data unchanged) and `dengue/data.py` |
| `Chaya_Dengue_Parameter_Table.csv`, `Chaya_SanJuan_DENV1_Cleaned.csv` | `archive/` |
| `e1_fit.py` (E1, Stage 2) | `dengue/calibration.py` |
| `dengue_calibration_E1.py` (E1, scaffold) | `archive/E1_calibration_scaffold_2nov.py` |
| `outputs/e1_fit_results.csv` | `tests/reference/` and `results/calibration_fits.csv` |
| `outputs/e1_fit_curves.csv`, `e1_m_phi_R0_Reff.csv`, `e1_fit_plot.png` | regenerated as `results/calibration_curves.csv`, `results/calibration_fits.csv`, `figures/fig1_calibration_fits.png` |
| `e2_observation 2/observation.py`, `diagnostics.py` | `dengue/observation.py`, `dengue/diagnostics.py` |
| `e2_observation 2/checks/` | `validation/` |
| `e2_observation 2/README.md` | `validation/README.md` (the checks) and `README.md` (the interfaces) |
| `e2_observation 2/data/` (copies of Section 4 files 01, 02, 03, 06) | replaced by the AB package: the same weekly counts, with unrounded DENV-1 shares and the cleaned OpenDengue flags |
| `e2_observation 2/outputs/` | regenerated in `results/validation/` (the same 22 files) |
| `e2_observation 2/tests/` | `tests/test_observation.py`, `tests/test_diagnostics.py` |
| `e2_observation 2/model_provisional.py` | `archive/E2_model_provisional.py` |
| `E2_Observation_and_Validation_Plan.pdf` | `docs/E2_observation_and_validation_plan.pdf` |
| F `f_fit.py` | `dengue/calibration.py` (same method; F's s0 values are in the fit grid) |
| F `f_model.py`, `f_scenarios.py` | `dengue/interventions.py` |
| F `f_contour_plots.py`, `f_figures.py` | `dengue/figures.py` |
| F private copies of `dengue_model.py`, `observation.py` | removed (the package has one copy) |
| F result tables | `tests/reference/` (regression targets) |
| F `fit tables/e1_fit_*.csv` | removed (identical copies of E1's outputs) |
| F documents (three PDFs) | `docs/F_*.pdf` |
| `r0_contour_plots.py` (F prototype) | `archive/F_r0_contour_prototype_3nov.py` |

## Checks that each member's numbers survive

- `tests/test_regression.py` refits E1's s0 = 0.2, 0.4, 0.6 and F's 0.39, 0.5, 0.64 and compares every column of their saved tables (agreement to their six saved decimals). From F's own fits, it rebuilds F's package table, weekly curves and equal-cut table; they agree to within one part in a million on any computer (the ODE solver itself works to about 1e-8).
- `tests/test_model.py`: the parameters equal the AB table; the disease-free equilibrium is a fixed point; populations are conserved; the next-generation and closed-form R0 agree; R0 scales with a and √m.
- `tests/test_interventions.py`: F's closed-form controlled R equals R0 of the changed model; F's hand-typed week boundaries equal the data's; the S5 rule gives 0.65 at s0 = 0.5; each single-lever threshold brings R to exactly 1.
- E2's 18 tests, now run against D's model and the AB data.

## What the project's settings change in F's answers

F's documents use s0 = 0.5 with a 0.39–0.64 band; the project uses 0.4 with 0.2–0.6. Answers about R_eff do not move. Answers about R0 do.

| Result | F's documents | This project |
|---|---|---|
| R_eff on 30 Apr 2010 | 1.44 | 1.44 |
| Single lever to stop 2010/11: κ, δ, p_w | 0.61, 0.52, 0.58 | 0.61, 0.52, 0.58 |
| S3 proposal's combination: R_eff, cases averted | 0.92, 86% | 0.92, 86% |
| RQ2 cases averted, 10/20/30% cut, a vs m | 15/40/72% vs 7/16/30% | 14/40/72% vs 6/16/30% |
| R0 at the headline s0 | 2.04 | 2.28 |
| δ alone for R0 < 1, headline (band) | 0.76 (0.69–0.81) | 0.81 (0.71–0.90) |
| p_w alone for R0 < 1, headline (band) | 0.84 (0.77–0.90) | 0.90 (0.79 to impossible) |
| S3: R0 across the band | 1.15–1.47 | 1.19–2.05 |
| S4 maximum feasible: R0 across the band | 0.49–0.63 | 0.51–0.87 |
| S5 balanced: δ, R0 across the band | 0.65, 0.88–1.13 | 0.75, 0.77–1.33 |

Two of F's conclusions change wording. "Larviciding alone needs more than 85%" was judged false at s0 = 0.39–0.64; it is true for s0 below about 0.31, which the project's band includes (0.90 at s0 = 0.2). And S4 is still the one package without *Wolbachia* that keeps R0 below 1 for every s0 in the band.

## Decisions for the team

1. **Headline s0.** The project uses 0.4 (AB's evidence, E1's fits, D's README). If the team prefers F's 0.5, change two lines in `config.py` and rerun; the table above shows what moves.
2. **F's setups document.** Correct the formula to the square-root form the code uses, rename the *Wolbachia* blocking from φ to q, and name the packages S0–S6.
3. **C's v0.2.** Add R_eff = R0 √s0 and the controlled-R formula, so the report can cite the theory for both.
4. **RQ1 wording.** With no cost data, "cheapest" cannot be answered (F's note). A version the results answer: "Which feasible levels of protection, mosquito reduction and *Wolbachia* bring R_eff and R0 below 1?"
5. **The proposal's third question** (whether dengue can persist when R0 < 1) has no code in the repository. The shared model has no importation or seasonality, so it cannot persist below threshold. Either drop the question from the report or add an importation term.
6. **One name for δ.** Suggested: "mosquito reduction", with source reduction and spraying as two ways to achieve it.
7. **η.** Report the main results at η = 0.5, with F's η table (0 to 0.9) as the sensitivity check.
