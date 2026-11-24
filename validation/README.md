# Validation checks (role E2)

Seven checks of the observation model and the calibration, written for E2's validation plan (`docs/E2_observation_and_validation_plan.pdf`). They run on the shared model (`dengue/model.py`) and the AB data.

```bash
python run_all.py                     # stage 4 runs all seven; tables and logs go to results/validation/
cd validation && python check_04_holdout.py    # or run one check on its own
```

`common.py` holds the shared loaders and E2's fitting helper. Every check writes its tables to `results/validation/`, and `run_all.py` keeps each check's printed output there as a `.log` file.

| Check | Question | Main output |
|---|---|---|
| `check_01_dispersion.py` | How noisy are the weekly counts, and which residual scale should the fit use? (no model involved) | `check01_*.csv` |
| `check_02_ridge.py` | What can one season pin down, and what do the plan's first bounds on m and φ do to the fit? | `check02_ridge_profile.csv`, `check02_plan_bounds_fit.csv` |
| `check_03_synthetic.py` | Does the start-up dip go away with the growth-mode start, and are known parameters recovered from simulated data? | `check03_startup.csv`, `check03_synthetic_recovery.csv` |
| `check_04_holdout.py` | Fitted to week k only, how well does the model forecast the rest of the season, compared with a naive baseline? | `check04_holdout_scores.csv`, `check04_*_forecast_curves.csv` |
| `check_05_cross_season_and_2024.py` | Does the 2010/11 fit carry over to 2012/13, and how does the 2024 season compare in shape? | `check05_*.csv` |
| `check_06_residual_diagnostics.py` | Are the full-season residuals free of pattern? | `check06_full_fit_residuals.csv` |
| `check_07_diagnostics.py` | What do the fit diagnostics say with s0 fixed, and with s0, m and φ all free? | `check07_diagnostics.csv` |

## What they show

- **Only m × s0 and φ × s0 are identified.** With s0 fixed, the 2010/11 fit has condition number 1.7. With s0, m and φ all free, the condition number is 46 and the three parameters are correlated at about ±1 (check 7). So the calibration fixes s0 and fits m and φ.
- **Holdout.** Fitted up to week 16 or later, the model beats the baseline easily: the forecast for the remaining weeks is 23% high from week 16, 16% high from week 20, and 8.5% low from week 26, against 330–500% errors for the baseline. Fitted only up to week 12, before the peak, it is 237% too high (check 4).
- **Next season.** Fitted to 2012/13 on its own, the model gives R_eff = 1.30 and m × s0 = 3.19, against 1.44 and 3.92 for 2010/11. But the 2010/11 fit, run forward, leaves R_eff at 0.64–0.85 by April 2012 and so predicts no 2012/13 epidemic (check 5). The model uses up susceptibles too fast: it has no seasonality, importation or heterogeneity.

## Numbers that depend on the computer

Two fits are deliberately ill-conditioned, to show what goes wrong. Their exact numbers change with the numpy and scipy versions, so compare only their conclusions:

- the holdout fitted to week 8 (its forecast runs away: about 90 times the observed total on one setup, about 160 times on another);
- check 7's fit with scipy's default finite-difference step (condition number 12.6 on one setup, over 900 on another).

Everything else agrees within 1% between setups, and the project's own results (fits, packages, thresholds, equal cuts) within 0.02%. `results/SUMMARY.md` comes out identical on Python 3.10 with pandas 2.3 and on Python 3.13 with pandas 3.1.
