# E2: observation model and validation tools (MATH 574, Team 2)

Owner: Fahad (E2). For D (shared ODE), E1 (fitting) and Stage 3 (validation).
Version of 31 Oct 2026. The plan itself is in the team doc "E2 Observation and Validation Plan".

## What is here

| File | What it does |
|---|---|
| `observation.py` | The E2 module: real week boundaries, starting state, weekly expected cases, residuals, likelihoods, holdout splits and scores, and the contract test for D's model |
| `diagnostics.py` | Residual diagnosis for E1's fits: RSS, AIC/AICc/BIC, standard errors, parameter correlation matrix, condition number, flags for parameters the data cannot separate |
| `model_provisional.py` | Mirrors D's `dengue_rhs(t, y, params)` (same equations, same keys a0/kappa/eta) so the tests run without D's file |
| `tests/` | 18 tests (run `python -m pytest -q` in this folder) |
| `checks/check_01_dispersion.py` | Model-free noise check: which residual to fit on |
| `checks/check_02_ridge.py` | What the 2010/11 data can pin down, and what happens with the plan's bounds |
| `checks/check_03_synthetic.py` | Start-up dip, and parameter recovery on fake data with a known answer |
| `checks/check_04_holdout.py` | Dry run of the holdout (synthetic data, plus a preview on the real season) |
| `checks/check_05_cross_season_and_2024.py` | 2012/13 cross-season check and 2024 shape comparison |
| `checks/check_06_residual_diagnostics.py` | Residual diagnostics on a full-season preview fit |
| `checks/check_07_diagnostics.py` | The diagnosis module on fits with s0 fixed vs s0, m, phi all free |
| `outputs/` | CSV results from the checks (used in the plan) |
| `data/` | Copies of package files 01, 02, 03 and 06 |

Requires Python 3.10+, numpy (1.x or 2.x), scipy, pandas, and pytest for the tests.
Each check runs in a few seconds: `cd checks && python check_02_ridge.py`.

## How D, E1 and E2 connect

```
D:   dengue_rhs(t, y, params) -> 7 derivatives
     state order  S_h, E_h, I_h, R_h, S_v, E_v, I_v
     param keys   N_h, m, a0, kappa, eta, b, c, eps_h, eps_v, gamma_h, d_h, d_v   (per day; D has no "a")
     contract     observation.validate_rhs(dengue_rhs, params) must return []

E2:  edges = week_edges(cal)                                   # real week boundaries, days
     mu    = expected_cases(params, s0, phi, edges, rhs=dengue_rhs, y_anchor=y_train)
     r     = residuals(y_train, mu[train_idx], kind="sqrt")    # data - model, square-root scale

E1:  least_squares on r with diff_step=1e-4; catch observation.ObservationError and
     return a large finite residual vector (and log the failure) instead of crashing.
     Then: rep = diagnostics.fit_report(res, param_names=[...], fun=resid)
```

The reporting fraction is called phi (as in D's README), because rho is the
spectral radius rho(K) in the theory derivation. To use D's model, import
`dengue_rhs` from `dengue_model.py` instead of `model_provisional.py`.

`expected_cases` returns expected reported cases per week:
phi x (new symptom onsets in that week), where onsets are integrated exactly
with an extra counter state. The starting state sits on the model's own growth
mode and is scaled so the first 3 weeks match the data, so the starting size is
not an extra parameter. Pass `onset_rate0=` instead of `y_anchor=` if E1 prefers
to fit the starting size.

## Two things E1 must know

1. Only m*s0 and phi*s0 are identified by one season of case data. Dividing s0
   by any factor and multiplying m and phi by the same factor leaves the weekly
   cases (almost) unchanged. Fix s0 at scenario values and fit m and phi.
2. phi*s0 is about 0.002 whatever a, b, c and m are, so phi in [0.035, 0.09]
   forces s0 <= 0.055 (almost nobody susceptible to DENV-1). Treat phi as the
   share of ALL infections that become a counted case (about 0.003-0.023 for
   DENV-1), not the share of symptomatic cases. With Table 1's a = 0.35, m <= 10
   also caps phi at about 0.005 (the cap scales with a^2*b*c).
