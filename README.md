# dengue_ode.py — README

MATH 574 Team 2 · Task D · Xindi Li

## Interface (agreed with E1 / E2)
 
```python
dengue_rhs(t, y, params)   # use: solve_ivp(dengue_rhs, (0, t_end), y0, args=(params,))
```
 
- **State order of `y`:** `[S_h, E_h, I_h, R_h, S_v, E_v, I_v]`
- **Units:** time in days; compartments are absolute counts.
- **`params` keys:** `N_h, m, a0, kappa, eta, b, c, eps_h, eps_v, gamma_h, d_h, d_v`
- There is no `"a"` key. The model computes `a = a0 * (1 - kappa * eta)` internally.
- **Weekly output:** `weekly_reported = phi * (7-day sum of eps_h * E_h)`. `phi` is the reporting fraction; ρ is reserved for ρ(K).
## Parameter values
 
Baselines come from the AB package: `data_and_literature/parameters/parameters.csv` (values, ranges, sources).
Ranges are for sensitivity analysis only.
 
| Name | Baseline | Source / status |
|---|---|---|
| `N_h` | 2,347,833 | 2010 Census, San Juan metro (fixed) |
| `d_h` | 3.4781e-5 | 1 / 78.717 years, PAHO 2010 (fixed) |
| `a0, b, c` | 0.35, 0.33, 0.33 | literature (fixed, sensitivity) |
| `eps_h, eps_v, gamma_h, d_v` | 0.169492, 0.093458, 0.2, 0.071429 | literature (fixed) |
| `m` | 2.0 | **fitted by E1** (range 1–20; about 9.8 at s0 = 0.4) |
| `kappa`, `eta` | 0, 0 | 0 in calibration; scenario inputs for F |
| `phi` | 0.005 | **fitted by E1** (range 0.003–0.023, per infection) |
| `s0` | 0.4 | scenario (0.2 and 0.6 as alternatives) |
 
`validate_rhs(dengue_rhs, PARAMS)` from `e2_observation/observation.py` returns `[]` (checked 9 Nov with both the old provisional and the AB parameters).
 
The demo in `__main__` still starts from 10 cases with s0 = 1 and m = 2 (R0 ≈ 1.03), so it shows almost no outbreak.
For real runs, use E2's `initial_state` / `expected_cases` and E1's fitted m and phi.