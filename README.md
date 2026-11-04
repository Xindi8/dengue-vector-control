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

## Hard-coded values (provisional)

| Name | Value | Replaced by |
|---|---|---|
| `N_h`, `d_h` | 2.35e6, 3.5e-5 | B01 |
| `a0, b, c, eps_h, eps_v, gamma_h, d_v` | 0.35, 0.33, 0.33, 1/5.9, 1/10.7, 1/5, 1/14 | C01 |
| `m` | 2.0 | E1 (fit) |
| `kappa`, `eta` | 0, 0 (no intervention) | scenario input |
| `I_h0`, `R_h0` | 10, 0 | E1 (fit) |
| `E_h0, E_v0, I_v0` | 0 (inside `y0`) | E1, if needed |
| `phi` | 0.005 (range 0.003–0.023, per infection, based on Fahad's E2 plan)  | E1 (fit) |
| `n_weeks` | 52 | length of the data series |

With the current values R0 ≈ 1.03, so almost no outbreak develops. This is expected until `m` is fitted.
