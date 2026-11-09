"""Check 2 - what can the 2010/11 data pin down?  (feasibility check, NOT a calibration)

For a grid of fixed s0 values we fit m and phi by least squares on the square-root
scale, holding every other parameter at its Table 1 baseline. If the fit quality
is the same for every s0, then s0, m and phi cannot be separated by these data:
only the combinations m*s0 and phi*s0 are identified.

Then we ask whether the team's current bounds (m in [1, 10] from Table 1,
phi in [0.035, 0.09] from the proposal and the AB brief) can both hold.
"""
import numpy as np
import pandas as pd

from common import OUT, fit_m_phi, load

cal, edges = load()
S0_GRID = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
SERIES = ["denv1_attributed_cases", "denv1_cases", "total_cases"]

rows = []
for col in SERIES:
    y = cal[col].to_numpy(float)
    prev = None
    for s0 in S0_GRID:
        # warm start along the ridge from the previous s0
        m0, phi0 = (4.0 / s0, 0.002 / s0) if prev is None else (prev["m"] * prev["s0"] / s0, prev["phi"] * prev["s0"] / s0)
        f = fit_m_phi(y, edges, s0, m0=m0, phi0=phi0)
        prev = f
        rows.append(dict(series=col, s0=s0, m=f["m"], phi=f["phi"], cost=f["cost"], R0=f["R0"],
                         Reff=f["Reff"], m_s0=f["m"] * s0, phi_s0=f["phi"] * s0,
                         phi_over_m=f["phi"] / f["m"], status=f["status"]))
tab = pd.DataFrame(rows)
pd.set_option("display.width", 200)
print(tab.round(5).to_string(index=False))
tab.to_csv(f"{OUT}/check02_ridge_profile.csv", index=False)

print("\nSpread of cost along the ridge (max - min) per series:")
print(tab.groupby("series").cost.agg(lambda c: c.max() - c.min()).round(4))

# Are the team's bounds compatible?  Along the ridge phi/m is constant.
print("\nIdentified ratio phi/m along the ridge vs the box implied by the bounds:")
box_lo, box_hi = 0.035 / 10, 0.09 / 1
for col, g in tab.groupby("series"):
    ratio = g.phi_over_m.median()
    m_needed = (0.035 / ratio, 0.09 / ratio)
    phi_if_m10 = ratio * 10
    s0_if_m10 = g.m_s0.median() / 10
    print(f"  {col:24s} phi/m = {ratio:.2e} (box needs {box_lo:.1e} to {box_hi:.2e}); "
          f"phi in [0.035, 0.09] needs m = {m_needed[0]:.0f}-{m_needed[1]:.0f}; "
          f"at m = 10: phi = {phi_if_m10:.4f}, s0 = {s0_if_m10:.2f}")

# What happens if E1 uses the bounds exactly as written in the plan?
# Fit s0, m, phi together with m in [1, 10], phi in [0.035, 0.09], s0 in [0.01, 1].
from scipy.optimize import least_squares as _ls
from common import baseline_params, dengue_rhs, obs, r0

y = cal["denv1_attributed_cases"].to_numpy(float)
base = baseline_params()


def _resid(theta):
    s0, m, phi = theta
    try:
        mu = obs.expected_cases(dict(base, m=m), s0, phi, edges, rhs=dengue_rhs, y_anchor=y)
    except obs.ObservationError:
        return 1e3 * np.ones(len(y))
    return obs.residuals(y, mu, kind="sqrt")


best = None
for x0 in ([0.5, 5, 0.05], [0.05, 9, 0.04], [0.9, 2, 0.08], [0.02, 10, 0.035]):
    r = _ls(_resid, x0, bounds=([0.01, 1, 0.035], [1, 10, 0.09]), x_scale="jac")
    if best is None or r.cost < best.cost:
        best = r
s0, m, phi = best.x
mu = obs.expected_cases(dict(base, m=m), s0, phi, edges, rhs=dengue_rhs, y_anchor=y)
print(f"\nWith the plan's bounds: s0={s0:.3f}, m={m:.2f}, phi={phi:.4f}, cost={best.cost:.1f} "
      f"(free ridge cost ~16-18); model peak {mu.max():.0f} vs data {y.max():.0f}; "
      f"model total {mu.sum():.0f} vs data {y.sum():.0f}; "
      f"R0={r0(dict(base, m=m)):.2f}, Reff={r0(dict(base, m=m)) * np.sqrt(s0):.2f}")
pd.DataFrame({"week": np.arange(1, len(y) + 1), "observed": y.round(1),
              "fit_with_plan_bounds": mu.round(1)}).to_csv(f"{OUT}/check02_plan_bounds_fit.csv", index=False)
