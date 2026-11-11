"""For s0 in {0.2, 0.4, 0.6}: fix s0, fit m and phi (log scale) to San Juan 2010/11
DENV-1-attributed weekly cases on the square-root residual scale, using
D's dengue_rhs and E2's expected_cases / fit_report.
"""
import os, sys, glob
import numpy as np, pandas as pd
# --- locate teammates' files (works from repo root or from this folder) ---
for d in ("e2_observation 2", "e2_observation", "."):
    if os.path.isdir(d): sys.path.insert(0, os.path.abspath(d))
_hits = glob.glob("**/02_calibration_san_juan_2010_11.csv", recursive=True)
DATA02 = _hits[0] if _hits else "data/02_calibration_san_juan_2010_11.csv"
os.makedirs("outputs", exist_ok=True)
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import observation as obs, diagnostics as dg
from dengue_model import dengue_rhs, PARAMS, r0

cal = obs.load_calibration(DATA02)
edges = obs.week_edges(cal)
y = cal["denv1_attributed_cases"].to_numpy(float)
base = dict(PARAMS)
assert obs.validate_rhs(dengue_rhs, base) == [], "contract test failed"

S0_LIST = [0.2, 0.4, 0.6]
LO, HI = np.log([0.05, 1e-5]), np.log([2000.0, 1.0])   # wide bounds on (m, phi)

def make_resid(s0):
    def resid(theta):
        p = dict(base, m=float(np.exp(theta[0])))
        try:
            mu = obs.expected_cases(p, s0, float(np.exp(theta[1])), edges, rhs=dengue_rhs, y_anchor=y)
        except (obs.ObservationError, ValueError) as e:
            return 1e3 * np.ones(len(y))
        return obs.residuals(y, mu, kind="sqrt")
    return resid

rows, curves, reports = [], {"week": np.arange(1, len(y)+1), "observed": y}, {}
for s0 in S0_LIST:
    resid = make_resid(s0)
    best = None
    for m0 in (4/s0, 2/s0, 8/s0):            # multi-start along the m*s0 ridge
        x0 = np.clip(np.log([m0, 0.002/s0]), LO+1e-9, HI-1e-9)
        r = least_squares(resid, x0, bounds=(LO, HI), x_scale="jac", diff_step=1e-4, xtol=1e-10, ftol=1e-10)
        if best is None or r.cost < best.cost: best = r
    m, phi = np.exp(best.x); p = dict(base, m=float(m))
    mu = obs.expected_cases(p, s0, phi, edges, rhs=dengue_rhs, y_anchor=y)
    rep = dg.fit_report(best, param_names=["log_m", "log_phi"], fun=resid)
    R0 = r0(p); Reff = R0 * np.sqrt(s0)
    rows.append(dict(s0=s0, m=m, phi=phi, R0=R0, R_eff=Reff, m_x_s0=m*s0, phi_x_s0=phi*s0,
                     RSS=rep["rss"], cost=best.cost, AIC=rep["aic"], AICc=rep["aicc"], BIC=rep["bic"],
                     cond_number=rep["condition_number"], corr_m_phi=rep["corr"][0,1],
                     se_log_m=rep["se"][0], se_log_phi=rep["se"][1],
                     lag1_autocorr=rep["lag1_autocorrelation"], runs_z=rep["runs_test_z"],
                     model_peak=mu.max(), data_peak=y.max(), model_total=mu.sum(), data_total=y.sum(),
                     peak_week_model=int(mu.argmax())+1, peak_week_data=int(y.argmax())+1,
                     solver_status=int(best.status)))
    curves[f"fit_s0_{s0}"] = mu
    reports[s0] = rep
    print(f"s0={s0}: m={m:.3f} phi={phi:.5f} R0={R0:.3f} Reff={Reff:.3f} RSS={rep['rss']:.3f} "
          f"cond={rep['condition_number']:.3g} status={best.status}")

res = pd.DataFrame(rows)
res.round(6).to_csv("outputs/e1_fit_results.csv", index=False)
res[["s0","m","phi","R0","R_eff"]].round(5).to_csv("outputs/e1_m_phi_R0_Reff.csv", index=False)
pd.DataFrame(curves).round(3).to_csv("outputs/e1_fit_curves.csv", index=False)

fig, ax = plt.subplots(figsize=(8,4.5))
ax.plot(curves["week"], y, "ko", ms=4, label="observed (DENV-1 attributed)")
for s0 in S0_LIST: ax.plot(curves["week"], curves[f"fit_s0_{s0}"], label=f"fit, s0={s0}")
ax.set_xlabel("epidemic week (from 30 Apr 2010)"); ax.set_ylabel("weekly cases"); ax.legend()
ax.set_title("San Juan 2010/11: ODE fit with s0 fixed, m and phi fitted"); fig.tight_layout()
fig.savefig("outputs/e1_fit_plot.png", dpi=150)
