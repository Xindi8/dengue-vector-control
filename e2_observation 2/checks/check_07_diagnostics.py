"""Check 7 - the residual diagnosis module (RSS, AIC, parameter correlations) on two fits:
  A. s0 fixed at 0.4, m and phi free (the plan's recommendation)
  B. s0, m and phi all free (E1's original brief)
on synthetic data with a known answer and on the real 2010/11 season (preview, provisional
parameters). B should be flagged: s0, m and phi are not separately identifiable."""
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from common import OUT, baseline_params, dengue_rhs, load, obs
import diagnostics as dg

cal, edges = load()
y_real = cal.denv1_attributed_cases.to_numpy(float)
base = baseline_params()
truth = obs.expected_cases(dict(base, m=10.0), 0.4, 0.005, edges, rhs=dengue_rhs, y_anchor=y_real)
y_syn = np.random.default_rng(7).negative_binomial(290, 290 / (290 + truth)).astype(float)


def fit(y, free_s0):
    def resid(th):
        if free_s0:
            s0, m, phi = np.exp(th)
        else:
            s0, (m, phi) = 0.4, np.exp(th)
        try:
            mu = obs.expected_cases(dict(base, m=m), s0, phi, edges, rhs=dengue_rhs, y_anchor=y)
        except obs.ObservationError:
            return 1e3 * np.ones(len(y))
        return obs.residuals(y, mu, kind="sqrt")

    if free_s0:
        x0, lo, hi = np.log([0.4, 10, 0.005]), np.log([0.02, 0.1, 1e-4]), np.log([1.0, 500, 0.5])
        names = ["log_s0", "log_m", "log_phi"]
    else:
        x0, lo, hi = np.log([10, 0.005]), np.log([0.1, 1e-4]), np.log([500, 0.5])
        names = ["log_m", "log_phi"]
    res = least_squares(resid, x0, bounds=(lo, hi), x_scale="jac", diff_step=1e-4,
                        xtol=1e-10, ftol=1e-10)
    return res, dg.fit_report(res, param_names=names, fun=resid)


# What goes wrong with the default finite-difference step (res.jac used as returned)
res_def = least_squares(lambda th: obs.residuals(y_syn, obs.expected_cases(
    dict(base, m=float(np.exp(th[1]))), float(np.exp(th[0])), float(np.exp(th[2])), edges,
    rhs=dengue_rhs, y_anchor=y_syn), kind="sqrt"),
    np.log([0.4, 10, 0.005]), bounds=(np.log([0.02, 0.1, 1e-4]), np.log([1.0, 500, 0.5])),
    x_scale="jac", xtol=1e-10, ftol=1e-10)
rep_def = dg.fit_report(res_def, param_names=["log_s0", "log_m", "log_phi"])
print("Default diff_step, res.jac as returned (synthetic, all three free): condition number "
      f"{rep_def['condition_number']:.3g}, corr(s0, m) {rep_def['corr'][0, 1]:+.3f}, "
      f"flagged pairs {len(rep_def['flagged_pairs'])}")

rows = []
for label, y in (("synthetic", y_syn), ("real 2010/11 (preview)", y_real)):
    for free_s0 in (False, True):
        res, rep = fit(y, free_s0)
        tag = "B: s0, m, phi free" if free_s0 else "A: s0 fixed at 0.4"
        print(f"\n=== {label} | {tag} | estimates {np.round(np.exp(res.x), 4)}")
        print(dg.format_report(rep))
        corr = rep["corr"]
        rows.append(dict(data=label, fit=tag, p=rep["p"], rss=rep["rss"], aic=rep["aic"], aicc=rep["aicc"],
                         condition_number=rep["condition_number"],
                         corr_m_phi=corr[-2, -1],
                         corr_s0_m=corr[0, 1] if free_s0 else np.nan,
                         corr_s0_phi=corr[0, 2] if free_s0 else np.nan,
                         flagged=len(rep["flagged_pairs"])))
tab = pd.DataFrame(rows)
pd.set_option("display.width", 220)
print("\n" + tab.round(4).to_string(index=False))
tab.round(5).to_csv(f"{OUT}/check07_diagnostics.csv", index=False)
