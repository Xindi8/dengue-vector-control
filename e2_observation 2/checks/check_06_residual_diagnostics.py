"""Check 6 - residual diagnostics on the full-season fit (preview, provisional parameters).
Shows what the Stage 3 diagnostics look like and how big the model error is
compared with the reporting noise measured in check 1."""
import numpy as np
import pandas as pd

from common import OUT, dengue_rhs, fit_m_phi, load, obs

cal, edges = load()
y = cal.denv1_attributed_cases.to_numpy(float)
f = fit_m_phi(y, edges, 0.4, m0=10, phi0=0.005)
mu = obs.expected_cases(f["params"], 0.4, f["phi"], edges, rhs=dengue_rhs, y_anchor=y)
r = obs.residuals(y, mu, kind="sqrt")
k = 293.0
pear = obs.residuals(y, mu, kind="pearson", k=k)
out = dict(
    sqrt_residual_variance=float(np.var(r)),
    noise_only_expectation="about 0.30-0.50 (check 1)",
    lag1_autocorrelation=obs.lag1_autocorrelation(r),
    runs_test_z=obs.runs_test_z(r),
    pearson_dispersion=float(np.sum(pear ** 2) / (len(y) - 2)),
    weeks_model_too_low_1_20=int(np.sum(r[:20] > 0)),
    weeks_model_too_low_21_52=int(np.sum(r[20:] > 0)),
)
for key, val in out.items():
    print(f"{key:28s} {val if isinstance(val, str) else round(val, 3)}")
pd.DataFrame({"week": np.arange(1, 53), "observed": y.round(1), "fitted": mu.round(1),
              "sqrt_residual": r.round(3)}).to_csv(f"{OUT}/check06_full_fit_residuals.csv", index=False)
