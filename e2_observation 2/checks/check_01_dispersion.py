"""Check 1 - which residual should E1 fit on?  (model-free, uses no fitted model)

Idea: if the weekly mean changes smoothly, the second difference
e_w = y_w - (y_{w-1} + y_{w+1}) / 2 is mostly noise, with variance ~1.5 * Var(y_w).
Comparing that variance with the local mean tells us how noisy the counts are:
Var = mean (Poisson) or Var = mean + mean^2 / k (negative binomial, overdispersed).
Then we check which transform makes the noise the same size at low and high counts.
"""
import numpy as np
import pandas as pd

from common import DATA, OUT, obs, load

cal, _ = load()
full = pd.read_csv(f"{DATA}/01_san_juan_weekly_cases_1990_2013.csv", parse_dates=["week_start_date"])

series = {
    "DENV-1 attributed, 2010/11": cal.denv1_attributed_cases.to_numpy(),
    "DENV-1 typed, 2010/11": cal.denv1_cases.to_numpy(float),
    "All lab-positive, 2010/11": cal.total_cases.to_numpy(float),
}

rows = []
for name, y in series.items():
    mu, v = obs.local_dispersion(y, min_mean=5)
    k = obs.estimate_nb_size(mu, v)
    rows.append(dict(series=name, weeks=len(mu), k=k, var_to_mean_median=float(np.median(v / mu))))

# Pool all 23 seasons, all serotypes and DENV-1, for a steadier estimate
pool_mu, pool_v = [], []
for _, g in full.groupby("season"):
    for col in ("total_cases", "denv1_cases"):
        m_, v_ = obs.local_dispersion(g[col].to_numpy(float), min_mean=5)
        pool_mu.append(m_)
        pool_v.append(v_)
pool_mu, pool_v = np.concatenate(pool_mu), np.concatenate(pool_v)
print(f"Distinct weeks in file 01: {len(full)}; usable week-points (local mean >= 5) "
      f"from the all-serotype and DENV-1 series together: {len(pool_mu)}")
# Lab-positive counts only (drop the source's estimate for untested specimens)
lab_mu, lab_v = [], []
for _, g in full.groupby("season"):
    m_, v_ = obs.local_dispersion((g.total_cases - g.estimated_untested_cases).to_numpy(float), min_mean=5)
    lab_mu.append(m_)
    lab_v.append(v_)
print(f"Lab-positive only (no estimated cases): k = {obs.estimate_nb_size(np.concatenate(lab_mu), np.concatenate(lab_v)):.0f}")
k_pool = obs.estimate_nb_size(pool_mu, pool_v)
rows.append(dict(series="All 23 seasons (total + DENV-1)", weeks=len(pool_mu), k=k_pool,
                 var_to_mean_median=float(np.median(pool_v / pool_mu))))
tab = pd.DataFrame(rows)
print(tab.round(3).to_string(index=False))

# Binned mean-variance table (pooled), for the plan's chart
bins = [5, 10, 20, 40, 80, 160, 400]
cut = pd.cut(pool_mu, bins)
binned = pd.DataFrame({"mu": pool_mu, "v": pool_v, "bin": cut}).groupby("bin", observed=True).agg(
    weeks=("mu", "size"), mean=("mu", "mean"), var=("v", "mean"))
binned["poisson_var"] = binned["mean"]
binned["nb_var"] = binned["mean"] + binned["mean"] ** 2 / k_pool
binned["var_to_mean"] = binned["var"] / binned["mean"]
print("\nPooled mean-variance by level (Var estimated from second differences):")
print(binned.round(2).to_string())
binned.round(3).to_csv(f"{OUT}/check01_mean_variance_bins.csv")

# Which transform evens out the noise?  Variance of the transformed second
# difference, low-count weeks versus high-count weeks (closer to 1 = better).
def transformed_noise(y, g):
    gy = g(np.asarray(y, float))
    e = gy[1:-1] - 0.5 * (gy[:-2] + gy[2:])
    m = (y[:-2] + y[1:-1] + y[2:]) / 3.0
    return m, e ** 2 / 1.5

transforms = {
    "raw": lambda x: x,
    "sqrt": np.sqrt,
    "log1p": np.log1p,
    f"nb_vst (k={k_pool:.0f})": lambda x: obs.nb_vst(x, k_pool),
}
out = []
for tname, g in transforms.items():
    ms, vs = [], []
    for _, gdf in full.groupby("season"):
        for col in ("total_cases", "denv1_cases"):
            y = gdf[col].to_numpy(float)
            m_, v_ = transformed_noise(y, g)
            ms.append(m_)
            vs.append(v_)
    ms, vs = np.concatenate(ms), np.concatenate(vs)
    low = vs[(ms >= 5) & (ms < 20)].mean()
    high = vs[ms >= 80].mean()
    out.append(dict(transform=tname, var_low_5_20=low, var_high_80plus=high, high_over_low=high / low))
tt = pd.DataFrame(out)
print("\nNoise size after each transform (pooled seasons):")
print(tt.round(3).to_string(index=False))
tab.round(3).to_csv(f"{OUT}/check01_dispersion_summary.csv", index=False)
tt.round(3).to_csv(f"{OUT}/check01_transform_balance.csv", index=False)
