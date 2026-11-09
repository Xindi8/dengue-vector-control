"""Check 4 - dry run of the proposed holdout (rolling-origin, within season).

For each cutoff K we fit on weeks 1..K only, forecast weeks K+1..52, and score
the forecast with the CDC challenge targets (peak week, peak size, season total)
plus weekly errors. A naive 'persistence' forecast is the baseline to beat.

Part A runs on synthetic data (tests the code: the true model should forecast
well once it has seen the peak). Part B is a PREVIEW on the real 2010/11 DENV-1
series with provisional Table 1 parameters - not a calibration; Stage 3 repeats
it with the agreed model.
Also checks that forecasts do not depend on which s0 is fixed (the ridge).
"""
import numpy as np
import pandas as pd

from common import OUT, baseline_params, dengue_rhs, fit_m_phi, load, obs

cal, edges = load()
y_real = cal.denv1_attributed_cases.to_numpy(float)
n = len(y_real)
CUTOFFS = (8, 12, 16, 20, 26)


def forecast_after(y, K, s0=0.4):
    split = obs.Split(K, np.arange(K), np.arange(K, n))
    f = fit_m_phi(y, edges, s0, train=split.train, m0=4.0 / s0, phi0=0.002 / s0)
    mu = obs.expected_cases(f["params"], s0, f["phi"], edges, rhs=dengue_rhs, y_anchor=y)
    return split, f, mu


def run(y, label):
    rows, curves = [], {"week": np.arange(1, n + 1), "observed": y}
    for K in CUTOFFS:
        split, f, mu = forecast_after(y, K)
        s = obs.score_split(y, mu, split)
        base = np.concatenate([y[split.train], obs.persistence_forecast(y[split.train], len(split.test))])
        b = obs.score_split(y, base, split)
        rows.append(dict(data=label, cutoff=K, Reff_fit=round(f["Reff"], 3),
                         peak_week_err=s["peak_week_error"], peak_rel_err=s["peak_rel_error"],
                         test_total_rel_err=s["test_total_rel_error"],
                         season_total_rel_err=s["season_total_rel_error"],
                         rmse_sqrt=s["test_rmse_sqrt"],
                         baseline_peak_week_err=b["peak_week_error"],
                         baseline_test_total_rel_err=b["test_total_rel_error"],
                         baseline_rmse_sqrt=b["test_rmse_sqrt"],
                         test_weeks_observed_total=float(y[split.test].sum())))
        curves[f"forecast_from_week_{K}"] = np.where(np.arange(n) >= K, mu, np.nan)
    return pd.DataFrame(rows), pd.DataFrame(curves)


# Part A: synthetic truth (s0=0.4, m=10, phi=0.005), NB noise k=290
p_true = dict(baseline_params(), m=10.0)
mu_true = obs.expected_cases(p_true, 0.4, 0.005, edges, rhs=dengue_rhs, y_anchor=y_real)
rng = np.random.default_rng(11)
y_syn = rng.negative_binomial(290, 290 / (290 + mu_true)).astype(float)
A, curvesA = run(y_syn, "synthetic")

# Part B: real DENV-1 attributed series (preview)
B, curvesB = run(y_real, "real 2010/11 (preview)")

pd.set_option("display.width", 220)
out = pd.concat([A, B])
print(out.round(3).to_string(index=False))
out.round(4).to_csv(f"{OUT}/check04_holdout_scores.csv", index=False)
curvesB.round(2).to_csv(f"{OUT}/check04_real_forecast_curves.csv", index=False)
curvesA.round(2).to_csv(f"{OUT}/check04_synthetic_forecast_curves.csv", index=False)

# Forecasts along the ridge: fit week 1-12 with three different fixed s0
print("\nForecast from week 12 with different fixed s0 (real data):")
fc = {}
for s0 in (0.2, 0.4, 0.8):
    split, f, mu = forecast_after(y_real, 12, s0=s0)
    fc[s0] = mu
    print(f"  s0={s0}: m={f['m']:.2f} phi={f['phi']:.5f} R0={f['R0']:.2f} Reff={f['Reff']:.3f} "
          f"forecast peak week={int(np.argmax(mu)) + 1}, peak={mu.max():.1f}, season total={mu.sum():.0f}")
for K in (12, 16, 26):
    fk = {s0: forecast_after(y_real, K, s0=s0)[2] for s0 in (0.2, 0.4, 0.8)}
    test = np.arange(K, n)
    worst = max(np.max(np.abs(fk[a][test] / fk[0.4][test] - 1)) for a in (0.2, 0.8))
    tot = max(abs(fk[a][test].sum() / fk[0.4][test].sum() - 1) for a in (0.2, 0.8))
    print(f"  fit to week {K}: largest weekly difference between s0 = 0.2/0.4/0.8 forecasts "
          f"{worst:.1%}; test-window totals differ by at most {tot:.1%}")
