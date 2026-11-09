"""Check 5 - the other two validation options.

(a) Cross-season holdout, 2012/13 DENV-1 (preview). Because only m*s0 and
    phi*s0 are identified in each season, the ratio phi/m (= phi*s0 / m*s0)
    is identified too. If reporting (phi) and mosquito density (m) were stable
    between 2010/11 and 2012/13, the two seasons should give the same phi/m,
    whatever s0 is. That makes a ridge-proof consistency test.
(b) Shape comparison with 2024 (OpenDengue, island-wide, all serotypes), to
    support the comparability checklist.
"""
import numpy as np
import pandas as pd

from common import DATA, OUT, fit_m_phi, obs

full = pd.read_csv(f"{DATA}/01_san_juan_weekly_cases_1990_2013.csv", parse_dates=["week_start_date"])


def season_frame(season):
    g = full[full.season == season].copy().reset_index(drop=True)
    g["days_since_start"] = (g.week_start_date - g.week_start_date.iloc[0]).dt.days
    share = g.denv1_cases.sum() / g.typed_cases.sum()
    wk_share = np.where(g.typed_cases > 0, g.denv1_cases / g.typed_cases.where(g.typed_cases > 0, 1), share)
    g["denv1_attributed_cases"] = g.total_cases * wk_share
    return g


rows = []
for season in ("2010/2011", "2012/2013"):
    g = season_frame(season)
    edges = obs.week_edges(g)
    y = g.denv1_attributed_cases.to_numpy(float)
    f = fit_m_phi(y, edges, 0.4, m0=10, phi0=0.005)
    rows.append(dict(season=season, weeks=len(y), denv1_attributed_total=round(y.sum()),
                     Reff=f["Reff"], m_s0=f["m"] * 0.4, phi_s0=f["phi"] * 0.4,
                     phi_over_m=f["phi"] / f["m"], cost=f["cost"]))
cross = pd.DataFrame(rows)
print("(a) Identified combinations by season (s0 fixed at 0.4 only to run the fit):")
print(cross.round(5).to_string(index=False))
r = cross.phi_over_m.to_numpy()
print(f"    phi/m ratio 2012/13 vs 2010/11: {r[1] / r[0]:.2f}")
ms0 = cross.m_s0.to_numpy()
print(f"    m*s0 ratio 2012/13 vs 2010/11: {ms0[1] / ms0[0]:.2f}  "
      "(= s0 ratio if m stayed the same; compare with depletion predicted by the model)")
cross.round(6).to_csv(f"{OUT}/check05_cross_season.csv", index=False)

# Forward run: what the 2010/11 fit predicts for the start of 2012/13 (29 Apr 2012)
from scipy.integrate import solve_ivp
from common import dengue_rhs, r0
g10 = season_frame("2010/2011")
e10 = obs.week_edges(g10)
y10 = g10.denv1_attributed_cases.to_numpy(float)
t_end = (pd.Timestamp("2012-04-29") - g10.week_start_date.iloc[0]).days
fw = []
for s0 in (0.1, 0.4, 1.0):
    f = fit_m_phi(y10, e10, s0, m0=3.92 / s0, phi0=0.0019 / s0)
    p = f["params"]
    y0 = obs.initial_state(p, s0, obs.anchored_onset_rate(dengue_rhs, p, s0, f["phi"], e10, y10))
    sol = solve_ivp(lambda t, z: dengue_rhs(t, z, p), (0, t_end), y0, method="LSODA", rtol=1e-8, atol=1e-6)
    s_end = sol.y[0, -1] / p["N_h"]
    fw.append(dict(s0=s0, m=f["m"], S_over_N_2012=s_end, ratio_to_s0=s_end / s0,
                   Reff_2012=r0(p) * np.sqrt(s_end)))
fw = pd.DataFrame(fw)
print("\n    Forward run of the 2010/11 fit to 29 Apr 2012:")
print(fw.round(3).to_string(index=False))
fw.round(4).to_csv(f"{OUT}/check05_forward_run_2012.csv", index=False)


# (b) shapes
def shape(y, dates, label, note):
    y = np.asarray(y, float)
    dates = np.asarray(dates)
    pk = int(np.argmax(y))
    return dict(series=label, weeks=len(y), total=round(y.sum()), peak=round(y.max(), 1),
                peak_week_start=str(pd.Timestamp(dates[pk]).date()),
                weeks_at_or_above_half_peak=int((y >= y.max() / 2).sum()), note=note)


cal = pd.read_csv(f"{DATA}/02_calibration_san_juan_2010_11.csv", parse_dates=["week_start_date"])
od = pd.read_csv(f"{DATA}/03_validation_puerto_rico_weekly_opendengue.csv", parse_dates=["week_start_date"])
# 2024 epi year = 52 weeks starting 31 Dec 2023 (the first is labelled year 2023 in file 03)
od24 = od[(od.week_start_date >= "2023-12-31") & (od.week_start_date <= "2024-12-22")]
g13 = season_frame("2012/2013")
shapes = pd.DataFrame([
    shape(cal.denv1_attributed_cases, cal.week_start_date, "San Juan 2010/11, DENV-1 attributed",
          "calibration target"),
    shape(cal.total_cases, cal.week_start_date, "San Juan 2010/11, all lab-positive", "all serotypes"),
    shape(g13.denv1_attributed_cases, g13.week_start_date, "San Juan 2012/13, DENV-1 attributed",
          "cross-season holdout"),
    shape(od24.dengue_total, od24.week_start_date, "Puerto Rico 2024, OpenDengue total",
          "island-wide, all serotypes, 52 epi weeks"),
])
print("\n(b) Epidemic shapes:")
print(shapes.to_string(index=False))
shapes.to_csv(f"{OUT}/check05_shapes.csv", index=False)
