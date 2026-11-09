#!/usr/bin/env python3
"""Build the AB data and literature package (MATH 574, Team 2) from the files in raw/.

Run from this folder:   python build_ab_package.py
Needs Python 3.10+, numpy, pandas and matplotlib. It rewrites data/, checks/, figures/
and parameters/parameters.json. Files in raw/ are only read.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.dates as mdates  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parent
RAW, DATA, CHECKS, FIG, PAR = (ROOT / p for p in ("raw", "data", "checks", "figures", "parameters"))
SEROTYPES = ["denv1_cases", "denv2_cases", "denv3_cases", "denv4_cases"]
MODEL_KEYS = ["N_h", "m", "a0", "kappa", "eta", "b", "c", "eps_h", "eps_v", "gamma_h", "d_h", "d_v"]
EIP_25C, EIP_30C = 15.0, 6.5  # Chan & Johansson 2012: mean extrinsic incubation period (days)
CALIBRATION, SECONDARY = "2010/2011", "2012/2013"
SOURCES = {
    "cdc": "https://github.com/cdcepi/dengue-forecasting-project-2015",
    "opendengue": "https://github.com/OpenDengue/master-repo (data release V1.3)",
    "paho": "https://opendata.paho.org/en/core-indicators (release of 13 Apr 2026)",
}


def eip_days(temp_c):
    """Mean extrinsic incubation period in days, log-linear in temperature through
    Chan & Johansson's means at 25 C and 30 C."""
    return EIP_25C * (EIP_30C / EIP_25C) ** ((np.asarray(temp_c, dtype=float) - 25.0) / 5.0)


def cdc_week_of_year(dates):
    """Week number used in the CDC San Juan files: 7-day steps from 1 January, with
    week 52 taking the extra 1-2 days at the end of December."""
    doy = pd.DatetimeIndex(dates).dayofyear.to_numpy()
    return np.minimum(52, (doy - 1) // 7 + 1)


def mmwr_week(start):
    """MMWR epidemiological year and week of a Sunday-start week."""
    mid = start + pd.Timedelta(days=3)
    year = mid.year
    jan4 = pd.Timestamp(year=year, month=1, day=4)
    week1 = jan4 - pd.Timedelta(days=(jan4.dayofweek + 1) % 7)
    return year, (start - week1).days // 7 + 1


# ---------------------------------------------------------------- cases
def load_cases() -> pd.DataFrame:
    c = pd.read_csv(RAW / "cdc" / "san_juan_testing_data.csv", parse_dates=["week_start_date"])
    c = c.sort_values("week_start_date").reset_index(drop=True)
    nxt = c.week_start_date.shift(-1)
    c["week_end_date"] = (nxt - pd.Timedelta(days=1)).fillna(c.week_start_date + pd.Timedelta(days=6))
    c["days_in_week"] = (c.week_end_date - c.week_start_date).dt.days + 1
    c["weekofyear"] = cdc_week_of_year(c.week_start_date)
    c["days_since_start"] = (c.week_start_date - c.groupby("season").week_start_date.transform("min")).dt.days
    typed = c[SEROTYPES].sum(axis=1)
    share = np.where(typed > 0, c.denv1_cases / typed.where(typed > 0, 1), np.nan)
    c["typed_cases"] = typed
    c["denv1_share_of_typed"] = np.round(share, 4)
    untyped = c.other_positive_cases + c.additional_cases
    c["denv1_attributed_cases"] = np.round(c.denv1_cases + untyped * np.nan_to_num(share), 3)
    return c


# ---------------------------------------------------------------- temperature
def load_daily_weather() -> pd.DataFrame:
    w = pd.read_csv(RAW / "cdc" / "SanJuanRQW00011641.csv", skipinitialspace=True)
    w.columns = [col.strip() for col in w.columns]
    w["date"] = pd.to_datetime(dict(year=w.YYYY, month=w.MM, day=w.DD))
    for col in ["TMAX", "TMIN", "TAVG", "TDTR", "PRCP"]:
        w[col] = pd.to_numeric(w[col], errors="coerce").where(lambda s: s > -999)
    return w


def weekly_weather(cases: pd.DataFrame, daily: pd.DataFrame) -> pd.DataFrame:
    span = (daily.date >= cases.week_start_date.min()) & (daily.date <= cases.week_end_date.max())
    d = daily[span].copy()
    d["row"] = np.searchsorted(cases.week_start_date.to_numpy(), d.date.to_numpy(), side="right") - 1
    g = d.groupby("row")
    out = pd.DataFrame({
        "temp_mean_c": g.TAVG.mean().round(2),
        "temp_min_mean_c": g.TMIN.mean().round(2),
        "temp_max_mean_c": g.TMAX.mean().round(2),
        "precip_total_mm": g.PRCP.sum(min_count=1).round(1),
        "days_with_temperature": g.TAVG.count(),
    }).reindex(range(len(cases)))
    out["eip_days"] = np.round(eip_days(out.temp_mean_c), 2)
    out["eps_v_per_day"] = np.round(1.0 / out.eip_days, 5)
    return out


# ---------------------------------------------------------------- population
def demography():
    pop = pd.read_csv(RAW / "cdc" / "san_juan_population_data.csv").rename(
        columns={"Year": "year", "Estimated_population": "san_juan_msa_population"})
    p = pd.read_csv(RAW / "paho" / "PAHO_Core_Indicators_2026_PRI_subset.csv")
    pick = {275: "puerto_rico_population_thousands", 74: "life_expectancy_years",
            414: "crude_birth_rate_per_1000", 415: "crude_death_rate_per_1000"}
    wide = (p[p.paho_indicator_id.isin(pick)]
            .pivot_table(index="time_dim", columns="paho_indicator_id", values="numeric_value", aggfunc="first")
            .rename(columns=pick).reset_index().rename(columns={"time_dim": "year"}))
    demo = wide.merge(pop, on="year", how="left").sort_values("year")
    demo = demo[(demo.year >= 1990) & (demo.year <= 2024)].reset_index(drop=True)
    demo["d_h_per_day"] = (1.0 / (demo.life_expectancy_years * 365.25)).round(9)
    demo["san_juan_msa_population"] = demo.san_juan_msa_population.astype("Int64")
    cols = ["year", "san_juan_msa_population", "puerto_rico_population_thousands", "life_expectancy_years",
            "crude_birth_rate_per_1000", "crude_death_rate_per_1000", "d_h_per_day"]
    # PAHO's "0-1" band covers ages 0 and 1, so it overlaps "1-4" by one birth cohort. The clean
    # table uses 0-4 = total population minus the bands from age 5 up.
    y10 = p[p.time_dim == 2010].drop_duplicates("paho_indicator_id").set_index("paho_indicator_id").numeric_value
    older = {358: "5-14", 359: "15-24", 360: "25-44", 361: "45-64", 362: "65-74", 363: "75+"}
    total = y10.loc[275]
    rows = [("0-4", total - sum(y10.loc[k] for k in older))] + [(v, y10.loc[k]) for k, v in older.items()]
    age = pd.DataFrame(rows, columns=["age_group", "population_thousands"])
    age["population_thousands"] = age.population_thousands.round(3)
    age["share"] = (age.population_thousands / total).round(4)
    age["cumulative_share"] = (age.population_thousands.cumsum() / total).round(4)
    overlap = (y10.loc[356] + y10.loc[357] - (total - sum(y10.loc[k] for k in older)), y10.loc[414] / 1000 * total)
    return demo[cols], age, overlap


# ---------------------------------------------------------------- OpenDengue
def puerto_rico_weekly():
    od = pd.read_csv(RAW / "opendengue" / "Temporal_extract_PAHO_V1_3_PUERTO_RICO.csv",
                     parse_dates=["calendar_start_date", "calendar_end_date"])
    wk = od[od.T_res == "Week"].sort_values("calendar_start_date").reset_index(drop=True)
    yw = [mmwr_week(s) for s in wk.calendar_start_date]
    out = pd.DataFrame({
        "week_start_date": wk.calendar_start_date.dt.date,
        "week_end_date": wk.calendar_end_date.dt.date,
        "epi_year": [y for y, _ in yw],
        "epi_week": [w for _, w in yw],
        "year_label_in_source": wk.Year,
        "dengue_total": wk.dengue_total,
        "case_definition": wk.case_definition_standardised,
        "spatial_resolution": wk.S_res,
        "source_id": wk.UUID,
    })
    # Years where every week is zero are reporting gaps, not years without dengue.
    out["dengue_total"] = out.dengue_total.round().astype("Int64")
    zeros = out.groupby("epi_year").dengue_total.transform(lambda s: int((s == 0).sum()))
    weeks = out.groupby("epi_year").dengue_total.transform("size")
    out.insert(6, "reporting_note", "")
    all_zero = zeros == weeks
    out.loc[all_zero, "reporting_note"] = "every week zero in source; treated as missing"
    out.loc[all_zero, "dengue_total"] = pd.NA
    out.loc[~all_zero & (zeros >= 13), "reporting_note"] = "13+ zero weeks this year; reporting likely incomplete"
    # OpenDengue marks values it filled in with a suffix on the record ID, e.g. "(Imputed)".
    out.insert(7, "opendengue_flag", wk.UUID.str.extract(r"\(([^)]*)\)\s*$")[0].str.lower().fillna("").to_numpy())
    nat = pd.read_csv(RAW / "opendengue" / "National_extract_V1_3_PUERTO_RICO.csv",
                      parse_dates=["calendar_start_date"])
    annual = nat[nat.T_res == "Year"].assign(year=lambda x: x.calendar_start_date.dt.year)
    return out, annual.set_index("year").dengue_total


# ---------------------------------------------------------------- figure
def case_plot(c: pd.DataFrame, path: Path):
    ink, ink2, grid, grey = "#0b0b0b", "#52514e", "#e6e5e0", "#a3a29c"
    blue, orange = "#2a78d6", "#eb6834"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.edgecolor": grid,
                         "axes.labelcolor": ink2, "xtick.color": ink2, "ytick.color": ink2})
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7.6), gridspec_kw={"height_ratios": [1, 1.05], "hspace": 0.55})
    fig.patch.set_facecolor("white")
    for ax in (ax1, ax2):
        ax.set_facecolor("white")
        ax.grid(axis="y", color=grid, linewidth=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(length=0)

    gap = c[(c.week_start_date >= "2000-04-30") & (c.week_start_date < "2007-04-29")]
    ax1.axvspan(gap.week_start_date.min(), gap.week_end_date.max(), color="#f1f0ec", zorder=0)
    ax1.text(gap.week_start_date.min() + pd.Timedelta(days=60), 395, "DENV-1 absent\n2000/01-2006/07",
             color=ink2, fontsize=8, va="top")
    ax1.plot(c.week_start_date, c.total_cases, color=grey, linewidth=1.0, label="All lab-positive cases (total_cases)")
    ax1.plot(c.week_start_date, c.denv1_attributed_cases, color=blue, linewidth=1.2,
             label="DENV-1 attributed (denv1_attributed_cases)")
    ax1.set_ylim(0, 470)
    ax1.set_ylabel("Cases per week")
    ax1.xaxis.set_major_locator(mdates.YearLocator(2))
    ax1.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax1.set_xlim(c.week_start_date.min(), c.week_end_date.max())
    ax1.legend(loc="upper right", frameon=False, fontsize=8, labelcolor=ink2)
    ax1.set_title("San Juan, 1990-2013: DENV-1 vanished for seven seasons, then drove the 2010/11 and 2012/13 seasons",
                  loc="left", fontsize=10.5, color=ink, fontweight="bold", pad=10)

    s = c[c.season == CALIBRATION]
    ax2.plot(s.week_start_date, s.total_cases, color=grey, linewidth=1.6, label="All lab-positive cases")
    ax2.plot(s.week_start_date, s.denv1_cases, color=orange, linewidth=1.6, linestyle=(0, (4, 2)),
             label="Typed as DENV-1 by RT-PCR")
    ax2.plot(s.week_start_date, s.denv1_attributed_cases, color=blue, linewidth=2.0,
             label="DENV-1 attributed: the fit target")
    pk = s.loc[s.denv1_attributed_cases.idxmax()]
    ax2.plot([pk.week_start_date], [pk.denv1_attributed_cases], "o", color=blue, markersize=6,
             markeredgecolor="white", markeredgewidth=1.5)
    ax2.annotate(f"Peak {pk.denv1_attributed_cases:.0f} (week of {pk.week_start_date:%d %b %Y})",
                 (pk.week_start_date, pk.denv1_attributed_cases), xytext=(14, 4), textcoords="offset points",
                 fontsize=8, color=ink)
    long = s[s.days_in_week > 7].iloc[0]
    ax2.annotate(f"{long.days_in_week}-day week ({long.week_start_date:%d}-{long.week_end_date:%d %b})",
                 (long.week_start_date, long.denv1_attributed_cases), xytext=(0, 38), textcoords="offset points",
                 fontsize=8, color=ink2, ha="center", arrowprops=dict(arrowstyle="-", color=ink2, linewidth=0.7))
    ax2.set_ylim(0, 300)
    ax2.set_ylabel("Cases per week")
    ax2.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[5, 7, 9, 11, 1, 3]))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax2.set_xlim(s.week_start_date.min(), s.week_end_date.max())
    ax2.legend(loc="upper right", frameon=False, fontsize=8, labelcolor=ink2)
    ax2.set_title(f"2010/11 season (30 Apr 2010 - 29 Apr 2011): {s.denv1_attributed_cases.sum():,.0f} DENV-1-attributed "
                  f"of {s.total_cases.sum():,} lab-positive cases", loc="left", fontsize=10.5, color=ink,
                  fontweight="bold", pad=10)
    fig.text(0.125, 0.015, "Source: CDC and Puerto Rico Department of Health weekly data, San Juan-Carolina-Caguas "
             "metro area (cdcepi/dengue-forecasting-project-2015). Weeks by symptom onset.", fontsize=7.5, color=ink2)
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- main
def main():
    for d in (DATA, CHECKS, FIG):
        d.mkdir(exist_ok=True)
    cases = load_cases()
    daily = load_daily_weather()
    wx = weekly_weather(cases, daily)
    full = pd.concat([cases, wx], axis=1)
    demo, age, overlap = demography()
    prw, pr_annual = puerto_rico_weekly()

    base = ["season", "season_week", "week_start_date", "week_end_date", "days_in_week", "days_since_start",
            "weekofyear", "total_cases", *SEROTYPES, "other_positive_cases", "additional_cases", "typed_cases",
            "denv1_share_of_typed", "denv1_attributed_cases"]
    climate = ["temp_mean_c", "eip_days", "eps_v_per_day"]
    date_fmt = {"week_start_date": lambda d: d.week_start_date.dt.date,
                "week_end_date": lambda d: d.week_end_date.dt.date}
    full[base].assign(**date_fmt).to_csv(DATA / "san_juan_weekly_cases_1990_2013.csv", index=False)
    for season, name in ((CALIBRATION, "calibration_san_juan_2010_11.csv"), (SECONDARY, "secondary_san_juan_2012_13.csv")):
        full.loc[full.season == season, base + climate].assign(**date_fmt).to_csv(DATA / name, index=False)
    wcols = ["season", "season_week", "week_start_date", "week_end_date", "days_in_week", "temp_mean_c",
             "temp_min_mean_c", "temp_max_mean_c", "precip_total_mm", "days_with_temperature", "eip_days", "eps_v_per_day"]
    full[wcols].assign(**date_fmt).to_csv(DATA / "temperature_san_juan_weekly_1990_2013.csv", index=False)
    demo.to_csv(DATA / "demography_puerto_rico_1990_2024.csv", index=False)
    age.to_csv(DATA / "age_structure_puerto_rico_2010.csv", index=False)
    prw.to_csv(DATA / "puerto_rico_weekly_opendengue.csv", index=False)
    prw[prw.epi_year == 2024].to_csv(DATA / "validation_puerto_rico_2024.csv", index=False)

    # parameters.json for D's dengue_model.py
    par = pd.read_csv(PAR / "parameters.csv").set_index("key")
    params = {k: float(par.loc[k, "baseline"]) for k in MODEL_KEYS}
    params["N_h"] = int(params["N_h"])
    rng = {k: [None if pd.isna(par.loc[k, "low"]) else float(par.loc[k, "low"]),
               None if pd.isna(par.loc[k, "high"]) else float(par.loc[k, "high"])] for k in par.index}
    (PAR / "parameters.json").write_text(json.dumps({
        "PARAMS": params,
        "fitted_or_scenario": {k: {"baseline": None if pd.isna(par.loc[k, "baseline"]) else float(par.loc[k, "baseline"]),
                                   "range": rng[k], "status": par.loc[k, "status"]}
                               for k in ["m", "phi", "s0", "kappa", "eta"]},
        "ranges": rng,
        "units": "time in days, rates per day, populations in persons, fractions for phi, s0, kappa and eta",
        "source": "parameters/parameters.csv",
    }, indent=2, ensure_ascii=False) + "\n")

    case_plot(cases, FIG / "cases_san_juan.png")
    write_checks(cases, wx, daily, demo, age, overlap, prw, pr_annual, par)
    write_manifest()
    print("built:", ", ".join(sorted(p.name for p in DATA.glob("*.csv"))))


def write_checks(c, wx, daily, demo, age, overlap, prw, pr_annual, par):
    L = ["# Data checks", "", "Written by `build_ab_package.py`; every number below is computed from `raw/`.", ""]

    def line(ok, text):
        L.append(f"- **{'PASS' if ok else 'NOTE'}**: {text}")

    L += ["## San Juan weekly cases (CDC/PRDH, 1990/91-2012/13)", ""]
    n_seasons, n_rows = c.season.nunique(), len(c)
    line(n_seasons == 23 and n_rows == 1196 and c.groupby("season").size().eq(52).all(),
         f"{n_seasons} seasons x 52 weeks = {n_rows} rows; first week starts {c.week_start_date.min():%d %b %Y}, "
         f"last week ends {c.week_end_date.max():%d %b %Y}.")
    line(not c.week_start_date.duplicated().any(), "No duplicate weeks.")
    cols = ["total_cases", *SEROTYPES, "other_positive_cases", "additional_cases"]
    line(not c[cols].isna().any().any() and (c[cols] >= 0).all().all(), "No missing or negative counts.")
    parts = c[SEROTYPES + ["other_positive_cases", "additional_cases"]].sum(axis=1)
    line((parts == c.total_cases).all(), "total_cases equals DENV-1 to -4 + untyped + estimated untested in every week.")
    lens = c.days_in_week.value_counts().sort_index()
    long = c[c.days_in_week > 7]
    line(set(lens.index) <= {7, 8, 9} and (long.week_end_date.dt.strftime("%m-%d") == "12-31").all(),
         "Week lengths: " + ", ".join(f"{n} weeks of {d} days" for d, n in lens.items())
         + ". Every long week ends on 31 Dec (weeks run from 1 Jan in 7-day steps), so fit on real dates, not week numbers.")
    blind = (c.typed_cases == 0) & (c.other_positive_cases + c.additional_cases > 0)
    late = c.week_start_date >= "2008-04-27"
    in_fit = blind & c.season.isin([CALIBRATION, SECONDARY])
    line(not in_fit.any(),
         f"Weeks with untyped cases but no typed case, where the DENV-1 share is unknown and attribution gives 0: "
         f"{int(blind.sum())} of {len(c)} in the full series, {int(blind[late].sum())} from 2008/09 on, "
         f"{int(in_fit.sum())} in 2010/11 and 2012/13. Before 2008 few positives were typed, so "
         f"`denv1_attributed_cases` is only meaningful from 2008/09.")
    line(True, "Geography: San Juan-Carolina-Caguas Metropolitan Statistical Area (raw/cdc/san_juan_metadata.md).")
    line(True, "Case definition: laboratory-positive (serotype by RT-PCR, or IgM-positive), plus `additional_cases`, "
               "which CDC estimated for weeks when not all specimens were tested. Dated by week of symptom onset.")

    for season in (CALIBRATION, SECONDARY):
        s = c[c.season == season]
        typed = s[SEROTYPES].sum()
        big = s[s.total_cases >= 20]
        lab = big.typed_cases + big.other_positive_cases
        pk = s.loc[s.denv1_attributed_cases.idxmax()]
        no_type = int(((s.typed_cases == 0) & (s.other_positive_cases + s.additional_cases > 0)).sum())
        L += ["", f"## Season {season}", ""]
        line(True, f"{s.week_start_date.min():%d %b %Y} to {s.week_end_date.max():%d %b %Y}, {len(s)} weeks, "
                   f"{s.days_in_week.sum()} days.")
        line(True, f"Totals: all lab-positive {s.total_cases.sum():,}; typed DENV-1 {s.denv1_cases.sum():,}; "
                   f"DENV-1 attributed {s.denv1_attributed_cases.sum():,.1f}; estimated untested {s.additional_cases.sum()}.")
        line(True, "Serotype mix of typed cases: " + ", ".join(
            f"DENV-{i + 1} {typed.iloc[i] / typed.sum():.1%}" for i in range(4)) + ".")
        line(True, f"Share of lab-positive cases with a serotype, weeks with 20+ cases: "
                   f"{(big.typed_cases / lab).min():.2f}-{(big.typed_cases / lab).max():.2f}.")
        line(no_type == 0, f"Weeks with untyped cases but no typed case (DENV-1 share unknown, attributed as 0): {no_type}.")
        line(True, f"DENV-1-attributed peak: {pk.denv1_attributed_cases:.1f} in the week of {pk.week_start_date:%d %b %Y} "
                   f"(season week {pk.season_week}); weeks at half the peak or more: "
                   f"{int((s.denv1_attributed_cases >= s.denv1_attributed_cases.max() / 2).sum())}.")

    hist = c.groupby("season")[SEROTYPES].sum()
    absent = hist.index[hist.denv1_cases <= 1].tolist()
    before = c[c.season.isin(["2007/2008", "2008/2009", "2009/2010"])]
    L += ["", "## DENV-1 history (used for the s0 evidence)", ""]
    line(True, f"Seasons with at most 1 typed DENV-1 case: {', '.join(absent)}.")
    line(True, "Typed DENV-1 cases before the calibration season: " + ", ".join(
        f"{k} {int(v)}" for k, v in before.groupby("season").denv1_cases.sum().items())
        + f"; DENV-1 attributed over those three seasons: {before.denv1_attributed_cases.sum():,.0f}.")

    L += ["", "## Population", ""]
    n10 = int(demo.loc[demo.year == 2010, "san_juan_msa_population"].iloc[0])
    pr10 = float(demo.loc[demo.year == 2010, "puerto_rico_population_thousands"].iloc[0])
    le10 = float(demo.loc[demo.year == 2010, "life_expectancy_years"].iloc[0])
    line(n10 == int(par.loc["N_h", "baseline"]), f"N_h = {n10:,} (San Juan metro, 2010, US Census estimate in the CDC repo) "
                                                 f"= {n10 / (pr10 * 1000):.0%} of Puerto Rico's {pr10 / 1000:.2f} million (PAHO).")
    dh = 1 / (le10 * 365.25)
    line(abs(dh - float(par.loc["d_h", "baseline"])) / dh < 0.01,
         f"d_h = 1 / ({le10:.3f} years x 365.25) = {dh:.4e} per day (PAHO life expectancy, Puerto Rico 2010).")
    u15 = age.loc[age.age_group.isin(["0-4", "5-14"]), "share"].sum()
    u25 = u15 + age.loc[age.age_group == "15-24", "share"].sum()
    ov, births = overlap
    line(abs(ov - births) / births < 0.1,
         f"PAHO's 0-1 and 1-4 age bands overlap at age 1: together they exceed the total by {ov:.1f} thousand, about one "
         f"year of births ({births:.1f} thousand). The clean table uses 0-4 = total minus the bands from age 5 up.")
    line(True, f"Puerto Rico 2010 age structure: {u15:.1%} under 15, {u25:.1%} under 25.")

    p = {k: float(par.loc[k, "baseline"]) for k in MODEL_KEYS}

    def r0(m):
        a = p["a0"] * (1 - p["kappa"] * p["eta"])
        return np.sqrt(a ** 2 * p["b"] * p["c"] * m * p["eps_h"] * p["eps_v"]
                       / ((p["d_h"] + p["eps_h"]) * (p["d_h"] + p["gamma_h"]) * (p["d_v"] + p["eps_v"]) * p["d_v"]))
    L += ["", "## Parameter table", ""]
    line(True, f"R0 from the closed form with the baseline values: {r0(p['m']):.2f} at m = {p['m']:g}, and "
               f"{r0(3.9 / 0.4):.2f} at m = 9.75 (the value implied at s0 = 0.4).")

    L += ["", "## Temperature (NOAA GHCN-Daily, station RQW00011641, San Juan airport)", ""]
    span = daily[(daily.date >= c.week_start_date.min()) & (daily.date <= c.week_end_date.max())]
    days = (c.week_end_date.max() - c.week_start_date.min()).days + 1
    line(len(span) == days and span.TAVG.notna().all(),
         f"{len(span):,} of {days:,} days present; days without a daily mean: {int(span.TAVG.isna().sum())}.")
    tmean = span.TAVG.mean()
    eip = float(eip_days(tmean))
    line(abs(1 / eip - float(par.loc["eps_v", "baseline"])) / (1 / eip) < 0.02,
         f"Mean daily temperature 30 Apr 1990-29 Apr 2013: {tmean:.2f} C, so mean extrinsic incubation "
         f"{eip:.2f} days and eps_v = {1 / eip:.4f} per day (parameter table: {float(par.loc['eps_v', 'baseline']):.4f}).")
    s = pd.concat([c, wx], axis=1)
    s = s[s.season == CALIBRATION]
    line(True, f"2010/11 weekly mean temperature {s.temp_mean_c.min():.1f}-{s.temp_mean_c.max():.1f} C, so extrinsic incubation "
               f"{s.eip_days.min():.1f}-{s.eip_days.max():.1f} days.")

    L += ["", "## Puerto Rico weekly cases, OpenDengue V1.3 (island-wide)", ""]
    p24 = prw[prw.epi_year == 2024]
    line(len(p24) == 52, f"Epi year 2024: {len(p24)} weeks, {p24.week_start_date.min()} to {p24.week_end_date.max()}, "
                         f"{int(p24.dengue_total.sum()):,} cases.")
    lab23 = int((p24.year_label_in_source != 2024).sum())
    line(lab23 == 0, f"Weeks of epi year 2024 labelled with another year in the source: {lab23} (filter by date or "
                     f"`epi_year`, not by the source's year label).")
    if 2024 in pr_annual.index:
        line(int(pr_annual.loc[2024]) == int(p24.dengue_total.sum()),
             f"Weekly 2024 total matches OpenDengue's annual 2024 figure ({int(pr_annual.loc[2024]):,}).")
    pk = p24.loc[p24.dengue_total.idxmax()]
    line(True, f"2024 peak {int(pk.dengue_total)} in the week of {pk.week_start_date}; weeks at half the peak or more: "
               f"{int((p24.dengue_total >= p24.dengue_total.max() / 2).sum())}.")
    gaps = int((pd.to_datetime(prw.week_start_date).diff().dt.days.dropna() != 7).sum())
    w53 = prw.groupby("epi_year").epi_week.max()
    line(gaps == 0, f"Weekly coverage {prw.week_start_date.min()} to {prw.week_end_date.max()}, {gaps} gaps; years with 53 "
                    f"epi weeks: {', '.join(map(str, w53.index[w53 == 53]))}.")
    miss = sorted(prw.loc[prw.reporting_note.str.startswith("every"), "epi_year"].unique())
    sparse = (prw[prw.reporting_note.str.startswith("13+")].groupby("epi_year").dengue_total
              .apply(lambda s: int((s == 0).sum())))
    line(False, f"Years where every week is zero in the source, set to missing: {', '.join(map(str, miss)) or 'none'}. "
                "Years flagged as likely incomplete: "
                + (", ".join(f"{y} ({n} zero weeks)" for y, n in sparse.items()) or "none") + ".")
    flagged = prw[prw.opendengue_flag != ""]
    by_year = flagged.groupby(["epi_year", "opendengue_flag"]).size()
    line(not (p24.opendengue_flag != "").any(),
         f"Weeks OpenDengue filled in itself (marked in the record ID): {len(flagged)}, by year and type: "
         + ", ".join(f"{y} {t} {n}" for (y, t), n in by_year.items())
         + f". Epi year 2024 has {int((p24.opendengue_flag != '').sum())}.")
    line(True, "Case definition 'Total' (suspected and confirmed as reported to PAHO), island-wide, no serotype, no "
               "municipality.")
    (CHECKS / "data_checks.md").write_text("\n".join(L) + "\n")


def write_manifest():
    rows = []
    for f in sorted(RAW.rglob("*")):
        if f.is_file():
            src = SOURCES[f.parent.name] if f.parent.name in SOURCES else ""
            rows.append({"file": str(f.relative_to(ROOT)), "bytes": f.stat().st_size,
                         "sha256": hashlib.sha256(f.read_bytes()).hexdigest(), "source": src})
    pd.DataFrame(rows).to_csv(CHECKS / "raw_manifest.csv", index=False)


if __name__ == "__main__":
    main()
