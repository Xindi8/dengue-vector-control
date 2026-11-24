"""Intervention scenarios, feasibility thresholds and the equal-cut comparison.

Role F (Deckar): scenario design, packages and analysis, from f_model.py, f_scenarios.py
and f_contour_plots.py. Each control changes one thing in the shared model:

    personal protection   a  -> a0 (1 - kappa eta)
    source reduction      m  -> m (1 - delta), and the starting mosquitoes shrink with it
    Wolbachia             b, c -> each times sqrt(1 - p_w q), so b c -> b c (1 - p_w q)

so the reproduction number scales as R (1 - kappa eta) sqrt(1 - delta) sqrt(1 - p_w q).
Simulations start from E2's anchored state on 30 Apr 2010, with E1's fitted m and phi,
and the controls in place all season.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config
from . import observation as obs
from .model import dengue_rhs, r0


def controlled_r(r_base, kappa=0.0, delta=0.0, p_w=0.0, eta=config.ETA, q=config.Q_WOLBACHIA):
    """R under controls, from the closed form (C's formula as used by F)."""
    return r_base * (1 - kappa * eta) * np.sqrt(1 - delta) * np.sqrt(1 - p_w * q)


def apply_controls(p, y0, kappa=0.0, delta=0.0, p_w=0.0, eta=config.ETA, q=config.Q_WOLBACHIA):
    pc = dict(p, kappa=kappa, eta=eta, m=p["m"] * (1 - delta))
    shrink = np.sqrt(1 - p_w * q)
    pc["b"], pc["c"] = p["b"] * shrink, p["c"] * shrink
    yc = np.array(y0, dtype=float)
    yc[4:7] *= (1 - delta)
    return pc, yc


class Season:
    """Observed season plus E1's fits: everything a scenario run needs."""

    def __init__(self, season_df: pd.DataFrame, fits: dict, base: dict):
        self.y = season_df["denv1_attributed_cases"].to_numpy(float)
        self.edges = obs.week_edges(season_df)
        self.edges_2y = np.concatenate([self.edges, self.edges[-1] + 7 * np.arange(1, 53)])
        self.fits = fits
        self.base = dict(base)

    def start(self, s0):
        m, phi = self.fits[round(s0, 2)]
        p = dict(self.base, m=m)
        rate = obs.anchored_onset_rate(dengue_rhs, p, s0, phi, self.edges, self.y)
        return p, obs.initial_state(p, s0, rate), phi

    def run(self, s0, kappa=0.0, delta=0.0, p_w=0.0, eta=config.ETA, q=config.Q_WOLBACHIA, two_years=False):
        p, y0, phi = self.start(s0)
        pc, yc = apply_controls(p, y0, kappa, delta, p_w, eta, q)
        edges = self.edges_2y if two_years else self.edges
        return pc, phi * obs.simulate_weekly_onsets(dengue_rhs, pc, yc, edges)

    def r0_at(self, s0):
        m, _ = self.fits[round(s0, 2)]
        return r0(dict(self.base, m=m))


def balanced_delta(r0_head, kappa=config.BALANCED_KAPPA, eta=config.ETA, step=config.BALANCED_STEP):
    """F's balanced package: keep kappa at its realistic value and raise delta in steps
    until R0 < 1 at the headline s0 (None if no delta up to 0.95 is enough)."""
    for i in range(1, int(round(0.95 / step)) + 1):
        delta = round(i * step, 4)
        if controlled_r(r0_head, kappa, delta, 0.0, eta) < 1:
            return delta
    return None


def packages(season: Season, headline=config.S0_HEADLINE):
    out = []
    for sid, name, k, d, w in config.PACKAGES:
        if d is None:
            d = balanced_delta(season.r0_at(headline))
        out.append((sid, name, k, d, w))
    return out


def scenario_table(season: Season, s0_levels, headline=config.S0_HEADLINE):
    """F's package table: R0 at each s0, and at the headline s0 R_eff, 52-week cases,
    share averted, peak, and cases in the following 52 weeks."""
    levels = sorted(set(s0_levels) | {headline})
    base = {s0: season.run(s0)[1] for s0 in levels}
    rows, curves = [], {}
    for sid, name, k, d, w in packages(season, headline):
        row = dict(ID=sid, Name=name, kappa=k, delta=d, p_w=w)
        for s0 in levels:
            pc, mu = season.run(s0, k, d, w)
            row[f"R0_s0={s0:g}"] = r0(pc)
            row[f"averted_s0={s0:g}"] = 100 * (1 - mu.sum() / base[s0].sum())
            if s0 == headline:
                row["R_eff"] = r0(pc) * np.sqrt(s0)
                row["cases_52wk"] = mu.sum()
                row["pct_averted"] = row[f"averted_s0={s0:g}"]
                row["peak"] = mu.max()
                row["peak_week"] = int(mu.argmax()) + 1
                curves[sid] = mu
                row["cases_wk53_104"] = season.run(s0, k, d, w, two_years=True)[1][52:].sum()
        row["feasible"] = k <= config.KAPPA_MAX and d <= config.DELTA_MAX
        rows.append(row)
    table = pd.DataFrame(rows)
    lead = ["ID", "Name", "kappa", "delta", "p_w", "R_eff", "cases_52wk", "pct_averted", "peak", "peak_week",
            "cases_wk53_104"]
    rest = [c for c in table.columns if c not in lead and c != "feasible"]
    curves = pd.DataFrame(curves, index=np.arange(1, len(season.y) + 1)).rename_axis("week")
    return table[lead + rest + ["feasible"]], curves


def rq2_table(season: Season, headline=config.S0_HEADLINE, cuts=config.RQ2_CUTS):
    """Same percentage cut to the biting rate (kappa = x at eta = 1) or to mosquito density."""
    _, base = season.run(headline)
    rows = []
    for x in cuts:
        for lever, (k, d) in (("biting rate a", (x, 0.0)), ("mosquito density m", (0.0, x))):
            pc, mu = season.run(headline, k, d, 0.0, eta=1.0)
            rows.append(dict(cut=x, lever=lever, R_eff=r0(pc) * np.sqrt(headline),
                             pct_averted=100 * (1 - mu.sum() / base.sum()),
                             peak_drop_pct=100 * (1 - mu.max() / base.max()),
                             peak_week=int(mu.argmax()) + 1))
    return pd.DataFrame(rows)


def single_lever_thresholds(r_eff, r0_by_s0: dict, eta=config.ETA, q=config.Q_WOLBACHIA):
    """Smallest level of each lever alone that brings R below 1: kappa = (1 - 1/R)/eta,
    delta = 1 - 1/R^2, p_w = (1 - 1/R^2)/q. Values above 1 are impossible."""
    targets = {"R_eff (stop 2010/11)": r_eff}
    targets.update({f"R0, s0 = {s0:g}": r for s0, r in sorted(r0_by_s0.items())})
    rows = []
    for lever, limit in (("Personal protection, kappa", config.KAPPA_MAX),
                         ("Source reduction, delta", config.DELTA_MAX),
                         ("Wolbachia coverage, p_w", 1.0)):
        row = dict(lever=lever, ceiling=limit)
        for label, R in targets.items():
            need = 1 - 1 / R if lever.startswith("Personal") else 1 - 1 / R ** 2
            row[label] = need / eta if lever.startswith("Personal") else (need / q if "Wolbachia" in lever else need)
        rows.append(row)
    return pd.DataFrame(rows)
