"""Calibration: fit m and phi for each fixed s0.

Role E1 (Zaid): fitting method, from e1_fit.py (Stage 2). It follows E2's plan: one
season pins down only m*s0 and phi*s0, so s0 is fixed at scenario values and m and phi
are fitted on log scales to the square-root residuals of DENV-1-attributed cases, with
E2's expected_cases (real week boundaries, growth-mode start) and fit_report diagnostics.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from . import config
from . import diagnostics as dg
from . import observation as obs
from .model import PARAMS, dengue_rhs, r0

FIT_COLUMNS = ["s0", "m", "phi", "R0", "R_eff", "m_x_s0", "phi_x_s0", "RSS", "cost", "AIC", "AICc",
               "BIC", "cond_number", "corr_m_phi", "se_log_m", "se_log_phi", "lag1_autocorr",
               "runs_z", "model_peak", "data_peak", "model_total", "data_total",
               "peak_week_model", "peak_week_data", "solver_status"]


def residual_function(y, edges, s0, base=None):
    """Square-root residuals as a function of theta = (log m, log phi)."""
    base = dict(base or PARAMS)
    big = 1e3 * np.ones(len(y))

    def resid(theta):
        p = dict(base, m=float(np.exp(theta[0])))
        try:
            mu = obs.expected_cases(p, s0, float(np.exp(theta[1])), edges, rhs=dengue_rhs, y_anchor=y)
        except (obs.ObservationError, ValueError):
            return big
        return obs.residuals(y, mu, kind="sqrt")
    return resid


def fit_s0(y, edges, s0, base=None):
    """Best fit of (m, phi) at a fixed s0. Returns (row dict, expected weekly cases)."""
    base = dict(base or PARAMS)
    y = np.asarray(y, dtype=float)
    lo = np.log([config.FIT_M_BOUNDS[0], config.FIT_PHI_BOUNDS[0]])
    hi = np.log([config.FIT_M_BOUNDS[1], config.FIT_PHI_BOUNDS[1]])
    resid = residual_function(y, edges, s0, base)
    best = None
    for k in config.FIT_STARTS_M_TIMES_S0:            # multi-start along the m*s0 ridge
        x0 = np.clip(np.log([k / s0, 0.002 / s0]), lo + 1e-9, hi - 1e-9)
        r = least_squares(resid, x0, bounds=(lo, hi), x_scale="jac", diff_step=config.FIT_DIFF_STEP,
                          xtol=1e-10, ftol=1e-10)
        if best is None or r.cost < best.cost:
            best = r
    m, phi = (float(v) for v in np.exp(best.x))
    p = dict(base, m=m)
    mu = obs.expected_cases(p, s0, phi, edges, rhs=dengue_rhs, y_anchor=y)
    rep = dg.fit_report(best, param_names=["log_m", "log_phi"], fun=resid)
    R0 = r0(p)
    row = dict(s0=s0, m=m, phi=phi, R0=R0, R_eff=R0 * np.sqrt(s0), m_x_s0=m * s0, phi_x_s0=phi * s0,
               RSS=rep["rss"], cost=float(best.cost), AIC=rep["aic"], AICc=rep["aicc"], BIC=rep["bic"],
               cond_number=rep["condition_number"], corr_m_phi=float(rep["corr"][0, 1]),
               se_log_m=float(rep["se"][0]), se_log_phi=float(rep["se"][1]),
               lag1_autocorr=rep["lag1_autocorrelation"], runs_z=rep["runs_test_z"],
               model_peak=float(mu.max()), data_peak=float(y.max()), model_total=float(mu.sum()),
               data_total=float(y.sum()), peak_week_model=int(mu.argmax()) + 1,
               peak_week_data=int(y.argmax()) + 1, solver_status=int(best.status))
    return row, mu


def calibrate(season: pd.DataFrame, s0_values=config.S0_FIT_GRID, base=None, verbose=True):
    """Fit every s0 in s0_values. Returns (fits table, weekly curves table)."""
    edges = obs.week_edges(season)
    y = season["denv1_attributed_cases"].to_numpy(float)
    rows, curves = [], {"week": np.arange(1, len(y) + 1), "week_start_date": season["week_start_date"].dt.date,
                        "observed": y}
    for s0 in s0_values:
        row, mu = fit_s0(y, edges, s0, base)
        rows.append(row)
        curves[f"fit_s0_{s0:g}"] = mu
        if verbose:
            print(f"  s0 = {s0:<4g}  m = {row['m']:8.3f}  phi = {row['phi']:.5f}  R0 = {row['R0']:.3f}  "
                  f"R_eff = {row['R_eff']:.3f}  RSS = {row['RSS']:.2f}")
    return pd.DataFrame(rows)[FIT_COLUMNS], pd.DataFrame(curves)


def fits_lookup(table: pd.DataFrame) -> dict:
    """{s0: (m, phi)} from a fits table."""
    return {round(float(r.s0), 2): (float(r.m), float(r.phi)) for r in table.itertuples()}
