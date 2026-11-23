"""Shared helpers for E2's validation checks (role E2, Fahad).

The checks use the project's shared model (dengue.model), parameters (AB's
parameters.json) and data (AB's package). Run them from this folder, or all at once
with ``python run_all.py``. Outputs go to results/validation/.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from dengue import data                                   # noqa: E402
from dengue import observation as obs                     # noqa: E402
from dengue.model import PARAMS, dengue_rhs, r0           # noqa: E402

OUT = os.path.join(ROOT, "results", "validation")
os.makedirs(OUT, exist_ok=True)


def baseline_params():
    return dict(PARAMS)


def load():
    cal = data.calibration()
    return cal, obs.week_edges(cal)


def load_full():
    """All 23 seasons, with the column names the checks were written for."""
    full = data.full_series()
    return full.rename(columns={"other_positive_cases": "untyped_positive_cases",
                                "additional_cases": "estimated_untested_cases"})


def load_opendengue():
    return data.opendengue_weekly()


def fit_m_phi(y, edges, s0, base=None, rhs=dengue_rhs, kind="sqrt", k=None,
              train=None, m0=4.0, phi0=0.01, m_bounds=(0.05, 2000.0),
              phi_bounds=(1e-5, 1.0), n_anchor=3):
    """Fit m and phi with s0 held fixed (log scale for both). Returns dict."""
    base = dict(base or baseline_params())
    y = np.asarray(y, dtype=float)
    idx = np.arange(len(y)) if train is None else np.asarray(train)
    lo = np.log([m_bounds[0], phi_bounds[0]])
    hi = np.log([m_bounds[1], phi_bounds[1]])
    big = 1e3 * np.ones(len(idx))

    def resid(theta):
        p = dict(base, m=float(np.exp(theta[0])))
        try:
            mu = obs.expected_cases(p, s0, float(np.exp(theta[1])), edges[: idx[-1] + 2],
                                    rhs=rhs, y_anchor=y, n_anchor=n_anchor)
        except obs.ObservationError:
            return big
        return obs.residuals(y[idx], mu[idx], kind=kind, k=k)

    x0 = np.clip(np.log([m0, phi0]), lo + 1e-9, hi - 1e-9)
    res = least_squares(resid, x0, bounds=(lo, hi), x_scale="jac", xtol=1e-10, ftol=1e-10)
    m, phi = np.exp(res.x)
    p = dict(base, m=float(m))
    return dict(s0=s0, m=float(m), phi=float(phi), cost=float(res.cost), status=int(res.status),
                nfev=int(res.nfev), R0=r0(p), Reff=r0(p) * np.sqrt(s0), params=p)
