"""Tests for E2's residual diagnosis module. Run from the e2 folder:  python -m pytest -q"""
import os
import sys

import numpy as np
from scipy.optimize import least_squares

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import diagnostics as dg                                  # noqa: E402
import observation as obs                                 # noqa: E402
from model_provisional import baseline_params, dengue_rhs  # noqa: E402


def _linear_fit(X, y):
    return least_squares(lambda b: X @ b - y, np.zeros(X.shape[1]))


def test_matches_ordinary_least_squares():
    rng = np.random.default_rng(3)
    X = np.column_stack([np.ones(60), rng.normal(size=60), rng.normal(size=60)])
    y = X @ np.array([1.0, 2.0, -0.5]) + rng.normal(scale=0.3, size=60)
    rep = dg.fit_report(_linear_fit(X, y), param_names=["b0", "b1", "b2"])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    rss = float(np.sum((y - X @ beta) ** 2))
    cov = rss / (60 - 3) * np.linalg.inv(X.T @ X)
    assert np.isclose(rep["rss"], rss, rtol=1e-8)
    assert np.allclose(rep["cov"], cov, rtol=1e-6)
    assert np.isclose(rep["aic"], 60 * np.log(rss / 60) + 6, rtol=1e-10)
    assert np.allclose(np.diag(rep["corr"]), 1) and np.allclose(rep["corr"], rep["corr"].T)
    assert rep["flagged_pairs"] == []


def test_flags_nearly_collinear_parameters():
    rng = np.random.default_rng(4)
    x = rng.normal(size=80)
    X = np.column_stack([x, x + 1e-4 * rng.normal(size=80)])     # two copies of one regressor
    y = 2 * x + rng.normal(scale=0.1, size=80)
    rep = dg.fit_report(_linear_fit(X, y), param_names=["u", "v"])
    assert rep["flagged_pairs"] and abs(rep["flagged_pairs"][0][2]) > 0.95
    assert rep["condition_number"] > 1e3
    assert "NOT SEPARATELY IDENTIFIABLE" in dg.format_report(rep)


def test_flags_the_s0_m_phi_ridge():
    """With s0, m and phi all free, the module must flag them (the identifiability ridge)."""
    cal = obs.load_calibration(os.path.join(ROOT, "data", "02_calibration_san_juan_2010_11.csv"))
    edges = obs.week_edges(cal)
    y_anchor = cal.denv1_attributed_cases.to_numpy(float)
    base = baseline_params()
    truth = obs.expected_cases(dict(base, m=10.0), 0.4, 0.005, edges, rhs=dengue_rhs, y_anchor=y_anchor)

    def resid(th):
        s0, m, phi = np.exp(th)
        mu = obs.expected_cases(dict(base, m=m), s0, phi, edges, rhs=dengue_rhs, y_anchor=y_anchor)
        return obs.residuals(truth, mu, kind="sqrt")

    res = least_squares(resid, np.log([0.35, 11.0, 0.0055]), diff_step=1e-4,
                        bounds=(np.log([0.01, 0.1, 1e-4]), np.log([1.0, 500.0, 0.5])))
    rep = dg.fit_report(res, param_names=["log_s0", "log_m", "log_phi"], fun=resid)
    pairs = {(a, b) for a, b, _ in rep["flagged_pairs"]}
    assert ("log_s0", "log_m") in pairs and ("log_s0", "log_phi") in pairs


def test_central_jacobian_matches_linear_design():
    X = np.column_stack([np.ones(10), np.arange(10.0)])
    J = dg.jacobian_central(lambda b: X @ b - 1.0, np.array([0.3, -0.2]))
    assert np.allclose(J, X)


def test_aic_negbin_counts_k():
    y = np.array([3.0, 5.0, 8.0])
    mu = np.array([3.5, 5.0, 7.0])
    a1 = dg.aic_negbin(y, mu, 50.0, n_params=2, k_estimated=True)
    a0 = dg.aic_negbin(y, mu, 50.0, n_params=2, k_estimated=False)
    assert np.isclose(a1 - a0, 2.0)
