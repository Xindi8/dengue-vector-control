"""Tests for E2's observation module. Run from the e2 folder:  python -m pytest -q"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import observation as obs                                   # noqa: E402
from model_provisional import baseline_params, dengue_rhs  # noqa: E402

CAL = os.path.join(ROOT, "data", "02_calibration_san_juan_2010_11.csv")


@pytest.fixture(scope="module")
def setup():
    cal = obs.load_calibration(CAL)
    edges = obs.week_edges(cal)
    p = dict(baseline_params(), m=10.0)
    y = cal.denv1_attributed_cases.to_numpy(float)
    return cal, edges, p, y


def test_week_edges_keep_the_8_day_week(setup):
    _, edges, _, _ = setup
    widths = np.diff(edges)
    assert len(widths) == 52 and edges[-1] == 365
    assert widths[34] == 8 and np.sum(widths == 8) == 1      # week 35: 24-31 Dec 2010


def test_week_edges_reject_bad_weeks(setup):
    cal, _, _, _ = setup
    bad = cal.copy()
    bad.loc[10, "days_since_start"] += 5
    with pytest.raises(ValueError):
        obs.week_edges(bad)


def test_provisional_model_passes_the_contract(setup):
    _, _, p, _ = setup
    assert obs.validate_rhs(dengue_rhs, p) == []


def test_contract_catches_a_wrong_state_order(setup):
    _, _, p, _ = setup

    def swapped(t, y, params):          # D accidentally orders vectors as (S_v, I_v, E_v)
        y2 = list(y)
        y2[5], y2[6] = y[6], y[5]
        d = dengue_rhs(t, y2, params)
        d[5], d[6] = d[6], d[5]
        return d

    assert obs.validate_rhs(swapped, p) != []


def test_initial_state_is_consistent(setup):
    _, _, p, _ = setup
    y0 = obs.initial_state(p, 0.4, 500.0)
    assert np.isclose(y0[0], 0.4 * p["N_h"])
    assert np.isclose(y0[:4].sum(), p["N_h"]) and np.isclose(y0[4:].sum(), p["m"] * p["N_h"])
    assert np.all(y0 >= 0) and np.isclose(p["eps_h"] * y0[1], 500.0)


def test_counter_equals_the_integral_of_onsets(setup):
    _, edges, p, _ = setup
    y0 = obs.initial_state(p, 0.4, 300.0)
    weekly, sol = obs.simulate_weekly_onsets(dengue_rhs, p, y0, edges, return_solution=True)
    from scipy.integrate import solve_ivp
    fine = solve_ivp(lambda t, y: dengue_rhs(t, y, p), (0, 365), y0, method="LSODA",
                     rtol=1e-10, atol=1e-8, dense_output=True)
    for w in (0, 15, 34, 51):          # includes the 8-day week
        t = np.linspace(edges[w], edges[w + 1], 4001)
        trap = getattr(np, "trapezoid", None) or np.trapz      # numpy < 2 has only trapz
        integral = trap(p["eps_h"] * fine.sol(t)[1], t)
        assert np.isclose(weekly[w], integral, rtol=1e-4)
    assert np.allclose(sol.y[:4].sum(0), p["N_h"], rtol=1e-9)


def test_growth_mode_start_grows_at_the_dominant_rate(setup):
    _, _, p, _ = setup
    r, _ = obs.growth_mode(p, 0.4)
    days = np.arange(0, 8, dtype=float)
    weekly = obs.simulate_weekly_onsets(dengue_rhs, p, obs.initial_state(p, 0.4, 1.0), days)
    rates = np.log(weekly[1:] / weekly[:-1])
    assert np.allclose(rates, r, atol=2e-4)                  # no start-up transient


def test_anchoring_matches_the_first_weeks(setup):
    _, edges, p, y = setup
    mu = obs.expected_cases(p, 0.4, 0.005, edges, rhs=dengue_rhs, y_anchor=y, n_anchor=3)
    assert np.isclose(mu[:3].sum(), y[:3].sum(), rtol=1e-4)


def test_ridge_symmetry(setup):
    """(s0, m, phi) and (s0/2, 2m, 2phi) give nearly the same weekly cases.
    Births and the ~2% of infected mosquitoes break the symmetry slightly,
    mostly in the low-count tail, so we check the main wave and the total."""
    _, edges, p, y = setup
    mu1 = obs.expected_cases(p, 0.4, 0.005, edges, rhs=dengue_rhs, y_anchor=y)
    mu2 = obs.expected_cases(dict(p, m=20.0), 0.2, 0.010, edges, rhs=dengue_rhs, y_anchor=y)
    assert np.max(np.abs(mu2[:30] / mu1[:30] - 1)) < 0.05
    assert abs(mu2.sum() / mu1.sum() - 1) < 0.05


def test_residual_definitions():
    y = np.array([0.0, 4.0, 100.0])
    for kind, k in (("sqrt", None), ("nb_vst", 50.0), ("log1p", None), ("raw", None), ("pearson", 50.0)):
        assert np.allclose(obs.residuals(y, y, kind=kind, k=k), 0)
    assert np.allclose(obs.nb_vst(y, 1e12), np.sqrt(y + 0.375), rtol=1e-6)
    assert obs.residuals([10.0], [4.0], kind="sqrt")[0] > 0   # positive = model too low


def test_nb_vst_stabilises_variance():
    rng = np.random.default_rng(1)
    k = 30.0
    for mu in (5.0, 50.0, 300.0):
        draws = rng.negative_binomial(k, k / (k + mu), 200_000)
        assert abs(np.var(obs.nb_vst(draws, k)) - 0.25) < 0.02


def test_scores_and_baseline():
    y = np.array([1, 3, 9, 20, 15, 8, 4, 2], dtype=float)
    split = next(obs.rolling_origin_splits(len(y), cutoffs=(4,)))
    s = obs.score_split(y, y, split)
    assert s["peak_week_error"] == 0 and abs(s["test_total_rel_error"]) < 1e-12 and s["test_mae"] == 0
    late = y.copy()
    late[5] = 30.0
    assert obs.score_split(y, late, split)["peak_week_error"] == 2
    assert len(obs.persistence_forecast(y[:4], 4)) == 4


def test_errors_are_raised_not_crashed(setup):
    _, edges, p, y = setup
    with pytest.raises(ValueError):
        obs.initial_state(p, 0.0, 10.0)
    with pytest.raises(obs.ObservationError):
        obs.initial_state(p, 0.4, 1e9)                       # more infected than people
    with pytest.raises(ValueError):
        obs.expected_cases(p, 0.4, 0.0, edges, rhs=dengue_rhs, y_anchor=y)
