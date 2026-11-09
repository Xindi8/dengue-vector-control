"""E2 - observation model, residuals and holdout tools (MATH 574, Team 2).

Owner: Fahad (E2). Interface to be agreed with D (model) and E1 (fitting).

What this module does
---------------------
1. Turns the ODE solution into the quantity the surveillance data count:
   new symptom onsets per reporting week, using the REAL week boundaries
   (one week a year is 8-9 days long in the CDC San Juan data).
2. Scales onsets by the reporting fraction phi to get expected weekly cases.
3. Starts the epidemic on the model's own dominant growth mode, scaled to the
   first observed weeks, so there is no artificial start-up dip.
4. Defines residuals and likelihoods for fitting (E1) and diagnostics (Stage 3).
5. Defines the holdout split and forecast scores for validation (Stage 3).

Interface with D
----------------
D provides ``dengue_rhs(t, y, params)`` with the 7-state order
(S_h, E_h, I_h, R_h, S_v, E_v, I_v) and the parameter keys listed in
``model_provisional.py``. Every function here takes ``rhs`` as an argument,
so swapping the provisional model for D's needs no other change.

Time-varying rates: if D makes eps_v follow temperature (file 02 carries a
weekly extrinsic incubation period), give growth_mode and validate_rhs the
week-1 value; everything else works unchanged.

Interface with E1
-----------------
E1's residual function for scipy.optimize.least_squares is

    def resid(theta):
        params, s0, phi = unpack(theta)             # E1 decides which are free
        mu = expected_cases(params, s0, phi, edges, y_anchor=y_train, rhs=rhs)
        return residuals(y_train, mu[train], kind="sqrt")

and E1 should catch ``ObservationError`` (solver failure, negative states)
and return a large finite residual vector instead of crashing.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.special import gammaln

STATE_NAMES = ("S_h", "E_h", "I_h", "R_h", "S_v", "E_v", "I_v")
E_H, I_H, E_V, I_V = 1, 2, 5, 6


class ObservationError(RuntimeError):
    """Raised when the ODE solve fails or returns impossible values."""


# ---------------------------------------------------------------------------
# 1. Data and week boundaries
# ---------------------------------------------------------------------------
CASE_COLUMNS = {
    # column in 02_calibration_san_juan_2010_11.csv : meaning
    "denv1_attributed_cases": "PRIMARY. Lab-positive cases attributed to DENV-1 "
                              "(untyped positives split by that week's serotype mix).",
    "denv1_cases": "SENSITIVITY. Lab-positive cases typed as DENV-1 by RT-PCR only.",
    "total_cases": "SENSITIVITY. All lab-positive cases, all serotypes.",
}


def load_calibration(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, parse_dates=["week_start_date"])
    if df["days_since_start"].iloc[0] != 0:
        raise ValueError("days_since_start must start at 0")
    return df


def week_edges(df: pd.DataFrame, last_week_days: int = 7) -> np.ndarray:
    """Week boundaries in days since the first week start, including the end
    of the last week. Uses the real dates, so the 8-9 day year-end week keeps
    its true length. Checks that every week is 7-9 days long."""
    starts = df["days_since_start"].to_numpy(dtype=float)
    edges = np.append(starts, starts[-1] + last_week_days)
    widths = np.diff(edges)
    if np.any(widths < 7) or np.any(widths > 9):
        raise ValueError(f"unexpected week lengths: {sorted(set(widths))}")
    return edges


# ---------------------------------------------------------------------------
# 2. Initial state on the dominant growth mode
# ---------------------------------------------------------------------------
def _biting_rate(params: dict) -> float:
    if "a0" in params:
        return params["a0"] * (1.0 - params.get("kappa", 0.0) * params.get("eta", 0.0))
    return params["a"]


def infected_jacobian(params: dict, s0: float) -> np.ndarray:
    """Linearised infected subsystem (E_h, I_h, E_v, I_v) at t = 0, when a
    fraction s0 of humans is susceptible and almost all mosquitoes are."""
    p, a = params, _biting_rate(params)
    return np.array([
        [-(p["eps_h"] + p["d_h"]), 0.0, 0.0, a * p["b"] * s0],
        [p["eps_h"], -(p["gamma_h"] + p["d_h"]), 0.0, 0.0],
        [0.0, a * p["c"] * p["m"], -(p["eps_v"] + p["d_v"]), 0.0],
        [0.0, 0.0, p["eps_v"], -p["d_v"]],
    ])


def growth_mode(params: dict, s0: float):
    """Dominant eigenvalue r (per day) and its eigenvector, scaled so that the
    symptom-onset rate eps_h * E_h equals 1 person per day."""
    vals, vecs = np.linalg.eig(infected_jacobian(params, s0))
    k = int(np.argmax(vals.real))
    r, v = float(vals[k].real), np.abs(vecs[:, k].real)   # Perron vector is positive
    v = v / (params["eps_h"] * v[0])
    return r, v


def initial_state(params: dict, s0: float, onset_rate0: float) -> np.ndarray:
    """7-state initial condition.

    s0           fraction of humans susceptible to the epidemic serotype at t = 0
                 (S_h(0) = s0 * N_h exactly)
    onset_rate0  symptom onsets per day at t = 0 (all infections, before phi)
    """
    if not 0.0 < s0 <= 1.0:
        raise ValueError("s0 must be in (0, 1]")
    N_h, N_v = params["N_h"], params["m"] * params["N_h"]
    _, v = growth_mode(params, s0)
    E_h, I_h, E_v, I_v = onset_rate0 * v
    S_h = s0 * N_h
    R_h = N_h - S_h - E_h - I_h
    if R_h < 0:                      # s0 close to 1: take the infected out of S_h
        S_h, R_h = N_h - E_h - I_h, 0.0
    S_v = N_v - E_v - I_v
    if S_h < 0 or S_v < 0:
        raise ObservationError("initial infected exceed the population")
    return np.array([S_h, E_h, I_h, R_h, S_v, E_v, I_v])


def validate_rhs(rhs, params: dict, s0: float = 0.5, rel_tol: float = 1e-6):
    """Contract test for D's ``dengue_rhs(t, y, params)``. Returns a list of
    problems; an empty list means the model plugs into E1/E2 unchanged.

    Checks: 7 derivatives; the disease-free state is an equilibrium; human and
    mosquito totals are conserved; and the linearised infected block matches
    the Jacobian E2 uses for the starting state (this catches a different
    state order or a different force of infection)."""
    problems = []
    N_h, N_v = params["N_h"], params["m"] * params["N_h"]
    dfe = np.array([N_h, 0, 0, 0, N_v, 0, 0], dtype=float)
    d = np.asarray(rhs(0.0, dfe, params), dtype=float)
    if d.shape != (7,):
        return ["rhs must return 7 derivatives in the order " + ", ".join(STATE_NAMES)]
    if np.max(np.abs(d)) > rel_tol * max(N_h, N_v):
        problems.append("the disease-free state is not an equilibrium")

    base = np.array([s0 * N_h, 0, 0, (1 - s0) * N_h, N_v, 0, 0], dtype=float)
    y = base + np.array([0, 50, 40, -90, -300, 200, 100], dtype=float)
    d = np.asarray(rhs(0.0, y, params), dtype=float)
    if abs(d[:4].sum()) > rel_tol * N_h:
        problems.append("human total is not conserved (births should equal deaths, b_h = d_h)")
    if abs(d[4:].sum()) > rel_tol * N_v:
        problems.append("mosquito total is not conserved (emergence should equal death, b_v = d_v)")

    # numerical Jacobian of the infected block at the partly immune DFE
    idx = [E_H, I_H, E_V, I_V]
    J = np.zeros((4, 4))
    f0 = np.asarray(rhs(0.0, base, params), dtype=float)
    for j, col in enumerate(idx):
        h = 1e-3
        yp = base.copy()
        yp[col] += h
        J[:, j] = (np.asarray(rhs(0.0, yp, params), dtype=float) - f0)[idx] / h
    J_ref = infected_jacobian(params, s0)
    scale = np.max(np.abs(J_ref))
    if np.max(np.abs(J - J_ref)) > 1e-4 * scale:
        problems.append("linearised infected block differs from E2's Jacobian: check the state "
                        "order and the force of infection (a*b*I_v/N_h*S_h, a*c*I_h/N_h*S_v)")
    return problems


# ---------------------------------------------------------------------------
# 3. Weekly onsets and expected reported cases
# ---------------------------------------------------------------------------
def simulate_weekly_onsets(rhs, params: dict, y0, edges, rtol=1e-8, atol=1e-6,
                           method="LSODA", return_solution=False):
    """New symptom onsets in each reporting week: integral of eps_h * E_h(t)
    over [edges[w], edges[w+1]). Uses an 8th 'counter' state C with
    dC/dt = eps_h * E_h, so the weekly totals are exact differences of C."""
    eps_h = params["eps_h"]

    def rhs8(t, y):
        dy = rhs(t, y[:7], params)
        return [*dy, eps_h * y[E_H]]

    y08 = np.append(np.asarray(y0, dtype=float), 0.0)
    try:
        sol = solve_ivp(rhs8, (edges[0], edges[-1]), y08, t_eval=edges,
                        method=method, rtol=rtol, atol=atol)
    except Exception as exc:                      # noqa: BLE001
        raise ObservationError(f"solver raised {exc!r}") from exc
    if not sol.success or sol.y.shape[1] != len(edges):
        raise ObservationError(f"solver failed: {sol.message}")
    if not np.all(np.isfinite(sol.y)):
        raise ObservationError("non-finite state")
    scale = max(params["N_h"], params["m"] * params["N_h"])
    if sol.y[:7].min() < -1e-6 * scale:
        raise ObservationError(f"negative state {sol.y[:7].min():.3g}")
    onsets = np.diff(sol.y[7])
    if onsets.min() < -1e-6:
        raise ObservationError("counter decreased")
    onsets = np.clip(onsets, 0.0, None)
    return (onsets, sol) if return_solution else onsets


def anchored_onset_rate(rhs, params, s0, phi, edges, y_anchor, n_anchor=3,
                        iterations=4, tol=1e-6):
    """Initial onset rate that makes the model's expected cases in the first
    n_anchor weeks equal the observed cases there.

    Early growth is almost linear in the starting size, so a unit run plus a
    rescale is close; a few rescaling passes on the short anchor window remove
    the small error from susceptible depletion (it matters when phi is small,
    because then many infections sit behind each reported case)."""
    y_anchor = np.asarray(y_anchor, dtype=float)[:n_anchor]
    target = y_anchor.sum()
    if target <= 0:
        raise ObservationError("anchor weeks have zero cases")
    short = edges[: n_anchor + 1]
    rate = 1.0
    for _ in range(iterations):
        model = phi * simulate_weekly_onsets(rhs, params, initial_state(params, s0, rate), short).sum()
        if model <= 0:
            raise ObservationError("model produces no onsets in the anchor weeks")
        ratio = target / model
        rate *= ratio
        if abs(ratio - 1.0) < tol:
            break
    return rate


def expected_cases(params: dict, s0: float, phi: float, edges, rhs,
                   y_anchor=None, n_anchor: int = 3, onset_rate0: float | None = None,
                   return_onsets=False):
    """Expected reported cases per week: mu_w = phi * onsets_w.

    Give either y_anchor (observed weekly cases; the first n_anchor weeks set the
    starting size, so no extra parameter is needed) or onset_rate0 (if E1 wants
    to fit the starting size as a free parameter instead)."""
    if not 0.0 < phi <= 1.0:
        raise ValueError("phi must be in (0, 1]")
    if onset_rate0 is None:
        if y_anchor is None:
            raise ValueError("give y_anchor or onset_rate0")
        onset_rate0 = anchored_onset_rate(rhs, params, s0, phi, edges, y_anchor, n_anchor)
    onsets = simulate_weekly_onsets(rhs, params, initial_state(params, s0, onset_rate0), edges)
    mu = phi * onsets
    return (mu, onsets) if return_onsets else mu


# ---------------------------------------------------------------------------
# 4. Residuals and likelihoods
# ---------------------------------------------------------------------------
def nb_vst(x, k):
    """Variance-stabilising (Anscombe) transform for negative-binomial counts
    with size k > 3/4:  g(x) = sqrt(k) * asinh(sqrt((x + 3/8) / (k - 3/4))).
    Var(g(Y)) is close to 1/4 at every mean, even for small counts.
    As k -> infinity it becomes sqrt(x + 3/8), the Poisson version."""
    if k <= 0.75:
        raise ValueError("nb_vst needs k > 0.75")
    x = np.asarray(x, dtype=float)
    return np.sqrt(k) * np.arcsinh(np.sqrt((np.clip(x, 0, None) + 0.375) / (k - 0.75)))


def residuals(y, mu, kind: str = "sqrt", k: float | None = None):
    """Residual = transformed data - transformed model (positive = model too low).

    kind
      "sqrt"     sqrt(y) - sqrt(mu): weights weeks evenly if noise is Poisson-like.
      "nb_vst"   nb_vst(y) - nb_vst(mu): weights weeks evenly under NB noise (needs k).
      "log1p"    log(y+1) - log(mu+1): suits purely multiplicative noise.
      "raw"      y - mu: lets the peak weeks dominate (not recommended for fitting).
      "pearson"  (y - mu) / sqrt(mu + mu^2/k): diagnostics only, not for fitting.
    """
    y = np.asarray(y, dtype=float)
    mu = np.asarray(mu, dtype=float)
    if kind == "sqrt":
        return np.sqrt(np.clip(y, 0, None)) - np.sqrt(np.clip(mu, 0, None))
    if kind == "nb_vst":
        if k is None:
            raise ValueError("nb_vst needs k")
        return nb_vst(y, k) - nb_vst(mu, k)
    if kind == "log1p":
        return np.log1p(y) - np.log1p(mu)
    if kind == "raw":
        return y - mu
    if kind == "pearson":
        if k is None:
            raise ValueError("pearson needs k (use np.inf for Poisson)")
        var = mu if np.isinf(k) else mu + mu ** 2 / k
        return (y - mu) / np.sqrt(np.clip(var, 1e-12, None))
    raise ValueError(f"unknown residual kind {kind!r}")


def nb_negloglik(y, mu, k):
    """Negative log-likelihood of NB(mean mu, size k); works for non-integer y."""
    y = np.asarray(y, dtype=float)
    mu = np.clip(np.asarray(mu, dtype=float), 1e-12, None)
    return -np.sum(gammaln(y + k) - gammaln(k) - gammaln(y + 1)
                   + k * np.log(k / (k + mu)) + y * np.log(mu / (k + mu)))


def poisson_negloglik(y, mu):
    y = np.asarray(y, dtype=float)
    mu = np.clip(np.asarray(mu, dtype=float), 1e-12, None)
    return -np.sum(y * np.log(mu) - mu - gammaln(y + 1))


def lag1_autocorrelation(r):
    r = np.asarray(r, dtype=float) - np.mean(r)
    return float(np.sum(r[1:] * r[:-1]) / np.sum(r * r))


def runs_test_z(r):
    """Wald-Wolfowitz runs test on residual signs. z < -2 means fewer sign
    changes than chance (long same-sign runs: the model misses a systematic
    pattern); z > +2 means the residuals flip sign more often than chance."""
    s = np.sign(np.asarray(r, dtype=float))
    s = s[s != 0]
    n1, n2 = np.sum(s > 0), np.sum(s < 0)
    runs = 1 + np.sum(s[1:] != s[:-1])
    n = n1 + n2
    mean = 2 * n1 * n2 / n + 1
    var = 2 * n1 * n2 * (2 * n1 * n2 - n) / (n ** 2 * (n - 1))
    return float((runs - mean) / np.sqrt(var))


# ---------------------------------------------------------------------------
# 5. Holdout split and forecast scores
# ---------------------------------------------------------------------------
@dataclass
class Split:
    cutoff: int          # number of training weeks (weeks 1..cutoff)
    train: np.ndarray    # indices 0..cutoff-1
    test: np.ndarray     # indices cutoff..n-1


def rolling_origin_splits(n_weeks: int, cutoffs=(12, 16, 20, 26)):
    """Expanding-window splits within one season: train on weeks 1..K, test on
    weeks K+1..n. The defaults sit before, at and after the 2010/11 peak
    (DENV-1 attributed cases peak in week 16); K = 26 is the headline split.
    Add K = 8 only as a stress test (growth phase only)."""
    for K in cutoffs:
        if not 3 <= K < n_weeks:
            raise ValueError(f"cutoff {K} out of range")
        yield Split(K, np.arange(K), np.arange(K, n_weeks))


def season_targets(y):
    """The three CDC 2015 forecasting-challenge targets (Johansson et al. 2019):
    peak week (1-based), peak weekly cases, and season total."""
    y = np.asarray(y, dtype=float)
    return {"peak_week": int(np.argmax(y)) + 1, "peak": float(y.max()), "total": float(y.sum())}


def score_split(y, mu, split: Split):
    """Scores for one split.

    Peak week and peak size combine the observed training weeks with the
    forecast test weeks, as in the CDC challenge; once the peak is inside the
    training weeks they are right by construction unless the forecast
    overshoots. The total is scored over the TEST weeks only, because a
    season total is mostly made of weeks the model has already seen."""
    y = np.asarray(y, dtype=float)
    mu = np.asarray(mu, dtype=float)
    combined = np.concatenate([y[split.train], mu[split.test]])
    truth, pred = season_targets(y), season_targets(combined)
    yt, mt = y[split.test], mu[split.test]
    return {
        "cutoff": split.cutoff,
        "peak_week_error": pred["peak_week"] - truth["peak_week"],
        "peak_rel_error": pred["peak"] / truth["peak"] - 1,
        "test_total_rel_error": mt.sum() / yt.sum() - 1,
        "season_total_rel_error": pred["total"] / truth["total"] - 1,
        "test_rmse_sqrt": float(np.sqrt(np.mean((np.sqrt(yt) - np.sqrt(mt)) ** 2))),
        "test_mae": float(np.mean(np.abs(yt - mt))),
    }


def persistence_forecast(y_train, horizon: int, window: int = 2):
    """Naive baseline: the mean of the last `window` observed weeks, held flat.
    The model is only worth its complexity if it beats this."""
    return np.full(horizon, float(np.mean(np.asarray(y_train, dtype=float)[-window:])))


# ---------------------------------------------------------------------------
# 6. Model-free dispersion check (choose the residual before calibration)
# ---------------------------------------------------------------------------
def local_dispersion(y, min_mean=5.0):
    """Model-free noise check. For counts with a smooth mean, the second
    difference e_w = y_w - (y_{w-1} + y_{w+1}) / 2 has variance about 1.5 * Var(y_w).
    Returns (local mean, e_w^2 / 1.5) pairs for weeks whose local mean >= min_mean."""
    y = np.asarray(y, dtype=float)
    e = y[1:-1] - 0.5 * (y[:-2] + y[2:])
    m = (y[:-2] + y[1:-1] + y[2:]) / 3.0
    keep = m >= min_mean
    return m[keep], (e[keep] ** 2) / 1.5


def estimate_nb_size(means, var_hat):
    """Fit Var = mu + mu^2 / k by least squares on (var - mu) = mu^2 * (1/k).
    Returns k (np.inf if no overdispersion)."""
    means, var_hat = np.asarray(means), np.asarray(var_hat)
    x, z = means ** 2, var_hat - means
    inv_k = float(np.sum(x * z) / np.sum(x * x))
    return np.inf if inv_k <= 0 else 1.0 / inv_k
