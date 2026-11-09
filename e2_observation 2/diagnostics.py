"""E2 - residual diagnosis module (MATH 574, Team 2).

The team plan asks E2 for a module that reports RSS, AIC and the parameter
correlation matrix for a fit. It works directly on E1's output:

    res = scipy.optimize.least_squares(resid, x0, ...)
    rep = fit_report(res, param_names=["log_m", "log_phi"])
    print(format_report(rep))

What it reports
- RSS, residual variance sigma2 = RSS / (n - p)
- AIC, AICc, BIC on the residual scale (Gaussian form). Compare these only
  between fits to the SAME weeks on the SAME residual scale (e.g. square root).
- covariance = sigma2 * (J^T J)^-1, standard errors, correlation matrix
- condition number of the Jacobian, and every parameter pair with |corr| >= 0.95,
  which means the data cannot tell those parameters apart
- lag-1 autocorrelation and runs test of the residuals

If E1 fits parameters on a log scale, the correlations are the same as on the
natural scale (to first order); a natural-scale standard error is
value * (log-scale standard error).

Important for E1: with an ODE inside the residuals, least_squares' default
finite-difference step (~1e-8) is smaller than the solver's own error, so
res.jac can be badly wrong. On a synthetic fit it hid the s0-m-phi ridge
completely (condition number 7 instead of about 950). Either fit with
diff_step=1e-4, or pass `fun=` here so the Jacobian is recomputed by central
differences with a safe step.
"""
from __future__ import annotations

import numpy as np

from observation import lag1_autocorrelation, nb_negloglik, runs_test_z


def jacobian_central(fun, x, step: float = 1e-4):
    """Central-difference Jacobian of the residual function at x (absolute step,
    suited to parameters on a log scale)."""
    x = np.asarray(x, dtype=float)
    cols = []
    for i in range(len(x)):
        e = np.zeros_like(x)
        e[i] = step
        cols.append((np.asarray(fun(x + e), float) - np.asarray(fun(x - e), float)) / (2 * step))
    return np.column_stack(cols)


def fit_report(res, param_names=None, fun=None, step: float = 1e-4,
               corr_flag: float = 0.95, rcond: float = 1e-12):
    """Residual diagnosis for a scipy.optimize.least_squares result.

    Give `fun` (the residual function E1 minimised) to recompute residuals and
    the Jacobian at res.x by central differences; otherwise res.fun and res.jac
    are used as they are (only safe if the fit used diff_step around 1e-4)."""
    if fun is not None:
        r = np.asarray(fun(res.x), dtype=float)
        J = jacobian_central(fun, res.x, step)
    else:
        r = np.asarray(res.fun, dtype=float)
        J = np.atleast_2d(np.asarray(res.jac, dtype=float))
    n, p = J.shape
    if n <= p:
        raise ValueError("need more residuals than parameters")
    names = list(param_names) if param_names is not None else [f"theta{i}" for i in range(p)]
    rss = float(r @ r)
    sigma2 = rss / (n - p)
    loglik_term = n * np.log(rss / n)
    aic = loglik_term + 2 * p
    aicc = aic + 2 * p * (p + 1) / (n - p - 1) if n - p - 1 > 0 else np.inf
    bic = loglik_term + p * np.log(n)

    U, s, Vt = np.linalg.svd(J, full_matrices=False)
    cond = float(s[0] / s[-1]) if s[-1] > 0 else np.inf
    singular = bool(s[-1] <= rcond * s[0])
    inv_s2 = np.where(s > rcond * s[0], 1.0 / np.maximum(s, 1e-300) ** 2, 0.0)   # pseudo-inverse
    cov = sigma2 * (Vt.T * inv_s2) @ Vt
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = cov / np.outer(se, se)
    np.fill_diagonal(corr, 1.0)
    flagged = [(names[i], names[j], float(corr[i, j]))
               for i in range(p) for j in range(i + 1, p)
               if np.isfinite(corr[i, j]) and abs(corr[i, j]) >= corr_flag]
    return dict(n=n, p=p, names=names, rss=rss, sigma2=sigma2, aic=float(aic), aicc=float(aicc),
                bic=float(bic), cov=cov, se=se, corr=corr, condition_number=cond,
                singular=singular, flagged_pairs=flagged,
                lag1_autocorrelation=lag1_autocorrelation(r), runs_test_z=runs_test_z(r))


def aic_negbin(y, mu, k, n_params: int, k_estimated: bool = True):
    """Likelihood-based AIC for a negative-binomial fit (for MCMC/likelihood work).
    n_params counts the model parameters; k adds one more if it was estimated."""
    p = n_params + (1 if k_estimated else 0)
    return 2 * p + 2 * float(nb_negloglik(y, mu, k))


def format_report(rep: dict, digits: int = 3) -> str:
    """Plain-text summary for run logs."""
    names = rep["names"]
    lines = [f"n = {rep['n']}, p = {rep['p']}, RSS = {rep['rss']:.{digits}f}, "
             f"sigma2 = {rep['sigma2']:.{digits}f}",
             f"AIC = {rep['aic']:.{digits}f}, AICc = {rep['aicc']:.{digits}f}, BIC = {rep['bic']:.{digits}f}",
             f"condition number = {rep['condition_number']:.3g}" + ("  (SINGULAR)" if rep["singular"] else ""),
             "standard errors: " + ", ".join(f"{nm} {s:.{digits}g}" for nm, s in zip(names, rep["se"])),
             "correlation matrix:"]
    width = max(len(nm) for nm in names) + 2
    lines.append(" " * width + "".join(f"{nm:>{width}}" for nm in names))
    for i, nm in enumerate(names):
        lines.append(f"{nm:<{width}}" + "".join(f"{rep['corr'][i, j]:>{width}.{digits}f}" for j in range(len(names))))
    if rep["flagged_pairs"]:
        lines.append("NOT SEPARATELY IDENTIFIABLE (|corr| >= 0.95): " +
                     "; ".join(f"{a}-{b} {c:+.3f}" for a, b, c in rep["flagged_pairs"]))
    lines.append(f"residuals: lag-1 autocorrelation {rep['lag1_autocorrelation']:+.3f}, "
                 f"runs test z {rep['runs_test_z']:+.2f}")
    return "\n".join(lines)
