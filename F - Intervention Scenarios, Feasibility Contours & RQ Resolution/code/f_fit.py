import numpy as np
import pandas as pd
from scipy.optimize import least_squares
import observation as obs
from dengue_model import dengue_rhs, PARAMS, basic_reproduction_number

S0_LIST = [0.39, 0.5, 0.64]
LO, HI = np.log([0.05, 1e-5]), np.log([2000.0, 1.0])
WEEK_EDGES = np.array([7 * k for k in range(35)] + [246 + 7 * k for k in range(18)], float)

y = pd.read_csv("../inputs/e1_fit_curves.csv")["observed"].to_numpy(float)
base = dict(PARAMS)
assert obs.validate_rhs(dengue_rhs, base) == []


def make_resid(s0):
    def resid(theta):
        p = dict(base, m=float(np.exp(theta[0])))
        try:
            mu = obs.expected_cases(p, s0, float(np.exp(theta[1])), WEEK_EDGES, rhs=dengue_rhs, y_anchor=y)
        except (obs.ObservationError, ValueError):
            return 1e3 * np.ones(len(y))
        return obs.residuals(y, mu, kind="sqrt")
    return resid


def fit(s0):
    resid = make_resid(s0)
    best = None
    for m0 in (4 / s0, 2 / s0, 8 / s0):
        x0 = np.clip(np.log([m0, 0.002 / s0]), LO + 1e-9, HI - 1e-9)
        r = least_squares(resid, x0, bounds=(LO, HI), x_scale="jac", diff_step=1e-4, xtol=1e-10, ftol=1e-10)
        if best is None or r.cost < best.cost:
            best = r
    m, phi = np.exp(best.x)
    p = dict(base, m=float(m))
    mu = obs.expected_cases(p, s0, phi, WEEK_EDGES, rhs=dengue_rhs, y_anchor=y)
    R0 = basic_reproduction_number(p)[1]
    return dict(s0=s0, m=m, phi=phi, R0=R0, R_eff=R0 * np.sqrt(s0), m_x_s0=m * s0, phi_x_s0=phi * s0,
                RSS=2 * best.cost, cost=best.cost, model_peak=mu.max(), model_total=mu.sum(),
                peak_week_model=int(mu.argmax()) + 1, solver_status=int(best.status))


if __name__ == "__main__":
    res = pd.DataFrame([fit(s0) for s0 in S0_LIST])
    res.round(6).to_csv("../outputs/f_fit_results.csv", index=False)
    print(res.round(6).to_string(index=False))
