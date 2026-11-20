import numpy as np
from dengue_model import PARAMS, dengue_rhs, basic_reproduction_number
import observation as obs

BASE = dict(PARAMS)

FIT_FILES = ("../inputs/e1_fit_results.csv", "../outputs/f_fit_results.csv")


def load_fit(paths=FIT_FILES):
    fit = {}
    for path in paths:
        rows = np.genfromtxt(path, delimiter=",", names=True)
        fit.update({round(float(r["s0"]), 2): (float(r["m"]), float(r["phi"])) for r in rows})
    return fit


FIT = load_fit()

WEEK_EDGES = np.array([7 * k for k in range(35)] + [246 + 7 * k for k in range(18)], float)


def r0(p):
    return basic_reproduction_number(p)[1]


def anchored_start(s0, observed):
    m, phi = FIT[s0]
    p = dict(BASE, m=m)
    rate = obs.anchored_onset_rate(dengue_rhs, p, s0, phi, WEEK_EDGES, observed)
    return p, obs.initial_state(p, s0, rate), phi


def weekly_cases(p, y0, phi, edges=WEEK_EDGES):
    return phi * obs.simulate_weekly_onsets(dengue_rhs, p, y0, edges)


def apply_controls(p, y0, kappa=0.0, delta=0.0, p_w=0.0, eta=0.5, q=0.9):
    pc = dict(p, kappa=kappa, eta=eta, m=p["m"] * (1 - delta))
    shrink = np.sqrt(1 - p_w * q)
    pc["b"], pc["c"] = p["b"] * shrink, p["c"] * shrink
    yc = np.array(y0, dtype=float)
    yc[4:7] *= (1 - delta)
    return pc, yc
