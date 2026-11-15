import numpy as np
from scipy.integrate import solve_ivp
from dengue_model import PARAMS, dengue_rhs, basic_reproduction_number

BASE = dict(PARAMS)

FIT_FILE = "../inputs/e1_fit_results.csv"


def load_fit(path=FIT_FILE):
    rows = np.genfromtxt(path, delimiter=",", names=True)
    return {round(float(r["s0"]), 2): (float(r["m"]), float(r["phi"])) for r in rows}


FIT = load_fit()

WEEK_EDGES = np.array([7 * k for k in range(35)] + [246 + 7 * k for k in range(18)], float)


def biting(p):
    return p["a0"] * (1 - p["kappa"] * p["eta"])


def rhs_with_onsets(t, y, p):
    return np.append(dengue_rhs(t, y[:7], p), p["eps_h"] * y[1])


def r0(p):
    return basic_reproduction_number(p)[1]


def growth_mode(p, s0):
    a = biting(p)
    J = np.array([
        [-(p["eps_h"] + p["d_h"]), 0, 0, a * p["b"] * s0],
        [p["eps_h"], -(p["gamma_h"] + p["d_h"]), 0, 0],
        [0, a * p["c"] * p["m"], -(p["eps_v"] + p["d_v"]), 0],
        [0, 0, p["eps_v"], -p["d_v"]]])
    vals, vecs = np.linalg.eig(J)
    k = np.argmax(vals.real)
    v = np.abs(vecs[:, k].real)
    return v / v[0], vals[k].real


def initial_state(p, s0, scale):
    v, _ = growth_mode(p, s0)
    E_h, I_h, E_v, I_v = scale * v * p["N_h"]
    N_h, N_v = p["N_h"], p["m"] * p["N_h"]
    return np.array([s0 * N_h - E_h - I_h, E_h, I_h, (1 - s0) * N_h,
                     N_v - E_v - I_v, E_v, I_v, 0.0])


def weekly_cases(p, y0, phi, edges=WEEK_EDGES):
    sol = solve_ivp(rhs_with_onsets, (edges[0], edges[-1]), y0, args=(p,), t_eval=edges,
                    method="LSODA", rtol=1e-9, atol=1e-6)
    if not sol.success:
        raise RuntimeError(sol.message)
    return phi * np.diff(sol.y[7]), sol


def anchored_start(s0, observed_first3):
    m, phi = FIT[s0]
    p = dict(BASE, m=m)
    scale = 1e-6
    for _ in range(20):
        got = weekly_cases(p, initial_state(p, s0, scale), phi)[0][:3].sum()
        if abs(got / observed_first3 - 1) < 1e-10:
            break
        scale *= observed_first3 / got
    return p, initial_state(p, s0, scale), phi


def apply_controls(p, y0, kappa=0.0, delta=0.0, p_w=0.0, eta=0.5, q=0.9):
    pc = dict(p, kappa=kappa, eta=eta, m=p["m"] * (1 - delta))
    shrink = np.sqrt(1 - p_w * q)
    pc["b"], pc["c"] = p["b"] * shrink, p["c"] * shrink
    yc = y0.copy()
    yc[4:7] *= (1 - delta)
    return pc, yc
