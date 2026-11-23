"""Shared 7-state SEIR_h-SEI_v dengue model.

Role D (Xindi): model code. Role C (Aileen): equations, disease-free equilibrium and
next-generation R0 (docs/C_theory_derivation_v0.1.pdf; section numbers in the comments).

State order y = [S_h, E_h, I_h, R_h, S_v, E_v, I_v]; time in days; compartments are counts.
Parameter keys: N_h, m, a0, kappa, eta, b, c, eps_h, eps_v, gamma_h, d_h, d_v (rates per day).
There is no "a" key: the biting rate is a = a0 * (1 - kappa * eta).

Baseline values are read from the AB package (data_and_literature/parameters/parameters.json),
so the parameter table there is the single source of truth.
"""
from __future__ import annotations

import json

import numpy as np

from .paths import PARAMETERS_JSON

STATE_NAMES = ["S_h", "E_h", "I_h", "R_h", "S_v", "E_v", "I_v"]
PARAM_KEYS = ("N_h", "m", "a0", "kappa", "eta", "b", "c", "eps_h", "eps_v", "gamma_h", "d_h", "d_v")


def load_params(path=PARAMETERS_JSON) -> dict:
    """Baseline parameters with exactly the model's keys."""
    with open(path, encoding="utf-8") as f:
        params = json.load(f)["PARAMS"]
    missing = [k for k in PARAM_KEYS if k not in params]
    if missing:
        raise KeyError(f"parameters.json is missing {missing}")
    return {k: float(params[k]) for k in PARAM_KEYS}


PARAMS = load_params()


def biting_rate(params: dict) -> float:
    """Effective biting rate a = a0 (1 - kappa eta)   [Theory 1.1]"""
    return params["a0"] * (1.0 - params["kappa"] * params["eta"])


def dengue_rhs(t, y, params):
    """Right-hand side for solve_ivp(dengue_rhs, (0, t_end), y0, args=(params,))."""
    S_h, E_h, I_h, R_h, S_v, E_v, I_v = y
    p = params
    b, c, m = p["b"], p["c"], p["m"]
    eps_h, eps_v, gamma_h = p["eps_h"], p["eps_v"], p["gamma_h"]
    d_h, d_v, N_h = p["d_h"], p["d_v"], p["N_h"]
    N_v = m * N_h                                    # N_v = m N_h   [Theory 2]
    a = biting_rate(p)                               # [Theory 1.1]

    # Governing equations   [Theory 1]
    dS_h = d_h * N_h - a * b * (I_v / N_h) * S_h - d_h * S_h
    dE_h = a * b * (I_v / N_h) * S_h - (eps_h + d_h) * E_h
    dI_h = eps_h * E_h - (gamma_h + d_h) * I_h
    dR_h = gamma_h * I_h - d_h * R_h
    dS_v = d_v * N_v - a * c * (I_h / N_h) * S_v - d_v * S_v
    dE_v = a * c * (I_h / N_h) * S_v - (eps_v + d_v) * E_v
    dI_v = eps_v * E_v - d_v * I_v
    return np.array([dS_h, dE_h, dI_h, dR_h, dS_v, dE_v, dI_v])


def disease_free_equilibrium(params: dict) -> np.ndarray:
    """E_0 = (N_h, 0, 0, 0, m N_h, 0, 0)   [Theory 3]"""
    N_h = params["N_h"]
    return np.array([N_h, 0.0, 0.0, 0.0, params["m"] * N_h, 0.0, 0.0])


def basic_reproduction_number(params: dict):
    """(spectral radius of the next-generation matrix, closed form)   [Theory 4-7]"""
    p = params
    a = biting_rate(p)
    b, c, m = p["b"], p["c"], p["m"]
    eps_h, eps_v, gamma_h = p["eps_h"], p["eps_v"], p["gamma_h"]
    d_h, d_v = p["d_h"], p["d_v"]
    # Jacobians at the DFE, infected order x = (E_h, I_h, E_v, I_v)   [Theory 4.3]
    F = np.array([[0, 0, 0, a * b],
                  [0, 0, 0, 0],
                  [0, a * c * m, 0, 0],
                  [0, 0, 0, 0]])
    V = np.array([[eps_h + d_h, 0, 0, 0],
                  [-eps_h, gamma_h + d_h, 0, 0],
                  [0, 0, eps_v + d_v, 0],
                  [0, 0, -eps_v, d_v]])
    K = F @ np.linalg.inv(V)                          # next-generation matrix   [Theory 5]
    r0_ngm = float(np.max(np.abs(np.linalg.eigvals(K))))
    return r0_ngm, r0(params)


def r0(params: dict) -> float:
    """Closed-form R0   [Theory 7]. Scales with a and with sqrt(m) and sqrt(b c)."""
    p = params
    a = biting_rate(p)
    return float(np.sqrt(a ** 2 * p["b"] * p["c"] * p["m"] * p["eps_h"] * p["eps_v"]
                         / ((p["d_h"] + p["eps_h"]) * (p["d_h"] + p["gamma_h"])
                            * (p["d_v"] + p["eps_v"]) * p["d_v"])))


def r_effective(params: dict, s0: float) -> float:
    """R with a share s0 of humans susceptible: R0 * sqrt(s0) (only the mosquito-to-human
    step meets the human susceptible fraction)."""
    return r0(params) * float(np.sqrt(s0))
