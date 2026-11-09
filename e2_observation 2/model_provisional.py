"""Copy of D's 7-state model for E2's tests (MATH 574, Team 2).

Same equations, state order and parameter keys as Xindi's dengue_model.py
and Aileen's derivation, so E2's tests and checks run on their own:

    dengue_rhs(t, y, params) -> 7 derivatives

For team runs, import ``dengue_rhs`` from dengue_model.py instead; nothing
else changes.

State order (index: name)
    0 S_h   susceptible humans
    1 E_h   exposed humans (infected, not yet infectious)
    2 I_h   infectious humans
    3 R_h   recovered humans
    4 S_v   susceptible female mosquitoes
    5 E_v   exposed mosquitoes
    6 I_v   infectious mosquitoes

Parameter keys (all rates per day; names match 06_parameters.csv)
    N_h      human population (persons)
    m        female mosquitoes per person, N_v = m * N_h
    a        effective biting rate per mosquito (1/day); if a0, kappa, eta
             are given instead, a = a0 * (1 - kappa * eta) as in the proposal
    b        mosquito-to-human transmission probability per bite
    c        human-to-mosquito transmission probability per bite
    eps_h    intrinsic incubation rate (1/latent period in humans)
    eps_v    extrinsic incubation rate (1/latent period in mosquitoes)
    gamma_h  human recovery rate (1/infectious period)
    d_h      human birth = death rate (b_h = d_h keeps N_h constant)
    d_v      mosquito emergence = death rate (b_v = d_v keeps N_v constant)
"""
from __future__ import annotations

import numpy as np

STATE_NAMES = ("S_h", "E_h", "I_h", "R_h", "S_v", "E_v", "I_v")
IDX = {name: i for i, name in enumerate(STATE_NAMES)}

REQUIRED_KEYS = ("N_h", "m", "b", "c", "eps_h", "eps_v", "gamma_h", "d_h", "d_v")


def effective_biting_rate(params: dict) -> float:
    """a = a0 * (1 - kappa * eta) if the compliance keys are present, else params['a']."""
    if "a0" in params:
        return params["a0"] * (1.0 - params.get("kappa", 0.0) * params.get("eta", 0.0))
    return params["a"]


def baseline_params() -> dict:
    """Table 1 baselines (06_parameters.csv), with D's keys (a0, kappa, eta; no "a").
    Provisional; owned by AB/C/D."""
    return dict(
        N_h=2_347_833.0,      # San Juan-Carolina-Caguas MSA, 2010 Census
        m=2.0,                # fitted in practice; prior 1-10 in Table 1
        a0=0.35,              # baseline biting rate, Zahid et al. 2023
        kappa=0.0,            # compliance with personal protection (0 = none)
        eta=0.0,              # protection efficacy per contact (0 = none); a = a0 * (1 - kappa * eta)
        b=0.33,
        c=0.33,
        eps_h=1 / 5.9,        # Chan & Johansson 2012
        eps_v=1 / 10.7,       # Chan & Johansson 2012, 27 C
        gamma_h=1 / 5.0,
        d_h=1 / (78.7 * 365),
        d_v=1 / 14.0,
    )


def dengue_rhs(t: float, y, params: dict):
    """Right-hand side of the 7-state host-vector model (proposal Section 3)."""
    S_h, E_h, I_h, R_h, S_v, E_v, I_v = y
    N_h = params["N_h"]
    N_v = params["m"] * N_h
    a = effective_biting_rate(params)
    d_h, d_v = params["d_h"], params["d_v"]

    inf_h = a * params["b"] * (I_v / N_h) * S_h   # new human infections per day
    inf_v = a * params["c"] * (I_h / N_h) * S_v   # new mosquito infections per day

    return [
        d_h * N_h - inf_h - d_h * S_h,
        inf_h - (params["eps_h"] + d_h) * E_h,
        params["eps_h"] * E_h - (params["gamma_h"] + d_h) * I_h,
        params["gamma_h"] * I_h - d_h * R_h,
        d_v * N_v - inf_v - d_v * S_v,
        inf_v - (params["eps_v"] + d_v) * E_v,
        params["eps_v"] * E_v - d_v * I_v,
    ]


def r0(params: dict) -> float:
    """Closed-form R0 from the proposal (Section 5, next-generation matrix, square-root form)."""
    a = effective_biting_rate(params)
    p = params
    num = a * a * p["b"] * p["c"] * p["m"] * p["eps_h"] * p["eps_v"]
    den = (p["d_h"] + p["eps_h"]) * (p["d_h"] + p["gamma_h"]) * (p["d_v"] + p["eps_v"]) * p["d_v"]
    return float(np.sqrt(num / den))
