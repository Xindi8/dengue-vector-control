"""
dengue_ode.py — shared 7-state ODE (SEIR_h – SEI_v), MATH 574 Team 2
State order y = [S_h, E_h, I_h, R_h, S_v, E_v, I_v]; time unit = days.
"""
import numpy as np

STATE_NAMES = ["S_h", "E_h", "I_h", "R_h", "S_v", "E_v", "I_v"]


# PROVISIONAL values (proposal Table 1 baselines) — replace after B01/C01
PARAMS = {
    "N_h": 2.35e6,       # human population
    "m": 2.0,            # mosquitoes per human, N_v = m * N_h
    "a0": 0.35,          # baseline biting rate (1/day)
    "kappa": 0.0,        # community compliance with personal protection, in [0,1] (0 = no intervention)
    "eta": 0.0,          # protection efficacy per contact, in [0,1] (0 = no intervention)
    "b": 0.33,           # mosquito -> human transmission prob per bite
    "c": 0.33,           # human -> mosquito transmission prob per bite
    "eps_h": 1 / 5.9,    # intrinsic incubation rate (1/day)
    "eps_v": 1 / 10.7,   # extrinsic incubation rate (1/day)
    "gamma_h": 1 / 5.0,  # human recovery rate (1/day)
    "d_h": 3.5e-5,       # human birth = death rate (1/day)
    "d_v": 1 / 14.0,     # mosquito birth = death rate (1/day)
}


def dengue_rhs(t, y, params):
    # State vector X = (S_h, E_h, I_h, R_h, S_v, E_v, I_v)^T   [Theory 3]
    S_h, E_h, I_h, R_h, S_v, E_v, I_v = y

    # Parameters (same symbols as Theory Derivation)
    p = params
    a0, kappa, eta = p["a0"], p["kappa"], p["eta"]
    b, c, m = p["b"], p["c"], p["m"]
    eps_h, eps_v = p["eps_h"], p["eps_v"]
    gamma_h = p["gamma_h"]
    d_h, d_v = p["d_h"], p["d_v"]
    N_h = p["N_h"]
    N_v = m * N_h                                   # N_v = m N_h   [Theory 2]

    # Effective biting rate   [Theory 1.1]
    a = a0 * (1 - kappa * eta)
    # Governing equations   [Theory 1]
    dS_h = d_h * N_h - a * b * (I_v / N_h) * S_h - d_h * S_h
    dE_h = a * b * (I_v / N_h) * S_h - (eps_h + d_h) * E_h
    dI_h = eps_h * E_h - (gamma_h + d_h) * I_h
    dR_h = gamma_h * I_h - d_h * R_h
    dS_v = d_v * N_v - a * c * (I_h / N_h) * S_v - d_v * S_v
    dE_v = a * c * (I_h / N_h) * S_v - (eps_v + d_v) * E_v
    dI_v = eps_v * E_v - d_v * I_v

    return np.array([dS_h, dE_h, dI_h, dR_h, dS_v, dE_v, dI_v])

def disease_free_equilibrium(params):
    # E_0 = (N_h, 0, 0, 0, m N_h, 0, 0)   [Theory 3]
    N_h = params["N_h"]
    return np.array([N_h, 0, 0, 0, params["m"] * N_h, 0, 0])

def basic_reproduction_number(params):
    p = params
    a = p["a0"] * (1 - p["kappa"] * p["eta"])          # [Theory 1.1]
    b, c, m = p["b"], p["c"], p["m"]
    eps_h, eps_v, gamma_h = p["eps_h"], p["eps_v"], p["gamma_h"]
    d_h, d_v = p["d_h"], p["d_v"]

    # Jacobians at DFE, infected order x = (E_h, I_h, E_v, I_v)   [Theory 4.3]
    F = np.array([[0, 0,         0, a * b],
                  [0, 0,         0, 0    ],
                  [0, a * c * m, 0, 0    ],
                  [0, 0,         0, 0    ]])
    V = np.array([[eps_h + d_h, 0,             0,           0  ],
                  [-eps_h,      gamma_h + d_h, 0,           0  ],
                  [0,           0,             eps_v + d_v, 0  ],
                  [0,           0,             -eps_v,      d_v]])

    K = F @ np.linalg.inv(V)                             # [Theory 5] next generation matrix
    R0_ngm = np.max(np.abs(np.linalg.eigvals(K)))        # [Theory 7] spectral radius

    # Closed form   [Theory 7]
    R0_closed = np.sqrt(a**2 * b * c * m * eps_h * eps_v
                        / ((d_h + eps_h) * (d_h + gamma_h) * (d_v + eps_v) * d_v))
    return R0_ngm, R0_closed

if __name__ == "__main__":
    from scipy.integrate import solve_ivp

    N_h = PARAMS["N_h"]        # human population
    N_v = PARAMS["m"] * N_h    # mosquito population
    #to be fitted by E1
    I_h0 = 10                  # Initial infected humans(hardcoded for now, to be fitted by E1)
    R_h0 = 0                  # Initial recovered humans (hardcoded for now, to be fitted by E1)
    S_h0 = N_h - I_h0 - R_h0
    y0 = [S_h0, 0, I_h0, R_h0, N_v, 0, 0]
    # for testing, print the initial derivative
    dydt = dengue_rhs(0, y0, PARAMS)
    print(dydt)
    n_weeks = 52
    t_end = 7 * n_weeks
    sol = solve_ivp(dengue_rhs, (0, t_end), y0, args=(PARAMS,),
                    t_eval=np.arange(0, t_end + 1), rtol=1e-8)
    last_day = sol.y[:, -1]       # state on the last day     # 
    last_day = last_day.round(1)   
    result = dict(zip(STATE_NAMES, last_day))
    print("final state:", result)

        # ---- checks: success, non-negativity, population balance ----
    print("solver success:", sol.success)
    y = sol.y
    print("all compartments >= 0:", (y >= 0).all())
    err_h = np.abs(y[:4].sum(axis=0) - N_h).max() / N_h   # human: S+E+I+R
    err_v = np.abs(y[4:].sum(axis=0) - N_v).max() / N_v   # mosquito: S+E+I
    print(f"max relative error  N_h: {err_h:.1e},  N_v: {err_v:.1e}")
        # ---- prototype weekly new cases ----
    phi = 0.06                                   # PROVISIONAL reporting fraction (phi; rho is reserved for spectral radius rho(K))
    daily_new = PARAMS["eps_h"] * y[1]           # y[1] = E_h, new symptomatic cases per day
    weekly_new = daily_new[:t_end].reshape(n_weeks, 7).sum(axis=1)   # n_weeks weeks, 7 days each
    weekly_reported = phi * weekly_new
    print("peak weekly reported cases:", weekly_reported.max().round(1),
          "in week", weekly_reported.argmax() + 1)
        # ---- theory checks ----
    E0 = disease_free_equilibrium(PARAMS)
    print("rhs at DFE (should be all 0):", dengue_rhs(0, E0, PARAMS))
    print("R0 (NGM numeric, closed form):", basic_reproduction_number(PARAMS))