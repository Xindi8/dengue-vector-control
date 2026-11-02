
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

# Import shared model components from Xindi's core file
try:
    from dengue_model import dengue_rhs, PARAMS, STATE_NAMES
except ImportError:
    from dengue_ode import dengue_rhs, PARAMS, STATE_NAMES



def simulate_weekly_cases(theta, base_params, n_weeks=52):
    """
    Runs the 7-ODE simulation and aggregates into weekly reported cases:
      weekly_reported = phi * (7-day sum of eps_h * E_h)
    """
    s_frac, phi, m = theta

    p = base_params.copy()
    p["m"] = m
    N_h = p["N_h"]
    N_v = m * N_h

    # Initial state seeding
    I_h0 = 10.0
    R_h0 = 0.0
    E_h0 = 0.0
    S_h0 = s_frac * N_h - I_h0
    y0 = [S_h0, E_h0, I_h0, R_h0, N_v, 0.0, 10.0]

    t_end = 7 * n_weeks
    sol = solve_ivp(
        dengue_rhs,
        (0, t_end),
        y0,
        args=(p,),
        t_eval=np.arange(0, t_end + 1),
        method="RK45",
        rtol=1e-7,
        atol=1e-9,
    )

    if not sol.success:
        return None

    # Observation mapping agreed with D and E2:
    daily_new = p["eps_h"] * sol.y[1]  # eps_h * E_h
    weekly_new = daily_new[:t_end].reshape(n_weeks, 7).sum(axis=1)
    weekly_reported = phi * weekly_new
    return weekly_reported


def generate_synthetic_truth(base_params, true_theta, n_weeks=52, seed=42):
    """
    Generates a synthetic outbreak curve with observation noise to verify
    the optimizer's capability to invert the true parameters.
    """
    np.random.seed(seed)
    clean_curve = simulate_weekly_cases(true_theta, base_params, n_weeks)
    # Poisson-like observation noise
    noise = np.random.normal(0, np.sqrt(np.maximum(1.0, clean_curve)))
    noisy_cases = np.maximum(0.0, clean_curve + noise)
    return clean_curve, noisy_cases


def run_calibration_scaffold(observed_data, base_params, bounds, starts):
    """
    Multi-start Trust Region Reflective (TRF) nonlinear least-squares.
    """
    n_weeks = len(observed_data)

    def residuals(theta):
        pred = simulate_weekly_cases(theta, base_params, n_weeks)
        if pred is None or np.any(np.isnan(pred)):
            return np.ones(n_weeks) * 1e5
        return pred - observed_data

    best_res = None
    best_cost = np.inf

   
    for i, x0 in enumerate(starts):
        try:
            res = least_squares(
                residuals,
                x0=x0,
                bounds=bounds,
                method="trf",
                ftol=1e-6,
                xtol=1e-6,
                max_nfev=350,
            )
            print(f"Pass {i+1}: Start={x0} -> Cost={res.cost:.2f}, Success={res.success}")
            if res.success and res.cost < best_cost:
                best_cost = res.cost
                best_res = res
        except Exception as err:
            print(f"Pass {i+1} failed with error: {err}")

    return best_res


if __name__ == "__main__":
    # Benchmark synthetic truth: [s_frac, phi, m]
    # m=5.5 produces an active outbreak season (R0 ~ 1.7)
    true_theta = [0.85, 0.06, 5.5]
    n_weeks = 52

    true_curve, noisy_cases = generate_synthetic_truth(PARAMS, true_theta, n_weeks)

    # Literature-based bounds:
    # s_frac in [0.4, 1.0], phi in [0.035, 0.090] (CDC PR range), m in [1.5, 10.0]
    param_bounds = (
        [0.40, 0.035, 1.5],
        [1.00, 0.090, 10.0],
    )

    # Multi-start initial guesses
    initial_guesses = [
        [0.70, 0.050, 4.0],
        [0.90, 0.075, 7.0],
        [0.55, 0.040, 3.0],
    ]

    fit = run_calibration_scaffold(noisy_cases, PARAMS, param_bounds, initial_guesses)

    print("\n---------------- CALIBRATION SUMMARY ----------------")
    print(f"True Benchmark    [s_frac, phi, m]: {true_theta}")
    print(f"Fitted Estimates  [s_frac, phi, m]: {np.round(fit.x, 4)}")

    fitted_curve = simulate_weekly_cases(fit.x, PARAMS, n_weeks)
    rss = np.sum((fitted_curve - noisy_cases) ** 2)
    print(f"Residual Sum of Squares (RSS): {rss:.2f}")
    print("-----------------------------------------------------")

    # Save verification plot
    plt.figure(figsize=(10, 5))
    weeks = np.arange(1, n_weeks + 1)
    plt.scatter(weeks, noisy_cases, color="black", alpha=0.6, label="Synthetic Observed Cases ($C_k$)")
    plt.plot(weeks, true_curve, "g--", linewidth=1.5, label="True Underlying Curve")
    plt.plot(weeks, fitted_curve, "r-", linewidth=2.0, label="E1 Calibrated Model Fit")
    plt.title("MATH 574 Team 2 — Role E1: Fitting Scaffold Verification (Zaid Uddin)")
    plt.xlabel("Epidemiological Week")
    plt.ylabel("Weekly Reported Cases")
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig("E1_synthetic_fit_verification.png", dpi=200)
    print("Saved plot as 'E1_synthetic_fit_verification.png'.")