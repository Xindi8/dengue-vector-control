"""Check 3 - synthetic data, where we know the truth.

(a) Start-up: starting the epidemic with no infected mosquitoes (the old y0)
    creates a fake dip in the first weeks; starting on the growth mode does not.
(b) Recovery: fake NB data from known parameters. With s0 fixed at the truth,
    m and phi come back. With s0 fixed at a wrong value, the fit is just as good
    and m*s0, phi*s0 still come back: even a perfect model cannot find s0.
"""
import numpy as np
import pandas as pd

from common import OUT, baseline_params, dengue_rhs, fit_m_phi, load, obs, r0

cal, edges = load()
TRUE = dict(s0=0.40, m=10.0, phi=0.005, k=290.0)
p_true = dict(baseline_params(), m=TRUE["m"])
y_real = cal.denv1_attributed_cases.to_numpy(float)

# ---------------- (a) start-up transient ----------------
rate = obs.anchored_onset_rate(dengue_rhs, p_true, TRUE["s0"], TRUE["phi"], edges, y_real)
N = p_true["N_h"]


def naive_start(E0):
    """The old start: E_h = I_h = E0, no infected mosquitoes."""
    return np.array([TRUE["s0"] * N, E0, E0, N - TRUE["s0"] * N - 2 * E0, p_true["m"] * N, 0.0, 0.0])


# anchor the naive start to the same first three weeks, so the comparison is fair
E0 = rate / p_true["eps_h"]
for _ in range(6):
    first3 = TRUE["phi"] * obs.simulate_weekly_onsets(dengue_rhs, p_true, naive_start(E0), edges[:4]).sum()
    E0 *= y_real[:3].sum() / first3
mu_naive = TRUE["phi"] * obs.simulate_weekly_onsets(dengue_rhs, p_true, naive_start(E0), edges)
mu_mode = TRUE["phi"] * obs.simulate_weekly_onsets(dengue_rhs, p_true, obs.initial_state(p_true, TRUE["s0"], rate), edges)
start = pd.DataFrame({"week": np.arange(1, 9), "data": y_real[:8].round(1),
                      "naive_start": mu_naive[:8].round(1), "growth_mode_start": mu_mode[:8].round(1)})
print("(a) First 8 weeks, expected cases:")
print(start.to_string(index=False))
start.to_csv(f"{OUT}/check03_startup.csv", index=False)

# ---------------- (b) recovery on synthetic data ----------------
mu_true = obs.expected_cases(p_true, TRUE["s0"], TRUE["phi"], edges, rhs=dengue_rhs, y_anchor=y_real)
rng = np.random.default_rng(574)
rows = []
for rep in range(20):
    k = TRUE["k"]
    y_syn = rng.negative_binomial(k, k / (k + mu_true)).astype(float)
    for s0_fixed in (0.40, 0.20, 0.80):
        f = fit_m_phi(y_syn, edges, s0_fixed, m0=4.0 / s0_fixed, phi0=0.002 / s0_fixed)
        rows.append(dict(rep=rep, s0_fixed=s0_fixed, m=f["m"], phi=f["phi"], cost=f["cost"],
                         m_s0=f["m"] * s0_fixed, phi_s0=f["phi"] * s0_fixed, R0=f["R0"], Reff=f["Reff"]))
res = pd.DataFrame(rows)
summ = res.groupby("s0_fixed").agg(
    m_median=("m", "median"), phi_median=("phi", "median"), m_s0=("m_s0", "median"),
    phi_s0=("phi_s0", "median"), Reff=("Reff", "median"), R0=("R0", "median"),
    cost_median=("cost", "median"))
print(f"\n(b) Truth: s0={TRUE['s0']}, m={TRUE['m']}, phi={TRUE['phi']}, "
      f"m*s0={TRUE['m'] * TRUE['s0']:.2f}, phi*s0={TRUE['phi'] * TRUE['s0']:.4f}, "
      f"R0={r0(p_true):.2f}, Reff={r0(p_true) * np.sqrt(TRUE['s0']):.3f}")
print(summ.round(4).to_string())
# paired cost difference: wrong s0 vs right s0, same synthetic data set
piv = res.pivot(index="rep", columns="s0_fixed", values="cost")
print("\nCost difference vs the true s0 (same data, 20 replicates):")
for s in (0.20, 0.80):
    d = piv[s] - piv[0.40]
    print(f"  s0 fixed at {s}: median {d.median():+.4f}, range {d.min():+.4f} to {d.max():+.4f}")
res.to_csv(f"{OUT}/check03_synthetic_recovery.csv", index=False)
