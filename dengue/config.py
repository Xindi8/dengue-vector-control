"""Every setting shared across roles, in one place. Change a value here and rerun
``python run_all.py``; every table and figure follows.

Two choices differ between the roles' own documents. Both are set here once:

* Headline s0. AB recommends 0.4 (evidence puts s0 at about 0.29-0.44), E1 fitted
  0.2 / 0.4 / 0.6, and D's README lists 0.4. F's documents used 0.5 with a 0.39-0.64
  band (the range E2 derived under Table 1's m <= 10 bound, which AB later raised to 20).
  The project uses AB's values; set S0_HEADLINE = 0.5 and S0_BAND = (0.39, 0.64) to
  reproduce F's documents exactly.
* Balanced package S5. F set delta = 0.65 by hand for s0 = 0.5. Here delta is found by
  F's rule (keep kappa at 0.35 and raise delta in 0.05 steps until R0 < 1 at the
  headline s0), which gives 0.65 at s0 = 0.5.
"""

# --- s0: share of the population susceptible to DENV-1 on 30 Apr 2010 (AB, E1, F) ---
S0_HEADLINE = 0.4
S0_BAND = (0.2, 0.6)
# Fitted once by E1. Covers every value any role reports, so all tables share one set of fits.
S0_FIT_GRID = (0.1, 0.2, 0.3, 0.39, 0.4, 0.5, 0.6, 0.64, 0.8, 1.0)

# --- calibration (E1, following E2's plan) ---
FIT_M_BOUNDS = (0.05, 2000.0)
FIT_PHI_BOUNDS = (1e-5, 1.0)
FIT_DIFF_STEP = 1e-4          # least_squares' default step hides the s0-m-phi ridge (E2, check 07)
FIT_STARTS_M_TIMES_S0 = (4.0, 2.0, 8.0)   # multi-start along the m*s0 ridge (E1)

# --- interventions (F) ---
ETA = 0.5                     # bites a protected person avoids (AB range 0.3-0.9)
Q_WOLBACHIA = 0.9             # transmission a Wolbachia carrier blocks (theoretical)
KAPPA_MAX = 0.6               # proposal's ceiling for compliance
DELTA_MAX = 0.85              # proposal's ceiling for source reduction
PROPOSAL_POINT = (0.35, 0.40)
P_W_LEVELS = (0.0, 0.4, 0.8)  # Wolbachia coverage for the contour panels
RQ2_CUTS = (0.1, 0.2, 0.3)
BALANCED_KAPPA = 0.35
BALANCED_STEP = 0.05

# (ID, name, kappa, delta, Wolbachia coverage). delta None = chosen by F's balanced rule.
PACKAGES = [
    ("S0", "Baseline", 0.0, 0.0, 0.0),
    ("S1", "Realistic protection", 0.35, 0.0, 0.0),
    ("S2", "Moderate spraying", 0.0, 0.4, 0.0),
    ("S3", "Proposal's combination", 0.35, 0.4, 0.0),
    ("S4", "Max feasible, no Wolbachia", 0.6, 0.85, 0.0),
    ("S5", "Balanced package", BALANCED_KAPPA, None, 0.0),
    ("S6", "Everything at its ceiling", 0.6, 0.85, 0.8),
]
