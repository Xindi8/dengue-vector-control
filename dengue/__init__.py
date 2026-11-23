"""Dengue vector control in San Juan, Puerto Rico (MATH 574, Team 2).

One package for the whole project:

    model          D   shared 7-state SEIR_h-SEI_v model and R0 (theory from C)
    data           AB  loaders for the data and literature package
    observation    E2  weekly expected cases, starting state, residuals, holdout tools
    diagnostics    E2  residual diagnosis for fits (RSS, AIC, correlations)
    calibration    E1  fits of m and phi for each fixed s0
    interventions  F   control scenarios, thresholds and equal-cut comparison
    figures            the project's figures
    config             every setting shared across roles, in one place

Run the whole pipeline with ``python run_all.py`` from the project folder.
"""

__version__ = "1.0.0"
