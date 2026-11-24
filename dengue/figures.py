"""The project's figures. Contours, trajectories and equal-cut bars follow F's designs
(f_contour_plots.py, f_figures.py); the calibration figure follows E1's fit plot."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from . import config  # noqa: E402
from .interventions import controlled_r  # noqa: E402

INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e0"


def _style():
    # DejaVu Sans ships with matplotlib, so the figures come out the same on every computer
    plt.rcParams.update({"font.family": "DejaVu Sans", "text.hinting": "default",
                         "font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11, "legend.fontsize": 9,
                         "axes.spines.top": False, "axes.spines.right": False})


def calibration_fits(curves, fits, path, headline=config.S0_HEADLINE, band=config.S0_BAND):
    _style()
    fig, ax = plt.subplots(figsize=(7.5, 4.6), constrained_layout=True)
    ax.plot(curves.week, curves.observed, "o", color="#a3a29c", ms=4, label="Observed, DENV-1 attributed")
    shown = sorted({band[0], headline, band[1]})
    shades = {s: c for s, c in zip(shown, ["#9cc3ef", "#2a78d6", "#14467f"])}
    for s0 in shown:
        row = fits[np.isclose(fits.s0, s0)].iloc[0]
        ax.plot(curves.week, curves[f"fit_s0_{s0:g}"], color=shades[s0], lw=2.4 if s0 == headline else 1.4,
                label=f"Fit, s0 = {s0:g}: m = {row.m:.1f}, φ = {row.phi:.4f}, R0 = {row.R0:.2f}")
    r_eff = float(fits.R_eff.mean())
    ax.set_title(f"San Juan 2010/11: every s0 fits the same curve (R_eff = {r_eff:.2f})", loc="left")
    ax.set_xlabel("Week of season (week 1 starts 30 Apr 2010)")
    ax.set_ylabel("Reported cases per week")
    ax.grid(axis="y", color=GRID)
    ax.legend(frameon=False)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def contours(r_eff, r0_head, r0_band, out_paths, headline=config.S0_HEADLINE, band=config.S0_BAND):
    """R = 1 thresholds over compliance kappa and source reduction delta, one panel per
    Wolbachia coverage level."""
    _style()
    kappa = np.linspace(0, 1, 401)
    delta = np.linspace(0, 0.99, 400)
    K, D = np.meshgrid(kappa, delta)

    def boundary(ax, R, **style):
        if R.min() < 1 < R.max():
            ax.contour(K, D, R, levels=[1.0], **style)
            return ""
        return ": already below 1" if R.max() <= 1 else ": above 1 on the whole grid"

    def delta_needed(R, p_w):
        """delta at which R (1 - kappa eta) sqrt(1 - delta) sqrt(1 - p_w q) = 1, for each kappa."""
        r_k = controlled_r(R, kappa, 0.0, p_w)
        return np.clip(1 - 1 / r_k ** 2, 0.0, 0.99)

    for p_w, path in zip(config.P_W_LEVELS, out_paths):
        fig, ax = plt.subplots(figsize=(6.5, 6.2), constrained_layout=True)
        ax.fill_between([0, config.KAPPA_MAX], 0, config.DELTA_MAX, color="0.93", zorder=0)
        ax.plot([0, config.KAPPA_MAX, config.KAPPA_MAX], [config.DELTA_MAX, config.DELTA_MAX, 0],
                color="0.35", lw=1.5, ls="--", zorder=2)
        # band between the R0 = 1 curves at the two ends of the s0 band (lower s0 -> higher R0)
        lower, upper = delta_needed(min(r0_band), p_w), delta_needed(max(r0_band), p_w)
        in_band = upper > lower
        if in_band.any():
            ax.fill_between(kappa, lower, upper, where=in_band, color="#c9b8e8", alpha=0.8, lw=0, zorder=1)
        handles = []
        note_r0 = boundary(ax, controlled_r(r0_head, K, D, p_w), colors="#6a3d9a", linewidths=2.6, zorder=3)
        note_reff = boundary(ax, controlled_r(r_eff, K, D, p_w), colors="black", linewidths=2.6, zorder=3)
        handles += [Line2D([], [], color="black", lw=2.6, label=f"R_eff = 1 (data, any s0){note_reff}"),
                    Line2D([], [], color="#6a3d9a", lw=2.6, label=f"R0 = 1, s0 = {headline:g} (headline){note_r0}"),
                    Patch(facecolor="#c9b8e8", label=f"R0 = 1 for s0 between {band[0]:g} and {band[1]:g}"
                          + ("" if in_band.any() else ": none on grid")),
                    Line2D([], [], color="0.35", ls="--", lw=1.5,
                           label=f"Feasible region (κ ≤ {config.KAPPA_MAX}, δ ≤ {config.DELTA_MAX})"),
                    Line2D([], [], marker="o", color="white", mec="black", ms=9, ls="",
                           label=f"Proposal: κ = {config.PROPOSAL_POINT[0]}, δ = {config.PROPOSAL_POINT[1]}")]
        ax.plot(*config.PROPOSAL_POINT, marker="o", color="white", mec="black", ms=9, zorder=4)
        ax.set_title(f"Threshold R = 1, Wolbachia coverage p_w = {p_w:.1f}\n"
                     f"(η = {config.ETA}, q = {config.Q_WOLBACHIA}; R0 and R_eff from the 2010/11 fits)")
        ax.set_xlabel("Compliance κ")
        ax.set_ylabel("Adult mosquito reduction δ")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 0.99)
        fig.legend(handles=handles, loc="outside lower center", ncol=1, frameon=False)
        fig.savefig(path, dpi=200)
        plt.close(fig)


STYLES = {"S0": ("black", "-", 2.6), "S1": ("#d62728", "-", 2.0), "S2": ("#ff7f0e", "--", 2.0),
          "S3": ("#9467bd", "-", 2.0), "S4": ("#1f77b4", "-.", 2.0), "S5": ("#2ca02c", "-", 2.0),
          "S6": ("#17becf", ":", 2.4)}


def trajectories(curves, table, observed, path, headline=config.S0_HEADLINE):
    _style()
    fig, ax = plt.subplots(figsize=(6.5, 5.6), constrained_layout=True)
    ax.plot(curves.index, observed, "o", color="0.55", ms=4, label="Observed 2010/11 (DENV-1)")
    t = table.set_index("ID")
    for sid, (color, ls, lw) in STYLES.items():
        row = t.loc[sid]
        ax.plot(curves.index, curves[sid], color=color, ls=ls, lw=lw,
                label=f"{sid} {row['Name']} (κ {row['kappa']:g}, δ {row['delta']:g}, p_w {row['p_w']:g}): "
                      f"R_eff {row['R_eff']:.2f}, {row['pct_averted']:.0f}% averted")
    ax.set_xlabel("Week of season (from 30 Apr 2010)")
    ax.set_ylabel("Expected reported cases per week")
    ax.set_title(f"Weekly reported cases under each package\n(controls start 30 Apr 2010; s0 = {headline:g} fit, "
                 f"η = {config.ETA}, q = {config.Q_WOLBACHIA})")
    ax.set_xlim(1, len(curves))
    ax.set_ylim(0, None)
    fig.legend(loc="outside lower center", ncol=1, frameon=False)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def equal_cuts(rq2, path, headline=config.S0_HEADLINE):
    _style()
    fig, axes = plt.subplots(1, 2, figsize=(6.8, 3.9), constrained_layout=True, sharey=True)
    cuts = sorted(rq2["cut"].unique())
    x = np.arange(len(cuts))
    w = 0.38
    for ax, col, title in zip(axes, ["pct_averted", "peak_drop_pct"],
                              ["Cases averted over 52 weeks", "Reduction in peak weekly cases"]):
        for i, (lever, color) in enumerate([("biting rate a", "#1f77b4"), ("mosquito density m", "#ff7f0e")]):
            vals = rq2[rq2["lever"] == lever].set_index("cut").loc[cuts, col].values
            bars = ax.bar(x + (i - 0.5) * w, vals, w, color=color, label=f"Cut {lever}")
            ax.bar_label(bars, fmt="%.0f%%", fontsize=8, padding=2)
        ax.set_xticks(x, [f"{int(round(c * 100))}% cut" for c in cuts])
        ax.set_title(title)
        ax.set_ylim(0, 100)
    axes[0].set_ylabel("Percent of baseline (S0)")
    fig.suptitle(f"Same percentage cut: biting rate vs mosquito density (s0 = {headline:g} fit)")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower center", ncol=2, frameon=False)
    fig.savefig(path, dpi=200)
    plt.close(fig)
