import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

FIT_RESULTS = "../inputs/e1_fit_results.csv"

ETA = 0.5
Q_BLOCK = 0.9
P_W_LEVELS = [0.0, 0.4, 0.8]
KAPPA_MAX = 0.6
DELTA_MAX = 0.85
PROPOSAL = (0.35, 0.40)
REF_S0 = 0.4
BAND_S0 = (0.29, 0.44)
SENS_S0 = (0.2, 0.6)


def load_fit(path=FIT_RESULTS):
    data = np.genfromtxt(path, delimiter=",", names=True)
    r_eff = float(np.mean(data["R_eff"]))
    r0 = {round(float(row["s0"]), 2): float(row["R0"]) for row in data}
    return r_eff, r0


def r0_at(s0, r0, r_eff):
    return r0.get(round(s0, 2), r_eff / np.sqrt(s0))


def r_controlled(r_base, kappa, delta, p_w, eta=ETA, q=Q_BLOCK):
    return r_base * (1 - kappa * eta) * np.sqrt(1 - delta) * np.sqrt(1 - p_w * q)


def boundary(ax, R, **style):
    if R.min() < 1 < R.max():
        ax.contour(K, D, R, levels=[1.0], **style)
        return ""
    return ": already below 1" if R.max() <= 1 else ": above 1 on whole grid"


r_eff, r0 = load_fit()
kappa = np.linspace(0, 1, 401)
delta = np.linspace(0, 0.99, 400)
K, D = np.meshgrid(kappa, delta)

plt.rcParams.update({"font.size": 12, "axes.titlesize": 13,
                     "axes.labelsize": 12, "legend.fontsize": 10})

for i, p_w in enumerate(P_W_LEVELS):
    fig, ax = plt.subplots(figsize=(6.5, 6.2), constrained_layout=True)
    ax.fill_between([0, KAPPA_MAX], 0, DELTA_MAX, color="0.93", zorder=0)
    ax.plot([0, KAPPA_MAX, KAPPA_MAX], [DELTA_MAX, DELTA_MAX, 0],
            color="0.35", lw=1.5, ls="--", zorder=1)

    lo, hi = (r_controlled(r0_at(s, r0, r_eff), K, D, p_w) for s in BAND_S0)
    band = (lo >= 1) & (hi < 1)
    if band.any():
        ax.contourf(K, D, band.astype(float), levels=[0.5, 1.5], colors=["#c9b8e8"], alpha=0.8, zorder=1)

    handles = []
    for s0 in SENS_S0:
        note = boundary(ax, r_controlled(r0_at(s0, r0, r_eff), K, D, p_w),
                        colors="0.55", linestyles=":", linewidths=1.6, zorder=2)
        handles.append(Line2D([], [], color="0.55", ls=":", lw=1.6,
                              label=f"R$_0$, s$_0$ = {s0} (sensitivity){note}"))

    note = boundary(ax, r_controlled(r0_at(REF_S0, r0, r_eff), K, D, p_w),
                    colors="#6a3d9a", linestyles="-", linewidths=2.6, zorder=3)
    handles.insert(0, Line2D([], [], color="#6a3d9a", lw=2.6,
                             label=f"R$_0$, s$_0$ = {REF_S0} (E1 reference){note}"))
    handles.insert(1, Patch(facecolor="#c9b8e8", label=f"R$_0$ = 1 range for s$_0$ = {BAND_S0[0]}–{BAND_S0[1]} (evidence)"
                            + ("" if band.any() else ": none on grid")))

    note = boundary(ax, r_controlled(r_eff, K, D, p_w),
                    colors="black", linestyles="-", linewidths=2.6, zorder=3)
    handles.insert(0, Line2D([], [], color="black", lw=2.6, label=f"R$_{{eff}}$ (data, any s$_0$){note}"))

    ax.plot(*PROPOSAL, marker="o", color="white", mec="black", ms=9, zorder=4)
    handles += [Line2D([], [], color="0.35", ls="--", lw=1.5,
                       label=f"Feasible region (κ ≤ {KAPPA_MAX}, δ ≤ {DELTA_MAX})"),
                Line2D([], [], marker="o", color="white", mec="black", ms=9, ls="",
                       label=f"Proposal: κ = {PROPOSAL[0]}, δ = {PROPOSAL[1]}")]

    ax.set_title(f"Threshold R = 1, Wolbachia coverage p$_w$ = {p_w:.1f}\n"
                 f"(η = {ETA}, q = {Q_BLOCK}; R$_0$ and R$_{{eff}}$ from E1's 2010/11 fit)")
    ax.set_xlabel("Compliance κ")
    ax.set_ylabel("Adult mosquito reduction δ")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.99)
    fig.legend(handles=handles, loc="outside lower center", ncol=1, frameon=False)

    fig.savefig(f"../figures/fig2{'abc'[i]}_contours_pw{int(round(p_w * 10)):02d}.png", dpi=300)
    plt.close(fig)
