import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R0_BASE = 2.0
ETA = 0.7
PHI = 0.9
P_LEVELS = [0.0, 0.5, 0.8]
KAPPA_MAX = 0.6
DELTA_MAX = 0.85
CONVENTION = "step"


def r0_controlled(kappa, delta, p, r0=R0_BASE, eta=ETA, phi=PHI,
                  convention=CONVENTION):
    step_factor = (1 - kappa * eta) * np.sqrt(1 - delta) * np.sqrt(1 - p * phi)
    if convention == "step":
        return r0 * step_factor
    if convention == "cycle":
        return r0 * step_factor**2
    raise ValueError("convention must be 'step' or 'cycle'")


def main():
    kappa = np.linspace(0, 1, 401)
    delta = np.linspace(0, 0.99, 400)
    K, D = np.meshgrid(kappa, delta)

    plt.rcParams.update({"font.size": 12, "axes.titlesize": 13,
                         "axes.labelsize": 12, "legend.fontsize": 10})

    handles = [Line2D([], [], color="black", lw=2.2, label="R$_0$ = 1 boundary"),
               Line2D([], [], color="navy", lw=2, ls="--",
                      label=f"Feasible limits (κ ≤ {KAPPA_MAX}, δ ≤ {DELTA_MAX})"),
               Line2D([], [], marker="o", color="white", mec="black", ms=8, ls="",
                      label="Proposal: κ = 0.35, δ = 0.4")]

    for p in P_LEVELS:
        fig, ax = plt.subplots(figsize=(6.5, 5.6), constrained_layout=True)
        R = r0_controlled(K, D, p)

        cf = ax.contourf(K, D, R, levels=np.linspace(0, R0_BASE, 21),
                         cmap="RdYlGn_r", extend="max")
        ax.contour(K, D, R, levels=[1.0], colors="black", linewidths=2.2)
        ax.plot([0, KAPPA_MAX, KAPPA_MAX], [DELTA_MAX, DELTA_MAX, 0],
                color="navy", lw=2, ls="--")
        ax.plot(0.35, 0.40, marker="o", color="white", mec="black", ms=8)

        ax.set_title(f"R$_0$ = 1 threshold, Wolbachia coverage p = {p:.1f}\n"
                     f"(baseline R$_0$ = {R0_BASE}, η = {ETA}, φ = {PHI})")
        ax.set_xlabel("Compliance κ")
        ax.set_ylabel("Adult mosquito reduction δ")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 0.99)
        fig.legend(handles=handles, loc="outside lower center", ncol=2, frameon=False)
        fig.colorbar(cf, ax=ax, label="Controlled R$_0$")

        fig.savefig(f"r0_contour_p{int(round(p * 10)):02d}.png", dpi=300)
        plt.close(fig)


if __name__ == "__main__":
    main()
