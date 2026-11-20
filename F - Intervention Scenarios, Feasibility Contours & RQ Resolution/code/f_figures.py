import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 11, "axes.titlesize": 12, "axes.labelsize": 11,
                     "legend.fontsize": 9, "axes.spines.top": False, "axes.spines.right": False})

curves = pd.read_csv("../outputs/f_scenario_curves.csv", index_col="week")
table = pd.read_csv("../outputs/f_scenario_results.csv").set_index("ID")
observed = pd.read_csv("../inputs/e1_fit_curves.csv")["observed"].values
rq2 = pd.read_csv("../outputs/f_rq2_equal_cuts.csv")

styles = {"S0": ("black", "-", 2.6), "S1": ("#d62728", "-", 2.0), "S2": ("#ff7f0e", "--", 2.0),
          "S3": ("#9467bd", "-", 2.0), "S4": ("#1f77b4", "-.", 2.0), "S5": ("#2ca02c", "-", 2.0),
          "S6": ("#17becf", ":", 2.4)}

fig, ax = plt.subplots(figsize=(6.5, 5.4), constrained_layout=True)
ax.plot(curves.index, observed, "o", color="0.55", ms=4, label="Observed 2010/11 (DENV-1)")
for sid, (color, ls, lw) in styles.items():
    row = table.loc[sid]
    ax.plot(curves.index, curves[sid], color=color, ls=ls, lw=lw,
            label=f"{sid} {row['Name']}: R$_{{eff}}$ {row['R_eff']:.2f}, {row['pct_averted']:.0f}% averted")
ax.set_xlabel("Week of season (from 30 Apr 2010)")
ax.set_ylabel("Expected reported cases per week")
ax.set_title("Weekly reported cases under each package\n(controls start 30 Apr 2010; s$_0$ = 0.5 fit, η = 0.5, q = 0.9)")
ax.set_xlim(1, 52)
ax.set_ylim(0, None)
fig.legend(loc="outside lower center", ncol=1, frameon=False)
fig.savefig("../figures/fig4_trajectories.png", dpi=300)
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(6.5, 3.8), constrained_layout=True, sharey=True)
cuts = sorted(rq2["cut"].unique())
x = np.arange(len(cuts))
w = 0.38
for ax, col, title in zip(axes, ["pct_averted", "peak_drop_pct"],
                          ["Cases averted over 52 weeks", "Reduction in peak weekly cases"]):
    for i, (lever, color) in enumerate([("biting rate a", "#1f77b4"), ("mosquito density m", "#ff7f0e")]):
        vals = rq2[rq2["lever"] == lever].set_index("cut").loc[cuts, col].values
        bars = ax.bar(x + (i - 0.5) * w, vals, w, color=color, label=f"Cut {lever}")
        ax.bar_label(bars, fmt="%.0f%%", fontsize=8, padding=2)
    ax.set_xticks(x, [f"{int(c * 100)}% cut" for c in cuts])
    ax.set_title(title)
    ax.set_ylim(0, 100)
axes[0].set_ylabel("Percent of baseline (S0)")
fig.suptitle("Same percentage cut, biting rate vs mosquito density (s$_0$ = 0.5 fit)")
fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower center", ncol=2, frameon=False)
fig.savefig("../figures/fig3_rq2_equal_cuts.png", dpi=300)
plt.close(fig)
