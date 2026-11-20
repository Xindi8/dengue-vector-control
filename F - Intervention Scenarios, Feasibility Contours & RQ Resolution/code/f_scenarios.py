import numpy as np
import pandas as pd
from f_model import anchored_start, apply_controls, weekly_cases, r0, WEEK_EDGES

OBSERVED = pd.read_csv("../inputs/e1_fit_curves.csv")["observed"].values
ETA, Q = 0.5, 0.9
S0_LEVELS = (0.39, 0.5, 0.64)
REF_S0 = 0.5

PACKAGES = [
    ("S0", "Baseline", 0.0, 0.0, 0.0),
    ("S1", "Realistic protection", 0.35, 0.0, 0.0),
    ("S2", "Moderate spraying", 0.0, 0.4, 0.0),
    ("S3", "Proposal's combination", 0.35, 0.4, 0.0),
    ("S4", "Max feasible, no Wolbachia", 0.6, 0.85, 0.0),
    ("S5", "Balanced package", 0.35, 0.65, 0.0),
    ("S6", "Everything at its ceiling", 0.6, 0.85, 0.8),
]

EDGES_2Y = np.concatenate([WEEK_EDGES, 365 + 7 * np.arange(1, 53)])


def run(s0, kappa, delta, p_w, eta=ETA, q=Q, edges=WEEK_EDGES):
    p, y0, phi = anchored_start(s0, OBSERVED)
    pc, yc = apply_controls(p, y0, kappa, delta, p_w, eta, q)
    mu = weekly_cases(pc, yc, phi, edges)
    return pc, mu


def scenario_table():
    rows, curves = [], {}
    base = {s0: run(s0, 0, 0, 0)[1] for s0 in S0_LEVELS}
    for sid, name, k, d, w in PACKAGES:
        row = dict(ID=sid, Name=name, kappa=k, delta=d, p_w=w)
        for s0 in S0_LEVELS:
            pc, mu = run(s0, k, d, w)
            row[f"R0_s0={s0}"] = r0(pc)
            if s0 == REF_S0:
                row["R_eff"] = r0(pc) * np.sqrt(s0)
                row["cases_52wk"] = mu.sum()
                row["pct_averted"] = 100 * (1 - mu.sum() / base[s0].sum())
                row["peak"] = mu.max()
                row["peak_week"] = int(mu.argmax()) + 1
                curves[sid] = mu
                _, mu2 = run(s0, k, d, w, edges=EDGES_2Y)
                row["cases_wk53_104"] = mu2[52:].sum()
            row[f"averted_s0={s0}"] = 100 * (1 - mu.sum() / base[s0].sum())
        row["feasible"] = k <= 0.6 and d <= 0.85
        rows.append(row)
    return pd.DataFrame(rows), curves


def rq2_table(cuts=(0.1, 0.2, 0.3)):
    _, base = run(REF_S0, 0, 0, 0)
    rows = []
    for x in cuts:
        for lever, args in (("biting rate a", (x, 0.0)), ("mosquito density m", (0.0, x))):
            pc, mu = run(REF_S0, args[0], args[1], 0.0, eta=1.0)
            rows.append(dict(cut=x, lever=lever, R_eff=r0(pc) * np.sqrt(REF_S0),
                             pct_averted=100 * (1 - mu.sum() / base.sum()),
                             peak_drop_pct=100 * (1 - mu.max() / base.max()),
                             peak_week=int(mu.argmax()) + 1))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    table, curves = scenario_table()
    table.to_csv("../outputs/f_scenario_results.csv", index=False)
    pd.DataFrame(curves, index=np.arange(1, 53)).rename_axis("week").to_csv("../outputs/f_scenario_curves.csv")
    rq2 = rq2_table()
    rq2.to_csv("../outputs/f_rq2_equal_cuts.csv", index=False)
    pd.set_option("display.width", 250)
    print(table.round(3).to_string(index=False))
    print(rq2.round(2).to_string(index=False))
