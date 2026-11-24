"""Run the whole project: data checks, model checks, calibration, validation, interventions,
figures and a summary of the key numbers.

    python run_all.py                    # everything (a few minutes)
    python run_all.py --skip-validation  # leave out E2's validation checks
    python run_all.py --rebuild-data     # also rebuild the AB data package from its raw files

Outputs go to results/ (tables, validation logs, SUMMARY.md) and figures/.
"""
from __future__ import annotations

import argparse
import hashlib
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd

from dengue import calibration, config, data, figures, interventions, observation as obs
from dengue.model import PARAMS, basic_reproduction_number, dengue_rhs, r0
from dengue.paths import AB, FIGURES, RESULTS, ROOT


def stage(title):
    print(f"\n=== {title} ===")
    return time.time()


def check_data(rebuild: bool):
    if rebuild:
        subprocess.run([sys.executable, str(AB / "build_ab_package.py")], check=True, cwd=AB)
    manifest = pd.read_csv(AB / "checks" / "raw_manifest.csv")
    bad = [r.file for r in manifest.itertuples()
           if hashlib.sha256((AB / r.file).read_bytes()).hexdigest() != r.sha256]
    if bad:
        raise SystemExit(f"Raw files differ from the manifest: {bad}")
    cal = data.calibration()
    print(f"  {len(manifest)} raw files match their checksums; calibration season has {len(cal)} weeks, "
          f"{cal.denv1_attributed_cases.sum():,.0f} DENV-1-attributed cases")
    return cal


def check_model():
    problems = obs.validate_rhs(dengue_rhs, PARAMS)
    ngm, closed = basic_reproduction_number(PARAMS)
    if problems or not np.isclose(ngm, closed, rtol=1e-10):
        raise SystemExit(f"Model checks failed: {problems}, R0 {ngm} vs {closed}")
    print(f"  E2's contract test passes; next-generation and closed-form R0 agree ({closed:.4f} at m = {PARAMS['m']:g})")


def run_validation():
    out = RESULTS / "validation"
    out.mkdir(parents=True, exist_ok=True)
    for script in sorted((ROOT / "validation").glob("check_0*.py")):
        log = out / f"{script.stem}.log"
        with open(log, "w", encoding="utf-8") as fh:   # UTF-8 on every OS, Windows included
            done = subprocess.run([sys.executable, str(script)], cwd=script.parent, stdout=fh,
                                  stderr=subprocess.STDOUT, env={**os.environ, "PYTHONIOENCODING": "utf-8"})
        print(f"  {script.name}: {'ok' if done.returncode == 0 else 'FAILED, see ' + str(log)}")
        if done.returncode:
            raise SystemExit(1)


def summary(cal, fits, table, rq2, thresholds, path):
    head, band = config.S0_HEADLINE, config.S0_BAND
    f = fits.set_index(fits.s0.round(2))
    t = table.set_index("ID")
    r_eff = float(fits.R_eff.mean())
    lines = [
        "# Key results",
        "",
        "Written by `run_all.py`. Change settings in `dengue/config.py` and rerun.",
        "",
        f"- Headline s0 = {head:g}; sensitivity band {band[0]:g}–{band[1]:g}; η = {config.ETA}, "
        f"Wolbachia blocking q = {config.Q_WOLBACHIA}.",
        f"- Calibration season: {len(cal)} weeks from 30 Apr 2010, {cal.denv1_attributed_cases.sum():,.0f} "
        "DENV-1-attributed cases.",
        "",
        "## Calibration (E1)",
        "",
        f"- R_eff on 30 Apr 2010 = {r_eff:.2f} at every s0 (range {fits.R_eff.min():.3f}–{fits.R_eff.max():.3f}): "
        f"m × s0 ≈ {fits.m_x_s0.mean():.2f} and φ × s0 ≈ {fits.phi_x_s0.mean():.4f}.",
        "",
        "| s0 | m | φ | R0 | RSS |",
        "|---|---|---|---|---|",
    ]
    shown = {round(s, 2) for s in (band[0], head, band[1])} | ({0.39, 0.5, 0.64} & set(f.index))   # F's values too
    for s0 in sorted(shown):
        r = f.loc[s0]
        lines.append(f"| {s0:g} | {r.m:.2f} | {r.phi:.5f} | {r.R0:.2f} | {r.RSS:.2f} |")
    lines += ["", "## Single levers needed to bring R below 1 (F)", "",
              "| Lever | Ceiling | " + " | ".join(c for c in thresholds.columns[2:]) + " |",
              "|---|---|" + "---|" * (len(thresholds.columns) - 2)]
    for r in thresholds.itertuples(index=False):
        vals = [f"{v:.3f}" + (" (impossible)" if v > 1 else " (above ceiling)" if v > r[1] else "") for v in r[2:]]
        lines.append(f"| {r[0]} | {r[1]:g} | " + " | ".join(vals) + " |")
    lines += ["", "R_eff stops the 2010/11 season as it started (the share susceptible on 30 Apr 2010). "
              "R0 prevents outbreaks in a fully susceptible population, so it depends on the s0 behind the fit."]
    lines += ["", f"## Packages (F), 52 weeks from 30 Apr 2010, s0 = {head:g}", "",
              "| ID | Package | κ | δ | p_w | R_eff | R0 at s0 = " + " / ".join(f"{s:g}" for s in sorted({band[0], head, band[1]}))
              + " | Cases | Averted | Peak (week) |", "|---|---|---|---|---|---|---|---|---|---|"]
    for sid, r in t.iterrows():
        r0s = " / ".join(f"{r[f'R0_s0={s:g}']:.2f}" for s in sorted({band[0], head, band[1]}))
        lines.append(f"| {sid} | {r.Name} | {r.kappa:g} | {r.delta:g} | {r.p_w:g} | {r.R_eff:.2f} | {r0s} | "
                     f"{r.cases_52wk:,.0f} | {r.pct_averted:.0f}% | {r.peak:.0f} ({int(r.peak_week)}) |")
    lines += ["", f"## Equal cuts, biting rate vs mosquito density (RQ2), s0 = {head:g}", "",
              "| Cut | Lever | R_eff | Cases averted | Peak reduction |", "|---|---|---|---|---|"]
    for r in rq2.itertuples():
        lines.append(f"| {r.cut:.0%} | {r.lever} | {r.R_eff:.2f} | {r.pct_averted:.0f}% | {r.peak_drop_pct:.0f}% |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skip-validation", action="store_true")
    ap.add_argument("--rebuild-data", action="store_true")
    args = ap.parse_args()
    RESULTS.mkdir(exist_ok=True)
    FIGURES.mkdir(exist_ok=True)
    start = time.time()

    stage("1. Data (AB)")
    cal = check_data(args.rebuild_data)

    stage("2. Model (D, C)")
    check_model()

    stage("3. Calibration (E1): fit m and phi for each fixed s0")
    grid = sorted({round(s, 2) for s in (*config.S0_FIT_GRID, config.S0_HEADLINE, *config.S0_BAND)})
    fits, curves = calibration.calibrate(cal, grid)
    fits.round(6).to_csv(RESULTS / "calibration_fits.csv", index=False)
    curves.round(3).to_csv(RESULTS / "calibration_curves.csv", index=False)
    figures.calibration_fits(curves, fits, FIGURES / "fig1_calibration_fits.png")

    if not args.skip_validation:
        stage("4. Validation (E2): holdout, cross-season, 2024 and residual checks")
        run_validation()

    stage("5. Interventions (F)")
    season = interventions.Season(cal, calibration.fits_lookup(fits), PARAMS)
    head, band = config.S0_HEADLINE, config.S0_BAND
    table, traj = interventions.scenario_table(season, band, head)
    rq2 = interventions.rq2_table(season, head)
    r_eff = float(fits.R_eff.mean())
    thresholds = interventions.single_lever_thresholds(
        r_eff, {s: season.r0_at(s) for s in sorted({band[0], head, band[1]})})
    table.to_csv(RESULTS / "intervention_packages.csv", index=False)
    traj.to_csv(RESULTS / "intervention_curves.csv")
    rq2.to_csv(RESULTS / "intervention_equal_cuts.csv", index=False)
    thresholds.to_csv(RESULTS / "intervention_single_lever_thresholds.csv", index=False)
    print(table[["ID", "Name", "kappa", "delta", "p_w", "R_eff", "pct_averted"]].round(3).to_string(index=False))
    figures.contours(r_eff, season.r0_at(head), (season.r0_at(band[0]), season.r0_at(band[1])),
                     [FIGURES / f"fig2{'abc'[i]}_contours_pw{int(round(p * 10)):02d}.png"
                      for i, p in enumerate(config.P_W_LEVELS)])
    figures.equal_cuts(rq2, FIGURES / "fig3_equal_cuts.png")
    figures.trajectories(traj, table, season.y, FIGURES / "fig4_package_trajectories.png")

    stage("6. Summary")
    summary(cal, fits, table, rq2, thresholds, RESULTS / "SUMMARY.md")
    print(f"  results/SUMMARY.md written. Total time {time.time() - start:.0f} s")


if __name__ == "__main__":
    main()
