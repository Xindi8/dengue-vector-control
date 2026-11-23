# Archive

Files from the team repository that the project no longer runs, kept so the history of each role's work stays visible.

| File | Role | Why it is here |
|---|---|---|
| `Chaya_Dengue_Parameter_Table.csv` | AB (Chaya) | First-pass parameter table, 4 Nov. Replaced by `data_and_literature/parameters/parameters.csv`, which has sourced baselines and ranges. |
| `Chaya_SanJuan_DENV1_Cleaned.csv` | AB (Chaya) | Typed DENV-1 cases for calendar years 2010–2011. The fit uses the 2010/11 season from 30 Apr 2010 with DENV-1-attributed cases (`data_and_literature/data/calibration_san_juan_2010_11.csv`). |
| `E1_calibration_scaffold_2nov.py` | E1 (Zaid) | First calibration scaffold, tested on synthetic data. Replaced by the Stage 2 method in `dengue/calibration.py`. |
| `E2_model_provisional.py` | E2 (Fahad) | Stand-in model used before D's model was final. The checks and tests now use `dengue/model.py`. |
| `F_r0_contour_prototype_3nov.py` | F (Deckar) | First R0 contour sketch with assumed values. Replaced by the contours in `dengue/figures.py`, drawn from the fitted R0 and R_eff. |

None of these files is imported by the package, and they may not run as they are.
