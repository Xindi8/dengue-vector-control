"""The integrated package must reproduce each member's own results.

tests/reference/ holds the result files E1 and F committed to the team repository
(outputs/e1_fit_results.csv and F's outputs/). Those files are rounded to 6 decimals,
so fitted values are compared to that precision. F's scenarios are recomputed from
F's own rounded fits. The ODE solver runs at rtol 1e-8, and its adaptive steps differ
slightly between computers (Linux, Mac, Windows), so these comparisons allow 1e-6:
loose enough for any computer, while a wrong formula would miss by 1e-3 or more.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from dengue import calibration, data, interventions as iv
from dengue.model import PARAMS

REF = Path(__file__).parent / "reference"


@pytest.fixture(scope="module")
def cal():
    return data.calibration()


def assert_frame_close(ours, ref, rtol, atol):
    """Every column of ref, compared by name."""
    for col in ref.columns:
        if ref[col].dtype.kind in "fi":
            assert np.allclose(ours[col].to_numpy(float), ref[col].to_numpy(float), rtol=rtol, atol=atol), col
        else:
            assert (ours[col].to_numpy() == ref[col].to_numpy()).all(), col


@pytest.mark.parametrize("role, file", [("E1", "E1_e1_fit_results.csv"), ("F", "F_f_fit_results.csv")])
def test_calibration_reproduces_member_fits(cal, role, file):
    ref = pd.read_csv(REF / file)
    fits, _ = calibration.calibrate(cal, tuple(ref.s0), verbose=False)
    cols = [c for c in ref.columns if c != "solver_status"]   # status codes 2 and 3 both mean converged
    assert_frame_close(fits.reset_index(drop=True), ref[cols], rtol=1e-5, atol=1e-6)


@pytest.fixture(scope="module")
def f_season(cal):
    ref = pd.read_csv(REF / "F_f_fit_results.csv")
    return iv.Season(cal, {round(r.s0, 2): (r.m, r.phi) for r in ref.itertuples()}, PARAMS)


def test_scenarios_reproduce_f(f_season):
    table, curves = iv.scenario_table(f_season, (0.39, 0.64), headline=0.5)
    ref = pd.read_csv(REF / "F_f_scenario_results.csv")
    assert sorted(table.columns) == sorted(ref.columns)
    assert_frame_close(table, ref, rtol=1e-6, atol=1e-6)
    ref_curves = pd.read_csv(REF / "F_f_scenario_curves.csv", index_col="week")
    assert np.allclose(curves[ref_curves.columns].to_numpy(), ref_curves.to_numpy(), rtol=1e-6, atol=1e-6)


def test_equal_cuts_reproduce_f(f_season):
    ref = pd.read_csv(REF / "F_f_rq2_equal_cuts.csv")
    assert_frame_close(iv.rq2_table(f_season, 0.5), ref, rtol=1e-6, atol=1e-6)
