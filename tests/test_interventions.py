"""Tests for the intervention formulas (role F)."""
import numpy as np
import pytest

from dengue import config, data, interventions as iv
from dengue.model import PARAMS, r0

# F's week boundaries, typed by hand in f_model.py (week 35 is 8 days long)
F_WEEK_EDGES = np.array([7 * k for k in range(35)] + [246 + 7 * k for k in range(18)], float)


@pytest.mark.parametrize("kappa, delta, p_w", [(0.35, 0.4, 0.0), (0.6, 0.85, 0.8), (0.0, 0.0, 0.5),
                                               (0.2, 0.1, 0.3)])
def test_closed_form_matches_the_modified_model(kappa, delta, p_w):
    p = dict(PARAMS, m=7.87)
    y0 = np.ones(7)
    pc, yc = iv.apply_controls(p, y0, kappa, delta, p_w)
    assert r0(pc) == pytest.approx(iv.controlled_r(r0(p), kappa, delta, p_w), rel=1e-12)
    assert np.allclose(yc[4:], 1 - delta) and np.allclose(yc[:4], 1.0)


def test_week_edges_match_f():
    season = iv.Season(data.calibration(), {}, PARAMS)
    assert np.array_equal(season.edges, F_WEEK_EDGES)
    assert season.edges_2y[-1] == 365 + 7 * 52


def test_balanced_rule_gives_f_choice_at_s0_half():
    # F's s0 = 0.5 fit: R0 = 2.040; F set delta = 0.65 by hand
    assert iv.balanced_delta(2.040018) == pytest.approx(0.65)


def test_balanced_rule_returns_none_when_impossible():
    assert iv.balanced_delta(50.0) is None


def test_single_lever_thresholds_bring_r_to_one():
    R = 1.44
    t = iv.single_lever_thresholds(R, {}).set_index("lever")["R_eff (stop 2010/11)"]
    assert iv.controlled_r(R, kappa=t["Personal protection, kappa"]) == pytest.approx(1.0)
    assert iv.controlled_r(R, delta=t["Source reduction, delta"]) == pytest.approx(1.0)
    assert iv.controlled_r(R, p_w=t["Wolbachia coverage, p_w"]) == pytest.approx(1.0)


def test_packages_use_the_config_ceilings():
    for sid, _, k, d, _ in config.PACKAGES:
        if d is not None and sid in ("S4", "S6"):
            assert (k, d) == (config.KAPPA_MAX, config.DELTA_MAX)
