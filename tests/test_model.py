"""Tests for the shared model (roles D and C)."""
import numpy as np
import pandas as pd
import pytest

from dengue import observation as obs
from dengue.model import (PARAM_KEYS, PARAMS, basic_reproduction_number, biting_rate,
                          dengue_rhs, disease_free_equilibrium, r0, r_effective)
from dengue.paths import PARAMETERS_CSV


def test_params_have_exactly_the_model_keys():
    assert set(PARAMS) == set(PARAM_KEYS)


def test_params_match_the_ab_parameter_table():
    table = pd.read_csv(PARAMETERS_CSV).set_index("key")["baseline"]
    for key in PARAM_KEYS:
        assert PARAMS[key] == pytest.approx(float(table[key]), rel=1e-12), key


def test_model_passes_e2_contract():
    assert obs.validate_rhs(dengue_rhs, PARAMS) == []


def test_disease_free_equilibrium_is_a_fixed_point():
    dfe = disease_free_equilibrium(PARAMS)
    assert np.allclose(dengue_rhs(0.0, dfe, PARAMS), 0.0, atol=1e-9)


def test_populations_are_conserved():
    p = dict(PARAMS, m=5.0)
    N_h, N_v = p["N_h"], p["m"] * p["N_h"]
    y = np.array([0.5 * N_h, 0.1 * N_h, 0.1 * N_h, 0.3 * N_h, 0.8 * N_v, 0.15 * N_v, 0.05 * N_v])
    d = dengue_rhs(0.0, y, p)
    assert abs(d[:4].sum()) < 1e-6 * N_h
    assert abs(d[4:].sum()) < 1e-6 * N_v


@pytest.mark.parametrize("seed", range(5))
def test_next_generation_matrix_equals_closed_form(seed):
    rng = np.random.default_rng(seed)
    p = dict(PARAMS, m=rng.uniform(0.5, 30), a0=rng.uniform(0.2, 1.0), kappa=rng.uniform(0, 0.6),
             eta=rng.uniform(0, 1), b=rng.uniform(0.1, 1), c=rng.uniform(0.1, 1))
    ngm, closed = basic_reproduction_number(p)
    assert ngm == pytest.approx(closed, rel=1e-10)


def test_r0_scales_with_a_and_with_root_m():
    base = dict(PARAMS, m=4.0)
    assert r0(dict(base, a0=2 * base["a0"])) == pytest.approx(2 * r0(base), rel=1e-12)
    assert r0(dict(base, m=4 * base["m"])) == pytest.approx(2 * r0(base), rel=1e-12)


def test_personal_protection_scales_the_biting_rate():
    p = dict(PARAMS, kappa=0.4, eta=0.5)
    assert biting_rate(p) == pytest.approx(PARAMS["a0"] * 0.8)
    assert r0(p) == pytest.approx(0.8 * r0(dict(PARAMS)), rel=1e-12)


def test_r_effective_is_r0_times_root_s0():
    p = dict(PARAMS, m=9.81)
    assert r_effective(p, 0.4) == pytest.approx(r0(p) * np.sqrt(0.4))
