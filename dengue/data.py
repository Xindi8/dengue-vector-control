"""Loaders for the AB data and literature package (data_and_literature/).

Role AB: Chaya (first-pass parameter table and cleaned series, now in archive/) and
Fahad (data and literature package). Column meanings are in
data_and_literature/data_dictionary.md.
"""
from __future__ import annotations

import pandas as pd

from . import observation as obs
from .paths import AB_DATA, PARAMETERS_CSV

CALIBRATION_SEASON = "2010/2011"
SECONDARY_SEASON = "2012/2013"


def calibration() -> pd.DataFrame:
    """2010/11 season, 52 weeks from 30 Apr 2010; fit target `denv1_attributed_cases`."""
    return obs.load_calibration(AB_DATA / "calibration_san_juan_2010_11.csv")


def secondary() -> pd.DataFrame:
    """2012/13 season, 52 weeks from 29 Apr 2012, for the cross-season check."""
    return obs.load_calibration(AB_DATA / "secondary_san_juan_2012_13.csv")


def full_series() -> pd.DataFrame:
    """All 23 seasons, 1990/91-2012/13 (1,196 weeks)."""
    return pd.read_csv(AB_DATA / "san_juan_weekly_cases_1990_2013.csv",
                       parse_dates=["week_start_date", "week_end_date"])


def opendengue_weekly() -> pd.DataFrame:
    """Island-wide weekly cases 2014-2024; 2017-18 blank (missing reports)."""
    return pd.read_csv(AB_DATA / "puerto_rico_weekly_opendengue.csv",
                       parse_dates=["week_start_date", "week_end_date"])


def validation_2024() -> pd.DataFrame:
    """The 52 epi weeks of 2024 (6,088 cases), for the qualitative comparison."""
    return pd.read_csv(AB_DATA / "validation_puerto_rico_2024.csv",
                       parse_dates=["week_start_date", "week_end_date"])


def demography() -> pd.DataFrame:
    return pd.read_csv(AB_DATA / "demography_puerto_rico_1990_2024.csv")


def metro_population(year: int) -> int:
    """San Juan metro population (US Census estimate in the CDC repository)."""
    d = demography().set_index("year")
    return int(d.loc[year, "san_juan_msa_population"])


def parameter_table() -> pd.DataFrame:
    return pd.read_csv(PARAMETERS_CSV)


def observed(season: pd.DataFrame):
    return season["denv1_attributed_cases"].to_numpy(float)
