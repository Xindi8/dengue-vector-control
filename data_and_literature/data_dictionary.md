# Data dictionary

Dates are ISO (YYYY-MM-DD) and counts are cases per week. In the San Juan files a week is the CDC week: 7-day steps from 1 January, with the week ending 31 December lasting 8 or 9 days. Blank cells mean the source has no value.

## `data/san_juan_weekly_cases_1990_2013.csv`

1,196 rows: 23 seasons (1990/91 to 2012/13) × 52 weeks. San Juan-Carolina-Caguas metro area, by week of symptom onset.

| Column | Type | Units | Description |
|---|---|---|---|
| `season` | text | – | Transmission season, e.g. `2010/2011`. Seasons start in late April, after the historically quietest week. |
| `season_week` | integer | 1–52 | Week of the season |
| `week_start_date` | date | – | First day of the week |
| `week_end_date` | date | – | Last day of the week (the day before the next week starts) |
| `days_in_week` | integer | days | 7, or 8–9 for the week ending 31 December |
| `days_since_start` | integer | days | Days from the season's first day to the start of this week. Use these as the week boundaries when integrating the model. |
| `weekofyear` | integer | 1–52 | CDC week of the calendar year (the plan's `weekofyear`) |
| `total_cases` | integer | cases | All lab-positive cases plus the estimated untested ones, all serotypes (CDC's `total_cases`) |
| `denv1_cases` … `denv4_cases` | integer | cases | Lab-confirmed cases of each serotype (RT-PCR) |
| `other_positive_cases` | integer | cases | Lab-positive cases without a serotype (IgM-positive or IgM conversion) |
| `additional_cases` | integer | cases | CDC's estimate of lab-positive cases among specimens that were not tested |
| `typed_cases` | integer | cases | `denv1_cases` + `denv2_cases` + `denv3_cases` + `denv4_cases` |
| `denv1_share_of_typed` | decimal | fraction | `denv1_cases` / `typed_cases`; blank when no case was typed that week |
| `denv1_attributed_cases` | decimal | cases | `denv1_cases` + (`other_positive_cases` + `additional_cases`) × `denv1_share_of_typed`. **The fit target.** Only meaningful from 2008/09 on, because few positives were typed earlier. |

## `data/calibration_san_juan_2010_11.csv` and `data/secondary_san_juan_2012_13.csv`

52 rows each: the 2010/11 calibration season (30 Apr 2010 to 29 Apr 2011) and the 2012/13 cross-season check (29 Apr 2012 to 29 Apr 2013). Same columns as above, plus three weather columns:

| Column | Type | Units | Description |
|---|---|---|---|
| `temp_mean_c` | decimal | °C | Mean over the week of the daily average temperature (NOAA station RQW00011641, San Juan airport). In this record the daily average is the midpoint of the daily maximum and minimum. |
| `eip_days` | decimal | days | Mean extrinsic incubation period at `temp_mean_c` (Chan & Johansson 2012, log-linear through their 25 °C and 30 °C means; weeks outside 25–30 °C follow the same curve) |
| `eps_v_per_day` | decimal | per day | 1 / `eip_days`, the weekly value of `eps_v` |

## `data/temperature_san_juan_weekly_1990_2013.csv`

1,196 rows, one per CDC week, matching the case file row for row.

| Column | Type | Units | Description |
|---|---|---|---|
| `season`, `season_week`, `week_start_date`, `week_end_date`, `days_in_week` | – | – | As in the case file |
| `temp_mean_c` | decimal | °C | Mean daily average temperature |
| `temp_min_mean_c` | decimal | °C | Mean daily minimum |
| `temp_max_mean_c` | decimal | °C | Mean daily maximum |
| `precip_total_mm` | decimal | mm | Total precipitation in the week |
| `days_with_temperature` | integer | days | Days with a daily average in the record (equal to `days_in_week` in every week) |
| `eip_days`, `eps_v_per_day` | decimal | days, per day | As in the calibration file |

## `data/demography_puerto_rico_1990_2024.csv`

One row per year, 1990–2024.

| Column | Type | Units | Description |
|---|---|---|---|
| `year` | integer | – | Calendar year |
| `san_juan_msa_population` | integer | persons | US Census estimate for the San Juan-Carolina-Caguas metro area (CDC repository); 1990 and 1999–2014 only |
| `puerto_rico_population_thousands` | decimal | thousands | Island-wide population (PAHO indicator 275) |
| `life_expectancy_years` | decimal | years | Life expectancy at birth (PAHO indicator 74), from 1995 |
| `crude_birth_rate_per_1000` | decimal | per 1,000 per year | PAHO indicator 414, from 1995 |
| `crude_death_rate_per_1000` | decimal | per 1,000 per year | PAHO indicator 415, from 1995 |
| `d_h_per_day` | decimal | per day | 1 / (`life_expectancy_years` × 365.25), the model's `d_h` |

## `data/age_structure_puerto_rico_2010.csv`

| Column | Type | Units | Description |
|---|---|---|---|
| `age_group` | text | years | Age band: 0-4, 5-14, 15-24, 25-44, 45-64, 65-74, 75+ |
| `population_thousands` | decimal | thousands | PAHO indicators 358–363, Puerto Rico 2010. PAHO's 0-1 and 1-4 bands overlap at age 1, so 0-4 is the total population (indicator 275) minus the bands from age 5 up. |
| `share` | decimal | fraction | Share of the total population |
| `cumulative_share` | decimal | fraction | Running share from the youngest band |

## `data/puerto_rico_weekly_opendengue.csv` and `data/validation_puerto_rico_2024.csv`

Island-wide weekly cases from OpenDengue V1.3: 574 weeks from 29 Dec 2013 to 28 Dec 2024, and the 52 weeks of epidemiological year 2024.

| Column | Type | Units | Description |
|---|---|---|---|
| `week_start_date`, `week_end_date` | date | – | Sunday-to-Saturday epidemiological week |
| `epi_year`, `epi_week` | integer | – | MMWR epidemiological year and week (week 1 contains 4 January) |
| `year_label_in_source` | integer | – | Year as labelled in OpenDengue. It differs from `epi_year` for weeks that start in late December, so filter by date or `epi_year`. |
| `dengue_total` | integer | cases | Dengue cases reported to PAHO, island-wide. Blank for 2017 and 2018, where the source reports zero for every week. |
| `reporting_note` | text | – | Why a value is blank or doubtful: every week zero in the source (2017, 2018), or 13+ zero weeks in the year (2016, 2019) |
| `opendengue_flag` | text | – | `imputed` or `zero filling` where OpenDengue filled in the value itself (marked in its record ID); blank otherwise. No 2024 week is flagged. |
| `case_definition` | text | – | OpenDengue's standardised case definition (`Total`: suspected and confirmed) |
| `spatial_resolution` | text | – | `Admin0` = all of Puerto Rico |
| `source_id` | text | – | OpenDengue record identifier, which traces each row to its source report |

## `parameters/parameters.csv`

| Column | Description |
|---|---|
| `key` | Parameter name. Model rows use the exact keys of `dengue_model.py`. |
| `symbol` | Symbol in the derivation |
| `description` | What the parameter is |
| `baseline`, `low`, `high` | Value to use and plausible range, in `units` |
| `units` | Rates are per day |
| `status` | `fixed`, `fitted (E1)`, `scenario` or `literature` (a factor that feeds another parameter) |
| `used_by` | Roles that need it |
| `source` | Where the value comes from (details in `literature/evidence.md`) |
| `notes` | Evidence gaps and how to use the value |
