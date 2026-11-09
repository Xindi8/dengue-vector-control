# Data checks

Written by `build_ab_package.py`; every number below is computed from `raw/`.

## San Juan weekly cases (CDC/PRDH, 1990/91-2012/13)

- **PASS**: 23 seasons x 52 weeks = 1196 rows; first week starts 30 Apr 1990, last week ends 29 Apr 2013.
- **PASS**: No duplicate weeks.
- **PASS**: No missing or negative counts.
- **PASS**: total_cases equals DENV-1 to -4 + untyped + estimated untested in every week.
- **PASS**: Week lengths: 1173 weeks of 7 days, 17 weeks of 8 days, 6 weeks of 9 days. Every long week ends on 31 Dec (weeks run from 1 Jan in 7-day steps), so fit on real dates, not week numbers.
- **PASS**: Weeks with untyped cases but no typed case, where the DENV-1 share is unknown and attribution gives 0: 96 of 1196 in the full series, 3 from 2008/09 on, 0 in 2010/11 and 2012/13. Before 2008 few positives were typed, so `denv1_attributed_cases` is only meaningful from 2008/09.
- **PASS**: Geography: San Juan-Carolina-Caguas Metropolitan Statistical Area (raw/cdc/san_juan_metadata.md).
- **PASS**: Case definition: laboratory-positive (serotype by RT-PCR, or IgM-positive), plus `additional_cases`, which CDC estimated for weeks when not all specimens were tested. Dated by week of symptom onset.

## Season 2010/2011

- **PASS**: 30 Apr 2010 to 29 Apr 2011, 52 weeks, 365 days.
- **PASS**: Totals: all lab-positive 4,696; typed DENV-1 2,594; DENV-1 attributed 3,863.9; estimated untested 22.
- **PASS**: Serotype mix of typed cases: DENV-1 82.2%, DENV-2 1.1%, DENV-3 0.0%, DENV-4 16.7%.
- **PASS**: Share of lab-positive cases with a serotype, weeks with 20+ cases: 0.55-0.94.
- **PASS**: Weeks with untyped cases but no typed case (DENV-1 share unknown, attributed as 0): 0.
- **PASS**: DENV-1-attributed peak: 228.6 in the week of 13 Aug 2010 (season week 16); weeks at half the peak or more: 15.

## Season 2012/2013

- **PASS**: 29 Apr 2012 to 29 Apr 2013, 52 weeks, 366 days.
- **PASS**: Totals: all lab-positive 5,283; typed DENV-1 2,692; DENV-1 attributed 4,091.2; estimated untested 38.
- **PASS**: Serotype mix of typed cases: DENV-1 77.5%, DENV-2 0.0%, DENV-3 0.0%, DENV-4 22.4%.
- **PASS**: Share of lab-positive cases with a serotype, weeks with 20+ cases: 0.52-0.92.
- **PASS**: Weeks with untyped cases but no typed case (DENV-1 share unknown, attributed as 0): 0.
- **PASS**: DENV-1-attributed peak: 185.8 in the week of 01 Jan 2013 (season week 36); weeks at half the peak or more: 18.

## DENV-1 history (used for the s0 evidence)

- **PASS**: Seasons with at most 1 typed DENV-1 case: 2000/2001, 2001/2002, 2002/2003, 2003/2004, 2004/2005, 2005/2006, 2006/2007.
- **PASS**: Typed DENV-1 cases before the calibration season: 2007/2008 125, 2008/2009 472, 2009/2010 1297; DENV-1 attributed over those three seasons: 2,892.

## Population

- **PASS**: N_h = 2,347,833 (San Juan metro, 2010, US Census estimate in the CDC repo) = 63% of Puerto Rico's 3.71 million (PAHO).
- **PASS**: d_h = 1 / (78.717 years x 365.25) = 3.4781e-05 per day (PAHO life expectancy, Puerto Rico 2010).
- **PASS**: PAHO's 0-1 and 1-4 age bands overlap at age 1: together they exceed the total by 44.1 thousand, about one year of births (42.3 thousand). The clean table uses 0-4 = total minus the bands from age 5 up.
- **PASS**: Puerto Rico 2010 age structure: 19.5% under 15, 34.1% under 25.

## Parameter table

- **PASS**: R0 from the closed form with the baseline values: 1.03 at m = 2, and 2.27 at m = 9.75 (the value implied at s0 = 0.4).

## Temperature (NOAA GHCN-Daily, station RQW00011641, San Juan airport)

- **PASS**: 8,401 of 8,401 days present; days without a daily mean: 0.
- **PASS**: Mean daily temperature 30 Apr 1990-29 Apr 2013: 27.05 C, so mean extrinsic incubation 10.65 days and eps_v = 0.0939 per day (parameter table: 0.0935).
- **PASS**: 2010/11 weekly mean temperature 24.3-29.3 C, so extrinsic incubation 7.3-17.0 days.

## Puerto Rico weekly cases, OpenDengue V1.3 (island-wide)

- **PASS**: Epi year 2024: 52 weeks, 2023-12-31 to 2024-12-28, 6,088 cases.
- **NOTE**: Weeks of epi year 2024 labelled with another year in the source: 1 (filter by date or `epi_year`, not by the source's year label).
- **PASS**: Weekly 2024 total matches OpenDengue's annual 2024 figure (6,088).
- **PASS**: 2024 peak 217 in the week of 2024-09-29; weeks at half the peak or more: 30.
- **PASS**: Weekly coverage 2013-12-29 to 2024-12-28, 0 gaps; years with 53 epi weeks: 2014, 2020.
- **NOTE**: Years where every week is zero in the source, set to missing: 2017, 2018. Years flagged as likely incomplete: 2016 (21 zero weeks), 2019 (43 zero weeks).
- **PASS**: Weeks OpenDengue filled in itself (marked in the record ID): 45, by year and type: 2015 imputed 16, 2016 imputed 3, 2017 imputed 12, 2017 zero filling 1, 2019 imputed 2, 2019 zero filling 5, 2021 imputed 2, 2022 imputed 4. Epi year 2024 has 0.
- **PASS**: Case definition 'Total' (suspected and confirmed as reported to PAHO), island-wide, no serotype, no municipality.
