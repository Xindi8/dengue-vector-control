# Literature evidence for the parameter table

MATH 574 Team 2, AB (data and literature). Version of 1 November 2026.

Every value in `parameters/parameters.csv` traces to a source below. Records for the Puerto Rico studies, the incubation periods and the intervention reviews were opened on PubMed on 1 November 2026 (abstracts, plus the full text of Adams et al. 2022). The Section 4 sources (Manore 2014, Gubler 1998, Zahid 2023) were checked against the papers on 3 October 2026. All rates are per day.

## 1. Model parameters (keys as in `dengue_model.py`)

| Key | Value used | Range | What the source reports | Source |
|---|---|---|---|---|
| `N_h` | 2,347,833 | – | San Juan-Carolina-Caguas metro area, 2010 | US Census Bureau estimate in the CDC/PRDH repository |
| `m` | 2 | 1–20 | Vector:host ratio 2:1, varied from 1:1 to 10:1 | Manore et al. 2014, Table 3 |
| `a0` | 0.35 | 0.25–1 | 0.35 bites per mosquito per day for Puerto Rico; the gonotrophic cycle gives 0.25–0.33 | Zahid et al. 2023 |
| `b` | 0.33 | 0.10–0.75 | Mosquito-to-human probability per bite | Manore et al. 2014, Table 3 |
| `c` | 0.33 | 0.10–1.0 | Human-to-mosquito probability per bite; 1.0 upper bound from the review | Manore et al. 2014; Andraud et al. 2012 |
| `eps_h` | 1/5.9 = 0.1695 | 1/10–1/3 | Mean intrinsic incubation 5.9 days, 95% between 3 and 10 days | Chan & Johansson 2012 |
| `eps_v` | 1/10.7 = 0.0935 | 1/15–1/8 | Mean extrinsic incubation 15 days at 25 °C and 6.5 days at 30 °C | Chan & Johansson 2012 |
| `gamma_h` | 1/5 = 0.2 | 1/12–1/2 | Viremia 4–5 days on average, range 2–12 (pp. 487, 490); infectious period 6 days (4–12) | Gubler 1998; Manore et al. 2014 |
| `d_h` | 3.478e-5 | – | Life expectancy at birth, Puerto Rico 2010: 78.717 years | PAHO Core Indicators 2026 (UN WPP 2024) |
| `d_v` | 1/14 = 0.0714 | 1/42–1/4 | Adult lifespan 14 days (8–42); field survival is shorter | Manore et al. 2014; Brady et al. 2013 |

**Extrinsic incubation and temperature.** Chan & Johansson give two means, so the package interpolates log-linearly: EIP(T) = 15 × (6.5/15)^((T − 25)/5) days. San Juan's mean daily temperature from 30 Apr 1990 to 29 Apr 2013 was 27.05 °C (NOAA station RQW00011641). That gives 10.65 days, which the table rounds to 1/10.7, the value D already uses. In 2010/11 the weekly means ran from 24.3 to 29.3 °C, so the extrinsic incubation period ran from 7.3 to 17.0 days. The weekly values are in `data/temperature_san_juan_weekly_1990_2013.csv` if D makes `eps_v` vary with time.

**Mosquito ratio.** No field estimate of `m` exists for San Juan. The 2010/11 fit fixes m × s0 at about 3.9 (E2 plan), so the upper bound has to reach 3.9 / s0. A bound of 20 covers every s0 down to 0.2.

## 2. Reporting fraction φ

The model has no asymptomatic class, so φ is the share of all DENV-1 infections that end up as a counted lab-positive case:

φ = P(symptoms | infection) × P(reported | symptoms) × P(lab-positive | reported)

| Factor | Evidence | Value |
|---|---|---|
| P(symptoms \| infection) | Worldwide, 96 million of 390 million yearly infections are apparent (Bhatt et al. 2013) | about 0.25 |
| | Zhongshan, China, 2013: 11.0 inapparent infections per symptomatic one in the towns with mostly DENV-1 (95% CI 3.7–∞), against 1.0 in the DENV-3 town (Wang et al. 2015) | about 0.083 for DENV-1 |
| | Patillas, Puerto Rico, 2007–08, ages 10–18: 61.1% of new infections were inapparent; symptoms came from a fever diary and were mostly mild (Argüello et al. 2015) | about 0.39 |
| P(reported \| symptoms) | 10–27 cases went unreported for every case reported to passive surveillance in the 1990s (Noyd & Sharp 2015). This is the team plan's p range. | 0.035–0.09 |
| P(lab-positive \| reported) | In 2010, 46.6% of 26,766 suspected cases island-wide were lab-positive (Sharp et al. 2013). Suspected cases include other febrile illnesses, so this bounds the factor only loosely. | up to 1 |

Multiplying the first two rows gives 0.083 × 0.035 = 0.003 up to 0.25 × 0.09 = 0.023, the E2 plan's range. The Puerto Rico cohort's 0.39 counts mild fevers that rarely reach a clinic, so multiplying it by the passive-surveillance factor would overstate φ.

## 3. s0: share susceptible to DENV-1 on 30 April 2010

One season of case data fixes only m × s0 (about 3.9) and φ × s0 (about 0.0019), so s0 must come from outside the case data. No serotype-specific serosurvey exists for San Juan around 2010. The evidence below bounds it.

| Evidence | What it shows | What it means for s0 |
|---|---|---|
| DENV-1 history in these data (`checks/data_checks.md`) | DENV-1 was 8–59% of typed cases in 1990/91–1999/00, then had at most 1 typed case a season from 2000/01 to 2006/07. It came back in 2007/08 (125 typed cases), 2008/09 (472) and 2009/10 (1,297). | Children under about 11 in 2010 had no DENV-1 exposure before 2007. The 2007–10 waves had already used up part of the pool before the fitted season starts. |
| Age structure, Puerto Rico 2010 (PAHO) | 19.5% were under 15 and 34.1% under 25 | Sets the size of the least-exposed cohorts |
| Argüello et al. 2015 (Patillas, 2007–08, ages 10–18) | 49.8% had antibodies to any dengue virus at enrolment | Half of adolescents had never had dengue; the rest were infected while DENV-1 was absent |
| Mohammed et al. 2012 (Red Cross donors, Feb–Mar 2006, mean age 44.6) | 92% had anti-dengue IgG. Of 32 specimens with one predominant serotype, 63% were DENV-3 or DENV-2. | Nearly all adults had been infected, mostly by other serotypes; DENV-1 immunity in adults was not measured |
| Sharp et al. 2013 (2010 epidemic, island-wide) | DENV-1 was 69.0% of RT-PCR positives. About 21% of cases were primary infections, and DENV-1 patients were 4.0 times more likely than DENV-4 patients to have a primary infection. | DENV-1 reached many people with no earlier dengue, as expected if its susceptible pool was large |
| Adams et al. 2023 (Ponce children, 2018–19) | 44% of 9–16-year-olds were seropositive; the force of infection peaked in 2007, 2010 and 2012–13 | 2010 was a high-transmission year island-wide |

**Consistency check.** If φ held steady from 2007 to 2011, the 2,892 DENV-1-attributed cases reported in San Juan over 2007/08–2009/10 mean about 2,892 / φ infections before 30 April 2010. With φ = 0.0019 / s0 and N_h = 2,347,833, that is 0.65 × s0 of the population. So the share that would have been susceptible in 2010 without the 2007–10 waves is about 1.65 × s0.

That share can be bounded from the 2010 age structure, with assumptions that are AB's own:
- 95–100% of under-15s, since DENV-1 circulated little after 1995 and not at all from 2000 to 2007.
- 70–90% of 15–24-year-olds, after some exposure in the early 1990s.
- 30–60% of adults. That range is an assumption: adults carried several infections each, but mostly not DENV-1.

Weighted by age, that gives 0.49–0.72, and so **s0 ≈ 0.29–0.44 on 30 April 2010**.

**Recommendation.** Use **s0 = 0.4** as the central scenario, with 0.2 and 0.6 as the main alternatives. Keep 0.1 and 1.0 only as extremes for sensitivity, since the evidence above argues against both. At the central value the fit implies m ≈ 9.8 (inside Table 1's cap of 10) and φ ≈ 0.0048. That φ matches the literature: for example 0.083 (DENV-1 symptomatic share) × 0.06 (mid reporting rate) × 1 ≈ 0.005.

| s0 | m implied (3.9 / s0) | φ implied (0.0019 / s0) |
|---|---|---|
| 0.2 | 19.5 | 0.0095 |
| 0.3 | 13.0 | 0.0063 |
| 0.4 | 9.8 | 0.0048 |
| 0.6 | 6.5 | 0.0032 |

**What would settle it:** a serotype-specific neutralization (PRNT or FRNT) survey of San Juan residents from around 2010, or age-specific DENV-1 case counts fitted with a catalytic model, as Kada et al. 2024 did for any-serotype dengue in US territories.

## 4. Interventions (for F)

| Lever in the model | Evidence | Value | Setting and caveat |
|---|---|---|---|
| κ, compliance | 44.8% of 4,035 COPA participants used repellent in the past 30 days; 63.4% had screens and 68.7% air conditioning | κ anchor about 0.45 | Ponce, 2018–19 (Adams et al. 2022). Repellent users were not less infected: 33.0% vs 30.2% chikungunya IgG. |
| η, per-bite efficacy | DEET 23.8% gave a mean of 301.5 minutes of complete protection in the laboratory | high while active | Fradin & Day 2002. Field studies found no significant effect of skin repellents on dengue (Bowman et al. 2016). The 0.3–0.9 range is an assumption. |
| a, through house screening | House screening: dengue OR 0.22 (95% CI 0.05–0.93) | about 78% lower risk | Meta-analysis (Bowman et al. 2016). In Ponce, screens gave aOR 0.51 for chikungunya and 0.66 for Zika; air conditioning gave 0.64 and 0.69 (Adams et al. 2022). |
| m, through traps | Three autocidal gravid ovitraps per home: female Ae. aegypti density fell 79% in the former reference area, and untreated areas had 88% more mosquitoes | m cut by about 80% | Southern Puerto Rico (Barrera et al. 2014) |
| m, through traps | The same traps in about 85% of homes: chikungunya infection 26.1% vs 43.8%, aPR 0.50 (95% CI 0.37–0.91); 0.21 for residents who spent most daytime hours in the community | infection about halved | Salinas and Guayama, Puerto Rico (Sharp et al. 2019) |
| m, through source reduction | Community environmental management with water-container covers: dengue OR 0.22 (95% CI 0.15–0.32) | large reduction | Meta-analysis (Bowman et al. 2016) |
| d_v, through spraying | No randomized trial has evaluated space spraying against dengue. Indoor residual spraying: OR 0.67 (95% CI 0.22–2.11), not significant. | no effect size | Bowman et al. 2016. Spraying scenarios rest on assumptions. |

Puerto Rico's 2024 response combined emergency larviciding, mass trapping and spraying (MMWR 74(5), 2025). That confounds any use of 2024 as validation.

## 5. Evidence gaps

1. **s0:** no serotype-specific serosurvey for San Juan near 2010. Section 3 bounds it, but it is still an estimate.
2. **m and mosquito lifespan:** no San Juan field estimates. m is fitted; d_v comes from reviews.
3. **η:** no per-bite field efficacy for repellents or screens. The κ anchor comes from Ponce in 2018–19, not San Juan in 2010.
4. **φ's factors:** the reporting rate dates from the 1990s, and the DENV-1 symptomatic share comes from one outbreak in China with a wide confidence interval.
5. **Spraying:** no effect size from trials.

## References

1. Adams LE et al. Risk factors for infection with chikungunya and Zika viruses in southern Puerto Rico: a community-based cross-sectional seroprevalence survey. PLoS Negl Trop Dis 16(6):e0010416, 2022. https://doi.org/10.1371/journal.pntd.0010416
2. Adams LE et al. Previous dengue infection among children in Puerto Rico and implications for dengue vaccine implementation. Am J Trop Med Hyg 109(2):413–419, 2023. https://doi.org/10.4269/ajtmh.23-0091
3. Andraud M, Hens N, Marais C, Beutels P. Dynamic epidemiological models for dengue transmission: a systematic review of structural approaches. PLoS ONE 7(11):e49085, 2012.
4. Argüello DF et al. Incidence of dengue virus infection in school-aged children in Puerto Rico: a prospective seroepidemiologic study. Am J Trop Med Hyg 92(3):486–491, 2015. https://doi.org/10.4269/ajtmh.14-0231
5. Barrera R et al. Sustained, area-wide control of Aedes aegypti using CDC autocidal gravid ovitraps. Am J Trop Med Hyg 91(6):1269–1276, 2014. https://doi.org/10.4269/ajtmh.14-0426
6. Bhatt S et al. The global distribution and burden of dengue. Nature 496:504–507, 2013. https://doi.org/10.1038/nature12060
7. Bowman LR, Donegan S, McCall PJ. Is dengue vector control deficient in effectiveness or evidence?: systematic review and meta-analysis. PLoS Negl Trop Dis 10(3):e0004551, 2016. https://doi.org/10.1371/journal.pntd.0004551
8. Brady OJ et al. Modelling adult Aedes aegypti and Aedes albopictus survival at different temperatures in laboratory and field settings. Parasit Vectors 6:351, 2013.
9. CDC. Dengue outbreak and response — Puerto Rico, 2024. MMWR 74(5):54–60, 2025. https://www.cdc.gov/mmwr/volumes/74/wr/mm7405a1.htm
10. Chan M, Johansson MA. The incubation periods of dengue viruses. PLoS ONE 7(11):e50972, 2012. https://doi.org/10.1371/journal.pone.0050972
11. Fradin MS, Day JF. Comparative efficacy of insect repellents against mosquito bites. N Engl J Med 347(1):13–18, 2002. https://doi.org/10.1056/NEJMoa011699
12. Gubler DJ. Dengue and dengue hemorrhagic fever. Clin Microbiol Rev 11(3):480–496, 1998.
13. Kada S et al. Age-specific case data reveal varying dengue transmission intensity in US states and territories. PLoS Negl Trop Dis 18(3):e0011143, 2024. https://doi.org/10.1371/journal.pntd.0011143
14. Manore CA, Hickmann KS, Xu S, Wearing HJ, Hyman JM. Comparing dengue and chikungunya emergence and endemic transmission in A. aegypti and A. albopictus. J Theor Biol 356:174–191, 2014.
15. Mohammed H et al. Prevalence of anti-dengue immunoglobulin G antibodies among American Red Cross blood donors in Puerto Rico, 2006. Transfusion 52(8):1652–1656, 2012. https://doi.org/10.1111/j.1537-2995.2011.03492.x
16. Noyd DH, Sharp TM. Recent advances in dengue: relevance to Puerto Rico. P R Health Sci J 34(2):65–70, 2015.
17. Pan American Health Organization. Core Indicators 2026 (release of 13 April 2026). https://opendata.paho.org/en/core-indicators
18. Sharp TM et al. Virus-specific differences in rates of disease during the 2010 dengue epidemic in Puerto Rico. PLoS Negl Trop Dis 7(4):e2159, 2013. https://doi.org/10.1371/journal.pntd.0002159
19. Sharp TM et al. Autocidal gravid ovitraps protect humans from chikungunya virus infection by reducing Aedes aegypti mosquito populations. PLoS Negl Trop Dis 13(7):e0007538, 2019. https://doi.org/10.1371/journal.pntd.0007538
20. Wang T et al. Evaluation of inapparent dengue infections during an outbreak in Southern China. PLoS Negl Trop Dis 9(3):e0003677, 2015. https://doi.org/10.1371/journal.pntd.0003677
21. Zahid MH et al. The biting rate of Aedes aegypti and its variability: a systematic review (1970–2022). PLoS Negl Trop Dis 17(8):e0010831, 2023.
