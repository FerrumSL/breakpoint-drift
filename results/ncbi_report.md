# Breakpoint drift: NCBI BioSample isolates (not surveillance data)

2011-2018 isolates read with CLSI 2018 or EUCAST 2018; 2019-2025 isolates with CLSI 2025 or EUCAST 2025. Changes are in percentage points of %R.

1783 isolates and 21555 MIC and disk tests (CLSI and EUCAST) from NCBI BioSample. Shown: the 34 organism x drug x test pairs with at least 10 readable isolates in both periods (all 455 are in `ncbi_drift.csv`). Submitters chose these isolates, often because they were resistant, so the rates describe no population. The edition in force is assumed from the collection year.

- **Observed**: %R in each period under its own guideline, as an antibiogram reports it.
- **Real**: change in the values, with the guideline held fixed (mean of holding either one).
- **Artefact**: the rest of the observed change, caused by the breakpoint revision alone.
- **Verdict**: Fisher's exact test, Benjamini-Hochberg corrected across 34 pairs.
- **Read**: 2011-2018 isolates that the old / new guideline could classify, of those tested. A capped MIC such as <=1 cannot be read against a breakpoint below it.
- **MIC shift**: change in the MIC distribution in doubling dilutions, with 95% CI, from an interval-censored model that uses no breakpoint. * = significant after Benjamini-Hochberg correction.

Verdicts: significant change 12, breakpoint artefact 3, noise 19, not interpretable 0.

## *Acinetobacter baumannii*

| Drug | Test | %R 2011-2018 | %R 2019-2025 | Observed | Real | Artefact | Verdict | Read 2011-2018 | MIC shift |
|---|---|--:|--:|--:|--:|--:|---|--:|--:|
| Ampicillin/sulbactam | MIC CLSI | 11.5 | 61.5 | +50.0 | +50.0 | +0.0 | significant change | 26 / 26 of 46 | +1.3 (-0.7 to +3.2) |
| Polymyxin B | MIC CLSI | 0.0 | 31.2 | +31.2 | +31.2 | +0.0 | significant change | 23 / 23 of 23 | +1.5 (+0.6 to +2.4) * |
| Meropenem | MIC CLSI | 45.0 | 75.0 | +30.0 | +30.0 | +0.0 | noise | 40 / 40 of 62 | +3.5 (-0.3 to +7.3) |
| Levofloxacin | MIC CLSI | 22.7 | 50.0 | +27.3 | +27.3 | +0.0 | noise | 22 / 22 of 61 | -1.8 (-6.0 to +2.3) |
| Ceftazidime | MIC CLSI | 39.5 | 65.0 | +25.5 | +25.5 | +0.0 | noise | 38 / 38 of 62 | +0.6 (-1.7 to +2.9) |
| Piperacillin/tazobactam | MIC CLSI | 85.2 | 64.3 | -20.9 | -20.9 | +0.0 | noise | 27 / 27 of 27 |  |
| Cefepime | MIC CLSI | 69.8 | 50.0 | -19.8 | -19.8 | +0.0 | noise | 53 / 53 of 76 | -1.5 (-4.4 to +1.5) |
| Tobramycin | MIC CLSI | 2.4 | 21.4 | +19.0 | +19.0 | +0.0 | noise | 41 / 41 of 60 | +0.6 (-1.7 to +2.8) |
| Trimethoprim/sulfamethoxazole | MIC CLSI | 54.5 | 37.5 | -17.0 | -17.0 | +0.0 | noise | 33 / 33 of 59 |  |
| Imipenem | MIC CLSI | 57.1 | 68.8 | +11.6 | +11.6 | +0.0 | noise | 49 / 49 of 71 | +0.4 (-6.6 to +7.4) |
| Ciprofloxacin | MIC CLSI | 73.3 | 64.3 | -9.0 | -9.0 | +0.0 | noise | 45 / 45 of 85 | -0.6 (-5.1 to +3.8) |
| Colistin | MIC CLSI | 4.8 | 9.1 | +4.3 | +4.3 | +0.0 | noise | 63 / 63 of 63 | +1.5 (+0.3 to +2.7) * |
| Gentamicin | MIC CLSI | 32.3 | 30.8 | -1.5 | -1.5 | +0.0 | noise | 31 / 31 of 62 | -0.8 (-3.1 to +1.5) |

## *Enterococcus faecium*

| Drug | Test | %R 2011-2018 | %R 2019-2025 | Observed | Real | Artefact | Verdict | Read 2011-2018 | MIC shift |
|---|---|--:|--:|--:|--:|--:|---|--:|--:|
| Daptomycin | MIC CLSI | 0.0 | 0.0 | +0.0 | -9.5 | +9.5 | noise | 77 / 95 of 95 | -0.1 (-0.9 to +0.7) |

## *Escherichia coli*

| Drug | Test | %R 2011-2018 | %R 2019-2025 | Observed | Real | Artefact | Verdict | Read 2011-2018 | MIC shift |
|---|---|--:|--:|--:|--:|--:|---|--:|--:|
| Amoxicillin/clavulanic acid | disk EUCAST | 40.0 | 100.0 | +60.0 | +60.0 | +0.0 | significant change | 20 / 20 of 20 |  |
| Gentamicin | MIC CLSI | 7.6 | 43.8 | +36.1 | +4.6 | +31.5 | breakpoint artefact | 171 / 236 of 244 | +2.2 (-7.3 to +11.8) |
| Ciprofloxacin | MIC CLSI | 47.4 | 81.2 | +33.8 | +3.8 | +30.0 | breakpoint artefact | 116 / 235 of 243 | -6.6 (-13.5 to +0.3) |
| Amikacin | MIC CLSI | 1.7 | 18.8 | +17.0 | +3.3 | +13.7 | breakpoint artefact | 234 / 232 of 243 | +1.8 (+0.4 to +3.1) * |
| Colistin | MIC CLSI | 3.0 | 15.4 | +12.4 | +10.8 | +1.5 | noise | 66 / 66 of 66 |  |
| Trimethoprim/sulfamethoxazole | MIC CLSI | 61.3 | 53.3 | -8.0 | -8.0 | +0.0 | noise | 243 / 243 of 244 |  |
| Cefoxitin | MIC CLSI |  | 11.8 |  |  |  | noise | 0 / 155 of 187 | -1.8 (-3.8 to +0.1) |

## *Klebsiella pneumoniae*

| Drug | Test | %R 2011-2018 | %R 2019-2025 | Observed | Real | Artefact | Verdict | Read 2011-2018 | MIC shift |
|---|---|--:|--:|--:|--:|--:|---|--:|--:|
| Ceftriaxone | MIC EUCAST | 0.0 | 58.8 | +58.8 | +58.8 | +0.0 | significant change | 24 / 24 of 24 |  |
| Ceftazidime | MIC EUCAST | 0.0 | 58.0 | +58.0 | +58.0 | +0.0 | significant change | 24 / 24 of 24 |  |
| Amoxicillin/clavulanic acid | MIC EUCAST | 8.3 | 55.9 | +47.6 | +47.6 | +0.0 | significant change | 24 / 24 of 24 | +8.8 (+5.7 to +12.0) * |
| Piperacillin/tazobactam | MIC EUCAST | 17.4 | 59.8 | +42.4 | +36.7 | +5.7 | significant change | 23 / 23 of 23 | +9.9 (+5.8 to +14.0) * |

## *Pseudomonas aeruginosa*

| Drug | Test | %R 2011-2018 | %R 2019-2025 | Observed | Real | Artefact | Verdict | Read 2011-2018 | MIC shift |
|---|---|--:|--:|--:|--:|--:|---|--:|--:|
| Meropenem | disk CLSI | 78.3 | 4.1 | -74.2 | -74.2 | +0.0 | significant change | 23 / 23 of 23 |  |
| Ceftazidime | disk CLSI | 56.5 | 7.6 | -48.9 | -48.9 | +0.0 | significant change | 23 / 23 of 23 |  |
| Ciprofloxacin | disk CLSI | 47.8 | 17.4 | -30.4 | -38.4 | +8.0 | significant change | 23 / 23 of 23 |  |
| Levofloxacin | disk CLSI | 52.2 | 21.7 | -30.4 | -40.2 | +9.8 | significant change | 23 / 23 of 23 |  |
| Amikacin | disk CLSI | 30.4 | 0.7 | -29.7 | -29.7 | +0.0 | significant change | 23 / 23 of 23 |  |
| Meropenem | MIC EUCAST | 87.5 | 68.5 | -19.0 | -19.0 | +0.0 | noise | 16 / 16 of 16 | -0.3 (-2.7 to +2.2) |
| Imipenem | MIC EUCAST | 93.8 | 77.5 | -16.3 | -22.4 | +6.1 | noise | 16 / 16 of 16 | -2.5 (-5.4 to +0.3) |
| Cefepime | MIC EUCAST | 80.0 | 66.7 | -13.3 | -13.3 | +0.0 | noise | 15 / 15 of 15 | -0.5 (-2.4 to +1.3) |
| Ceftazidime | MIC EUCAST | 68.8 | 68.6 | -0.1 | -0.1 | +0.0 | noise | 16 / 16 of 16 | -0.2 (-2.3 to +2.0) |

