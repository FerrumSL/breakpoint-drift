# Breakpoint drift: NCBI BioSample isolates (not surveillance data)

2016-2020 isolates read with CLSI 2020 or EUCAST 2020; 2021-2025 isolates with CLSI 2025 or EUCAST 2025. Changes are in percentage points of %R.

1496 isolates and 16832 MIC and disk tests (CLSI and EUCAST) from NCBI BioSample. Shown: the 47 organism x drug x test pairs with at least 10 readable isolates in both periods (all 363 are in `ncbi_drift.csv`). Submitters chose these isolates, often because they were resistant, so the rates describe no population. The edition in force is assumed from the collection year.

- **Observed**: %R in each period under its own guideline, as an antibiogram reports it.
- **Real**: change in the values, with the guideline held fixed (mean of holding either one).
- **Artefact**: the rest of the observed change, caused by the breakpoint revision alone.
- **Verdict**: Fisher's exact test, Benjamini-Hochberg corrected across 47 pairs.
- **Read**: 2016-2020 isolates that the old / new guideline could classify, of those tested. A capped MIC such as <=1 cannot be read against a breakpoint below it.
- **MIC shift**: change in the MIC distribution in doubling dilutions, with 95% CI, from an interval-censored model that uses no breakpoint. * = significant after Benjamini-Hochberg correction.

Verdicts: significant change 7, breakpoint artefact 0, noise 40, not interpretable 0.

## *Acinetobacter baumannii*

| Drug | Test | %R 2016-2020 | %R 2021-2025 | Observed | Real | Artefact | Verdict | Read 2016-2020 | MIC shift |
|---|---|--:|--:|--:|--:|--:|---|--:|--:|
| Imipenem | MIC EUCAST | 62.5 | 100.0 | +37.5 | +37.5 | +0.0 | noise | 16 / 16 of 16 |  |
| Trimethoprim/sulfamethoxazole | MIC CLSI | 56.2 | 28.6 | -27.7 | -27.7 | +0.0 | noise | 16 / 16 of 18 |  |
| Ceftazidime | MIC CLSI | 77.3 | 50.0 | -27.3 | -27.3 | +0.0 | noise | 22 / 22 of 25 | -2.9 (-5.9 to +0.1) |
| Meropenem | MIC EUCAST | 62.5 | 89.5 | +27.0 | +27.0 | +0.0 | noise | 16 / 16 of 16 |  |
| Cefepime | MIC CLSI | 66.7 | 41.7 | -25.0 | -25.0 | +0.0 | noise | 18 / 18 of 20 | +0.1 (-2.0 to +2.3) |
| Ciprofloxacin | MIC EUCAST | 81.2 | 100.0 | +18.8 | +18.8 | +0.0 | noise | 16 / 16 of 16 |  |
| Meropenem | MIC CLSI | 77.3 | 64.3 | -13.0 | -13.0 | +0.0 | noise | 22 / 22 of 25 |  |
| Tobramycin | MIC CLSI | 13.3 | 8.3 | -5.0 | -5.0 | +0.0 | noise | 15 / 15 of 19 |  |

## *Enterococcus faecium*

| Drug | Test | %R 2016-2020 | %R 2021-2025 | Observed | Real | Artefact | Verdict | Read 2016-2020 | MIC shift |
|---|---|--:|--:|--:|--:|--:|---|--:|--:|
| Daptomycin | MIC CLSI | 18.9 | 0.0 | -18.9 | -18.9 | +0.0 | noise | 95 / 95 of 95 | -0.1 (-0.9 to +0.7) |

## *Klebsiella pneumoniae*

| Drug | Test | %R 2016-2020 | %R 2021-2025 | Observed | Real | Artefact | Verdict | Read 2016-2020 | MIC shift |
|---|---|--:|--:|--:|--:|--:|---|--:|--:|
| Ceftriaxone | MIC EUCAST | 18.8 | 96.9 | +78.0 | +78.0 | +0.0 | significant change | 154 / 154 of 154 |  |
| Ceftazidime | MIC EUCAST | 18.7 | 95.0 | +76.3 | +76.3 | +0.0 | significant change | 155 / 155 of 155 |  |
| Amoxicillin/clavulanic acid | MIC EUCAST | 20.1 | 89.9 | +69.8 | +69.8 | +0.0 | significant change | 154 / 154 of 154 |  |
| Piperacillin/tazobactam | MIC EUCAST | 23.5 | 89.9 | +66.4 | +59.4 | +7.0 | significant change | 153 / 153 of 153 |  |
| Polymyxin B | MIC CLSI | 80.0 | 50.3 | -29.7 | -29.7 | +0.0 | significant change | 30 / 30 of 30 |  |
| Cefoxitin | MIC CLSI | 80.6 | 68.1 | -12.5 | -12.5 | +0.0 | noise | 31 / 31 of 31 |  |
| Aztreonam | MIC EUCAST | 75.0 | 87.0 | +12.0 | +12.0 | +0.0 | noise | 12 / 12 of 12 |  |
| Ertapenem | MIC EUCAST | 84.4 | 73.3 | -11.1 | -11.1 | +0.0 | noise | 32 / 32 of 32 |  |
| Ceftolozane/tazobactam | MIC EUCAST | 100.0 | 90.0 | -10.0 | -10.0 | +0.0 | noise | 19 / 19 of 19 |  |
| Meropenem | MIC EUCAST | 75.0 | 65.6 | -9.4 | -9.4 | +0.0 | noise | 32 / 32 of 32 |  |
| Ciprofloxacin | MIC EUCAST | 87.5 | 96.3 | +8.8 | +8.8 | +0.0 | noise | 32 / 32 of 32 |  |
| Ticarcillin/clavulanic acid | MIC EUCAST | 87.1 | 95.0 | +7.9 | +7.9 | +0.0 | noise | 31 / 31 of 31 |  |
| Cefoperazone | MIC CLSI | 87.1 | 94.3 | +7.2 | +7.2 | +0.0 | noise | 31 / 31 of 31 |  |
| Cefepime | MIC EUCAST | 87.5 | 94.4 | +6.9 | +6.9 | +0.0 | noise | 32 / 32 of 32 |  |
| Imipenem | MIC EUCAST | 67.7 | 63.5 | -4.2 | -4.2 | +0.0 | noise | 31 / 31 of 31 |  |
| Ampicillin | MIC EUCAST | 96.9 | 100.0 | +3.1 | +3.1 | +0.0 | noise | 32 / 32 of 32 |  |
| Trimethoprim/sulfamethoxazole | MIC EUCAST | 83.9 | 84.3 | +0.4 | +0.4 | +0.0 | noise | 31 / 31 of 31 |  |

## *Pseudomonas aeruginosa*

| Drug | Test | %R 2016-2020 | %R 2021-2025 | Observed | Real | Artefact | Verdict | Read 2016-2020 | MIC shift |
|---|---|--:|--:|--:|--:|--:|---|--:|--:|
| Ciprofloxacin | disk CLSI | 28.9 | 17.6 | -11.2 | -11.2 | +0.0 | noise | 45 / 45 of 45 |  |
| Ciprofloxacin | MIC EUCAST | 66.7 | 77.5 | +10.8 | +10.8 | +0.0 | noise | 30 / 30 of 30 |  |
| Levofloxacin | disk CLSI | 31.1 | 21.6 | -9.5 | -9.5 | +0.0 | noise | 45 / 45 of 45 |  |
| Cefepime | MIC EUCAST | 72.9 | 63.9 | -8.9 | -8.9 | +0.0 | noise | 59 / 59 of 59 | -0.6 (-1.9 to +0.7) |
| Ceftazidime | disk CLSI | 15.6 | 7.1 | -8.4 | -8.4 | +0.0 | noise | 45 / 45 of 45 |  |
| Ceftazidime | MIC EUCAST | 72.4 | 64.4 | -8.0 | -8.0 | +0.0 | noise | 58 / 58 of 58 | -0.8 (-2.3 to +0.7) |
| Imipenem | MIC EUCAST | 82.8 | 78.0 | -4.8 | -4.8 | +0.0 | noise | 58 / 58 of 58 |  |
| Amikacin | disk CLSI | 4.4 | 1.0 | -3.5 | -3.5 | +0.0 | noise | 45 / 45 of 45 |  |
| Aztreonam | MIC EUCAST | 53.3 | 50.0 | -3.3 | -3.3 | +0.0 | noise | 30 / 30 of 30 | -0.3 (-2.4 to +1.8) |
| Meropenem | MIC EUCAST | 72.4 | 69.2 | -3.2 | -3.2 | +0.0 | noise | 58 / 58 of 58 | +0.1 (-1.6 to +1.8) |
| Ticarcillin/clavulanic acid | MIC EUCAST | 83.3 | 84.6 | +1.3 | +1.3 | +0.0 | noise | 30 / 30 of 30 |  |
| Piperacillin/tazobactam | MIC EUCAST | 71.0 | 71.8 | +0.8 | +0.8 | +0.0 | noise | 31 / 31 of 31 |  |
| Polymyxin B | MIC CLSI | 0.0 | 0.0 | +0.0 | +0.0 | +0.0 | noise | 30 / 30 of 30 | -0.5 (-0.8 to -0.2) * |

## *Staphylococcus aureus*

| Drug | Test | %R 2016-2020 | %R 2021-2025 | Observed | Real | Artefact | Verdict | Read 2016-2020 | MIC shift |
|---|---|--:|--:|--:|--:|--:|---|--:|--:|
| Oxacillin | MIC CLSI | 91.7 | 17.1 | -74.5 | -74.5 | +0.0 | significant change | 12 / 12 of 12 |  |
| Erythromycin | MIC CLSI | 61.5 | 18.0 | -43.5 | -43.5 | +0.0 | significant change | 13 / 13 of 13 |  |
| Trimethoprim/sulfamethoxazole | MIC CLSI | 16.7 | 1.2 | -15.4 | -15.4 | +0.0 | noise | 12 / 12 of 12 |  |
| Gentamicin | MIC CLSI | 0.0 | 2.9 | +2.9 | +2.9 | +0.0 | noise | 12 / 12 of 13 |  |
| Clindamycin | MIC CLSI | 7.7 | 6.3 | -1.4 | -1.4 | +0.0 | noise | 13 / 13 of 13 |  |
| Linezolid | MIC CLSI | 0.0 | 0.0 | +0.0 | +0.0 | +0.0 | noise | 12 / 12 of 12 | -0.5 (-0.7 to -0.2) * |
| Rifampicin | MIC CLSI | 0.0 | 0.0 | +0.0 | +0.0 | +0.0 | noise | 13 / 13 of 13 |  |
| Vancomycin | MIC CLSI | 0.0 | 0.0 | +0.0 | +0.0 | +0.0 | noise | 13 / 13 of 13 | -0.2 (-0.5 to +0.1) |

