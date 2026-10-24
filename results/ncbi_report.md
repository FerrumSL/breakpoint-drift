# Breakpoint drift: NCBI BioSample isolates (not surveillance data)

2011-2018 isolates read with CLSI 2018, 2019-2025 isolates with CLSI 2025. Changes are in percentage points of %R.

1249 isolates and 12147 CLSI MIC tests from NCBI BioSample. Shown: the 20 organism x drug pairs with at least 10 readable isolates in both periods (all 189 are in `ncbi_drift.csv`). Submitters chose these isolates, often because they were resistant, so the rates describe no population. The edition in force is assumed from the collection year.

- **Observed**: %R in each year under its own guideline, as an antibiogram reports it.
- **Real**: change in the MICs, with the guideline held fixed (mean of holding either one).
- **Artefact**: the rest of the observed change, caused by the breakpoint revision alone.
- **Verdict**: Fisher's exact test, Benjamini-Hochberg corrected across 20 pairs.
- **Read**: 2011-2018 isolates that CLSI 2018 / CLSI 2025 could classify, of those tested. A capped MIC such as <=1 cannot be read against a breakpoint below it.

Verdicts: significant change 1, breakpoint artefact 2, noise 17, not interpretable 0.

## *Acinetobacter baumannii*

| Drug | %R 2011-2018 | %R 2019-2025 | Observed | Real | Artefact | Verdict | Read 2011-2018 |
|---|--:|--:|--:|--:|--:|---|--:|
| Ampicillin/sulbactam | 11.5 | 61.5 | +50.0 | +50.0 | +0.0 | significant change | 26 / 26 of 46 |
| Polymyxin B | 0.0 | 31.2 | +31.2 | +31.2 | +0.0 | noise | 23 / 23 of 23 |
| Meropenem | 45.0 | 75.0 | +30.0 | +30.0 | +0.0 | noise | 40 / 40 of 62 |
| Levofloxacin | 22.7 | 50.0 | +27.3 | +27.3 | +0.0 | noise | 22 / 22 of 61 |
| Ceftazidime | 39.5 | 65.0 | +25.5 | +25.5 | +0.0 | noise | 38 / 38 of 62 |
| Piperacillin/tazobactam | 85.2 | 64.3 | -20.9 | -20.9 | +0.0 | noise | 27 / 27 of 27 |
| Cefepime | 69.8 | 50.0 | -19.8 | -19.8 | +0.0 | noise | 53 / 53 of 76 |
| Tobramycin | 2.4 | 21.4 | +19.0 | +19.0 | +0.0 | noise | 41 / 41 of 60 |
| Trimethoprim/sulfamethoxazole | 54.5 | 37.5 | -17.0 | -17.0 | +0.0 | noise | 33 / 33 of 59 |
| Imipenem | 57.1 | 68.8 | +11.6 | +11.6 | +0.0 | noise | 49 / 49 of 71 |
| Ciprofloxacin | 73.3 | 64.3 | -9.0 | -9.0 | +0.0 | noise | 45 / 45 of 85 |
| Colistin | 4.8 | 9.1 | +4.3 | +4.3 | +0.0 | noise | 63 / 63 of 63 |
| Gentamicin | 32.3 | 30.8 | -1.5 | -1.5 | +0.0 | noise | 31 / 31 of 62 |

## *Enterococcus faecium*

| Drug | %R 2011-2018 | %R 2019-2025 | Observed | Real | Artefact | Verdict | Read 2011-2018 |
|---|--:|--:|--:|--:|--:|---|--:|
| Daptomycin | 0.0 | 0.0 | +0.0 | -9.5 | +9.5 | noise | 77 / 95 of 95 |

## *Escherichia coli*

| Drug | %R 2011-2018 | %R 2019-2025 | Observed | Real | Artefact | Verdict | Read 2011-2018 |
|---|--:|--:|--:|--:|--:|---|--:|
| Gentamicin | 7.6 | 43.8 | +36.1 | +4.6 | +31.5 | breakpoint artefact | 171 / 236 of 244 |
| Ciprofloxacin | 47.4 | 81.2 | +33.8 | +3.8 | +30.0 | noise | 116 / 235 of 243 |
| Amikacin | 1.7 | 18.8 | +17.0 | +3.3 | +13.7 | breakpoint artefact | 234 / 232 of 243 |
| Colistin | 3.0 | 15.4 | +12.4 | +10.8 | +1.5 | noise | 66 / 66 of 66 |
| Trimethoprim/sulfamethoxazole | 61.3 | 53.3 | -8.0 | -8.0 | +0.0 | noise | 243 / 243 of 244 |
| Cefoxitin |  | 11.8 |  |  |  | noise | 0 / 155 of 187 |

