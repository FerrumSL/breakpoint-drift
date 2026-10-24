# NARMS cecal isolates: CLSI 2020 against CLSI 2023

2019-2020 isolates read with CLSI 2020; 2022-2023 isolates with CLSI 2023. Changes are in percentage points of %R.

11744 Salmonella and E. coli isolates from food-animal guts at slaughter.

- **Observed**: %R in each period under its own guideline, as an antibiogram reports it.
- **Real**: change in the values, with the guideline held fixed (mean of holding either one).
- **Artefact**: the rest of the observed change, caused by the breakpoint revision alone.
- **Verdict**: Fisher's exact test, Benjamini-Hochberg corrected across 24 pairs.
- **Read**: 2019-2020 isolates that the old / new guideline could classify, of those tested. A capped MIC such as <=1 cannot be read against a breakpoint below it.
- **MIC shift**: change in the MIC distribution in doubling dilutions, with 95% CI, from an interval-censored model that uses no breakpoint. * = significant after Benjamini-Hochberg correction.

Verdicts: significant change 5, breakpoint artefact 1, noise 18, not interpretable 0.

## *Escherichia coli*

| Drug | Test | %R 2019-2020 | %R 2022-2023 | Observed | Real | Artefact | Verdict | Read 2019-2020 | MIC shift |
|---|---|--:|--:|--:|--:|--:|---|--:|--:|
| Tetracycline | MIC CLSI | 43.5 | 48.7 | +5.2 | +5.2 | +0.0 | significant change | 3951 / 3951 of 3951 |  |
| Ampicillin | MIC CLSI | 17.7 | 22.0 | +4.3 | +4.3 | +0.0 | significant change | 3952 / 3952 of 3952 | +0.4 (+0.3 to +0.5) * |
| Trimethoprim/sulfamethoxazole | MIC CLSI | 4.8 | 6.4 | +1.6 | +1.6 | +0.0 | significant change | 3951 / 3951 of 3951 |  |
| Chloramphenicol | MIC CLSI | 5.7 | 7.2 | +1.5 | +1.5 | +0.0 | noise | 3952 / 3952 of 3952 | +0.1 (+0.0 to +0.1) * |
| Gentamicin | MIC CLSI | 5.3 | 6.7 | +1.4 | +0.3 | +1.1 | noise | 3952 / 3952 of 3952 | -0.1 (-0.2 to -0.0) * |
| Nalidixic acid | MIC CLSI | 2.9 | 4.0 | +1.1 | +1.1 | +0.0 | noise | 3937 / 3937 of 3937 | +0.2 (+0.1 to +0.3) * |
| Ceftriaxone | MIC CLSI | 3.3 | 4.4 | +1.1 | +1.1 | +0.0 | noise | 3952 / 3952 of 3952 |  |
| Amoxicillin/clavulanic acid | MIC CLSI | 2.6 | 3.2 | +0.6 | +0.6 | +0.0 | noise | 3952 / 3952 of 3952 | +0.2 (+0.1 to +0.2) * |
| Cefoxitin | MIC CLSI | 2.6 | 3.1 | +0.6 | +0.6 | +0.0 | noise | 3950 / 3950 of 3950 | -0.0 (-0.1 to +0.0) |
| Ciprofloxacin | MIC CLSI | 1.0 | 1.4 | +0.4 | +0.4 | +0.0 | noise | 3952 / 3952 of 3952 |  |
| Colistin | MIC CLSI | 0.1 | 0.0 | -0.1 | -0.1 | +0.0 | noise | 1702 / 1702 of 1702 |  |
| Meropenem | MIC CLSI | 0.0 | 0.0 | +0.0 | +0.0 | +0.0 | noise | 3952 / 3952 of 3952 |  |

## *Salmonella enterica*

| Drug | Test | %R 2019-2020 | %R 2022-2023 | Observed | Real | Artefact | Verdict | Read 2019-2020 | MIC shift |
|---|---|--:|--:|--:|--:|--:|---|--:|--:|
| Nalidixic acid | MIC CLSI | 8.4 | 12.6 | +4.2 | +4.2 | +0.0 | significant change | 2732 / 2732 of 2732 | +0.3 (+0.2 to +0.3) * |
| Gentamicin | MIC CLSI | 1.7 | 5.8 | +4.1 | +0.7 | +3.4 | breakpoint artefact | 2734 / 2734 of 2734 | -0.2 (-0.3 to -0.1) * |
| Tetracycline | MIC CLSI | 26.9 | 30.9 | +4.0 | +4.0 | +0.0 | significant change | 2734 / 2734 of 2734 |  |
| Chloramphenicol | MIC CLSI | 5.2 | 6.7 | +1.5 | +1.5 | +0.0 | noise | 2734 / 2734 of 2734 | -0.0 (-0.1 to +0.1) |
| Amoxicillin/clavulanic acid | MIC CLSI | 1.9 | 2.5 | +0.6 | +0.6 | +0.0 | noise | 2734 / 2734 of 2734 |  |
| Trimethoprim/sulfamethoxazole | MIC CLSI | 3.1 | 2.5 | -0.5 | -0.5 | +0.0 | noise | 2734 / 2734 of 2734 |  |
| Ampicillin | MIC CLSI | 9.0 | 8.7 | -0.3 | -0.3 | +0.0 | noise | 2734 / 2734 of 2734 |  |
| Cefoxitin | MIC CLSI | 1.9 | 2.2 | +0.3 | +0.3 | +0.0 | noise | 2734 / 2734 of 2734 | -0.0 (-0.1 to +0.0) |
| Ceftriaxone | MIC CLSI | 4.7 | 4.9 | +0.2 | +0.2 | +0.0 | noise | 2734 / 2734 of 2734 |  |
| Ciprofloxacin | MIC CLSI | 0.3 | 0.4 | +0.2 | +0.2 | +0.0 | noise | 2734 / 2734 of 2734 |  |
| Colistin | MIC CLSI | 0.3 | 0.5 | +0.2 | +0.2 | +0.0 | noise | 1157 / 1157 of 1157 |  |
| Meropenem | MIC CLSI | 0.0 | 0.0 | +0.0 | +0.0 | +0.0 | noise | 2734 / 2734 of 2734 |  |

