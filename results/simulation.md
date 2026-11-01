# Simulation

Written by scripts/simulation.R (one seed per cell from 2026, 200 datasets per cell, B = 999, 186 cells, 49.9 minutes).
Setup in the script header. Bias in percentage points, the rest in percent. Ranges are over the
cells with that n; per-cell results, with the true components, are in simulation.csv.

## %R, base cells

| | n = 200 | n = 1000 | n = 5000 |
|---|--:|--:|--:|
| Bias, period, largest absolute | 0.41 | 0.31 | 0.17 |
| Bias, guideline, largest absolute | 0.16 | 0.09 | 0.03 |
| Coverage, period, no capping | 93.0-98.0 | 91.5-97.0 | 88.0-98.5 |
| Coverage, period, 30% capped | 92.0-97.0 | 92.5-98.0 | 93.5-96.5 |
| Coverage, guideline, shift > 0, no capping | 92.5-97.5 | 92.5-96.0 | 93.0-96.5 |
| Coverage, guideline, shift > 0, 30% capped | 90.0-96.0 | 93.0-96.5 | 94.0-97.0 |
| Period shown, no period change | 2.0-6.5 | 2.0-8.5 | 1.5-12.0 |
| Period shown, period change | 16.5-40.0 | 64.0-99.0 | 100.0 |
| Guideline shown, shift 0 | 0.0 | 0.0 | 0.0 |
| Guideline shown, shift > 0 | 100.0 | 100.0 | 100.0 |
| Larger named correctly, shift > 0 and period change | 1.0-99.5 | 2.0-100.0 | 2.0-100.0 |
| Larger named wrongly, any cell | 0.0-5.0 | 0.0-6.5 | 0.0-12.0 |
| Values excluded, 30% capped | 0.0 | 0.0 | 0.0 |

## Other cells

| | n = 200 | n = 1000 | n = 5000 |
|---|--:|--:|--:|
| NS, bias of guideline, no capping | 0.15 | 0.05 | 0.04 |
| NS, bias of guideline, 30% capped (mean) | 6.0 | 5.8 | 5.9 |
| NS, coverage of guideline, 30% capped | 3.5-26.0 | 0.0 | 0.0 |
| NS, values excluded, 30% capped | 18.1-22.5 | 18.1-22.5 | 18.1-22.5 |
| Near, coverage of period | 92.5-97.0 | 92.0-97.0 | 92.0-96.5 |
| Near, coverage of guideline, shift > 0 | 93.0-97.5 | 92.5-97.0 | 90.5-96.5 |
| Near, guideline shown, shift > 0 | 100.0 | 100.0 | 100.0 |
| Reported, compliant, coverage of guideline, shift > 0 | 93.0-96.0 | 93.5-95.5 | 90.0-98.0 |
| Reported, compliant, guideline shown, shift 0 | 0.0 | 0.0 | 0.0 |
| Reported, own rule, coverage of guideline | 94.0-96.5 | 95.0-96.0 | 90.5-96.0 |
| Reported, own rule, guideline shown | 3.0-5.0 | 3.0-4.0 | 0.0 |
| Clustered, coverage of period, isolate bootstrap | 85.5-93.0 | 52.5-65.5 | 33.0-39.5 |
| Clustered, coverage of period, site bootstrap | 91.5-96.5 | 91.5-97.0 | 93.5-95.5 |
| Clustered, period shown, no period change, isolate bootstrap | 7.0-8.5 | 34.5-36.0 | 66.0-66.5 |
| Clustered, period shown, no period change, site bootstrap | 3.5-5.5 | 3.0-8.5 | 6.0-6.5 |
| Clustered, coverage of guideline, shift 1, isolate bootstrap | 93.0-93.5 | 91.0-94.5 | 91.5-94.5 |
| Clustered, coverage of guideline, shift 1, site bootstrap | 92.0-93.0 | 89.0-92.0 | 93.5-94.0 |
| Clustered, larger named wrongly, isolate bootstrap | 0.0-8.5 | 0.0-34.5 | 0.0-66.0 |
| Clustered, larger named wrongly, site bootstrap | 0.0-5.5 | 0.0-3.0 | 0.0-6.5 |
| Low, coverage of period | 91.5-97.5 | 93.0-97.0 | 92.5-95.5 |
| Low, coverage of guideline, shift 1 | 93.5-98.0 | 92.0-94.5 | 91.0-97.5 |
| Low, period shown, no period change | 2.5-5.5 | 2.5-6.5 | 0.0-1.5 |
| Low, period shown, +1 point | 5.0-8.5 | 19.0-31.5 | 44.5-49.5 |
| Low, guideline shown, shift 0 | 0.0 | 0.0 | 0.0 |
| Low, guideline shown, shift 1 | 100.0 | 100.0 | 100.0 |
| Low, larger named wrongly | 0.0-5.5 | 0.0-6.5 | 0.0-1.0 |
| Straddle, not interpretable | 100.0 | 100.0 | 100.0 |

True components: period 4.4-4.8 points with a period change (base cells); guideline 4.5-5.0 points (shift 1) and 16.8-20.1 (shift 2); reported own-rule guideline -0.1-0.0.

Base cell with the lowest period coverage (n = 5000, prev 0.20, shift 0, delta 0.00, capped 0.0), rerun with 600 new datasets: coverage 95.8, period shown 3.8.

## Twenty pairs per dataset, Benjamini-Hochberg over them

| | % |
|---|--:|
| guideline shown, revised pairs | 100.0 |
| guideline shown, unrevised pairs | 0.0 |
| period shown, +5 points | 51.2 |
| period shown, no period change | 1.1 |
| false discovery rate, guideline | 0.0 |
| false discovery rate, period | 1.6 |
| pairs not interpretable | 0.0 |

## MIC-shift model

| Case | Target | Estimated | Mean shift | Coverage | p < 0.05 |
|---|--:|--:|--:|--:|--:|
| panel | 0.50 | 100.0 | 0.50 | 94.5 | 100.0 |
| bimodal | 0.00 | 100.0 | 0.27 | 37.0 | 63.0 |
