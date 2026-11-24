# Simulation

Written by scripts/simulation.R (one seed per cell from 2026, 1000 datasets per cell, B = 999, 186 cells, 416.8 minutes).
Setup in the script header. Bias in percentage points, the rest in percent. Ranges are over the
cells with that n; per-cell results, with the true components, are in simulation.csv.

## %R, base cells

| | n = 200 | n = 1000 | n = 5000 |
|---|--:|--:|--:|
| Bias, fixed criteria, largest absolute | 0.24 | 0.11 | 0.05 |
| Bias, reclassification, largest absolute | 0.11 | 0.04 | 0.02 |
| Coverage, fixed criteria, no capping | 93.7-96.2 | 93.4-95.4 | 93.0-95.8 |
| Coverage, fixed criteria, 30% capped | 93.9-95.6 | 93.9-96.4 | 94.1-96.1 |
| Coverage, reclassification, shift > 0, no capping | 93.4-96.2 | 93.5-94.8 | 93.8-95.6 |
| Coverage, reclassification, shift > 0, 30% capped | 93.3-95.3 | 93.2-96.0 | 94.2-95.5 |
| Fixed criteria shown, no change at fixed criteria | 3.2-5.3 | 3.2-6.3 | 3.3-5.9 |
| Fixed criteria shown, change at fixed criteria | 17.6-39.1 | 63.6-97.1 | 99.9-100.0 |
| Reclassification shown, shift 0 | 0.0 | 0.0 | 0.0 |
| Reclassification shown, shift > 0 | 100.0 | 100.0 | 100.0 |
| Larger named correctly, shift > 0 and change at fixed criteria | 1.0-97.7 | 2.5-100.0 | 4.1-100.0 |
| Larger named wrongly, any cell | 0.0-4.5 | 0.0-4.7 | 0.0-5.9 |
| Values excluded, 30% capped | 0.0 | 0.0 | 0.0 |

## Other cells

| | n = 200 | n = 1000 | n = 5000 |
|---|--:|--:|--:|
| NS, bias of reclassification, no capping | 0.08 | 0.04 | 0.01 |
| NS, bias of reclassification, 30% capped (mean) | 5.9 | 5.9 | 5.9 |
| NS, coverage of reclassification, 30% capped | 5.2-24.9 | 0.0 | 0.0 |
| NS, values excluded, 30% capped | 18.2-22.5 | 18.1-22.5 | 18.1-22.5 |
| Near, coverage of fixed criteria | 93.6-95.8 | 93.9-96.3 | 93.9-95.7 |
| Near, coverage of reclassification, shift > 0 | 93.3-96.0 | 94.4-95.6 | 93.4-96.5 |
| Near, reclassification shown, shift > 0 | 100.0 | 100.0 | 100.0 |
| Reported, compliant, coverage of reclassification, shift > 0 | 93.8-95.2 | 94.2-95.5 | 92.9-95.9 |
| Reported, compliant, reclassification shown, shift 0 | 0.0 | 0.0 | 0.0 |
| Reported, own rule, coverage of reclassification | 93.5-95.8 | 95.6-95.8 | 95.1-95.2 |
| Reported, own rule, reclassification shown | 3.7-4.4 | 2.8-3.2 | 0.0 |
| Clustered, coverage of fixed criteria, isolate bootstrap | 85.0-89.3 | 62.3-65.5 | 35.4-39.1 |
| Clustered, coverage of fixed criteria, site bootstrap | 92.9-94.9 | 93.5-94.4 | 93.1-93.8 |
| Clustered, fixed criteria shown, no change at fixed criteria, isolate bootstrap | 10.1 | 32.9-34.3 | 60.4-63.1 |
| Clustered, fixed criteria shown, no change at fixed criteria, site bootstrap | 4.8-5.0 | 5.9 | 5.6-5.7 |
| Clustered, coverage of reclassification, shift 1, isolate bootstrap | 94.2 | 93.7-94.9 | 92.2-94.3 |
| Clustered, coverage of reclassification, shift 1, site bootstrap | 93.1-94.1 | 92.2-93.5 | 92.0-94.0 |
| Clustered, larger named wrongly, isolate bootstrap | 0.0-10.1 | 0.0-34.3 | 0.0-63.1 |
| Clustered, larger named wrongly, site bootstrap | 0.0-4.8 | 0.0-5.9 | 0.0-5.6 |
| Low, coverage of fixed criteria | 92.7-95.9 | 94.1-95.9 | 93.4-95.6 |
| Low, coverage of reclassification, shift 1 | 92.7-94.6 | 92.9-95.0 | 94.5-95.0 |
| Low, fixed criteria shown, no change at fixed criteria | 3.6-5.3 | 3.7-5.1 | 0.0-1.4 |
| Low, fixed criteria shown, +1 point | 6.4-9.2 | 15.7-33.2 | 44.0-48.5 |
| Low, reclassification shown, shift 0 | 0.0 | 0.0 | 0.0 |
| Low, reclassification shown, shift 1 | 100.0 | 100.0 | 100.0 |
| Low, larger named wrongly | 0.0-5.2 | 0.0-5.1 | 0.0-0.6 |
| Straddle, not interpretable | 100.0 | 100.0 | 100.0 |

True components: fixed criteria 4.4-4.8 points with a change at fixed criteria (base cells); reclassification 4.5-5.0 points (shift 1) and 16.8-20.1 (shift 2); reported own-rule reclassification -0.1-0.0.

Base cell with the lowest fixed-criteria coverage (n = 5000, prev 0.05, shift 1, delta 0.05, capped 0.0), rerun with 600 new datasets: coverage 95.8, fixed criteria shown 100.0.

## Twenty pairs per dataset, Benjamini-Hochberg over them

| | % |
|---|--:|
| reclassification shown, revised pairs | 100.0 |
| reclassification shown, unrevised pairs | 0.0 |
| fixed criteria shown, +5 points | 52.5 |
| fixed criteria shown, no change at fixed criteria | 1.6 |
| false discovery rate, reclassification | 0.0 |
| false discovery rate, fixed criteria | 2.3 |
| pairs not interpretable | 0.0 |

## MIC-shift model

| Case | Target | Estimated | Mean shift | Coverage | p < 0.05 |
|---|--:|--:|--:|--:|--:|
| panel | 0.50 | 100.0 | 0.50 | 94.9 | 100.0 |
| bimodal | 0.00 | 100.0 | 0.26 | 39.6 | 60.4 |
