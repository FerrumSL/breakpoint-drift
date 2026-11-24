# breakpointdrift

[![test](https://github.com/FerrumSL/breakpoint-drift/actions/workflows/test.yml/badge.svg)](https://github.com/FerrumSL/breakpoint-drift/actions/workflows/test.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23770271.svg)](https://doi.org/10.5281/zenodo.23770271)
[![Project Status: Active](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License: GPL-2](https://img.shields.io/badge/license-GPL--2-blue.svg)](LICENSE)

`breakpointdrift` shows how much of a change in an antibiogram comes from new breakpoints and how
much from the isolates themselves.

When CLSI or EUCAST lowers a breakpoint, an isolate with the same MIC can move from S to R. The
share of resistant isolates (%R) then rises with no change in the bacteria. A 2025 review calls
this *breakpoint drift* and notes that there is no standard method for separating it from other
change in historical susceptibility data
([Sallam 2025](https://doi.org/10.3389/fmicb.2025.1674968)). The same happens when a laboratory
or a surveillance programme changes its own reporting rule.

The package reads raw MICs or zone diameters from two periods under the criteria in force in
each, and splits the change in %R for every organism and drug into two parts:

- **reclassification**: the part that comes from reading the two periods with different
  criteria;
- **change at fixed criteria**: the change when both periods are read with the same criteria.

The two parts add up to the change the antibiogram shows.

## What it does

- Reads raw MICs, including capped values such as `<=0.25`, and zone diameters, with CLSI or
  EUCAST breakpoints from the [AMR package](https://amr-for-r.org).
- Splits the change in %R, or in the share non-susceptible (%NS), for each organism and drug.
- Gives bootstrap intervals and tests, corrected for the number of pairs, and a verdict for each
  pair. Resamples whole clusters, such as sampling sites, when isolates are not independent.
- Uses each laboratory's adoption dates and its originally reported S/I/R calls when you have
  them.
- Says what it could not read: missing breakpoints, and capped values that fall on both sides of
  a breakpoint.
- Writes a Markdown report for a stewardship committee.

## Installation

```r
# install.packages("remotes")
remotes::install_github("FerrumSL/breakpoint-drift")
```

It installs `AMR` (3.1.0 or later) and `survival`. The tests run on R 4.6.

## Quick start

Gentamicin MICs of *E. coli* from 2022 and 2023. CLSI lowered the gentamicin breakpoints for
Enterobacterales in 2023 (S ≤ 4 / R ≥ 16 became S ≤ 2 / R ≥ 8).

```r
library(breakpointdrift)

values <- c("<=1", "2", "4", "8", "16", ">=32")
isolates <- data.frame(
  year = rep(c(2022, 2023), each = 200),
  mo = "E. coli",
  ab = "gentamicin",
  mic = c(rep(values, c(110, 40, 20, 14, 10, 6)),   # 2022
          rep(values, c(100, 42, 22, 18, 11, 7))),  # 2023
  isolate = 1:400
)

set.seed(1)
result <- drift_test(drift_bootstrap(isolates, B = 2000))
result[c("observed", "fixed", "reclass", "verdict")]
#>   observed fixed reclass          verdict
#> 1      0.1  0.02    0.08 reclassification
```

Read with the edition in force each year, %R rose from 8% to 18%, a change of 10 percentage
points. Reclassification accounts for 8 of them (95% interval 5.2 to 10.8). The change at fixed
criteria is 2 points (−3.8 to 7.8), which is not significant.

`drift_report(result, "report.md")` writes the same result as a table for each organism:

| Drug | Test | %R 2022 | %R 2023 | Observed | Fixed criteria | Reclassification | Verdict | Larger | Read 2022 |
|---|---|--:|--:|--:|--:|--:|---|---|--:|
| Gentamicin | MIC CLSI | 8.0 | 18.0 | +10.0 | +2.0 | +8.0 | reclassification | neither | 200 / 200 of 200 |

## Your own data

One row per isolate and drug, after removing duplicate isolates (for example with
`AMR::first_isolate()`).

| Column | Required | Content |
|---|---|---|
| `year` | yes | Year the isolate was collected |
| `mo` | yes | Organism, any name `AMR::as.mo()` reads |
| `ab` | yes | Drug, any name `AMR::as.ab()` reads |
| `mic` or `disk` | one per row | Raw MIC as stored (`"<=0.25"`, `"8"`, `">=16"`), or zone diameter in mm |
| `isolate` | no | Isolate ID: keeps an isolate's drugs together in the bootstrap |
| `reported` | no | S/I/R as originally reported: the observed change is then the reported one |
| `org` | no | `"CLSI"` (default) or `"EUCAST"` |

Periods can be several years: `drift_bootstrap(isolates, from = 2019:2020, to = 2022:2023)`. A
period is read with the edition in force in its last year. By default that is the edition named
for the year. If you know when each laboratory adopted each edition, pass a table instead:

```r
adoption <- data.frame(lab = c("north", "north", "south"),
                       year = c(2014, 2021, 2020),
                       guideline = c("CLSI 2014", "CLSI 2020", "EUCAST 2020"))
b <- drift_bootstrap(isolates, from = 2019:2020, to = 2022:2023, adoption = adoption)
```

`isolates` then needs a `lab` column. With a column of sampling sites, `cluster = "site"`
resamples whole sites instead of isolates.

## How it works

Let p(y, g) be the share resistant among isolates from period y read with criteria g. The two
shares an antibiogram reports are p(y₁, g₁) and p(y₂, g₂).

| | Older criteria g₁ | Newer criteria g₂ |
|---|---|---|
| **First period y₁** | p(y₁, g₁), reported | p(y₁, g₂) |
| **Second period y₂** | p(y₂, g₁) | p(y₂, g₂), reported |

The change at fixed criteria compares the two periods within a column. It is the mean of the two
columns (the Shapley value), so the result does not depend on which criteria are held fixed.
Reclassification is the rest of the observed change. This is the two-factor decomposition of
Kitagawa (1955) applied to a change of interpretive criteria.

A component is *detected* when its bootstrap p-value, corrected by Benjamini-Hochberg across the
pairs, is below 0.05 and the component is at least 1 percentage point. The verdict names the
components detected:

| Verdict | Meaning |
|---|---|
| reclassification | Only reclassification is detected |
| change at fixed criteria | Only the change at fixed criteria is detected |
| both | Both are detected |
| unattributed | The observed change is significant, but neither part is detected on its own |
| no change detected | No evidence of change after correction (not evidence of no change) |
| not interpretable | No isolate could be read under one edition, or capped values leave a part unknown over 1 point or more; `verdict_note` says which |

[`docs/method.md`](docs/method.md) explains how capped values and missing breakpoints are
handled, how the tests work, the breakpoint-free MIC shift, and what the automated tests check.

## Validation

- Read with the CLSI edition of their year, 27,060 MICs from US national surveillance (NARMS,
  2023) agree with NARMS's own S/I/R calls in 99.5% of tests.
- In 2023 NARMS changed how it reports ciprofloxacin. Reported resistance in *Salmonella* from
  cecal samples fell from 11.4% to 0.6%. Reclassification accounts for −11.1 points (95%
  interval −12.7 to −9.5), the change at fixed criteria for +0.3 (−0.2 to +0.8).
- In simulated data with independent isolates, both parts were close to unbiased, and 95%
  intervals covered the true value in 93.0–96.4% of datasets in every base scenario (1,000
  datasets each). With clustered isolates the cluster bootstrap was needed.

Every number comes from a script in `scripts/` and is written to `results/`; see
[`docs/method.md`](docs/method.md#validation-on-us-national-surveillance-narms).

## Limitations

- The method needs raw MICs or zone diameters. Results stored only as S/I/R cannot be read again.
- The AMR package's CLSI editions before 2020 lack Enterobacterales rows for many drugs. Such
  comparisons come out *not interpretable*. You can supply the missing rows with the
  `breakpoints` argument.
- The change at fixed criteria is not a causal estimate of change in the bacteria: sampling,
  test panels and case mix can change between periods too.
- The default assumes a laboratory adopts each edition in the year it is published. Many adopt
  later; pass an adoption table when you know the dates.

## Reproducing the results

```sh
Rscript -e 'renv::restore()'           # packages as pinned in renv.lock
Rscript scripts/narms.R                # NARMS validation, about 10 minutes
Rscript scripts/ncbi_fluoroquinolones.R  # CLSI's 2019 fluoroquinolone revision in NCBI data
Rscript scripts/simulation.R           # simulation, about 2 hours on 16 cores
Rscript scripts/whonet_audit.R         # audit of the historical CLSI tables
```

The scripts `source("R/drift.R")`, so they run from a clone without installing the package. The
tests run with `Rscript -e 'testthat::test_local()'`.

## Related software

- The [AMR package](https://amr-for-r.org) interprets MICs and zones and holds the breakpoint
  tables this package reads. `breakpointdrift` adds the decomposition of trends on top.
- [WHONET](https://whonet.org) and its interpretation engine
  [AMRIE](https://github.com/WHONET-Development/AMRIE) maintain the breakpoint file the AMR
  package imports.
- [series-break-backcast](https://github.com/FerrumSL/series-break-backcast) deals with the same
  problem in official statistics, where a change of method looks like a change in the economy.

## Citation

Cite the software with its DOI, [10.5281/zenodo.23770271](https://doi.org/10.5281/zenodo.23770271),
which always points to the latest version. GitHub's *Cite this repository* button gives the full
reference from `CITATION.cff`.

## Licence and credits

GPL-2. The breakpoint tables and their reading come from the AMR package
([Berends et al., J Stat Softw 2022](https://doi.org/10.18637/jss.v104.i03), GPL-2). NARMS data:
Food and Drug Administration, NARMS Now, U.S. Department of Health and Human Services. NCBI data:
BioSample antibiograms, National Center for Biotechnology Information. Questions and bug reports:
[GitHub issues](https://github.com/FerrumSL/breakpoint-drift/issues).
