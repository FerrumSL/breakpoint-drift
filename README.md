# Breakpoint drift auditor

When CLSI or EUCAST lowers a breakpoint, isolates with the same MIC move from S to R, and %R
in the antibiogram jumps with no change in the bacteria. A 2025 review in *Frontiers in
Microbiology* calls this breakpoint drift and notes that there is no standard method for
separating it from real change in historical AST data
([doi:10.3389/fmicb.2025.1674968](https://doi.org/10.3389/fmicb.2025.1674968)). The effect is
not small: across 13,101 Enterobacterales and *P. aeruginosa* isolates over a decade, every
β-lactam and fluoroquinolone looked less active after CLSI updates
([Antimicrob Steward Healthc Epidemiol](https://www.cambridge.org/core/journals/antimicrobial-stewardship-and-healthcare-epidemiology/article/changing-times-the-impact-of-gramnegative-breakpoint-changes-over-the-previous-decade/7FD634795CF29C39B441B75125B9843D)).

This tool takes raw MICs or zone diameters from two periods and splits the change in %R for
each organism and drug into the part caused by the breakpoint revision and the part caused by
the values themselves. It then tests each pair, estimates how far the MIC distribution moved
without using any breakpoint, and writes a report for a stewardship committee. CLSI and EUCAST
are both supported. The breakpoint tables, and the reading of each value, come from the
[AMR package](https://amr-for-r.org).

## Method

Let p(y, g) be the share of isolates from period y that are R under guideline g, and g1, g2
the editions in force in the two periods.

- observed change = p(y2, g2) − p(y1, g1), what two antibiograms side by side show
- real change = p(y2, g2) − p(y1, g2): the same guideline, values from two periods
- artefact = p(y1, g2) − p(y1, g1): the same values, two guidelines

The split depends on which step comes first: holding g1 instead gives
real change = p(y2, g1) − p(y1, g1). The tool reports both, and uses their mean (the Shapley
value) as the real change, so real change and artefact always add up to the observed change.

**Capped MICs.** A stored "≤1" cannot be read against a new breakpoint of S ≤ 0.25. AMR returns
NI (not interpretable) for such values, they are left out of that p(y, g), and the report shows
how many isolates each guideline could read. The Dutch ISIS-AR study reported the same coverage
figure when it reinterpreted old *E. coli* MICs: 70-90% of isolates, depending on the drug
([Clin Microbiol Infect](https://www.sciencedirect.com/science/article/pii/S1198743X1460759X)).

**Verdicts.** Two Fisher's exact tests per pair, each with a Benjamini-Hochberg correction across
pairs: the observed comparison and the real one (both periods under g2). A pair is a
*significant change* if the real test is significant, a *breakpoint artefact* if only the
observed test is and the artefact is larger than the real change, and *noise* otherwise.

**MICs and zones, CLSI and EUCAST.** Each organism and drug gets one row per test method and
organisation. MIC and disk breakpoints are revised separately, and a lab moving from disk to
MIC is a change of its own, so the two are never pooled. An `org` column says which
organisation each row was read with.

**MIC shift.** %R only moves when MICs cross a breakpoint. `drift_mic_shift()` measures the
whole MIC distribution instead. A stored MIC is an interval: "8" means above 4 and at most 8,
"≤0.25" means at most 0.25. A normal model of log2 MIC with a period effect is fitted to those
intervals by maximum likelihood (`survival::survreg`), and the period effect is the shift in
doubling dilutions, with a 95% CI. It uses no breakpoint, so it picks up MIC creep while %R is
still flat. A period needs at least 10 values, 5 of them not capped; a period of only "≤1"
values has no estimable mean, and the shift is left blank.

## Use

```r
source("R/drift.R")
# one row per first isolate and drug: year, mo, ab, and mic ("<=0.25", "8", ">=16")
# or disk (mm); optionally org ("CLSI" or "EUCAST", default CLSI)
d <- drift_decompose(isolates, from = 2013:2016, to = 2021:2024)
t <- drift_test(d)
shift <- drift_mic_shift(isolates, from = 2013:2016, to = 2021:2024)
drift_report(merge(t, shift, by = c("mo", "ab"), all.x = TRUE), "report.md")
```

`guideline_in_force()` decides which edition was in force: the one named for the year, or
for the last year of a period, since that is when a cumulative antibiogram for it would have
been written. Labs often adopt a new edition late; `in_force = function(y, org) guideline_in_force(y, org, lag = 1)`
models a one-year lag. Any function of year and organisation that returns a guideline name works.

## Tests

`tests/testthat/test_drift.R` builds MICs and zones by hand, or draws them from a known
distribution, so the true answer is known. It uses three revisions in the AMR package's
tables: CLSI 2020 cut the *E. coli* ciprofloxacin MIC breakpoints from S ≤ 1 / R ≥ 4 to
S ≤ 0.25 / R ≥ 1, CLSI 2023 raised the gentamicin disk breakpoints (14 mm went from I to R),
and EUCAST 2017 cut ciprofloxacin (an MIC of 1 went from I to R). The tests check that:

- the same 100 MICs read under both editions give real change 0 and artefact +30 points
- a shifted MIC distribution under unchanged breakpoints gives artefact 0
- real change and artefact always sum to the observed change
- "≤1" values are excluded under the new breakpoint and counted
- zones and MICs for one drug are decomposed apart, and EUCAST rows get EUCAST editions
- verdicts separate an artefact, a real shift, and a three-isolate pair (noise)
- the censored model recovers a known one-dilution shift from capped panel values, and finds
  a shift while %R stays at 0 (creep below the breakpoint)
- it returns blank rather than a number when a period is all capped or too small

```sh
Rscript -e 'testthat::test_dir("tests/testthat")'
```

## Example on real data

[`results/ncbi_report.md`](results/ncbi_report.md) runs the tool on 21,555 MIC and disk tests
(CLSI and EUCAST) from 1,783 human isolates in NCBI BioSample, 2011–2018 against 2019–2025
(`Rscript scripts/ncbi.R`, about 90 seconds). It reports the 34 pairs with at least 10 readable
isolates in both periods. The clearest case is *E. coli* gentamicin: %R rose from 7.6% to 43.8%, but
31.5 of those 36.1 points come from CLSI's 2023 cut to the aminoglycoside breakpoints, and the
real change (+4.6) is not significant. Amikacin, revised in the same edition, shows the same
pattern, and there the MIC shift (+1.8 dilutions, 95% CI +0.4 to +3.1) says part of the rise
is real too.

This run shows the tool works on real, messy exports. It does not show real resistance trends:
submitters chose which isolates to sequence, often because they were resistant, and which labs
submitted changes from year to year. The 2019–2025 *E. coli* group has only 16 isolates. The
proof that the split is right is the tests above.

## Limits

- The AMR package holds CLSI and EUCAST editions from 2011 onward; earlier years have no
  breakpoints to read them with. Older editions are also sparser: AMR 3.1.0 has CLSI disk
  breakpoints for Enterobacterales only from the 2020 edition for most drugs, so older CLSI
  zones for them come out *not interpretable*.
- NCBI records the organisation but not the edition a lab used, so the edition is assumed from
  the collection year.
- The MIC shift fits one normal distribution per period. Many drugs have a wild-type and a
  resistant population, and then the size of the shift is overstated (two shifts in the
  NCBI run exceed 8 dilutions); its direction and significance are more reliable than its size.
  A mixture model with censoring
  ([PLOS ONE 2019](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0220427))
  would fix this. Zone diameters get no shift estimate.
- Deduplication is left to the caller: run `AMR::first_isolate()` (CLSI M39) before the
  decomposition.

## Data and credits

`data/raw/ncbi_ast_long.csv.gz` is the NCBI BioSample snapshot made by `R/00_fetch_ncbi.R` in
[amr-surveillance-r](https://github.com/FerrumSL/amr-surveillance-r), copied unchanged
(SHA-256 `a4ed673b…0b498`). The AMR package is by Berends et al.
([J Stat Softw 2022](https://doi.org/10.18637/jss.v104.i03)), GPL-2.

R 4.6.1; packages pinned in `renv.lock` (`renv::restore()`). Licence: GPL-2.
