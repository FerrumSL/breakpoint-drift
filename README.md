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

On 25,136 MICs from US national surveillance (NARMS, 2023) the tool's reading agrees with
NARMS's own calls in 99.5% of tests. When NARMS changed its own ciprofloxacin rules in 2023,
reported *Salmonella* resistance fell from 11.4% to 0.6%; the tool attributes −11.1 points to
the relabelling and +0.3 to the MICs ([below](#validation-on-us-national-surveillance-narms)).

The same problem, a change in measurement mistaken for a change in reality, comes up in
official statistics when an office changes its method; see
[series-break-backcast](https://github.com/FerrumSL/series-break-backcast).

## Method

Let p(y, g) be the share of isolates from period y that are R under guideline g, and g1, g2
the editions in force in the two periods.

- observed change = p(y2, g2) − p(y1, g1), what two antibiograms side by side show
- real change = p(y2, g2) − p(y1, g2): the same guideline, values from two periods
- artefact = p(y1, g2) − p(y1, g1): the same values, two guidelines

The split depends on which step comes first: holding g1 instead gives
real change = p(y2, g1) − p(y1, g1). The tool reports both, and uses their mean (the Shapley
value) as the real change, so real change and artefact always add up to the observed change.

**Reported calls.** If the data also holds the S/I/R each isolate was originally reported with
(a `reported` column), the observed change is the change in reported %R, as the published
antibiograms showed it. The artefact then also covers a lab's own reporting rules, not only
guideline editions.

**Capped MICs.** A stored "≤1" cannot be read against a new breakpoint of S ≤ 0.25. AMR returns
NI (not interpretable) for such values, they are left out of that p(y, g), and the report shows
how many isolates each guideline could read. The Dutch ISIS-AR study reported the same coverage
figure when it reinterpreted old *E. coli* MICs: 70-90% of isolates, depending on the drug
([Clin Microbiol Infect](https://www.sciencedirect.com/science/article/pii/S1198743X1460759X)).

**Missing breakpoints.** When an edition has no breakpoint for an organism's own group, AMR
falls back to a broader one. For Enterobacterales that is the generic Gram-negative row (CLSI
Table 2B-5, "other non-Enterobacterales"), which does not apply to them. In AMR 3.1.0 the CLSI
editions before 2020 lack Enterobacterales rows for many drugs, so 2015 *E. coli* ceftriaxone
would be read with S ≤ 8 instead of S ≤ 1. The tool drops those readings, and such pairs come
out *not interpretable* rather than as a false artefact.

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
still flat. The shift is left blank when it cannot be estimated: a period with fewer than 10
values or fewer than 5 uncapped ones, or a fit whose mean lies outside the range the panel
tested (mostly capped data let the fit run off).

## Use

```r
source("R/drift.R")
# one row per first isolate and drug: year, mo, ab, and mic ("<=0.25", "8", ">=16")
# or disk (mm); optionally org ("CLSI" or "EUCAST", default CLSI) and reported (S/I/R)
d <- drift_decompose(isolates, from = 2020:2021, to = 2023:2024)
t <- drift_test(d)
shift <- drift_mic_shift(isolates, from = 2020:2021, to = 2023:2024)
drift_report(merge(t, shift, by = c("mo", "ab"), all.x = TRUE), "report.md")
```

**Which edition was in force.** By default, the one named for the year, or for the last year of
a period, since that is when a cumulative antibiogram for it would have been written.
`in_force = function(y, org) guideline_in_force(y, org, lag = 1)` models a one-year lag. A lab
that knows when it adopted each edition should pass that instead, as a table:

```r
adoption <- read.csv(text = "
lab,year,guideline
north,2014,CLSI 2014
north,2021,CLSI 2020
south,2014,CLSI 2014
south,2020,EUCAST 2020")
d <- drift_decompose(isolates, from = 2020:2021, to = 2023:2024, adoption = adoption)
```

From `year` on, that lab read with `guideline`. `isolates` then needs a `lab` column; leave
`lab` out of the table if all sites changed together. Labs on different editions of the same
organisation are pooled, each isolate read with its own lab's edition, and `guideline_from` /
`guideline_to` list them all ("CLSI 2020 + CLSI 2023"). A lab that switched from CLSI to EUCAST
between the periods, like `south` here, is a source of drift in itself and gets its own row
("CLSI to EUCAST").

## Tests

`tests/testthat/test_drift.R` builds MICs and zones by hand, or draws them from a known
distribution, so the true answer is known. Most tests use CLSI's 2023 cut to the
Enterobacterales gentamicin breakpoints (MIC S ≤ 4 / R ≥ 16 to S ≤ 2 / R ≥ 8; disk 14 mm from I
to R), and one uses EUCAST's 2021 cut to piperacillin/tazobactam. The tests check that:

- the same 100 MICs read under both editions give real change 0 and artefact +30 points
- a shifted MIC distribution under unchanged breakpoints gives artefact 0
- real change and artefact always sum to the observed change
- "≤4" values are excluded under the new breakpoint and counted
- an *E. coli* reading from the generic Gram-negative fallback counts as not readable
- zones and MICs for one drug are decomposed apart, and EUCAST rows get EUCAST editions
- verdicts separate an artefact, a real shift, and a three-isolate pair (noise)
- reported calls replace the guideline reading as the observed change
- the censored model recovers a known one-dilution shift from capped panel values, and finds
  a shift while %R stays at 0 (creep below the breakpoint)
- it returns blank rather than a number when a period is all capped, too small, or fitted
  outside the tested range
- an adoption table overrides the year rule, pools labs on different editions correctly, and
  separates a lab that moved from CLSI to EUCAST

```sh
Rscript -e 'testthat::test_dir("tests/testthat")'
```

## Validation on US national surveillance (NARMS)

NARMS, the US National Antimicrobial Resistance Monitoring System, samples retail meat and
food-animal guts on a fixed plan and tests every isolate by broth microdilution in one lab
network. FDA publishes the isolate-level MICs, with NARMS's own S/I/R call, free and without
registration. `Rscript scripts/narms.R` (about 90 seconds) uses the *Salmonella* and
*E. coli* isolates: 229,543 isolates and 744,155 MIC tests from 2002 to 2023.

**Agreement.** Read with the CLSI edition of their year, 2020–2023 cecal isolates agree with
NARMS's calls in 99.2–99.5% of tests per year (25,011 of 25,136 in 2023,
[`narms_agreement.csv`](results/narms_agreement.csv)). Every 2023 disagreement is a rule NARMS
applies differently from CLSI 2023, none a misread value: gentamicin MICs of 4 and 8 (NARMS
still uses the pre-2023 breakpoint, S ≤ 4 / R ≥ 16) and *E. coli* ciprofloxacin MICs of 0.12
and 0.25 (NARMS calls them I, as for *Salmonella*; CLSI's Enterobacterales breakpoint makes
them S). The reading itself is a table lookup, so it can only be as wrong as the table.

**NARMS's 2023 relabelling** ([`narms_relabel.md`](results/narms_relabel.md)). Through 2022,
NARMS called a cecal *Salmonella* or *E. coli* ciprofloxacin MIC of 0.12 mg/L or more R; in 2023
it moved to CLSI, where 0.12–0.5 is I. With the `reported` calls as the observed change:

| 2022 to 2023, cecal | Reported %R | Observed | Real | Artefact | Verdict |
|---|--:|--:|--:|--:|---|
| *Salmonella*, ciprofloxacin | 11.4 to 0.6 | −10.8 | +0.3 | −11.1 | breakpoint artefact |
| *E. coli*, ciprofloxacin | 7.6 to 1.7 | −5.9 | +0.4 | −6.4 | breakpoint artefact |

Over the same year the table separates real rises from the relabelling: *Salmonella*
tetracycline +5.7 and nalidixic acid +5.5 points, both significant and all real.

**A revision NARMS has not applied yet** ([`narms_editions.md`](results/narms_editions.md)).
In 2023 NARMS still reported cecal gentamicin with the old CLSI breakpoint (an MIC of 8 is I),
while its 2024 refresh of the retail file applies the new one backwards (8 is R). Reading
2019–2020 and 2022–2023 cecal isolates with CLSI 2020 and CLSI 2023, *Salmonella* gentamicin
%R rises from 1.7% to 5.8%: +3.4 points artefact, +0.7 real, and the MIC distribution moved
−0.2 dilutions (95% CI −0.3 to −0.1). Read with CLSI 2023 instead of CLSI 2020, the same
2022–2023 isolates go from 2.3% to 5.8% R: once NARMS adopts the 2023 breakpoint, its
*Salmonella* gentamicin resistance will more than double with no change in the bacteria.

**Twenty years of MIC shift** ([`narms_mic_shift.csv`](results/narms_mic_shift.csv)). On retail
isolates, 2002–2006 against 2017–2021, the shift model needs no breakpoint, so the AMR
package's gaps before 2020 do not matter. *Salmonella* nalidixic acid MICs rose 0.6 dilutions
(0.5 to 0.7) and *E. coli* nalidixic acid fell 0.5; the other eight estimated shifts are under
0.25 dilutions. Ten of the 20 pairs tested in both periods were left blank, because most of
their values sat at the panel's lowest dilution and the fit fell outside the tested range.

NARMS isolates come from food animals and meat, not patients, and the breakpoints are human
clinical ones; that is fine for testing the method, which does not depend on where an isolate
came from.

## Example on NCBI isolates

[`results/ncbi_report.md`](results/ncbi_report.md) runs the tool on 16,832 MIC and disk tests
(CLSI and EUCAST) from 1,496 human isolates in NCBI BioSample, 2016–2020 against 2021–2025
(`Rscript scripts/ncbi.R`, about a minute). It reports the 47 pairs with at least 10 readable
isolates in both periods: 7 significant changes, no breakpoint artefacts.

This run shows the tool works on real, messy exports of human isolates. It does not show real
resistance trends: submitters chose which isolates to sequence, often because they were
resistant, and which labs submitted changes from year to year.

## Limits

- The AMR package holds CLSI and EUCAST editions from 2011 onward, but its CLSI editions before
  2020 lack Enterobacterales rows for many drugs, MIC and disk (see *Missing breakpoints*). For
  Enterobacterales, comparisons that need a pre-2020 CLSI edition mostly come out not
  interpretable. Adding the missing rows from the CLSI M100 editions would fix this; the tool
  passes the AMR package's tables through unchanged.
- NCBI and NARMS record the organisation but not the edition a lab used, so both runs assume
  the edition from the collection year. With a lab's own adoption dates this guess is not
  needed.
- The MIC shift fits one normal distribution per period. Many drugs have a wild-type and a
  resistant population, and then the size of the shift is overstated; its direction and
  significance are more reliable than its size. A mixture model with censoring
  ([PLOS ONE 2019](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0220427))
  would fix this. With tens of thousands of isolates, a shift of 0.1 dilution is significant
  without mattering. Zone diameters get no shift estimate.
- Deduplication is left to the caller: run `AMR::first_isolate()` (CLSI M39) before the
  decomposition.

## Data and credits

- `data/raw/ncbi_ast_long.csv.gz` is the NCBI BioSample snapshot made by `R/00_fetch_ncbi.R` in
  [amr-surveillance-r](https://github.com/FerrumSL/amr-surveillance-r), copied unchanged
  (SHA-256 `a4ed673b…0b498`).
- `data/raw/narms_extract.csv.gz` holds the *Salmonella* and *E. coli* columns of two FDA NARMS
  files, *Retail Meats* and *Cecal* (updated 27 September 2024), from the
  [NARMS download page](https://www.fda.gov/animal-veterinary/national-antimicrobial-resistance-monitoring-system/integrated-reportssummaries).
  `scripts/narms.R` lists their SHA-256 and rebuilds the extract if they are put in
  `data/narms/`. Source: Food and Drug Administration (FDA), NARMS Now, Rockville, MD:
  U.S. Department of Health and Human Services.
- The AMR package is by Berends et al.
  ([J Stat Softw 2022](https://doi.org/10.18637/jss.v104.i03)), GPL-2.

R 4.6.1; packages pinned in `renv.lock` (`renv::restore()`). Licence: GPL-2.
