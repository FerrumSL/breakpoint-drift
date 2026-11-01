# Breakpoint drift auditor

When CLSI or EUCAST lowers a breakpoint, isolates with the same MIC move from S to R, and %R
in the antibiogram jumps with no change in the bacteria. A 2025 review in *Frontiers in
Microbiology* calls this breakpoint drift and notes that there is no standard method for
separating it from other change in historical AST data
([doi:10.3389/fmicb.2025.1674968](https://doi.org/10.3389/fmicb.2025.1674968)). The effect is
not small: across 13,101 Enterobacterales and *P. aeruginosa* isolates over a decade, every
β-lactam and fluoroquinolone looked less active after CLSI updates
([Antimicrob Steward Healthc Epidemiol](https://www.cambridge.org/core/journals/antimicrobial-stewardship-and-healthcare-epidemiology/article/changing-times-the-impact-of-gramnegative-breakpoint-changes-over-the-previous-decade/7FD634795CF29C39B441B75125B9843D)).

This tool takes raw MICs or zone diameters from two periods and splits the change in %R for
each organism and drug into a guideline component, the part due to reading the two periods with
different interpretive criteria, and a period component, the change with the guideline held
fixed. The period component is not a causal estimate: sampling, panels and case mix can change
too. It then tests each pair, estimates how far the MIC distribution moved
without using any breakpoint, and writes a report for a stewardship committee. CLSI and EUCAST
are both supported. The breakpoint tables, and the reading of each value, come from the
[AMR package](https://amr-for-r.org).

On 27,060 MICs from US national surveillance (NARMS, 2023) the tool's reading agrees with
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
- period component = p(y2, g2) − p(y1, g2): the same guideline, values from two periods
- guideline component = p(y1, g2) − p(y1, g1): the same values, two guidelines

The split depends on which step comes first: holding g1 instead gives
period component = p(y2, g1) − p(y1, g1). The tool reports both (`real_at_new`, `real_at_old`)
and their difference (`interaction`), and uses their mean (the Shapley value) as the period
component, so the two components always add up to the observed change. This is the symmetric
two-factor decomposition of Kitagawa and Das Gupta applied to reinterpretation; the split of the
interaction is a convention. The columns are `real` (period) and `artefact` (guideline). The paper calls these the
*fixed-criteria* and the *criteria* component, and its verdicts "criteria change" and
"fixed-criteria change" are the software's "guideline change" and "period change"; "no change
detected" is "noise" and "change not attributed" is "unresolved".

**Reported calls.** If the data also holds the S/I/R each isolate was originally reported with
(a `reported` column), the observed change is the change in reported %R, as the published
antibiograms showed it. The guideline component then also covers a lab's own reporting rules and
any disagreement between the lab's reading and the tool's, not only guideline editions.

**Capped MICs.** A stored value is an interval of possible MICs: "≤4" is anything up to 4. The
tool reads both ends of the interval and counts the value if they fall on the same side of the
outcome's boundary: "≤4" is not R under R ≥ 8, although AMR calls it NI because it may be S or
I. A value that straddles the boundary ("≤1" against R ≥ 0.5, or "≤4" for %NS under S ≤ 2) is
left out of that p(y, g), and the report shows how many isolates each guideline could read.
A value that straddles the boundary under one guideline only ("≤8" against R ≥ 8 but not
R ≥ 16) leaves the two shares of a period on different isolates. The tool then gives the range
each component can lie in (`real_min`, `real_max`, `artefact_min`, `artefact_max`), counting
every such value as counted and as not, and a pair whose range is `min_effect` or wider is *not
interpretable*.

**Missing breakpoints.** When an edition has no breakpoint for an organism's own group, AMR
falls back to a broader one. For Enterobacterales that is the generic Gram-negative row (CLSI
Table 2B-5, "other non-Enterobacterales"), which does not apply to them. In AMR 3.1.0 the CLSI
editions before 2020 lack Enterobacterales rows for many drugs, so 2015 *E. coli* ceftriaxone
would be read with S ≤ 8 instead of S ≤ 1. The tool drops those readings, and such pairs come
out *not interpretable* rather than with a false guideline component.

**Tests and verdicts.** `drift_bootstrap()` resamples isolates within each period (or whole
clusters, such as sampling sites, with `cluster =`, drawing G − 1 of the G clusters each time;
at least 5 clusters per period, with a warning below 10)
and gives percentile intervals and bootstrap p-values for the observed change and both
components. `drift_test()` corrects them (Benjamini-Hochberg over the interpretable pairs) and
adds Fisher's exact test of the observed change. A component is *shown* when its corrected
p-value is below 0.05 and it is at least 1 percentage point (`min_effect`). The verdict is
*both*, *guideline change* or *period change* by which components are shown, *unresolved* if
neither is but the observed change is and both components are non-zero, and *noise* otherwise.
Column `larger` names the larger component when the bootstrap test of |guideline| − |period| is
significant and that component is shown. Neither label depends on which test is checked first,
and `verdict_note` says why a pair is *not interpretable*.

**%R or %NS.** `outcome = "NS"` counts non-susceptible isolates (I, SDD or R) instead of
resistant ones. A revision that only moves the I/R boundary then leaves no guideline component.

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
# install.packages("remotes")
remotes::install_github("FerrumSL/breakpoint-drift")
library(breakpointdrift)

# one row per first isolate and drug: year, mo, ab, and mic ("<=0.25", "8", ">=16")
# or disk (mm); optionally org ("CLSI" or "EUCAST", default CLSI) and reported (S/I/R)
set.seed(1)
b <- drift_bootstrap(isolates, from = 2020:2021, to = 2023:2024)  # decomposition, intervals, tests
t <- drift_test(b)
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

- the same 100 MICs read under both editions give period component 0 and guideline +30 points
- a shifted MIC distribution under unchanged breakpoints gives guideline component 0
- the two components always sum to the observed change
- "≤4" counts as not R under R ≥ 8, is left out only for %NS, and a mostly capped panel gives no
  false guideline component; ">4" is counted under R ≥ 8 only; values straddling one guideline's
  boundary only ("≤8", ">4" against R ≥ 16) make the pair *not interpretable*, one such value
  among 300 does not
- an *E. coli* reading from the generic Gram-negative fallback counts as not readable
- zones and MICs for one drug are decomposed apart, and EUCAST rows get EUCAST editions
- verdicts separate a guideline change, a period change, both, and a three-isolate pair (noise),
  and a guideline component four times the period one is named larger with 1,000, 5,000 and 20,000 isolates per period
- reported calls replace the guideline reading as the observed change, and a lab that followed
  each year's edition has the revision as its whole guideline component
- %NS ignores a revision that only moves the I/R boundary
- bootstrap intervals cover the estimate and the true values and are reproducible with a seed;
  the bootstrap p-values follow the Shapley components and are 1 for a component that is 0 in
  every replicate; the cluster bootstrap reports its clusters, refuses fewer than 5 per period and
  warns below 10
- a significant observed change that neither component shows is *unresolved*, and a component
  under one point does not count; with one component exactly 0 the pair is never *unresolved*,
  and `larger` never names a component the verdict does not show
- pairs that are not interpretable stay out of the corrections
- the censored model recovers a known one-dilution shift from capped panel values, finds a shift
  while %R stays at 0, and returns blank, with the reason, when a period is all capped, too
  small, or fitted outside the tested range
- an adoption table overrides the year rule, pools labs on different editions correctly, and
  separates a lab that moved from CLSI to EUCAST

```sh
R CMD build . && R CMD check --no-manual breakpointdrift_1.1.0.tar.gz
```

## Validation on US national surveillance (NARMS)

NARMS, the US National Antimicrobial Resistance Monitoring System, samples retail meat and
food-animal guts on a fixed plan and tests every isolate by broth microdilution in one lab
network. FDA publishes the isolate-level MICs, with NARMS's own S/I/R call, free and without
registration. `Rscript scripts/narms.R` (about 10 minutes, most of it the bootstrap) uses the *Salmonella*
and *E. coli* records: 229,543 sample records, of which 163,241 had no growth, leaving 65,273
isolates with an MIC and 770,135 MIC tests from 2002 to 2023. NARMS's 2023 calls the tool cannot
read (azithromycin, no CLSI Enterobacterales breakpoint) are counted in `narms_agreement.csv`.

**Agreement.** Read with the CLSI edition of their year, 2020–2023 cecal isolates agree with
NARMS's calls in 99.2–99.5% of tests per year (26,925 of 27,060 in 2023,
[`narms_agreement.csv`](results/narms_agreement.csv)). Every 2023 disagreement is a rule NARMS
applies differently from CLSI 2023, none a misread value: gentamicin MICs of 4 and 8 (NARMS
still uses the pre-2023 breakpoint, S ≤ 4 / R ≥ 16) and *E. coli* ciprofloxacin MICs of 0.12
and 0.25 (NARMS calls them I, as for *Salmonella*; CLSI's Enterobacterales breakpoint makes
them S). The reading itself is a table lookup, so it can only be as wrong as the table.

**NARMS's 2023 relabelling** ([`narms_relabel.md`](results/narms_relabel.md)). Through 2022,
NARMS called a cecal *Salmonella* or *E. coli* ciprofloxacin MIC of 0.12 mg/L or more R; in 2023
it moved to CLSI, where 0.12–0.5 is I. With the `reported` calls as the observed change:

| 2022 to 2023, cecal | Reported %R | Observed | Period | Guideline | Verdict |
|---|--:|--:|--:|--:|---|
| *Salmonella*, ciprofloxacin | 11.4 to 0.6 | −10.8 | +0.3 (−0.2 to +0.8) | −11.1 (−12.7 to −9.5) | guideline change |
| *E. coli*, ciprofloxacin | 7.6 to 1.8 | −5.8 | +0.6 (−0.4 to +1.6) | −6.4 (−7.8 to −5.1) | guideline change |

Bootstrap 95% intervals, 5,000 replicates over isolates. The cecal file has no sampling site;
the sensitivity analyses below resample whole collection months instead.

Over the same year the table separates period changes from the relabelling: *Salmonella*
tetracycline +5.7 and nalidixic acid +5.5 points, both period changes with no guideline component.

The relabelling also hid a rise. Counting non-susceptible isolates (I or R) instead
([`narms_relabel_ns.csv`](results/narms_relabel_ns.csv)), *Salmonella* ciprofloxacin went from
11.4% to 17.0% between 2022 and 2023: +5.6 points (2.9 to 8.2), all period component. Reported
%R fell to 0.6% in the same year
([`fig_cip_salmonella.pdf`](results/fig_cip_salmonella.pdf),
[`narms_cip_series.csv`](results/narms_cip_series.csv)).

**Checks** ([`narms_ecoff.csv`](results/narms_ecoff.csv), [`narms_strata.csv`](results/narms_strata.csv),
[`narms_flow.csv`](results/narms_flow.csv), [`narms_drugs.csv`](results/narms_drugs.csv)). Read
against the EUCAST wild-type cut-off (0.064 mg/L), which no breakpoint revision moves,
non-wild-type *Salmonella* ciprofloxacin also rose from 11.4% to 17.0%. Part of that rise is
composition: standardised to the pooled mix of host species it is +2.1 points, to the pooled
mix of serotypes +2.8, and to host and serotype jointly +0.8 (−0.6 to +2.3); within *S.* Infantis, mostly from chickens, it rose from 74% to
89%. The flow table counts records, de-duplication and readable results by file and year.

**Sensitivity analyses** ([`narms_cluster.csv`](results/narms_cluster.csv),
[`narms_shift_state.csv`](results/narms_shift_state.csv),
[`narms_mic_model.csv`](results/narms_mic_model.csv), [`narms_panels.csv`](results/narms_panels.csv)).
The cecal file records no sampling site and one isolate per organism per sample, so the cecal
intervals are re-run resampling whole collection months; the retail MIC shift is refitted with
a bootstrap over the 10 states sampled in both periods. The MIC-shift model is refitted with a separate standard
deviation per period and checked against the observed distribution, and the panel table lists
every drug whose tested range changed between years. `results/sessionInfo.txt` records the R
session of the last run.

**Breakpoint-table audit** (`Rscript scripts/whonet_audit.R`; [`whonet_audit.csv`](results/whonet_audit.csv),
[`whonet_example.csv`](results/whonet_example.csv)). Counts the CLSI Table 2A (Enterobacterales)
rows per edition 2011-2020 in WHONET's `Breakpoints.txt`, pinned to AMRIE commit 69de201 and
checked by SHA-256, and in the AMR package, and records which breakpoint row `as.sir()` uses for
an *E. coli* MIC of 1 mg/L under CLSI 2019 and 2020.

**A revision NARMS has not applied yet** ([`narms_editions.md`](results/narms_editions.md)).
In 2023 NARMS still reported cecal gentamicin with the old CLSI breakpoint (an MIC of 8 is I),
while its 2024 refresh of the retail file applies the new one backwards (8 is R). Reading
2019–2020 and 2022–2023 cecal isolates with CLSI 2020 and CLSI 2023, *Salmonella* gentamicin
%R rises from 1.7% to 5.8%: +3.4 points guideline component, +0.7 period component, and the MIC distribution moved
−0.2 dilutions (95% CI −0.3 to −0.1). The bulk of the distribution moved down (MICs ≤ 0.25
went from 27.8% to 38.9%) while the tail at 8 or more grew slightly, and %R follows the tail. Read with CLSI 2023 instead of CLSI 2020, the same
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
(`Rscript scripts/ncbi.R`, about a minute). It reports the 50 pairs with at least 10 readable
isolates in both periods: 9 period changes, 1 pair with both components, no guideline change alone.

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
  resistant population; in the simulation, a rise in the resistant share alone, with the wild
  type unchanged, showed as a shift of about 0.3 dilutions, significant in most datasets. Read
  it as a summary of the bulk of the distribution, not as a breakpoint-free measure of change. A mixture model with censoring
  ([PLOS ONE 2020](https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0220427))
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

R 4.6.2 (CI checks on R 4.6.1); packages pinned in `renv.lock` (`renv::restore()`). The scripts `source("R/drift.R")`
so they run from a clone without installing the package. Licence: GPL-2.
