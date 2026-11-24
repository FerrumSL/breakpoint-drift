# How breakpointdrift works

This page holds the detail behind the [README](../README.md): the decomposition, how values are
read, the tests, what the automated tests check, and the validation runs.

## Decomposition

Let p(y, g) be the share of isolates from period y that are R under criteria g, and g₁, g₂ the
editions in force in the two periods.

- observed change = p(y₂, g₂) − p(y₁, g₁), what two antibiograms side by side show
- change at fixed criteria = p(y₂, g₂) − p(y₁, g₂): the same criteria, values from two periods
- reclassification = p(y₁, g₂) − p(y₁, g₁): the same values, two sets of criteria

The split depends on which step comes first: holding g₁ instead gives a change at fixed criteria
of p(y₂, g₁) − p(y₁, g₁). The package reports both orderings (`fixed_at_new`, `fixed_at_old`) and
their difference (`interaction`), and uses their mean (the Shapley value) as the change at fixed
criteria, so the two parts always add up to the observed change. This is the symmetric
two-factor decomposition of Kitagawa and Das Gupta applied to reinterpretation; the equal split
of the interaction is a convention. The columns are `fixed` and `reclass`.

**Reported calls.** If the data also hold the S/I/R each isolate was originally reported with (a
`reported` column), the observed change is the change in reported %R, as the published
antibiograms showed it. Reclassification then also covers a laboratory's own reporting rules and
any disagreement between the laboratory's reading and the package's, not only guideline
editions.

**%R or %NS.** `outcome = "NS"` counts non-susceptible isolates (I, SDD or R) instead of
resistant ones. A revision that only moves the I/R boundary then causes no reclassification.

**MICs and zones, CLSI and EUCAST.** Each organism and drug gets one row per test method and
organisation. MIC and disk breakpoints are revised separately, and a laboratory moving from disk
to MIC is a change of its own, so the two are never pooled. An `org` column says which
organisation each row was read with.

## Reading the values

**Capped MICs.** A stored value is an interval of possible MICs: "≤4" is anything up to 4. The
package reads both ends of the interval and counts the value if they fall on the same side of the
outcome's boundary: "≤4" is not R under R ≥ 8, although the AMR package calls it NI because it
may be S or I. A value that straddles the boundary ("≤1" against R ≥ 0.5, or "≤4" for %NS under
S ≤ 2) is left out of that p(y, g), and the report shows how many isolates each edition could
read. A value that straddles the boundary under one edition only ("≤8" against R ≥ 8 but not
R ≥ 16) leaves the two shares of a period on different isolates. The package then gives the
range each part can lie in (`fixed_min`, `fixed_max`, `reclass_min`, `reclass_max`), counting
every such value as counted and as not, and a pair whose range is `min_effect` or wider is *not
interpretable*.

**Missing breakpoints.** When an edition has no breakpoint for an organism's own group, the AMR
package falls back to a broader one. For Enterobacterales that is the generic Gram-negative row
(CLSI Table 2B-5, "other non-Enterobacterales"), which does not apply to them. In AMR 3.1.0 the
CLSI editions before 2020 lack Enterobacterales rows for many drugs, so 2015 *E. coli*
ceftriaxone would be read with S ≤ 8 instead of S ≤ 1. The package drops those readings, and
such pairs come out *not interpretable* rather than with a false reclassification. Missing MIC
rows can be supplied with the `breakpoints` argument.

**Which edition was in force.** By default, the one named for the year, or for the last year of a
period, since that is when a cumulative antibiogram for it would have been written.
`in_force = function(y, org) guideline_in_force(y, org, lag = 1)` models a one-year lag. A
laboratory that knows when it adopted each edition should pass an adoption table instead. From
`year` on, that laboratory read with `guideline`. Laboratories on different editions of the same
organisation are pooled, each isolate read with its own laboratory's edition, and
`guideline_from` / `guideline_to` list them all ("CLSI 2020 + CLSI 2023"). A laboratory that
switched from CLSI to EUCAST between the periods is a source of drift in itself and gets its own
row ("CLSI to EUCAST").

## Tests and verdicts

`drift_bootstrap()` resamples isolates within each period (or whole clusters, such as sampling
sites, with `cluster =`, drawing G − 1 of the G clusters each time; at least 5 clusters per
period, with a warning below 10) and gives percentile intervals and bootstrap p-values for the
observed change and both parts. `drift_test()` corrects them (Benjamini-Hochberg over the
interpretable pairs) and adds Fisher's exact test of the observed change. A part is *detected*
when its corrected p-value is below 0.05 and it is at least 1 percentage point (`min_effect`).
The verdict is *both*, *reclassification* or *change at fixed criteria* by which parts are
detected, *unattributed* if neither is but the observed change is and both parts are non-zero,
and *no change detected* otherwise. Column `larger` names the larger part when the bootstrap test
of |reclass| − |fixed| is significant and that part is detected. Neither label depends on which
test is checked first, and `verdict_note` says why a pair is *not interpretable*.

## MIC shift

%R only moves when MICs cross a breakpoint. `drift_mic_shift()` measures the whole MIC
distribution instead. A stored MIC is an interval: "8" means above 4 and at most 8, "≤0.25" means
at most 0.25. A normal model of log2 MIC with a period effect is fitted to those intervals by
maximum likelihood (`survival::survreg`), and the period effect is the shift in doubling
dilutions, with a 95% CI. It uses no breakpoint, so it picks up MIC creep while %R is still flat.
The shift is left blank when it cannot be estimated: a period with fewer than 10 values or fewer
than 5 uncapped ones, or a fit whose mean lies outside the range the panel tested.

The model fits one normal distribution per period. Many drugs have a wild-type and a resistant
population; in the simulation, a rise in the resistant share alone, with the wild type unchanged,
showed as a shift of about 0.3 dilutions, significant in most datasets. Read it as a summary of
the bulk of the distribution, not as a breakpoint-free measure of change. A mixture model with
censoring ([PLOS ONE 2020](https://doi.org/10.1371/journal.pone.0220427)) would fix this. With
tens of thousands of isolates, a shift of 0.1 dilution is significant without mattering. Zone
diameters get no shift estimate.

## What the automated tests check

`tests/testthat/test_drift.R` builds MICs and zones by hand, or draws them from a known
distribution, so the true answer is known. Most tests use CLSI's 2023 cut to the Enterobacterales
gentamicin breakpoints (MIC S ≤ 4 / R ≥ 16 to S ≤ 2 / R ≥ 8; disk 14 mm from I to R), and one
uses EUCAST's 2021 cut to piperacillin/tazobactam. The tests check that:

- the same 100 MICs read under both editions give a change at fixed criteria of 0 and
  reclassification of +30 points
- a shifted MIC distribution under unchanged breakpoints gives no reclassification
- the two parts always sum to the observed change
- "≤4" counts as not R under R ≥ 8, is left out only for %NS, and a mostly capped panel gives no
  false reclassification; ">4" is counted under R ≥ 8 only; values straddling one edition's
  boundary only ("≤8", ">4" against R ≥ 16) make the pair *not interpretable*, one such value
  among 300 does not
- an *E. coli* reading from the generic Gram-negative fallback counts as not readable
- zones and MICs for one drug are decomposed apart, and EUCAST rows get EUCAST editions
- verdicts separate reclassification, change at fixed criteria, both, and a three-isolate pair
  (no change detected), and reclassification four times the change at fixed criteria is named
  larger with 1,000, 5,000 and 20,000 isolates per period
- reported calls replace the edition's reading as the observed change, and a laboratory that
  followed each year's edition has the revision as its whole reclassification
- %NS ignores a revision that only moves the I/R boundary
- bootstrap intervals cover the estimate and the true values and are reproducible with a seed;
  the bootstrap p-values follow the Shapley parts and are 1 for a part that is 0 in every
  replicate; the cluster bootstrap reports its clusters, refuses fewer than 5 per period and warns
  below 10
- a significant observed change that neither part shows is *unattributed*, and a part under one
  point does not count; with one part exactly 0 the pair is never *unattributed*, and `larger`
  never names a part the verdict does not show
- pairs that are not interpretable stay out of the corrections
- the censored model recovers a known one-dilution shift from capped panel values, finds a shift
  while %R stays at 0, and returns blank, with the reason, when a period is all capped, too small,
  or fitted outside the tested range
- an adoption table overrides the year rule, pools laboratories on different editions correctly,
  and separates a laboratory that moved from CLSI to EUCAST
- extra breakpoint rows replace the AMR package's row for the same edition, organism and drug,
  and keep CLSI's "R at or above the breakpoint" rule

## Validation on US national surveillance (NARMS)

NARMS, the US National Antimicrobial Resistance Monitoring System, samples retail meat and
food-animal guts on a fixed plan and tests every isolate by broth microdilution in one laboratory
network. FDA publishes the isolate-level MICs, with NARMS's own S/I/R call, free and without
registration. `Rscript scripts/narms.R` (about 10 minutes, most of it the bootstrap) uses the
*Salmonella* and *E. coli* records: 229,543 sample records, of which 163,241 had no growth,
leaving 65,273 isolates with an MIC and 770,135 MIC tests from 2002 to 2023. NARMS's 2023 calls
the package cannot read (azithromycin, no CLSI Enterobacterales breakpoint) are counted in
`narms_agreement.csv`.

**Agreement.** Read with the CLSI edition of their year, 2020–2023 cecal isolates agree with
NARMS's calls in 99.2–99.5% of tests per year (26,925 of 27,060 in 2023,
[`narms_agreement.csv`](../results/narms_agreement.csv)). Every 2023 disagreement is a rule NARMS
applies differently from CLSI 2023, none a misread value: gentamicin MICs of 4 and 8 (NARMS still
uses the pre-2023 breakpoint, S ≤ 4 / R ≥ 16) and *E. coli* ciprofloxacin MICs of 0.12 and 0.25
(NARMS calls them I, as for *Salmonella*; CLSI's Enterobacterales breakpoint makes them S). The
reading itself is a table lookup, so it can only be as wrong as the table.

**NARMS's 2023 change of reporting rule** ([`narms_relabel.md`](../results/narms_relabel.md)).
Through 2022, NARMS called a cecal *Salmonella* or *E. coli* ciprofloxacin MIC of 0.12 mg/L or
more R; in 2023 it moved to CLSI, where 0.12–0.5 is I. With the `reported` calls as the observed
change:

| 2022 to 2023, cecal | Reported %R | Observed | Fixed criteria | Reclassification | Verdict |
|---|--:|--:|--:|--:|---|
| *Salmonella*, ciprofloxacin | 11.4 to 0.6 | −10.8 | +0.3 (−0.2 to +0.8) | −11.1 (−12.7 to −9.5) | reclassification |
| *E. coli*, ciprofloxacin | 7.6 to 1.8 | −5.8 | +0.6 (−0.4 to +1.6) | −6.4 (−7.8 to −5.1) | reclassification |

Bootstrap 95% intervals, 5,000 replicates over isolates. The cecal file has no sampling site;
the sensitivity analyses below resample whole collection months instead.

Over the same year the table separates changes at fixed criteria from the change of rule:
*Salmonella* tetracycline +5.7 and nalidixic acid +5.5 points, both changes at fixed criteria with
no reclassification.

The change of rule also hid a rise. Counting non-susceptible isolates (I or R) instead
([`narms_relabel_ns.csv`](../results/narms_relabel_ns.csv)), *Salmonella* ciprofloxacin went from
11.4% to 17.0% between 2022 and 2023: +5.6 points (2.9 to 8.2), all change at fixed criteria.
Reported %R fell to 0.6% in the same year
([`fig_cip_salmonella.pdf`](../results/fig_cip_salmonella.pdf),
[`narms_cip_series.csv`](../results/narms_cip_series.csv)).

**Checks** ([`narms_ecoff.csv`](../results/narms_ecoff.csv),
[`narms_strata.csv`](../results/narms_strata.csv), [`narms_flow.csv`](../results/narms_flow.csv),
[`narms_drugs.csv`](../results/narms_drugs.csv)). Read against the EUCAST wild-type cut-off
(0.064 mg/L), which no breakpoint revision moves, non-wild-type *Salmonella* ciprofloxacin also
rose from 11.4% to 17.0%. Part of that rise is composition: standardised to the pooled mix of host
species it is +2.1 points, to the pooled mix of serotypes +2.8, and to host and serotype jointly
+0.8 (−0.6 to +2.3); within *S.* Infantis, mostly from chickens, it rose from 74% to 89%. The flow
table counts records, de-duplication and readable results by file and year.

**Sensitivity analyses** ([`narms_cluster.csv`](../results/narms_cluster.csv),
[`narms_shift_state.csv`](../results/narms_shift_state.csv),
[`narms_mic_model.csv`](../results/narms_mic_model.csv),
[`narms_panels.csv`](../results/narms_panels.csv)). The cecal file records no sampling site and
one isolate per organism per sample, so the cecal intervals are re-run resampling whole collection
months; the retail MIC shift is refitted with a bootstrap over the 10 states sampled in both
periods. The MIC-shift model is refitted with a separate standard deviation per period and
checked against the observed distribution, and the panel table lists every drug whose tested range
changed between years. `results/sessionInfo.txt` records the R session of the last run.

**Breakpoint-table audit** (`Rscript scripts/whonet_audit.R`;
[`whonet_audit.csv`](../results/whonet_audit.csv), [`whonet_example.csv`](../results/whonet_example.csv)).
Counts the CLSI Table 2A (Enterobacterales) rows per edition 2011-2020 in WHONET's
`Breakpoints.txt`, pinned to AMRIE commit 69de201 and checked by SHA-256, and in the AMR package,
and records which breakpoint row `as.sir()` uses for an *E. coli* MIC of 1 mg/L under CLSI 2019
and 2020.

**A revision NARMS has not applied yet** ([`narms_editions.md`](../results/narms_editions.md)).
In 2023 NARMS still reported cecal gentamicin with the old CLSI breakpoint (an MIC of 8 is I),
while its 2024 refresh of the retail file applies the new one backwards (8 is R). Reading
2019–2020 and 2022–2023 cecal isolates with CLSI 2020 and CLSI 2023, *Salmonella* gentamicin %R
rises from 1.7% to 5.8%: +3.4 points reclassification, +0.7 change at fixed criteria, and the MIC
distribution moved −0.2 dilutions (95% CI −0.3 to −0.1). The bulk of the distribution moved down
(MICs ≤ 0.25 went from 27.8% to 38.9%) while the tail at 8 or more grew slightly, and %R follows
the tail. Read with CLSI 2023 instead of CLSI 2020, the same 2022–2023 isolates go from 2.3% to
5.8% R: once NARMS adopts the 2023 breakpoint, its *Salmonella* gentamicin resistance will more
than double with no change in the bacteria.

**Twenty years of MIC shift** ([`narms_mic_shift.csv`](../results/narms_mic_shift.csv)). On retail
isolates, 2002–2006 against 2017–2021, the shift model needs no breakpoint, so the AMR package's
gaps before 2020 do not matter. *Salmonella* nalidixic acid MICs rose 0.6 dilutions (0.5 to 0.7)
and *E. coli* nalidixic acid fell 0.5; the other eight estimated shifts are under 0.25 dilutions.
Ten of the 20 pairs tested in both periods were left blank, because most of their values sat at
the panel's lowest dilution and the fit fell outside the tested range.

NARMS isolates come from food animals and meat, not patients, and the breakpoints are human
clinical ones; that is fine for testing the method, which does not depend on where an isolate
came from.

## A revision applied by laboratories (NCBI)

`Rscript scripts/ncbi_fluoroquinolones.R` reads CLSI's 2019 ciprofloxacin and levofloxacin
revision for Enterobacterales other than *Salmonella* in the calls of NCBI BioSample antibiograms
of human isolates, 2014–2018 against 2023–2025
([`ncbi_fluoroquinolones.csv`](../results/ncbi_fluoroquinolones.csv)). The AMR package has no
Enterobacterales fluoroquinolone rows before CLSI 2020, so the earlier breakpoints come from
`data/breakpoints_clsi_fluoroquinolones.csv` through the `breakpoints` argument. Submitters chose
which isolates to sequence, so the change at fixed criteria describes no population;
reclassification reads the same isolates twice and does not depend on how they were chosen.

[`results/ncbi_report.md`](../results/ncbi_report.md) is a second example on NCBI BioSample data
(`Rscript scripts/ncbi.R`, about a minute): 16,832 MIC and disk tests (CLSI and EUCAST) from 1,496
human isolates, 2016–2020 against 2021–2025. It shows that the package runs on real, messy exports
of human isolates. It does not show real resistance trends: submitters chose which isolates to
sequence, often because they were resistant, and which laboratories submitted changed from year
to year.

## Data

- `data/raw/ncbi_ast_long.csv.gz` is the NCBI BioSample snapshot made by `R/00_fetch_ncbi.R` in
  [amr-surveillance-r](https://github.com/FerrumSL/amr-surveillance-r), copied unchanged
  (SHA-256 `a4ed673b…0b498`).
- `data/raw/ncbi_enterobacterales_ast.csv.gz` holds every Enterobacterales BioSample with an
  antibiogram, fetched by `scripts/fetch_ncbi_enterobacterales.py`.
- `data/raw/narms_extract.csv.gz` holds the *Salmonella* and *E. coli* columns of two FDA NARMS
  files, *Retail Meats* and *Cecal* (updated 27 September 2024), from the
  [NARMS download page](https://www.fda.gov/animal-veterinary/national-antimicrobial-resistance-monitoring-system/integrated-reportssummaries).
  `scripts/narms.R` lists their SHA-256 and rebuilds the extract if they are put in `data/narms/`.
  Source: Food and Drug Administration (FDA), NARMS Now, Rockville, MD: U.S. Department of Health
  and Human Services.

R 4.6.2, the same version CI uses; packages pinned in `renv.lock` (`renv::restore()`).
