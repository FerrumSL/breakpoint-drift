# breakpointdrift 1.2.0

## Renamed (breaking)

The paper, the README and the code now use one set of names. Results are unchanged; only names
differ.

- The two components are **reclassification** (column `reclass`, was `artefact`) and the
  **change at fixed criteria** (column `fixed`, was `real`). Every derived column follows:
  `fixed_low`, `fixed_high`, `fixed_min`, `fixed_max`, `fixed_at_new`, `fixed_at_old`,
  `boot_p_fixed`, `p_fixed`, and the same with `reclass`.
- Verdicts: "guideline change" is now "reclassification", "period change" "change at fixed
  criteria", "unresolved" "unattributed" and "noise" "no change detected". "both" and "not
  interpretable" are unchanged.
- Column `larger` says "reclassification" or "fixed criteria" (was "guideline" or "period").
- `drift_report()` heads its columns "Fixed criteria" and "Reclassification".
- The supplementary tables, the simulation outputs and the keys of `results/narms_numbers.csv`
  use the same words.

## Since 1.1.0

Apart from the names, no real-data result in `results/` changes. New files are listed below, and the simulation
(`results/simulation*`) was rerun with 1,000 datasets per cell instead of 200.

## New

- `drift_decompose()` and `drift_bootstrap()` take `breakpoints`: MIC rows the AMR package
  lacks, which replace its row for the same guideline, organism and drug. Editions without such
  rows are read as before. The AMR package reads any table but its own with EUCAST's "R if
  MIC > breakpoint" rule; the CLSI rows of an extended table are adjusted so that CLSI's
  "R if MIC >= breakpoint" still holds (all 1,578 of its own MIC rows of CLSI 2018 and 2023 read
  the same; the tests check every tenth).
- `scripts/ncbi_fluoroquinolones.R`: CLSI's 2019 ciprofloxacin and levofloxacin revision in
  the calls of NCBI BioSample antibiograms (all Enterobacterales BioSamples with one, fetched
  2026-11-17 by `scripts/fetch_ncbi_enterobacterales.py`), with the earlier breakpoints from
  `data/breakpoints_clsi_fluoroquinolones.csv`. Writes `results/ncbi_fluoroquinolones*.csv`.

## Faster

- `drift_bootstrap()` counts the draws per kind of unit (units whose rows are alike) instead of
  per unit, and `mic_ends()` parses each distinct value once. The random draws and every result
  are unchanged; a simulated dataset with 5,000 isolates per period runs about 2.7 times as fast.

## Fixed

- Optional columns (`isolate`, `reported`, `org`, `lab`) were looked up with `$`, which also
  takes a column whose name only starts that way. A column `isolate_source` with no `isolate`
  column became the bootstrap unit, so two sources gave two units per period and a zero-width
  interval; `reported_date` made the observed change NA; `organism` stopped with an error.
  They are now matched by their exact name.
- `drift_bootstrap()` put every row with a missing `isolate` ID into one unit. Such a row is now
  its own unit.
- A capped value off the dilution grid, such as ">1.5" from a gradient strip, was read as
  above 2 (">1.5" started at 4). It now starts at the next grid value, 2.
- `drift_report()` counted pairs made not interpretable by straddling capped values among the
  interpretable pairs in its header.

## Changed

- `guideline_in_force()` warns when the AMR package lacks the edition of a year between two
  editions it has, and says that the latest earlier edition is used.

# breakpointdrift 1.1.0

## Fixed

- Capped MICs were left out whenever AMR called them NI, even where the outcome was determined:
  "<=4" is not R under R >= 8 although it may be S or I. A capped value is now an interval of
  possible MICs, read at both ends, and left out only if it falls on both sides of the outcome's
  boundary (R, or NS). With short panels the old rule inflated the guideline component (a +80
  point component from 80 values at "<=4" and 20 at "16", read the same in both years). NARMS
  results are not affected: every pair was read on the same isolates under both guidelines.
- A capped value that straddles the boundary under one guideline only ("<=8" against R >= 8 but
  not R >= 16) is still left out of that guideline's share only, so the two shares of a period
  rest on different isolates and the components are not identified. `drift_decompose()` now
  counts such values (`straddled_y1`, `straddled_y2`) and gives the range each component can lie
  in (`real_min`, `real_max`, `artefact_min`, `artefact_max`); `drift_test()` makes the pair
  "not interpretable" when that range is `min_effect` or wider, and says why in `verdict_note`. Before,
  80 values at "<=8" and 20 at "16" in both years gave a +80 point "guideline change".
- `drift_bootstrap()` with `cluster` keeps an isolate with no cluster whole (all its drugs) instead
  of drawing each of its rows alone.
- `drift_report()` takes `outcome`, so a %NS table is no longer labelled %R.
- `drift_bootstrap()` with `cluster` stops when a period has fewer than `min_clusters` (default
  5) clusters and warns below 10. With one cluster per period it gave a zero-width interval and
  the smallest p-value.
- `larger` names a component only when that component is itself shown, and "unresolved" needs
  both components non-zero. The test of |guideline| - |period| and the component tests are
  corrected in separate families, so with a guideline component of exactly 0 a pair could be
  "noise" with the period named larger.
- `guideline_in_force()` takes the latest edition on or before the year, so a gap in the
  editions is never filled by a name the AMR package lacks.
- Periods that are not runs of years are labelled as sets ("2015, 2018"), not as a range.
- Pairs that are not interpretable no longer enter the Benjamini-Hochberg corrections.

## New

- `drift_bootstrap()`: percentile intervals for the observed change and both components,
  resampling isolates (or clusters, `cluster =`, drawing G - 1 of the G clusters of a period)
  within each period, and bootstrap p-values for the Shapley period component (`boot_p_real`),
  the guideline component (`boot_p_artefact`) and |guideline| - |period| (`boot_p_larger`).
  `boot_na` counts replicates with an empty share; `units_y1`, `units_y2` the units resampled.
- `drift_decompose()` returns both orderings of the period component (`real_at_new`,
  `real_at_old`) and their difference (`interaction`). `up_y*` and `down_y*` count the
  reclassified isolates descriptively.
- `outcome = "NS"` counts non-susceptible isolates (I, SDD or R) instead of resistant ones.
- `drift_mic_shift()` says why it withheld a shift (`note`).

## Changed tests and verdicts

- `drift_test()` needs a `drift_bootstrap()` table. Fisher's exact test is kept for the observed
  change only; the period and guideline components are tested by the bootstrap, so each test is
  of the quantity reported. The components are not tested with an
  exact McNemar test: with a one-way revision its "down" count is always 0, it weights periods by
  n, and in reported mode it counts half- and full-weight isolates alike.
- Verdicts no longer depend on which test is checked first: "both" if both components are shown,
  "guideline change" or "period change" if only one is, "unresolved" if neither is but the
  observed change is, "noise", "not interpretable". A component is shown when its corrected
  p-value is below `alpha` and it is at least `min_effect` (default 1 percentage point). Column
  `larger` names the larger component where the test of |guideline| - |period| is significant.

## Internal

- Each distinct MIC or zone is read once, so large data sets read much faster. Readings dropped
  as a generic Gram-negative fallback are marked (`attr(, "fallback")`).

## NARMS validation (scripts/narms.R)

- De-duplication key is now file + organism + serotype + isolate + drug, and records without an isolate ID
  are kept as separate isolates. The old key (file + isolate + drug) merged all records without
  an ID into one and dropped the Salmonella that shares an ID with an E. coli from the same
  retail sample.
- B = 5,000 bootstrap replicates (seed 2026).
- New outputs: `narms_cip_series.csv` and `fig_cip_salmonella.pdf` (with n per year, the 2023
  rule change and a generated caption), `narms_relabel_ns.csv` (%NS), `supplement_tables.tex`
  (n per period, `p_artefact`, structural zeros printed as 0), `narms_ecoff.csv` (%non-wild-type
  by EUCAST ECOFF), `narms_strata.csv` (Salmonella ciprofloxacin %NS by host and serotype, with
  a directly standardised change), `narms_flow.csv` (records, de-duplication and readability by
  file and year) and `narms_drugs.csv` (which drug reaches which table).
- `scripts/simulation.R` checks bias, interval coverage and verdict rates on simulated MICs with
  known components: %R and %NS, short panels, a wild type near the breakpoint, reported calls
  (a compliant lab and one with its own rule), clustered sites, 20 pairs with Benjamini-Hochberg
  correction, and the MIC-shift model under a panel change and a bimodal change
  (`results/simulation.csv`, `.md`). Each cell has its own seed.
- The extract now carries the collection month (cecal) and state (retail). New sensitivity
  outputs: `narms_cluster.csv` (cecal bootstrap intervals resampling whole collection months),
  `narms_shift_state.csv` (retail MIC shift refitted with a bootstrap over states),
  `narms_mic_model.csv` (MIC-shift model with a separate SD per period, likelihood-ratio test of
  equal SD, largest gap between fitted and observed CDF), `narms_panels.csv` (tested range by
  file, drug and year) and `sessionInfo.txt`. `supplement_tables.tex` adds p-values for the
  observed and period tests and Tables S4-S7.
- `scripts/whonet_audit.R` counts the CLSI Table 2A rows per edition in WHONET's
  `Breakpoints.txt` (pinned to AMRIE commit 69de201, SHA-256 checked) and in the AMR package,
  and shows which row `as.sir()` uses for an E. coli MIC under CLSI 2019 and 2020
  (`whonet_audit.csv`, `whonet_example.csv`, `supplement_whonet.tex`).
- The extract carries the GROWTH flag. Of 229,543 Salmonella and E. coli sample records, 163,241
  had no growth; 65,273 isolates have at least one MIC. Table S4 counts records, no-growth
  records, isolates and isolates with an MIC. `narms_codes.csv` lists every drug code in the two
  files and why it is or is not analysed. De-duplication also keys on serotype (two retail 2019
  isolate IDs are shared by Salmonella of different serotypes).
- `narms_shift_state.csv` now refits the retail shift on the states sampled in both periods
  with a bootstrap over whole states, replacing the sandwich standard errors from 10 and 22
  clusters. `narms_cluster.csv` gives the number of months per period.
- New: `narms_windows.csv` (other windows for both comparisons), `narms_fallback_check.csv`
  (the false guideline components the generic Gram-negative fallback would produce with CLSI
  2019), `narms_standardised.csv` (host, serotype and joint standardisation with bootstrap
  intervals), `narms_numbers.csv` (every number the paper quotes that no table holds),
  `supplement_extra.tex` (Tables S9-S12), Clopper-Pearson intervals in `narms_cip_series.csv`
  and the figure, and unread NARMS calls (`calls` - `tests`) in `narms_agreement.csv`.
- `scripts/ncbi.R` uses the bootstrap tests (B = 2,000, seed 2026).
- `scripts/simulation.R` adds cells for rare resistance (0.5% and 2%) and for capped values
  straddling the revised R breakpoint, and reports the larger-component label and guideline
  coverage under the site bootstrap. `scripts/simulation_table.R` writes the summary table
  (`results/simulation_table.tex`).
- `scripts/narms.R` adds the gentamicin MIC distributions (`fig_gen_mic.pdf`,
  `narms_gen_distribution.csv`), self-contained table captions, and labels the supplementary
  tables with the paper's terms: criteria component (`artefact`) and fixed-criteria component
  (`real`). The software's column names and verdicts are unchanged.
