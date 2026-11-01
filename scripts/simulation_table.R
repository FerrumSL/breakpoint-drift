# Summary table of the simulation for the main text, from results/simulation.csv.
# Run from the project root after scripts/simulation.R: Rscript scripts/simulation_table.R
# Each row pools the cells of one scenario over n = 200, 1,000 and 5,000; ranges are over cells.

out <- read.csv("results/simulation.csv")
rng <- function(v, f = function(x) sprintf("%.1f", 100 * x)) {
  v <- v[!is.na(v)]
  if (length(v) == 0) return("--")
  if (diff(range(v)) < 5e-4) f(v[1]) else paste0(f(min(v)), "--", f(max(v)))
}
big <- function(v) if (all(is.na(v))) "--" else sprintf("%.1f", 100 * max(abs(v), na.rm = TRUE))
row <- function(label, x, site = FALSE) {
  g <- x[x$shift > 0, ]
  cover_r <- if (site) x$cover_real_site else x$cover_real
  cover_g <- if (site) g$cover_artefact_site else g$cover_artefact
  false_p <- if (site) x$period_shown_site[x$delta == 0] else x$period_shown[x$delta == 0]
  # Every pair not interpretable: no estimate is reported, so bias and coverage are not shown
  if (all(x$not_interpretable == 1)) return(sprintf("%s & -- & -- & -- & -- & %s & %s \\\\", label, rng(false_p), rng(x$not_interpretable)))
  sprintf("%s & %s & %s & %s & %s & %s & %s \\\\", label, big(x$bias_real), big(x$bias_artefact),
          rng(cover_r), rng(cover_g), rng(false_p), rng(x$not_interpretable))
}
b <- out$set == "base"
rows <- c(
  row("\\%R, no capped values", out[b & out$capped == 0, ]),
  row("\\%R, 30\\% capped, one side of the breakpoints", out[b & out$capped > 0, ]),
  row("\\%R, 30\\% capped, straddling the revised R breakpoint", out[out$set == "straddle", ]),
  row("\\%NS, 30\\% capped, straddling the revised S breakpoint", out[out$set == "NS" & out$capped > 0, ]),
  row("Rare resistance (0.5\\% or 2\\%)", out[out$set == "low", ]),
  row("Wild type near the revised breakpoint", out[out$set == "near", ]),
  row("Reported calls, laboratory following each edition", out[out$reported == "compliant", ]),
  row("20 sites per period, isolate bootstrap", out[out$set == "clustered", ]),
  row("20 sites per period, site bootstrap", out[out$set == "clustered", ], site = TRUE))
writeLines(c(
  "\\begin{table}[h]", "\\centering", "\\scriptsize", "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{lrrllll}", "\\toprule",
  "Scenario & \\multicolumn{2}{c}{Largest bias} & \\multicolumn{2}{c}{Coverage (\\%)} & False fixed- & Not \\\\",
  " & Fixed & Criteria & Fixed & Criteria & criteria change (\\%) & interpretable (\\%) \\\\", "\\midrule",
  rows, "\\bottomrule", "\\end{tabular}",
  paste("\\caption{Simulation summary. Each scenario pools its cells over $n$ = 200, 1,000 and 5,000 isolates per period,",
        "200 datasets per cell, 999 bootstrap replicates each; ranges are over cells. Largest bias: largest absolute mean",
        "error over the cells, percentage points. Coverage of the 95\\% intervals; criteria coverage only where a breakpoint",
        "moved. Fixed: fixed-criteria component. False fixed-criteria change: share of datasets with no true change under",
        "fixed criteria in which one was detected.",
        "Not interpretable: share of datasets so labelled; where every dataset was, no estimate is reported. Every cell is in Supplementary file \\texttt{simulation.csv}.}"),
  "\\label{tab:sim}", "\\end{table}"), "results/simulation_table.tex")
