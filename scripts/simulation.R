# Simulation: does the decomposition recover a known period component and a known guideline
# component, and how often does each verdict appear?
# Run from the project root: Rscript scripts/simulation.R  (about 50 minutes on 16 cores)
#
# Each period draws n true log2 MICs from a mixture: a wild type N(r - 4, 1.3) (or N(r - 2, 1.3),
# "near": its tail reaches the revised breakpoint) and a resistant share N(r + 1, 1), where r is
# the old guideline's log2 R breakpoint. The resistant share is `prev` in period 1 and
# `prev + delta` in period 2. Values are stored as a panel does, rounded up to the next dilution
# and capped at the ends; a share `capped` of isolates is on a short panel that starts at the old
# S breakpoint ("<=4" for gentamicin), or, with `cap_at` "R", at the revised R breakpoint ("<=8"),
# so that the capped value straddles it. The guidelines are real revisions in the AMR package:
#   shift 0  E. coli gentamicin, CLSI 2023 to 2024 (unchanged: S <= 2, R >= 8)
#   shift 1  E. coli gentamicin, CLSI 2022 to 2023 (S <= 4 / R >= 16 to S <= 2 / R >= 8)
#   shift 2  E. coli piperacillin/tazobactam, CLSI 2021 to 2022 (S <= 16 / R >= 128 to S <= 8 / R >= 32)
# The true components are those of the latent MICs with nothing capped. Cells:
#   base       %R, every combination of n, prev, delta, shift and capping
#   NS         %NS, where a short panel's "<=4" straddles the revised S breakpoint
#   near       wild type near the breakpoint
#   reported   reported calls: a lab following each year's edition ("compliant"), or one whose own
#              rule already called the revised breakpoint's R in period 1 ("own rule", shift 1)
#   clustered  20 sites per period whose resistant share varies (log-normal, SD 0.5, mean prev),
#              analysed with the isolate and with the site bootstrap
#   low        rare resistance: resistant share 0.5% or 2%, period change 0 or +1 point
#   straddle   %R with 30% of values capped at "<=8", which straddles the revised R breakpoint (shift 1);
#              these pairs should be "not interpretable"
# Then two further checks: 20 pairs in one dataset with Benjamini-Hochberg correction, and the
# MIC-shift model under a panel change and under a bimodal change.

source("R/drift.R")
library(parallel)

# Every dataset holds the same few stored values, so each is read with AMR once per edition
read_sir_amr <- read_sir
read_sir <- local({
  cache <- list()
  function(x, guideline) {
    key <- paste(guideline, x$mo, x$ab, x$mic)
    new <- !duplicated(key) & !key %in% names(cache)
    if (any(new)) cache[key[new]] <<- as.list(read_sir_amr(x[new, ], guideline))
    unlist(cache[key], use.names = FALSE)
  }
})

revisions <- list(
  `0` = list(ab = "GEN", years = c(2023, 2024), r = c(3, 3), s = c(1, 1), panel = -2:5),
  `1` = list(ab = "GEN", years = c(2022, 2023), r = c(4, 3), s = c(2, 1), panel = -2:5),
  `2` = list(ab = "TZP", years = c(2021, 2022), r = c(7, 5), s = c(4, 3), panel = 0:8))

# P(stored MIC >= 2^t): the stored value is 2^ceiling(z), so z > t - 1
p_at <- function(prev, t, rev, wt) {
  (1 - prev) * stats::pnorm(t - 1, rev$r[1] - wt, 1.3, lower.tail = FALSE) +
    prev * stats::pnorm(t - 1, rev$r[1] + 1, 1, lower.tail = FALSE)
}

store <- function(z, rev, short, cap_at = "S") {
  k <- pmin(pmax(ceiling(z), min(rev$panel)), max(rev$panel) + 1)
  low <- ifelse(short, if (cap_at == "R") rev$r[2] else rev$s[1], min(rev$panel))
  ifelse(k <= low, paste0("<=", 2^low), ifelse(k > max(rev$panel), paste0(">", 2^max(rev$panel)), as.character(2^k)))
}

draw <- function(n, prev, rev, wt, capped, sites = 0, cap_at = "S") {
  site <- NA
  if (sites > 0) {  # each site's resistant share, mean prev
    site <- rep(seq_len(sites), length.out = n)
    prev <- pmin(1, prev * exp(stats::rnorm(sites, -0.125, 0.5)))[site]
  }
  res <- stats::runif(n) < prev
  z <- ifelse(res, stats::rnorm(n, rev$r[1] + 1, 1), stats::rnorm(n, rev$r[1] - wt, 1.3))
  data.frame(mic = store(z, rev, stats::runif(n) < capped, cap_at), site = site)
}

# The true components of a cell
truth <- function(cell) {
  rev <- revisions[[as.character(cell$shift)]]
  wt <- if (cell$wt == "near") 2 else 4
  t <- if (cell$outcome == "R") rev$r else rev$s + 1
  p <- outer(c(cell$prev, cell$prev + cell$delta), t, Vectorize(function(q, k) p_at(q, k, rev, wt)))
  real <- ((p[2, 2] - p[1, 2]) + (p[2, 1] - p[1, 1])) / 2
  observed <- p[2, 2] - p[1, 1]
  if (cell$reported == "own rule") observed <- p[2, 2] - p_at(cell$prev, rev$r[1] - 1, rev, wt)
  c(real = real, artefact = observed - real)
}

one_run <- function(cell, B) {
  rev <- revisions[[as.character(cell$shift)]]
  wt <- if (cell$wt == "near") 2 else 4
  sites <- if (cell$clustered) 20 else 0
  data <- rbind(draw(cell$n, cell$prev, rev, wt, cell$capped, sites, cell$cap_at),
                draw(cell$n, cell$prev + cell$delta, rev, wt, cell$capped, sites, cell$cap_at))
  data <- data.frame(year = rep(rev$years, each = cell$n), mo = "B_ESCHR_COLI", ab = rev$ab, mic = data$mic,
                     isolate = seq_len(2 * cell$n), site = paste(rep(1:2, each = cell$n), data$site))
  if (cell$reported == "compliant") {
    x <- transform(data, method = "MIC", disk = NA)
    data$reported <- ifelse(data$year == rev$years[1], read_sir(x, paste("CLSI", rev$years[1])),
                            read_sir(x, paste("CLSI", rev$years[2])))
  }
  if (cell$reported == "own rule") {
    value <- as.numeric(sub("^[<>=]+", "", data$mic))
    data$reported <- ifelse(value >= 2^(rev$r[1] - 1), "R", "S")
  }
  t <- drift_test(drift_bootstrap(data, B = B, outcome = cell$outcome))
  out <- c(real = t$real, artefact = t$artefact, real_low = t$real_low, real_high = t$real_high,
           artefact_low = t$artefact_low, artefact_high = t$artefact_high,
           verdict = match(t$verdict, verdict_levels), larger = match(t$larger, c("guideline", "period", "neither")),
           excluded = 1 - (t$n_y1g2 + t$n_y2g2) / (t$tested_y1 + t$tested_y2))
  if (cell$clustered) {
    s <- drift_test(drift_bootstrap(data, B = B, outcome = cell$outcome, cluster = "site"))
    out <- c(out, real_low_site = s$real_low, real_high_site = s$real_high,
             artefact_low_site = s$artefact_low, artefact_high_site = s$artefact_high,
             verdict_site = match(s$verdict, verdict_levels),
             larger_site = match(s$larger, c("guideline", "period", "neither")))
  }
  out
}

cells <- function(set, ...) {
  g <- expand.grid(n = c(200, 1000, 5000), ..., stringsAsFactors = FALSE)
  defaults <- list(prev = 0.20, delta = c(0, 0.05), shift = 1, capped = 0, outcome = "R", wt = "far",
                   reported = "none", clustered = FALSE, cap_at = "S")
  for (k in names(defaults)) if (is.null(g[[k]])) g <- merge(g, setNames(data.frame(defaults[[k]], stringsAsFactors = FALSE), k))
  cbind(set = set, g[c("n", names(defaults))])
}
grid <- rbind(
  cells("base", prev = c(0.05, 0.20), delta = c(0, 0.05), shift = 0:2, capped = c(0, 0.3)),
  cells("NS", shift = 1:2, capped = c(0, 0.3), outcome = "NS"),
  cells("near", prev = c(0.05, 0.20), shift = 0:2, wt = "near"),
  cells("reported", shift = 0:2, reported = "compliant"),
  cells("reported", shift = 1, reported = "own rule"),
  cells("clustered", shift = 0:1, clustered = TRUE),
  cells("low", prev = c(0.005, 0.02), delta = c(0, 0.01), shift = 0:1),
  cells("straddle", capped = 0.3, cap_at = "R"))
grid <- grid[!(grid$shift == 0 & grid$capped > 0), ]  # capping at S <= 2 changes nothing when S stays 2
grid <- grid[order(grid$set != "base"), ]
rownames(grid) <- NULL
reps <- as.integer(Sys.getenv("SIM_REPS", "200"))  # SIM_REPS=2 for a quick check
B <- as.integer(Sys.getenv("SIM_B", "999"))
# 200 datasets per cell; Monte Carlo error on a 95% coverage is about 1.5 points

cores <- max(1, detectCores() - 1)
started <- Sys.time()
# Each cell has its own seed, so the results do not depend on the number of cores or the platform
runs <- mclapply(seq_len(nrow(grid)), function(j) {
  set.seed(2026 + j)
  t(replicate(reps, one_run(grid[j, ], B)))
}, mc.cores = cores, mc.preschedule = FALSE)

out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(j) {
  cell <- grid[j, ]
  tr <- truth(cell)
  x <- as.data.frame(runs[[j]])
  v <- verdict_levels[x$verdict]
  has_g <- abs(tr[["artefact"]]) > 1e-9
  has_p <- abs(tr[["real"]]) > 1e-9
  true_larger <- if (!has_g && !has_p) NA else if (abs(tr[["artefact"]]) > abs(tr[["real"]])) "guideline" else "period"
  larger <- c("guideline", "period", "neither")[x$larger]
  wrong <- function(l) if (is.na(true_larger)) mean(l != "neither") else mean(l != true_larger & l != "neither")
  site <- if (cell$clustered) c(
    cover_real_site = mean(x$real_low_site <= tr[["real"]] & tr[["real"]] <= x$real_high_site),
    cover_artefact_site = mean(x$artefact_low_site <= tr[["artefact"]] & tr[["artefact"]] <= x$artefact_high_site),
    period_shown_site = mean(verdict_levels[x$verdict_site] %in% c("period change", "both")),
    larger_wrong_site = wrong(c("guideline", "period", "neither")[x$larger_site])) else
    c(cover_real_site = NA, cover_artefact_site = NA, period_shown_site = NA, larger_wrong_site = NA)
  data.frame(cell, true_real = tr[["real"]], true_artefact = tr[["artefact"]],
             bias_real = mean(x$real) - tr[["real"]], bias_artefact = mean(x$artefact) - tr[["artefact"]],
             cover_real = mean(x$real_low <= tr[["real"]] & tr[["real"]] <= x$real_high),
             cover_artefact = mean(x$artefact_low <= tr[["artefact"]] & tr[["artefact"]] <= x$artefact_high),
             t(site), excluded = mean(x$excluded),
             period_change = mean(v == "period change"), guideline_change = mean(v == "guideline change"),
             both = mean(v == "both"), unresolved = mean(v == "unresolved"), noise = mean(v == "noise"),
             guideline_shown = mean(v %in% c("guideline change", "both")),
             period_shown = mean(v %in% c("period change", "both")),
             larger_correct = if (is.na(true_larger)) NA else mean(larger == true_larger),
             larger_wrong = wrong(larger), not_interpretable = mean(v == "not interpretable"))
}))
dir.create("results", showWarnings = FALSE)
num <- vapply(out, is.numeric, logical(1))
out[num] <- lapply(out[num], signif, 4)
write.csv(out, "results/simulation.csv", row.names = FALSE)

# ---- The base cell with the lowest period coverage, again with 600 new datasets -----------------
# With 156 cells, some fall a few Monte Carlo errors from 95% by chance; this checks the worst.
worst <- which.min(ifelse(out$set == "base", out$cover_real, Inf))
recheck_runs <- do.call(rbind, mclapply(seq_len(600), function(k) {
  set.seed(90000 + k)
  one_run(grid[worst, ], B)
}, mc.cores = cores))
recheck <- data.frame(grid[worst, ], datasets = 600,
                      cover_real = mean(recheck_runs[, "real_low"] <= out$true_real[worst] & out$true_real[worst] <= recheck_runs[, "real_high"]),
                      period_shown = mean(verdict_levels[recheck_runs[, "verdict"]] %in% c("period change", "both")))
write.csv(recheck, "results/simulation_recheck.csv", row.names = FALSE)

# ---- Several pairs in one dataset, Benjamini-Hochberg over them ------------------------------------
# 2022 against 2023: gentamicin (revised, shift 1) and piperacillin/tazobactam (unchanged between
# these editions) for each of 10 Enterobacterales species, n = 1,000 per period and pair, resistant
# share 20%; half the species get a period change of +5 points in both drugs. The true guideline
# component is 0 for every piperacillin/tazobactam pair, and the period component 0 for half.
species <- as.character(AMR::as.mo(c("Escherichia coli", "Klebsiella pneumoniae", "Klebsiella oxytoca",
                                     "Enterobacter cloacae", "Proteus mirabilis", "Citrobacter freundii",
                                     "Serratia marcescens", "Klebsiella aerogenes", "Morganella morganii",
                                     "Providencia stuartii"), info = FALSE))
multi_run <- function(B) {
  data <- do.call(rbind, lapply(seq_along(species), function(i) do.call(rbind, lapply(c("GEN", "TZP"), function(ab) {
    rev <- if (ab == "GEN") revisions[["1"]] else modifyList(revisions[["2"]], list(r = c(5, 5), s = c(3, 3)))
    delta <- if (i <= 5) 0.05 else 0
    mics <- c(draw(1000, 0.2, rev, 4, 0)$mic, draw(1000, 0.2 + delta, rev, 4, 0)$mic)
    data.frame(year = rep(2022:2023, each = 1000), mo = species[i], ab = ab, mic = mics)
  }))))
  data$isolate <- seq_len(nrow(data))
  t <- drift_test(drift_bootstrap(data, B = B))
  t$true_period <- t$mo %in% species[1:5]
  t$true_guideline <- t$ab == "GEN"
  t[c("mo", "ab", "verdict", "true_period", "true_guideline")]
}
multi <- mclapply(seq_len(reps), function(k) {
  set.seed(4026 + k)
  cbind(dataset = k, multi_run(B))
}, mc.cores = cores, mc.preschedule = FALSE)
multi <- do.call(rbind, multi)
g_shown <- multi$verdict %in% c("guideline change", "both")
p_shown <- multi$verdict %in% c("period change", "both")
fdr <- function(shown, true) mean(tapply(shown & !true, multi$dataset, sum) / pmax(1, tapply(shown, multi$dataset, sum)))
multi_summary <- data.frame(
  measure = c("guideline shown, revised pairs", "guideline shown, unrevised pairs", "period shown, +5 points",
              "period shown, no period change", "false discovery rate, guideline", "false discovery rate, period",
              "pairs not interpretable"),
  value = c(mean(g_shown[multi$true_guideline]), mean(g_shown[!multi$true_guideline]), mean(p_shown[multi$true_period]),
            mean(p_shown[!multi$true_period]), fdr(g_shown, multi$true_guideline), fdr(p_shown, multi$true_period),
            mean(multi$verdict == "not interpretable")))

# ---- The MIC-shift model --------------------------------------------------------------------------
# (a) one normal distribution, shifted by 0.5 dilutions, on a panel whose tested range moves up
# two dilutions between the periods; (b) the base mixture (gentamicin, wild type far) with the
# resistant share rising from 20% to 25% and the wild type unchanged, which no single location
# shift describes. n = 1,000 per period.
shift_run <- function(case) {
  if (case == "panel") {
    mics <- c(store(stats::rnorm(1000, -1, 1.3), list(panel = -4:3), FALSE),
              store(stats::rnorm(1000, -0.5, 1.3), list(panel = -2:5), FALSE))
  } else {
    mics <- c(draw(1000, 0.20, revisions[["1"]], 4, 0)$mic, draw(1000, 0.25, revisions[["1"]], 4, 0)$mic)
  }
  s <- drift_mic_shift(data.frame(year = rep(c(2015, 2024), each = 1000), mo = "B_ESCHR_COLI", ab = "GEN", mic = mics))
  c(shift = s$shift, low = s$shift_low, high = s$shift_high, p = s$p_shift)
}
shift_sim <- do.call(rbind, lapply(c("panel", "bimodal"), function(case) {
  set.seed(if (case == "panel") 5026 else 6026)
  x <- as.data.frame(t(replicate(reps, shift_run(case))))
  target <- if (case == "panel") 0.5 else 0
  data.frame(case = case, target = target, estimated = mean(!is.na(x$shift)), mean_shift = mean(x$shift, na.rm = TRUE),
             cover = mean(x$low <= target & target <= x$high, na.rm = TRUE), rejected = mean(x$p < 0.05, na.rm = TRUE))
}))
minutes <- as.numeric(difftime(Sys.time(), started, units = "mins"))

# ---- Summary ----------------------------------------------------------------------------------------
pct <- function(v) sprintf("%.1f", 100 * v)
rng <- function(v) if (length(v) == 0) "" else if (diff(range(v)) < 5e-4) pct(v[1]) else paste0(pct(min(v)), "-", pct(max(v)))
by_n <- function(f, rows) paste(vapply(c(200, 1000, 5000), function(n) f(out[rows & out$n == n, ]), character(1)), collapse = " | ")
range_pct <- function(col) function(x) rng(x[[col]])
points <- function(col, f = mean, d = 2) function(x) sprintf(paste0("%.", d, "f"), 100 * f(x[[col]]))
b <- out$set == "base"
uncapped <- b & out$capped == 0
capped <- b & out$capped > 0
row <- function(label, f, rows) paste("|", label, "|", by_n(f, rows), "|")
lines <- c(
  "# Simulation", "",
  sprintf("Written by scripts/simulation.R (one seed per cell from 2026, %d datasets per cell, B = %d, %d cells, %.1f minutes).",
          reps, B, nrow(out), minutes),
  "Setup in the script header. Bias in percentage points, the rest in percent. Ranges are over the",
  "cells with that n; per-cell results, with the true components, are in simulation.csv.", "",
  "## %R, base cells", "",
  "| | n = 200 | n = 1000 | n = 5000 |", "|---|--:|--:|--:|",
  row("Bias, period, largest absolute", points("bias_real", function(v) max(abs(v))), b),
  row("Bias, guideline, largest absolute", points("bias_artefact", function(v) max(abs(v))), b),
  row("Coverage, period, no capping", range_pct("cover_real"), uncapped),
  row("Coverage, period, 30% capped", range_pct("cover_real"), capped),
  row("Coverage, guideline, shift > 0, no capping", range_pct("cover_artefact"), uncapped & out$shift > 0),
  row("Coverage, guideline, shift > 0, 30% capped", range_pct("cover_artefact"), capped),
  row("Period shown, no period change", range_pct("period_shown"), b & out$delta == 0),
  row("Period shown, period change", range_pct("period_shown"), b & out$delta > 0),
  row("Guideline shown, shift 0", range_pct("guideline_shown"), b & out$shift == 0),
  row("Guideline shown, shift > 0", range_pct("guideline_shown"), b & out$shift > 0),
  row("Larger named correctly, shift > 0 and period change", range_pct("larger_correct"), b & out$shift > 0 & out$delta > 0),
  row("Larger named wrongly, any cell", range_pct("larger_wrong"), b),
  row("Values excluded, 30% capped", range_pct("excluded"), capped),
  "",
  "## Other cells", "",
  "| | n = 200 | n = 1000 | n = 5000 |", "|---|--:|--:|--:|",
  row("NS, bias of guideline, no capping", points("bias_artefact", function(v) max(abs(v))), out$set == "NS" & out$capped == 0),
  row("NS, bias of guideline, 30% capped (mean)", points("bias_artefact", d = 1), out$set == "NS" & out$capped > 0),
  row("NS, coverage of guideline, 30% capped", range_pct("cover_artefact"), out$set == "NS" & out$capped > 0),
  row("NS, values excluded, 30% capped", range_pct("excluded"), out$set == "NS" & out$capped > 0),
  row("Near, coverage of period", range_pct("cover_real"), out$set == "near"),
  row("Near, coverage of guideline, shift > 0", range_pct("cover_artefact"), out$set == "near" & out$shift > 0),
  row("Near, guideline shown, shift > 0", range_pct("guideline_shown"), out$set == "near" & out$shift > 0),
  row("Reported, compliant, coverage of guideline, shift > 0", range_pct("cover_artefact"), out$reported == "compliant" & out$shift > 0),
  row("Reported, compliant, guideline shown, shift 0", range_pct("guideline_shown"), out$reported == "compliant" & out$shift == 0),
  row("Reported, own rule, coverage of guideline", range_pct("cover_artefact"), out$reported == "own rule"),
  row("Reported, own rule, guideline shown", range_pct("guideline_shown"), out$reported == "own rule"),
  row("Clustered, coverage of period, isolate bootstrap", range_pct("cover_real"), out$set == "clustered"),
  row("Clustered, coverage of period, site bootstrap", range_pct("cover_real_site"), out$set == "clustered"),
  row("Clustered, period shown, no period change, isolate bootstrap", range_pct("period_shown"), out$set == "clustered" & out$delta == 0),
  row("Clustered, period shown, no period change, site bootstrap", range_pct("period_shown_site"), out$set == "clustered" & out$delta == 0),
  row("Clustered, coverage of guideline, shift 1, isolate bootstrap", range_pct("cover_artefact"), out$set == "clustered" & out$shift == 1),
  row("Clustered, coverage of guideline, shift 1, site bootstrap", range_pct("cover_artefact_site"), out$set == "clustered" & out$shift == 1),
  row("Clustered, larger named wrongly, isolate bootstrap", range_pct("larger_wrong"), out$set == "clustered"),
  row("Clustered, larger named wrongly, site bootstrap", range_pct("larger_wrong_site"), out$set == "clustered"),
  row("Low, coverage of period", range_pct("cover_real"), out$set == "low"),
  row("Low, coverage of guideline, shift 1", range_pct("cover_artefact"), out$set == "low" & out$shift == 1),
  row("Low, period shown, no period change", range_pct("period_shown"), out$set == "low" & out$delta == 0),
  row("Low, period shown, +1 point", range_pct("period_shown"), out$set == "low" & out$delta > 0),
  row("Low, guideline shown, shift 0", range_pct("guideline_shown"), out$set == "low" & out$shift == 0),
  row("Low, guideline shown, shift 1", range_pct("guideline_shown"), out$set == "low" & out$shift == 1),
  row("Low, larger named wrongly", range_pct("larger_wrong"), out$set == "low"),
  row("Straddle, not interpretable", range_pct("not_interpretable"), out$set == "straddle"),
  "",
  sprintf(paste("True components: period %s points with a period change (base cells); guideline %s points (shift 1)",
                "and %s (shift 2); reported own-rule guideline %s."),
          rng(out$true_real[b & out$delta > 0]), rng(out$true_artefact[b & out$shift == 1]),
          rng(out$true_artefact[b & out$shift == 2]), rng(out$true_artefact[out$reported == "own rule"])),
  "",
  sprintf("Base cell with the lowest period coverage (n = %d, prev %.2f, shift %d, delta %.2f, capped %.1f), rerun with 600 new datasets: coverage %s, period shown %s.",
          recheck$n, recheck$prev, recheck$shift, recheck$delta, recheck$capped, pct(recheck$cover_real), pct(recheck$period_shown)),
  "",
  "## Twenty pairs per dataset, Benjamini-Hochberg over them", "",
  "| | % |", "|---|--:|",
  paste0("| ", multi_summary$measure, " | ", pct(multi_summary$value), " |"),
  "",
  "## MIC-shift model", "",
  "| Case | Target | Estimated | Mean shift | Coverage | p < 0.05 |", "|---|--:|--:|--:|--:|--:|",
  sprintf("| %s | %.2f | %s | %.2f | %s | %s |", shift_sim$case, shift_sim$target, pct(shift_sim$estimated), shift_sim$mean_shift,
          pct(shift_sim$cover), pct(shift_sim$rejected)))
writeLines(lines, "results/simulation.md")
message(paste(lines, collapse = "\n"))
