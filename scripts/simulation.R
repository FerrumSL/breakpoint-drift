# Simulation: does the decomposition recover a known change at fixed criteria and a known
# reclassification, and how often does each verdict appear?
# Run from the project root: Rscript scripts/simulation.R  (about 2 hours on 16 cores)
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
#   low        rare resistance: resistant share 0.5% or 2%, change at fixed criteria 0 or +1 point
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
  fixed <- ((p[2, 2] - p[1, 2]) + (p[2, 1] - p[1, 1])) / 2
  observed <- p[2, 2] - p[1, 1]
  if (cell$reported == "own rule") observed <- p[2, 2] - p_at(cell$prev, rev$r[1] - 1, rev, wt)
  c(fixed = fixed, reclass = observed - fixed)
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
  out <- c(fixed = t$fixed, reclass = t$reclass, fixed_low = t$fixed_low, fixed_high = t$fixed_high,
           reclass_low = t$reclass_low, reclass_high = t$reclass_high,
           verdict = match(t$verdict, verdict_levels), larger = match(t$larger, c("reclassification", "fixed criteria", "neither")),
           excluded = 1 - (t$n_y1g2 + t$n_y2g2) / (t$tested_y1 + t$tested_y2))
  if (cell$clustered) {
    s <- drift_test(drift_bootstrap(data, B = B, outcome = cell$outcome, cluster = "site"))
    out <- c(out, fixed_low_site = s$fixed_low, fixed_high_site = s$fixed_high,
             reclass_low_site = s$reclass_low, reclass_high_site = s$reclass_high,
             verdict_site = match(s$verdict, verdict_levels),
             larger_site = match(s$larger, c("reclassification", "fixed criteria", "neither")))
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
reps <- as.integer(Sys.getenv("SIM_REPS", "1000"))  # SIM_REPS=2 for a quick check
B <- as.integer(Sys.getenv("SIM_B", "999"))
# 1,000 datasets per cell; Monte Carlo error on a 95% coverage is about 0.7 points

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
  has_g <- abs(tr[["reclass"]]) > 1e-9
  has_p <- abs(tr[["fixed"]]) > 1e-9
  true_larger <- if (!has_g && !has_p) NA else if (abs(tr[["reclass"]]) > abs(tr[["fixed"]])) "reclassification" else "fixed criteria"
  larger <- c("reclassification", "fixed criteria", "neither")[x$larger]
  wrong <- function(l) if (is.na(true_larger)) mean(l != "neither") else mean(l != true_larger & l != "neither")
  site <- if (cell$clustered) c(
    cover_fixed_site = mean(x$fixed_low_site <= tr[["fixed"]] & tr[["fixed"]] <= x$fixed_high_site),
    cover_reclass_site = mean(x$reclass_low_site <= tr[["reclass"]] & tr[["reclass"]] <= x$reclass_high_site),
    fixed_shown_site = mean(verdict_levels[x$verdict_site] %in% c("change at fixed criteria", "both")),
    larger_wrong_site = wrong(c("reclassification", "fixed criteria", "neither")[x$larger_site])) else
    c(cover_fixed_site = NA, cover_reclass_site = NA, fixed_shown_site = NA, larger_wrong_site = NA)
  data.frame(cell, true_fixed = tr[["fixed"]], true_reclass = tr[["reclass"]],
             bias_fixed = mean(x$fixed) - tr[["fixed"]], bias_reclass = mean(x$reclass) - tr[["reclass"]],
             cover_fixed = mean(x$fixed_low <= tr[["fixed"]] & tr[["fixed"]] <= x$fixed_high),
             cover_reclass = mean(x$reclass_low <= tr[["reclass"]] & tr[["reclass"]] <= x$reclass_high),
             t(site), excluded = mean(x$excluded),
             fixed_change = mean(v == "change at fixed criteria"), reclassification = mean(v == "reclassification"),
             both = mean(v == "both"), unattributed = mean(v == "unattributed"), no_change = mean(v == "no change detected"),
             reclass_shown = mean(v %in% c("reclassification", "both")),
             fixed_shown = mean(v %in% c("change at fixed criteria", "both")),
             larger_correct = if (is.na(true_larger)) NA else mean(larger == true_larger),
             larger_wrong = wrong(larger), not_interpretable = mean(v == "not interpretable"))
}))
dir.create("results", showWarnings = FALSE)
num <- vapply(out, is.numeric, logical(1))
out[num] <- lapply(out[num], signif, 4)
write.csv(out, "results/simulation.csv", row.names = FALSE)

# ---- The base cell with the lowest fixed-criteria coverage, again with 600 new datasets -----------------
# With 156 cells, some fall a few Monte Carlo errors from 95% by chance; this checks the worst.
worst <- which.min(ifelse(out$set == "base", out$cover_fixed, Inf))
recheck_runs <- do.call(rbind, mclapply(seq_len(600), function(k) {
  set.seed(90000 + k)
  one_run(grid[worst, ], B)
}, mc.cores = cores))
recheck <- data.frame(grid[worst, ], datasets = 600,
                      cover_fixed = mean(recheck_runs[, "fixed_low"] <= out$true_fixed[worst] & out$true_fixed[worst] <= recheck_runs[, "fixed_high"]),
                      fixed_shown = mean(verdict_levels[recheck_runs[, "verdict"]] %in% c("change at fixed criteria", "both")))
write.csv(recheck, "results/simulation_recheck.csv", row.names = FALSE)

# ---- Several pairs in one dataset, Benjamini-Hochberg over them ------------------------------------
# 2022 against 2023: gentamicin (revised, shift 1) and piperacillin/tazobactam (unchanged between
# these editions) for each of 10 Enterobacterales species, n = 1,000 per period and pair, resistant
# share 20%; half the species get a change at fixed criteria of +5 points in both drugs. True
# reclassification is 0 for every piperacillin/tazobactam pair, and the change at fixed criteria 0 for half.
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
  t$true_fixed <- t$mo %in% species[1:5]
  t$true_reclass <- t$ab == "GEN"
  t[c("mo", "ab", "verdict", "true_fixed", "true_reclass")]
}
multi <- mclapply(seq_len(reps), function(k) {
  set.seed(4026 + k)
  cbind(dataset = k, multi_run(B))
}, mc.cores = cores, mc.preschedule = FALSE)
multi <- do.call(rbind, multi)
r_shown <- multi$verdict %in% c("reclassification", "both")
f_shown <- multi$verdict %in% c("change at fixed criteria", "both")
fdr <- function(shown, true) mean(tapply(shown & !true, multi$dataset, sum) / pmax(1, tapply(shown, multi$dataset, sum)))
multi_summary <- data.frame(
  measure = c("reclassification shown, revised pairs", "reclassification shown, unrevised pairs", "fixed criteria shown, +5 points",
              "fixed criteria shown, no change at fixed criteria", "false discovery rate, reclassification", "false discovery rate, fixed criteria",
              "pairs not interpretable"),
  value = c(mean(r_shown[multi$true_reclass]), mean(r_shown[!multi$true_reclass]), mean(f_shown[multi$true_fixed]),
            mean(f_shown[!multi$true_fixed]), fdr(r_shown, multi$true_reclass), fdr(f_shown, multi$true_fixed),
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
  row("Bias, fixed criteria, largest absolute", points("bias_fixed", function(v) max(abs(v))), b),
  row("Bias, reclassification, largest absolute", points("bias_reclass", function(v) max(abs(v))), b),
  row("Coverage, fixed criteria, no capping", range_pct("cover_fixed"), uncapped),
  row("Coverage, fixed criteria, 30% capped", range_pct("cover_fixed"), capped),
  row("Coverage, reclassification, shift > 0, no capping", range_pct("cover_reclass"), uncapped & out$shift > 0),
  row("Coverage, reclassification, shift > 0, 30% capped", range_pct("cover_reclass"), capped),
  row("Fixed criteria shown, no change at fixed criteria", range_pct("fixed_shown"), b & out$delta == 0),
  row("Fixed criteria shown, change at fixed criteria", range_pct("fixed_shown"), b & out$delta > 0),
  row("Reclassification shown, shift 0", range_pct("reclass_shown"), b & out$shift == 0),
  row("Reclassification shown, shift > 0", range_pct("reclass_shown"), b & out$shift > 0),
  row("Larger named correctly, shift > 0 and change at fixed criteria", range_pct("larger_correct"), b & out$shift > 0 & out$delta > 0),
  row("Larger named wrongly, any cell", range_pct("larger_wrong"), b),
  row("Values excluded, 30% capped", range_pct("excluded"), capped),
  "",
  "## Other cells", "",
  "| | n = 200 | n = 1000 | n = 5000 |", "|---|--:|--:|--:|",
  row("NS, bias of reclassification, no capping", points("bias_reclass", function(v) max(abs(v))), out$set == "NS" & out$capped == 0),
  row("NS, bias of reclassification, 30% capped (mean)", points("bias_reclass", d = 1), out$set == "NS" & out$capped > 0),
  row("NS, coverage of reclassification, 30% capped", range_pct("cover_reclass"), out$set == "NS" & out$capped > 0),
  row("NS, values excluded, 30% capped", range_pct("excluded"), out$set == "NS" & out$capped > 0),
  row("Near, coverage of fixed criteria", range_pct("cover_fixed"), out$set == "near"),
  row("Near, coverage of reclassification, shift > 0", range_pct("cover_reclass"), out$set == "near" & out$shift > 0),
  row("Near, reclassification shown, shift > 0", range_pct("reclass_shown"), out$set == "near" & out$shift > 0),
  row("Reported, compliant, coverage of reclassification, shift > 0", range_pct("cover_reclass"), out$reported == "compliant" & out$shift > 0),
  row("Reported, compliant, reclassification shown, shift 0", range_pct("reclass_shown"), out$reported == "compliant" & out$shift == 0),
  row("Reported, own rule, coverage of reclassification", range_pct("cover_reclass"), out$reported == "own rule"),
  row("Reported, own rule, reclassification shown", range_pct("reclass_shown"), out$reported == "own rule"),
  row("Clustered, coverage of fixed criteria, isolate bootstrap", range_pct("cover_fixed"), out$set == "clustered"),
  row("Clustered, coverage of fixed criteria, site bootstrap", range_pct("cover_fixed_site"), out$set == "clustered"),
  row("Clustered, fixed criteria shown, no change at fixed criteria, isolate bootstrap", range_pct("fixed_shown"), out$set == "clustered" & out$delta == 0),
  row("Clustered, fixed criteria shown, no change at fixed criteria, site bootstrap", range_pct("fixed_shown_site"), out$set == "clustered" & out$delta == 0),
  row("Clustered, coverage of reclassification, shift 1, isolate bootstrap", range_pct("cover_reclass"), out$set == "clustered" & out$shift == 1),
  row("Clustered, coverage of reclassification, shift 1, site bootstrap", range_pct("cover_reclass_site"), out$set == "clustered" & out$shift == 1),
  row("Clustered, larger named wrongly, isolate bootstrap", range_pct("larger_wrong"), out$set == "clustered"),
  row("Clustered, larger named wrongly, site bootstrap", range_pct("larger_wrong_site"), out$set == "clustered"),
  row("Low, coverage of fixed criteria", range_pct("cover_fixed"), out$set == "low"),
  row("Low, coverage of reclassification, shift 1", range_pct("cover_reclass"), out$set == "low" & out$shift == 1),
  row("Low, fixed criteria shown, no change at fixed criteria", range_pct("fixed_shown"), out$set == "low" & out$delta == 0),
  row("Low, fixed criteria shown, +1 point", range_pct("fixed_shown"), out$set == "low" & out$delta > 0),
  row("Low, reclassification shown, shift 0", range_pct("reclass_shown"), out$set == "low" & out$shift == 0),
  row("Low, reclassification shown, shift 1", range_pct("reclass_shown"), out$set == "low" & out$shift == 1),
  row("Low, larger named wrongly", range_pct("larger_wrong"), out$set == "low"),
  row("Straddle, not interpretable", range_pct("not_interpretable"), out$set == "straddle"),
  "",
  sprintf(paste("True components: fixed criteria %s points with a change at fixed criteria (base cells); reclassification %s points (shift 1)",
                "and %s (shift 2); reported own-rule reclassification %s."),
          rng(out$true_fixed[b & out$delta > 0]), rng(out$true_reclass[b & out$shift == 1]),
          rng(out$true_reclass[b & out$shift == 2]), rng(out$true_reclass[out$reported == "own rule"])),
  "",
  sprintf("Base cell with the lowest fixed-criteria coverage (n = %d, prev %.2f, shift %d, delta %.2f, capped %.1f), rerun with 600 new datasets: coverage %s, fixed criteria shown %s.",
          recheck$n, recheck$prev, recheck$shift, recheck$delta, recheck$capped, pct(recheck$cover_fixed), pct(recheck$fixed_shown)),
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
