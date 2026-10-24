# Split a change in %R between two years into the part caused by a breakpoint revision and the
# part caused by the MICs themselves.
#
# p(y, g) is the share of isolates from year y that are R under guideline g. With y1 < y2 and
# g1, g2 the guidelines in force in those years:
#   observed change = p(y2, g2) - p(y1, g1)
#   real change     = p(y2, g2) - p(y1, g2)   (same guideline, MICs from two years)
#   artefact        = p(y1, g2) - p(y1, g1)   (same MICs, two guidelines)
# The split depends on the order of the two steps. Holding the old guideline instead gives
#   real change     = p(y2, g1) - p(y1, g1)
# The reported real change is the mean of the two (the Shapley value), and the artefact is the
# rest, so the two parts always add up to the observed change.
#
# A capped MIC can fall on both sides of a new breakpoint: "<=1" cannot be read against S <= 0.25.
# AMR::as.sir() returns NI (not interpretable) for those, and they are left out of p(y, g). The
# coverage columns say how many tested isolates each p(y, g) still rests on.


# Guideline edition in force in a given year
#
# Assumes a lab uses the edition named for the year, `lag` years late. Years before the oldest
# edition in the AMR package give NA; years after the newest give the newest.
guideline_in_force <- function(year, org = "CLSI", lag = 0) {
  editions <- unique(AMR::clinical_breakpoints$guideline)
  editions <- as.integer(sub(".* ", "", editions[startsWith(editions, paste0(org, " "))]))
  if (length(editions) == 0) stop("The AMR package has no ", org, " breakpoints.")
  year <- year - lag
  ifelse(year < min(editions), NA_character_, paste(org, pmin(year, max(editions))))
}


# Decompose the change in %R between two years, or two periods of years
#
# `data` has one row per isolate and drug, with columns `year`, `mo` (organism), `ab` (drug) and
# `mic` (raw MIC, e.g. "<=0.25"). Deduplicate it first, e.g. with AMR::first_isolate().
# `from` and `to` are a year each or a range such as 2013:2016. A period is read with the edition
# in force in its last year, as a cumulative antibiogram for that period would have been.
# `in_force` maps a year to a guideline name such as "CLSI 2019".
drift_decompose <- function(data, from = min(data$year), to = max(data$year),
                            in_force = guideline_in_force) {
  missing <- setdiff(c("year", "mo", "ab", "mic"), names(data))
  if (length(missing) > 0) stop("`data` lacks column(s): ", paste(missing, collapse = ", "))
  stopifnot(max(from) < min(to))
  g <- c(in_force(max(from)), in_force(max(to)))
  if (anyNA(g)) stop("No breakpoints in the AMR package for ", c(max(from), max(to))[is.na(g)][1], ".")
  label <- function(w) if (length(w) == 1) as.character(w) else paste0(min(w), "-", max(w))
  periods <- list(from, to)

  data <- data[data$year %in% c(from, to), ]
  data$mo <- as.character(AMR::as.mo(data$mo, info = FALSE))
  data$ab <- as.character(AMR::as.ab(data$ab, info = FALSE))
  data$mic <- AMR::as.mic(data$mic)
  out <- unique(data[c("mo", "ab")])
  out <- out[order(out$mo, out$ab), ]
  rownames(out) <- NULL
  key <- paste(out$mo, out$ab)
  out$year_from <- label(from)
  out$year_to <- label(to)
  out$guideline_from <- g[1]
  out$guideline_to <- g[2]

  for (y in 1:2) {
    x <- data[data$year %in% periods[[y]], ]
    sum_by_pair <- function(v) as.vector(tapply(v, factor(paste(x$mo, x$ab), levels = key), sum, default = 0))
    out[[paste0("tested_y", y)]] <- sum_by_pair(rep(1L, nrow(x)))
    for (k in 1:2) {
      sir <- AMR::as.sir(x$mic, mo = x$mo, ab = x$ab, guideline = g[k],
                         capped_mic_handling = "standard", info = FALSE)
      cell <- paste0("y", y, "g", k)
      n <- sum_by_pair(!is.na(sir) & sir != "NI")
      r <- sum_by_pair(sir %in% "R")
      out[[paste0("n_", cell)]] <- n
      out[[paste0("r_", cell)]] <- r
      out[[paste0("p_", cell)]] <- ifelse(n > 0, r / n, NA_real_)
    }
  }

  out$observed <- out$p_y2g2 - out$p_y1g1
  out$real_at_new <- out$p_y2g2 - out$p_y1g2
  out$real_at_old <- out$p_y2g1 - out$p_y1g1
  out$real <- (out$real_at_new + out$real_at_old) / 2
  out$artefact <- out$observed - out$real
  out
}


# Test each organism x drug pair and label it
#
# Two Fisher's exact tests per pair, each with a Benjamini-Hochberg correction across pairs:
# the observed comparison (each year under its own guideline, what a naive antibiogram shows)
# and the real one (both years under the newer guideline). A pair is a "significant change" if
# the real test is significant, a "breakpoint artefact" if only the observed one is and the
# artefact is larger than the real change, and "noise" otherwise. Clopper-Pearson 95% intervals are added for both years under the newer guideline.
drift_test <- function(x, alpha = 0.05) {
  fisher <- function(r1, n1, r2, n2) mapply(function(r1, n1, r2, n2) {
    if (n1 == 0 || n2 == 0) return(NA_real_)
    stats::fisher.test(matrix(c(r1, n1 - r1, r2, n2 - r2), 2))$p.value
  }, r1, n1, r2, n2)
  x$p_observed <- stats::p.adjust(fisher(x$r_y1g1, x$n_y1g1, x$r_y2g2, x$n_y2g2), "BH")
  x$p_real <- stats::p.adjust(fisher(x$r_y1g2, x$n_y1g2, x$r_y2g2, x$n_y2g2), "BH")

  for (cell in c("y1g2", "y2g2")) {
    ci <- mapply(function(r, n) if (n == 0) c(NA, NA) else stats::binom.test(r, n)$conf.int,
                 x[[paste0("r_", cell)]], x[[paste0("n_", cell)]])
    x[[paste0("ci_low_", cell)]] <- ci[1, ]
    x[[paste0("ci_high_", cell)]] <- ci[2, ]
  }

  sig <- function(p) !is.na(p) & p < alpha
  x$verdict <- ifelse(is.na(x$p_real), "not interpretable",
               ifelse(sig(x$p_real), "significant change",
               ifelse(sig(x$p_observed) & abs(x$artefact) > abs(x$real), "breakpoint artefact", "noise")))
  x
}


# Write the stewardship report as Markdown, one section per organism
drift_report <- function(x, path, title = "Breakpoint drift report", note = NULL) {
  pct <- function(p) ifelse(is.na(p), "", sprintf("%.1f", 100 * p))
  pp <- function(p) ifelse(is.na(p), "", sprintf("%+.1f", 100 * p))
  g1 <- x$guideline_from[1]
  g2 <- x$guideline_to[1]
  y1 <- x$year_from[1]
  y2 <- x$year_to[1]
  counts <- table(factor(x$verdict, c("significant change", "breakpoint artefact", "noise", "not interpretable")))

  lines <- c(
    paste("#", title), "",
    sprintf("%s isolates read with %s, %s isolates with %s. Changes are in percentage points of %%R.",
            y1, g1, y2, g2),
    if (!is.null(note)) c("", note),
    "",
    "- **Observed**: %R in each year under its own guideline, as an antibiogram reports it.",
    "- **Real**: change in the MICs, with the guideline held fixed (mean of holding either one).",
    "- **Artefact**: the rest of the observed change, caused by the breakpoint revision alone.",
    sprintf("- **Verdict**: Fisher's exact test, Benjamini-Hochberg corrected across %d pairs.", nrow(x)),
    sprintf("- **Read**: %s isolates that %s / %s could classify, of those tested. A capped MIC such as <=1 cannot be read against a breakpoint below it.", y1, g1, g2),
    "",
    paste0("Verdicts: ", paste(sprintf("%s %d", names(counts), counts), collapse = ", "), "."),
    ""
  )
  for (mo in unique(x$mo)) {
    s <- x[x$mo == mo, ]
    s <- s[order(-abs(s$observed)), ]
    lines <- c(lines, paste0("## *", AMR::mo_name(mo, info = FALSE), "*"), "",
      sprintf("| Drug | %%R %s | %%R %s | Observed | Real | Artefact | Verdict | Read %s |", y1, y2, y1),
      "|---|--:|--:|--:|--:|--:|---|--:|",
      sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |",
              AMR::ab_name(s$ab, info = FALSE), pct(s$p_y1g1), pct(s$p_y2g2), pp(s$observed),
              pp(s$real), pp(s$artefact), s$verdict,
              sprintf("%d / %d of %d", s$n_y1g1, s$n_y1g2, s$tested_y1)),
      "")
  }
  writeLines(lines, path)
  invisible(path)
}
