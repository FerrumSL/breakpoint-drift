# Split a change in %R between two years into the part caused by a breakpoint revision and the
# part caused by the MICs or zone diameters themselves.
#
# p(y, g) is the share of isolates from year y that are R under guideline g. With y1 < y2 and
# g1, g2 the guidelines in force in those years:
#   observed change = p(y2, g2) - p(y1, g1)
#   real change     = p(y2, g2) - p(y1, g2)   (same guideline, values from two years)
#   artefact        = p(y1, g2) - p(y1, g1)   (same values, two guidelines)
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
# a raw value in `mic` (e.g. "<=0.25") or `disk` (zone diameter in mm), or both columns with one
# filled per row. An optional `org` column ("CLSI" or "EUCAST") says which organisation's
# guideline each row was read with; without it every row is CLSI. Deduplicate first, e.g. with
# AMR::first_isolate().
# `from` and `to` are a year each or a range such as 2013:2016. A period is read with the edition
# in force in its last year, as a cumulative antibiogram for that period would have been.
# `in_force(year, org)` returns a guideline name such as "CLSI 2019".
#
# One row per organism, drug, method (MIC or disk) and organisation: breakpoints for MICs and
# zones are revised separately, and a lab switching method is a change of its own.
drift_decompose <- function(data, from = min(data$year), to = max(data$year),
                            in_force = guideline_in_force) {
  missing <- setdiff(c("year", "mo", "ab"), names(data))
  if (!any(c("mic", "disk") %in% names(data))) missing <- c(missing, "mic or disk")
  if (length(missing) > 0) stop("`data` lacks column(s): ", paste(missing, collapse = ", "))
  stopifnot(max(from) < min(to))
  for (col in c("mic", "disk")) if (is.null(data[[col]])) data[[col]] <- NA
  if (is.null(data$org)) data$org <- "CLSI"
  if (any(!is.na(data$mic) & !is.na(data$disk))) stop("A row has both an MIC and a zone diameter.")
  data$method <- ifelse(!is.na(data$mic), "MIC", "disk")

  data <- data[data$year %in% c(from, to) & !(is.na(data$mic) & is.na(data$disk)), ]
  data$mo <- as.character(AMR::as.mo(data$mo, info = FALSE))
  data$ab <- as.character(AMR::as.ab(data$ab, info = FALSE))
  out <- do.call(rbind, lapply(split(data, data$org), decompose_org, from, to, in_force))
  rownames(out) <- NULL

  out$observed <- out$p_y2g2 - out$p_y1g1
  out$real_at_new <- out$p_y2g2 - out$p_y1g2
  out$real_at_old <- out$p_y2g1 - out$p_y1g1
  out$real <- (out$real_at_new + out$real_at_old) / 2
  out$artefact <- out$observed - out$real
  out
}

decompose_org <- function(data, from, to, in_force) {
  org <- data$org[1]
  g <- c(in_force(max(from), org), in_force(max(to), org))
  if (anyNA(g)) stop("No ", org, " breakpoints in the AMR package for ",
                     c(max(from), max(to))[is.na(g)][1], ".")
  label <- function(w) if (length(w) == 1) as.character(w) else paste0(min(w), "-", max(w))
  periods <- list(from, to)

  out <- unique(data[c("mo", "ab", "method", "org")])
  out <- out[order(out$mo, out$ab, out$method), ]
  key <- paste(out$mo, out$ab, out$method)
  out$year_from <- label(from)
  out$year_to <- label(to)
  out$guideline_from <- g[1]
  out$guideline_to <- g[2]

  for (y in 1:2) {
    x <- data[data$year %in% periods[[y]], ]
    sum_by_pair <- function(v) {
      as.vector(tapply(v, factor(paste(x$mo, x$ab, x$method), levels = key), sum, default = 0))
    }
    out[[paste0("tested_y", y)]] <- sum_by_pair(rep(1L, nrow(x)))
    for (k in 1:2) {
      sir <- read_sir(x, g[k])
      cell <- paste0("y", y, "g", k)
      n <- sum_by_pair(!is.na(sir) & sir != "NI")
      r <- sum_by_pair(sir %in% "R")
      out[[paste0("n_", cell)]] <- n
      out[[paste0("r_", cell)]] <- r
      out[[paste0("p_", cell)]] <- ifelse(n > 0, r / n, NA_real_)
    }
  }
  out
}

# S/I/R/NI as text for each row, from its MIC or zone diameter
read_sir <- function(x, guideline) {
  sir <- rep(NA_character_, nrow(x))
  for (method in c("MIC", "disk")) {
    i <- x$method == method
    if (!any(i)) next
    sir[i] <- as.character(if (method == "MIC") {
      AMR::as.sir(AMR::as.mic(x$mic[i]), mo = x$mo[i], ab = x$ab[i], guideline = guideline,
                  capped_mic_handling = "standard", info = FALSE)
    } else {
      AMR::as.sir(AMR::as.disk(x$disk[i]), mo = x$mo[i], ab = x$ab[i], guideline = guideline, info = FALSE)
    })
  }
  sir
}


# Test each organism x drug pair and label it
#
# Two Fisher's exact tests per pair, each with a Benjamini-Hochberg correction across pairs:
# the observed comparison (each year under its own guideline, what a naive antibiogram shows)
# and the real one (both years under the newer guideline). A pair is a "significant change" if
# the real test is significant, a "breakpoint artefact" if only the observed one is and the
# artefact is larger than the real change, and "noise" otherwise. Clopper-Pearson 95% intervals
# are added for both years under the newer guideline.
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


# Shift in the MIC distribution between two periods, in doubling dilutions
#
# Each stored MIC is an interval on the log2 scale: "8" means the true MIC is above 4 and at most
# 8, "<=0.25" at most 0.25, ">=16" (or ">8") above 8. For each organism and drug, a normal model
# of log2 MIC with a period effect is fitted to these intervals by maximum likelihood
# (survival::survreg). The shift does not use any breakpoint, so it also catches MICs creeping
# up below one, before %R moves. p-values are Benjamini-Hochberg corrected across pairs.
# A period whose values are nearly all capped ("<=1" and nothing else) has no estimable mean, so
# each period needs `min_n` values, of which `min_exact` are not capped; otherwise the shift is NA.
# ponytail: one normal per period. A wild-type / non-wild-type mixture (PLOS ONE 2019,
# doi:10.1371/journal.pone.0220427) fits bimodal populations better; upgrade if a resistant
# subpopulation shifting on its own matters.
drift_mic_shift <- function(data, from = min(data$year), to = max(data$year), min_n = 10,
                            min_exact = 5) {
  data <- data[data$year %in% c(from, to) & !is.na(data$mic), ]
  data$mo <- as.character(AMR::as.mo(data$mo, info = FALSE))
  data$ab <- as.character(AMR::as.ab(data$ab, info = FALSE))
  data$late <- data$year %in% to
  data <- cbind(data, mic_interval(data$mic))
  data <- data[!(is.na(data$lo) & is.na(data$hi)), ]

  out <- unique(data[c("mo", "ab")])
  out <- out[order(out$mo, out$ab), ]
  rownames(out) <- NULL
  fits <- lapply(seq_len(nrow(out)), function(i) {
    fit_shift(data[data$mo == out$mo[i] & data$ab == out$ab[i], ], min_n, min_exact)
  })
  out <- cbind(out, do.call(rbind, fits))
  out$p_shift <- stats::p.adjust(out$p_shift, "BH")
  out
}

# log2 bounds of the true MIC behind each stored value; NA is an open end
mic_interval <- function(mic) {
  s <- as.character(AMR::as.mic(mic))
  v <- round(log2(as.numeric(sub("^[<>=]+", "", s))))  # 2^-5 is written 0.03 or 0.032
  sign <- sub("^([<>=]*).*", "\\1", s)
  data.frame(lo = ifelse(sign %in% c("<", "<="), NA, ifelse(sign == ">", v, v - 1)),
             hi = ifelse(sign %in% c(">", ">="), NA, ifelse(sign == "<", v - 1, v)))
}

fit_shift <- function(x, min_n, min_exact) {
  n <- c(sum(!x$late), sum(x$late))
  exact <- !is.na(x$lo) & !is.na(x$hi)
  n_exact <- c(sum(exact & !x$late), sum(exact & x$late))
  res <- data.frame(n_shift_from = n[1], n_shift_to = n[2], gm_from = NA_real_, gm_to = NA_real_,
                    shift = NA_real_, shift_low = NA_real_, shift_high = NA_real_, p_shift = NA_real_)
  if (min(n) < min_n || min(n_exact) < min_exact) return(res)
  # A warning here means no finite estimate, e.g. every value in a period censored the same way
  fit <- tryCatch(survival::survreg(survival::Surv(lo, hi, type = "interval2") ~ late, data = x,
                                    dist = "gaussian"),
                  error = function(e) NULL, warning = function(w) NULL)
  if (is.null(fit)) return(res)
  b <- stats::coef(fit)
  se <- sqrt(stats::vcov(fit)["lateTRUE", "lateTRUE"])
  if (!is.finite(se)) return(res)
  res$gm_from <- 2^b[[1]]
  res$gm_to <- 2^(b[[1]] + b[[2]])
  res$shift <- b[[2]]
  res$shift_low <- b[[2]] - stats::qnorm(0.975) * se
  res$shift_high <- b[[2]] + stats::qnorm(0.975) * se
  res$p_shift <- 2 * stats::pnorm(-abs(b[[2]] / se))
  res
}


# Write the stewardship report as Markdown, one section per organism
#
# `x` comes from drift_test(). If it also has the columns of drift_mic_shift() (merged by `mo`
# and `ab`, MIC rows only), the report adds the MIC shift.
drift_report <- function(x, path, title = "Breakpoint drift report", note = NULL) {
  pct <- function(p) ifelse(is.na(p), "", sprintf("%.1f", 100 * p))
  pp <- function(p) ifelse(is.na(p), "", sprintf("%+.1f", 100 * p))
  y1 <- x$year_from[1]
  y2 <- x$year_to[1]
  editions <- unique(x[c("guideline_from", "guideline_to")])
  counts <- table(factor(x$verdict, c("significant change", "breakpoint artefact", "noise", "not interpretable")))
  has_shift <- !is.null(x$shift)

  lines <- c(
    paste("#", title), "",
    sprintf("%s isolates read with %s; %s isolates with %s. Changes are in percentage points of %%R.",
            y1, paste(editions$guideline_from, collapse = " or "),
            y2, paste(editions$guideline_to, collapse = " or ")),
    if (!is.null(note)) c("", note),
    "",
    "- **Observed**: %R in each period under its own guideline, as an antibiogram reports it.",
    "- **Real**: change in the values, with the guideline held fixed (mean of holding either one).",
    "- **Artefact**: the rest of the observed change, caused by the breakpoint revision alone.",
    sprintf("- **Verdict**: Fisher's exact test, Benjamini-Hochberg corrected across %d pairs.", nrow(x)),
    sprintf("- **Read**: %s isolates that the old / new guideline could classify, of those tested. A capped MIC such as <=1 cannot be read against a breakpoint below it.", y1),
    if (has_shift) "- **MIC shift**: change in the MIC distribution in doubling dilutions, with 95% CI, from an interval-censored model that uses no breakpoint. * = significant after Benjamini-Hochberg correction.",
    "",
    paste0("Verdicts: ", paste(sprintf("%s %d", names(counts), counts), collapse = ", "), "."),
    ""
  )
  for (mo in unique(x$mo)) {
    s <- x[x$mo == mo, ]
    s <- s[order(is.na(s$observed), -abs(s$observed)), ]
    shift <- if (!has_shift) NULL else
      ifelse(is.na(s$shift), "", sprintf("%+.1f (%+.1f to %+.1f)%s", s$shift, s$shift_low, s$shift_high,
                                         ifelse(!is.na(s$p_shift) & s$p_shift < 0.05, " *", "")))
    lines <- c(lines, paste0("## *", AMR::mo_name(mo, info = FALSE), "*"), "",
      paste0(sprintf("| Drug | Test | %%R %s | %%R %s | Observed | Real | Artefact | Verdict | Read %s |",
                     y1, y2, y1), if (has_shift) " MIC shift |"),
      paste0("|---|---|--:|--:|--:|--:|--:|---|--:|", if (has_shift) "--:|"),
      paste0(sprintf("| %s | %s %s | %s | %s | %s | %s | %s | %s | %s |",
                     AMR::ab_name(s$ab, info = FALSE), s$method, s$org, pct(s$p_y1g1), pct(s$p_y2g2),
                     pp(s$observed), pp(s$real), pp(s$artefact), s$verdict,
                     sprintf("%d / %d of %d", s$n_y1g1, s$n_y1g2, s$tested_y1)),
             if (has_shift) paste0(" ", shift, " |")),
      "")
  }
  writeLines(lines, path)
  invisible(path)
}
