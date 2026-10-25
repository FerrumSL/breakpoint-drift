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


#' Guideline edition in force in a given year
#'
#' Assumes a lab uses the edition named for the year, `lag` years late. Years before the oldest
#' edition in the AMR package give NA; years after the newest give the newest.
#'
#' @param year Year(s).
#' @param org `"CLSI"` or `"EUCAST"`.
#' @param lag Years a lab takes to adopt a new edition.
#' @return Character vector of editions, such as `"CLSI 2023"`.
#' @examples
#' guideline_in_force(2021:2024)
#' guideline_in_force(2024, "EUCAST", lag = 1)
#' @export
guideline_in_force <- function(year, org = "CLSI", lag = 0) {
  editions <- unique(AMR::clinical_breakpoints$guideline)
  editions <- as.integer(sub(".* ", "", editions[startsWith(editions, paste0(org, " "))]))
  if (length(editions) == 0) stop("The AMR package has no ", org, " breakpoints.")
  year <- year - lag
  ifelse(year < min(editions), NA_character_, paste(org, pmin(year, max(editions))))
}


#' Edition each lab had adopted by a given year
#'
#' Looks up a table of adoption dates. A lab may switch organisation, e.g. CLSI to EUCAST.
#'
#' @param adoption Data frame with columns `year` and `guideline` (`"CLSI 2019"`), and `lab` if
#'   labs differ: from `year` on, that lab read with `guideline`.
#' @param year Year(s) to look up.
#' @param lab Lab of each year, if `adoption` has a `lab` column.
#' @return Character vector of editions; NA for a lab or year the table does not cover.
#' @examples
#' adoption <- data.frame(year = c(2015, 2021), guideline = c("CLSI 2015", "EUCAST 2021"))
#' adopted_edition(adoption, c(2014, 2018, 2023))
#' @export
adopted_edition <- function(adoption, year, lab = NA) {
  vapply(seq_along(year), function(i) {
    a <- adoption[adoption$year <= year[i], ]
    if (!is.null(adoption$lab)) a <- a[a$lab %in% rep_len(lab, length(year))[i], ]
    if (nrow(a) == 0) NA_character_ else a$guideline[which.max(a$year)]
  }, character(1))
}


#' Decompose the change in %R between two years, or two periods of years
#'
#' Splits the observed change in %R for each organism and drug into the real change (the values
#' moved, guideline held fixed) and the breakpoint artefact (the rest). The real change is the
#' mean of holding either guideline fixed (the Shapley value), so the two parts add up to the
#' observed change. Capped MICs that a guideline cannot read (NI) are left out; the `n_*`
#' columns say how many isolates each share rests on.
#'
#' Which edition was in force comes from one of two places:
#' - `adoption`, a table of when each lab adopted which edition (see [adopted_edition()]);
#'   `data` then needs a `lab` column if the table has one. Labs on different editions are
#'   pooled.
#' - otherwise the rule `in_force(year, org)`, by default the edition named for the year, with
#'   an optional `org` column (`"CLSI"` or `"EUCAST"`, default CLSI) per row.
#'
#' @param data One row per isolate and drug, with columns `year`, `mo` (organism), `ab` (drug)
#'   and a raw value in `mic` (e.g. `"<=0.25"`) or `disk` (zone diameter in mm), or both columns
#'   with one filled per row. Deduplicate first, e.g. with [AMR::first_isolate()]. An optional
#'   `reported` column holds the S/I/R each isolate was originally reported with; the observed
#'   change is then the change in reported %R, as the published antibiograms showed it, and the
#'   artefact also covers a lab's own reporting rules.
#' @param from,to A year each or a range such as `2013:2016`. A period is read with the edition
#'   in force at its last year, as a cumulative antibiogram for that period would have been.
#' @param in_force Function of (year, org) giving the edition in force.
#' @param adoption Optional table of adoption dates; see [adopted_edition()].
#' @return A data frame with one row per organism, drug, method (MIC or disk) and organisation:
#'   breakpoints for MICs and zones are revised separately, and a lab that switched
#'   organisation between the periods gets its own row (`"CLSI to EUCAST"`). Columns
#'   `p_y{1,2}g{1,2}` hold %R for each period under each guideline, then `observed`, `real`
#'   (with `real_at_new` and `real_at_old`) and `artefact`, as shares.
#' @examples
#' # CLSI 2023 lowered the Enterobacterales gentamicin breakpoints; the MICs did not change
#' mics <- rep(c("1", "2", "8", "16"), c(40, 20, 30, 10))
#' isolates <- data.frame(year = rep(c(2022, 2023), each = 100), mo = "B_ESCHR_COLI",
#'                        ab = "GEN", mic = mics)
#' d <- drift_decompose(isolates)
#' d[c("observed", "real", "artefact")]
#' @export
drift_decompose <- function(data, from = min(data$year), to = max(data$year),
                            in_force = guideline_in_force, adoption = NULL) {
  missing <- setdiff(c("year", "mo", "ab", if (!is.null(adoption$lab)) "lab"), names(data))
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

  # The two editions each row is read with
  if (is.null(adoption)) {
    for (k in 1:2) {
      year <- max(list(from, to)[[k]])
      g <- vapply(unique(data$org), function(o) in_force(year, o), character(1))
      if (anyNA(g)) stop("No ", names(g)[is.na(g)][1], " breakpoints in the AMR package for ", year, ".")
      data[[paste0("g", k)]] <- unname(g[data$org])
    }
  } else {
    labs <- if (is.null(adoption$lab)) NA else data$lab
    key <- unique(data.frame(lab = labs))
    for (k in 1:2) {
      year <- max(list(from, to)[[k]])
      g <- adopted_edition(adoption, rep(year, nrow(key)), key$lab)
      if (anyNA(g)) stop("The adoption table has no edition for lab ", key$lab[is.na(g)][1],
                         " by ", year, ".")
      data[[paste0("g", k)]] <- g[match(labs, key$lab)]
    }
    org <- function(g) sub(" .*", "", g)
    data$org <- ifelse(org(data$g1) == org(data$g2), org(data$g2), paste(org(data$g1), "to", org(data$g2)))
  }

  out <- unique(data[c("mo", "ab", "method", "org")])
  out <- out[order(out$mo, out$ab, out$method, out$org), ]
  rownames(out) <- NULL
  key <- paste(out$mo, out$ab, out$method, out$org)
  pair <- factor(paste(data$mo, data$ab, data$method, data$org), levels = key)
  label <- function(w) if (length(w) == 1) as.character(w) else paste0(min(w), "-", max(w))
  out$year_from <- label(from)
  out$year_to <- label(to)
  editions <- function(g) as.vector(tapply(g, pair, function(x) paste(sort(unique(x)), collapse = " + ")))
  out$guideline_from <- editions(data$g1)
  out$guideline_to <- editions(data$g2)

  for (y in 1:2) {
    in_period <- data$year %in% list(from, to)[[y]]
    x <- data[in_period, ]
    sum_by_pair <- function(v) as.vector(tapply(v, pair[in_period], sum, default = 0))
    out[[paste0("tested_y", y)]] <- sum_by_pair(rep(1L, nrow(x)))
    for (k in 1:2) {
      g <- x[[paste0("g", k)]]
      sir <- rep(NA_character_, nrow(x))
      for (edition in unique(g)) sir[g == edition] <- read_sir(x[g == edition, ], edition)
      cell <- paste0("y", y, "g", k)
      n <- sum_by_pair(!is.na(sir) & sir != "NI")
      r <- sum_by_pair(sir %in% "R")
      out[[paste0("n_", cell)]] <- n
      out[[paste0("r_", cell)]] <- r
      out[[paste0("p_", cell)]] <- ifelse(n > 0, r / n, NA_real_)
    }
    if (!is.null(x$reported)) {
      n <- sum_by_pair(x$reported %in% c("S", "SDD", "I", "R"))
      r <- sum_by_pair(x$reported %in% "R")
      out[[paste0("n_rep_y", y)]] <- n
      out[[paste0("r_rep_y", y)]] <- r
      out[[paste0("p_rep_y", y)]] <- ifelse(n > 0, r / n, NA_real_)
    }
  }

  out$observed <- if (is.null(data$reported)) out$p_y2g2 - out$p_y1g1 else out$p_rep_y2 - out$p_rep_y1
  out$real_at_new <- out$p_y2g2 - out$p_y1g2
  out$real_at_old <- out$p_y2g1 - out$p_y1g1
  out$real <- (out$real_at_new + out$real_at_old) / 2
  out$artefact <- out$observed - out$real
  out
}

# S/I/R/NI as text for each row, from its MIC or zone diameter
#
# Where a guideline has no breakpoint for an organism's own group, AMR::as.sir() falls back to a
# broader one. For Enterobacterales that broader one is the generic Gram-negative row (CLSI
# Table 2B-5, "other non-Enterobacterales"), which is wrong for them: in AMR 3.1.0, CLSI 2011-2019
# lack Enterobacterales rows for many drugs, and E. coli ceftriaxone would be read with S <= 8
# instead of S <= 1. Such readings are dropped (NA), so they count as not readable.
read_sir <- function(x, guideline) {
  sir <- rep(NA_character_, nrow(x))
  for (method in c("MIC", "disk")) {
    i <- which(x$method == method)
    if (length(i) == 0) next
    invisible(AMR::sir_interpretation_history(clean = TRUE))
    sir[i] <- as.character(if (method == "MIC") {
      AMR::as.sir(AMR::as.mic(x$mic[i]), mo = x$mo[i], ab = x$ab[i], guideline = guideline,
                  capped_mic_handling = "standard", info = FALSE)
    } else {
      AMR::as.sir(AMR::as.disk(x$disk[i]), mo = x$mo[i], ab = x$ab[i], guideline = guideline, info = FALSE)
    })
    log <- AMR::sir_interpretation_history(clean = TRUE)
    generic <- log$index[as.character(log$mo) %in% c("B_GRAMN", "B_GRAMP")]
    mos <- unique(x$mo[i])
    enterobacterales <- mos[AMR::mo_order(mos, info = FALSE) %in% "Enterobacterales"]
    sir[i[generic[x$mo[i][generic] %in% enterobacterales]]] <- NA
  }
  sir
}


#' Test each organism x drug pair and label it
#'
#' Two Fisher's exact tests per pair, each with a Benjamini-Hochberg correction across pairs:
#' the observed comparison (each year under its own guideline, what a naive antibiogram shows)
#' and the real one (both years under the newer guideline). A pair is a "significant change" if
#' the real test is significant, a "breakpoint artefact" if only the observed one is and the
#' artefact is larger than the real change, and "noise" otherwise.
#'
#' @param x A [drift_decompose()] table.
#' @param alpha Significance level after correction.
#' @return `x` with `p_observed`, `p_real`, Clopper-Pearson 95% intervals for both years under
#'   the newer guideline (`ci_low_*`, `ci_high_*`) and `verdict`.
#' @export
drift_test <- function(x, alpha = 0.05) {
  fisher <- function(r1, n1, r2, n2) mapply(function(r1, n1, r2, n2) {
    if (n1 == 0 || n2 == 0) return(NA_real_)
    stats::fisher.test(matrix(c(r1, n1 - r1, r2, n2 - r2), 2))$p.value
  }, r1, n1, r2, n2)
  x$p_observed <- stats::p.adjust(if (is.null(x$r_rep_y1)) {
    fisher(x$r_y1g1, x$n_y1g1, x$r_y2g2, x$n_y2g2)
  } else {
    fisher(x$r_rep_y1, x$n_rep_y1, x$r_rep_y2, x$n_rep_y2)
  }, "BH")
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


#' Shift in the MIC distribution between two periods, in doubling dilutions
#'
#' Each stored MIC is an interval on the log2 scale: "8" means the true MIC is above 4 and at
#' most 8, "<=0.25" at most 0.25, ">=16" (or ">8") above 8. For each organism and drug, a normal
#' model of log2 MIC with a period effect is fitted to these intervals by maximum likelihood
#' ([survival::survreg()]). The shift does not use any breakpoint, so it also catches MICs
#' creeping up below one, before %R moves. p-values are Benjamini-Hochberg corrected across
#' pairs.
#'
#' A period whose values are nearly all capped ("<=1" and nothing else) has no estimable mean.
#' The shift is NA unless each period has `min_n` values, of which `min_exact` are not capped,
#' and also when either period's fitted geometric mean lies outside the range the panel tested:
#' mostly capped values let the fit run below or above any tested dilution.
#'
#' @inheritParams drift_decompose
#' @param data As for [drift_decompose()]; only rows with an `mic` are used.
#' @param min_n,min_exact Minimum values, and uncapped values, per period.
#' @return A data frame with one row per organism and drug: counts, geometric mean MIC in each
#'   period (`gm_from`, `gm_to`), `shift` in doubling dilutions with a 95% interval, and
#'   `p_shift`.
#' @export
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
  tested <- range(c(x$lo, x$hi), na.rm = TRUE)
  means <- c(b[[1]], b[[1]] + b[[2]])
  if (any(means < tested[1] | means > tested[2])) return(res)
  res$gm_from <- 2^b[[1]]
  res$gm_to <- 2^(b[[1]] + b[[2]])
  res$shift <- b[[2]]
  res$shift_low <- b[[2]] - stats::qnorm(0.975) * se
  res$shift_high <- b[[2]] + stats::qnorm(0.975) * se
  res$p_shift <- 2 * stats::pnorm(-abs(b[[2]] / se))
  res
}


#' Write the stewardship report as Markdown, one section per organism
#'
#' @param x A [drift_test()] table. If it also has the columns of [drift_mic_shift()] (merged by
#'   `mo` and `ab`, MIC rows only), the report adds the MIC shift.
#' @param path File to write.
#' @param title Report title.
#' @param note Optional paragraph under the title.
#' @return `path`, invisibly.
#' @export
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
    if (is.null(x$p_rep_y1)) "- **Observed**: %R in each period under its own guideline, as an antibiogram reports it."
    else "- **Observed**: %R as originally reported in each period.",
    "- **Real**: change in the values, with the guideline held fixed (mean of holding either one).",
    "- **Artefact**: the rest of the observed change, caused by the breakpoint revision alone.",
    sprintf("- **Verdict**: Fisher's exact test, Benjamini-Hochberg corrected across %d pairs.", nrow(x)),
    sprintf("- **Read**: %s isolates that the old / new guideline could classify, of those tested. A capped MIC such as <=1 cannot be read against a breakpoint below it.", y1),
    if (has_shift) "- **MIC shift**: change in the MIC distribution in doubling dilutions, with 95% CI, from an interval-censored model that uses no breakpoint. * = significant after Benjamini-Hochberg correction.",
    "",
    paste0("Verdicts: ", paste(sprintf("%s %d", names(counts), counts), collapse = ", "), "."),
    ""
  )
  p1 <- if (is.null(x$p_rep_y1)) x$p_y1g1 else x$p_rep_y1
  p2 <- if (is.null(x$p_rep_y2)) x$p_y2g2 else x$p_rep_y2
  for (mo in unique(x$mo)) {
    i <- which(x$mo == mo)
    i <- i[order(is.na(x$observed[i]), -abs(x$observed[i]))]
    s <- x[i, ]
    shift <- if (!has_shift) NULL else
      ifelse(is.na(s$shift), "", sprintf("%+.1f (%+.1f to %+.1f)%s", s$shift, s$shift_low, s$shift_high,
                                         ifelse(!is.na(s$p_shift) & s$p_shift < 0.05, " *", "")))
    lines <- c(lines, paste0("## *", AMR::mo_name(mo, info = FALSE), "*"), "",
      paste0(sprintf("| Drug | Test | %%R %s | %%R %s | Observed | Real | Artefact | Verdict | Read %s |",
                     y1, y2, y1), if (has_shift) " MIC shift |"),
      paste0("|---|---|--:|--:|--:|--:|--:|---|--:|", if (has_shift) "--:|"),
      paste0(sprintf("| %s | %s %s | %s | %s | %s | %s | %s | %s | %s |",
                     AMR::ab_name(s$ab, info = FALSE), s$method, s$org, pct(p1[i]), pct(p2[i]),
                     pp(s$observed), pp(s$real), pp(s$artefact), s$verdict,
                     sprintf("%d / %d of %d", s$n_y1g1, s$n_y1g2, s$tested_y1)),
             if (has_shift) paste0(" ", shift, " |")),
      "")
  }
  writeLines(lines, path)
  invisible(path)
}
