# Split a change in %R between two years into a period component (the change with the guideline
# held fixed) and a guideline component (the part due to reading the two years with different
# interpretive criteria). The software calls them `real` and `artefact`.
#
# p(y, g) is the share of isolates from year y that are R under guideline g. With y1 < y2 and
# g1, g2 the guidelines in force in those years:
#   observed change     = p(y2, g2) - p(y1, g1)
#   period component    = p(y2, g2) - p(y1, g2)   (same guideline, values from two years)
#   guideline component = p(y1, g2) - p(y1, g1)   (same values, two guidelines)
# The split depends on the order of the two steps. Holding the old guideline instead gives
#   period component    = p(y2, g1) - p(y1, g1)
# The reported period component is the mean of the two (the Shapley value), and the guideline
# component is the rest, so the two parts always add up to the observed change. The two orderings
# differ by the interaction p(y2, g2) - p(y2, g1) - p(y1, g2) + p(y1, g1), which the Shapley value
# splits equally by convention.
#
# A capped MIC is an interval of possible true MICs. It is counted when every MIC in that interval
# falls on the same side of the outcome's boundary: "<=4" is not R under R >= 8, although AMR calls
# it NI because it may be S or I. Only a value that straddles the boundary is left out of p(y, g);
# the coverage columns say how many tested isolates each p(y, g) still rests on.


#' Guideline edition in force in a given year
#'
#' Assumes a lab uses the latest edition named for the year or before it, `lag` years late.
#' Years before the oldest edition in the AMR package give NA; years after the newest give the
#' newest.
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
  editions <- sort(as.integer(sub(".* ", "", editions[startsWith(editions, paste0(org, " "))])))
  if (length(editions) == 0) stop("The AMR package has no ", org, " breakpoints.")
  # The latest edition on or before the year, so a gap in the editions is never filled by a name
  # the AMR package does not have
  k <- findInterval(year - lag, editions)
  ifelse(is.na(k) | k == 0, NA_character_, paste(org, editions[pmax(k, 1)]))
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
#' Splits the observed change in %R for each organism and drug into the period component (the
#' change with the guideline held fixed, column `real`) and the guideline component (the rest,
#' column `artefact`). The period component is the mean of holding either guideline fixed (the
#' Shapley value), so the two parts add up to the observed change. A capped MIC that falls on
#' both sides of the outcome's boundary under a guideline is left out of that guideline's share;
#' the `n_*` columns say how many isolates each share rests on.
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
#'   guideline component also covers a lab's own reporting rules and any disagreement between the
#'   lab's reading and this one.
#' @param from,to A year each or a set of years such as `2013:2016`. A period is read with the
#'   edition in force at its last year, as a cumulative antibiogram for that period would have
#'   been.
#' @param in_force Function of (year, org) giving the edition in force.
#' @param adoption Optional table of adoption dates; see [adopted_edition()].
#' @param outcome `"R"` to count resistant isolates (%R), or `"NS"` to count non-susceptible ones
#'   (I, SDD or R). A breakpoint revision that only moves the I/R boundary changes %R but not %NS.
#' @return A data frame with one row per organism, drug, method (MIC or disk) and organisation:
#'   breakpoints for MICs and zones are revised separately, and a lab that switched
#'   organisation between the periods gets its own row (`"CLSI to EUCAST"`). Columns
#'   `p_y{1,2}g{1,2}` hold %R for each period under each guideline, then `observed`, `real`
#'   (the period component, with its two orderings `real_at_new` and `real_at_old`),
#'   `interaction` (`real_at_new - real_at_old`) and `artefact` (the guideline component), as
#'   shares. `up_y{1,2}` and `down_y{1,2}` count, descriptively, the isolates of each period,
#'   readable both ways, whose own two readings differ: counted under the new guideline but not
#'   the old one (`up`) or the reverse (`down`). With a `reported` column, an isolate's reading as
#'   reported is compared with the mean of its two guideline readings instead.
#'   `straddled_y{1,2}` count the capped values read under one guideline but straddling the
#'   outcome's boundary under the other (such as "<=8" against R >= 8 but not R >= 16); they are
#'   left out of the one share only, so the components are not identified exactly.
#'   `real_min`, `real_max`, `artefact_min` and `artefact_max` give the range each component can
#'   lie in, with every straddled value counted as counted or not; [drift_test()] makes a pair
#'   not interpretable when that range is `min_effect` or wider.
#' @examples
#' # CLSI 2023 lowered the Enterobacterales gentamicin breakpoints; the MICs did not change
#' mics <- rep(c("1", "2", "8", "16"), c(40, 20, 30, 10))
#' isolates <- data.frame(year = rep(c(2022, 2023), each = 100), mo = "B_ESCHR_COLI",
#'                        ab = "GEN", mic = mics)
#' d <- drift_decompose(isolates)
#' d[c("observed", "real", "artefact")]
#' @export
drift_decompose <- function(data, from = min(data$year), to = max(data$year),
                            in_force = guideline_in_force, adoption = NULL, outcome = "R") {
  rows <- interpret_rows(data, from, to, in_force, adoption, outcome)
  reclassified(rows, tabulate_rows(rows))
}

# Isolates whose own readings differ, signed so that `up` raises the guideline component: the
# change from the old-side to the new-side reading (period 1: as observed to the guidelines' mean;
# period 2: the guidelines' mean to as observed). Descriptive only; the guideline component is
# tested by the bootstrap.
reclassified <- function(rows, out) {
  d <- rows$data
  both <- d$read_g1 & d$read_g2
  obs <- if (is.null(d$read_rep)) ifelse(d$period == 1, d$hit_g1, d$hit_g2) else d$hit_rep
  if (!is.null(d$read_rep)) both <- both & d$read_rep
  s <- sign(obs - (d$hit_g1 + d$hit_g2) / 2) * ifelse(d$period == 1, -1, 1)
  for (y in 1:2) {
    out[[paste0("up_y", y)]] <- as.vector(table(rows$pair[both & d$period == y & s > 0]))
    out[[paste0("down_y", y)]] <- as.vector(table(rows$pair[both & d$period == y & s < 0]))
  }
  out
}

# Reads every row once under both editions; cell_sums() then sums them, with weights so
# drift_bootstrap() can resample without reading again.
interpret_rows <- function(data, from, to, in_force, adoption, outcome) {
  outcome <- match.arg(outcome, c("R", "NS"))
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
  out$year_from <- year_label(from)
  out$year_to <- year_label(to)
  editions <- function(g) as.vector(tapply(g, pair, function(x) paste(sort(unique(x)), collapse = " + ")))
  out$guideline_from <- editions(data$g1)
  out$guideline_to <- editions(data$g2)

  counted <- if (outcome == "R") "R" else c("I", "SDD", "R")
  data$period <- ifelse(data$year %in% from, 1L, 2L)
  ends <- mic_ends(data$mic)
  for (k in 1:2) {
    g <- data[[paste0("g", k)]]
    hit <- list()
    for (end in c("lo", "hi")) {
      x <- data
      x$mic <- ends[[end]]
      sir <- rep(NA_character_, nrow(data))
      for (edition in unique(g)) sir[g == edition] <- read_sir(x[g == edition, ], edition)
      hit[[end]] <- ifelse(is.na(sir) | sir == "NI", NA, sir %in% counted)
    }
    # Read when both ends of the value's interval are read and fall on the same side
    read <- !is.na(hit$lo) & !is.na(hit$hi) & hit$lo == hit$hi
    data[[paste0("read_g", k)]] <- read
    data[[paste0("hit_g", k)]] <- read & hit$hi
  }
  if (!is.null(data$reported)) {
    data$read_rep <- data$reported %in% c("S", "SDD", "I", "R")
    data$hit_rep <- data$reported %in% counted
  }
  list(data = data, out = out, pair = pair,
       group = as.integer(pair) + nlevels(pair) * (data$period - 1L))
}

# "2019-2020" for a run of years, "2015, 2018" otherwise
year_label <- function(w) {
  w <- sort(unique(w))
  if (length(w) == 1) as.character(w) else if (all(diff(w) == 1)) paste0(min(w), "-", max(w)) else paste(w, collapse = ", ")
}

# The two ends of each stored MIC's interval as plain MICs on the dilution grid: "<=4" spans the
# lowest MIC the AMR package knows to 4, ">8" spans 16 to the highest. A zone diameter (NA here)
# is read as it is.
mic_ends <- function(mic) {
  s <- as.character(AMR::as.mic(mic))
  number <- sub("^[<>=]+", "", s)
  sign <- sub("^([<>=]*).*", "\\1", s)
  step <- function(k) as.character(AMR::as.mic(2^(round(log2(as.numeric(number))) + k)))
  list(lo = ifelse(sign %in% c("<", "<="), "0.0001", ifelse(sign == ">", step(1), number)),
       hi = ifelse(sign %in% c(">", ">="), "4096", ifelse(sign == "<", step(-1), number)))
}

# Per pair and period, the sums every share is built from: isolates tested, read and counted under
# each guideline, and as reported. Each column of `w` is one set of row weights (one replicate);
# each sum is a matrix with a row per pair and period (period 1 first) and a column per replicate.
cell_sums <- function(rows, w = matrix(1, nrow(rows$data), 1)) {
  d <- rows$data
  v <- list(tested = rep(1, nrow(d)), n_g1 = d$read_g1, r_g1 = d$hit_g1, n_g2 = d$read_g2, r_g2 = d$hit_g2,
            u_g1 = d$read_g2 & !d$read_g1, u_g2 = d$read_g1 & !d$read_g2)
  if (!is.null(d$read_rep)) v <- c(v, list(n_rep = d$read_rep, r_rep = d$hit_rep))
  groups <- 2 * nlevels(rows$pair)
  lapply(v, function(x) {
    m <- matrix(0, groups, ncol(w))
    s <- rowsum(x * w, rows$group)
    m[as.integer(rownames(s)), ] <- s
    m
  })
}

# The decomposition from cell_sums(), one column per replicate
components <- function(s) {
  y1 <- seq_len(nrow(s$tested) / 2)
  y2 <- y1 + length(y1)
  share <- function(sum, i) s[[paste0("r_", sum)]][i, , drop = FALSE] / ifelse(s[[paste0("n_", sum)]][i, , drop = FALSE] > 0,
                                                                            s[[paste0("n_", sum)]][i, , drop = FALSE], NA)
  p11 <- share("g1", y1)
  p12 <- share("g2", y1)
  p21 <- share("g1", y2)
  p22 <- share("g2", y2)
  observed <- if (is.null(s$n_rep)) p22 - p11 else share("rep", y2) - share("rep", y1)
  at_new <- p22 - p12
  at_old <- p21 - p11
  real <- (at_new + at_old) / 2
  # Rounded so that a component that is 0 in exact arithmetic is 0, not 1e-17, in every replicate
  lapply(list(observed = observed, real = real, artefact = observed - real, real_at_new = at_new,
              real_at_old = at_old), round, 12)
}

# %R (or %NS) per pair, period and edition, and the decomposition
tabulate_rows <- function(rows) {
  out <- rows$out
  s <- cell_sums(rows)
  for (y in 1:2) {
    i <- seq_len(nrow(out)) + (y - 1) * nrow(out)
    out[[paste0("tested_y", y)]] <- s$tested[i, 1]
    for (k in c(paste0("g", 1:2), if (!is.null(s$n_rep)) "rep")) {
      cell <- if (k == "rep") paste0("rep_y", y) else paste0("y", y, k)
      n <- s[[paste0("n_", k)]][i, 1]
      r <- s[[paste0("r_", k)]][i, 1]
      out[[paste0("n_", cell)]] <- n
      out[[paste0("r_", cell)]] <- r
      out[[paste0("p_", cell)]] <- ifelse(n > 0, r / n, NA_real_)
    }
  }
  parts <- components(s)
  for (name in c("observed", "real_at_new", "real_at_old", "real", "artefact")) out[[name]] <- parts[[name]][, 1]
  out$interaction <- out$real_at_new - out$real_at_old
  range <- straddle_range(s)
  for (name in names(range)) {
    out[[name]] <- if (grepl("^straddled", name)) range[[name]] else ifelse(is.na(out$real), NA, range[[name]])
  }
  out
}

# A capped value read under one edition but straddling the outcome's boundary under the other is
# left out of that edition's share only, so the two shares of a period rest on different isolates.
# Counting every such value as hit, or as not hit, in each of the four shares gives the range the
# components can lie in; it collapses to the estimate when no value straddles one edition only.
straddle_range <- function(s) {
  y <- list(seq_len(nrow(s$tested) / 2))
  y[[2]] <- y[[1]] + length(y[[1]])
  corners <- as.matrix(expand.grid(rep(list(0:1), 4)))
  parts <- lapply(seq_len(nrow(corners)), function(j) {
    f <- s
    bit <- 0
    for (p in 1:2) for (k in c("g1", "g2")) {
      bit <- bit + 1
      u <- s[[paste0("u_", k)]][y[[p]], , drop = FALSE]
      f[[paste0("n_", k)]][y[[p]], ] <- f[[paste0("n_", k)]][y[[p]], ] + u
      f[[paste0("r_", k)]][y[[p]], ] <- f[[paste0("r_", k)]][y[[p]], ] + corners[j, bit] * u
    }
    components(f)
  })
  bound <- function(name, f) apply(do.call(cbind, lapply(parts, function(x) x[[name]][, 1])), 1, f)
  list(straddled_y1 = s$u_g1[y[[1]], 1] + s$u_g2[y[[1]], 1],
       straddled_y2 = s$u_g1[y[[2]], 1] + s$u_g2[y[[2]], 1],
       real_min = bound("real", min), real_max = bound("real", max),
       artefact_min = bound("artefact", min), artefact_max = bound("artefact", max))
}


#' Bootstrap intervals and tests for the observed change and its two components
#'
#' Resamples isolates (or clusters of isolates) with replacement within each period, reusing the
#' S/I/R readings of [drift_decompose()], and recomputes the decomposition each time. The
#' intervals and p-values are for the quantities reported: the Shapley period component and the
#' guideline component.
#'
#' With `cluster`, whole clusters are drawn, G - 1 of the G clusters of a period each time
#' (McCarthy and Snowden 1985; Rao and Wu 1988). Drawing G - 1 makes the bootstrap variance of a
#' linear statistic, such as a mean of cluster totals, unbiased; the shares here are ratios, so
#' this holds only approximately, and percentile intervals from few clusters remain approximate.
#' A period with fewer than `min_clusters` clusters is an error, and fewer than 10 a warning.
#' Clusters are resampled within each period, so a cluster sampled in both periods is treated as
#' two independent ones.
#'
#' @inheritParams drift_decompose
#' @param B Bootstrap replicates.
#' @param cluster Optional name of a column in `data` (such as a sampling site): whole clusters are
#'   resampled within each period instead of isolates. An isolate with no cluster is its own.
#' @param min_clusters Fewest clusters a period may have with `cluster`.
#' @param level Interval coverage.
#' @return The [drift_decompose()] table with percentile intervals `observed_low`,
#'   `observed_high`, `real_low`, `real_high`, `artefact_low` and `artefact_high`, and
#'   uncorrected two-sided bootstrap p-values (twice the smaller share of replicates at or beyond
#'   0, plus one, so a component that is 0 in every replicate gets 1): `boot_p_real` for the
#'   period component, `boot_p_artefact` for the guideline component, and `boot_p_larger` for
#'   |guideline| - |period|, the test of which component is larger (each component taken with the
#'   sign of its estimate, so the difference is linear and 0 when both estimates are 0). [drift_test()] corrects
#'   them. `boot_na` counts replicates in which a share had no isolate (left out of intervals and
#'   p-values), and `units_y1`, `units_y2` the isolates or clusters resampled in each period. Set
#'   a seed first for reproducible results.
#' @examples
#' mics <- rep(c("1", "2", "8", "16"), c(40, 20, 30, 10))
#' isolates <- data.frame(year = rep(c(2022, 2023), each = 100), mo = "B_ESCHR_COLI",
#'                        ab = "GEN", mic = mics, isolate = 1:200)
#' set.seed(1)
#' drift_bootstrap(isolates, B = 200)[c("artefact", "artefact_low", "artefact_high")]
#' @export
drift_bootstrap <- function(data, from = min(data$year), to = max(data$year),
                            in_force = guideline_in_force, adoption = NULL, outcome = "R",
                            B = 1000, cluster = NULL, min_clusters = 5, level = 0.95) {
  if (!is.null(cluster) && is.null(data[[cluster]])) stop("`data` has no column `", cluster, "`.")
  rows <- interpret_rows(data, from, to, in_force, adoption, outcome)
  d <- rows$data
  # The unit drawn: a cluster, else an isolate (all its drugs together), else a row
  own <- if (!is.null(d$isolate)) paste("isolate", d$isolate) else paste("row", seq_len(nrow(d)))
  unit <- if (!is.null(cluster)) ifelse(is.na(d[[cluster]]), own, paste("cluster", d[[cluster]])) else own
  unit <- paste(d$period, unit)
  first <- !duplicated(unit)
  units <- split(unit[first], d$period[first])
  index <- match(unit, unlist(units, use.names = FALSE))
  if (!is.null(cluster)) {
    g <- lengths(units)
    if (min(g) < min_clusters) stop("A period has ", min(g), " clusters; the cluster bootstrap needs at least ",
                                    min_clusters, ".")
    if (min(g) < 10) warning("A period has only ", min(g), " clusters; the cluster bootstrap is approximate.")
  }
  size <- vapply(units, function(u) if (is.null(cluster)) length(u) else length(u) - 1, numeric(1))
  est <- reclassified(rows, tabulate_rows(rows))

  draws <- lapply(c(observed = 1, real = 2, artefact = 3, larger = 4), function(i) matrix(NA_real_, nrow(est), B))
  chunk <- max(1, floor(2e7 / max(1, nrow(d))))  # replicates per pass, about 160 MB of weights
  for (start in seq(1, B, by = chunk)) {
    b <- min(chunk, B - start + 1)
    counts <- do.call(rbind, lapply(seq_along(units), function(p) {
      g <- length(units[[p]])
      draw <- matrix(sample.int(g, size[p] * b, replace = TRUE), size[p], b)
      matrix(apply(draw, 2, tabulate, g), g, b)
    }))
    parts <- components(cell_sums(rows, counts[index, , drop = FALSE]))
    cols <- start:(start + b - 1)
    draws$observed[, cols] <- parts$observed
    draws$real[, cols] <- parts$real
    draws$artefact[, cols] <- parts$artefact
    # |guideline| - |period| with the signs of the estimates, so the difference stays linear in
    # the draws and a component near 0 adds no upward bias of its absolute value
    draws$larger[, cols] <- sign(est$artefact) * parts$artefact - sign(est$real) * parts$real
  }

  q <- c((1 - level) / 2, 1 - (1 - level) / 2)
  for (name in c("observed", "real", "artefact")) {
    ci <- matrix(apply(draws[[name]], 1, stats::quantile, q, na.rm = TRUE, names = FALSE), 2)
    est[[paste0(name, "_low")]] <- ci[1, ]
    est[[paste0(name, "_high")]] <- ci[2, ]
  }
  boot_p <- function(x) {
    tail <- pmin(rowSums(x <= 0, na.rm = TRUE), rowSums(x >= 0, na.rm = TRUE))
    ifelse(is.na(est$real), NA, pmin(1, 2 * (tail + 1) / (rowSums(!is.na(x)) + 1)))
  }
  est$boot_p_real <- boot_p(draws$real)
  est$boot_p_artefact <- boot_p(draws$artefact)
  est$boot_p_larger <- boot_p(draws$larger)
  est$boot_na <- rowSums(is.na(draws$real))
  est$units_y1 <- length(units[["1"]])
  est$units_y2 <- length(units[["2"]])
  est
}

# S/I/R/NI as text for each row, from its MIC or zone diameter
#
# Where a guideline has no breakpoint for an organism's own group, AMR::as.sir() falls back to a
# broader one. For Enterobacterales that broader one is the generic Gram-negative row (CLSI
# Table 2B-5, "other non-Enterobacterales"), which is wrong for them: in AMR 3.1.0, CLSI 2011-2019
# lack Enterobacterales rows for many drugs, and E. coli ceftriaxone would be read with S <= 8
# instead of S <= 1. Such readings are dropped (NA), so they count as not readable; attribute
# `fallback` marks them. Each distinct value is read once.
read_sir <- function(x, guideline) {
  key <- paste(x$mo, x$ab, x$method, x$mic, x$disk)
  first <- !duplicated(key)
  if (!all(first)) {
    sir <- read_sir(x[first, ], guideline)
    m <- match(key, key[first])
    return(structure(as.vector(sir)[m], fallback = attr(sir, "fallback")[m]))
  }
  sir <- rep(NA_character_, nrow(x))
  fallback <- rep(FALSE, nrow(x))
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
    dropped <- i[generic[x$mo[i][generic] %in% enterobacterales]]
    sir[dropped] <- NA
    fallback[dropped] <- TRUE
  }
  structure(sir, fallback = fallback)
}


#' Test each organism x drug pair and label it
#'
#' Tests the observed change, the period component and the guideline component of each pair, and
#' which of the two components is larger, each with a Benjamini-Hochberg correction across the
#' interpretable pairs. The observed change gets Fisher's exact test (each period under its own
#' guideline or as reported, what a naive antibiogram shows). The two components, and the
#' difference |guideline| - |period|, get the bootstrap tests of [drift_bootstrap()], so each test
#' is of the quantity reported.
#'
#' A component is shown when its corrected p-value is below `alpha` and its size is at least
#' `min_effect`. Verdicts: "not interpretable" if the period component cannot be estimated, or if
#' capped values straddling one guideline's boundary leave a component unknown over a range of
#' `min_effect` or more (see [drift_decompose()]); "both" if both components are shown;
#' "guideline change" or "period change" if only that one is; "unresolved" if neither is but the
#' observed change is and both components are non-zero (it is real, but neither component can be
#' shown to drive it); otherwise "noise", which means no evidence of change, not evidence of none.
#' Column `larger` says which component is larger, "guideline" or "period", when the test of
#' |guideline| - |period| is significant and that component is shown, and "neither" otherwise.
#' Neither label depends on which component is tested first, and a large sample does not turn a
#' small component into the label. Column `verdict_note` gives the reason for "not interpretable".
#'
#' @param x A [drift_bootstrap()] table.
#' @param alpha Significance level after correction.
#' @param min_effect Smallest component, as a share (0.01 is one percentage point), that a
#'   verdict counts. 0 counts any significant component.
#' @return `x` with `p_observed`, `p_real` (period component), `p_artefact` (guideline
#'   component), `p_larger`, Clopper-Pearson 95% intervals for both years under the newer
#'   guideline (`ci_low_*`, `ci_high_*`), `larger`, `verdict` and `verdict_note`.
#' @export
drift_test <- function(x, alpha = 0.05, min_effect = 0.01) {
  if (is.null(x$boot_p_artefact)) {
    stop("drift_test() needs a drift_bootstrap() table: the components are tested by the bootstrap.")
  }
  fisher <- function(r1, n1, r2, n2) mapply(function(r1, n1, r2, n2) {
    if (n1 == 0 || n2 == 0) return(NA_real_)
    stats::fisher.test(matrix(c(r1, n1 - r1, r2, n2 - r2), 2))$p.value
  }, r1, n1, r2, n2)
  # Width of the range the straddled capped values leave each component in (0 for older tables)
  width <- if (is.null(x$artefact_max)) 0 else
    pmax(x$real_max - x$real_min, x$artefact_max - x$artefact_min)
  straddled <- !is.na(width) & width > 1e-9 & width >= min_effect - 1e-9
  ok <- !is.na(x$real) & !straddled
  x$verdict_note <- ifelse(is.na(x$real), "no reading under one guideline in a period",
            ifelse(straddled, sprintf("capped values leave a component unknown over %.1f points", 100 * width), ""))
  adjust <- function(p) stats::p.adjust(ifelse(ok, p, NA), "BH")  # over the interpretable pairs
  x$p_observed <- adjust(if (is.null(x$r_rep_y1)) {
    fisher(x$r_y1g1, x$n_y1g1, x$r_y2g2, x$n_y2g2)
  } else {
    fisher(x$r_rep_y1, x$n_rep_y1, x$r_rep_y2, x$n_rep_y2)
  })
  x$p_real <- adjust(x$boot_p_real)
  x$p_artefact <- adjust(x$boot_p_artefact)
  x$p_larger <- adjust(x$boot_p_larger)

  for (cell in c("y1g2", "y2g2")) {
    ci <- mapply(function(r, n) if (n == 0) c(NA, NA) else stats::binom.test(r, n)$conf.int,
                 x[[paste0("r_", cell)]], x[[paste0("n_", cell)]])
    x[[paste0("ci_low_", cell)]] <- ci[1, ]
    x[[paste0("ci_high_", cell)]] <- ci[2, ]
  }

  shown <- function(p, est) !is.na(p) & p < alpha & !is.na(est) & abs(est) >= min_effect - 1e-9
  period <- shown(x$p_real, x$real)
  guideline <- shown(x$p_artefact, x$artefact)
  # The larger component is named only when it is itself shown: the two tests are corrected in
  # separate families, so the test of |guideline| - |period| can pass where the component's own
  # test does not
  bigger <- ifelse(abs(x$artefact) > abs(x$real), "guideline", "period")
  x$larger <- ifelse(!ok, NA, ifelse(!is.na(x$p_larger) & x$p_larger < alpha &
                                     ifelse(bigger == "guideline", guideline, period), bigger, "neither"))
  # With one component exactly 0 the observed change is the other one, so it is judged by that
  # component's test alone, never "unresolved"
  both_present <- !is.na(x$real) & !is.na(x$artefact) & x$real != 0 & x$artefact != 0
  x$verdict <- ifelse(!ok, "not interpretable",
               ifelse(period & guideline, "both",
               ifelse(guideline, "guideline change",
               ifelse(period, "period change",
               ifelse(both_present & shown(x$p_observed, x$observed), "unresolved", "noise")))))
  x
}

verdict_levels <- c("period change", "guideline change", "both", "unresolved", "noise", "not interpretable")


#' Shift in the MIC distribution between two periods, in doubling dilutions
#'
#' Each stored MIC is an interval on the log2 scale: "8" means the true MIC is above 4 and at
#' most 8, "<=0.25" at most 0.25, ">=16" (or ">8") above 8. For each organism and drug, a normal
#' model of log2 MIC with a period effect and one standard deviation is fitted to these intervals
#' by maximum likelihood ([survival::survreg()]). The shift does not use any breakpoint, so it
#' also catches MICs creeping up below one, before %R moves. It summarises the bulk of the
#' distribution: where the two periods' spreads differ, or a resistant subpopulation moves on its
#' own, one location shift does not describe the change. p-values are Benjamini-Hochberg
#' corrected across pairs.
#'
#' A period whose values are nearly all capped ("<=1" and nothing else) has no estimable mean.
#' The shift is NA unless each period has `min_n` values, of which `min_exact` are not capped,
#' and also when either period's fitted geometric mean lies outside the range the panel tested:
#' mostly capped values let the fit run below or above any tested dilution. Column `note` gives
#' the reason for every NA.
#'
#' @inheritParams drift_decompose
#' @param data As for [drift_decompose()]; only rows with an `mic` are used.
#' @param min_n,min_exact Minimum values, and uncapped values, per period.
#' @return A data frame with one row per organism and drug: counts, geometric mean MIC in each
#'   period (`gm_from`, `gm_to`), `shift` in doubling dilutions with a 95% interval, `p_shift`,
#'   and `note`.
#' @export
# Limitation: one normal per period. A wild-type / non-wild-type mixture (PLOS ONE 2020,
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
                    shift = NA_real_, shift_low = NA_real_, shift_high = NA_real_, p_shift = NA_real_,
                    note = "")
  withheld <- function(why) {
    res$note <- why
    res
  }
  if (min(n) < min_n) return(withheld("too few values"))
  if (min(n_exact) < min_exact) return(withheld("too few uncapped values"))
  # A warning here means no finite estimate, e.g. every value in a period censored the same way
  fit <- tryCatch(survival::survreg(survival::Surv(lo, hi, type = "interval2") ~ late, data = x,
                                    dist = "gaussian"),
                  error = function(e) paste("fit failed:", conditionMessage(e)),
                  warning = function(w) paste("fit warning:", conditionMessage(w)))
  if (is.character(fit)) return(withheld(fit))
  b <- stats::coef(fit)
  se <- sqrt(stats::vcov(fit)["lateTRUE", "lateTRUE"])
  if (!is.finite(se)) return(withheld("no finite standard error"))
  tested <- range(c(x$lo, x$hi), na.rm = TRUE)
  means <- c(b[[1]], b[[1]] + b[[2]])
  if (any(means < tested[1] | means > tested[2])) return(withheld("fitted mean outside the tested range"))
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
#' @param outcome The `outcome` that `x` was computed with, `"R"` or `"NS"`; sets the labels.
#' @return `path`, invisibly.
#' @export
drift_report <- function(x, path, title = "Breakpoint drift report", note = NULL, outcome = "R") {
  outcome <- match.arg(outcome, c("R", "NS"))
  pct <- function(p) ifelse(is.na(p), "", sprintf("%.1f", 100 * p))
  pp <- function(p) ifelse(is.na(p), "", sprintf("%+.1f", 100 * p))
  y1 <- x$year_from[1]
  y2 <- x$year_to[1]
  editions <- unique(x[c("guideline_from", "guideline_to")])
  counts <- table(factor(x$verdict, verdict_levels))
  has_shift <- !is.null(x$shift)
  reported <- !is.null(x$p_rep_y1)

  lines <- c(
    paste("#", title), "",
    sprintf("%s isolates read with %s; %s isolates with %s. Changes are in percentage points of %%R.",
            y1, paste(editions$guideline_from, collapse = " or "),
            y2, paste(editions$guideline_to, collapse = " or ")),
    if (!is.null(note)) c("", note),
    "",
    if (!reported) "- **Observed**: %R in each period under its own guideline, as an antibiogram reports it."
    else "- **Observed**: %R as originally reported in each period.",
    "- **Period**: change in %R with the guideline held fixed (mean of holding either one). Not a causal estimate: sampling, panels and case mix can change too.",
    if (!reported) "- **Guideline**: the rest of the observed change, from reading the periods with different guidelines."
    else "- **Guideline**: the rest of the observed change: the guideline revision, the lab's own reporting rules, and any disagreement between the lab's reading and this one.",
    sprintf("- **Verdict**: which component is shown by its bootstrap test, Benjamini-Hochberg corrected across the %d interpretable pairs; **Larger**: which component is larger, where that is shown. Fisher's exact test of the observed change decides *unresolved*.", sum(!is.na(x$real))),
    sprintf("- **Read**: %s isolates that the old / new guideline could classify, of those tested. A capped MIC such as <=1 is left out only if it can fall on both sides of the breakpoint.", y1),
    if (has_shift) "- **MIC shift**: change in the MIC distribution in doubling dilutions, with 95% CI, from an interval-censored normal model that uses no breakpoint. * = significant after Benjamini-Hochberg correction.",
    "",
    paste0("Verdicts: ", paste(sprintf("%s %d", names(counts), counts), collapse = ", "), "."),
    ""
  )
  p1 <- if (!reported) x$p_y1g1 else x$p_rep_y1
  p2 <- if (!reported) x$p_y2g2 else x$p_rep_y2
  for (mo in unique(x$mo)) {
    i <- which(x$mo == mo)
    i <- i[order(is.na(x$observed[i]), -abs(x$observed[i]))]
    s <- x[i, ]
    shift <- if (!has_shift) NULL else
      ifelse(is.na(s$shift), "", sprintf("%+.1f (%+.1f to %+.1f)%s", s$shift, s$shift_low, s$shift_high,
                                         ifelse(!is.na(s$p_shift) & s$p_shift < 0.05, " *", "")))
    lines <- c(lines, paste0("## *", AMR::mo_name(mo, info = FALSE), "*"), "",
      paste0(sprintf("| Drug | Test | %%R %s | %%R %s | Observed | Period | Guideline | Verdict | Larger | Read %s |",
                     y1, y2, y1), if (has_shift) " MIC shift |"),
      paste0("|---|---|--:|--:|--:|--:|--:|---|---|--:|", if (has_shift) "--:|"),
      paste0(sprintf("| %s | %s %s | %s | %s | %s | %s | %s | %s | %s | %s |",
                     AMR::ab_name(s$ab, info = FALSE), s$method, s$org, pct(p1[i]), pct(p2[i]),
                     pp(s$observed), pp(s$real), pp(s$artefact), s$verdict, ifelse(is.na(s$larger), "", s$larger),
                     sprintf("%d / %d of %d", s$n_y1g1, s$n_y1g2, s$tested_y1)),
             if (has_shift) paste0(" ", shift, " |")),
      "")
  }
  if (outcome == "NS") lines <- gsub("%R", "%NS", lines, fixed = TRUE)
  writeLines(lines, path)
  invisible(path)
}
