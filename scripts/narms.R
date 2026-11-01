# Validate the drift audit on NARMS, the US National Antimicrobial Resistance Monitoring System.
# Run from the project root: Rscript scripts/narms.R
#
# NARMS samples retail meat and food-animal guts on a fixed plan and tests every isolate in
# one network of labs by broth microdilution, so unlike NCBI its isolates are not picked by
# submitters. FDA publishes the isolate-level MICs with NARMS's own S/I/R call, free and without
# registration ("not confidential"; cite NARMS):
#   https://www.fda.gov/animal-veterinary/national-antimicrobial-resistance-monitoring-system/integrated-reportssummaries
#   Retail Meats  cvm-narms-retail-meats_0.xlsx  SHA-256 a989afc9daa0a41b116425524b1e541fdccb09498d7897423f2fec3d3afeb336
#   Cecal         CVM-NARMS-Cecal.xlsx           SHA-256 116831649876a47623040f76d50f806d2ac9940fde3a9983187159452696d7af
# Put them in data/narms/ (not committed: 74 MB). The first run writes the Salmonella and E. coli
# columns used here to data/raw/narms_extract.csv.gz, which is committed, and later runs use that.
# It also writes results/narms_codes.csv, every drug code in the two files with its result count
# and why it is or is not analysed, since the extract keeps only the mapped codes.
#
# A record is a sample in which Salmonella or E. coli was sought; GROWTH says whether it was
# isolated. Only records with growth are isolates, and nearly all of them have MICs.
#
# Four results:
#   1. NARMS's own relabelling: up to 2022 NARMS called a Salmonella ciprofloxacin MIC of 0.12
#      or more R; from 2023 it uses CLSI (0.12-0.5 is I). Each drug's reported %R in 2022 and
#      2023 is split into real change (both years read with CLSI 2023) and artefact.
#   2. The tool's own decomposition, CLSI 2020 against CLSI 2023, on cecal isolates.
#   3. Agreement of the tool's reading with NARMS's calls, by drug and year.
#   4. The MIC shift on retail isolates over 20 years, which needs no breakpoint at all.
# Then further checks: wild-type cut-offs (ECOFF) instead of clinical
# breakpoints, the ciprofloxacin rise by host and serotype, and the isolate flow.

source("R/drift.R")

extract <- "data/raw/narms_extract.csv.gz"
# NARMS code -> AMR code. Spelled out because AMR::as.ab() reads some NARMS codes as other
# drugs (PTZ as pentizidone, TET as tetroxoprim).
drugs <- c(AMC = "AMC", AMI = "AMK", AMP = "AMP", ATM = "ATM", AXO = "CRO", AZI = "AZM",
           CAZ = "CAZ", CHL = "CHL", CIP = "CIP", COL = "COL", COT = "SXT", CTX = "CTX",
           FEP = "FEP", FOX = "FOX", GEN = "GEN", IMI = "IPM", MER = "MEM", NAL = "NAL",
           PTZ = "TZP", TET = "TCY", TGC = "TGC")

if (!file.exists(extract)) local({  # local: nothing from the build reaches the analysis below
  codes_seen <- list()
  read_source <- function(path, source) {
    d <- readxl::read_excel(path, col_types = "text")
    names(d)[names(d) %in% c("Year", "YEAR")] <- "year"
    names(d)[names(d) %in% c("Month", "MONTH")] <- "month"
    d <- d[d$GENUS_NAME %in% c("Salmonella", "Escherichia"), ]
    all_codes <- sub(" Sign$", "", grep(" Sign$", names(d), value = TRUE))
    all_codes <- all_codes[all_codes %in% names(d)]
    codes_seen[[source]] <<- data.frame(source = source, code = all_codes,
      results = vapply(all_codes, function(k) sum(!is.na(d[[k]]) & d$GROWTH %in% "YES"), numeric(1)),
      first_year = vapply(all_codes, function(k) suppressWarnings(min(as.integer(d$year[!is.na(d[[k]])]))), numeric(1)),
      last_year = vapply(all_codes, function(k) suppressWarnings(max(as.integer(d$year[!is.na(d[[k]])]))), numeric(1)))
    out <- data.frame(source = source, isolate = d$ISOLATE_ID, year = as.integer(d$year),
                      genus = d$GENUS_NAME, host = d$HOST_SPECIES, serotype = d$SEROTYPE,
                      month = d$month, state = d$STATE, growth = d$GROWTH)
    for (code in intersect(names(drugs), names(d))) {
      sign <- d[[paste(code, "Sign")]]
      value <- format(signif(as.numeric(d[[code]]), 3), scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
      out[[code]] <- ifelse(is.na(d[[code]]), NA, paste0(ifelse(sign %in% c("<", "<=", ">", ">="), sign, ""), value))
      out[[paste0(code, "_sir")]] <- d[[paste(code, "SIR")]]
    }
    out
  }
  raw <- merge(read_source("data/narms/cvm-narms-retail-meats_0.xlsx", "retail"),
               read_source("data/narms/CVM-NARMS-Cecal.xlsx", "cecal"), all = TRUE)
  # Rows in the order of the extract before month, state and growth were added (sorted on the
  # other columns)
  keep <- setdiff(names(raw), c("month", "state", "growth"))
  raw <- raw[c(keep, "month", "state", "growth")][do.call(order, unname(as.list(raw[keep]))), ]
  write.csv(raw, gzfile(extract), row.names = FALSE, na = "")

  # Every drug code, and why it is not analysed if it is not. Duplicate columns: the ones whose
  # values equal another code's in every Salmonella and E. coli record of both files.
  codes <- do.call(rbind, codes_seen)
  codes <- codes[codes$results > 0, ]
  ab <- suppressWarnings(as.character(AMR::as.ab(codes$code, info = FALSE)))
  ab[codes$code %in% names(drugs)] <- drugs[codes$code[codes$code %in% names(drugs)]]
  ab[codes$code %in% c("FIS", "SUF", "SMX")] <- "SSS"  # sulfisoxazole and other sulfonamides
  bp <- AMR::clinical_breakpoints
  has_bp <- ab %in% bp$ab[bp$guideline == "CLSI 2023" & bp$type == "human" & bp$method == "MIC" &
                          as.character(bp$mo) %in% c("B_[ORD]_ENTRBCTR", "B_SLMNL", "B_ESCHR_COLI")]
  codes$analysed <- codes$code %in% names(drugs)
  codes$clsi_2023_enterobacterales_breakpoint <- has_bp
  codes$reason <- ifelse(codes$analysed, "analysed",
                  ifelse(codes$code %in% c("CIP2", "SUF"), paste("duplicate of", ifelse(codes$code == "CIP2", "CIP", "FIS")),
                  ifelse(!has_bp, "no CLSI 2023 Enterobacterales MIC breakpoint in AMR 3.1.0",
                  ifelse(codes$last_year < 2019, paste("last tested in", codes$last_year), "highest concentration below the R breakpoint"))))
  write.csv(codes, "results/narms_codes.csv", row.names = FALSE)
})

raw <- read.csv(extract, colClasses = "character", na.strings = "")
raw$year <- as.integer(raw$year)
records <- raw  # every sample record, for the flow table
raw <- raw[raw$growth %in% "YES", ]  # isolates
raw$mo <- ifelse(raw$genus == "Salmonella", "B_SLMNL_ENTR", "B_ESCHR_COLI")
# 3,189 isolates have no isolate ID, 2,162 of them with MICs (cecal 2013 and 2023, retail 2015
# and 2020): each is its own isolate, so it is neither merged with the others as a duplicate nor resampled with them as
# one unit
no_id <- is.na(raw$isolate)
raw$isolate[no_id] <- paste0("record", which(no_id))

# One row per isolate and drug. NARMS panels print 2^-6 as 0.015 and 2^-4 as 0.06; the AMR
# package writes them 0.016 and 0.064.
amr_level <- c("0.015" = "0.016", "0.03" = "0.032", "0.06" = "0.064", "0.12" = "0.125")
codes <- intersect(names(drugs), names(raw))
long <- do.call(rbind, lapply(codes, function(code) {
  x <- raw[!is.na(raw[[code]]), ]
  if (nrow(x) == 0) return(NULL)
  number <- sub("^[<>=]+", "", x[[code]])
  # month and state: sampling batch (cecal, which records no state) and state (retail), for
  # the cluster sensitivity analyses at the end
  data.frame(source = x$source, isolate = x$isolate, year = x$year, mo = x$mo, ab = drugs[[code]],
             batch = ifelse(is.na(x$month), NA, paste(x$year, x$month)), state = x$state, serotype = x$serotype,
             mic = paste0(sub("[0-9.].*$", "", x[[code]]), ifelse(number %in% names(amr_level), amr_level[number], number)),
             narms = x[[paste0(code, "_sir")]])
}))
# Retail reuses an isolate ID for an E. coli and a Salmonella from the same sample (2018-2019),
# and in 2019 for two Salmonella of different serotypes, so organism and serotype are part of the
# key
tested <- long
long <- long[!duplicated(long[c("source", "mo", "serotype", "isolate", "ab")]), ]
dir.create("results", showWarnings = FALSE)
round_numbers <- function(x) {  # as in scripts/ncbi.R: the same file on every platform
  num <- vapply(x, is.numeric, logical(1))
  x[num] <- lapply(x[num], signif, 6)
  x
}

# ---- 3. Agreement with NARMS's calls ------------------------------------------------------------
# Read each isolate with the CLSI edition named for its year, from 2020 (the first edition with
# Enterobacterales rows for most drugs in the AMR package).
recent <- long[long$year >= 2020 & !is.na(long$narms), ]
recent$method <- "MIC"
recent$disk <- NA
recent$ours <- NA_character_
for (y in unique(recent$year)) {
  i <- recent$year == y
  recent$ours[i] <- read_sir(recent[i, ], guideline_in_force(y))
}
# calls: NARMS calls; tests: those the tool could also read (not NI, no generic fallback, a
# breakpoint exists), on which agreement is computed
recent$readable <- !is.na(recent$ours) & recent$ours != "NI"
recent$agree <- recent$readable & recent$ours == recent$narms
agreement <- aggregate(cbind(calls = 1, tests = readable, agree = agree) ~ source + mo + ab + year,
                       data = recent, FUN = sum)
agreement$agree_pct <- ifelse(agreement$tests > 0, 100 * agreement$agree / agreement$tests, NA)
write.csv(round_numbers(agreement), "results/narms_agreement.csv", row.names = FALSE)

# ---- 1. NARMS's own relabelling, 2022 to 2023 (cecal) --------------------------------------------
# Observed: %R as NARMS reported it each year. Period component: the MICs, read with CLSI 2022 and
# CLSI 2023 (Shapley mean of holding either). The guideline component is everything else: NARMS's
# own rules and the CLSI revisions. Pairs need 10 NARMS calls in each year (a working minimum). NARMS cecal calls through 2022 treat a Salmonella or E. coli ciprofloxacin MIC
# of 0.12 or more as R; in 2023 they follow CLSI, except that they keep the gentamicin breakpoint
# CLSI cut in 2023 (the agreement table shows both).
relabel_data <- long[long$source == "cecal" & long$year %in% 2022:2023, ]
relabel_data$reported <- relabel_data$narms
set.seed(2026)  # bootstrap intervals and tests (isolate level; the cecal file has no sampling site)
relabel <- drift_bootstrap(relabel_data, from = 2022, to = 2023, B = 5000)
relabel <- drift_test(relabel[relabel$n_rep_y1 >= 10 & relabel$n_rep_y2 >= 10, ])
write.csv(round_numbers(relabel), "results/narms_relabel.csv", row.names = FALSE)

# The same comparison for the share non-susceptible (I or R). NARMS's 2023 change moved
# ciprofloxacin 0.12-0.5 from R to I, so %R fell while %NS need not.
relabel_ns <- drift_bootstrap(relabel_data, from = 2022, to = 2023, outcome = "NS", B = 5000)
relabel_ns <- drift_test(relabel_ns[relabel_ns$n_rep_y1 >= 10 & relabel_ns$n_rep_y2 >= 10, ])
write.csv(round_numbers(relabel_ns), "results/narms_relabel_ns.csv", row.names = FALSE)

# ---- 2. CLSI 2020 against CLSI 2023 (cecal) ------------------------------------------------------
# Between these editions CLSI cut the Enterobacterales aminoglycoside breakpoints (2023) and
# piperacillin/tazobactam (2022).
cecal <- long[long$source == "cecal", ]
editions <- drift_bootstrap(cecal, from = 2019:2020, to = 2022:2023, B = 5000)
editions <- drift_test(editions[editions$n_y1g2 >= 10 & editions$n_y2g2 >= 10, ])
shift_cecal <- drift_mic_shift(cecal, from = 2019:2020, to = 2022:2023)
editions <- merge(editions, shift_cecal, by = c("mo", "ab"), all.x = TRUE)
write.csv(round_numbers(editions), "results/narms_editions.csv", row.names = FALSE)

# ---- 2b. Other windows -----------------------------------------------------------------------------
# The windows above were chosen after looking at the data. The same two comparisons over other
# windows, each period read with the edition in force in its last year, for the four pairs the
# paper discusses, and the number of pairs with each verdict.
window_runs <- list(
  list("relabel", 2022, 2023), list("relabel", 2021:2022, 2023), list("relabel", 2021, 2023),
  list("editions", 2019:2020, 2022:2023), list("editions", 2020, 2023), list("editions", 2019, 2023),
  list("editions", 2020:2021, 2022:2023), list("editions", 2019:2020, 2023))
set.seed(2026)
windows <- do.call(rbind, lapply(window_runs, function(w) {
  x <- if (w[[1]] == "relabel") transform(long[long$source == "cecal" & long$year %in% c(w[[2]], w[[3]]), ], reported = narms) else cecal
  b <- drift_bootstrap(x, from = w[[2]], to = w[[3]], B = 2000)
  n1 <- if (w[[1]] == "relabel") b$n_rep_y1 else b$n_y1g2
  n2 <- if (w[[1]] == "relabel") b$n_rep_y2 else b$n_y2g2
  t <- drift_test(b[n1 >= 10 & n2 >= 10, ])
  counts <- table(factor(t$verdict, verdict_levels))
  key <- t[t$ab %in% if (w[[1]] == "relabel") "CIP" else "GEN", ]
  cbind(comparison = w[[1]], key[c("year_from", "year_to", "guideline_from", "guideline_to", "mo", "ab", "observed",
                                   "real", "real_low", "real_high", "artefact", "artefact_low", "artefact_high", "verdict")],
        t(as.matrix(counts)), pairs = nrow(t))
}))
write.csv(round_numbers(windows), "results/narms_windows.csv", row.names = FALSE)

# ---- 2c. The generic Gram-negative fallback, if it were not dropped ------------------------------
# Cecal 2018-2019 read with CLSI 2019 against 2022-2023 read with CLSI 2023, with AMR's readings
# as they are (no fallback filter), and the same comparison with CLSI 2020, which has the
# Enterobacterales rows, as the older edition. The difference is the false guideline component
# the missing rows would produce.
read_sir_filtered <- read_sir
read_sir <- function(x, guideline) {
  key <- paste(x$mo, x$ab, x$mic)
  first <- !duplicated(key)
  sir <- as.character(AMR::as.sir(AMR::as.mic(x$mic[first]), mo = x$mo[first], ab = x$ab[first],
                                  guideline = guideline, capped_mic_handling = "standard", info = FALSE))
  sir[match(key, key[first])]
}
unfiltered <- drift_decompose(cecal, from = 2018:2019, to = 2022:2023)
read_sir <- read_sir_filtered
filtered <- drift_decompose(cecal, from = 2018:2019, to = 2022:2023)
with_2020 <- drift_decompose(cecal, from = 2018:2019, to = 2022:2023,
                             in_force = function(y, o) if (y == 2019) "CLSI 2020" else guideline_in_force(y, o))
fallback_check <- data.frame(unfiltered[c("mo", "ab")], n = unfiltered$n_y1g1,
                             artefact_clsi2019_unfiltered = unfiltered$artefact,
                             interpretable_clsi2019_filtered = !is.na(filtered$real),
                             artefact_clsi2020 = with_2020$artefact[match(paste(unfiltered$mo, unfiltered$ab), paste(with_2020$mo, with_2020$ab))])
fallback_check$false_part <- fallback_check$artefact_clsi2019_unfiltered - fallback_check$artefact_clsi2020
fallback_check <- fallback_check[fallback_check$n >= 10 & !is.na(fallback_check$artefact_clsi2019_unfiltered), ]
write.csv(round_numbers(fallback_check), "results/narms_fallback_check.csv", row.names = FALSE)

# ---- 4. Twenty years of MIC shift (retail) -------------------------------------------------------
retail <- long[long$source == "retail", ]
shift_retail <- drift_mic_shift(retail, from = 2002:2006, to = 2017:2021)
write.csv(round_numbers(shift_retail), "results/narms_mic_shift.csv", row.names = FALSE)

# ---- 5. Salmonella ciprofloxacin by year (cecal), for the figure -------------------------------
cip <- long[long$source == "cecal" & long$mo == "B_SLMNL_ENTR" & long$ab == "CIP", ]
cip$method <- "MIC"
cip$disk <- NA
sir <- read_sir(cip, "CLSI 2023")
cip_series <- do.call(rbind, lapply(sort(unique(cip$year)), function(y) {
  i <- cip$year == y
  read <- i & !is.na(sir) & sir != "NI"
  ci <- function(k, n) stats::binom.test(k, n)$conf.int
  rep_r <- sum(cip$narms[i] %in% "R")
  r <- sum(sir[read] == "R")
  ns <- sum(sir[read] %in% c("I", "R"))
  data.frame(year = y, n = sum(i), reported_r = rep_r / sum(i), clsi2023_r = r / sum(read), clsi2023_ns = ns / sum(read),
             reported_r_low = ci(rep_r, sum(i))[1], reported_r_high = ci(rep_r, sum(i))[2],
             clsi2023_r_low = ci(r, sum(read))[1], clsi2023_r_high = ci(r, sum(read))[2],
             clsi2023_ns_low = ci(ns, sum(read))[1], clsi2023_ns_high = ci(ns, sum(read))[2])
}))
write.csv(round_numbers(cip_series), "results/narms_cip_series.csv", row.names = FALSE)

# Colours: categorical slots 1-3 of the dataviz default palette, validated for colour-vision
# deficiency; each line also has its own marker and line type, so colour is never alone. Up to
# 2022 NARMS's R was CLSI's I or R, so the reported line is drawn dashed and offset by a fifth of
# a year to keep both visible. Vertical bars: Clopper-Pearson 95% intervals.
pdf("results/fig_cip_salmonella.pdf", width = 6.5, height = 3.6, pointsize = 9)
par(mar = c(4.2, 3.6, 0.8, 0.8), mgp = c(2.1, 0.6, 0), las = 1, bty = "n")
series <- list(list("clsi2023_r", "#1baf7a", 15, 1, "Read with CLSI 2023: R", -0.1),
               list("clsi2023_ns", "#2a78d6", 17, 1, "Read with CLSI 2023: I or R", -0.1),
               list("reported_r", "#eb6834", 16, 2, "Reported by NARMS: R", 0.1))
plot(NA, xlim = range(cip_series$year) + c(-0.3, 0.3), ylim = c(0, 20), xlab = "", ylab = "Isolates (%)",
     xaxt = "n", yaxt = "n")
title(xlab = "Year", line = 2.9)
axis(1, at = cip_series$year, col = "grey60", col.axis = "grey25", lwd = 0.5)
axis(1, at = cip_series$year, labels = paste0("n=", cip_series$n), line = 0.9, tick = FALSE,
     col.axis = "grey40", cex.axis = 0.8)
axis(2, at = seq(0, 20, 5), col = "grey60", col.axis = "grey25", lwd = 0.5)
abline(h = seq(0, 20, 5), col = "grey90", lwd = 0.5)
abline(v = 2022.5, col = "grey40", lty = 3)  # NARMS moves to CLSI's ciprofloxacin rules
text(2022.5, 19.5, "NARMS adopts\nCLSI rules (2023)", pos = 2, cex = 0.85, col = "grey25")
for (s in series) {
  x <- cip_series$year + s[[6]]
  segments(x, 100 * cip_series[[paste0(s[[1]], "_low")]], x, 100 * cip_series[[paste0(s[[1]], "_high")]], col = s[[2]], lwd = 1)
  lines(x, 100 * cip_series[[s[[1]]]], col = s[[2]], lwd = 2, lty = s[[4]])
  points(x, 100 * cip_series[[s[[1]]]], col = s[[2]], pch = s[[3]], cex = 1.1)
}
legend("topleft", legend = rev(sapply(series, `[[`, 5)), col = rev(sapply(series, `[[`, 2)),
       pch = rev(sapply(series, `[[`, 3)), lty = rev(sapply(series, `[[`, 4)), lwd = 2,
       bty = "n", text.col = "grey15", seg.len = 2.5)
invisible(dev.off())
writeLines(sprintf(paste(
  "\\textit{Salmonella} ciprofloxacin, cecal isolates %d--%d (%s--%s isolates a year; n under each",
  "year): share R as reported by NARMS, and share R and I or R with the same MICs read with CLSI",
  "2023, with Clopper-Pearson 95\\%% intervals; the reported series is drawn slightly to the right. Dotted",
  "line: from 2023 NARMS reports ciprofloxacin by CLSI's rules (0.12--0.5~mg/L is I); through 2022 it",
  "called 0.12~mg/L or more R, so its series matched the CLSI share of I or R."), min(cip_series$year), max(cip_series$year),
  format(min(cip_series$n), big.mark = ","), format(max(cip_series$n), big.mark = ",")),
  "results/fig_cip_salmonella_caption.tex")

# ---- Supplementary tables (LaTeX) ----------------------------------------------------------------
pct <- function(p) ifelse(is.na(p), "", sprintf("%.1f", 100 * p))
# sprintf rounds -0.0004 to "-0.0"; a part that is 0 in every replicate (no isolate between the
# two breakpoints) is printed as 0
pp <- function(p) ifelse(is.na(p), "", sub("^[-+](0\\.0)$", "\\1", sprintf("%+.1f", 100 * p)))
ci <- function(est, lo, hi) ifelse(!is.na(lo) & est == 0 & lo == 0 & hi == 0, "0",
                                   ifelse(is.na(lo), pp(est), sprintf("%s (%s to %s)", pp(est), pp(lo), pp(hi))))
pval <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "$<$0.001", sprintf("%.3f", p)))
org <- function(mo) ifelse(mo == "B_SLMNL_ENTR", "\\textit{Salmonella}", "\\textit{E.~coli}")
# The paper's terms for the software's labels: the guideline component is the criteria component,
# the period component the fixed-criteria component
term <- function(v) unname(c(`guideline change` = "criteria change", `period change` = "fixed-criteria change",
                             both = "both", unresolved = "change not attributed", noise = "no change detected",
                             `not interpretable` = "not interpretable", guideline = "criteria", period = "fixed criteria",
                             neither = "neither")[v])
# A long table breaks across pages (longtable, caption on top, header repeated)
table_tex <- function(spec, header, rows, caption, long = FALSE) {
  if (long) return(c("\\clearpage", "{\\scriptsize", paste0("\\begin{longtable}{", spec, "}"), paste0("\\caption{", caption, "}\\\\"),
                     "\\toprule", header, "\\midrule", "\\endfirsthead", "\\toprule", header, "\\midrule", "\\endhead",
                     rows, "\\bottomrule", "\\end{longtable}}", ""))
  c("\\begin{table}[h]", "\\centering", "\\scriptsize", paste0("\\begin{tabular}{", spec, "}"), "\\toprule",
    header, "\\midrule", rows, "\\bottomrule", "\\end{tabular}", paste0("\\caption{", caption, "}"), "\\end{table}", "")
}
supplement <- function(x, caption) {
  rep <- !is.null(x$p_rep_y1)
  p1 <- if (rep) x$p_rep_y1 else x$p_y1g1
  p2 <- if (rep) x$p_rep_y2 else x$p_y2g2
  n1 <- if (rep) x$n_rep_y1 else x$n_y1g1
  n2 <- if (rep) x$n_rep_y2 else x$n_y2g2
  rows <- sprintf("%s & %s & %d & %d & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", org(x$mo),
                  AMR::ab_name(x$ab, info = FALSE), n1, n2, pct(p1), pct(p2), pp(x$observed),
                  ci(x$real, x$real_low, x$real_high), ci(x$artefact, x$artefact_low, x$artefact_high),
                  pval(x$p_observed), pval(x$p_real), pval(x$p_artefact), pval(x$p_larger), term(x$verdict),
                  ifelse(is.na(x$larger), "", term(x$larger)))
  table_tex("llrrrrrllrrrrll",
            "Organism & Drug & n before & n after & Before & After & Observed & Fixed criteria (95\\% CI) & Criteria (95\\% CI) & $p$ obs. & $p$ fixed & $p$ criteria & $p$ larger & Verdict & Larger \\\\",
            rows, caption)
}
columns <- paste(
  "n: isolates read (or reported) in each period. Before, After: shares in percent; Observed and the components in",
  "percentage points, with percentile intervals from 5,000 bootstrap replicates over isolates (not adjusted); 0: zero",
  "in every replicate. $p$ obs.: Fisher's exact test of the observed change; $p$ fixed, $p$ criteria: bootstrap tests",
  "of the components; $p$ larger: bootstrap test of |criteria| $-$ |fixed criteria|; each column Benjamini-Hochberg",
  "corrected over the interpretable pairs of this table. A component is detected when its corrected $p$ is below 0.05",
  "and it is at least one point; Verdict says which components are detected, Larger which is larger where that test",
  "is significant and that component is detected. Pairs without a breakpoint are not interpretable.")
writeLines(c(
  supplement(relabel, paste("Cecal isolates, 2022 against 2023, \\%R as reported by NARMS. Fixed-criteria component from MICs read",
                            "with CLSI 2022 and CLSI 2023; criteria component: the rest of the observed change.", columns)),
  supplement(relabel_ns, paste("Cecal isolates, 2022 against 2023, share non-susceptible (I or R) as reported by NARMS.",
                               "Fixed-criteria component from MICs read with CLSI 2022 and CLSI 2023; criteria component: the",
                               "rest of the observed change.", columns)),
  supplement(editions, paste("Cecal isolates, 2019-2020 read with CLSI 2020 against 2022-2023 read with CLSI 2023, \\%R.",
                             "Fixed-criteria component: the change with one edition applied to both periods; criteria",
                             "component: the rest.", columns))),
  "results/supplement_tables.tex")

# ---- 6. Wild-type cut-offs (ECOFF) ---------------------------------------------------------------
# %non-wild-type uses no clinical breakpoint at all, so no revision can move it. EUCAST ECOFFs
# from the AMR package (EUCAST 2026 table; earlier tables have none for these pairs):
# ciprofloxacin 0.064 mg/L and gentamicin 2 mg/L for both organisms. AMR keeps the Salmonella
# ECOFFs under S. enterica subsp. enterica, not under S. enterica.
ecoff_mo <- function(mo) ifelse(mo == "B_SLMNL_ENTR", "B_SLMNL_ENTR_ENTR", mo)
ecoff <- do.call(rbind, lapply(list(list("B_SLMNL_ENTR", "CIP", 2022, 2023), list("B_ESCHR_COLI", "CIP", 2022, 2023),
                                    list("B_SLMNL_ENTR", "GEN", 2019:2020, 2022:2023),
                                    list("B_ESCHR_COLI", "GEN", 2019:2020, 2022:2023)), function(e) {
  x <- long[long$source == "cecal" & long$mo == e[[1]] & long$ab == e[[2]], ]
  sir <- as.character(AMR::as.sir(AMR::as.mic(x$mic), mo = ecoff_mo(x$mo), ab = x$ab, guideline = "EUCAST 2026",
                                  breakpoint_type = "ECOFF", capped_mic_handling = "standard", info = FALSE))
  bp <- AMR::clinical_breakpoints
  cut <- bp$breakpoint_S[bp$type == "ECOFF" & bp$guideline == "EUCAST 2026" & bp$method == "MIC" &
                         bp$ab == e[[2]] & as.character(bp$mo) == ecoff_mo(e[[1]])]
  count <- function(years) {
    i <- x$year %in% years
    c(sum(i), sum(i & sir %in% c("WT", "NWT")), sum(i & sir %in% "NWT"))
  }
  a <- count(e[[3]])
  b <- count(e[[4]])
  label <- function(w) if (length(w) == 1) as.character(w) else paste0(min(w), "-", max(w))
  data.frame(mo = e[[1]], ab = e[[2]], ecoff = cut, year_from = label(e[[3]]), year_to = label(e[[4]]),
             tested_from = a[1], read_from = a[2], nwt_from = a[3], pct_nwt_from = 100 * a[3] / a[2],
             tested_to = b[1], read_to = b[2], nwt_to = b[3], pct_nwt_to = 100 * b[3] / b[2],
             change = 100 * (b[3] / b[2] - a[3] / a[2]),
             p_fisher = stats::fisher.test(matrix(c(a[3], a[2] - a[3], b[3], b[2] - b[3]), 2))$p.value)
}))
write.csv(round_numbers(ecoff), "results/narms_ecoff.csv", row.names = FALSE)

# ---- 6b. Gentamicin MIC distributions (cecal), for the figure ----------------------------------
# The guideline component of Table S3 made visible: the same MICs, the CLSI 2020 R breakpoint
# (>= 16), the CLSI 2023 one (>= 8) and the ECOFF (2 mg/L). Share of isolates at each stored value.
gen <- long[long$source == "cecal" & long$ab == "GEN" & long$year %in% c(2019:2020, 2022:2023), ]
gen$period <- ifelse(gen$year <= 2020, "2019-2020", "2022-2023")
levels_gen <- as.character(sort(unique(AMR::as.mic(gen$mic))))
gen_dist <- as.data.frame(with(gen, table(mo, period, mic = factor(as.character(AMR::as.mic(mic)), levels_gen))))
gen_dist$share <- gen_dist$Freq / ave(gen_dist$Freq, gen_dist$mo, gen_dist$period, FUN = sum)
write.csv(round_numbers(gen_dist), "results/narms_gen_distribution.csv", row.names = FALSE)
# Log scale, so the few isolates near the breakpoints are visible next to the wild-type peak
pdf("results/fig_gen_mic.pdf", width = 6.5, height = 3.6, pointsize = 9)
par(mfrow = c(1, 2), mar = c(4, 3.6, 2, 0.5), oma = c(2.6, 0, 0, 0), mgp = c(2.1, 0.6, 0), las = 1, bty = "n")
for (mo in c("B_SLMNL_ENTR", "B_ESCHR_COLI")) {
  d <- gen_dist[gen_dist$mo == mo, ]
  m <- 100 * rbind(d$share[d$period == "2019-2020"], d$share[d$period == "2022-2023"])
  x <- barplot(m, beside = TRUE, names.arg = levels_gen, col = c("#9fb7d9", "#2a78d6"), border = NA,
               log = "y", ylim = c(0.1, 100), ylab = "Isolates (%, log scale)", xlab = "Gentamicin MIC (mg/L)", cex.names = 0.8,
               main = if (mo == "B_SLMNL_ENTR") "Salmonella" else "E. coli", font.main = 3)
  # Boundaries between bars: just below the first value counted resistant, and just above the ECOFF
  edge <- function(v) mean(x[, match(v, levels_gen) - c(1, 0)])
  abline(v = edge("16"), lty = 2, col = "grey30")
  abline(v = edge("8"), lty = 1, col = "grey30")
  abline(v = edge("4"), lty = 3, col = "#eb6834")
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot.new()
legend("bottom", c("2019-2020", "2022-2023", "R from 16, CLSI 2020", "R from 8, CLSI 2023", "above ECOFF 2"),
       fill = c("#9fb7d9", "#2a78d6", NA, NA, NA), border = NA, lty = c(NA, NA, 2, 1, 3),
       col = c(NA, NA, "grey30", "grey30", "#eb6834"), bty = "n", ncol = 3, cex = 0.85, seg.len = 1.8)
invisible(dev.off())

# ---- 7. Salmonella ciprofloxacin %NS by host and serotype (cecal, 2022 to 2023) ------------------
# Read with CLSI 2023 (its ciprofloxacin breakpoint for Salmonella is the 2022 one). Direct
# standardisation: each year's %NS per stratum, weighted by the stratum's share of both years
# pooled, over strata seen in both years. Serotypes with under 20 isolates pooled are "Other".
ns <- long[long$source == "cecal" & long$mo == "B_SLMNL_ENTR" & long$ab == "CIP" & long$year %in% 2022:2023, ]
ns$method <- "MIC"
ns$disk <- NA
sir <- read_sir(ns, "CLSI 2023")
ns <- ns[!is.na(sir) & sir != "NI", ]
ns$hit <- sir[!is.na(sir) & sir != "NI"] %in% c("I", "R")
m <- match(paste(ns$source, ns$mo, ns$isolate), paste(raw$source, raw$mo, raw$isolate))
ns$host <- ifelse(is.na(raw$host[m]), "Unknown", raw$host[m])
ns$serotype <- ifelse(is.na(raw$serotype[m]), "Unknown", sub("^S\\. ", "", trimws(raw$serotype[m])))
ns$serotype[ns$serotype %in% names(which(table(ns$serotype) < 20))] <- "Other"
strata <- do.call(rbind, lapply(c("host", "serotype"), function(by) {
  n <- table(factor(ns[[by]]), ns$year)
  r <- table(factor(ns[[by]])[ns$hit], ns$year[ns$hit])
  s <- data.frame(by = by, stratum = rownames(n), n_2022 = n[, "2022"], ns_2022 = r[, "2022"],
                  n_2023 = n[, "2023"], ns_2023 = r[, "2023"])
  s$pct_2022 <- 100 * s$ns_2022 / s$n_2022
  s$pct_2023 <- 100 * s$ns_2023 / s$n_2023
  both <- s$n_2022 > 0 & s$n_2023 > 0
  s$weight <- ifelse(both, (s$n_2022 + s$n_2023) / sum((s$n_2022 + s$n_2023)[both]), 0)
  s$change <- s$pct_2023 - s$pct_2022
  total <- function(name, p1, p2) data.frame(by = by, stratum = name, n_2022 = sum(s$n_2022), ns_2022 = sum(s$ns_2022),
                                             n_2023 = sum(s$n_2023), ns_2023 = sum(s$ns_2023), pct_2022 = p1,
                                             pct_2023 = p2, weight = 1, change = p2 - p1)
  rbind(s[order(-s$n_2022 - s$n_2023), ],
        total("All (crude)", 100 * sum(s$ns_2022) / sum(s$n_2022), 100 * sum(s$ns_2023) / sum(s$n_2023)),
        total("All (standardised)", sum(s$weight * s$pct_2022, na.rm = TRUE), sum(s$weight * s$pct_2023, na.rm = TRUE)))
}))
write.csv(round_numbers(strata), "results/narms_strata.csv", row.names = FALSE)

# Standardised change with bootstrap intervals, by host, by serotype and by host and serotype
# jointly (Infantis is concentrated in chickens, so the two are confounded). Isolates resampled
# within year, 2,000 replicates; the strata present in both years are fixed by the data.
ns$joint <- paste(ns$host, ns$serotype)
standardise <- function(x, by) {
  keep <- x[[by]] %in% intersect(x[[by]][x$year == 2022], x[[by]][x$year == 2023])
  x <- x[keep, ]
  w <- table(x[[by]]) / nrow(x)
  rate <- function(y) tapply(x$hit[x$year == y], factor(x[[by]][x$year == y], names(w)), mean)
  c(100 * sum(w * rate(2022), na.rm = TRUE), 100 * sum(w * rate(2023), na.rm = TRUE), sum(keep) / length(keep))
}
set.seed(2026)
standardised <- do.call(rbind, lapply(c("host", "serotype", "joint"), function(by) {
  est <- standardise(ns, by)
  draws <- replicate(2000, {
    b <- do.call(rbind, lapply(split(ns, ns$year), function(d) d[sample.int(nrow(d), replace = TRUE), ]))
    e <- standardise(b, by)
    e[2] - e[1]
  })
  ci <- stats::quantile(draws, c(0.025, 0.975), names = FALSE)
  data.frame(by = by, pct_2022 = est[1], pct_2023 = est[2], change = est[2] - est[1], change_low = ci[1],
             change_high = ci[2], share_in_common_strata = est[3])
}))
crude <- c(100 * mean(ns$hit[ns$year == 2022]), 100 * mean(ns$hit[ns$year == 2023]))
standardised <- rbind(data.frame(by = "crude", pct_2022 = crude[1], pct_2023 = crude[2], change = diff(crude),
                                 change_low = NA, change_high = NA, share_in_common_strata = 1), standardised)
write.csv(round_numbers(standardised), "results/narms_standardised.csv", row.names = FALSE)

# ---- 8. Isolate flow and drugs -------------------------------------------------------------------
# Per file and year: sample records in the extract (Salmonella and E. coli sought), those with no
# growth, isolates (growth), isolates with at least one MIC, MIC results
# before and after de-duplication (old key source + isolate + drug, which merged records without
# an ID and the E. coli and Salmonella sharing one; new key with the organism), then how each
# result reads with the CLSI edition named for its year: interpretable, capped and not
# interpretable (NI), dropped as a generic Gram-negative fallback, no breakpoint, or a year
# before the AMR package's first CLSI edition (2011).
flow <- long
flow$method <- "MIC"
flow$disk <- NA
flow$guideline <- guideline_in_force(flow$year)
sir <- rep(NA_character_, nrow(flow))
fallback <- rep(FALSE, nrow(flow))
for (g in unique(na.omit(flow$guideline))) {
  i <- which(flow$guideline %in% g)
  s <- read_sir(flow[i, ], g)
  sir[i] <- s
  fallback[i] <- attr(s, "fallback")
}
old_id <- ifelse(startsWith(tested$isolate, "record"), NA, tested$isolate)
by_year <- function(x, v) tapply(v, factor(paste(x$source, x$year), keys), sum, default = 0)
keys <- sort(unique(paste(records$source, records$year)))
has_mic <- paste(raw$source, raw$mo, raw$isolate) %in% paste(long$source, long$mo, long$isolate)
counts <- data.frame(
  records = by_year(records, rep(1, nrow(records))),
  no_growth = by_year(records, !records$growth %in% "YES"),
  isolates = by_year(raw, rep(1, nrow(raw))),
  isolates_mic = by_year(raw, has_mic),
  isolates_no_id = by_year(raw, no_id),
  results = by_year(tested, rep(1, nrow(tested))),
  results_old_key = by_year(tested, !duplicated(data.frame(tested$source, old_id, tested$ab))),
  results_dedup = by_year(flow, rep(1, nrow(flow))),
  interpretable = by_year(flow, !is.na(sir) & sir != "NI"),
  capped_ni = by_year(flow, sir %in% "NI"),
  fallback_only = by_year(flow, fallback),
  no_breakpoint = by_year(flow, !is.na(flow$guideline) & is.na(sir) & !fallback),
  no_edition = by_year(flow, is.na(flow$guideline)))
flow_table <- cbind(source = sub(" .*", "", keys), year = as.integer(sub(".* ", "", keys)), counts)
flow_table <- rbind(flow_table, cbind(source = "total", year = NA, as.data.frame(t(colSums(counts)))))
rownames(flow_table) <- NULL
write.csv(flow_table, "results/narms_flow.csv", row.names = FALSE, na = "")
# Table S4, appended to the supplementary tables written above. No result was capped and not
# interpretable with the edition of its own year, so that column is left out, and de-duplication
# removed none, so that one too.
stopifnot(all(flow_table$capped_ni == 0), all(flow_table$results == flow_table$results_dedup))
big <- function(v) format(v, big.mark = ",", trim = TRUE)
cat(c("\\begin{table}[h]", "\\centering", "\\scriptsize", "\\begin{tabular}{llrrrrrrrrr}", "\\toprule",
      "File & Year & Records & No growth & Isolates & With MIC & Results & Interpretable & Generic fallback & No breakpoint & Before 2011 \\\\",
      "\\midrule",
      with(flow_table, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", ifelse(source == "total", "Total", source),
                               ifelse(is.na(year), "", year), big(records), big(no_growth), big(isolates), big(isolates_mic),
                               big(results),
                               big(interpretable), big(fallback_only), big(no_breakpoint), big(no_edition))),
      "\\bottomrule", "\\end{tabular}",
      paste("\\caption{Sample records, isolates and MIC results by file and year, \\textit{Salmonella} and \\textit{E.~coli}.",
            "A record is a sample cultured for the organism; no growth: none isolated (GROWTH = NO). With MIC: isolates",
            "with at least one MIC result for the 21 drug codes analysed. Each result is one file, organism, serotype,",
            "isolate and drug; de-duplication on that key removed none. The last four columns split",
            "these results by how they read with the CLSI edition named for their year; none was capped and",
            "not interpretable.}"),
      "\\end{table}", ""), file = "results/supplement_tables.tex", sep = "\n", append = TRUE)

# Which drugs reach which table. The map above lists 21 NARMS codes; tigecycline has no
# Salmonella or E. coli result in either file, so the extract holds 20.
organisms <- function(x, ab) paste(sort(unique(ifelse(x$mo[x$ab == ab] == "B_SLMNL_ENTR", "Salmonella", "E. coli"))), collapse = "; ")
drug_table <- data.frame(code = names(drugs), ab = unname(drugs),
                         results_retail = vapply(drugs, function(a) sum(long$ab == a & long$source == "retail"), numeric(1)),
                         results_cecal = vapply(drugs, function(a) sum(long$ab == a & long$source == "cecal"), numeric(1)),
                         results_cecal_2019_2023 = vapply(drugs, function(a) sum(long$ab == a & long$source == "cecal" & long$year >= 2019), numeric(1)),
                         reported_cecal_2022_2023 = vapply(drugs, function(a) sum(relabel_data$ab == a & !is.na(relabel_data$reported)), numeric(1)))
for (t in c("relabel", "relabel_ns", "editions", "shift_retail", "agreement")) {
  x <- get(t)
  if (t == "shift_retail") x <- x[!is.na(x$shift), ]
  drug_table[[t]] <- vapply(drug_table$ab, function(a) organisms(x, a), character(1))
}
write.csv(drug_table, "results/narms_drugs.csv", row.names = FALSE)

# ---- Report --------------------------------------------------------------------------------------
note <- function(...) paste(...)
drift_report(relabel, "results/narms_relabel.md",
             title = "NARMS cecal isolates: NARMS's 2023 relabelling",
             note = note(
               "Up to 2022 NARMS reported cecal isolates with its own rules, e.g. a Salmonella",
               "ciprofloxacin MIC of 0.12 or more as R; from 2023 it follows CLSI, but keeps the",
               "gentamicin breakpoint CLSI cut in 2023. **Observed** is %R as NARMS reported it;",
               "**period** reads the MICs with CLSI 2022 and CLSI 2023."))
drift_report(editions, "results/narms_editions.md",
             title = "NARMS cecal isolates: CLSI 2020 against CLSI 2023",
             note = note(sprintf("%d Salmonella and E. coli isolates from food-animal guts at slaughter.",
                                 length(unique(cecal$isolate[cecal$year %in% c(2019:2020, 2022:2023)])))))

agree_total <- aggregate(cbind(tests, agree) ~ mo + ab, data = agreement[agreement$source == "cecal" & agreement$year == 2023, ], FUN = sum)
message("Sample records: ", nrow(records), "; isolates (growth): ", nrow(raw), "; with an MIC: ", sum(has_mic),
        "; isolate x drug tests: ", nrow(long))
message("Agreement with NARMS calls, cecal 2023: ", sum(agree_total$agree), " / ", sum(agree_total$tests))
message("Relabelling verdicts: ", paste(names(table(relabel$verdict)), table(relabel$verdict), collapse = ", "))
message("Non-susceptible relabelling verdicts: ", paste(names(table(relabel_ns$verdict)), table(relabel_ns$verdict), collapse = ", "))
message("Edition verdicts: ", paste(names(table(editions$verdict)), table(editions$verdict), collapse = ", "))
message("Isolates: ", nrow(raw), " (", sum(no_id), " without an ID); MIC results: ",
        nrow(tested), ", ", nrow(long), " after de-duplication (", flow_table$results_old_key[flow_table$source == "total"],
        " with the old key)")
message("Drugs with results: ", sum(drug_table$results_retail + drug_table$results_cecal > 0), " of ", nrow(drug_table))

# ---- 9. Panels: tested range by file, drug and year ----------------------------------------------
# A panel's limits show only in censored values: "<=x" is its lowest concentration and ">x" its
# highest. Per file, organism, drug and year, the distinct lower and upper limits seen (empty if
# no value was censored that way). A pair whose known limits differ between years changed panel.
limits <- function(m, signs) {
  v <- as.numeric(sub("^[<>=]+", "", m[sub("[0-9.].*$", "", m) %in% signs]))
  paste(format(sort(unique(v)), trim = TRUE, drop0trailing = TRUE), collapse = "/")
}
panel <- aggregate(mic ~ source + mo + ab + year, data = long, FUN = function(m) limits(m, c("<", "<=")))
names(panel)[names(panel) == "mic"] <- "lower"
panel$upper <- aggregate(mic ~ source + mo + ab + year, data = long, FUN = function(m) limits(m, c(">", ">=")))$mic
panel <- panel[order(panel$source, panel$mo, panel$ab, panel$year), ]
write.csv(panel, "results/narms_panels.csv", row.names = FALSE)
years_text <- function(y) {  # 2013 2014 2015 2018 -> "2013-2015, 2018"
  run <- cumsum(c(TRUE, diff(y) != 1))
  paste(tapply(y, run, function(r) if (length(r) == 1) r else paste0(min(r), "-", max(r))), collapse = ", ")
}
history <- function(v, y, sign) {
  k <- v != ""
  if (length(unique(v[k])) < 2) return(NA_character_)
  paste(vapply(unique(v[k]), function(u) sprintf("%s%s (%s)", sign, u, years_text(y[k & v == u])), character(1)), collapse = "; ")
}
pairs <- split(panel, paste(panel$source, panel$mo, panel$ab), drop = TRUE)
changed <- do.call(rbind, lapply(pairs, function(p)
  data.frame(source = p$source[1], mo = p$mo[1], ab = p$ab[1], lower = history(p$lower, p$year, "$\\leq$"),
             upper = history(p$upper, p$year, "$>$"))))
changed <- changed[!is.na(changed$lower) | !is.na(changed$upper), ]
org <- function(mo) ifelse(mo == "B_SLMNL_ENTR", "\\textit{Salmonella}", "\\textit{E.~coli}")
cat(c("\\begin{table}[h]", "\\centering", "\\scriptsize", "\\begin{tabular}{lllp{6.5cm}p{6.5cm}}", "\\toprule",
      "File & Organism & Drug & Lowest concentration, mg/L (years) & Highest concentration, mg/L (years) \\\\", "\\midrule",
      sprintf("%s & %s & %s & %s & %s \\\\", changed$source, org(changed$mo), AMR::ab_name(changed$ab, info = FALSE),
              ifelse(is.na(changed$lower), "unchanged", changed$lower), ifelse(is.na(changed$upper), "unchanged", changed$upper)),
      "\\bottomrule", "\\end{tabular}",
      sprintf(paste("\\caption{Panel changes, read from censored values: ``$\\leq$x'' marks the lowest concentration",
                    "tested and ``$>$x'' the highest. File, organism and drug pairs whose limits differed between years",
                    "(%d of %d pairs); a/b: both limits in one year. Years in which no value was censored at a limit",
                    "are not listed for it.}"), nrow(changed), length(pairs)),
      "\\end{table}", ""), file = "results/supplement_tables.tex", sep = "\n", append = TRUE)

# ---- 10. Clustering --------------------------------------------------------------------------
# The cecal file records one isolate per organism per sample and no establishment or state, so
# the only grouping is the month of collection (about a dozen per period), which is not the
# likely unit of dependence (plant, flock or herd). Resampling whole year-months, G - 1 of G each
# time, shows how much the intervals widen if isolates collected together are correlated.
# Retail records the state: the retail MIC shift is refitted on the states sampled in both
# periods, with a bootstrap that resamples whole states (G - 1 of G).
suppressPackageStartupMessages(library(survival))  # survreg() finds strata() by name
set.seed(2026)
relabel_cl <- drift_bootstrap(relabel_data, from = 2022, to = 2023, B = 5000, cluster = "batch")
editions_cl <- drift_bootstrap(cecal, from = 2019:2020, to = 2022:2023, B = 5000, cluster = "batch")
cluster_table <- function(iso, cl, comparison) {
  m <- merge(iso[!is.na(iso$real), c("mo", "ab", "real", "real_low", "real_high", "artefact", "artefact_low", "artefact_high", "verdict")],
             cl[c("mo", "ab", "real_low", "real_high", "artefact_low", "artefact_high", "boot_p_real", "units_y1", "units_y2")],
             by = c("mo", "ab"), suffixes = c("", "_month"))
  cbind(comparison = comparison, m)
}
clustering <- rbind(cluster_table(relabel, relabel_cl, "2022 v 2023"),
                    cluster_table(editions, editions_cl, "2019-2020 v 2022-2023"))
names(clustering)[names(clustering) %in% c("units_y1", "units_y2")] <- c("months_y1", "months_y2")
write.csv(round_numbers(clustering), "results/narms_cluster.csv", row.names = FALSE)
widen <- function(lo, hi, lo2, hi2) ifelse(hi - lo > 0, (hi2 - lo2) / (hi - lo), NA)
clustering$ratio_real <- widen(clustering$real_low, clustering$real_high, clustering$real_low_month, clustering$real_high_month)
clustering$ratio_artefact <- widen(clustering$artefact_low, clustering$artefact_high, clustering$artefact_low_month, clustering$artefact_high_month)
excl0 <- function(lo, hi) !is.na(lo) & (lo > 0 | hi < 0)

shift_data <- function(d, from, to) {
  d <- d[d$year %in% c(from, to) & !is.na(d$mic), ]
  d$late <- d$year %in% to
  d <- cbind(d, mic_interval(d$mic))
  d[!(is.na(d$lo) & is.na(d$hi)), ]
}
retail_shift <- shift_data(retail, 2002:2006, 2017:2021)
common_states <- intersect(retail_shift$state[!retail_shift$late], retail_shift$state[retail_shift$late])
fit_late <- function(x) {
  f <- tryCatch(survival::survreg(survival::Surv(lo, hi, type = "interval2") ~ late, data = x, dist = "gaussian"),
                error = function(e) NULL, warning = function(w) NULL)
  if (is.null(f)) NA else stats::coef(f)[["lateTRUE"]]
}
state_boot <- function(x, B = 999) {
  by_period <- split(split(x, x$state), vapply(split(x, x$state), function(d) d$late[1], logical(1)))
  draws <- replicate(B, fit_late(do.call(rbind, unlist(lapply(by_period, function(st)
    st[sample.int(length(st), max(1, length(st) - 1), replace = TRUE)]), recursive = FALSE))))
  stats::quantile(draws, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
}
est <- shift_retail[!is.na(shift_retail$shift), ]
set.seed(2026)
state_ci <- t(vapply(seq_len(nrow(est)), function(i) {
  x <- retail_shift[retail_shift$mo == est$mo[i] & retail_shift$ab == est$ab[i] & retail_shift$state %in% common_states, ]
  # each state appears in both periods; split it into its early and late part so that a state is
  # resampled within each period
  x$state <- paste(x$state, x$late)
  c(fit_late(x), state_boot(x))
}, numeric(3)))
est$shift_common_states <- state_ci[, 1]
est$shift_low_state <- state_ci[, 2]
est$shift_high_state <- state_ci[, 3]
est$states_from <- length(unique(retail_shift$state[!retail_shift$late]))
est$states_to <- length(unique(retail_shift$state[retail_shift$late]))
est$states_common <- length(common_states)
write.csv(round_numbers(est[c("mo", "ab", "shift", "shift_low", "shift_high", "shift_common_states", "shift_low_state",
                              "shift_high_state", "states_from", "states_to", "states_common")]),
          "results/narms_shift_state.csv", row.names = FALSE)

ci1 <- function(lo, hi) ifelse(is.na(lo), "", ifelse(lo == 0 & hi == 0, "0", sprintf("%s to %s", pp(lo), pp(hi))))
rows <- with(clustering, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %d / %d \\\\", comparison, org(mo),
                                 AMR::ab_name(ab, info = FALSE), pp(real), ci1(real_low, real_high),
                                 ci1(real_low_month, real_high_month), pp(artefact), ci1(artefact_low, artefact_high),
                                 ci1(artefact_low_month, artefact_high_month), months_y1, months_y2))
cat(table_tex("lllrllrllr",
      "Comparison & Organism & Drug & Fixed criteria & Isolates & Months & Criteria & Isolates & Months & Months (G) \\\\",
      rows,
      paste("Cluster sensitivity, cecal isolates, comparisons of Tables S1 and S3. Fixed-criteria and criteria components (percentage points) with",
            "95\\% bootstrap intervals that resample isolates (as in Tables S1 and S3) or whole collection months",
            "(G $-$ 1 of the G months of each period per replicate), 5,000 replicates each. G: months per period.",
            "Interpretable pairs only."), long = TRUE), file = "results/supplement_tables.tex", sep = "\n", append = TRUE)

# ---- 11. Fit of the MIC-shift model ------------------------------------------------------------
# For every estimated shift: the same model with a separate standard deviation per period
# (likelihood-ratio test against one shared SD), and the largest gap between the fitted normal
# CDF and the observed share at or below each stored dilution, per period.
model_check <- function(d, shifts, label) {
  est <- shifts[!is.na(shifts$shift), ]
  do.call(rbind, lapply(seq_len(nrow(est)), function(i) {
    x <- d[d$mo == est$mo[i] & d$ab == est$ab[i], ]
    f1 <- tryCatch(survival::survreg(survival::Surv(lo, hi, type = "interval2") ~ late, data = x, dist = "gaussian"),
                   error = function(e) NULL, warning = function(w) NULL)
    if (is.null(f1)) return(NULL)  # drift_mic_shift() fitted the same model, so this does not happen
    f2 <- tryCatch(survival::survreg(survival::Surv(lo, hi, type = "interval2") ~ late + strata(late), data = x,
                                     dist = "gaussian"), error = function(e) NULL, warning = function(w) NULL)
    gap <- vapply(c(FALSE, TRUE), function(l) {
      xi <- x[x$late == l, ]
      mu <- sum(stats::coef(f1) * c(1, l))
      h <- sort(unique(xi$hi[!is.na(xi$hi)]))
      # values whose interval straddles h (a capped value from another panel) are left out at h
      max(vapply(h, function(k) {
        known <- (!is.na(xi$hi) & xi$hi <= k) | (!is.na(xi$lo) & xi$lo >= k)
        abs(mean(!is.na(xi$hi[known]) & xi$hi[known] <= k) - stats::pnorm((k - mu) / f1$scale))
      }, numeric(1)))
    }, numeric(1))
    data.frame(data = label, mo = est$mo[i], ab = est$ab[i], shift = est$shift[i],
               shift_unequal = if (is.null(f2)) NA else stats::coef(f2)[["lateTRUE"]],
               sd = f1$scale, sd_from = if (is.null(f2)) NA else f2$scale[1], sd_to = if (is.null(f2)) NA else f2$scale[2],
               p_equal_sd = if (is.null(f2)) NA else stats::pchisq(2 * (f2$loglik[2] - f1$loglik[2]), 1, lower.tail = FALSE),
               max_gap_from = gap[1], max_gap_to = gap[2])
  }))
}
cecal_shift <- shift_data(cecal, 2019:2020, 2022:2023)
mic_model <- rbind(model_check(cecal_shift, shift_cecal, "cecal 2019-2020 v 2022-2023"),
                   model_check(retail_shift, shift_retail, "retail 2002-2006 v 2017-2021"))
write.csv(round_numbers(mic_model), "results/narms_mic_model.csv", row.names = FALSE)
f2 <- function(v) ifelse(is.na(v), "", sprintf("%.2f", v))
cat(c("\\begin{table}[h]", "\\centering", "\\scriptsize", "\\begin{tabular}{lllrrrrrrr}", "\\toprule",
      "Data & Organism & Drug & Shift & Shift, own SD & SD & SD before & SD after & $p$ equal SD & Largest CDF gap \\\\",
      "\\midrule",
      with(mic_model, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", data, org(mo), AMR::ab_name(ab, info = FALSE),
                              f2(shift), f2(shift_unequal), f2(sd), f2(sd_from), f2(sd_to),
                              ifelse(is.na(p_equal_sd), "", ifelse(p_equal_sd < 0.001, "$<$0.001", sprintf("%.3f", p_equal_sd))),
                              f2(pmax(max_gap_from, max_gap_to)))),
      "\\bottomrule", "\\end{tabular}",
      paste("\\caption{Fit of the MIC-shift model, for every estimated shift (log$_2$ dilutions). Shift, own SD: the",
            "period effect when each period has its own standard deviation; $p$ equal SD: likelihood-ratio test of",
            "one shared SD (unadjusted). Largest CDF gap: the largest absolute difference, in either period,",
            "between the share of values at or below a stored dilution and the fitted normal share.}"),
      "\\end{table}", ""), file = "results/supplement_tables.tex", sep = "\n", append = TRUE)

# ---- 12. Further supplementary tables (S9-S12) and the numbers quoted in the text --------------
extra <- character(0)
orderings <- function(x, comparison) {
  x <- x[!is.na(x$real), ]
  rep <- !is.null(x$p_rep_y1)
  # Equation 2: the guideline component as the revision term plus the two reporting gaps
  revision <- x$p_y2g2 - x$p_y1g1 - x$real
  gap1 <- if (rep) x$p_rep_y1 - x$p_y1g1 else NA
  gap2 <- if (rep) x$p_rep_y2 - x$p_y2g2 else NA
  sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", comparison, org(x$mo), AMR::ab_name(x$ab, info = FALSE),
          pp(x$real_at_new), pp(x$real_at_old), pp(x$interaction), pp(x$artefact), pp(revision),
          if (rep) pp(gap1) else "", if (rep) pp(gap2) else "")
}
extra <- c(extra, table_tex("lllrrrrrrr",
  "Comparison & Organism & Drug & Fixed at $g_2$ & Fixed at $g_1$ & Interaction & Criteria & Revision & Gap before & Gap after \\\\",
  c(orderings(relabel, "S1"), orderings(relabel_ns, "S2"), orderings(editions, "S3")),
  paste("Comparisons of Tables S1-S3, interpretable pairs. The two orderings of the fixed-criteria component, $p(y_2, g_2) - p(y_1, g_2)$ and $p(y_2, g_1) - p(y_1, g_1)$,",
        "their difference (the interaction, which the Shapley mean splits equally by convention), and the",
        "criteria component. For reported calls (S1, S2) the criteria component is split as in Equation 2 of the main text:",
        "Revision, $p(y_2, g_2) - p(y_1, g_1) -$ fixed-criteria component; Gap before and after, $r(y) - p(y, g)$, the reported",
        "share minus the tool's reading, so criteria = revision + gap after $-$ gap before. Percentage points."), long = TRUE))
w <- windows
extra <- c(extra, table_tex("llllllrlll",
  "Comparison & Before & After & Organism & Drug & Editions & Observed & Fixed criteria (95\\% CI) & Criteria (95\\% CI) & Verdict \\\\",
  sprintf("%s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\", ifelse(w$comparison == "relabel", "Reported", "Editions"),
          w$year_from, w$year_to, org(w$mo), AMR::ab_name(w$ab, info = FALSE),
          paste(sub("CLSI ", "", w$guideline_from), "/", sub("CLSI ", "", w$guideline_to)), pp(w$observed),
          ci(w$real, w$real_low, w$real_high), ci(w$artefact, w$artefact_low, w$artefact_high), term(w$verdict)),
  paste("Sensitivity to the choice of windows, cecal isolates. Reported: \\%R as reported by NARMS (as Table S1),",
        "ciprofloxacin; Editions: each period read with the CLSI edition in force in its last year (as Table S3),",
        "gentamicin. Bootstrap intervals and tests, 2,000 replicates; verdicts with the correction over all",
        "interpretable pairs of that window. The first row of each comparison repeats the main-text analysis with",
        "its own 2,000 replicates, so its intervals can differ from the main text by a tenth of a point.")))
f <- fallback_check
extra <- c(extra, table_tex("llrrrl",
  "Organism & Drug & n & Criteria, CLSI 2019 unfiltered & Criteria, CLSI 2020 as older edition & Interpretable with filter \\\\",
  sprintf("%s & %s & %s & %s & %s & %s \\\\", org(f$mo), AMR::ab_name(f$ab, info = FALSE), format(f$n, big.mark = ","),
          pp(f$artefact_clsi2019_unfiltered), pp(f$artefact_clsi2020), ifelse(f$interpretable_clsi2019_filtered, "yes", "no")),
  paste("The generic Gram-negative fallback, cecal isolates 2018-2019 against 2022-2023 (CLSI 2023). Criteria",
        "component (percentage points) with the older period read under CLSI 2019 as AMR 3.1.0 reads it, without",
        "the fallback filter, and with CLSI 2020, which has the Enterobacterales rows, as the older edition. n:",
        "isolates tested in 2018-2019. Last column: whether the pair is interpretable under CLSI 2019 with the",
        "filter the tool applies.")))
st <- standardised
extra <- c(extra, table_tex("lrrrl",
  "Standardisation & \\%NS 2022 & \\%NS 2023 & Change (95\\% CI) & Isolates in common strata \\\\",
  sprintf("%s & %.1f & %.1f & %s & %s \\\\", c(crude = "Crude", host = "Host", serotype = "Serotype", joint = "Host and serotype")[st$by],
          st$pct_2022, st$pct_2023,
          ifelse(is.na(st$change_low), sprintf("%+.1f", st$change), sprintf("%+.1f (%+.1f to %+.1f)", st$change, st$change_low, st$change_high)),
          sprintf("%.1f\\%%", 100 * st$share_in_common_strata)),
  paste("Cecal \\textit{Salmonella} ciprofloxacin, share non-susceptible under CLSI 2023, 2022 and 2023, directly",
        "standardised to the pooled mix of both years over the strata present in both years. Serotypes with under",
        "20 isolates pooled as Other. Bootstrap intervals, isolates resampled within year, 2,000 replicates.")))
writeLines(extra, "results/supplement_extra.tex")

# Numbers quoted in the paper that no table above holds
sal_gen <- cecal_shift[cecal_shift$mo == "B_SLMNL_ENTR" & cecal_shift$ab == "GEN", ]
share_at <- function(x, at_most = NULL, at_least = NULL) {
  if (!is.null(at_most)) mean(!is.na(x$hi) & x$hi <= log2(at_most)) else mean(!is.na(x$lo) & x$lo >= log2(at_least) - 1)
}
sal <- ns[ns$year == 2022, ]
retail_gen <- long[long$source == "retail" & long$ab == "GEN" & long$year %in% 2019:2021 & long$mic %in% c("4", "8"), ]
state_ratio <- (est$shift_high_state - est$shift_low_state) / (est$shift_high - est$shift_low)
zero <- function(x) sum(!is.na(x$real) & abs(x$artefact) < 5e-4)
numbers <- c(
  sample_records = nrow(records), no_growth_records = sum(!records$growth %in% "YES"), isolates = nrow(raw),
  isolates_with_mic = sum(has_mic), isolates_with_mic_cecal = sum(has_mic & raw$source == "cecal"),
  isolates_with_mic_retail = sum(has_mic & raw$source == "retail"), isolates_no_id = sum(no_id),
  mic_results_dedup = nrow(long),
  sal_gen_le025_2019_2020 = 100 * share_at(sal_gen[!sal_gen$late, ], at_most = 0.25),
  sal_gen_le025_2022_2023 = 100 * share_at(sal_gen[sal_gen$late, ], at_most = 0.25),
  sal_gen_ge8_2019_2020 = 100 * share_at(sal_gen[!sal_gen$late, ], at_least = 8),
  sal_gen_ge8_2022_2023 = 100 * share_at(sal_gen[sal_gen$late, ], at_least = 8),
  sal_cip_chicken_share_2022 = 100 * mean(ns$host[ns$year == 2022] == "Chickens"),
  sal_cip_chicken_share_2023 = 100 * mean(ns$host[ns$year == 2023] == "Chickens"),
  sal_cip_chicken_ns_2022 = 100 * mean(sal$hit[sal$host == "Chickens"]),
  sal_cip_infantis_share_2022 = 100 * mean(ns$serotype[ns$year == 2022] == "Infantis"),
  sal_cip_infantis_share_2023 = 100 * mean(ns$serotype[ns$year == 2023] == "Infantis"),
  sal_cip_infantis_ns_2022 = 100 * mean(sal$hit[sal$serotype == "Infantis"]),
  sal_cip_infantis_ns_2023 = 100 * mean(ns$hit[ns$year == 2023 & ns$serotype == "Infantis"]),
  retail_gen_2019_2021_mic4_called_i = mean(retail_gen$narms[retail_gen$mic == "4"] %in% "I"),
  retail_gen_2019_2021_mic8_called_r = mean(retail_gen$narms[retail_gen$mic == "8"] %in% "R"),
  state_ci_ratio_median = stats::median(state_ratio, na.rm = TRUE),
  state_ci_ratio_min = min(state_ratio, na.rm = TRUE), state_ci_ratio_max = max(state_ratio, na.rm = TRUE),
  month_ratio_period_median = stats::median(clustering$ratio_real, na.rm = TRUE),
  month_ratio_period_min = min(clustering$ratio_real, na.rm = TRUE), month_ratio_period_max = max(clustering$ratio_real, na.rm = TRUE),
  month_ratio_guideline_median = stats::median(clustering$ratio_artefact, na.rm = TRUE),
  month_ratio_guideline_min = min(clustering$ratio_artefact, na.rm = TRUE),
  month_ratio_guideline_max = max(clustering$ratio_artefact, na.rm = TRUE),
  period_excl0_isolates = sum(excl0(clustering$real_low, clustering$real_high)),
  period_excl0_months = sum(excl0(clustering$real_low_month, clustering$real_high_month)),
  period_excl0_lost = sum(excl0(clustering$real_low, clustering$real_high) & !excl0(clustering$real_low_month, clustering$real_high_month)),
  period_excl0_gained = sum(!excl0(clustering$real_low, clustering$real_high) & excl0(clustering$real_low_month, clustering$real_high_month)),
  guideline_excl0_isolates = sum(excl0(clustering$artefact_low, clustering$artefact_high)),
  guideline_excl0_months = sum(excl0(clustering$artefact_low_month, clustering$artefact_high_month)),
  pairs_interpretable_s1 = sum(!is.na(relabel$real)), pairs_zero_guideline_s1 = zero(relabel),
  pairs_interpretable_s2 = sum(!is.na(relabel_ns$real)), pairs_zero_guideline_s2 = zero(relabel_ns),
  pairs_interpretable_s3 = sum(!is.na(editions$real)), pairs_zero_guideline_s3 = zero(editions),
  pairs_same_isolates_both_guidelines = sum(with(rbind(relabel[c("n_y1g1", "n_y1g2", "n_y2g1", "n_y2g2", "real")],
                                                         relabel_ns[c("n_y1g1", "n_y1g2", "n_y2g1", "n_y2g2", "real")],
                                                         editions[c("n_y1g1", "n_y1g2", "n_y2g1", "n_y2g2", "real")]),
                                                   !is.na(real) & n_y1g1 == n_y1g2 & n_y2g1 == n_y2g2)),
  retail_unread_calls_2020_2021 = sum(agreement$calls[agreement$source == "retail"] - agreement$tests[agreement$source == "retail"]),
  cecal_unread_calls_2023 = sum(agreement$calls[agreement$source == "cecal" & agreement$year == 2023] -
                                agreement$tests[agreement$source == "cecal" & agreement$year == 2023]))
for (v in verdict_levels) numbers[paste0("verdicts_s1_", gsub(" ", "_", v))] <- sum(relabel$verdict == v)
for (v in verdict_levels) numbers[paste0("verdicts_s3_", gsub(" ", "_", v))] <- sum(editions$verdict == v)
write.csv(data.frame(name = names(numbers), value = signif(unname(numbers), 6)), "results/narms_numbers.csv", row.names = FALSE)

message("Cluster sensitivity: median interval ratio (months / isolates), period ",
        round(stats::median(clustering$ratio_real, na.rm = TRUE), 2), ", guideline ",
        round(stats::median(clustering$ratio_artefact, na.rm = TRUE), 2), "; period intervals excluding 0: ",
        sum(excl0(clustering$real_low, clustering$real_high)), " with isolates, ",
        sum(excl0(clustering$real_low_month, clustering$real_high_month)), " with months; guideline: ",
        sum(excl0(clustering$artefact_low, clustering$artefact_high)), " and ",
        sum(excl0(clustering$artefact_low_month, clustering$artefact_high_month)))
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
