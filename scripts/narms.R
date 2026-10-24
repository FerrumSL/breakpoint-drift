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
#
# Four results:
#   1. NARMS's own relabelling: up to 2022 NARMS called a Salmonella ciprofloxacin MIC of 0.12
#      or more R; from 2023 it uses CLSI (0.12-0.5 is I). Each drug's reported %R in 2022 and
#      2023 is split into real change (both years read with CLSI 2023) and artefact.
#   2. The tool's own decomposition, CLSI 2020 against CLSI 2023, on cecal isolates.
#   3. Agreement of the tool's reading with NARMS's calls, by drug and year.
#   4. The MIC shift on retail isolates over 20 years, which needs no breakpoint at all.

source("R/drift.R")

extract <- "data/raw/narms_extract.csv.gz"
# NARMS code -> AMR code. Spelled out because AMR::as.ab() reads some NARMS codes as other
# drugs (PTZ as pentizidone, TET as tetroxoprim).
drugs <- c(AMC = "AMC", AMI = "AMK", AMP = "AMP", ATM = "ATM", AXO = "CRO", AZI = "AZM",
           CAZ = "CAZ", CHL = "CHL", CIP = "CIP", COL = "COL", COT = "SXT", CTX = "CTX",
           FEP = "FEP", FOX = "FOX", GEN = "GEN", IMI = "IPM", MER = "MEM", NAL = "NAL",
           PTZ = "TZP", TET = "TCY", TGC = "TGC")

if (!file.exists(extract)) {
  read_source <- function(path, source) {
    d <- readxl::read_excel(path, col_types = "text")
    names(d)[names(d) %in% c("Year", "YEAR")] <- "year"
    d <- d[d$GENUS_NAME %in% c("Salmonella", "Escherichia"), ]
    out <- data.frame(source = source, isolate = d$ISOLATE_ID, year = as.integer(d$year),
                      genus = d$GENUS_NAME, host = d$HOST_SPECIES)
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
  write.csv(raw, gzfile(extract), row.names = FALSE, na = "")
}

raw <- read.csv(extract, colClasses = "character", na.strings = "")
raw$year <- as.integer(raw$year)
raw$mo <- ifelse(raw$genus == "Salmonella", "B_SLMNL_ENTR", "B_ESCHR_COLI")

# One row per isolate and drug. NARMS panels print 2^-6 as 0.015 and 2^-4 as 0.06; the AMR
# package writes them 0.016 and 0.064.
amr_level <- c("0.015" = "0.016", "0.03" = "0.032", "0.06" = "0.064", "0.12" = "0.125")
codes <- intersect(names(drugs), names(raw))
long <- do.call(rbind, lapply(codes, function(code) {
  x <- raw[!is.na(raw[[code]]), ]
  if (nrow(x) == 0) return(NULL)
  number <- sub("^[<>=]+", "", x[[code]])
  data.frame(source = x$source, isolate = x$isolate, year = x$year, mo = x$mo, ab = drugs[[code]],
             mic = paste0(sub("[0-9.].*$", "", x[[code]]), ifelse(number %in% names(amr_level), amr_level[number], number)),
             narms = x[[paste0(code, "_sir")]])
}))
long <- long[!duplicated(long[c("source", "isolate", "ab")]), ]
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
read <- recent[!is.na(recent$ours) & recent$ours != "NI", ]
agreement <- aggregate(cbind(tests = 1, agree = read$ours == read$narms) ~ source + mo + ab + year,
                       data = read, FUN = sum)
agreement$agree_pct <- 100 * agreement$agree / agreement$tests
write.csv(round_numbers(agreement), "results/narms_agreement.csv", row.names = FALSE)

# ---- 1. NARMS's own relabelling, 2022 to 2023 (cecal) --------------------------------------------
# Observed: %R as NARMS reported it each year. Real: the MICs, read with CLSI 2022 and CLSI 2023
# (Shapley mean of holding either). The artefact is everything else: NARMS's own rules and the
# CLSI revisions. NARMS cecal calls through 2022 treat a Salmonella or E. coli ciprofloxacin MIC
# of 0.12 or more as R; in 2023 they follow CLSI, except that they keep the gentamicin breakpoint
# CLSI cut in 2023 (the agreement table shows both).
relabel_data <- long[long$source == "cecal" & long$year %in% 2022:2023, ]
relabel_data$reported <- relabel_data$narms
relabel <- drift_decompose(relabel_data, from = 2022, to = 2023)
relabel <- drift_test(relabel[relabel$n_rep_y1 >= 10 & relabel$n_rep_y2 >= 10, ])
write.csv(round_numbers(relabel), "results/narms_relabel.csv", row.names = FALSE)

# ---- 2. CLSI 2020 against CLSI 2023 (cecal) ------------------------------------------------------
# Between these editions CLSI cut the Enterobacterales aminoglycoside breakpoints (2023) and
# piperacillin/tazobactam (2022).
cecal <- long[long$source == "cecal", ]
editions <- drift_decompose(cecal, from = 2019:2020, to = 2022:2023)
editions <- drift_test(editions[editions$n_y1g2 >= 10 & editions$n_y2g2 >= 10, ])
shift_cecal <- drift_mic_shift(cecal, from = 2019:2020, to = 2022:2023)
editions <- merge(editions, shift_cecal, by = c("mo", "ab"), all.x = TRUE)
write.csv(round_numbers(editions), "results/narms_editions.csv", row.names = FALSE)

# ---- 4. Twenty years of MIC shift (retail) -------------------------------------------------------
retail <- long[long$source == "retail", ]
shift_retail <- drift_mic_shift(retail, from = 2002:2006, to = 2017:2021)
write.csv(round_numbers(shift_retail), "results/narms_mic_shift.csv", row.names = FALSE)

# ---- Report --------------------------------------------------------------------------------------
note <- function(...) paste(...)
drift_report(relabel, "results/narms_relabel.md",
             title = "NARMS cecal isolates: NARMS's 2023 relabelling",
             note = note(
               "Up to 2022 NARMS reported cecal isolates with its own rules, e.g. a Salmonella",
               "ciprofloxacin MIC of 0.12 or more as R; from 2023 it follows CLSI, but keeps the",
               "gentamicin breakpoint CLSI cut in 2023. **Observed** is %R as NARMS reported it;",
               "**real** reads the MICs with CLSI 2022 and CLSI 2023."))
drift_report(editions, "results/narms_editions.md",
             title = "NARMS cecal isolates: CLSI 2020 against CLSI 2023",
             note = note(sprintf("%d Salmonella and E. coli isolates from food-animal guts at slaughter.",
                                 length(unique(cecal$isolate[cecal$year %in% c(2019:2020, 2022:2023)])))))

agree_total <- aggregate(cbind(tests, agree) ~ mo + ab, data = agreement[agreement$source == "cecal" & agreement$year == 2023, ], FUN = sum)
message("Isolates in extract: ", nrow(raw), "; isolate x drug tests: ", nrow(long))
message("Agreement with NARMS calls, cecal 2023: ", sum(agree_total$agree), " / ", sum(agree_total$tests))
message("Relabelling verdicts: ", paste(names(table(relabel$verdict)), table(relabel$verdict), collapse = ", "))
message("Edition verdicts: ", paste(names(table(editions$verdict)), table(editions$verdict), collapse = ", "))
