# A real breakpoint revision in reported calls: CLSI's 2019 fluoroquinolone revision for
# Enterobacterales (M100 29th edition), in NCBI BioSample antibiograms.
# Run from the project root: Rscript scripts/ncbi_fluoroquinolones.R
#
# data/raw/ncbi_enterobacterales_ast.csv.gz: every Enterobacterales BioSample with an antibiogram,
# fetched on 2026-11-17 by scripts/fetch_ncbi_enterobacterales.py (one row per isolate and test,
# as published, with the submitting lab's own call). Unlike data/raw/ncbi_ast_long.csv.gz, which
# holds 400 isolates per species, it holds all of them.
#
# The AMR package has no Enterobacterales ciprofloxacin or levofloxacin row before its CLSI 2020,
# so the pre-revision breakpoints come from data/breakpoints_clsi_fluoroquinolones.csv. The calls
# the labs reported are the observed series. Submitters chose what to sequence, so the change at
# fixed criteria describes no population; reclassification compares the same isolates read twice
# and does not depend on how they were chosen.

source("R/drift.R")

long <- read.csv("data/raw/ncbi_enterobacterales_ast.csv.gz", colClasses = "character", na.strings = "")
breakpoints <- read.csv("data/breakpoints_clsi_fluoroquinolones.csv")

d <- long[tolower(long$antibiotic) %in% c("ciprofloxacin", "levofloxacin") & long$standard %in% "CLSI" &
            long$method %in% "MIC" & long$units %in% "mg/L" & !is.na(long$measurement) &
            long$host %in% "Homo sapiens" & !startsWith(long$organism, "Salmonella"), ]
# Salmonella has its own fluoroquinolone breakpoints (revised in 2012), so it is left out.
# The year is the submission year: the call is the one the lab reported then. Labs adopted the
# 2019 edition late (calls in 2019-2022 mix both), so the periods are 2014-2018 and 2023-2025.
d$year <- as.integer(substr(d$submission_date, 1, 4))
d$mo <- as.character(AMR::as.mo(d$organism, info = FALSE))
d$ab <- as.character(AMR::as.ab(d$antibiotic, info = FALSE))
d$mic <- paste0(ifelse(d$sign %in% c("==", "="), "", d$sign), d$measurement)
d$reported <- c(susceptible = "S", intermediate = "I", resistant = "R")[d$phenotype]  # others NA
d <- d[d$mo != "UNKNOWN" & !duplicated(d[c("accession", "ab")]), ]
from <- 2014:2018
to <- 2023:2025
d <- d[d$year %in% c(from, to), ]

cols <- c("year", "mo", "ab", "mic", "reported")
set.seed(2026)
by_species <- drift_bootstrap(cbind(d[cols], isolate = d$accession), from = from, to = to, B = 2000,
                              breakpoints = breakpoints)
# All Enterobacterales together, as a cumulative antibiogram would pool them
pooled <- d
pooled$mo <- "B_[ORD]_ENTRBCTR"
set.seed(2027)
pooled <- drift_bootstrap(cbind(pooled[cols], isolate = pooled$accession), from = from, to = to, B = 2000,
                          breakpoints = breakpoints)
result <- rbind(pooled, by_species)
enough <- result$n_y1g2 >= 10 & result$n_y2g2 >= 10
result <- merge(result, drift_test(result[enough, ]), all.x = TRUE, sort = FALSE)

dir.create("results", showWarnings = FALSE)
num <- vapply(result, is.numeric, logical(1))
result[num] <- lapply(result[num], signif, 6)
write.csv(result, "results/ncbi_fluoroquinolones.csv", row.names = FALSE)

# The calls at the MIC values the revision moved, as reported in each period
moved <- d[d$mic %in% c("0.5", "1", "2", "4"), ]
moved$period <- ifelse(moved$year %in% from, "2014-2018", "2023-2025")
calls <- as.data.frame(table(ab = moved$ab, mic = moved$mic, period = moved$period, reported = moved$reported))
calls <- calls[calls$Freq > 0, ]
write.csv(calls, "results/ncbi_fluoroquinolones_calls.csv", row.names = FALSE)

pct <- function(p) sprintf("%.1f", 100 * p)
ci <- function(x, lo, hi) sprintf("%+.1f (%+.1f to %+.1f)", 100 * x, 100 * lo, 100 * hi)
shown <- result[!is.na(result$verdict), ]
message("Isolates: ", length(unique(d$accession)), "; tests: ", nrow(d))
print(data.frame(mo = AMR::mo_name(shown$mo, info = FALSE), ab = shown$ab, n1 = shown$tested_y1, n2 = shown$tested_y2,
                 reported_1 = pct(shown$p_rep_y1), reported_2 = pct(shown$p_rep_y2),
                 observed = ci(shown$observed, shown$observed_low, shown$observed_high),
                 fixed = ci(shown$fixed, shown$fixed_low, shown$fixed_high),
                 reclass = ci(shown$reclass, shown$reclass_low, shown$reclass_high), verdict = shown$verdict),
      row.names = FALSE)
