# Run the drift audit on real MICs from NCBI BioSample.
# Run from the project root: Rscript scripts/ncbi.R
#
# data/raw/ncbi_ast_long.csv.gz is the snapshot made by R/00_fetch_ncbi.R in
# github.com/FerrumSL/amr-surveillance-r (human isolates with an NCBI antibiogram table, one row
# per isolate and test, as published). It is copied here unchanged.
#
# These isolates show that the method runs on real, messy data. They cannot show real trends:
# submitters chose what to sequence (often because it was resistant), and which labs submitted
# changes from year to year. The proof that the split is right is in tests/testthat.

source("R/drift.R")

long <- read.csv("data/raw/ncbi_ast_long.csv.gz", colClasses = "character", na.strings = "")

# CLSI MICs only: the guideline in force is looked up by year within one organisation. NCBI
# gives only "CLSI", not the edition, so the edition is assumed from the collection year.
d <- long[long$method %in% "MIC" & tolower(long$units) %in% "mg/l" & long$standard %in% "CLSI" &
            !is.na(long$measurement), ]
# Dates come as a day, a month, a year or a range; a range keeps its start.
d$year <- as.integer(ifelse(grepl("^\\d{4}", d$collection_date), substr(d$collection_date, 1, 4), NA))
d$mo <- as.character(AMR::as.mo(d$organism, info = FALSE))
d$ab <- as.character(suppressWarnings(AMR::as.ab(d$antibiotic, info = FALSE)))
d$mic <- paste0(ifelse(d$sign %in% c("==", "="), "", d$sign), d$measurement)
d <- d[!is.na(d$year) & !is.na(d$ab), ]
d <- d[!duplicated(d[c("accession", "ab")]), ]  # a few isolates repeat a test

# Each BioSample is one isolate from its own patient, so there is nothing for
# AMR::first_isolate() to remove. Two periods either side of CLSI's 2020 edition, which cut the
# Enterobacterales fluoroquinolone breakpoints; 2011 is the oldest edition in the AMR package.
from <- 2011:2018
to <- 2019:2025
d <- d[d$year %in% c(from, to), ]

result <- drift_decompose(d[c("year", "mo", "ab", "mic")], from = from, to = to)
# Pairs with fewer than 10 readable isolates in either period are kept in the CSV but not tested
# or reported, so the multiple-comparison correction runs over the reported pairs only.
enough <- result$n_y1g2 >= 10 & result$n_y2g2 >= 10
reported <- drift_test(result[enough, ])
result <- merge(result, reported, all.x = TRUE)

dir.create("results", showWarnings = FALSE)
write.csv(result, "results/ncbi_drift.csv", row.names = FALSE)
drift_report(reported, "results/ncbi_report.md",
             title = "Breakpoint drift: NCBI BioSample isolates (not surveillance data)",
             note = sprintf(paste(
               "%d isolates and %d CLSI MIC tests from NCBI BioSample. Shown: the %d organism x drug",
               "pairs with at least 10 readable isolates in both periods (all %d are in",
               "`ncbi_drift.csv`). Submitters chose these isolates, often because they were",
               "resistant, so the rates describe no population. The edition in force is assumed from",
               "the collection year."),
               length(unique(d$accession)), nrow(d), nrow(reported), nrow(result)))

message("Isolates: ", length(unique(d$accession)), "; tests: ", nrow(d), "; pairs: ", nrow(result),
        "; reported: ", nrow(reported))
print(table(reported$verdict))
