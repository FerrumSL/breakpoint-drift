# Audit of the CLSI Enterobacterales rows in WHONET's Breakpoints.txt and in the AMR package.
# Run from the project root: Rscript scripts/whonet_audit.R
#
# AMR imports its CLSI breakpoints from WHONET's AMRIE file (AMR issue #334, AMRIE issue #27).
# This script counts, per CLSI edition 2011-2020, the Table 2A (Enterobacterales) rows in both,
# and shows with as.sir() what an E. coli MIC is read with when the family rows are missing.
# The WHONET file is pinned to AMRIE commit 69de201 and checked against its SHA-256; it is
# downloaded to data/whonet/ (not committed: 5 MB) on the first run.

file <- "data/whonet/Breakpoints.txt"
url <- "https://raw.githubusercontent.com/WHONET-Development/AMRIE/69de201/Interpretation%20Engine/Resources/Breakpoints.txt"
sha <- "c9c6e72844d144941051e0fa1919ae9698846993d21cce043693d0a25092d202"
if (!file.exists(file)) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  utils::download.file(url, file, mode = "wb", quiet = TRUE)
}
got <- unname(tools::sha256sum(file))
if (got != sha) stop(file, " has SHA-256 ", got, ", expected ", sha, ".")

editions <- 2011:2020
who <- utils::read.delim(file, colClasses = "character", quote = "")
who <- who[who$GUIDELINES == "CLSI" & who$YEAR %in% editions & who$BREAKPOINT_TYPE == "Human", ]
amr <- AMR::clinical_breakpoints
amr <- amr[amr$guideline %in% paste("CLSI", editions) & amr$type == "human", ]
amr$year <- as.integer(sub("CLSI ", "", amr$guideline))

count <- function(v) as.vector(table(factor(v, levels = editions)))
audit <- data.frame(
  edition = paste("CLSI", editions),
  whonet_2a_family = count(who$YEAR[who$REFERENCE_TABLE == "Table 2A" & who$ORGANISM_CODE == "EBC"]),
  whonet_2a_all = count(who$YEAR[who$REFERENCE_TABLE == "Table 2A"]),
  whonet_other = count(who$YEAR[who$REFERENCE_TABLE != "Table 2A"]),
  amr_2a_family = count(amr$year[amr$ref_tbl == "Table 2A" & as.character(amr$mo) == "B_[ORD]_ENTRBCTR"]),
  amr_2a_all = count(amr$year[amr$ref_tbl == "Table 2A"]),
  amr_other = count(amr$year[amr$ref_tbl != "Table 2A"]))

# One E. coli MIC of 1 mg/L, every drug with an Enterobacterales MIC row in CLSI 2020, read with
# CLSI 2019 and CLSI 2020: which reading comes from the generic Gram-negative row (B_GRAMN)?
drugs <- sort(unique(as.character(amr$ab[amr$guideline == "CLSI 2020" & amr$ref_tbl == "Table 2A" &
                                         amr$method == "MIC" & as.character(amr$mo) == "B_[ORD]_ENTRBCTR"])))
source_row <- function(guideline) vapply(drugs, function(ab) {
  AMR::sir_interpretation_history(clean = TRUE)
  AMR::as.sir(AMR::as.mic(1), mo = "B_ESCHR_COLI", ab = ab, guideline = guideline, info = FALSE)
  h <- AMR::sir_interpretation_history()
  # no breakpoint at all: outcome NA, and the row names the input organism
  if (nrow(h) == 0 || is.na(h$outcome[nrow(h)])) NA_character_ else as.character(h$mo[nrow(h)])
}, character(1))
rows_2019 <- source_row("CLSI 2019")
rows_2020 <- source_row("CLSI 2020")
example <- data.frame(ab = drugs, name = AMR::ab_name(drugs, info = FALSE), mo_2019 = rows_2019, mo_2020 = rows_2020)

dir.create("results", showWarnings = FALSE)
write.csv(audit, "results/whonet_audit.csv", row.names = FALSE)
write.csv(example, "results/whonet_example.csv", row.names = FALSE, na = "")
writeLines(c(
  "\\begin{table}[h]", "\\centering", "\\scriptsize", "\\begin{tabular}{lrrrrrr}", "\\toprule",
  " & \\multicolumn{3}{c}{WHONET \\texttt{Breakpoints.txt}} & \\multicolumn{3}{c}{AMR 3.1.0} \\\\",
  "Edition & 2A family & 2A all & Other tables & 2A family & 2A all & Other tables \\\\", "\\midrule",
  with(audit, sprintf("%s & %d & %d & %d & %d & %d & %d \\\\", edition, whonet_2a_family, whonet_2a_all, whonet_other,
                      amr_2a_family, amr_2a_all, amr_other)),
  "\\bottomrule", "\\end{tabular}",
  sprintf(paste("\\caption{Human CLSI breakpoint rows per edition in WHONET's \\texttt{Breakpoints.txt} (AMRIE commit",
                "69de201, SHA-256 %s\\ldots) and in the AMR package 3.1.0. 2A family: Table 2A rows coded for the",
                "whole Enterobacterales family (WHONET code EBC, AMR \\texttt{B\\_[ORD]\\_ENTRBCTR}); 2A all: every",
                "Table 2A row, including those for \\textit{Salmonella} and single species. For an \\textit{E.~coli}",
                "MIC of 1~mg/L and the %d drugs with an Enterobacterales MIC row in CLSI 2020, \\texttt{as.sir()}",
                "used the generic Gram-negative row (\\texttt{B\\_GRAMN}, Table 2B-5) for %d under CLSI 2019 and %d",
                "under CLSI 2020; %d had no CLSI 2019 breakpoint at all (\\texttt{results/whonet\\_example.csv}).}"),
          substr(sha, 1, 12), length(drugs), sum(rows_2019 %in% "B_GRAMN"), sum(rows_2020 %in% "B_GRAMN"),
          sum(is.na(rows_2019))),
  "\\end{table}", ""), "results/supplement_whonet.tex")
message("E. coli MIC 1, ", length(drugs), " drugs: B_GRAMN under CLSI 2019 ", sum(rows_2019 %in% "B_GRAMN"),
        ", under CLSI 2020 ", sum(rows_2020 %in% "B_GRAMN"), "; no reading ", sum(is.na(rows_2019)), " and ", sum(is.na(rows_2020)))
