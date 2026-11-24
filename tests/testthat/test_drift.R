# Every test builds MICs by hand, so the true change at fixed criteria and reclassification are known exactly.
# Most use CLSI's 2023 revision of the Enterobacterales gentamicin MIC breakpoints, from
# S <= 4 / R >= 16 (2022 edition) to S <= 2 / R >= 8. Both editions have their own
# Enterobacterales rows in the AMR package; editions before 2020 often do not (see read_sir()).

isolates <- function(year, mics, mo = "B_ESCHR_COLI", ab = "GEN") {
  data.frame(year = year, mo = mo, ab = ab, mic = mics)
}
# 100 isolates: 10% R under CLSI 2022 (only the 16s), 40% R under CLSI 2023 (8 and 16)
fixed_mics <- rep(c("1", "2", "8", "16"), c(40, 20, 30, 10))

test_that("an unchanged MIC distribution under a revised breakpoint is all reclass", {
  d <- drift_decompose(rbind(isolates(2022, fixed_mics), isolates(2023, fixed_mics)))
  expect_equal(c(d$guideline_from, d$guideline_to), c("CLSI 2022", "CLSI 2023"))
  expect_equal(c(d$p_y1g1, d$p_y2g2), c(0.10, 0.40))
  expect_equal(d$observed, 0.30)
  expect_equal(d$fixed, 0)
  expect_equal(d$reclass, 0.30)
})

test_that("a shifted MIC distribution under fixed breakpoints is all fixed change", {
  # Ciprofloxacin: S <= 0.25 / R >= 1 in both editions
  later <- rep(c("0.25", "0.5", "1", "8"), c(50, 10, 20, 20))  # 40% R
  d <- drift_decompose(rbind(isolates(2021, rep(c("0.25", "1"), c(90, 10)), ab = "CIP"),
                             isolates(2023, later, ab = "CIP")))
  expect_equal(d$observed, 0.30)
  expect_equal(d$fixed, 0.30)
  expect_equal(d$reclass, 0)
})

test_that("fixed change and reclass together add up to the observed change", {
  set.seed(1)
  dilutions <- c("0.5", "1", "2", "4", "8", "16", "32")
  early <- sample(dilutions, 300, replace = TRUE, prob = 7:1)
  late <- sample(dilutions, 300, replace = TRUE, prob = c(4:2, 1:4))
  d <- drift_decompose(rbind(isolates(2020, early), isolates(2024, late)))
  expect_equal(d$fixed + d$reclass, d$observed)
  expect_equal(d$fixed, (d$fixed_at_new + d$fixed_at_old) / 2)
  expect_gt(d$fixed, 0)
  expect_gt(d$reclass, 0)
})

test_that("a capped MIC counts wherever it falls on one side of the outcome's boundary", {
  capped <- c(rep("<=4", 20), fixed_mics[-(1:20)])  # 20 of the 1s were only recorded as <=4
  d <- drift_decompose(rbind(isolates(2022, capped), isolates(2023, fixed_mics)))
  expect_equal(d$tested_y1, 100)
  expect_equal(c(d$n_y1g1, d$n_y1g2), c(100, 100))  # <=4 is not R under R >= 16 or R >= 8
  expect_equal(d$p_y1g2, 0.40)
  expect_equal(c(d$fixed, d$reclass), c(0, 0.30))
  # For %NS, <=4 is S or I under CLSI 2023 (S <= 2, I = 4), so it is left out there only
  ns <- drift_decompose(rbind(isolates(2022, capped), isolates(2023, fixed_mics)), outcome = "NS")
  expect_equal(c(ns$n_y1g1, ns$n_y1g2), c(100, 80))
  expect_equal(ns$p_y1g2, 40 / 80)
  # Mostly capped panels give no false reclassification (80 at <=4 and 20 at 16, both years)
  short <- rep(c("<=4", "16"), c(80, 20))
  d <- drift_decompose(rbind(isolates(2022, short), isolates(2023, short)))
  expect_equal(c(d$p_y1g1, d$p_y1g2, d$reclass), c(0.20, 0.20, 0))
  # A value capped above: >4 is R under R >= 8 but may be I or R under R >= 16
  high <- rep(c("1", ">4"), c(50, 50))
  d <- drift_decompose(rbind(isolates(2022, high), isolates(2023, high)))
  expect_equal(c(d$n_y1g1, d$n_y1g2), c(50, 100))
})

test_that("capped values straddling one guideline's boundary only make a pair not interpretable", {
  # <=8 is not R under R >= 16 but straddles R >= 8, so the 2023 share rests on the 16s alone.
  # Reclassification then lies anywhere from 0 to +80, and no verdict may be given.
  short <- rep(c("<=8", "16"), c(80, 20))
  data <- rbind(isolates(2022, short), isolates(2023, short))
  data$isolate <- seq_len(nrow(data))
  d <- drift_decompose(data)
  expect_equal(c(d$straddled_y1, d$straddled_y2), c(80, 80))
  expect_equal(c(d$reclass_min, d$reclass_max), c(0, 0.80))
  set.seed(1)
  t <- drift_test(drift_bootstrap(data, B = 99))
  expect_equal(t$verdict, "not interpretable")
  expect_true(is.na(t$larger) && is.na(t$p_reclass))
  expect_match(t$verdict_note, "capped values")
  # Same with a value capped above (>4 straddles R >= 16 only): range 0 to +50
  high <- rep(c("1", ">4"), c(50, 50))
  d <- drift_decompose(rbind(isolates(2022, high), isolates(2023, high)))
  expect_equal(c(d$reclass_min, d$reclass_max), c(0, 0.50))
  # Without straddled values the range is the estimate itself
  d <- drift_decompose(rbind(isolates(2022, fixed_mics), isolates(2023, fixed_mics)))
  expect_equal(c(d$straddled_y1, d$reclass_min, d$reclass_max, d$fixed_min, d$fixed_max), c(0, 0.30, 0.30, 0, 0))
  # A single straddled value among many moves a component less than one point and is tolerated
  few <- c(rep(fixed_mics, 3)[-1], "<=8")
  data <- rbind(isolates(2022, few), isolates(2023, rep(fixed_mics, 3)))
  data$isolate <- seq_len(nrow(data))
  set.seed(1)
  t <- drift_test(drift_bootstrap(data, B = 99))
  expect_lt(t$reclass_max - t$reclass_min, 0.01)
  expect_equal(t$verdict, "reclassification")
})

test_that("periods pool their years and use the edition in force in their last year", {
  d <- drift_decompose(rbind(isolates(2020, fixed_mics), isolates(2021, fixed_mics),
                             isolates(2023, fixed_mics), isolates(2024, fixed_mics),
                             isolates(2025, rep("16", 100))),  # outside both periods
                       from = 2020:2021, to = 2023:2024)
  expect_equal(c(d$year_from, d$year_to), c("2020-2021", "2023-2024"))
  expect_equal(c(d$guideline_from, d$guideline_to), c("CLSI 2021", "CLSI 2024"))
  expect_equal(c(d$tested_y1, d$tested_y2), c(200, 200))
  expect_equal(c(d$fixed, d$reclass), c(0, 0.30))
})

test_that("pairs are kept apart", {
  d <- drift_decompose(rbind(isolates(2022, fixed_mics), isolates(2023, fixed_mics),
                             isolates(2022, fixed_mics, mo = "B_KLBSL_PNMN"),
                             isolates(2023, rep("16", 100), mo = "B_KLBSL_PNMN")))
  expect_equal(nrow(d), 2)
  expect_equal(d$observed[d$mo == "B_ESCHR_COLI"], 0.30)
  expect_equal(d$observed[d$mo == "B_KLBSL_PNMN"], 0.90)
})

test_that("verdicts separate reclassification, change at fixed criteria, both, and no change", {
  big <- rep(fixed_mics, 3)
  shifted <- rep(c("1", "16"), c(60, 240))  # 40% R becomes 80% R under either edition
  data <- rbind(
    isolates(2022, big), isolates(2023, big),                                 # guideline only
    isolates(2022, big, mo = "B_KLBSL_PNMN"), isolates(2023, shifted, mo = "B_KLBSL_PNMN"),
    isolates(2022, c("0.25", "0.25", "8"), mo = "B_ESCHR_COLI", ab = "LVX"),  # too few to tell
    isolates(2023, c("0.25", "8", "8"), mo = "B_ESCHR_COLI", ab = "LVX")
  )
  data$isolate <- seq_len(nrow(data))
  set.seed(1)
  t <- drift_test(drift_bootstrap(data, B = 499))
  verdict <- setNames(t$verdict, paste(t$mo, t$ab))
  expect_equal(verdict[["B_ESCHR_COLI GEN"]], "reclassification")
  expect_equal(verdict[["B_KLBSL_PNMN GEN"]], "both")  # +40 at fixed criteria, and 8s read R in 2022 under CLSI 2023
  expect_equal(verdict[["B_ESCHR_COLI LVX"]], "no change detected")
  expect_equal(t$larger[t$mo == "B_KLBSL_PNMN"], "fixed criteria")
  ci <- t[t$mo == "B_ESCHR_COLI" & t$ab == "GEN", c("ci_low_y1g2", "ci_high_y1g2")]
  expect_lt(ci[[1]], 0.40)
  expect_gt(ci[[2]], 0.40)
})

test_that("a drug with no breakpoint is not interpretable and stays out of the corrections", {
  data <- rbind(isolates(2022, fixed_mics, ab = "VAN"), isolates(2023, fixed_mics, ab = "VAN"),
                isolates(2022, fixed_mics), isolates(2023, fixed_mics))
  t <- drift_test(drift_bootstrap(data, B = 99))
  expect_equal(t$verdict[t$ab == "VAN"], "not interpretable")
  expect_true(all(is.na(unlist(t[t$ab == "VAN", c("p_observed", "p_fixed", "p_reclass", "p_larger")]))))
  expect_error(drift_test(drift_decompose(data)), "drift_bootstrap")
})

test_that("guideline in force follows the year, a lag, and the AMR package's coverage", {
  expect_equal(guideline_in_force(c(2015, 2020)), c("CLSI 2015", "CLSI 2020"))
  expect_equal(guideline_in_force(2020, lag = 1), "CLSI 2019")
  expect_equal(guideline_in_force(2020, "EUCAST"), "EUCAST 2020")
  expect_true(is.na(guideline_in_force(2005)))
  expect_equal(guideline_in_force(2099), guideline_in_force(2026))
  expect_equal(guideline_in_force(c(2005, 2024)), c(NA, "CLSI 2024"))
  expect_error(drift_decompose(rbind(isolates(2005, "1"), isolates(2020, "1"))), "2005")
})

test_that("zone diameters are decomposed the same way, apart from MICs", {
  # CLSI 2023 raised the E. coli gentamicin disk breakpoints: 14 mm went from I to R
  zones <- rep(c(20, 16, 14, 10), c(70, 10, 10, 10))
  disk <- function(year) data.frame(year = year, mo = "B_ESCHR_COLI", ab = "GEN", mic = NA, disk = zones)
  mic <- function(year) cbind(isolates(year, fixed_mics, ab = "GEN"), disk = NA)
  data <- rbind(disk(2022), disk(2023), mic(2022), mic(2023))
  d <- drift_decompose(data)
  expect_setequal(d$method, c("disk", "MIC"))
  z <- d[d$method == "disk", ]
  expect_equal(c(z$p_y1g1, z$p_y2g2), c(0.10, 0.20))
  expect_equal(c(z$fixed, z$reclass), c(0, 0.10))
  expect_error(drift_decompose(data.frame(year = 2020:2021, mo = "E. coli", ab = "GEN", mic = "1", disk = 20)),
               "both")
})

test_that("EUCAST rows are read with EUCAST editions, CLSI rows with CLSI", {
  # EUCAST 2021 cut the E. coli piperacillin/tazobactam breakpoint: an MIC of 16 went from I to R
  eucast <- rep(c("2", "8", "16", "32"), c(60, 20, 10, 10))
  data <- rbind(cbind(isolates(2020, eucast, ab = "TZP"), org = "EUCAST"),
                cbind(isolates(2021, eucast, ab = "TZP"), org = "EUCAST"),
                cbind(isolates(2020, fixed_mics), org = "CLSI"), cbind(isolates(2021, fixed_mics), org = "CLSI"))
  d <- drift_decompose(data)
  e <- d[d$org == "EUCAST", ]
  expect_equal(c(e$guideline_from, e$guideline_to), c("EUCAST 2020", "EUCAST 2021"))
  expect_equal(c(e$fixed, e$reclass), c(0, 0.10))
  expect_equal(d$reclass[d$org == "CLSI"], 0)  # no CLSI change between these editions
})

# Draw true log2 MICs, then store them the way a dilution panel does: rounded up to the next
# dilution, capped at the panel ends.
panel <- function(z) {
  dilutions <- c("0.008", "0.016", "0.032", "0.064", "0.125", "0.25", "0.5", "1", "2", "4", "8")  # 2^-7 to 2^3
  k <- pmin(pmax(ceiling(z), -7), 4) + 8
  ifelse(k == 1, "<=0.008", ifelse(k == 12, ">8", dilutions[pmin(k, 11)]))
}

test_that("the censored model recovers a known one-dilution shift", {
  set.seed(42)
  data <- rbind(isolates(2015, panel(rnorm(500, -3, 1.5))), isolates(2024, panel(rnorm(500, -2, 1.5))))
  s <- drift_mic_shift(data, from = 2015, to = 2024)
  expect_lt(abs(s$shift - 1), 0.25)
  expect_true(s$shift_low < 1 && s$shift_high > 1)
  expect_lt(s$p_shift, 0.001)
  expect_equal(c(s$n_shift_from, s$n_shift_to), c(500, 500))
})

test_that("MIC creep below the breakpoint shows as a shift while %R stays at zero", {
  set.seed(7)
  data <- rbind(isolates(2021, panel(rnorm(300, -5, 0.5))), isolates(2024, panel(rnorm(300, -4, 0.5))))
  d <- drift_decompose(data)
  s <- drift_mic_shift(data)
  expect_equal(c(d$p_y1g2, d$p_y2g2), c(0, 0))
  expect_lt(abs(s$shift - 1), 0.2)
  expect_lt(s$p_shift, 0.001)
})

test_that("the censored model gives NA, not a number, when it cannot estimate", {
  all_capped <- rbind(isolates(2015, rep("<=0.25", 20)), isolates(2024, rep("<=0.25", 20)))
  expect_true(is.na(drift_mic_shift(all_capped)$shift))
  one_side <- rbind(isolates(2015, rep("<=1", 40)), isolates(2024, rep(c("2", "8", ">16"), 10)))
  expect_true(is.na(drift_mic_shift(one_side)$shift))
  few <- rbind(isolates(2015, c("1", "2")), isolates(2024, c("2", "4")))
  expect_true(is.na(drift_mic_shift(few)$shift))
})

test_that("stored MICs become the right log2 intervals", {
  b <- mic_interval(c("8", "<=0.25", ">=16", ">8", "<0.5"))
  expect_equal(b$lo, c(2, NA, 3, 3, NA))
  expect_equal(b$hi, c(3, -2, NA, NA, -2))
})

test_that("the report lists every pair under its organism", {
  data <- rbind(isolates(2022, fixed_mics), isolates(2023, fixed_mics),
                isolates(2022, fixed_mics, mo = "B_KLBSL_PNMN"), isolates(2023, fixed_mics, mo = "B_KLBSL_PNMN"))
  t <- merge(drift_test(drift_bootstrap(data, B = 99)), drift_mic_shift(data), by = c("mo", "ab"))
  path <- withr::local_tempfile(fileext = ".md")
  drift_report(t, path)
  md <- readLines(path)
  expect_equal(sum(startsWith(md, "## ")), 2)
  expect_equal(sum(startsWith(md, "| Gentamicin | MIC CLSI | 10.0 | 40.0 | +30.0 | +0.0 | +30.0 |")), 2)
  expect_true(any(grepl("MIC shift", md)))
  drift_report(t, path, outcome = "NS")
  md <- readLines(path)
  expect_true(any(grepl("%NS 2022", md, fixed = TRUE)))
  expect_false(any(grepl("%R", md, fixed = TRUE)))
})

test_that("an adoption table replaces the year rule", {
  data <- rbind(isolates(2022, fixed_mics), isolates(2024, fixed_mics))
  expect_equal(drift_decompose(data)$reclass, 0.30)  # the rule: CLSI 2022, then CLSI 2024
  # This lab only moved past CLSI 2022 in 2025, so both periods were read the same way
  late <- data.frame(year = c(2015, 2025), guideline = c("CLSI 2022", "CLSI 2025"))
  d <- drift_decompose(data, adoption = late)
  expect_equal(c(d$guideline_from, d$guideline_to), c("CLSI 2022", "CLSI 2022"))
  expect_equal(c(d$observed, d$reclass), c(0, 0))
})

test_that("labs on different editions are pooled, each read with its own", {
  data <- rbind(cbind(isolates(2022, fixed_mics), lab = "A"), cbind(isolates(2023, fixed_mics), lab = "A"),
                cbind(isolates(2022, fixed_mics), lab = "B"), cbind(isolates(2023, fixed_mics), lab = "B"))
  adoption <- data.frame(lab = c("A", "A", "B"), year = c(2020, 2023, 2020),
                         guideline = c("CLSI 2020", "CLSI 2023", "CLSI 2020"))
  d <- drift_decompose(data, adoption = adoption)
  expect_equal(d$guideline_to, "CLSI 2020 + CLSI 2023")
  expect_equal(c(d$p_y1g1, d$p_y2g2), c(0.10, 0.25))  # only lab A's half moved to 40% R
  expect_equal(c(d$fixed, d$reclass), c(0, 0.15))
  expect_error(drift_decompose(rbind(data, cbind(isolates(2023, "1"), lab = "C")), adoption = adoption),
               "lab C")
  expect_error(drift_decompose(data[names(data) != "lab"], adoption = adoption), "lab")
})

test_that("a lab that switched from CLSI to EUCAST gets its own row", {
  data <- rbind(cbind(isolates(2021, fixed_mics), lab = "A"), cbind(isolates(2024, fixed_mics), lab = "A"),
                cbind(isolates(2021, fixed_mics), lab = "B"), cbind(isolates(2024, fixed_mics), lab = "B"))
  adoption <- data.frame(lab = c("A", "B", "B"), year = c(2020, 2020, 2022),
                         guideline = c("CLSI 2020", "CLSI 2020", "EUCAST 2022"))
  d <- drift_decompose(data, adoption = adoption)
  expect_equal(sort(d$org), c("CLSI", "CLSI to EUCAST"))
  expect_equal(d$guideline_to[d$org == "CLSI to EUCAST"], "EUCAST 2022")
  expect_equal(d$tested_y1, c(100, 100))
})

test_that("adopted_edition takes the latest entry on or before the year", {
  a <- data.frame(year = c(2014, 2018), guideline = c("CLSI 2014", "EUCAST 2018"))
  expect_equal(adopted_edition(a, c(2013, 2014, 2017, 2018, 2030)),
               c(NA, "CLSI 2014", "CLSI 2014", "EUCAST 2018", "EUCAST 2018"))
})

test_that("the fallback check follows the order, so genera moved out of Enterobacteriaceae in 2016 count too", {
  # Proteus (now Morganellaceae) is still Enterobacterales, which CLSI Table 2A covers; a
  # non-fermenter such as Chryseobacterium is not, and Table 2B-5 is its proper table
  data <- rbind(isolates(2015, rep("8", 50), mo = "B_PROTS_MRBL", ab = "CRO"),
                isolates(2021, rep("8", 50), mo = "B_PROTS_MRBL", ab = "CRO"),
                isolates(2015, rep("8", 50), mo = "B_CHRYS_INDL", ab = "CRO"),
                isolates(2021, rep("8", 50), mo = "B_CHRYS_INDL", ab = "CRO"))
  d <- drift_decompose(data)
  expect_equal(d$n_y1g1[d$mo == "B_PROTS_MRBL"], 0)
  expect_equal(d$n_y1g1[d$mo == "B_CHRYS_INDL"], 50)
})

test_that("Enterobacterales read with the generic Gram-negative fallback count as not readable", {
  # AMR 3.1.0 has no Enterobacterales ceftriaxone row in CLSI 2015 and falls back to S <= 8
  data <- rbind(isolates(2015, rep("8", 50), ab = "CRO"), isolates(2021, rep("8", 50), ab = "CRO"),
                isolates(2015, rep("8", 50), mo = "B_PSDMN_AERG", ab = "CAZ"),
                isolates(2021, rep("8", 50), mo = "B_PSDMN_AERG", ab = "CAZ"))
  d <- drift_decompose(data)
  cro <- d[d$ab == "CRO", ]
  expect_equal(c(cro$n_y1g1, cro$n_y1g2), c(0, 50))  # CLSI 2021 has the Enterobacterales row
  expect_equal(cro$p_y1g2, 1)                          # 8 is R under S <= 1, R >= 4
  expect_equal(d$n_y1g1[d$ab == "CAZ"], 50)            # P. aeruginosa has its own row in 2015
})

test_that("the censored model gives NA when the fit lies outside the tested range", {
  # 95% of each period at the bottom of the panel: the fitted mean falls below any tested dilution
  low <- function(year, top) isolates(year, c(rep("<=0.25", 190), rep(top, 10)))
  expect_true(is.na(drift_mic_shift(rbind(low(2015, "0.5"), low(2024, "1")))$shift))
})

test_that("reported calls replace the guideline reading as the observed change", {
  # The lab reported 8 as R in 2022 (its own stricter rule) and switched to CLSI in 2023, where
  # 8 is R as well: reported %R stays 40%, while the CLSI editions alone would show +30 points
  rule <- function(year, mics, r) cbind(isolates(year, rep(mics, 3)), reported = ifelse(rep(mics, 3) %in% r, "R", "S"))
  data <- rbind(rule(2022, fixed_mics, c("8", "16")), rule(2023, fixed_mics, c("8", "16")))
  t <- drift_test(drift_bootstrap(data, B = 199))
  expect_equal(c(t$p_rep_y1, t$p_rep_y2), c(0.40, 0.40))
  expect_equal(c(t$observed, t$fixed, t$reclass), c(0, 0, 0))
  expect_equal(t$verdict, "no change detected")
  # Now the lab relaxes its rule instead: reported %R falls 30 points with the MICs unchanged
  data <- rbind(rule(2022, fixed_mics, c("8", "16")), rule(2023, fixed_mics, "16"))
  t <- drift_test(drift_bootstrap(data, B = 199))
  expect_equal(c(t$observed, t$fixed, t$reclass), c(-0.30, 0, -0.30))
  expect_equal(t$verdict, "reclassification")
  path <- withr::local_tempfile(fileext = ".md")
  drift_report(t, path)
  expect_true(any(startsWith(readLines(path), "| Gentamicin | MIC CLSI | 40.0 | 10.0 | -30.0 |")))
})

test_that("non-susceptible share ignores a revision that only moves the I/R boundary", {
  # CLSI 2022: 8 is I and 16 is R; CLSI 2023: both are R. Either way 8 and 16 are non-susceptible
  d <- drift_decompose(rbind(isolates(2022, fixed_mics), isolates(2023, fixed_mics)), outcome = "NS")
  expect_equal(c(d$p_y1g1, d$p_y2g2), c(0.40, 0.40))
  expect_equal(c(d$observed, d$fixed, d$reclass), c(0, 0, 0))
  expect_error(drift_decompose(isolates(2022, fixed_mics), outcome = "S"))
})

test_that("bootstrap intervals of one dataset contain its estimate and true values, and repeat with the seed", {
  # One dataset only: coverage over many datasets is checked by scripts/simulation.R
  data <- rbind(isolates(2022, rep(fixed_mics, 3)), isolates(2023, rep(fixed_mics, 3)))
  data$isolate <- seq_len(nrow(data))
  set.seed(1)
  b <- drift_bootstrap(data, B = 300)
  d <- drift_decompose(data)
  expect_equal(b[c("observed", "fixed", "reclass")], d[c("observed", "fixed", "reclass")])
  expect_true(b$reclass_low < 0.30 && b$reclass_high > 0.30)
  expect_true(b$fixed_low < 0 && b$fixed_high > 0)
  set.seed(1)
  expect_equal(drift_bootstrap(data, B = 300)$reclass_low, b$reclass_low)
})

test_that("a cluster bootstrap with too few clusters per period stops or warns", {
  data <- rbind(isolates(2022, fixed_mics), isolates(2023, fixed_mics))
  data$site <- rep(c("a", "b"), each = 100)  # one cluster per period: no spread, so no test
  expect_error(drift_bootstrap(data, B = 50, cluster = "site"), "at least 5")
  data$site <- paste(data$year, rep(1:2, 50))
  expect_error(drift_bootstrap(data, B = 50, cluster = "site"), "2 clusters")
  data$site <- paste(data$year, rep(1:5, 20))
  expect_warning(drift_bootstrap(data, B = 20, cluster = "site"), "approximate")
  expect_error(drift_bootstrap(data, B = 10, cluster = "lab"))
})

test_that("the larger component is named only where the verdict shows it", {
  # Reclassification exactly 0, so |reclass| - |fixed| has the fixed component's raw p-value; the
  # separate corrections must not name fixed criteria larger while the verdict says no change detected
  x <- data.frame(r_y1g1 = c(50, 50, 50), n_y1g1 = 500, r_y2g2 = 60, n_y2g2 = 500, r_y1g2 = 50,
                  n_y1g2 = 500, observed = 0.02, fixed = 0.02, reclass = 0,
                  boot_p_fixed = c(0.03, 0.04, 0.9), boot_p_reclass = 1, boot_p_larger = c(0.03, 0.04, 0.001))
  t <- drift_test(x)
  expect_true(all(t$larger == "neither" | t$verdict %in% c("change at fixed criteria", "reclassification", "both")))
  expect_equal(t$larger[3], "neither")
  # Observed change equal to the change at fixed criteria: never "unattributed"
  x <- data.frame(r_y1g1 = 50, n_y1g1 = 1000, r_y2g2 = 75, n_y2g2 = 1000, r_y1g2 = 50, n_y1g2 = 1000,
                  observed = 0.025, fixed = 0.025, reclass = 0, boot_p_fixed = 0.06,
                  boot_p_reclass = 1, boot_p_larger = 0.06)
  expect_equal(drift_test(x)$verdict, "no change detected")
})

test_that("a lab following each year's guideline exactly has no reclassification beyond the revision", {
  # Reported calls equal CLSI 2022 in 2022 and CLSI 2023 in 2023: reclassification is the
  # revision alone, as without reported calls
  read <- function(year, edition) {
    x <- isolates(year, rep(fixed_mics, 3))
    x$method <- "MIC"
    x$disk <- NA
    cbind(x[c("year", "mo", "ab", "mic")], reported = read_sir(x, edition))
  }
  data <- rbind(read(2022, "CLSI 2022"), read(2023, "CLSI 2023"))
  set.seed(1)
  t <- drift_test(drift_bootstrap(data, B = 499))
  expect_equal(c(t$observed, t$fixed, t$reclass), c(0.30, 0, 0.30))
  plain <- drift_decompose(data[names(data) != "reported"])
  expect_equal(t$reclass, plain$reclass)
  expect_equal(t$verdict, "reclassification")
})

test_that("reclassification much larger than the change at fixed criteria is not hidden by sample size", {
  # Two-dilution revision of piperacillin/tazobactam (CLSI 2021 to 2022: R >= 128 to R >= 32):
  # 20% of isolates sit at 32-64 in both years, and 5 points more are R by any reading in 2022
  early <- rep(c("4", "32", "128"), c(750, 200, 50))
  late <- rep(c("4", "32", "128"), c(700, 200, 100))
  for (k in c(1, 5, 20)) {  # 1,000, 5,000 and 20,000 isolates per period
    data <- rbind(isolates(2021, rep(early, k), ab = "TZP"), isolates(2022, rep(late, k), ab = "TZP"))
    set.seed(1)
    t <- drift_test(drift_bootstrap(data, B = 499))
    expect_equal(c(t$fixed, t$reclass), c(0.05, 0.20))
    expect_equal(t$verdict, "both")
    expect_equal(t$larger, "reclassification")
  }
})

test_that("drift_test calls a significant observed change that neither component shows unattributed", {
  # Given p-values, not bootstrap output: this checks the verdict rule alone
  x <- data.frame(r_y1g1 = 50, n_y1g1 = 500, r_y2g2 = 100, n_y2g2 = 500, r_y1g2 = 90, n_y1g2 = 500,
                  observed = 0.10, fixed = 0.07, reclass = 0.03, boot_p_fixed = 0.2,
                  boot_p_reclass = 0.2, boot_p_larger = 0.5)
  expect_equal(drift_test(x)$verdict, "unattributed")
  x$boot_p_reclass <- 0.001
  expect_equal(drift_test(x)$verdict, "reclassification")
  x$reclass <- 0.005  # significant but under one point
  expect_equal(drift_test(x)$verdict, "unattributed")
  expect_equal(drift_test(x, min_effect = 0)$verdict, "reclassification")
})

test_that("bootstrap p-values follow the Shapley components and are 1 for a component that is always 0", {
  # Ciprofloxacin: R >= 1 in both editions; 10% R becomes 50% R
  data <- rbind(isolates(2022, rep(fixed_mics, 3)), isolates(2023, rep(fixed_mics, 3)),
                isolates(2022, rep(c("0.25", "1"), c(270, 30)), ab = "CIP"),
                isolates(2023, rep(c("0.25", "1"), c(150, 150)), ab = "CIP"))
  set.seed(1)
  b <- drift_bootstrap(data, B = 199)
  p <- setNames(b$boot_p_fixed, b$ab)
  expect_equal(p[["CIP"]], 2 / 200)  # no replicate reaches 0
  expect_gt(p[["GEN"]], 0.2)
  expect_equal(b$boot_p_reclass[b$ab == "CIP"], 1)  # ciprofloxacin did not change in 2023
  expect_equal(b$boot_p_reclass[b$ab == "GEN"], 2 / 200)
  t <- drift_test(b)
  expect_equal(t$p_fixed, p.adjust(b$boot_p_fixed, "BH"))
})

test_that("the cluster bootstrap draws G - 1 clusters and says how many there were", {
  data <- rbind(isolates(2022, rep(fixed_mics, 10)), isolates(2023, rep(fixed_mics, 10)))
  data$site <- paste(data$year, rep(1:10, each = 100))
  b <- drift_bootstrap(data, B = 20, cluster = "site")
  expect_equal(c(b$units_y1, b$units_y2), c(10, 10))
  expect_equal(c(b$reclass_low, b$reclass_high), c(0.30, 0.30))  # identical sites
  # Five sites of 20 per period, one all resistant (gentamicin, unchanged from CLSI 2023 to 2024).
  # Drawing 4 of the 5 sites, each period's %R is a multiple of 1/4; drawing 5, of 1/5
  site_mics <- rep(c("32", "1"), c(20, 80))
  data <- rbind(isolates(2023, site_mics), isolates(2024, site_mics))
  data$site <- paste(data$year, rep(rep(1:5, each = 20), 2))
  b <- suppressWarnings(drift_bootstrap(data, B = 999, cluster = "site"))
  ends <- 4 * c(b$observed_low, b$observed_high)
  expect_equal(ends, round(ends))
  expect_true(all(ends != 0))
})

test_that("with a cluster column, an isolate with no cluster keeps its drugs together", {
  data <- rbind(isolates(2022, fixed_mics), isolates(2022, fixed_mics, ab = "CIP"),
                isolates(2023, fixed_mics), isolates(2023, fixed_mics, ab = "CIP"))
  data$isolate <- rep(rep(1:100, 2), 2) + rep(c(0, 100), each = 200)
  data$site <- ifelse((data$isolate - 1) %% 100 < 50, (data$isolate - 1) %/% 10, NA)  # 5 sites, 50 without
  b <- suppressWarnings(drift_bootstrap(data, B = 20, cluster = "site"))
  expect_equal(unique(c(b$units_y1, b$units_y2)), 5 + 50)
})

test_that("periods that are not runs of years are labelled as sets", {
  d <- drift_decompose(rbind(isolates(2015, fixed_mics), isolates(2018, fixed_mics), isolates(2023, fixed_mics)),
                       from = c(2015, 2018), to = 2023)
  expect_equal(d$year_from, "2015, 2018")
})

test_that("the censored model says why it withheld a shift", {
  few <- rbind(isolates(2015, c("1", "2")), isolates(2024, c("2", "4")))
  expect_equal(drift_mic_shift(few)$note, "too few values")
  all_capped <- rbind(isolates(2015, rep("<=0.25", 20)), isolates(2024, rep("<=0.25", 20)))
  expect_equal(drift_mic_shift(all_capped)$note, "too few uncapped values")
})

test_that("optional columns are matched by their exact name only", {
  data <- rbind(isolates(2022, fixed_mics), isolates(2023, fixed_mics))
  data$isolate_source <- "chicken"
  data$reported_date <- "2024-01-01"
  data$organism <- "E. coli"
  data$laboratory <- "A"
  b <- drift_bootstrap(data, B = 20)
  expect_equal(c(b$units_y1, b$units_y2), c(100, 100))
  expect_equal(b$observed, 0.30)
  expect_null(b$n_rep_y1)
  adoption <- data.frame(year = c(2022, 2023), guideline = c("CLSI 2022", "CLSI 2023"), label = "x")
  expect_equal(drift_decompose(data, adoption = adoption)$reclass, 0.30)
})

test_that("an isolate with no ID is its own unit", {
  data <- rbind(isolates(2022, fixed_mics), isolates(2023, fixed_mics))
  data$isolate <- c(rep(NA, 100), 1:100)
  b <- drift_bootstrap(data, B = 20)
  expect_equal(c(b$units_y1, b$units_y2), c(100, 100))
})

test_that("a capped value off the dilution grid moves to the next grid value", {
  e <- mic_ends(c(">1.5", "<1.5", ">0.12", "<0.5", "<=4", ">8"))
  expect_equal(e$lo, c("2", "0.0001", "0.25", "0.0001", "0.0001", "16"))
  expect_equal(e$hi, c("4096", "1", "4096", "0.25", "4", "4096"))
})

# CLSI ciprofloxacin, Enterobacteriaceae, M100 28th edition (2018): S <= 1, R >= 4. The AMR
# package has no Enterobacterales row for it, so these readings are otherwise dropped.
cip_2018 <- data.frame(guideline = "CLSI 2018", mo = "Enterobacterales", ab = "CIP",
                       breakpoint_S = 1, breakpoint_R = 4)

test_that("extra breakpoint rows are read with CLSI's R >= rule", {
  x <- data.frame(mo = "B_ESCHR_COLI", ab = "CIP", method = "MIC", mic = c("0.5", "1", "2", "4", "8"), disk = NA)
  expect_equal(as.vector(read_sir(x, "CLSI 2018", breakpoint_rows(cip_2018))), c("S", "S", "I", "R", "R"))
  expect_true(all(is.na(read_sir(x, "CLSI 2018"))))
})

test_that("extra rows leave every other row of the AMR package read as before", {
  b <- AMR::clinical_breakpoints
  rows <- unique(as.data.frame(b[b$guideline %in% c("CLSI 2018", "CLSI 2023") & b$method == "MIC" &
                                   b$type == "human" & !b$ab %in% "CIP", c("guideline", "mo", "ab")]))
  rows <- rows[seq(1, nrow(rows), by = 10), ]  # every 10th row; all 1,578 matched on 2026-11-17
  extra <- breakpoint_rows(rbind(cip_2018, transform(cip_2018, guideline = "CLSI 2023")))
  mics <- c("<=0.25", "0.5", "1", "2", "4", "8", "16", "32", ">64")
  for (g in unique(rows$guideline)) {
    r <- rows[rows$guideline == g, ]
    x <- data.frame(mo = rep(as.character(r$mo), each = length(mics)), ab = rep(as.character(r$ab), each = length(mics)),
                    method = "MIC", mic = mics, disk = NA)
    # AMR warns that it moves some MICs to its log2 levels; it does so in both readings
    expect_identical(suppressWarnings(as.vector(read_sir(x, g, extra))), suppressWarnings(as.vector(read_sir(x, g))))
  }
})

test_that("breakpoints reach the decomposition and leave other editions alone", {
  data <- isolates(rep(c(2018, 2023), each = 4), rep(c("0.5", "1", "2", "4"), 2), ab = "CIP")
  d <- drift_decompose(data, breakpoints = cip_2018)
  # 2018: only 4 is R (25%); 2023 (S <= 0.25, R >= 1): 1, 2 and 4 are R (75%)
  expect_equal(c(d$p_y1g1, d$p_y2g2, d$reclass), c(0.25, 0.75, 0.5))
  expect_equal(drift_decompose(data[data$year == 2023, ], from = 2022, to = 2023)$p_y2g2,
               drift_decompose(data[data$year == 2023, ], from = 2022, to = 2023, breakpoints = cip_2018)$p_y2g2)
})
