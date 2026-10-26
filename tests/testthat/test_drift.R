# Every test builds MICs by hand, so the true real change and artefact are known exactly.
# Most use CLSI's 2023 revision of the Enterobacterales gentamicin MIC breakpoints, from
# S <= 4 / R >= 16 (2022 edition) to S <= 2 / R >= 8. Both editions have their own
# Enterobacterales rows in the AMR package; editions before 2020 often do not (see read_sir()).

isolates <- function(year, mics, mo = "B_ESCHR_COLI", ab = "GEN") {
  data.frame(year = year, mo = mo, ab = ab, mic = mics)
}
# 100 isolates: 10% R under CLSI 2022 (only the 16s), 40% R under CLSI 2023 (8 and 16)
fixed_mics <- rep(c("1", "2", "8", "16"), c(40, 20, 30, 10))

test_that("an unchanged MIC distribution under a revised breakpoint is all artefact", {
  d <- drift_decompose(rbind(isolates(2022, fixed_mics), isolates(2023, fixed_mics)))
  expect_equal(c(d$guideline_from, d$guideline_to), c("CLSI 2022", "CLSI 2023"))
  expect_equal(c(d$p_y1g1, d$p_y2g2), c(0.10, 0.40))
  expect_equal(d$observed, 0.30)
  expect_equal(d$real, 0)
  expect_equal(d$artefact, 0.30)
})

test_that("a shifted MIC distribution under fixed breakpoints is all real change", {
  # Ciprofloxacin: S <= 0.25 / R >= 1 in both editions
  later <- rep(c("0.25", "0.5", "1", "8"), c(50, 10, 20, 20))  # 40% R
  d <- drift_decompose(rbind(isolates(2021, rep(c("0.25", "1"), c(90, 10)), ab = "CIP"),
                             isolates(2023, later, ab = "CIP")))
  expect_equal(d$observed, 0.30)
  expect_equal(d$real, 0.30)
  expect_equal(d$artefact, 0)
})

test_that("real change and artefact together add up to the observed change", {
  set.seed(1)
  dilutions <- c("0.5", "1", "2", "4", "8", "16", "32")
  early <- sample(dilutions, 300, replace = TRUE, prob = 7:1)
  late <- sample(dilutions, 300, replace = TRUE, prob = c(4:2, 1:4))
  d <- drift_decompose(rbind(isolates(2020, early), isolates(2024, late)))
  expect_equal(d$real + d$artefact, d$observed)
  expect_equal(d$real, (d$real_at_new + d$real_at_old) / 2)
  expect_gt(d$real, 0)
  expect_gt(d$artefact, 0)
})

test_that("a capped MIC the new breakpoint cannot read is left out and counted", {
  capped <- c(rep("<=4", 20), fixed_mics[-(1:20)])  # 20 of the 1s were only recorded as <=4
  d <- drift_decompose(rbind(isolates(2022, capped), isolates(2023, fixed_mics)))
  expect_equal(d$tested_y1, 100)
  expect_equal(d$n_y1g1, 100)      # <=4 is S under S <= 4
  expect_equal(d$n_y1g2, 80)       # but S or I under S <= 2, R >= 8
  expect_equal(d$p_y1g2, 40 / 80)  # the 8s and 16s among the 80 readable isolates
})

test_that("periods pool their years and use the edition in force in their last year", {
  d <- drift_decompose(rbind(isolates(2020, fixed_mics), isolates(2021, fixed_mics),
                             isolates(2023, fixed_mics), isolates(2024, fixed_mics),
                             isolates(2025, rep("16", 100))),  # outside both periods
                       from = 2020:2021, to = 2023:2024)
  expect_equal(c(d$year_from, d$year_to), c("2020-2021", "2023-2024"))
  expect_equal(c(d$guideline_from, d$guideline_to), c("CLSI 2021", "CLSI 2024"))
  expect_equal(c(d$tested_y1, d$tested_y2), c(200, 200))
  expect_equal(c(d$real, d$artefact), c(0, 0.30))
})

test_that("pairs are kept apart", {
  d <- drift_decompose(rbind(isolates(2022, fixed_mics), isolates(2023, fixed_mics),
                             isolates(2022, fixed_mics, mo = "B_KLBSL_PNMN"),
                             isolates(2023, rep("16", 100), mo = "B_KLBSL_PNMN")))
  expect_equal(nrow(d), 2)
  expect_equal(d$observed[d$mo == "B_ESCHR_COLI"], 0.30)
  expect_equal(d$observed[d$mo == "B_KLBSL_PNMN"], 0.90)
})

test_that("verdicts separate artefact, real change and noise", {
  big <- rep(fixed_mics, 3)
  shifted <- rep(c("1", "16"), c(60, 240))  # 40% R becomes 80% R
  data <- rbind(
    isolates(2022, big), isolates(2023, big),                                 # artefact only
    isolates(2022, big, mo = "B_KLBSL_PNMN"), isolates(2023, shifted, mo = "B_KLBSL_PNMN"),
    isolates(2022, c("0.25", "0.25", "8"), mo = "B_ESCHR_COLI", ab = "LVX"),  # too few to tell
    isolates(2023, c("0.25", "8", "8"), mo = "B_ESCHR_COLI", ab = "LVX")
  )
  t <- drift_test(drift_decompose(data))
  verdict <- setNames(t$verdict, paste(t$mo, t$ab))
  expect_equal(verdict[["B_ESCHR_COLI GEN"]], "breakpoint artefact")
  expect_equal(verdict[["B_KLBSL_PNMN GEN"]], "significant change")
  expect_equal(verdict[["B_ESCHR_COLI LVX"]], "noise")
  ci <- t[t$mo == "B_ESCHR_COLI" & t$ab == "GEN", c("ci_low_y1g2", "ci_high_y1g2")]
  expect_lt(ci[[1]], 0.40)
  expect_gt(ci[[2]], 0.40)
})

test_that("a drug with no breakpoint is not interpretable rather than noise", {
  t <- drift_test(drift_decompose(rbind(isolates(2022, fixed_mics, ab = "VAN"),
                                        isolates(2023, fixed_mics, ab = "VAN"))))
  expect_equal(t$verdict, "not interpretable")
})

test_that("guideline in force follows the year, a lag, and the AMR package's coverage", {
  expect_equal(guideline_in_force(c(2015, 2020)), c("CLSI 2015", "CLSI 2020"))
  expect_equal(guideline_in_force(2020, lag = 1), "CLSI 2019")
  expect_equal(guideline_in_force(2020, "EUCAST"), "EUCAST 2020")
  expect_true(is.na(guideline_in_force(2005)))
  expect_equal(guideline_in_force(2099), guideline_in_force(2026))
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
  expect_equal(c(z$real, z$artefact), c(0, 0.10))
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
  expect_equal(c(e$real, e$artefact), c(0, 0.10))
  expect_equal(d$artefact[d$org == "CLSI"], 0)  # no CLSI change between these editions
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
  t <- merge(drift_test(drift_decompose(data)), drift_mic_shift(data), by = c("mo", "ab"))
  path <- withr::local_tempfile(fileext = ".md")
  drift_report(t, path)
  md <- readLines(path)
  expect_equal(sum(startsWith(md, "## ")), 2)
  expect_equal(sum(startsWith(md, "| Gentamicin | MIC CLSI | 10.0 | 40.0 | +30.0 | +0.0 | +30.0 |")), 2)
  expect_true(any(grepl("MIC shift", md)))
})

test_that("an adoption table replaces the year rule", {
  data <- rbind(isolates(2022, fixed_mics), isolates(2024, fixed_mics))
  expect_equal(drift_decompose(data)$artefact, 0.30)  # the rule: CLSI 2022, then CLSI 2024
  # This lab only moved past CLSI 2022 in 2025, so both periods were read the same way
  late <- data.frame(year = c(2015, 2025), guideline = c("CLSI 2022", "CLSI 2025"))
  d <- drift_decompose(data, adoption = late)
  expect_equal(c(d$guideline_from, d$guideline_to), c("CLSI 2022", "CLSI 2022"))
  expect_equal(c(d$observed, d$artefact), c(0, 0))
})

test_that("labs on different editions are pooled, each read with its own", {
  data <- rbind(cbind(isolates(2022, fixed_mics), lab = "A"), cbind(isolates(2023, fixed_mics), lab = "A"),
                cbind(isolates(2022, fixed_mics), lab = "B"), cbind(isolates(2023, fixed_mics), lab = "B"))
  adoption <- data.frame(lab = c("A", "A", "B"), year = c(2020, 2023, 2020),
                         guideline = c("CLSI 2020", "CLSI 2023", "CLSI 2020"))
  d <- drift_decompose(data, adoption = adoption)
  expect_equal(d$guideline_to, "CLSI 2020 + CLSI 2023")
  expect_equal(c(d$p_y1g1, d$p_y2g2), c(0.10, 0.25))  # only lab A's half moved to 40% R
  expect_equal(c(d$real, d$artefact), c(0, 0.15))
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
  t <- drift_test(drift_decompose(data))
  expect_equal(c(t$p_rep_y1, t$p_rep_y2), c(0.40, 0.40))
  expect_equal(c(t$observed, t$real, t$artefact), c(0, 0, 0))
  expect_equal(t$verdict, "noise")
  # Now the lab relaxes its rule instead: reported %R falls 30 points with the MICs unchanged
  data <- rbind(rule(2022, fixed_mics, c("8", "16")), rule(2023, fixed_mics, "16"))
  t <- drift_test(drift_decompose(data))
  expect_equal(c(t$observed, t$real, t$artefact), c(-0.30, 0, -0.30))
  expect_equal(t$verdict, "breakpoint artefact")
  path <- withr::local_tempfile(fileext = ".md")
  drift_report(t, path)
  expect_true(any(startsWith(readLines(path), "| Gentamicin | MIC CLSI | 40.0 | 10.0 | -30.0 |")))
})
