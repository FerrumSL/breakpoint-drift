# Run from the project root: Rscript -e 'testthat::test_dir("tests/testthat")'
#
# Every test builds MICs by hand, so the true real change and artefact are known exactly.
# In the AMR package's tables, CLSI moved the E. coli ciprofloxacin breakpoints from S <= 1 /
# R >= 4 (2019 edition) to S <= 0.25 / R >= 1 (2020 edition), and left them there since.
root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
source(file.path(root, "R/drift.R"), local = TRUE)

isolates <- function(year, mics, mo = "B_ESCHR_COLI", ab = "CIP") {
  data.frame(year = year, mo = mo, ab = ab, mic = mics)
}
# 100 isolates: 10% R under CLSI 2019 (only the 8s), 40% R under CLSI 2020 (1, 2 and 8)
fixed_mics <- rep(c("0.25", "0.5", "1", "2", "8"), c(40, 20, 20, 10, 10))

test_that("an unchanged MIC distribution under a revised breakpoint is all artefact", {
  d <- drift_decompose(rbind(isolates(2019, fixed_mics), isolates(2020, fixed_mics)))
  expect_equal(c(d$guideline_from, d$guideline_to), c("CLSI 2019", "CLSI 2020"))
  expect_equal(c(d$p_y1g1, d$p_y2g2), c(0.10, 0.40))
  expect_equal(d$observed, 0.30)
  expect_equal(d$real, 0)
  expect_equal(d$artefact, 0.30)
})

test_that("a shifted MIC distribution under fixed breakpoints is all real change", {
  later <- rep(c("0.25", "0.5", "1", "8"), c(50, 10, 20, 20))  # 40% R under CLSI 2020
  d <- drift_decompose(rbind(isolates(2021, rep(c("0.25", "1"), c(90, 10))), isolates(2023, later)))
  expect_equal(d$observed, 0.30)
  expect_equal(d$real, 0.30)
  expect_equal(d$artefact, 0)
})

test_that("real change and artefact together add up to the observed change", {
  set.seed(1)
  dilutions <- c("0.125", "0.25", "0.5", "1", "2", "4", "8")
  early <- sample(dilutions, 300, replace = TRUE, prob = 7:1)
  late <- sample(dilutions, 300, replace = TRUE, prob = c(4:2, 1:4))
  d <- drift_decompose(rbind(isolates(2016, early), isolates(2024, late)))
  expect_equal(d$real + d$artefact, d$observed)
  expect_equal(d$real, (d$real_at_new + d$real_at_old) / 2)
  expect_gt(d$real, 0)
  expect_gt(d$artefact, 0)
})

test_that("a capped MIC the new breakpoint cannot read is left out and counted", {
  capped <- c(rep("<=1", 20), fixed_mics[-(1:20)])  # 20 of the 0.25s were only recorded as <=1
  d <- drift_decompose(rbind(isolates(2019, capped), isolates(2020, fixed_mics)))
  expect_equal(d$tested_y1, 100)
  expect_equal(d$n_y1g1, 100)      # <=1 is S under S <= 1
  expect_equal(d$n_y1g2, 80)       # but neither S nor R under S <= 0.25, R >= 1
  expect_equal(d$p_y1g2, 40 / 80)  # the 1, 2 and 8 among the 80 readable isolates
})

test_that("periods pool their years and use the edition in force in their last year", {
  d <- drift_decompose(rbind(isolates(2017, fixed_mics), isolates(2019, fixed_mics),
                             isolates(2020, fixed_mics), isolates(2024, fixed_mics),
                             isolates(2025, rep("8", 100))),  # outside both periods
                       from = 2017:2019, to = 2020:2024)
  expect_equal(c(d$year_from, d$year_to), c("2017-2019", "2020-2024"))
  expect_equal(c(d$guideline_from, d$guideline_to), c("CLSI 2019", "CLSI 2024"))
  expect_equal(c(d$tested_y1, d$tested_y2), c(200, 200))
  expect_equal(c(d$real, d$artefact), c(0, 0.30))
})

test_that("pairs are kept apart", {
  d <- drift_decompose(rbind(isolates(2019, fixed_mics), isolates(2020, fixed_mics),
                             isolates(2019, fixed_mics, mo = "B_KLBSL_PNMN"),
                             isolates(2020, rep("8", 100), mo = "B_KLBSL_PNMN")))
  expect_equal(nrow(d), 2)
  expect_equal(d$observed[d$mo == "B_ESCHR_COLI"], 0.30)
  expect_equal(d$observed[d$mo == "B_KLBSL_PNMN"], 0.90)
})

test_that("verdicts separate artefact, real change and noise", {
  big <- rep(fixed_mics, 3)
  shifted <- rep(c("0.25", "8"), c(60, 240))  # 40% R becomes 80% R
  data <- rbind(
    isolates(2019, big), isolates(2020, big),                                 # artefact only
    isolates(2019, big, mo = "B_KLBSL_PNMN"), isolates(2020, shifted, mo = "B_KLBSL_PNMN"),
    isolates(2019, c("0.25", "0.25", "8"), mo = "B_ESCHR_COLI", ab = "LVX"),  # too few to tell
    isolates(2020, c("0.25", "8", "8"), mo = "B_ESCHR_COLI", ab = "LVX")
  )
  t <- drift_test(drift_decompose(data))
  verdict <- setNames(t$verdict, paste(t$mo, t$ab))
  expect_equal(verdict[["B_ESCHR_COLI CIP"]], "breakpoint artefact")
  expect_equal(verdict[["B_KLBSL_PNMN CIP"]], "significant change")
  expect_equal(verdict[["B_ESCHR_COLI LVX"]], "noise")
  ci <- t[t$mo == "B_ESCHR_COLI" & t$ab == "CIP", c("ci_low_y1g2", "ci_high_y1g2")]
  expect_lt(ci[[1]], 0.40)
  expect_gt(ci[[2]], 0.40)
})

test_that("a drug with no breakpoint is not interpretable rather than noise", {
  t <- drift_test(drift_decompose(rbind(isolates(2019, fixed_mics, ab = "VAN"),
                                        isolates(2020, fixed_mics, ab = "VAN"))))
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
