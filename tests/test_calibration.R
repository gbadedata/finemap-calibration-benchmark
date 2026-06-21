# test_calibration.R
# Tests for the calibration/recovery logic. These mirror the Python prototype
# that was verified before porting, and they do NOT require susieR, so they run
# fast and prove the measurement arithmetic independently of the fine-mapper.
#
# Run with:  Rscript tests/test_calibration.R

source("R/calibration.R")

ok <- function(cond, msg) {
  if (!isTRUE(cond)) stop(paste("FAIL:", msg))
  cat("  pass:", msg, "\n")
}

cat("cs_contains_causal\n")
ok(cs_contains_causal(list(c(5L, 6L, 7L)), c(6L)), "causal inside a set is found")
ok(!cs_contains_causal(list(c(1L, 2L)), c(9L)), "causal absent is not found")
ok(!cs_contains_causal(list(), c(3L)), "no sets means not found")

cat("mean_cs_size\n")
ok(mean_cs_size(list(c(1L,2L,3L), c(4L,5L))) == 2.5, "mean set size is correct")
ok(is.na(mean_cs_size(list())), "empty sets give NA size")

cat("empirical_coverage\n")
ok(abs(empirical_coverage(c(TRUE, TRUE, TRUE, FALSE)) - 0.75) < 1e-9,
   "coverage is fraction containing causal")
ok(empirical_coverage(c(NA, TRUE, TRUE)) == 1.0, "NA loci are dropped")

cat("calibration_gap\n")
ok(abs(calibration_gap(0.80, 0.95) + 0.15) < 1e-9, "anti-conservative gap is negative")
ok(calibration_gap(0.96, 0.95) > 0, "conservative gap is positive")

cat("cs_false_discovery_rate\n")
# two loci: locus1 has one good set + one false set; locus2 has one good set
sets <- list(list(c(5L,6L), c(50L,51L)), list(c(10L)))
causal <- list(c(6L), c(10L))
fdr <- cs_false_discovery_rate(sets, causal)
ok(abs(fdr - (1/3)) < 1e-9, "FDR counts false sets over all sets")

cat("\nSimulated calibration recovery (the headline property)\n")
# Mock fine-mapper with KNOWN coverage: validates the metric end to end,
# mirroring the verified Python prototype.
set.seed(0)
mock <- function(n, true_cov) as.logical(runif(n) < true_cov)
emp_cal <- empirical_coverage(mock(5000, 0.95))
ok(abs(emp_cal - 0.95) < 0.02, "calibrated mock -> ~0.95 empirical coverage")
emp_bad <- empirical_coverage(mock(5000, 0.80))
ok(emp_bad < 0.85 && calibration_gap(emp_bad, 0.95) < -0.10,
   "anti-conservative mock is detected")

cat("\nAll calibration-logic tests passed.\n")
