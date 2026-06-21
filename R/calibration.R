# calibration.R
# Recovery and calibration metrics. This is the scientific heart of the
# project and the part that mirrors the calibration work in the companion
# variant-calling and ClinVar benchmarks: it asks whether a stated confidence
# level is honest. Here the stated confidence is the credible-set coverage
# level, and the question is whether a 95% credible set actually contains the
# true causal variant 95% of the time.
#
# The measurement logic was prototyped and unit-tested in Python before being
# ported here, so the arithmetic is known-correct independently of susieR.

# Does ANY credible set contain ANY of the true causal variants?
# (Per-locus recovery: did fine-mapping capture the truth at all.)
cs_contains_causal <- function(credible_sets, causal) {
  if (length(credible_sets) == 0) return(FALSE)
  for (cs in credible_sets) {
    if (any(causal %in% cs)) return(TRUE)
  }
  FALSE
}

# Mean credible-set size (resolution): smaller is sharper. Averaged over the
# sets returned for a locus.
mean_cs_size <- function(credible_sets) {
  if (length(credible_sets) == 0) return(NA_real_)
  mean(vapply(credible_sets, length, numeric(1)))
}

# Empirical coverage across many loci: the fraction whose credible sets
# contained a true causal variant. For a calibrated method at coverage=0.95
# this should be ~0.95.
empirical_coverage <- function(contains_flags) {
  flags <- contains_flags[!is.na(contains_flags)]
  if (length(flags) == 0) return(0.0)
  mean(flags)
}

# Signed calibration gap: empirical minus nominal. Negative means
# anti-conservative (credible sets miss the causal more often than promised),
# which is the dangerous direction; positive means conservative.
calibration_gap <- function(empirical, nominal) {
  empirical - nominal
}

# A false discovery at the credible-set level: a returned credible set that
# contains NONE of the true causal variants. FDR = such sets / all sets.
cs_false_discovery_rate <- function(all_locus_sets, all_locus_causal) {
  total_sets <- 0L
  false_sets <- 0L
  for (i in seq_along(all_locus_sets)) {
    sets <- all_locus_sets[[i]]
    causal <- all_locus_causal[[i]]
    for (cs in sets) {
      total_sets <- total_sets + 1L
      if (!any(causal %in% cs)) false_sets <- false_sets + 1L
    }
  }
  if (total_sets == 0L) return(NA_real_)
  false_sets / total_sets
}
