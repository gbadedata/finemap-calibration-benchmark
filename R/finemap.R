# finemap.R
# Wrapper around SuSiE-RSS (the field-standard fine-mapping method, susieR).
#
# Given simulated summary statistics (z, R, n) it runs susie_rss and extracts
# the credible sets and per-SNP posterior inclusion probabilities (PIPs).
# Keeping this thin and separate means the calibration logic never depends on
# susieR internals: it only sees credible sets and a truth label.

suppressMessages({
  library(susieR)
})

# Run SuSiE-RSS on one locus.
#   z   z-scores (length p)
#   R   LD matrix (p x p)
#   n   GWAS sample size
#   L   maximum number of single effects (credible sets) SuSiE may return
# Returns the fitted susie object.
run_susie <- function(z, R, n, L = 10) {
  fit <- susie_rss(z = z, R = R, n = n, L = L,
                   estimate_residual_variance = FALSE)
  fit
}

# Extract credible sets as a list of integer vectors (SNP indices), using
# susieR's own credible-set construction at the requested coverage level.
#
# Passing Xcorr = R is essential: it activates SuSiE's purity filter, which
# prunes credible sets to SNPs that are genuinely correlated with each other
# (minimum absolute correlation, default 0.5). Without it, sets balloon with
# uncorrelated noise SNPs and credible-set size is meaningless.
get_credible_sets <- function(fit, R, coverage = 0.95) {
  cs <- susie_get_cs(fit, Xcorr = R, coverage = coverage)
  if (is.null(cs$cs)) return(list())
  cs$cs
}

get_pip <- function(fit) {
  susie_get_pip(fit)
}
