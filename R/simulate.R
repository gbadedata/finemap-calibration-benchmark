# simulate.R
# Simulation harness for fine-mapping calibration.
#
# Fine-mapping calibration can ONLY be measured against known ground truth,
# and on real GWAS data the causal variant is unknown (that is why fine-mapping
# exists). Simulation is therefore the scientifically correct method, the same
# approach used to validate SuSiE, FiniMOM, SuSiEx, and SuSiE 2.0. We plant
# known causal variants in a locus with realistic linkage disequilibrium and
# generate GWAS summary statistics (z-scores) consistent with that truth.
#
# Output per locus: z-scores, the LD matrix R, the indices of the true causal
# variants, and the simulation parameters. This is exactly the input SuSiE-RSS
# consumes (susie_rss(z, R, n)).

suppressMessages({
  library(MASS)  # mvrnorm for correlated genotypes
})

# Build the LD matrix: local correlation with a hard neighbourhood cut.
# Correlation is strong within a small window around each SNP and ~0 beyond it,
# which is a clean, standard way to simulate a fine-mapping locus where the
# causal variant has a small set of genuine proxies. block_w maps to the
# within-neighbourhood correlation strength (higher = stronger local LD).
#   p             number of SNPs
#   block_w       within-neighbourhood correlation (0..1)
#   neighbourhood window half-width of non-zero correlation
make_ld_matrix <- function(p, block_w = 0.9, neighbourhood = 6) {
  idx <- seq_len(p)
  d <- abs(outer(idx, idx, "-"))
  R <- block_w ^ d
  R[d > neighbourhood] <- 0
  diag(R) <- 1
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < 1e-6) {
    R <- R + (1e-6 - min(ev)) * diag(p)
    s <- sqrt(diag(R)); R <- R / outer(s, s)
  }
  R
}

# Simulate one locus.
#   p             number of SNPs in the locus
#   n             GWAS sample size
#   n_causal      number of true causal variants
#   pve           per-causal proportion of variance explained (effect size)
#   block_w       LD block strength (higher = stronger LD, harder to localise)
#   seed          reproducibility
# Returns a list with z (z-scores), R (LD matrix), causal (true indices),
# and the parameters used.
simulate_locus <- function(p = 600, n = 20000, n_causal = 1,
                           pve = 0.002, block_w = 0.9, seed = 1) {
  set.seed(seed)
  R <- make_ld_matrix(p, block_w = block_w, neighbourhood = 6)

  causal <- sort(sample.int(p, n_causal))

  ncp <- sqrt(n * pve)
  b <- numeric(p)
  b[causal] <- ncp / sqrt(n_causal)

  mu <- as.vector(R %*% b)
  z <- as.vector(mvrnorm(1, mu = mu, Sigma = R))

  list(z = z, R = R, causal = causal,
       params = list(p = p, n = n, n_causal = n_causal,
                     pve = pve, block_w = block_w, seed = seed))
}
