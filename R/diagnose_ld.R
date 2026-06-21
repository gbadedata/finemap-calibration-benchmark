# diagnose_ld.R
# Fast diagnostic: run ONE locus through SuSiE for several LD settings and
# print the resulting credible-set size, so we can tune LD realism against
# real SuSiE output in seconds instead of running the full benchmark.
#
# Run:  Rscript R/diagnose_ld.R

source("R/finemap.R")
suppressMessages(library(MASS))

# Local-decay LD with a hard neighbourhood cut: strong correlation only within
# a small window around each SNP, then ~0. This gives the causal SNP few
# proxies, which should produce realistically small credible sets.
ld_local <- function(p, neighbourhood = 8, within = 0.9, decay = 0.6) {
  idx <- seq_len(p)
  d <- abs(outer(idx, idx, "-"))
  R <- within * (decay ^ d)
  R[d > neighbourhood] <- 0
  diag(R) <- 1
  ev <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < 1e-6) { R <- R + (1e-6 - min(ev)) * diag(p); s <- sqrt(diag(R)); R <- R / outer(s, s) }
  R
}

sim_one <- function(p, n, pve, ld_fun, seed = 1) {
  set.seed(seed)
  R <- ld_fun(p)
  causal <- p %/% 2
  b <- numeric(p); b[causal] <- sqrt(n * pve)
  z <- as.vector(mvrnorm(1, mu = as.vector(R %*% b), Sigma = R))
  list(z = z, R = R, causal = causal, n = n)
}

cat("Diagnostic: real SuSiE credible-set sizes for different LD settings\n")
cat("(20 loci each, single causal, p=600)\n\n")
cat(sprintf("%-38s %12s %10s\n", "setting (LD + signal)", "n_cred_sets", "set_size"))

settings <- list(
  list("r99 nb=5, pve=0.002", function(p) ld_local(p, 5, 0.995, 0.99), 0.002),
  list("r99 nb=8, pve=0.002", function(p) ld_local(p, 8, 0.995, 0.99), 0.002),
  list("r99 nb=5, pve=0.001", function(p) ld_local(p, 5, 0.995, 0.99), 0.001),
  list("r999 nb=8, pve=0.0015", function(p) ld_local(p, 8, 0.999, 0.995), 0.0015)
)

for (s in settings) {
  sizes <- c(); nsets <- c()
  for (seed in 1:20) {
    sim <- sim_one(600, 20000, s[[3]], s[[2]], seed = seed)
    fit <- run_susie(sim$z, sim$R, sim$n, L = 5)
    cs <- get_credible_sets(fit, sim$R, coverage = 0.95)
    nsets <- c(nsets, length(cs))
    if (length(cs) > 0) sizes <- c(sizes, mean(vapply(cs, length, numeric(1))))
  }
  cat(sprintf("%-38s %12.1f %10.1f\n", s[[1]], mean(nsets),
              if (length(sizes)) mean(sizes) else NA))
}
cat("\nTarget: set_size in the realistic ~3-15 range with ~1 credible set.\n")
