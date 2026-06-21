# run_benchmark.R
# The full fine-mapping calibration experiment.
#
# Across many simulated loci, run SuSiE-RSS and measure whether its 95%
# credible sets are calibrated (contain the true causal variant 95% of the
# time), how sharp they are (resolution), and how calibration and power
# degrade as the genetics gets harder: stronger LD, more causal variants,
# weaker effects.
#
# Run with:  Rscript R/run_benchmark.R
# Requires:  susieR, MASS  (install.packages(c("susieR","MASS")))

source("R/simulate.R")
source("R/finemap.R")
source("R/calibration.R")

suppressMessages({ library(jsonlite) })

# Evaluate a single condition: n_rep loci with fixed difficulty parameters.
evaluate_condition <- function(label, n_rep, p, n, n_causal, pve, block_w,
                               coverage = 0.95, base_seed = 1000) {
  contains <- rep(NA, n_rep)   # NA where SuSiE abstained (no credible set)
  sizes <- numeric(n_rep)
  abstained <- logical(n_rep)
  all_sets <- vector("list", n_rep)
  all_causal <- vector("list", n_rep)

  for (i in seq_len(n_rep)) {
    sim <- simulate_locus(p = p, n = n, n_causal = n_causal,
                          pve = pve, block_w = block_w, seed = base_seed + i)
    fit <- run_susie(sim$z, sim$R, sim$params$n, L = max(5, n_causal + 2))
    sets <- get_credible_sets(fit, sim$R, coverage = coverage)
    if (length(sets) == 0) {
      # SuSiE returned no credible set: it declined to commit. This is an
      # abstention, recorded separately, NOT counted as a coverage failure.
      abstained[i] <- TRUE
      contains[i] <- NA
      sizes[i] <- NA
    } else {
      contains[i] <- cs_contains_causal(sets, sim$causal)
      sizes[i] <- mean_cs_size(sets)
    }
    all_sets[[i]] <- sets
    all_causal[[i]] <- sim$causal
  }

  # Coverage is measured AMONG loci where a credible set was reported: when
  # SuSiE makes a 95% claim, does it hold? Abstentions are reported separately.
  emp <- empirical_coverage(contains)
  list(
    label = label,
    n_rep = n_rep,
    params = list(p = p, n = n, n_causal = n_causal, pve = pve, block_w = block_w),
    nominal_coverage = coverage,
    empirical_coverage = round(emp, 4),
    calibration_gap = round(calibration_gap(emp, coverage), 4),
    mean_cs_size = round(mean(sizes, na.rm = TRUE), 2),
    abstention_rate = round(mean(abstained), 4),
    cs_fdr = round(cs_false_discovery_rate(all_sets, all_causal), 4)
  )
}

main <- function() {
  n_rep <- as.integer(Sys.getenv("FM_NREP", "1000"))
  cat(sprintf("Fine-mapping calibration benchmark (n_rep=%d per condition)\n\n", n_rep))

  conditions <- list(
    # Baseline: single causal, moderate local LD, marginal-but-clear effect.
    list("baseline_1causal", n_rep, 600, 20000, 1, 0.002, 0.9),
    # Stronger LD: tighter proxy correlation around the causal.
    list("strong_LD",        n_rep, 600, 20000, 1, 0.002, 0.97),
    # Multiple causal variants.
    list("two_causal",       n_rep, 600, 20000, 2, 0.002, 0.9),
    list("three_causal",     n_rep, 600, 20000, 3, 0.002, 0.9),
    # Weaker effect (lower power): the regime where SuSiE abstains more and
    # coverage can slip.
    list("weak_effect",      n_rep, 600, 20000, 1, 0.001, 0.9),
    # Very weak effect: deep in the abstention regime.
    list("very_weak_effect", n_rep, 600, 20000, 1, 0.0007, 0.9)
  )

  results <- lapply(conditions, function(c) do.call(evaluate_condition, c))

  cat(sprintf("%-18s %8s %10s %8s %10s %8s\n",
              "condition", "emp.cov", "gap", "size", "abstain", "cs_fdr"))
  for (r in results) {
    cat(sprintf("%-18s %8.3f %10.3f %8.2f %10.3f %8.3f\n",
                r$label, r$empirical_coverage, r$calibration_gap,
                r$mean_cs_size, r$abstention_rate, r$cs_fdr))
  }

  dir.create("evidence/reports", recursive = TRUE, showWarnings = FALSE)
  write_json(results, "evidence/reports/calibration_benchmark.json",
             auto_unbox = TRUE, pretty = TRUE)
  cat("\nReport written to evidence/reports/calibration_benchmark.json\n")
}

if (sys.nframe() == 0) main()
