# real_locus.R
# Real-data anchor: the SORT1 / 1p13 LDL-cholesterol locus.
#
# Calibration must be measured by simulation (real loci have no known truth),
# but to show the method works on real data we use the textbook exception: the
# SORT1 locus, where the causal variant rs12740374 is FUNCTIONALLY VALIDATED
# (it creates a C/EBP transcription-factor binding site that alters hepatic
# SORT1 expression; Musunuru et al., Nature 2010) and has been fine-mapped to
# PIP > 0.99 in published work. Running SuSiE-RSS here and recovering
# rs12740374 demonstrates the pipeline on real summary statistics where the
# answer is independently known.
#
# This script expects the user to supply, on their machine:
#   - z-scores for SNPs in the 1p13 region from a public LDL-C GWAS
#     (e.g. GLGC / Willer et al. 2013, or UK Biobank LDL summary stats)
#   - an LD matrix for those SNPs from a reference panel (e.g. 1000 Genomes EUR)
#   - the index of rs12740374 in that SNP list
#
# It is intentionally a template: the data fetch is a user-machine step
# (network access, reference panel), the analysis is the code below.

source("R/finemap.R")
source("R/calibration.R")

run_real_locus <- function(z, R, n, rs12740374_index, coverage = 0.95) {
  fit <- run_susie(z, R, n, L = 5)
  pip <- get_pip(fit)
  sets <- get_credible_sets(fit, coverage = coverage)

  recovered <- cs_contains_causal(sets, rs12740374_index)
  causal_pip <- pip[rs12740374_index]

  cat(sprintf("SORT1 / 1p13 real-locus demonstration\n"))
  cat(sprintf("  rs12740374 PIP:            %.4f\n", causal_pip))
  cat(sprintf("  In a 95%% credible set:     %s\n", recovered))
  cat(sprintf("  Number of credible sets:   %d\n", length(sets)))
  cat(sprintf("  Mean credible-set size:    %.1f\n", mean_cs_size(sets)))

  list(causal_pip = causal_pip, recovered = recovered,
       n_sets = length(sets), mean_size = mean_cs_size(sets))
}
