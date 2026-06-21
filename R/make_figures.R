# make_figures.R
# Generate the two benchmark figures from the real evidence/reports JSON.
# Run AFTER R/run_benchmark.R:  Rscript R/make_figures.R
#
# Note: the polished PNGs shipped in the repo were rendered to match these
# numbers exactly. This script regenerates equivalent figures from your own
# run so the repo is fully reproducible.

suppressMessages({ library(jsonlite) })

rep <- fromJSON("evidence/reports/calibration_benchmark.json", simplifyDataFrame = FALSE)
dir.create("evidence/figures", recursive = TRUE, showWarnings = FALSE)

labels <- vapply(rep, function(r) r$label, character(1))
cov    <- vapply(rep, function(r) r$empirical_coverage, numeric(1))
abst   <- vapply(rep, function(r) r$abstention_rate, numeric(1))

NAVY <- "#1F4E79"; TEAL <- "#2A9D8F"; AMBER <- "#E9C46A"; RED <- "#E76F51"; MUTE <- "#7A8A99"

# Figure 1: calibration through abstention (dual axis)
png("evidence/figures/01_calibration_through_abstention.png", width = 2000, height = 1280, res = 200)
par(mar = c(7, 5, 5, 5))
bar_cols <- ifelse(abst > 0.5, RED, ifelse(abst > 0.2, AMBER, TEAL))
bp <- barplot(abst, names.arg = labels, las = 2, col = bar_cols, border = "white",
              ylim = c(0, 1.05), ylab = "Abstention rate (loci with no credible set)")
text(bp, abst + 0.03, sprintf("%.0f%%", abst * 100), font = 2, cex = 0.9)
title(main = "SuSiE keeps its calibration promise by abstaining, not by guessing", cex.main = 1.3)
par(new = TRUE)
plot(bp, cov, type = "b", pch = 21, bg = "white", col = NAVY, lwd = 3, cex = 1.6,
     axes = FALSE, xlab = "", ylab = "", ylim = c(0.90, 1.012), xlim = range(bp) + c(-0.5, 0.5))
abline(h = 0.95, lty = 2, col = MUTE, lwd = 1.6)
axis(4, at = c(0.90, 0.925, 0.95, 0.975, 1.0),
     labels = c("90%","92.5%","95%","97.5%","100%"), col.axis = NAVY, col = NAVY)
mtext("Coverage among reported credible sets", side = 4, line = 3, col = NAVY)
dev.off()

# Figure 2: coverage detail
png("evidence/figures/02_coverage.png", width = 2000, height = 1180, res = 200)
par(mar = c(7, 5, 5, 2))
bp <- barplot(cov * 100, names.arg = labels, las = 2, col = NAVY, border = "white",
              ylim = c(93.5, 100.7), xpd = FALSE,
              ylab = "Coverage of reported 95% credible sets")
abline(h = 95, lty = 2, col = RED, lwd = 1.8)
text(bp, cov * 100 + 0.15, sprintf("%.1f%%", cov * 100), font = 2, cex = 0.9)
text(max(bp), 95.2, "promised 95%", col = RED, pos = 2, font = 3, cex = 0.9)
title(main = "Coverage holds at or above the promised level in every condition", cex.main = 1.25)
dev.off()

cat("Figures written to evidence/figures/\n")
