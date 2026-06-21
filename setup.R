# setup.R  -- install dependencies. Run once: Rscript setup.R
pkgs <- c("susieR", "MASS", "jsonlite")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) {
  install.packages(to_install, repos = "https://cloud.r-project.org")
}
cat("Dependencies ready:\n")
for (p in pkgs) cat(sprintf("  %s: %s\n", p,
    as.character(packageVersion(p))))
