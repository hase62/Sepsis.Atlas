#!/usr/bin/env Rscript
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1800)
if (getRversion() < "4.3" || getRversion() >= "4.4") {
  stop("This installer is intended for R 4.3.x; found ", getRversion())
}
lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (!nzchar(lib)) lib <- .libPaths()[1]
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(unique(c(lib, .libPaths())))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib = lib)
if (as.character(BiocManager::version()) != "3.18") {
  BiocManager::install(version = "3.18", ask = FALSE, update = FALSE, lib = lib)
}
needed <- c("SingleCellExperiment", "scDblFinder", "BiocParallel", "scater", "miQC")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) BiocManager::install(missing, ask = FALSE, update = FALSE, lib = lib)
cran <- c("flexmix")
missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran, lib = lib)
all_pkgs <- c(needed, cran)
for (p in all_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) stop("Package still missing: ", p)
  cat(sprintf("%-24s %s\n", p, as.character(packageVersion(p))))
}
cat("PASS: final-QC packages are available\n")
