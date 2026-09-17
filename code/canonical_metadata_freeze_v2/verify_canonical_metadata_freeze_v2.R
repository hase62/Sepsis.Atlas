#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) args[[1]] else "."
project_root <- normalizePath(project_root, mustWork = TRUE)
out_dir <- file.path(project_root, "pre_integration", "canonical_metadata")

required <- c(
  "canonical_library_metadata.tsv",
  "canonical_field_audit_by_library.tsv",
  "canonical_field_coverage_by_project.tsv",
  "canonical_field_classes_by_project.tsv",
  "canonical_metadata_issues.tsv",
  "canonical_metadata_freeze_summary.tsv"
)

paths <- file.path(out_dir, required)
missing <- paths[!file.exists(paths)]
if (length(missing)) {
  stop("Missing outputs:\n", paste(missing, collapse = "\n"))
}

summary <- read.delim(
  file.path(out_dir, "canonical_metadata_freeze_summary.tsv"),
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

library_metadata <- read.delim(
  file.path(out_dir, "canonical_library_metadata.tsv"),
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

issues <- read.delim(
  file.path(out_dir, "canonical_metadata_issues.tsv"),
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

value <- function(metric) {
  x <- summary$value[summary$metric == metric]
  if (!length(x)) NA_character_ else x[[1]]
}

stopifnot(as.integer(value("n_included_libraries")) == 158L)
stopifnot(as.integer(value("n_included_cells")) == 665816L)
stopifnot(nrow(library_metadata) == 158L)
stopifnot(!anyDuplicated(library_metadata$atlas_library_id))

fatal_issue <- issues$issue %in% c(
  "manifest_identity_mismatch",
  "database_accession_not_project_id",
  "multiple_values_within_library"
)

if (any(fatal_issue)) {
  print(issues[fatal_issue, , drop = FALSE], row.names = FALSE)
  stop("Fatal canonical metadata integrity issues detected")
}

cat(
  "PASS: canonical metadata freeze is internally consistent\n",
  "Included libraries: ", nrow(library_metadata), "\n",
  "Included cells    : ", value("n_included_cells"), "\n",
  "Missing canonical columns reported: ",
  value("n_missing_canonical_column_issues"), "\n",
  sep = ""
)
