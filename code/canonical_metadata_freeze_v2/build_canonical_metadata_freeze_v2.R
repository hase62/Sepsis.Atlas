#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) args[[1]] else "."
project_root <- normalizePath(project_root, mustWork = TRUE)

suppressPackageStartupMessages(library(SeuratObject))

truthy <- function(x) {
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
}

is_missing <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    is.na(x) | !nzchar(trimws(x)) | tolower(trimws(x)) %in% c("na", "nan", "null", "none")
  } else {
    is.na(x)
  }
}

normalise_text <- function(x) {
  x <- as.character(x)
  x[is_missing(x)] <- NA_character_
  trimws(x)
}

canonical_fields <- c(
  "database_accession",
  "sample_id",
  "library_id",
  "library_unit_id",
  "donor_id",
  "patient_id",
  "study_group",
  "is_sepsis",
  "cov19",
  "ards",
  "septic_shock",
  "severity",
  "days",
  "sex",
  "age",
  "agemin",
  "agemax",
  "tissue_label",
  "source_blood_fraction",
  "frozen_or_fresh",
  "chemistry_version",
  "cellranger_or_equivalent_version",
  "reference_genome",
  "sofa",
  "peep"
)

manifest_path <- file.path(
  project_root,
  "pre_integration",
  "final_qc_library_manifest.tsv"
)

if (!file.exists(manifest_path)) {
  stop("Missing manifest: ", manifest_path)
}

manifest <- read.delim(
  manifest_path,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_manifest <- c(
  "project_id",
  "library_key",
  "include_in_atlas",
  "final_rds"
)

missing_manifest <- setdiff(required_manifest, names(manifest))
if (length(missing_manifest)) {
  stop(
    "Manifest is missing columns: ",
    paste(missing_manifest, collapse = ", ")
  )
}

manifest$include_in_atlas <- truthy(manifest$include_in_atlas)
manifest <- manifest[
  !is.na(manifest$include_in_atlas) & manifest$include_in_atlas,
  ,
  drop = FALSE
]

if (nrow(manifest) != 158L) {
  stop("Expected 158 included libraries; found ", nrow(manifest))
}

out_dir <- file.path(
  project_root,
  "pre_integration",
  "canonical_metadata"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library_rows <- vector("list", nrow(manifest))
field_rows <- list()
issue_rows <- list()
total_cells <- 0L

message("Freezing exact canonical metadata for ", nrow(manifest), " libraries...")

for (i in seq_len(nrow(manifest))) {
  m <- manifest[i, , drop = FALSE]
  path <- m$final_rds

  if (!file.exists(path)) {
    stop("Missing final RDS: ", path)
  }

  obj <- readRDS(path)
  md <- obj[[]]

  if (!identical(rownames(md), colnames(obj))) {
    stop(
      "Cell metadata order mismatch: ",
      m$project_id, " / ", m$library_key
    )
  }

  n_cells <- nrow(md)
  total_cells <- total_cells + n_cells

  row <- data.frame(
    project_id = m$project_id,
    library_key = m$library_key,
    atlas_library_id = paste(m$project_id, m$library_key, sep = "::"),
    n_cells = n_cells,
    final_rds = path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if ("project_id" %in% names(md)) {
    values <- unique(normalise_text(md$project_id))
    values <- values[!is.na(values)]
    if (length(values) != 1L || values[[1]] != m$project_id) {
      issue_rows[[length(issue_rows) + 1L]] <- data.frame(
        project_id = m$project_id,
        library_key = m$library_key,
        field = "project_id",
        issue = "manifest_identity_mismatch",
        detail = paste(values, collapse = " | "),
        stringsAsFactors = FALSE
      )
    }
  } else {
    issue_rows[[length(issue_rows) + 1L]] <- data.frame(
      project_id = m$project_id,
      library_key = m$library_key,
      field = "project_id",
      issue = "column_missing",
      detail = "",
      stringsAsFactors = FALSE
    )
  }

  if ("library_key" %in% names(md)) {
    values <- unique(normalise_text(md$library_key))
    values <- values[!is.na(values)]
    if (length(values) != 1L || values[[1]] != m$library_key) {
      issue_rows[[length(issue_rows) + 1L]] <- data.frame(
        project_id = m$project_id,
        library_key = m$library_key,
        field = "library_key",
        issue = "manifest_identity_mismatch",
        detail = paste(values, collapse = " | "),
        stringsAsFactors = FALSE
      )
    }
  } else {
    issue_rows[[length(issue_rows) + 1L]] <- data.frame(
      project_id = m$project_id,
      library_key = m$library_key,
      field = "library_key",
      issue = "column_missing",
      detail = "",
      stringsAsFactors = FALSE
    )
  }

  for (field in canonical_fields) {
    present <- field %in% names(md)

    if (!present) {
      row[[field]] <- NA_character_
      field_rows[[length(field_rows) + 1L]] <- data.frame(
        project_id = m$project_id,
        library_key = m$library_key,
        field = field,
        present = FALSE,
        column_class = NA_character_,
        n_cells = n_cells,
        n_nonmissing = 0L,
        nonmissing_fraction = 0,
        n_unique_nonmissing = 0L,
        library_constant = TRUE,
        value = NA_character_,
        stringsAsFactors = FALSE
      )
      issue_rows[[length(issue_rows) + 1L]] <- data.frame(
        project_id = m$project_id,
        library_key = m$library_key,
        field = field,
        issue = "canonical_column_missing",
        detail = "",
        stringsAsFactors = FALSE
      )
      next
    }

    x <- md[[field]]
    x_text <- normalise_text(x)
    values <- unique(x_text[!is.na(x_text)])
    n_nonmissing <- sum(!is.na(x_text))
    n_unique <- length(values)
    constant <- n_unique <= 1L

    row[[field]] <- if (n_unique == 1L) values[[1]] else NA_character_

    field_rows[[length(field_rows) + 1L]] <- data.frame(
      project_id = m$project_id,
      library_key = m$library_key,
      field = field,
      present = TRUE,
      column_class = paste(class(x), collapse = ";"),
      n_cells = n_cells,
      n_nonmissing = n_nonmissing,
      nonmissing_fraction = n_nonmissing / n_cells,
      n_unique_nonmissing = n_unique,
      library_constant = constant,
      value = if (n_unique == 1L) values[[1]] else paste(values, collapse = " | "),
      stringsAsFactors = FALSE
    )

    if (!constant) {
      issue_rows[[length(issue_rows) + 1L]] <- data.frame(
        project_id = m$project_id,
        library_key = m$library_key,
        field = field,
        issue = "multiple_values_within_library",
        detail = paste(head(values, 20L), collapse = " | "),
        stringsAsFactors = FALSE
      )
    }
  }

  if (!is.na(row$database_accession) && row$database_accession != m$project_id) {
    issue_rows[[length(issue_rows) + 1L]] <- data.frame(
      project_id = m$project_id,
      library_key = m$library_key,
      field = "database_accession",
      issue = "database_accession_not_project_id",
      detail = row$database_accession,
      stringsAsFactors = FALSE
    )
  }

  library_rows[[i]] <- row

  rm(obj, md)
  if (i %% 10L == 0L || i == nrow(manifest)) {
    message("  frozen ", i, "/", nrow(manifest))
    gc(verbose = FALSE)
  }
}

library_metadata <- do.call(rbind, library_rows)
field_audit <- do.call(rbind, field_rows)

if (length(issue_rows)) {
  issues <- do.call(rbind, issue_rows)
} else {
  issues <- data.frame(
    project_id = character(),
    library_key = character(),
    field = character(),
    issue = character(),
    detail = character(),
    stringsAsFactors = FALSE
  )
}

coverage <- aggregate(
  cbind(
    n_libraries_present = as.integer(field_audit$present),
    n_libraries_nonmissing = as.integer(field_audit$n_nonmissing > 0L),
    n_libraries_complete = as.integer(field_audit$nonmissing_fraction == 1),
    n_libraries_nonconstant = as.integer(!field_audit$library_constant),
    n_cells_nonmissing = field_audit$n_nonmissing,
    n_cells_total = field_audit$n_cells
  ) ~ project_id + field,
  data = field_audit,
  FUN = sum
)

project_library_n <- aggregate(
  library_key ~ project_id,
  data = library_metadata,
  FUN = length
)
names(project_library_n)[2] <- "n_libraries"

coverage <- merge(
  coverage,
  project_library_n,
  by = "project_id",
  all.x = TRUE,
  sort = FALSE
)
coverage$library_present_fraction <- coverage$n_libraries_present / coverage$n_libraries
coverage$library_nonmissing_fraction <- coverage$n_libraries_nonmissing / coverage$n_libraries
coverage$cell_nonmissing_fraction <- coverage$n_cells_nonmissing / coverage$n_cells_total

field_class <- aggregate(
  library_key ~ project_id + field + column_class,
  data = field_audit[field_audit$present, , drop = FALSE],
  FUN = length
)
names(field_class)[4] <- "n_libraries"

library_path <- file.path(out_dir, "canonical_library_metadata.tsv")
field_path <- file.path(out_dir, "canonical_field_audit_by_library.tsv")
coverage_path <- file.path(out_dir, "canonical_field_coverage_by_project.tsv")
class_path <- file.path(out_dir, "canonical_field_classes_by_project.tsv")
issue_path <- file.path(out_dir, "canonical_metadata_issues.tsv")
summary_path <- file.path(out_dir, "canonical_metadata_freeze_summary.tsv")

write.table(
  library_metadata,
  library_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
write.table(
  field_audit,
  field_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
write.table(
  coverage,
  coverage_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
write.table(
  field_class,
  class_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)
write.table(
  issues,
  issue_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

summary <- data.frame(
  metric = c(
    "n_included_libraries",
    "n_included_cells",
    "n_canonical_fields",
    "n_identity_or_semantic_issues",
    "n_missing_canonical_column_issues",
    "n_nonconstant_canonical_field_issues"
  ),
  value = c(
    nrow(library_metadata),
    total_cells,
    length(canonical_fields),
    nrow(issues),
    sum(issues$issue == "canonical_column_missing"),
    sum(issues$issue == "multiple_values_within_library")
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary,
  summary_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n=== CANONICAL METADATA FREEZE SUMMARY ===\n")
print(summary, row.names = FALSE)

cat("\n=== ISSUES ===\n")
if (nrow(issues)) {
  print(issues, row.names = FALSE)
} else {
  cat("None\n")
}

cat(
  "\nPASS: exact canonical metadata frozen without modifying RDS files\n",
  "Library table: ", library_path, "\n",
  "Coverage     : ", coverage_path, "\n",
  "Issues       : ", issue_path, "\n",
  sep = ""
)
