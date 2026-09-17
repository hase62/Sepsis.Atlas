#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".", mustWork=TRUE)

suppressPackageStartupMessages({
  library(SeuratObject)
  library(data.table)
})

tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
out <- file.path(
  root, "publication", "scientific_data",
  paste0("reuse_reference_preflight_v1__", tag)
)
dir.create(out, recursive=TRUE)

p05b <- file.path(
  root, "publication", "scientific_data",
  "processed_data_release_candidate_v1__20260819_194718"
)

cell_file <- file.path(
  p05b, "metadata",
  "SepsisAtlas_cell_metadata_v1__processed_data_release_candidate_v1__20260819_194718.tsv.gz"
)

lib_file <- file.path(
  p05b, "metadata",
  "SepsisAtlas_library_metadata_with_source_accessions_v1__processed_data_release_candidate_v1__20260819_194718.tsv"
)

stopifnot(file.exists(cell_file), file.exists(lib_file))

cells <- fread(cell_file)
libs  <- fread(lib_file)

# ------------------------------------------------------------
# 1. Final Atlas support by compartment/project/core identity
# ------------------------------------------------------------

fwrite(
  cells[, .N, by=.(final_compartment_v1, project_id)],
  file.path(out, "compartment_by_project_cell_counts.tsv"),
  sep="\t"
)

fwrite(
  cells[, .N, by=.(
    final_compartment_v1,
    atlas_core_identity_v1,
    project_id
  )],
  file.path(out, "core_identity_by_project_cell_counts.tsv"),
  sep="\t"
)

# ------------------------------------------------------------
# 2. Inventory metadata carried by all 158 final-QC RDS
# ------------------------------------------------------------

rds <- Sys.glob(file.path(
  root, "pre_integration", "GSE*",
  "rds_final_qc", "*.rds"
))

stopifnot(length(rds) == 158L)

inv <- list()
lib_const <- list()

for(i in seq_along(rds)) {
  f <- rds[[i]]
  project <- basename(dirname(dirname(dirname(f))))
  library_key <- sub("__final_qc\\.rds$", "", basename(f))

  z <- readRDS(f)
  md <- z@meta.data

  inv[[i]] <- rbindlist(lapply(names(md), function(nm) {
    x <- md[[nm]]
    nonmissing <- !(is.na(x) | (is.character(x) & !nzchar(trimws(x))))
    data.table(
      project_id=project,
      library_key=library_key,
      column=nm,
      class=paste(class(x), collapse=";"),
      n_cells=nrow(md),
      n_nonmissing=sum(nonmissing),
      n_unique=length(unique(x[nonmissing]))
    )
  }))

  vals <- lapply(names(md), function(nm) {
    x <- md[[nm]]
    x <- x[!(is.na(x) | (is.character(x) & !nzchar(trimws(x))))]
    u <- unique(as.character(x))
    if(length(u) == 1L) u else NA_character_
  })

  lib_const[[i]] <- data.table(
    project_id=project,
    library_key=library_key,
    setNames(vals, names(md))
  )

  rm(z, md)
  invisible(gc())
}

inventory <- rbindlist(inv, fill=TRUE)

summary_inv <- inventory[, .(
  n_libraries_present=.N,
  n_libraries_with_nonmissing=sum(n_nonmissing > 0),
  total_nonmissing_cells=sum(n_nonmissing),
  max_unique_within_library=max(n_unique, na.rm=TRUE),
  classes=paste(sort(unique(class)), collapse=";")
), by=column][order(column)]

fwrite(
  summary_inv,
  file.path(out, "final_qc_metadata_column_inventory.tsv"),
  sep="\t"
)

full_library_md <- rbindlist(lib_const, fill=TRUE)

fwrite(
  full_library_md,
  file.path(out, "library_constant_metadata_full.tsv"),
  sep="\t",
  quote=TRUE,
  na="NA"
)

# ------------------------------------------------------------
# 3. Join library-level condition/clinical fields to support
# ------------------------------------------------------------

preferred <- intersect(
  c(
    "project_id","library_key",
    "condition_binary","clinical_domain_std",
    "study_group","is_sepsis","cov19",
    "ARDS","ards","pneumonia",
    "mortality","shock_status",
    "days","timepoint_label","severity",
    "sex","age","SOFA","ventilation"
  ),
  names(full_library_md)
)

fwrite(
  unique(full_library_md[, ..preferred]),
  file.path(out, "library_condition_clinical_support.tsv"),
  sep="\t",
  quote=TRUE,
  na="NA"
)

# ------------------------------------------------------------
# 4. Existing scMerge2 model audit
# ------------------------------------------------------------

model_dir <- file.path(
  root, "pre_integration",
  "full_atlas_primary_v1",
  "scmerge2_primary_runs"
)

models <- Sys.glob(file.path(model_dir, "*__full_v1.rds"))

model_rows <- lapply(models, function(f) {
  z <- readRDS(f)
  data.table(
    file=basename(f),
    compartment=z$compartment,
    n_cells=length(z$cells),
    ruvK=z$ruvK,
    k_pseudoBulk=z$k_pseudoBulk,
    n_chosen_hvg=length(z$chosen_hvg),
    n_controls=length(z$controls),
    has_fullalpha=!is.null(z$fullalpha),
    has_M=!is.null(z$M),
    scMerge_version=z$scMerge_version
  )
})

model_audit <- rbindlist(model_rows, fill=TRUE)

fwrite(
  model_audit,
  file.path(out, "existing_scmerge2_model_audit.tsv"),
  sep="\t"
)

# ------------------------------------------------------------
# 5. Summary
# ------------------------------------------------------------

cat(
  "===== P05E0 REUSE REFERENCE PREFLIGHT =====\n",
  "status=PASS\n",
  "n_cells=", nrow(cells), "\n",
  "n_libraries=", nrow(libs), "\n",
  "n_final_qc_RDS=", length(rds), "\n",
  "n_cell_metadata_columns_P05B=", ncol(cells), "\n",
  "n_library_metadata_columns_P05B=", ncol(libs), "\n",
  "n_final_qc_metadata_columns_union=", nrow(summary_inv), "\n",
  "n_existing_scMerge2_models=", nrow(model_audit), "\n",
  "cellranger_count=NOT_ACCESSED\n",
  "expression_modified=NO\n",
  "annotation_modified=NO\n",
  sep="",
  file=file.path(out, "P05E0_SUMMARY.txt")
)

cat(readLines(file.path(out, "P05E0_SUMMARY.txt")), sep="\n")
cat("\nOUT_DIR=", out, "\n", sep="")
