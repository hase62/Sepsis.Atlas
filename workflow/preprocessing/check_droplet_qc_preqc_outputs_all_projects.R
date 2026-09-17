#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(readr)
  library(dplyr)
  library(purrr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg[[1]]), winslash = "/", mustWork = FALSE))
repo_root_default <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = FALSE)
project_root <- normalizePath(if (length(args)) args[[1]] else repo_root_default, winslash = "/", mustWork = TRUE)
default_projects <- c("GSE151263","GSE163668","GSE167363","GSE216007","GSE216020","GSE217906","GSE220189","GSE242127","GSE252331")
projects <- if (length(args) >= 2L) args[-1] else default_projects
required_meta <- c(
  "database_accession","sample_id","geo_accession","patient_id","study_group","factor",
  "is_sepsis","ards","days","timepoint_label","source_name","tissue_label",
  "cell_calling_method","cell_call_policy","emptydrops_fdr_threshold","emptydrops_lower",
  "validrops_label_dead","validrops_drop_dead","soupx_method","soupx_status",
  "soupx_contamination_fraction","n_raw_barcodes","n_nonzero_barcodes",
  "n_emptydrops_pass","n_validrops_pass","n_called_after_droplet_qc","preqc_stage"
)

missing_like <- function(x) {
  is.na(x) | (is.character(x) & !nzchar(trimws(x)))
}

one_numeric_value <- function(x) {
  z <- unique(x[is.finite(suppressWarnings(as.numeric(x)))])
  suppressWarnings(as.numeric(z))
}

rows <- list()
for (project_id in projects) {
  rdir <- file.path(project_root, "pre_integration", project_id, "rds_preqc_raw")
  files <- if (dir.exists(rdir)) list.files(rdir, pattern = "\\.rds$", full.names = TRUE) else character()
  if (!length(files)) {
    rows[[length(rows)+1]] <- tibble(project_id, file=NA_character_, pass=FALSE, problem="no preQC RDS")
    next
  }

  for (f in files) {
    rec <- tryCatch({
      obj <- readRDS(f)
      stopifnot(inherits(obj, "Seurat"))
      assays <- names(obj@assays)
      missing <- setdiff(required_meta, colnames(obj@meta.data))
      problems <- character()

      if (!all(c("RNA", "RAW") %in% assays)) {
        problems <- c(problems, "RNA/RAW assays missing")
      } else {
        if (!identical(colnames(obj[["RNA"]]), colnames(obj[["RAW"]]))) {
          problems <- c(problems, "RNA/RAW barcode order differs")
        }
        if (!identical(rownames(obj[["RNA"]]), rownames(obj[["RAW"]]))) {
          problems <- c(problems, "RNA/RAW feature order differs")
        }
      }

      if (length(missing)) {
        problems <- c(problems, paste0("missing_meta:", paste(missing, collapse=";")))
      }
      if (ncol(obj) == 0) problems <- c(problems, "zero cells")

      if (!length(missing)) {
        must_be_populated <- c(
          "cell_calling_method", "cell_call_policy", "emptydrops_fdr_threshold",
          "emptydrops_lower", "validrops_label_dead", "validrops_drop_dead",
          "soupx_method", "soupx_status", "n_raw_barcodes", "n_nonzero_barcodes",
          "n_emptydrops_pass", "n_validrops_pass", "n_called_after_droplet_qc",
          "preqc_stage"
        )
        all_na <- must_be_populated[vapply(
          must_be_populated,
          function(nm) all(missing_like(obj@meta.data[[nm]])),
          logical(1)
        )]
        if (length(all_na)) {
          problems <- c(problems, paste0("all_NA_meta:", paste(all_na, collapse=";")))
        }

        if (all(missing_like(obj$soupx_contamination_fraction))) {
          problems <- c(problems, "all_NA_meta:soupx_contamination_fraction")
        }
        if (!"emptydrops_fdr" %in% colnames(obj@meta.data) || all(is.na(obj$emptydrops_fdr))) {
          problems <- c(problems, "emptyDrops per-barcode results missing/all NA")
        }
        if (!"validrops_qc.pass" %in% colnames(obj@meta.data) || all(is.na(obj$validrops_qc.pass))) {
          problems <- c(problems, "valiDrops per-barcode results missing/all NA")
        }
        if (!isTRUE(all(!is.na(obj$preqc_stage) & obj$preqc_stage == "post_emptyDrops_validrops_SoupX_pre_scDblFinder_miQC"))) {
          problems <- c(problems, "unexpected preqc_stage")
        }
        if (!isTRUE(all(!is.na(obj$cell_calling_method) & obj$cell_calling_method == "valiDrops_with_emptyDrops_audit"))) {
          problems <- c(problems, "unexpected cell_calling_method")
        }
        if (!isTRUE(all(!is.na(obj$soupx_method) & obj$soupx_method %in% c("autoEstCont", "manual_fallback")))) {
          problems <- c(problems, "unexpected soupx_method")
        }

        n_raw <- one_numeric_value(obj$n_raw_barcodes)
        n_nonzero <- one_numeric_value(obj$n_nonzero_barcodes)
        n_ed <- one_numeric_value(obj$n_emptydrops_pass)
        n_vd <- one_numeric_value(obj$n_validrops_pass)
        n_called <- one_numeric_value(obj$n_called_after_droplet_qc)

        if (length(n_raw) != 1L || length(n_nonzero) != 1L ||
            length(n_ed) != 1L || length(n_vd) != 1L || length(n_called) != 1L) {
          problems <- c(problems, "droplet summary counts are not scalar/non-missing")
        } else {
          if (n_raw < n_called || n_nonzero < n_called) problems <- c(problems, "raw barcode count smaller than called count")
          if (n_called != ncol(obj)) problems <- c(problems, "n_called_after_droplet_qc differs from object cell count")
          if (n_ed <= 0 || n_vd <= 0) problems <- c(problems, "emptyDrops/valiDrops retained zero barcodes")
        }
      }

      tibble(
        project_id,
        file=basename(f),
        pass=!length(problems),
        problem=paste(unique(problems), collapse=" | "),
        n_cells=ncol(obj),
        n_genes=nrow(obj)
      )
    }, error=function(e) tibble(
      project_id,
      file=basename(f),
      pass=FALSE,
      problem=conditionMessage(e),
      n_cells=NA_integer_,
      n_genes=NA_integer_
    ))
    rows[[length(rows)+1]] <- rec
  }
}

out <- bind_rows(rows)
out_dir <- file.path(project_root,"pre_integration","_audit")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
write_csv(out,file.path(out_dir,"droplet_qc_preqc_audit.csv"))
print(out %>% count(project_id,pass,problem),n=Inf)
if (any(is.na(out$pass)) || any(!out$pass, na.rm = TRUE)) quit(status=1)
cat("PASS: all preQC RDS contain populated emptyDrops/valiDrops/SoupX results, corrected RNA, uncorrected RAW, and canonical metadata\n")
