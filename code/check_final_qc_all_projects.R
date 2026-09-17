#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", winslash = "/", mustWork = TRUE)
projects <- if (length(args) > 1L) args[-1L] else c(
  "GSE151263", "GSE163668", "GSE167363", "GSE216007", "GSE216020",
  "GSE217906", "GSE220189", "GSE242127", "GSE252331"
)
source(file.path(root, "R", "atlas_final_qc_pipeline.R"))
suppressPackageStartupMessages(library(Seurat))

read_tsv1 <- function(path) {
  utils::read.delim(path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA")
}
rbind_fill <- function(xs) {
  if (!length(xs)) return(data.frame())
  all_names <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    miss <- setdiff(all_names, names(x))
    for (nm in miss) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, xs)
}
rows <- list()

for (p in projects) {
  pdir <- file.path(root, "pre_integration", p)
  inputs <- sort(list.files(file.path(pdir, "rds_preqc_raw"), pattern = "\\.rds$", full.names = TRUE))
  sums <- sort(list.files(file.path(pdir, "qc", "final_qc"), pattern = "__final_qc_summary\\.tsv$", full.names = TRUE))
  sumdfs <- lapply(sums, read_tsv1)
  sm <- rbind_fill(sumdfs)
  missing_summary <- if (nrow(sm)) setdiff(normalizePath(inputs, winslash = "/", mustWork = FALSE), sm$input_rds) else inputs

  invalid <- character()
  status_counts <- table(if (nrow(sm)) sm$status else character())
  if (length(missing_summary)) invalid <- c(invalid, paste0("missing_summary=", length(missing_summary)))

  if (nrow(sm)) {
    for (i in seq_len(nrow(sm))) {
      st <- sm$status[[i]]
      out <- sm$output_rds[[i]]
      if (st %in% c("pass", "review") && !is.na(out) && nzchar(out)) {
        if (!file.exists(out)) {
          invalid <- c(invalid, paste0("missing_output:", basename(out)))
          next
        }
        obj <- tryCatch(readRDS(out), error = function(e) e)
        if (inherits(obj, "error") || !inherits(obj, "Seurat")) {
          invalid <- c(invalid, paste0("unreadable_output:", basename(out)))
          next
        }
        if (!all(c("RAW", "RNA") %in% names(obj@assays))) invalid <- c(invalid, paste0("missing_assay:", basename(out)))
        raw <- atlas_get_assay_counts(obj, "RAW")
        rna <- atlas_get_assay_counts(obj, "RNA")
        if (!identical(dim(raw), dim(rna)) || !identical(colnames(raw), colnames(rna)) || !identical(rownames(raw), rownames(rna))) invalid <- c(invalid, paste0("assay_mismatch:", basename(out)))
        req <- c("scDblFinder.class", "scDblFinder.score", "scDblFinder.keep", "miQC.prob_compromised", "miQC.keep", "general_qc_keep", "final_qc_keep", "final_qc_exclusion_reason", "final_qc_stage")
        miss <- setdiff(req, colnames(obj@meta.data))
        if (length(miss)) invalid <- c(invalid, paste0("missing_metadata:", basename(out), ":", paste(miss, collapse = ",")))
        if ("final_qc_keep" %in% colnames(obj@meta.data) && !all(obj$final_qc_keep %in% TRUE)) invalid <- c(invalid, paste0("nonkeep_cell_in_final:", basename(out)))
        if ("n_final" %in% names(sm) && is.finite(sm$n_final[[i]]) && ncol(obj) != sm$n_final[[i]]) invalid <- c(invalid, paste0("n_final_mismatch:", basename(out)))
        rm(obj, raw, rna); invisible(gc())
      }
    }
  }

  n_pass <- if ("pass" %in% names(status_counts)) unname(status_counts[["pass"]]) else 0L
  n_review <- if ("review" %in% names(status_counts)) unname(status_counts[["review"]]) else 0L
  n_failed <- if ("failed" %in% names(status_counts)) unname(status_counts[["failed"]]) else 0L
  overall <- if (length(invalid) || n_failed > 0L || nrow(sm) < length(inputs)) "fail" else if (n_review > 0L) "needs_review" else "pass"
  rows[[p]] <- data.frame(
    project_id = p,
    n_input_rds = length(inputs),
    n_summary = nrow(sm),
    n_pass = n_pass,
    n_review = n_review,
    n_failed = n_failed,
    n_final_rds = length(list.files(file.path(pdir, "rds_final_qc"), pattern = "\\.rds$")),
    status = overall,
    problem = paste(unique(invalid), collapse = " | "),
    stringsAsFactors = FALSE
  )
}

result <- do.call(rbind, rows)
print(result, row.names = FALSE)
out <- file.path(root, "pre_integration", "final_qc_all_projects_summary.tsv")
atlas_write_tsv(result, out)
cat("Summary:", out, "\n")
if (any(result$status == "fail")) quit(save = "no", status = 1L)
if (any(result$status == "needs_review")) quit(save = "no", status = 2L)
cat("PASS: all final-QC outputs validated\n")
