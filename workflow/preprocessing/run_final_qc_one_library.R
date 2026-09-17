#!/usr/bin/env Rscript

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--") || i == length(args)) stop("Invalid arguments")
    out[[sub("^--", "", key)]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- normalizePath(args[["project-root"]], winslash = "/", mustWork = TRUE)
input_rds <- normalizePath(args[["input"]], winslash = "/", mustWork = TRUE)
source(file.path(project_root, "R", "atlas_final_qc_pipeline.R"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

required <- c("SingleCellExperiment", "scDblFinder", "BiocParallel", "scater", "miQC", "flexmix")
for (p in required) atlas_require_namespace(p)

parts <- strsplit(input_rds, "/", fixed = TRUE)[[1]]
idx <- match("pre_integration", parts)
if (is.na(idx) || idx == length(parts)) stop("Input is not under pre_integration/<project>")
project_id <- parts[[idx + 1L]]
input_name <- basename(input_rds)
library_key <- sub("__preqc_raw\\.rds$", "", input_name)
if (identical(library_key, input_name)) library_key <- sub("\\.rds$", "", input_name)

project_dir <- file.path(project_root, "pre_integration", project_id)
out_dir <- file.path(project_dir, "rds_final_qc")
qc_dir <- file.path(project_dir, "qc", "final_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

output_rds <- file.path(out_dir, paste0(library_key, "__final_qc.rds"))
cell_qc_path <- file.path(qc_dir, paste0(library_key, "__final_qc_cells.csv.gz"))
summary_path <- file.path(qc_dir, paste0(library_key, "__final_qc_summary.tsv"))

write_summary_and_quit <- function(status, reason, summary_extra = list(), exit_status = 2L) {
  base <- list(
    project_id = project_id,
    library_key = library_key,
    input_rds = input_rds,
    output_rds = if (file.exists(output_rds)) output_rds else NA_character_,
    status = status,
    reason = reason,
    pipeline_version = atlas_final_qc_pipeline_version,
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  all <- c(base, summary_extra)
  df <- as.data.frame(all, stringsAsFactors = FALSE, optional = TRUE)
  atlas_write_tsv(df, summary_path)
  cat("ATLAS_FINAL_QC_RESULT\t", status, "\t", project_id, "\t", library_key, "\t", reason, "\n", sep = "")
  quit(save = "no", status = exit_status)
}

message("Project: ", project_id)
message("Library: ", library_key)
message("Input  : ", input_rds)
message("Output : ", output_rds)

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object")
if (!all(c("RAW", "RNA") %in% names(obj@assays))) stop("Input must contain RAW and RNA assays")
raw_counts <- atlas_get_assay_counts(obj, "RAW")
rna_counts <- atlas_get_assay_counts(obj, "RNA")
if (!identical(rownames(raw_counts), rownames(rna_counts))) stop("RAW and RNA gene order differs")
if (!identical(colnames(raw_counts), colnames(rna_counts))) stop("RAW and RNA barcode order differs")
if (!identical(colnames(raw_counts), rownames(obj@meta.data))) stop("Seurat metadata barcode order differs")
if (ncol(raw_counts) < 1L) stop("Input has zero cells")

md <- obj@meta.data
qc <- atlas_compute_raw_qc(raw_counts)
for (nm in setdiff(names(qc), "cell_barcode")) md[[nm]] <- qc[rownames(md), nm]
md <- atlas_add_general_qc(md)
n_input <- nrow(md)

message("Running scDblFinder on ", n_input, " called cells")
dbl <- atlas_run_scdblfinder(raw_counts, library_key)
if (!dbl$ok) {
  atlas_write_csv_gz(cbind(cell_barcode = rownames(md), md), cell_qc_path)
  write_summary_and_quit(
    "review", paste0("scDblFinder_failed: ", dbl$error),
    list(n_input = n_input, n_singlet = NA, n_doublet = NA, n_miqc_keep = NA, n_final = NA),
    exit_status = 2L
  )
}
for (nm in names(dbl$metadata)) md[[nm]] <- dbl$metadata[rownames(md), nm]
md$scDblFinder.keep <- tolower(as.character(md$scDblFinder.class)) == "singlet"
singlets <- rownames(md)[md$scDblFinder.keep %in% TRUE]
n_singlet <- length(singlets)
n_doublet <- n_input - n_singlet
message("scDblFinder: singlets=", n_singlet, ", doublets=", n_doublet)
if (n_singlet < 1L) {
  atlas_write_csv_gz(cbind(cell_barcode = rownames(md), md), cell_qc_path)
  write_summary_and_quit(
    "review", "scDblFinder_retained_zero_singlets",
    list(n_input = n_input, n_singlet = 0, n_doublet = n_doublet, n_miqc_keep = NA, n_final = 0),
    exit_status = 2L
  )
}

message("Running miQC on ", n_singlet, " singlets")
mi <- atlas_run_miqc(raw_counts, singlets, library_key)
md$miQC.prob_compromised <- NA_real_
md$miQC.keep_posterior <- FALSE
md$miQC.keep_below_boundary <- FALSE
md$miQC.keep_model <- FALSE
md$miQC.keep <- FALSE
md$miQC.auto_applied <- FALSE
md$miQC.auto_review_reason <- NA_character_
if (!mi$ok) {
  # miQC model failure is not evidence that all singlets are poor quality.
  # Retain scDblFinder singlets that pass conservative general QC, record
  # every miQC-derived field as unavailable, and force library-level review.
  miqc_failure_reason <- paste0("miQC_failed: ", mi$error)
  md[singlets, "miQC.keep_posterior"] <- NA
  md[singlets, "miQC.keep_below_boundary"] <- NA
  md[singlets, "miQC.keep_model"] <- NA
  md[singlets, "miQC.keep"] <- TRUE
  md[singlets, "miQC.auto_applied"] <- FALSE
  md[singlets, "miQC.auto_review_reason"] <- miqc_failure_reason
  md$miQC.status <- ifelse(
    md$scDblFinder.keep,
    "model_failed_singlet_retained_for_review",
    "not_evaluated_doublet"
  )
  md$miQC.posterior_cutoff <- atlas_env_num(
    "ATLAS_MIQC_POSTERIOR_CUTOFF", 0.90, min = 0, max = 1
  )

  md$final_qc_keep <- md$scDblFinder.keep & md$general_qc_keep
  md$final_qc_exclusion_reason <- atlas_make_exclusion_reasons(md)
  md$final_qc_stage <- "post_scDblFinder_miQC_failed_safe_generalQC"
  md$final_qc_pipeline_version <- atlas_final_qc_pipeline_version
  md$scDblFinder.package_version <- dbl$package_version
  md$miQC.package_version <- as.character(utils::packageVersion("miQC"))
  md$scDblFinder.seed <- dbl$seed
  md$scDblFinder.dbr_input <- dbl$dbr
  md$scDblFinder.dbr_per1k <- dbl$dbr_per1k
  md$scDblFinder.dbr_sd <- dbl$dbr_sd
  md$miQC.model_type <- NA_character_
  md$miQC.model_logLik <- NA_real_
  md$miQC.nstarts_success <- 0L
  md$miQC.posterior_removed_fraction <- NA_real_
  md$miQC.boundary_removed_fraction <- NA_real_
  md$miQC.model_removed_fraction <- NA_real_
  md$miQC.max_auto_removed_fraction <- atlas_env_num(
    "ATLAS_MIQC_MAX_AUTO_REMOVED_FRACTION", 0.50, min = 0, max = 1
  )

  keep_cells <- rownames(md)[md$final_qc_keep %in% TRUE]
  n_final <- length(keep_cells)
  doublet_fraction <- n_doublet / n_input
  final_fraction <- n_final / n_input

  atlas_write_csv_gz(cbind(cell_barcode = rownames(md), md), cell_qc_path)

  if (n_final < 1L) {
    write_summary_and_quit(
      "review",
      paste0(miqc_failure_reason, ";final_qc_retained_zero_cells"),
      list(
        n_input = n_input,
        n_singlet = n_singlet,
        n_doublet = n_doublet,
        doublet_fraction = doublet_fraction,
        n_miqc_keep_model = NA,
        n_miqc_keep = 0L,
        miqc_model_removed_fraction = NA,
        miqc_removed_fraction = NA,
        miQC_auto_applied = FALSE,
        miQC_auto_review_reason = miqc_failure_reason,
        n_general_qc_keep = sum(md$general_qc_keep, na.rm = TRUE),
        n_final = 0L,
        final_fraction = 0,
        scDblFinder_version = dbl$package_version,
        scDblFinder_dbr = dbl$dbr,
        scDblFinder_dbr_per1k = dbl$dbr_per1k,
        miQC_version = as.character(utils::packageVersion("miQC")),
        miQC_model_type = NA,
        miQC_posterior_cutoff = md$miQC.posterior_cutoff[[1]],
        miQC_logLik = NA,
        miQC_nstarts_success = 0L,
        cell_qc_csv = cell_qc_path
      ),
      exit_status = 2L
    )
  }

  obj@meta.data <- md
  obj_final <- subset(obj, cells = keep_cells)
  obj_final$final_qc_keep <- TRUE
  SeuratObject::DefaultAssay(obj_final) <- "RNA"

  partial <- paste0(output_rds, ".partial.", Sys.getpid())
  saveRDS(obj_final, partial, compress = "gzip")
  if (!file.rename(partial, output_rds)) {
    unlink(partial)
    stop("Could not atomically move final RDS to ", output_rds)
  }

  review_reasons <- miqc_failure_reason
  min_final <- atlas_env_int("ATLAS_FINAL_MIN_CELLS_REVIEW", 200L, min = 1L)
  if (n_final < min_final) {
    review_reasons <- paste0(
      review_reasons,
      ";final_cells_below_",
      min_final
    )
  }

  write_summary_and_quit(
    "review",
    review_reasons,
    list(
      n_input = n_input,
      n_singlet = n_singlet,
      n_doublet = n_doublet,
      doublet_fraction = doublet_fraction,
      n_miqc_keep_model = NA,
      n_miqc_keep = n_singlet,
      miqc_model_removed_fraction = NA,
      miqc_removed_fraction = 0,
      miQC_auto_applied = FALSE,
      miQC_auto_review_reason = miqc_failure_reason,
      miQC_posterior_removed_fraction = NA,
      miQC_boundary_removed_fraction = NA,
      n_general_qc_keep = sum(md$general_qc_keep, na.rm = TRUE),
      n_final = n_final,
      final_fraction = final_fraction,
      scDblFinder_version = dbl$package_version,
      scDblFinder_dbr = dbl$dbr,
      scDblFinder_dbr_per1k = dbl$dbr_per1k,
      miQC_version = as.character(utils::packageVersion("miQC")),
      miQC_model_type = NA,
      miQC_posterior_cutoff = md$miQC.posterior_cutoff[[1]],
      miQC_logLik = NA,
      miQC_nstarts_success = 0L,
      cell_qc_csv = cell_qc_path
    ),
    exit_status = 2L
  )
}
md[singlets, "miQC.prob_compromised"] <- mi$metadata[singlets, "miQC.prob_compromised"]
md[singlets, "miQC.keep_posterior"] <- mi$metadata[singlets, "miQC.keep_posterior"]
md[singlets, "miQC.keep_below_boundary"] <- mi$metadata[singlets, "miQC.keep_below_boundary"]
md[singlets, "miQC.keep_model"] <- mi$metadata[singlets, "miQC.keep_model"]

if (isTRUE(mi$auto_apply)) {
  md[singlets, "miQC.keep"] <- md[singlets, "miQC.keep_model"]
  md[singlets, "miQC.auto_applied"] <- TRUE
  md[singlets, "miQC.auto_review_reason"] <- "ok"
  md$miQC.status <- ifelse(md$scDblFinder.keep, "evaluated_auto_applied", "not_evaluated_doublet")
} else {
  # Conservative safeguard: retain singlets rather than discarding a majority
  # automatically. The model recommendation remains in miQC.keep_model.
  md[singlets, "miQC.keep"] <- TRUE
  md[singlets, "miQC.auto_applied"] <- FALSE
  md[singlets, "miQC.auto_review_reason"] <- mi$auto_review_reason
  md$miQC.status <- ifelse(md$scDblFinder.keep, "evaluated_auto_suppressed_review", "not_evaluated_doublet")
}
md$miQC.posterior_cutoff <- mi$posterior_cutoff

md$final_qc_keep <- md$scDblFinder.keep & md$miQC.keep & md$general_qc_keep
md$final_qc_exclusion_reason <- atlas_make_exclusion_reasons(md)
md$final_qc_stage <- "post_scDblFinder_miQC_safe_generalQC"
md$final_qc_pipeline_version <- atlas_final_qc_pipeline_version
md$scDblFinder.package_version <- dbl$package_version
md$miQC.package_version <- mi$package_version
md$scDblFinder.seed <- dbl$seed
md$scDblFinder.dbr_input <- dbl$dbr
md$scDblFinder.dbr_per1k <- dbl$dbr_per1k
md$scDblFinder.dbr_sd <- dbl$dbr_sd
md$miQC.model_type <- mi$model_type
md$miQC.model_logLik <- mi$logLik
md$miQC.nstarts_success <- mi$nstarts_success
md$miQC.posterior_removed_fraction <- mi$posterior_removed_fraction
md$miQC.boundary_removed_fraction <- mi$boundary_removed_fraction
md$miQC.model_removed_fraction <- mi$model_removed_fraction
md$miQC.max_auto_removed_fraction <- mi$max_auto_removed_fraction

n_miqc_keep_model <- sum(md$miQC.keep_model & md$scDblFinder.keep, na.rm = TRUE)
n_miqc_keep <- sum(md$miQC.keep & md$scDblFinder.keep, na.rm = TRUE)
keep_cells <- rownames(md)[md$final_qc_keep %in% TRUE]
n_final <- length(keep_cells)
message(
  "miQC model keep among singlets=", n_miqc_keep_model,
  "; applied keep=", n_miqc_keep,
  "; auto_apply=", mi$auto_apply,
  "; final cells=", n_final
)

atlas_write_csv_gz(cbind(cell_barcode = rownames(md), md), cell_qc_path)
if (n_final < 1L) {
  write_summary_and_quit(
    "review", "final_qc_retained_zero_cells",
    list(n_input = n_input, n_singlet = n_singlet, n_doublet = n_doublet,
         n_miqc_keep = n_miqc_keep, n_final = 0),
    exit_status = 2L
  )
}

obj@meta.data <- md
obj_final <- subset(obj, cells = keep_cells)
obj_final$final_qc_keep <- TRUE
SeuratObject::DefaultAssay(obj_final) <- "RNA"

partial <- paste0(output_rds, ".partial.", Sys.getpid())
saveRDS(obj_final, partial, compress = "gzip")
if (!file.rename(partial, output_rds)) {
  unlink(partial)
  stop("Could not atomically move final RDS to ", output_rds)
}

min_final <- atlas_env_int("ATLAS_FINAL_MIN_CELLS_REVIEW", 200L, min = 1L)
doublet_fraction <- n_doublet / n_input
miqc_removed_fraction <- if (n_singlet > 0) 1 - n_miqc_keep / n_singlet else NA_real_
miqc_model_removed_fraction <- if (n_singlet > 0) 1 - n_miqc_keep_model / n_singlet else NA_real_
final_fraction <- n_final / n_input
review_reasons <- character()
if (!isTRUE(mi$auto_apply)) review_reasons <- c(review_reasons, paste0("miQC_auto_filter_suppressed:", mi$auto_review_reason))
if (n_final < min_final) review_reasons <- c(review_reasons, paste0("final_cells_below_", min_final))
if (doublet_fraction > atlas_env_num("ATLAS_REVIEW_MAX_DOUBLET_FRACTION", 0.30, min = 0, max = 1)) review_reasons <- c(review_reasons, "high_doublet_fraction")
if (is.finite(miqc_model_removed_fraction) && miqc_model_removed_fraction > atlas_env_num("ATLAS_REVIEW_MAX_MIQC_REMOVED_FRACTION", 0.50, min = 0, max = 1)) review_reasons <- c(review_reasons, "high_miQC_model_removed_fraction")
if (final_fraction < atlas_env_num("ATLAS_REVIEW_MIN_FINAL_FRACTION", 0.20, min = 0, max = 1)) review_reasons <- c(review_reasons, "low_final_fraction")

status <- if (length(review_reasons)) "review" else "pass"
reason <- if (length(review_reasons)) paste(review_reasons, collapse = ";") else "ok"
summary_extra <- list(
  n_input = n_input,
  n_singlet = n_singlet,
  n_doublet = n_doublet,
  doublet_fraction = doublet_fraction,
  n_miqc_keep_model = n_miqc_keep_model,
  n_miqc_keep = n_miqc_keep,
  miqc_model_removed_fraction = miqc_model_removed_fraction,
  miqc_removed_fraction = miqc_removed_fraction,
  miQC_auto_applied = mi$auto_apply,
  miQC_auto_review_reason = mi$auto_review_reason,
  miQC_posterior_removed_fraction = mi$posterior_removed_fraction,
  miQC_boundary_removed_fraction = mi$boundary_removed_fraction,
  n_general_qc_keep = sum(md$general_qc_keep, na.rm = TRUE),
  n_final = n_final,
  final_fraction = final_fraction,
  scDblFinder_version = dbl$package_version,
  scDblFinder_dbr = dbl$dbr,
  scDblFinder_dbr_per1k = dbl$dbr_per1k,
  miQC_version = mi$package_version,
  miQC_model_type = mi$model_type,
  miQC_posterior_cutoff = mi$posterior_cutoff,
  miQC_logLik = mi$logLik,
  miQC_nstarts_success = mi$nstarts_success,
  cell_qc_csv = cell_qc_path
)
write_summary_and_quit(status, reason, summary_extra, exit_status = if (status == "pass") 0L else 2L)
