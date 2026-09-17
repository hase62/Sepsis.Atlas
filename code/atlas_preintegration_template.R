#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(Matrix)
  library(yaml)
})

options(stringsAsFactors = FALSE)

# Shared Atlas QC and metadata helpers -----------------------------------------
atlas_script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
atlas_script_dir <- if (length(atlas_script_arg)) dirname(normalizePath(sub("^--file=", "", atlas_script_arg[[1]]), winslash = "/", mustWork = FALSE)) else normalizePath(getwd(), winslash = "/", mustWork = FALSE)
atlas_helper_candidates <- unique(c(file.path(atlas_script_dir, "R"), file.path(getwd(), "R")))
atlas_helper_root <- atlas_helper_candidates[file.exists(file.path(atlas_helper_candidates, "atlas_droplet_qc_pipeline.R")) & file.exists(file.path(atlas_helper_candidates, "atlas_metadata_schema.R"))]
if (length(atlas_helper_root) == 0L) stop("Cannot find shared Atlas helper scripts under R/")
source(file.path(atlas_helper_root[[1]], "atlas_droplet_qc_pipeline.R"))
source(file.path(atlas_helper_root[[1]], "atlas_metadata_schema.R"))

# -----------------------------------------------------------------------------
# Config loading
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
config_arg <- if (length(args) >= 1) args[[1]] else Sys.getenv("ATLAS_CONFIG", unset = "")
if (!nzchar(config_arg)) {
  stop(
    "Usage: Rscript atlas_preintegration_template.R /path/to/study.yaml\n",
    "or set ATLAS_CONFIG=/path/to/study.yaml"
  )
}

config_path <- normalizePath(config_arg, winslash = "/", mustWork = TRUE)
config_dir  <- dirname(config_path)
cfg <- yaml::read_yaml(config_path)

cfg_scalar <- function(x, default = NA) {
  if (is.null(x) || length(x) == 0) return(default)
  if (length(x) > 1) return(x[[1]])
  x
}

cfg_first_nonempty <- function(...) {
  vals <- list(...)
  for (v in vals) {
    if (is.null(v) || length(v) == 0) next
    if (is.character(v) && length(v) == 1 && !nzchar(v)) next
    return(v)
  }
  NA
}

resolve_rel_path <- function(base_dir, path_value) {
  if (is.null(path_value) || length(path_value) == 0 || is.na(path_value) || !nzchar(path_value)) {
    return(NA_character_)
  }
  if (grepl("^/", path_value)) return(path_value)
  normalizePath(file.path(base_dir, path_value), winslash = "/", mustWork = FALSE)
}

require_cols <- function(df, required, name) {
  miss <- setdiff(required, names(df))
  if (length(miss) > 0) {
    stop(name, " is missing required columns: ", paste(miss, collapse = ", "))
  }
}

set_scalar_columns <- function(df, values, overwrite = TRUE) {
  for (nm in names(values)) {
    val <- values[[nm]]
    if (is.null(val)) val <- NA
    if (overwrite || !(nm %in% names(df))) {
      df[[nm]] <- val
    }
  }
  df
}

ensure_columns <- function(df, columns) {
  for (nm in setdiff(columns, names(df))) df[[nm]] <- NA
  df[, columns, drop = FALSE]
}

safe_save_rds <- function(object, file) {
  if (exists("SaveSeuratRds", mode = "function")) {
    SaveSeuratRds(object, file = file)
  } else {
    saveRDS(object, file = file)
  }
}

write_session_info <- function(outfile) {
  utils::capture.output(sessionInfo(), file = outfile)
}

# -----------------------------------------------------------------------------
# Project config
# -----------------------------------------------------------------------------
project_id  <- cfg_scalar(cfg$project$project_id)
sra_study   <- cfg_scalar(cfg$project$sra_study)
bioproject  <- cfg_scalar(cfg$project$bioproject)
organism    <- cfg_scalar(cfg$project$organism, "Homo sapiens")
atlas_pipeline_version <- cfg_scalar(cfg$project$atlas_pipeline_version, "atlas_preintegration_template_v1")

if (is.na(project_id) || is.na(sra_study) || is.na(bioproject)) {
  stop("project.project_id, project.sra_study, and project.bioproject are required in YAML.")
}

project_root_cfg <- cfg_scalar(cfg$paths$project_root, ".")
project_root <- normalizePath(resolve_rel_path(config_dir, project_root_cfg), winslash = "/", mustWork = FALSE)
output_base  <- cfg_scalar(cfg$paths$output_base, "pre_integration")
input_root_cfg <- resolve_rel_path(config_dir, cfg_scalar(cfg$paths$input_root, NA_character_))

resolve_input_root <- function(project_root, sra_study, input_root_cfg = NA_character_) {
  candidates <- c(
    Sys.getenv("CELLRANGER_COUNT_ROOT", unset = ""),
    input_root_cfg,
    file.path(project_root, "cellranger_count", "output", sra_study)
  )
  candidates <- unique(candidates[nzchar(candidates) & !is.na(candidates)])
  hits <- candidates[dir.exists(candidates)]
  if (length(hits) == 0) {
    stop(
      "Could not find Cell Ranger output directory. Tried:\n",
      paste0(" - ", candidates, collapse = "\n"),
      "\nSet CELLRANGER_COUNT_ROOT or paths.input_root in the YAML."
    )
  }
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

input_root  <- resolve_input_root(project_root, sra_study, input_root_cfg)
output_root <- file.path(project_root, output_base, project_id)

rds_preqc_dir   <- file.path(output_root, "rds_preqc_raw")
rds_qc_dir      <- file.path(output_root, "rds_qcfiltered_raw")
rds_working_dir <- file.path(output_root, "rds_working")
meta_dir        <- file.path(output_root, "metadata")
qc_dir          <- file.path(output_root, "qc")
manifest_dir    <- file.path(output_root, "manifest")

for (d in c(output_root, rds_preqc_dir, rds_qc_dir, rds_working_dir, meta_dir, qc_dir, manifest_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("PROJECT_ROOT: ", project_root)
message("input_root  : ", input_root)
message("output_root : ", output_root)

# -----------------------------------------------------------------------------
# Study-specific tables
# -----------------------------------------------------------------------------
sample_map_path <- resolve_rel_path(config_dir, cfg$tables$sample_map_csv)
clinical_path   <- resolve_rel_path(config_dir, cfg$tables$clinical_csv)
run_map_path    <- resolve_rel_path(config_dir, cfg$tables$run_map_csv)

sample_map <- readr::read_csv(sample_map_path, show_col_types = FALSE)
clinical_tbl <- readr::read_csv(clinical_path, show_col_types = FALSE)
run_tbl <- readr::read_csv(run_map_path, show_col_types = FALSE)

require_cols(sample_map,
             c("biosample_accession", "experiment_accession", "geo_accession", "paper_sample"),
             "sample_map_csv")
require_cols(clinical_tbl,
             c("geo_accession", "paper_sample"),
             "clinical_csv")
require_cols(run_tbl,
             c("run_accession", "biosample_accession", "experiment_accession", "geo_accession"),
             "run_map_csv")

# -----------------------------------------------------------------------------
# Runtime settings
# -----------------------------------------------------------------------------
run_emptydrops          <- as.logical(cfg_scalar(cfg$settings$cell_calling$run_emptydrops, TRUE))
emptydrops_lower        <- as.integer(cfg_scalar(cfg$settings$cell_calling$emptydrops_lower, 100L))
emptydrops_fdr          <- as.numeric(cfg_scalar(cfg$settings$cell_calling$emptydrops_fdr, 0.001))
fallback_min_counts     <- as.integer(cfg_scalar(cfg$settings$cell_calling$fallback_min_counts, 500L))

atlas_min_features      <- as.integer(cfg_scalar(cfg$settings$qc$atlas_min_features, 200L))
atlas_max_percent_mt    <- as.numeric(cfg_scalar(cfg$settings$qc$atlas_max_percent_mt, 10))
atlas_min_cells_saved   <- as.integer(cfg_scalar(cfg$settings$qc$atlas_min_cells_saved, 500L))
apply_paper_qc_filter   <- as.logical(cfg_scalar(cfg$settings$qc$apply_paper_qc_filter, FALSE))

run_scDblFinder         <- as.logical(cfg_scalar(cfg$settings$doublet$run_scDblFinder, FALSE))
drop_doublets           <- as.logical(cfg_scalar(cfg$settings$doublet$drop_doublets, FALSE))
scDblFinder_clusters    <- as.logical(cfg_scalar(cfg$settings$doublet$scDblFinder_clusters, TRUE))

run_working_normalize   <- as.logical(cfg_scalar(cfg$settings$working_object$run_working_normalize, TRUE))
run_preliminary_annotation <- as.logical(cfg_scalar(cfg$settings$working_object$run_preliminary_annotation, FALSE))
save_objlists <- as.logical(cfg_scalar(cfg$settings$save_objlists, FALSE))

# -----------------------------------------------------------------------------
# Provenance and constant metadata
# -----------------------------------------------------------------------------
paper_md <- cfg$paper_metadata
reprocess_md <- cfg$reprocess_metadata
study_defaults <- cfg$study_defaults

reprocess_cellranger_version <- cfg_first_nonempty(
  cfg_scalar(reprocess_md$cellranger_or_equivalent_version, NA),
  Sys.getenv("REPROCESS_CELLRANGER_VERSION", unset = "")
)
reprocess_reference_genome <- cfg_first_nonempty(
  cfg_scalar(reprocess_md$reference_genome, NA),
  Sys.getenv("REPROCESS_REFERENCE_GENOME", unset = "")
)
reprocess_gene_annotation <- cfg_first_nonempty(
  cfg_scalar(reprocess_md$gene_annotation, NA),
  Sys.getenv("REPROCESS_GENE_ANNOTATION", unset = "")
)
reprocess_chemistry_version <- cfg_first_nonempty(
  cfg_scalar(reprocess_md$chemistry_version, NA),
  Sys.getenv("REPROCESS_CHEMISTRY", unset = "")
)
reprocess_command_hash <- cfg_first_nonempty(
  cfg_scalar(reprocess_md$command_hash, NA),
  Sys.getenv("REPROCESS_COMMAND_HASH", unset = "")
)
reprocess_date <- cfg_first_nonempty(
  cfg_scalar(reprocess_md$date, NA),
  format(Sys.time(), "%Y-%m-%d")
)

common_meta_columns <- c(
  "database_accession", "sra_study", "bioproject", "organism",
  "biosample_accession", "experiment_accession", "run_accession", "geo_accession",
  "sample_id", "sample_alias", "paper_sample", "patient_id", "donor_id",
  "study_group", "factor", "cov19", "is_sepsis", "ARDS", "pneumonia",
  "infection_site", "pathogen_etiology", "mortality", "shock_status", "icu_admission",
  "days", "timepoint_label", "severity", "severityatday", "ventilation",
  "sex", "age", "SOFA", "tidal_cc_perkg", "PEEP", "plateau_pressure",
  "source_name", "tissue_label", "source_blood_fraction",
  "blood_collection_tube", "pbmc_isolation_method", "ficoll", "rbc_lysis",
  "sample_storage", "frozen_or_fresh",
  "platform", "sequencer_model", "library_type", "chemistry_version",
  "feature_barcoding", "vdj_capture", "read_len",
  "aligner", "cellranger_or_equivalent_version",
  "reference_genome", "gene_annotation", "seq_depth",
  "n_loaded_cells", "paper_loaded_cells_per_sample",
  "paper_qc_min_genes_per_cell", "paper_qc_max_percent_mt", "paper_qc_max_percent_hb",
  "enrollment_window", "batch_unit_recommended", "source_note",
  "predicted.cluster",
  "atlas_pipeline_version",
  "paper_platform", "paper_sequencer_model", "paper_library_type", "paper_chemistry_version",
  "paper_cellranger_or_equivalent_version", "paper_reference_genome", "paper_gene_annotation",
  "reprocess_platform", "reprocess_sequencer_model", "reprocess_library_type",
  "reprocess_chemistry_version", "reprocess_cellranger_version",
  "reprocess_reference_genome", "reprocess_gene_annotation",
  "reprocess_date", "reprocess_command_hash",
  "cell_calling_method", "ambient_correction_method", "doublet_detection_method"
)

sample_meta <- sample_map %>%
  left_join(clinical_tbl, by = c("geo_accession", "paper_sample"))

sample_meta$database_accession <- project_id
sample_meta$sra_study <- sra_study
sample_meta$bioproject <- bioproject
sample_meta$organism <- organism
sample_meta$sample_id <- sample_meta$geo_accession
if (!("donor_id" %in% names(sample_meta)) || all(is.na(sample_meta$donor_id))) {
  sample_meta$donor_id <- paste0(project_id, "-", sample_meta$paper_sample)
}
sample_meta$atlas_pipeline_version <- atlas_pipeline_version
sample_meta$predicted.cluster <- NA_character_

# Study defaults that are constant within the study but not specific to the template
sample_meta <- set_scalar_columns(sample_meta, list(
  source_name = cfg_scalar(study_defaults$source_name, NA),
  tissue_label = cfg_scalar(study_defaults$tissue_label, NA),
  source_blood_fraction = cfg_scalar(study_defaults$source_blood_fraction, NA),
  blood_collection_tube = cfg_scalar(study_defaults$blood_collection_tube, NA),
  pbmc_isolation_method = cfg_scalar(study_defaults$pbmc_isolation_method, NA),
  ficoll = cfg_scalar(study_defaults$ficoll, NA),
  rbc_lysis = cfg_scalar(study_defaults$rbc_lysis, NA),
  sample_storage = cfg_scalar(study_defaults$sample_storage, NA),
  frozen_or_fresh = cfg_scalar(study_defaults$frozen_or_fresh, NA),
  feature_barcoding = cfg_scalar(study_defaults$feature_barcoding, NA),
  vdj_capture = cfg_scalar(study_defaults$vdj_capture, NA),
  read_len = cfg_scalar(study_defaults$read_len, NA),
  seq_depth = cfg_scalar(study_defaults$seq_depth, NA),
  n_loaded_cells = cfg_scalar(study_defaults$n_loaded_cells, NA),
  paper_loaded_cells_per_sample = cfg_scalar(study_defaults$paper_loaded_cells_per_sample, NA),
  paper_qc_min_genes_per_cell = cfg_scalar(study_defaults$paper_qc_min_genes_per_cell, NA),
  paper_qc_max_percent_mt = cfg_scalar(study_defaults$paper_qc_max_percent_mt, NA),
  paper_qc_max_percent_hb = cfg_scalar(study_defaults$paper_qc_max_percent_hb, NA),
  enrollment_window = cfg_scalar(study_defaults$enrollment_window, NA),
  batch_unit_recommended = cfg_scalar(study_defaults$batch_unit_recommended, "sample"),
  cov19 = cfg_scalar(study_defaults$cov19, NA),
  pneumonia = cfg_scalar(study_defaults$pneumonia, NA)
), overwrite = FALSE)

# Explicit provenance split
sample_meta <- set_scalar_columns(sample_meta, list(
  paper_platform = cfg_scalar(paper_md$platform, NA),
  paper_sequencer_model = cfg_scalar(paper_md$sequencer_model, NA),
  paper_library_type = cfg_scalar(paper_md$library_type, NA),
  paper_chemistry_version = cfg_scalar(paper_md$chemistry_version, NA),
  paper_cellranger_or_equivalent_version = cfg_scalar(paper_md$cellranger_or_equivalent_version, NA),
  paper_reference_genome = cfg_scalar(paper_md$reference_genome, NA),
  paper_gene_annotation = cfg_scalar(paper_md$gene_annotation, NA),
  reprocess_platform = cfg_scalar(reprocess_md$platform, cfg_scalar(paper_md$platform, NA)),
  reprocess_sequencer_model = cfg_scalar(reprocess_md$sequencer_model, cfg_scalar(paper_md$sequencer_model, NA)),
  reprocess_library_type = cfg_scalar(reprocess_md$library_type, cfg_scalar(paper_md$library_type, NA)),
  reprocess_chemistry_version = reprocess_chemistry_version,
  reprocess_cellranger_version = reprocess_cellranger_version,
  reprocess_reference_genome = reprocess_reference_genome,
  reprocess_gene_annotation = reprocess_gene_annotation,
  reprocess_date = reprocess_date,
  reprocess_command_hash = reprocess_command_hash
))

# Legacy fields describe the actual matrices being used for downstream analysis
sample_meta <- set_scalar_columns(sample_meta, list(
  platform = cfg_scalar(reprocess_md$platform, cfg_scalar(paper_md$platform, NA)),
  sequencer_model = cfg_scalar(reprocess_md$sequencer_model, cfg_scalar(paper_md$sequencer_model, NA)),
  library_type = cfg_scalar(reprocess_md$library_type, cfg_scalar(paper_md$library_type, NA)),
  chemistry_version = cfg_first_nonempty(reprocess_chemistry_version, cfg_scalar(paper_md$chemistry_version, NA)),
  aligner = cfg_scalar(reprocess_md$aligner, cfg_scalar(paper_md$aligner, "Cell Ranger")),
  cellranger_or_equivalent_version = cfg_first_nonempty(reprocess_cellranger_version, cfg_scalar(paper_md$cellranger_or_equivalent_version, NA)),
  reference_genome = cfg_first_nonempty(reprocess_reference_genome, cfg_scalar(paper_md$reference_genome, NA)),
  gene_annotation = cfg_first_nonempty(reprocess_gene_annotation, cfg_scalar(paper_md$gene_annotation, NA)),
  cell_calling_method = if (run_emptydrops) "emptyDrops" else "fallback_feature_count",
  ambient_correction_method = cfg_scalar(study_defaults$ambient_correction_method, "none"),
  doublet_detection_method = if (run_scDblFinder) "scDblFinder" else "none"
))

sample_meta <- ensure_columns(sample_meta, common_meta_columns)

run_meta <- run_tbl %>%
  left_join(sample_meta %>% select(-run_accession),
            by = c("biosample_accession", "experiment_accession", "geo_accession")) %>%
  ensure_columns(common_meta_columns)

write_csv(sample_meta, file.path(meta_dir, paste0(project_id, "_sample_metadata.csv")))
write_csv(run_meta, file.path(meta_dir, paste0(project_id, "_run_metadata.csv")))
write_session_info(file.path(manifest_dir, "sessionInfo.txt"))

config_manifest <- tibble(
  key = c(
    "project_id", "sra_study", "bioproject", "organism",
    "atlas_pipeline_version", "PROJECT_ROOT", "input_root", "output_root",
    "run_emptydrops", "emptydrops_lower", "emptydrops_fdr", "fallback_min_counts",
    "atlas_min_features", "atlas_max_percent_mt", "atlas_min_cells_saved", "apply_paper_qc_filter",
    "run_scDblFinder", "drop_doublets", "run_working_normalize",
    "reprocess_cellranger_version", "reprocess_reference_genome",
    "reprocess_gene_annotation", "reprocess_chemistry_version",
    "reprocess_command_hash", "reprocess_date", "config_path"
  ),
  value = c(
    project_id, sra_study, bioproject, organism,
    atlas_pipeline_version, project_root, input_root, output_root,
    as.character(run_emptydrops), as.character(emptydrops_lower), as.character(emptydrops_fdr), as.character(fallback_min_counts),
    as.character(atlas_min_features), as.character(atlas_max_percent_mt), as.character(atlas_min_cells_saved), as.character(apply_paper_qc_filter),
    as.character(run_scDblFinder), as.character(drop_doublets), as.character(run_working_normalize),
    reprocess_cellranger_version, reprocess_reference_genome,
    reprocess_gene_annotation, reprocess_chemistry_version,
    reprocess_command_hash, reprocess_date, config_path
  )
)
write_csv(config_manifest, file.path(manifest_dir, "config_manifest.csv"))

# -----------------------------------------------------------------------------
# Data helpers
# -----------------------------------------------------------------------------
get_assay_matrix <- function(obj, assay = "RNA", layer = "counts") {
  tryCatch({
    SeuratObject::GetAssayData(obj, assay = assay, layer = layer)
  }, error = function(e) {
    SeuratObject::GetAssayData(obj, assay = assay, slot = layer)
  })
}

compute_percent_hb <- function(obj) {
  hb_genes <- intersect(
    rownames(obj),
    c("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
  )
  if (length(hb_genes) == 0) return(rep(0, ncol(obj)))

  counts <- get_assay_matrix(obj, assay = "RNA", layer = "counts")
  denom <- Matrix::colSums(counts)
  denom[denom == 0] <- NA_real_
  out <- 100 * Matrix::colSums(counts[hb_genes, , drop = FALSE]) / denom
  out[is.na(out)] <- 0
  as.numeric(out)
}

attach_sample_metadata <- function(obj, meta_row) {
  stopifnot(nrow(meta_row) == 1)
  for (nm in colnames(meta_row)) {
    obj[[nm]] <- meta_row[[nm]][1]
  }
  obj
}

read_gene_expression <- function(raw_dir) {
  mat <- Read10X(raw_dir)
  if (is.list(mat)) {
    if ("Gene Expression" %in% names(mat)) {
      mat <- mat[["Gene Expression"]]
    } else {
      mat <- mat[[1]]
    }
  }
  mat
}

call_cells_from_raw <- function(counts, sample_name, lower = 100L, fdr = 0.001,
                                fallback_min_features = 200L, fallback_min_counts = 500L) {
  totals <- Matrix::colSums(counts)
  nfeat  <- Matrix::colSums(counts > 0)

  if (run_emptydrops &&
      requireNamespace("DropletUtils", quietly = TRUE) &&
      requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    set.seed(123)
    ed <- DropletUtils::emptyDrops(counts, lower = lower)
    retain <- !is.na(ed$FDR) & ed$FDR <= fdr
    method <- "emptyDrops"

    call_tbl <- tibble(
      barcode = colnames(counts),
      total_umi = as.numeric(totals),
      detected_features = as.numeric(nfeat),
      emptydrops_fdr = ed$FDR,
      emptydrops_logprob = ed$LogProb,
      emptydrops_limited = ed$Limited,
      cell_call_retain = retain,
      cell_calling_method = method
    )

    if (sum(retain) == 0) {
      warning(sample_name, ": emptyDrops retained zero barcodes. Falling back to feature/count thresholds.")
      retain <- nfeat >= fallback_min_features & totals >= fallback_min_counts
      call_tbl$cell_call_retain <- retain
      call_tbl$cell_calling_method <- "fallback_feature_count"
    }
  } else {
    retain <- nfeat >= fallback_min_features & totals >= fallback_min_counts
    call_tbl <- tibble(
      barcode = colnames(counts),
      total_umi = as.numeric(totals),
      detected_features = as.numeric(nfeat),
      emptydrops_fdr = NA_real_,
      emptydrops_logprob = NA_real_,
      emptydrops_limited = NA,
      cell_call_retain = retain,
      cell_calling_method = "fallback_feature_count"
    )
  }

  list(
    counts = counts[, retain, drop = FALSE],
    barcode_calls = call_tbl,
    n_input_barcodes = length(retain),
    n_retained_barcodes = sum(retain)
  )
}

run_fixed_scDblFinder <- function(obj, clusters = TRUE) {
  obj$doublet_method <- "not_run"
  obj$doublet_score <- NA_real_
  obj$doublet_call <- NA_character_

  if (!run_scDblFinder) return(obj)
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) return(obj)
  if (!requireNamespace("scDblFinder", quietly = TRUE)) return(obj)

  counts <- get_assay_matrix(obj, assay = "RNA", layer = "counts")
  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = counts))

  set.seed(123)
  sce <- scDblFinder::scDblFinder(sce, clusters = clusters)

  obj$doublet_method <- "scDblFinder"
  obj$doublet_score  <- sce$scDblFinder.score
  obj$doublet_call   <- as.character(sce$scDblFinder.class)
  obj
}

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------
sample_meta <- atlas_standardize_metadata(sample_meta, project_id = project_id)

input_dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
input_dirs <- input_dirs[basename(input_dirs) %in% sample_meta$biosample_accession]
input_dirs <- input_dirs[order(match(basename(input_dirs), sample_meta$biosample_accession))]

qc_registry <- list()
barcode_metric_registry <- list()
obj_list_qc <- list()
obj_list_working <- list()

for (samn_dir in input_dirs) {
  biosample_id <- basename(samn_dir)
  raw_dir <- file.path(samn_dir, "outs", "raw_feature_bc_matrix")
  if (!dir.exists(raw_dir)) {
    message("Skipping missing raw_feature_bc_matrix: ", raw_dir)
    next
  }

  meta_row <- sample_meta %>% filter(biosample_accession == biosample_id)
  if (nrow(meta_row) != 1) {
    stop("Sample metadata row is not unique for ", biosample_id)
  }

  geo_id <- meta_row$geo_accession[[1]]
  paper_sample <- meta_row$paper_sample[[1]]
  message("Reading ", biosample_id, " -> ", geo_id)

  res <- atlas_preprocess_10x_raw(
    raw_dir = raw_dir,
    sample_id = biosample_id,
    qc_output_dir = file.path(qc_dir, "droplet_qc"),
    project_id = project_id
  )

  obj <- CreateSeuratObject(
    counts = res$counts,
    project = geo_id,
    meta.data = res$meta,
    min.cells = 0,
    min.features = 0
  )
  obj <- atlas_add_droplet_qc_assays(obj, res)
  obj <- attach_sample_metadata(obj, meta_row)

  obj$percent.mt <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj$percent.hb <- compute_percent_hb(obj)

  obj$paper_qc_pass <- with(obj@meta.data,
    nFeature_RNA >= paper_qc_min_genes_per_cell &
      percent.mt <= paper_qc_max_percent_mt &
      percent.hb <= paper_qc_max_percent_hb
  )

  obj$atlas_qc_pass <- with(obj@meta.data,
    nFeature_RNA >= atlas_min_features &
      percent.mt < atlas_max_percent_mt
  )

  n_input_barcodes <- unique(obj$n_raw_barcodes)[1]
  n_called_cells   <- ncol(obj)

  obj <- atlas_standardize_seurat_metadata(obj, project_id = project_id)
  preqc_path <- file.path(rds_preqc_dir, paste0(geo_id, "__", biosample_id, "__preqc_raw.rds"))
  safe_save_rds(obj, preqc_path)

  qc_subset <- rep(TRUE, ncol(obj))
  if (apply_paper_qc_filter) qc_subset <- qc_subset & obj$paper_qc_pass
  qc_subset <- qc_subset & obj$atlas_qc_pass
  obj <- subset(obj, cells = colnames(obj)[qc_subset])

  n_after_qc <- ncol(obj)
  if (n_after_qc < atlas_min_cells_saved) {
    warning(geo_id, " dropped because only ", n_after_qc, " cells remained after QC.")
    qc_registry[[length(qc_registry) + 1]] <- tibble(
      biosample_accession = biosample_id,
      geo_accession = geo_id,
      paper_sample = paper_sample,
      n_input_barcodes = n_input_barcodes,
      n_called_cells = n_called_cells,
      n_after_qc = n_after_qc,
      n_doublet = NA_integer_,
      saved = FALSE,
      qcfiltered_rds = NA_character_,
      working_rds = NA_character_
    )
    next
  }

  obj <- run_fixed_scDblFinder(obj, clusters = scDblFinder_clusters)

  if (drop_doublets && "doublet_call" %in% colnames(obj@meta.data)) {
    obj <- subset(obj, cells = colnames(obj)[obj$doublet_call != "doublet"])
  }

  n_doublet <- if ("doublet_call" %in% colnames(obj@meta.data)) {
    sum(obj$doublet_call == "doublet", na.rm = TRUE)
  } else {
    NA_integer_
  }

  qcfiltered_path <- file.path(rds_qc_dir, paste0(geo_id, "__", biosample_id, "__qcfiltered_raw.rds"))
  safe_save_rds(obj, qcfiltered_path)

  working_obj <- obj
  if (run_working_normalize) {
    working_obj <- NormalizeData(working_obj, verbose = FALSE)
    working_obj <- FindVariableFeatures(working_obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  }

  if (run_preliminary_annotation) {
    working_obj$predicted.cluster <- NA_character_
  }

  working_path <- file.path(rds_working_dir, paste0(geo_id, "__", biosample_id, "__working.rds"))
  safe_save_rds(working_obj, working_path)

  cell_metrics_path <- file.path(qc_dir, paste0(geo_id, "__", biosample_id, "__cell_metrics.csv.gz"))
  write_csv(as_tibble(working_obj@meta.data, rownames = "barcode"), cell_metrics_path)

  barcode_metric_registry[[length(barcode_metric_registry) + 1]] <- tibble(
    biosample_accession = biosample_id,
    geo_accession = geo_id,
    cell_metrics_path = cell_metrics_path
  )

  qc_registry[[length(qc_registry) + 1]] <- tibble(
    biosample_accession = biosample_id,
    geo_accession = geo_id,
    paper_sample = paper_sample,
    n_input_barcodes = n_input_barcodes,
    n_called_cells = n_called_cells,
    n_after_qc = ncol(obj),
    n_doublet = n_doublet,
    saved = TRUE,
    qcfiltered_rds = qcfiltered_path,
    working_rds = working_path
  )

  if (save_objlists) {
    obj_list_qc[[geo_id]] <- obj
    obj_list_working[[geo_id]] <- working_obj
  }
}

qc_registry <- bind_rows(qc_registry)
barcode_metric_registry <- bind_rows(barcode_metric_registry)

write_csv(qc_registry, file.path(qc_dir, paste0(project_id, "_qc_summary.csv")))
write_csv(barcode_metric_registry, file.path(qc_dir, paste0(project_id, "_cell_metrics_registry.csv")))
write_csv(qc_registry %>% filter(saved), file.path(output_root, paste0(project_id, "_preintegration_registry.csv")))

if (save_objlists) {
  saveRDS(obj_list_qc, file = file.path(output_root, paste0(project_id, "_qcfiltered_objlist.rds")))
  saveRDS(obj_list_working, file = file.path(output_root, paste0(project_id, "_working_objlist.rds")))
}

message("Done.")
message("QC-filtered raw objects : ", rds_qc_dir)
message("Working objects         : ", rds_working_dir)
message("Registry                : ", file.path(output_root, paste0(project_id, "_preintegration_registry.csv")))
