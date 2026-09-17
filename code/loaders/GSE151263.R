# GSE151263 / SRP264928 / PRJNA635353
# Build pre-integration Seurat objects from existing Cell Ranger outputs (BioSample-level)
# Input layout expected:
#   ${ATLAS_CELLRANGER_OUTPUT_ROOT}/SRP264928/SAMN*/outs/raw_feature_bc_matrix
#
# Design goals
# 1) respect the actual processed location (SAMN-level Cell Ranger outputs)
# 2) carry forward the harmonized metadata/annotations used across the shared GSE loader scripts
# 3) produce QC-filtered, unintegrated objects suitable as input to downstream integration (e.g. scMerge2)
# 4) optionally use your existing helper functions/ref object when they are available in the session
#
# Notes
# - GSM<->paper sample mapping should follow GEO, and run-level grouping should follow SRA Run Selector.
# - Here we read the already-merged BioSample-level Cell Ranger outputs, so SRR-level merging is not needed.
# - This script keeps metadata fields broadly compatible with the other shared scripts by defining a common column set.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(Matrix)
})

options(stringsAsFactors = FALSE)

# Shared Atlas QC and metadata helpers -----------------------------------------
atlas_script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
atlas_script_dir <- if (length(atlas_script_arg)) {
  dirname(normalizePath(sub("^--file=", "", atlas_script_arg[[1]]), winslash = "/", mustWork = FALSE))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}
atlas_helper_candidates <- unique(c(
  file.path(atlas_script_dir, "R"),
  file.path(dirname(atlas_script_dir), "R"),
  file.path(getwd(), "R")
))
atlas_helper_root <- atlas_helper_candidates[
  file.exists(file.path(atlas_helper_candidates, "atlas_droplet_qc_pipeline.R")) &
    file.exists(file.path(atlas_helper_candidates, "atlas_metadata_schema.R"))
]
if (length(atlas_helper_root) == 0L) {
  stop("Cannot find R/atlas_droplet_qc_pipeline.R and R/atlas_metadata_schema.R")
}
source(file.path(atlas_helper_root[[1]], "atlas_droplet_qc_pipeline.R"))
source(file.path(atlas_helper_root[[1]], "atlas_metadata_schema.R"))

# -----------------------------------------------------------------------------
# Config
# -----------------------------------------------------------------------------
project_id  <- "GSE151263"
sra_study   <- "SRP264928"
bioproject  <- "PRJNA635353"
organism    <- "Homo sapiens"

# Portable path handling -------------------------------------------------------
# Default behavior:
# - PROJECT_ROOT: current working directory
# - output_root : <PROJECT_ROOT>/pre_integration/GSE151263
# - input_root  : tries, in order:
#       1) Sys.getenv("CELLRANGER_COUNT_ROOT")
#       2) <PROJECT_ROOT>/cellranger_count/output/SRP264928
#       3) legacy absolute path used previously
#
# Recommended usage:
#   setwd("/path/to/Sepsis.Atlas")
#
# Or, if Cell Ranger outputs are elsewhere:
#   Sys.setenv(CELLRANGER_COUNT_ROOT = "/path/to/cellranger_count/output/SRP264928")

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

resolve_input_root <- function(project_root, sra_study) {
  candidates <- c(
    Sys.getenv("CELLRANGER_COUNT_ROOT", unset = ""),
    file.path(project_root, "cellranger_count", "output", sra_study),
    file.path(Sys.getenv("ATLAS_CELLRANGER_OUTPUT_ROOT", unset = ""), sra_study)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[dir.exists(candidates)]
  if (length(hits) == 0) {
    stop(
      "Could not find Cell Ranger output directory. Tried:
",
      paste0(" - ", candidates, collapse = "
"),
      "
Set Sys.setenv(CELLRANGER_COUNT_ROOT = '/path/to/.../", sra_study, "') and rerun."
    )
  }
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

input_root  <- resolve_input_root(PROJECT_ROOT, sra_study)
output_root <- file.path(PROJECT_ROOT, "pre_integration", project_id)
rds_preqc_dir <- file.path(output_root, "rds_preqc_raw")
rds_qc_dir    <- file.path(output_root, "rds_qcfiltered_raw")
rds_work_dir  <- file.path(output_root, "rds_working")
meta_dir      <- file.path(output_root, "metadata")
qc_dir        <- file.path(output_root, "qc")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_preqc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
message("PROJECT_ROOT: ", PROJECT_ROOT)
message("input_root  : ", input_root)
message("output_root : ", output_root)

# If TRUE, use the same mitochondrial cutoff pattern as your shared scripts after preprocessing.
apply_atlas_mt10_filter <- TRUE
# If TRUE, retain only cells that pass the paper-level QC thresholds as well.
apply_paper_qc_filter <- FALSE
# Minimum cells kept after QC to save the object.
min_cells_after_qc <- 500L

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
ensure_columns <- function(df, columns) {
  for (nm in setdiff(columns, names(df))) df[[nm]] <- NA
  df[, columns, drop = FALSE]
}

set_all_local <- function(df, col, value) {
  df[[col]] <- value
  df
}

safe_save_rds <- function(object, file) {
  if (exists("SaveSeuratRds", mode = "function")) {
    SaveSeuratRds(object, file = file)
  } else {
    saveRDS(object, file = file)
  }
}

read_gene_expression <- function(raw_dir, sample_id) {
  atlas_preprocess_10x_raw(
    raw_dir = raw_dir,
    sample_id = sample_id,
    qc_output_dir = file.path(qc_dir, "droplet_qc"),
    project_id = project_id
  )
}


get_assay_matrix <- function(obj, assay = "RNA", layer = "counts") {
  # SeuratObject v5 uses `layer`; v4 used `slot`.
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
  100 * Matrix::colSums(counts[hb_genes, , drop = FALSE]) / Matrix::colSums(counts)
}

compute_percent_ribo <- function(obj) {
  ribo_genes <- grep("^RP[SL]", rownames(obj), value = TRUE)
  if (length(ribo_genes) == 0) return(rep(0, ncol(obj)))
  counts <- get_assay_matrix(obj, assay = "RNA", layer = "counts")
  totals <- Matrix::colSums(counts)
  as.numeric(100 * Matrix::colSums(counts[ribo_genes, , drop = FALSE]) / pmax(totals, 1))
}

attach_sample_metadata <- function(obj, meta_row) {
  atlas_attach_sample_metadata(obj, meta_row)
}

run_optional_helpers <- function(obj) {
  if (exists("rough_filter", mode = "function")) obj <- rough_filter(obj)
  if (exists("run_cellfilter", mode = "function")) obj <- run_cellfilter(obj)
  if (exists("run_miQC", mode = "function")) obj <- run_miQC(obj)
  obj
}

run_optional_singler <- function(obj) {
  if (!exists("ref", inherits = TRUE)) {
    obj$predicted.cluster <- NA_character_
    return(obj)
  }
  if (!requireNamespace("SingleR", quietly = TRUE)) {
    obj$predicted.cluster <- NA_character_
    return(obj)
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    obj$predicted.cluster <- NA_character_
    return(obj)
  }

  obj_norm <- NormalizeData(obj, verbose = FALSE)
  sce <- as.SingleCellExperiment(obj_norm, assay = "RNA")
  pred <- SingleR::SingleR(
    test = sce,
    ref = ref,
    labels = ref$label.fine,
    assay.type.test = "logcounts",
    assay.type.ref = "logcounts"
  )
  obj$predicted.cluster <- pred$labels
  obj
}

# broad harmonized metadata columns seen across the shared scripts
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
  "paper_qc_min_genes_per_cell", "paper_qc_max_percent_mt", "paper_qc_max_percent_ribo", "paper_qc_max_percent_hb",
  "enrollment_window", "batch_unit_recommended", "source_note",
  "read_mode", "predicted.cluster"
)

# -----------------------------------------------------------------------------
# Sample-level mapping
# -----------------------------------------------------------------------------
# Cell Ranger outputs are already BioSample-level, so these are the true sample units.
sample_map <- tibble::tribble(
  ~biosample_accession, ~experiment_accession, ~geo_accession, ~paper_sample, ~study_group,       ~group_label_runselector,
  "SAMN15031958",      "SRX8405330",         "GSM4569780",  "ARDS_1",     "sepsis_plus_ARDS", "sepsis+ARDS",
  "SAMN15031957",      "SRX8405331",         "GSM4569781",  "ARDS_2",     "sepsis_plus_ARDS", "sepsis+ARDS",
  "SAMN15031956",      "SRX8405332",         "GSM4569782",  "ARDS_3",     "sepsis_plus_ARDS", "sepsis+ARDS",
  "SAMN15031955",      "SRX8405333",         "GSM4569783",  "Sepsis_1",   "sepsis_only",      "sepsis",
  "SAMN15031954",      "SRX8405334",         "GSM4569784",  "Sepsis_2",   "sepsis_only",      "sepsis",
  "SAMN15031953",      "SRX8405335",         "GSM4569785",  "Sepsis_3",   "sepsis_only",      "sepsis",
  "SAMN15031952",      "SRX8405336",         "GSM4569786",  "Sepsis_4",   "sepsis_only",      "sepsis"
)

# Clinical table from the paper (Table 1), mapped to GEO according to GEO titles / SRA mapping.
clinical_tbl <- tibble::tribble(
  ~geo_accession, ~paper_sample, ~sample_alias, ~patient_id, ~days, ~timepoint_label,
  ~severity, ~severityatday, ~ventilation, ~sex, ~age, ~is_sepsis, ~ARDS,
  ~factor, ~infection_site, ~pathogen_etiology, ~source_blood_fraction,
  ~SOFA, ~tidal_cc_perkg, ~PEEP, ~plateau_pressure, ~mortality, ~source_note,

  "GSM4569780", "ARDS_1",   "PBMC from patient 1 sepsis+ARDS", "ARDS_1",   1L, "within_24h_of_ventilation",
  "severe", "severe", "Y", "F", 73L, "YES", "YES",
  "Influenza", "lung", "Influenza", "PBMC",
  7L, 8.05, 8L, 21L, NA,
  "Mapped using GEO/SRA; clinical values from paper Table 1.",

  "GSM4569781", "ARDS_2",   "PBMC from patient 2 sepsis+ARDS", "ARDS_2",   1L, "within_24h_of_ventilation",
  "severe", "severe", "Y", "F", 51L, "YES", "YES",
  NA, "lung", NA, "PBMC",
  4L, 7.32, 8L, 22L, NA,
  "Mapped using GEO/SRA; clinical values from paper Table 1.",

  "GSM4569782", "ARDS_3",   "PBMC from patient 3 sepsis+ARDS", "ARDS_3",   1L, "within_24h_of_ventilation",
  "severe", "severe", "Y", "F", 36L, "YES", "YES",
  NA, "lung", NA, "PBMC",
  4L, 5.84, 18L, 27L, NA,
  "Mapped using GEO/SRA; clinical values from paper Table 1.",

  "GSM4569783", "Sepsis_1", "PBMC from patient 1 sepsis only",  "Sepsis_1", 1L, "within_24h_of_ventilation",
  "severe", "severe", "Y", "F", 65L, "YES", "NO",
  "bacteria", "lung", "Klebsiella pneumoniae", "PBMC",
  6L, 7.09, 5L, 16L, NA,
  "Mapped using GEO/SRA; clinical values from paper Table 1.",

  "GSM4569784", "Sepsis_2", "PBMC from patient 2 sepsis only",  "Sepsis_2", 1L, "within_24h_of_ventilation",
  "severe", "severe", "Y", "M", 62L, "YES", "NO",
  "bacteria", "lung", "Klebsiella oxytoca; MSSA", "PBMC",
  7L, 7.28, 14L, 20L, NA,
  "Mapped using GEO/SRA; clinical values from paper Table 1.",

  "GSM4569785", "Sepsis_3", "PBMC from patient 3 sepsis only",  "Sepsis_3", 1L, "within_24h_of_ventilation",
  "severe", "severe", "Y", "M", 51L, "YES", "NO",
  "bacteria", "lung", "Enterobacter species", "PBMC",
  6L, 7.16, 5L, 14L, NA,
  "Mapped using GEO/SRA; clinical values from paper Table 1.",

  "GSM4569786", "Sepsis_4", "PBMC from patient 4 sepsis only",  "Sepsis_4", 1L, "within_24h_of_ventilation",
  "severe", "severe", "Y", "M", 78L, "YES", "NO",
  "bacteria", "lung", "Enterobacter aerogenes; Escherichia coli; MRSA", "PBMC",
  9L, 7.10, 5L, 16L, NA,
  "Mapped using GEO/SRA; clinical values from paper Table 1."
)

# Run-level bookkeeping for completeness (these are NOT the object units here).
run_tbl <- tibble::tribble(
  ~run_accession, ~biosample_accession, ~experiment_accession, ~geo_accession,
  "SRR11855222", "SAMN15031958",      "SRX8405330",         "GSM4569780",
  "SRR11855223", "SAMN15031958",      "SRX8405330",         "GSM4569780",
  "SRR11855224", "SAMN15031958",      "SRX8405330",         "GSM4569780",
  "SRR11855225", "SAMN15031957",      "SRX8405331",         "GSM4569781",
  "SRR11855226", "SAMN15031956",      "SRX8405332",         "GSM4569782",
  "SRR11855227", "SAMN15031956",      "SRX8405332",         "GSM4569782",
  "SRR11855228", "SAMN15031955",      "SRX8405333",         "GSM4569783",
  "SRR11855229", "SAMN15031955",      "SRX8405333",         "GSM4569783",
  "SRR11855230", "SAMN15031954",      "SRX8405334",         "GSM4569784",
  "SRR11855231", "SAMN15031954",      "SRX8405334",         "GSM4569784",
  "SRR11855232", "SAMN15031953",      "SRX8405335",         "GSM4569785",
  "SRR11855233", "SAMN15031953",      "SRX8405335",         "GSM4569785",
  "SRR11855234", "SAMN15031952",      "SRX8405336",         "GSM4569786",
  "SRR11855235", "SAMN15031952",      "SRX8405336",         "GSM4569786"
)

sample_meta <- sample_map %>%
  left_join(clinical_tbl, by = c("geo_accession", "paper_sample")) %>%
  mutate(
    database_accession = project_id,
    sra_study = sra_study,
    bioproject = bioproject,
    organism = organism,
    sample_id = geo_accession,
    donor_id = paste0(project_id, "-", paper_sample),
    source_name = "Blood",
    tissue_label = "PBMC",
    blood_collection_tube = "ACD",
    pbmc_isolation_method = "Ficoll-Paque PLUS",
    ficoll = "yes (Ficoll-Paque PLUS)",
    rbc_lysis = NA,
    sample_storage = "Frozen PBMC",
    frozen_or_fresh = "Frozen (cryopreserved)",
    platform = "10x Genomics",
    sequencer_model = "Illumina HiSeq 4000",
    library_type = "10x 3prime GEX",
    chemistry_version = "Chromium Single Cell 3' v2",
    feature_barcoding = "none",
    vdj_capture = "no",
    read_len = NA,
    aligner = "Cell Ranger",
    cellranger_or_equivalent_version = "3.0.2",
    reference_genome = "hg19",
    gene_annotation = NA,
    seq_depth = NA,
    n_loaded_cells = 7000L,
    paper_loaded_cells_per_sample = 7000L,
    paper_qc_min_genes_per_cell = 3000L,
    paper_qc_max_percent_mt = 20,
    paper_qc_max_percent_ribo = NA_real_,
    paper_qc_max_percent_hb = 10,
    enrollment_window = "within 24 hours of initiation of mechanical ventilation",
    batch_unit_recommended = "sample",
    cov19 = "Unknown",
    pneumonia = "yes"
  ) %>%
  ensure_columns(common_meta_columns)

run_meta <- run_tbl %>%
  left_join(sample_meta %>% select(-run_accession), by = c("biosample_accession", "experiment_accession", "geo_accession")) %>%
  ensure_columns(common_meta_columns)

sample_meta <- atlas_standardize_metadata(sample_meta, project_id = project_id)
run_meta <- atlas_standardize_metadata(run_meta, project_id = project_id)
write_csv(sample_meta, file.path(meta_dir, "GSE151263_sample_metadata.csv"))
write_csv(run_meta,    file.path(meta_dir, "GSE151263_run_metadata.csv"))

# -----------------------------------------------------------------------------
# Read each BioSample-level matrix and build preQC / QC-filtered / working objects
# -----------------------------------------------------------------------------
input_dirs_all <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
input_dirs <- input_dirs_all[basename(input_dirs_all) %in% sample_meta$biosample_accession]
input_dirs <- input_dirs[order(match(basename(input_dirs), sample_meta$biosample_accession))]

if (nrow(sample_meta) == 0) {
  stop("No sample metadata rows were generated.")
}
if (length(input_dirs) == 0) {
  stop(
    "No BioSample directories under input_root matched sample metadata.\n",
    "Expected examples: ", paste(head(sample_meta$biosample_accession, 10), collapse = ", "), "\n",
    "Available under input_root: ", paste(head(basename(input_dirs_all), 20), collapse = ", ")
  )
}

message("Sample metadata rows: ", nrow(sample_meta))
message("Matched BioSample dirs: ", length(input_dirs))

sample_meta <- atlas_standardize_metadata(sample_meta, project_id = project_id)

obj_list <- list()
qc_summary <- tibble(
  biosample_accession = character(),
  geo_accession = character(),
  paper_sample = character(),
  n_raw = integer(),
  n_after_qc = integer(),
  saved = logical(),
  preqc_rds_path = character(),
  qc_rds_path = character(),
  working_rds_path = character()
)

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

  message("Reading ", biosample_id, " -> ", meta_row$geo_accession)
  res <- read_gene_expression(raw_dir, sample_id = biosample_id)

  obj <- CreateSeuratObject(
    counts = res$counts,
    project = meta_row$geo_accession,
    meta.data = res$meta,
    min.cells = 0,
    min.features = 0
  )
  obj <- atlas_add_droplet_qc_assays(obj, res)
  obj <- attach_sample_metadata(obj, mutate(meta_row, read_mode = res$read_mode))

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- compute_percent_ribo(obj)
  obj[["percent.hb"]] <- compute_percent_hb(obj)

  obj[["paper_qc_pass"]] <- with(obj@meta.data,
    nFeature_RNA >= paper_qc_min_genes_per_cell &
    percent.mt <= paper_qc_max_percent_mt &
    (is.na(paper_qc_max_percent_ribo) | percent.ribo <= paper_qc_max_percent_ribo) &
    percent.hb <= paper_qc_max_percent_hb
  )
  obj[["atlas_qc_pass"]] <- with(obj@meta.data,
    nFeature_RNA >= 200 &
    percent.mt < if (apply_atlas_mt10_filter) 10 else 20
  )

  n_raw <- ncol(obj)

  # Save pre-QC raw object (counts + metadata + QC metrics, before any filtering)
  obj <- atlas_standardize_seurat_metadata(obj, project_id = project_id)
  preqc_rds <- file.path(rds_preqc_dir, paste0(meta_row$geo_accession, "__", biosample_id, "__preqc_raw.rds"))
  safe_save_rds(obj, preqc_rds)

  # optional custom helper QC from your existing pipeline
  obj <- run_optional_helpers(obj)

  # enforce explicit QC for QC-filtered tier
  qc_subset <- rep(TRUE, ncol(obj))
  if (apply_paper_qc_filter) qc_subset <- qc_subset & obj$paper_qc_pass
  qc_subset <- qc_subset & obj$atlas_qc_pass
  obj_qc <- subset(obj, cells = colnames(obj)[qc_subset])

  n_after_qc <- ncol(obj_qc)
  if (n_after_qc < min_cells_after_qc) {
    warning(meta_row$geo_accession, " dropped because only ", n_after_qc, " cells remained after QC.")
    qc_summary <- bind_rows(
      qc_summary,
      tibble(
        biosample_accession = biosample_id,
        geo_accession = meta_row$geo_accession,
        paper_sample = meta_row$paper_sample,
        n_raw = as.integer(n_raw),
        n_after_qc = as.integer(n_after_qc),
        saved = FALSE,
        preqc_rds_path = preqc_rds,
        qc_rds_path = NA_character_,
        working_rds_path = NA_character_
      )
    )
    next
  }

  qc_rds <- file.path(rds_qc_dir, paste0(meta_row$geo_accession, "__", biosample_id, "__qcfiltered_raw.rds"))
  safe_save_rds(obj_qc, qc_rds)

  # working object only; no integration or doublet filtering here
  obj_work <- NormalizeData(obj_qc, verbose = FALSE)
  obj_work <- FindVariableFeatures(obj_work, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  obj_work <- run_optional_singler(obj_work)

  work_rds <- file.path(rds_work_dir, paste0(meta_row$geo_accession, "__", biosample_id, "__working.rds"))
  safe_save_rds(obj_work, work_rds)
  obj_list[[meta_row$geo_accession]] <- obj_work

  qc_summary <- bind_rows(
    qc_summary,
    tibble(
      biosample_accession = biosample_id,
      geo_accession = meta_row$geo_accession,
      paper_sample = meta_row$paper_sample,
      n_raw = as.integer(n_raw),
      n_after_qc = as.integer(n_after_qc),
      saved = TRUE,
      preqc_rds_path = preqc_rds,
      qc_rds_path = qc_rds,
      working_rds_path = work_rds
    )
  )
}

write_csv(qc_summary, file.path(qc_dir, "GSE151263_qc_summary.csv"))

sample_registry <- qc_summary %>% dplyr::filter(saved %in% TRUE)
write_csv(sample_registry, file.path(output_root, "GSE151263_preintegration_registry.csv"))
saveRDS(obj_list, file = file.path(output_root, "GSE151263_preintegration_objlist.rds"))

message("Done. Saved pre-QC raw objects to   : ", rds_preqc_dir)
message("Done. Saved QC-filtered raw objects to: ", rds_qc_dir)
message("Done. Saved working objects to      : ", rds_work_dir)
message("Registry: ", file.path(output_root, "GSE151263_preintegration_registry.csv"))
