# GSE163668 / SRP299788 / PRJNA690170
# Build pre-integration Seurat objects from existing Cell Ranger outputs
# for the UNPOOLED libraries only.
#
# External inputs kept as tables:
#   configs/GSE163668/GSE163668_SraRunTable.csv
#   configs/GSE163668/GSE163668_COMET_10X_CLINICAL_SCORES_PAPER.csv
#
# Embedded in this script:
#   - GEO title / pooled-vs-unpooled structure
#   - Paper / supplement-derived experimental metadata
#   - Paper QC thresholds
#
# Design notes
# - object unit here is BioSample (SAMN*) because Cell Ranger outputs already exist
#   at the BioSample level under SRP299788/SAMN*/outs/raw_feature_bc_matrix.
# - pooled GSMs are NOT turned into subject-level objects here.
#   They are written to a pooled manifest only.
# - this script is intentionally resilient to header differences in the external CSVs.
#
# Paper context:
# - Whole blood preserved workflow, RBC lysis, fresh samples, EDTA tube
# - 15,000 cells per individual loaded; some samples pooled before GEM generation
# - Chromium Single Cell 5' v5.1; NovaSeq 6000
# - Paper processing: Cell Ranger 3.0.2, GRCh38, Ensembl v85
# - Paper QC: genes >= 100, percent.mt <= 20, percent.ribo <= 50

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

# ------------------------------------------------------------------------------
# Config
# ------------------------------------------------------------------------------
project_id  <- "GSE163668"
sra_study   <- "SRP299788"
bioproject  <- "PRJNA690170"
organism    <- "Homo sapiens"

apply_atlas_mt10_filter <- TRUE
apply_paper_qc_filter   <- FALSE
min_cells_after_qc      <- 500L

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
      "Could not find Cell Ranger output directory. Tried:\n",
      paste0(" - ", candidates, collapse = "\n"),
      "\nSet Sys.setenv(CELLRANGER_COUNT_ROOT = '/path/to/.../", sra_study, "') and rerun."
    )
  }
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

find_first_existing <- function(paths) {
  hits <- paths[file.exists(paths)]
  if (length(hits) == 0) return(NA_character_)
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

input_root <- resolve_input_root(PROJECT_ROOT, sra_study)

sra_run_table <- find_first_existing(c(
  file.path(PROJECT_ROOT, "configs", "GSE163668", "GSE163668_SraRunTable.csv"),
  file.path(PROJECT_ROOT, "configs", "GSE163668_unpooled", "GSE163668_SraRunTable.csv"),
  file.path(PROJECT_ROOT, "GSE163668_SraRunTable.csv")
))
clinical_csv <- find_first_existing(c(
  file.path(PROJECT_ROOT, "configs", "GSE163668", "GSE163668_COMET_10X_CLINICAL_SCORES_PAPER.csv"),
  file.path(PROJECT_ROOT, "configs", "GSE163668_unpooled", "GSE163668_COMET_10X_CLINICAL_SCORES_PAPER.csv"),
  file.path(PROJECT_ROOT, "GSE163668_COMET_10X_CLINICAL_SCORES_PAPER.csv")
))

if (is.na(sra_run_table)) stop("Could not find GSE163668_SraRunTable.csv under configs/GSE163668 or project root.")
if (is.na(clinical_csv))  stop("Could not find GSE163668_COMET_10X_CLINICAL_SCORES_PAPER.csv under configs/GSE163668 or project root.")

output_root <- file.path(PROJECT_ROOT, "pre_integration", project_id)
rds_preqc_dir <- file.path(output_root, "rds_preqc_raw")
rds_qc_dir    <- file.path(output_root, "rds_qcfiltered_raw")
rds_work_dir  <- file.path(output_root, "rds_working")
meta_dir      <- file.path(output_root, "metadata")
qc_dir        <- file.path(output_root, "qc")
manifest_dir  <- file.path(output_root, "manifest")

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_preqc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)

message("PROJECT_ROOT      : ", PROJECT_ROOT)
message("input_root        : ", input_root)
message("SRA run table     : ", sra_run_table)
message("COMET clinical csv: ", clinical_csv)
message("output_root       : ", output_root)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
normalize_names <- function(x) {
  x <- enc2utf8(x)
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- trimws(x)
  x <- gsub("[[:space:]]+", "_", x)
  x <- gsub("[^[:alnum:]_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

safe_save_rds <- function(object, file) {
  if (exists("SaveSeuratRds", mode = "function")) {
    SaveSeuratRds(object, file = file)
  } else {
    saveRDS(object, file = file)
  }
}

ensure_columns <- function(df, columns) {
  for (nm in setdiff(columns, names(df))) df[[nm]] <- NA
  df[, columns, drop = FALSE]
}

find_col <- function(nms, exact = NULL, regex = NULL) {
  if (!is.null(exact)) {
    for (nm in exact) if (nm %in% nms) return(nm)
  }
  if (!is.null(regex)) {
    hit <- grep(regex, nms, perl = TRUE, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

pick_column <- function(df, exact = NULL, regex = NULL, required = FALSE, label = NULL) {
  nms <- names(df)
  hit <- find_col(nms, exact = exact, regex = regex)
  if (is.na(hit) && required) {
    stop(
      sprintf(
        "Could not infer required column%s. Available normalized headers are:\n%s",
        if (!is.null(label)) paste0(" for ", label) else "",
        paste(nms, collapse = ", ")
      )
    )
  }
  hit
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
  counts <- get_counts_matrix(obj, assay = "RNA")
  totals <- Matrix::colSums(counts)
  as.numeric(100 * Matrix::colSums(counts[hb_genes, , drop = FALSE]) / pmax(totals, 1))
}

get_counts_matrix <- function(obj, assay = "RNA") {
  tryCatch({
    SeuratObject::GetAssayData(obj, assay = assay, layer = "counts")
  }, error = function(e) {
    SeuratObject::GetAssayData(obj, assay = assay, slot = "counts")
  })
}

compute_percent_ribo <- function(obj) {
  ribo_genes <- grep("^RP[SL]", rownames(obj), value = TRUE)
  if (length(ribo_genes) == 0) return(rep(0, ncol(obj)))
  counts <- get_counts_matrix(obj, assay = "RNA")
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

common_meta_columns <- c(
  "database_accession", "sra_study", "bioproject", "organism",
  "biosample_accession", "experiment_accession", "run_accession", "geo_accession",
  "sample_id", "sample_alias", "paper_sample", "patient_id", "donor_id",
  "study_group", "factor", "cov19", "is_sepsis", "ards", "pneumonia",
  "infection_site", "pathogen_etiology", "mortality", "shock_status", "icu_admission",
  "days", "timepoint_label", "severity", "severityatday", "ventilation",
  "sex", "age", "sofa", "rbc_count", "hgb", "wbc_count", "plt_count",
  "source_name", "tissue_label", "source_blood_fraction",
  "blood_collection_tube", "pbmc_isolation_method", "ficoll", "rbc_lysis",
  "sample_storage", "frozen_or_fresh",
  "platform", "sequencer_model", "library_type", "chemistry_version",
  "feature_barcoding", "vdj_capture", "read_len",
  "aligner", "cellranger_or_equivalent_version",
  "reference_genome", "gene_annotation", "seq_depth",
  "n_loaded_cells", "paper_loaded_cells_per_sample",
  "paper_qc_min_genes_per_cell", "paper_qc_max_percent_mt",
  "paper_qc_max_percent_ribo", "paper_qc_max_percent_hb",
  "enrollment_window", "batch_unit_recommended", "source_note",
  "pool_id", "pool_members", "pooled_flag", "demux_status",
  "read_mode", "predicted.cluster"
)

# ------------------------------------------------------------------------------
# GEO / paper embedded mapping
# ------------------------------------------------------------------------------
geo_title_tbl <- tibble::tribble(
  ~geo_accession, ~geo_title, ~paper_sample, ~pooled_flag,
  "GSM4995425", "Pooled 10X GEX libraries for Patients 1, 2, and 3", "pool_patients_1_2_3", TRUE,
  "GSM4995426", "Pooled 10X GEX libraries for Patients 5 and 6", "pool_patients_5_6", TRUE,
  "GSM4995427", "Pooled 10X GEX libraries for Patients 7, 8, and 9", "pool_patients_7_8_9", TRUE,
  "GSM4995428", "Pooled 10X GEX libraries for Patient 10", "pool_patient_10", TRUE,
  "GSM4995429", "Pooled 10X GEX libraries for Patients 17, 20, and 21", "pool_patients_17_20_21", TRUE,
  "GSM4995430", "Pooled 10X GEX libraries for Patients 50 and 51", "pool_patients_50_51", TRUE,
  "GSM4995431", "10X GEX library for Patient 12", "Patient_12", FALSE,
  "GSM4995432", "10X GEX library for Patient 14", "Patient_14", FALSE,
  "GSM4995433", "10X GEX library for Patient 16", "Patient_16", FALSE,
  "GSM4995434", "10X GEX library for Patient 25", "Patient_25", FALSE,
  "GSM4995435", "10X GEX library for Patient 26", "Patient_26", FALSE,
  "GSM4995436", "10X GEX library for Patient 29", "Patient_29", FALSE,
  "GSM4995437", "10X GEX library for Patient 31", "Patient_31", FALSE,
  "GSM4995438", "10X GEX library for Patient 38", "Patient_38", FALSE,
  "GSM4995439", "10X GEX library for Patient 44", "Patient_44", FALSE,
  "GSM4995440", "10X GEX library for Patient 45", "Patient_45", FALSE,
  "GSM4995441", "10X GEX library for Patient 46", "Patient_46", FALSE,
  "GSM4995442", "10X GEX library for Patient 47", "Patient_47", FALSE,
  "GSM4995443", "10X GEX library for Patient 49", "Patient_49", FALSE,
  "GSM4995444", "10X GEX library for Patient 52", "Patient_52", FALSE,
  "GSM4995445", "10X GEX library for Patient 54", "Patient_54", FALSE,
  "GSM4995446", "10X GEX library for Patient 55", "Patient_55", FALSE,
  "GSM4995447", "10X GEX library for Patient 60", "Patient_60", FALSE,
  "GSM4995448", "10X GEX library for Patient 62", "Patient_62", FALSE,
  "GSM4995449", "10X GEX library for Donor 3", "Donor_3", FALSE,
  "GSM4995450", "10X GEX library for Donor 4", "Donor_4", FALSE,
  "GSM4995451", "10X GEX library for Donor 5", "Donor_5", FALSE,
  "GSM4995452", "10X GEX library for Donor 6", "Donor_6", FALSE,
  "GSM4995453", "10X GEX library for Donor 7", "Donor_7", FALSE,
  "GSM4995454", "10X GEX library for Donor 14", "Donor_14", FALSE,
  "GSM4995455", "10X GEX library for Donor 19", "Donor_19", FALSE,
  "GSM4995456", "10X GEX library for Donor 20", "Donor_20", FALSE,
  "GSM4995457", "10X GEX library for Donor 23", "Donor_23", FALSE,
  "GSM4995458", "10X GEX library for Donor 29", "Donor_29", FALSE,
  "GSM4995459", "10X GEX library for Donor 30", "Donor_30", FALSE,
  "GSM4995460", "10X GEX library for Donor 31", "Donor_31", FALSE,
  "GSM4995461", "10X GEX library for Donor 32", "Donor_32", FALSE,
  "GSM4995462", "10X GEX library for Donor 33", "Donor_33", FALSE
)

paper_clinical_tbl <- tibble::tribble(
  ~comet_id, ~covid_status_paper, ~severity_paper, ~age_paper, ~sex_paper,
  ~days_from_onset_paper, ~other_infection_paper, ~icu_admission_paper,
  ~days_ventilated_paper, ~discharged_paper,
  "1008", "NEG", "Mild/Moderate", 35.7, "M", 20, "Bacterial", "No", 0, "Yes",
  "1016", "NEG", "Mild/Moderate", 63.5, "F", 4, "None", "No", 0, "Yes",
  "1017", "NEG", "Mild/Moderate", 65.2, "F", 9, "None", "No", 0, "Yes",
  "1045", "NEG", "Mild/Moderate", 82.6, "M", NA, "None", "Yes", 0, "Yes",
  "1049", "NEG", "Mild/Moderate", 68.8, "F", 1, "Bacterial", "No", 0, "Yes",
  "1062", "NEG", "Mild/Moderate", 75.8, "F", NA, "Bacterial", "No", 0, "Yes",
  "1009", "NEG", "Severe", 67.0, "M", NA, "Bacterial", "Yes", 5, "Yes",
  "1010", "NEG", "Severe", 50.7, "F", 9, "Bacterial", "Yes", 5, "Yes",
  "1020", "NEG", "Severe", 78.6, "M", 5, "Bacterial", "Yes", 2, "Yes",
  "1021", "NEG", "Severe", 77.0, "F", 2, "Bacterial", "Yes", 5, "Yes",
  "1046", "NEG", "Severe", 49.6, "M", 3, "Bacterial", "Yes", 1, "Yes",
  "1005", "POS", "Mild/Moderate", 83.4, "M", 9, "None", "No", 0, "Yes",
  "1006", "POS", "Mild/Moderate", 71.9, "F", 16, "None", "No", 0, "Yes",
  "1007", "POS", "Mild/Moderate", 55.3, "M", 13, "Bacterial", "No", 0, "Yes",
  "1012", "POS", "Mild/Moderate", 41.1, "M", 3, "None", "No", 0, "Yes",
  "1014", "POS", "Mild/Moderate", 58.6, "M", 4, "None", "Yes", 0, "Yes",
  "1025", "POS", "Mild/Moderate", 73.2, "M", NA, "None", "No", 0, "Yes",
  "1026", "POS", "Mild/Moderate", 45.5, "F", 6, "None", "No", 0, "Yes",
  "1029", "POS", "Mild/Moderate", 85.0, "F", 2, "Viral", "No", 0, "Yes",
  "1044", "POS", "Mild/Moderate", 73.9, "M", NA, "None", "No", 0, "Yes",
  "1052", "POS", "Mild/Moderate", 37.9, "F", 4, "Viral", "No", 0, "Yes",
  "1054", "POS", "Mild/Moderate", 25.8, "M", NA, "None", "No", 0, "Yes",
  "1001", "POS", "Severe", 34.8, "F", 12, "Viral", "Yes", 10, "Yes",
  "1002", "POS", "Severe", 79.6, "M", 4, "Bacterial", "Yes", 7, "Yes",
  "1003", "POS", "Severe", 43.6, "F", 13, "Viral", "Yes", 5, "Yes",
  "1031", "POS", "Severe", 44.2, "M", 20, "Viral", "Yes", 28, "Yes",
  "1038", "POS", "Severe", 61.4, "F", 14, "None", "Yes", NA, "Yes",
  "1047", "POS", "Severe", 48.3, "F", 6, "Viral + Bacterial", "Yes", 28, "Yes",
  "1050", "POS", "Severe", 55.0, "M", 13, "Bacterial", "Yes", 29, "Yes",
  "1051", "POS", "Severe", 63.3, "M", 19, "Bacterial", "Yes", 33, "Yes",
  "1055", "POS", "Severe", 62.4, "M", 4, "Viral", "Yes", 31, "Yes",
  "1060", "POS", "Severe", 44.1, "M", 6, "Bacterial", "Yes", 5, "Yes",
  "1064", "POS", "Mild/Moderate", 44.0, "M", NA, "None", "No", 0, "Yes",
  "1066", "POS", "Mild/Moderate", 38.0, "M", NA, "None", "No", 0, "Yes",
  "1070", "POS", "Mild/Moderate", 40.0, "M", NA, "Bacterial", "No", 0, "Yes",
  "1071", "POS", "Mild/Moderate", 35.0, "F", NA, "None", "No", 0, "Yes",
  "1076", "POS", "Mild/Moderate", 46.0, "M", NA, "None", "No", 0, "Yes",
  "1080", "POS", "Mild/Moderate", 73.2, "M", NA, "None", "No", 0, "Yes",
  "1096", "POS", "Mild/Moderate", 48.0, "M", NA, "None", "No", 0, "Yes",
  "1098", "POS", "Mild/Moderate", 42.0, "M", NA, "None", "No", 0, "Yes",
  "1069", "POS", "Severe", 88.0, "F", NA, "none", "Yes", 0, "No",
  "1072", "POS", "Severe", 47.0, "M", NA, "none", "Yes", 21, "Yes",
  "1077", "POS", "Severe", 54.0, "M", NA, "Bacterial", "Yes", 8, "Yes",
  "1078", "POS", "Severe", 43.0, "M", NA, "none", "Yes", 15, "Yes",
  "1089", "POS", "Severe", 30.0, "F", NA, "Bacterial", "Yes", 16, "No",
  "1099", "POS", "Severe", 36.0, "M", NA, "none", "Yes", NA, "Yes"
)

# ------------------------------------------------------------------------------
# Read external tables and infer columns
# ------------------------------------------------------------------------------
sra_raw <- readr::read_csv(sra_run_table, show_col_types = FALSE)
names(sra_raw) <- normalize_names(names(sra_raw))

run_col <- pick_column(sra_raw,
  exact = c("run", "run_accession"),
  regex = "(^|_)run(_|$)",
  required = TRUE, label = "SRA run table"
)
exp_col <- pick_column(sra_raw,
  exact = c("experiment", "experiment_accession"),
  regex = "(^|_)experiment(_|$)",
  required = TRUE, label = "SRA run table"
)
biosample_col <- pick_column(sra_raw,
  exact = c("biosample", "biosample_accession", "bio_sample"),
  regex = "(^|_)biosample(_|$)|(^|_)bio_sample(_|$)",
  required = TRUE, label = "SRA run table"
)
geo_col <- pick_column(sra_raw,
  exact = c("geo_accession_exp", "geo_accession", "geo_accession_gsm"),
  regex = "geo.*accession",
  required = TRUE, label = "SRA run table"
)
title_col <- pick_column(sra_raw,
  exact = c("libraryname", "sample_name", "title"),
  regex = "libraryname|sample_name|title",
  required = FALSE
)

sra_tbl <- tibble(
  run_accession        = as.character(sra_raw[[run_col]]),
  experiment_accession = as.character(sra_raw[[exp_col]]),
  biosample_accession  = as.character(sra_raw[[biosample_col]]),
  geo_accession        = as.character(sra_raw[[geo_col]])
)
if (!is.na(title_col)) sra_tbl$run_title <- as.character(sra_raw[[title_col]]) else sra_tbl$run_title <- NA_character_

clinical_raw <- readr::read_csv(clinical_csv, show_col_types = FALSE)
names(clinical_raw) <- normalize_names(names(clinical_raw))

comet_col <- pick_column(clinical_raw,
  exact = c("master_record_id", "comet_id"),
  regex = "master_record_id|comet_id",
  required = FALSE
)
clinical_geo_col <- pick_column(clinical_raw,
  exact = c("geo_accession", "geo_accession_exp"),
  regex = "^geo_accession$|geo_accession_exp",
  required = FALSE
)
covid_col <- pick_column(clinical_raw,
  exact = c("sars_cov_2_status", "covid_status"),
  regex = "sars.*cov.*2.*status|covid.*status",
  required = FALSE
)
severity_col <- pick_column(clinical_raw,
  exact = c("disease_severity", "severity"),
  regex = "disease.*severity|^severity$",
  required = FALSE
)
age_col <- pick_column(clinical_raw, exact = c("age"), regex = "^age$", required = FALSE)
sex_col <- pick_column(clinical_raw,
  exact = c("gender", "sex"),
  regex = "gender|^sex$",
  required = FALSE
)
days_col <- pick_column(clinical_raw,
  exact = c("days_between_symptoms_and_sampling", "days_from_onset"),
  regex = "days.*symptoms.*sampling|days.*onset",
  required = FALSE
)
other_inf_col <- pick_column(clinical_raw,
  exact = c("other_infection_s", "other_infections"),
  regex = "other.*infection",
  required = FALSE
)
icu_col <- pick_column(clinical_raw,
  exact = c("icu_during_hospital_stay", "icu"),
  regex = "^icu$|icu.*hospital",
  required = FALSE
)
vent_days_col <- pick_column(clinical_raw,
  exact = c("days_under_mechanical_ventilation", "days_ventilated"),
  regex = "days.*mechanical.*ventilation|days.*ventilat",
  required = FALSE
)
discharged_col <- pick_column(clinical_raw,
  exact = c("discharged", "death"),
  regex = "^discharged$|death",
  required = FALSE
)
rbc_col <- pick_column(clinical_raw,
  exact = c("rbc_count"),
  regex = "^rbc(_count)?$|rbc_count",
  required = FALSE
)
hgb_col <- pick_column(clinical_raw,
  exact = c("hemoglobin", "hgb"),
  regex = "hemoglobin|^hgb$",
  required = FALSE
)
wbc_col <- pick_column(clinical_raw,
  exact = c("wbc", "wbc_count", "white_blood_cell_count"),
  regex = "^wbc$|wbc_count|white.*blood.*cell",
  required = FALSE
)
plt_col <- pick_column(clinical_raw,
  exact = c("platelet", "platelet_count", "plt", "plt_count"),
  regex = "platelet|^plt$|plt_count",
  required = FALSE
)

clinical_tbl_ext <- tibble(
  comet_id            = if (!is.na(comet_col)) as.character(clinical_raw[[comet_col]]) else NA_character_,
  geo_accession       = if (!is.na(clinical_geo_col)) as.character(clinical_raw[[clinical_geo_col]]) else NA_character_,
  covid_status        = if (!is.na(covid_col)) as.character(clinical_raw[[covid_col]]) else NA_character_,
  disease_severity    = if (!is.na(severity_col)) as.character(clinical_raw[[severity_col]]) else NA_character_,
  age                 = if (!is.na(age_col)) safe_numeric(clinical_raw[[age_col]]) else NA_real_,
  sex                 = if (!is.na(sex_col)) as.character(clinical_raw[[sex_col]]) else NA_character_,
  days_from_onset     = if (!is.na(days_col)) safe_numeric(clinical_raw[[days_col]]) else NA_real_,
  other_infection     = if (!is.na(other_inf_col)) as.character(clinical_raw[[other_inf_col]]) else NA_character_,
  icu_admission       = if (!is.na(icu_col)) as.character(clinical_raw[[icu_col]]) else NA_character_,
  days_ventilated     = if (!is.na(vent_days_col)) safe_numeric(clinical_raw[[vent_days_col]]) else NA_real_,
  discharged          = if (!is.na(discharged_col)) as.character(clinical_raw[[discharged_col]]) else NA_character_,
  rbc_count           = if (!is.na(rbc_col)) safe_numeric(clinical_raw[[rbc_col]]) else NA_real_,
  hgb                 = if (!is.na(hgb_col)) safe_numeric(clinical_raw[[hgb_col]]) else NA_real_,
  wbc_count           = if (!is.na(wbc_col)) safe_numeric(clinical_raw[[wbc_col]]) else NA_real_,
  plt_count           = if (!is.na(plt_col)) safe_numeric(clinical_raw[[plt_col]]) else NA_real_
)

# pad COMET IDs for joins
fix_comet_id <- function(x) {
  x <- as.character(x)
  x <- gsub("^ICC_", "", x)
  x <- trimws(x)
  x
}

clinical_tbl_ext <- clinical_tbl_ext %>% mutate(comet_id = fix_comet_id(comet_id))
paper_clinical_tbl <- paper_clinical_tbl %>% mutate(comet_id = fix_comet_id(comet_id))

# external clinical takes precedence for fields it actually contains;
# paper supplement fills gaps and provides a stable baseline
clinical_tbl <- paper_clinical_tbl %>%
  full_join(clinical_tbl_ext, by = "comet_id") %>%
  mutate(
    geo_accession        = coalesce(geo_accession, NA_character_),
    covid_status_final   = coalesce(covid_status, covid_status_paper),
    severity_final       = coalesce(disease_severity, severity_paper),
    age_final            = coalesce(age, age_paper),
    sex_final            = coalesce(sex, sex_paper),
    days_final           = coalesce(days_from_onset, days_from_onset_paper),
    other_infection_final= coalesce(other_infection, other_infection_paper),
    icu_final            = coalesce(icu_admission, icu_admission_paper),
    vent_days_final      = coalesce(days_ventilated, days_ventilated_paper),
    discharged_final     = coalesce(discharged, discharged_paper)
  )

# ------------------------------------------------------------------------------
# Build manifests
# ------------------------------------------------------------------------------
run_tbl <- sra_tbl %>%
  left_join(geo_title_tbl, by = "geo_accession") %>%
  mutate(
    pooled_flag = as.logical(pooled_flag),
    demux_status = ifelse(pooled_flag, "required", "not_needed")
  )

# attach clinical info via GEO when present
run_tbl <- run_tbl %>%
  left_join(
    clinical_tbl %>%
      select(comet_id, geo_accession, covid_status_final, severity_final,
             age_final, sex_final, days_final, other_infection_final,
             icu_final, vent_days_final, discharged_final,
             rbc_count, hgb, wbc_count, plt_count),
    by = "geo_accession"
  )

pooled_manifest <- run_tbl %>%
  filter(pooled_flag) %>%
  group_by(geo_accession, paper_sample, geo_title, biosample_accession, experiment_accession) %>%
  summarise(
    run_accessions = paste(sort(unique(run_accession)), collapse = ";"),
    pool_members = paste(sort(unique(na.omit(comet_id))), collapse = ";"),
    covid_status_members = paste(sort(unique(na.omit(covid_status_final))), collapse = ";"),
    severity_members = paste(sort(unique(na.omit(severity_final))), collapse = ";"),
    demux_required_for_subject_assignment = TRUE,
    .groups = "drop"
  )

write_csv(pooled_manifest, file.path(manifest_dir, "GSE163668_pooled_manifest.csv"))
run_tbl <- atlas_standardize_metadata(run_tbl, project_id = project_id)
write_csv(run_tbl, file.path(meta_dir, "GSE163668_run_metadata_full.csv"))
write_csv(sra_tbl, file.path(manifest_dir, "GSE163668_sra_table_normalized.csv"))
write_csv(clinical_tbl, file.path(manifest_dir, "GSE163668_clinical_join_table.csv"))

# ------------------------------------------------------------------------------
# Unpooled sample table (BioSample-level object units)
# ------------------------------------------------------------------------------
unpooled_samples <- run_tbl %>%
  filter(!pooled_flag) %>%
  group_by(geo_accession, paper_sample, geo_title, biosample_accession, experiment_accession) %>%
  summarise(
    run_accession = paste(sort(unique(run_accession)), collapse = ";"),
    comet_id = paste(sort(unique(na.omit(comet_id))), collapse = ";"),
    covid_status = paste(sort(unique(na.omit(covid_status_final))), collapse = ";"),
    severity = paste(sort(unique(na.omit(severity_final))), collapse = ";"),
    age = suppressWarnings(first(na.omit(age_final))),
    sex = paste(sort(unique(na.omit(sex_final))), collapse = ";"),
    days = suppressWarnings(first(na.omit(days_final))),
    other_infection = paste(sort(unique(na.omit(other_infection_final))), collapse = ";"),
    icu_admission = paste(sort(unique(na.omit(icu_final))), collapse = ";"),
    days_ventilated = suppressWarnings(first(na.omit(vent_days_final))),
    discharged = paste(sort(unique(na.omit(discharged_final))), collapse = ";"),
    rbc_count = suppressWarnings(first(na.omit(rbc_count))),
    hgb = suppressWarnings(first(na.omit(hgb))),
    wbc_count = suppressWarnings(first(na.omit(wbc_count))),
    plt_count = suppressWarnings(first(na.omit(plt_count))),
    .groups = "drop"
  ) %>%
  mutate(
    database_accession = project_id,
    sra_study = sra_study,
    bioproject = bioproject,
    organism = organism,
    sample_id = geo_accession,
    sample_alias = geo_title,
    patient_id = ifelse(grepl("^Patient", paper_sample), paper_sample, NA_character_),
    donor_id = ifelse(grepl("^Donor", paper_sample), paper_sample, NA_character_),
    study_group = case_when(
      grepl("^Donor", paper_sample) ~ "healthy_control",
      grepl("^Patient", paper_sample) & covid_status == "POS" ~ "covid19",
      grepl("^Patient", paper_sample) & covid_status == "NEG" ~ "non_covid_respiratory",
      TRUE ~ NA_character_
    ),
    factor = other_infection,
    cov19 = covid_status,
    is_sepsis = NA_character_,
    ards = case_when(
      severity == "Severe" ~ "possible",
      TRUE ~ NA_character_
    ),
    pneumonia = NA_character_,
    infection_site = "blood/systemic",
    pathogen_etiology = other_infection,
    mortality = case_when(
      toupper(discharged) %in% c("NO", "DEAD", "DIED") ~ "yes",
      toupper(discharged) %in% c("YES", "DISCHARGED") ~ "no",
      TRUE ~ NA_character_
    ),
    shock_status = NA_character_,
    timepoint_label = "first_day_after_admission",
    severityatday = severity,
    ventilation = case_when(
      !is.na(days_ventilated) & days_ventilated > 0 ~ "Y",
      !is.na(days_ventilated) & days_ventilated == 0 ~ "N",
      TRUE ~ NA_character_
    ),
    source_name = "Whole blood",
    tissue_label = "whole_blood",
    source_blood_fraction = "whole blood",
    blood_collection_tube = "EDTA",
    pbmc_isolation_method = NA_character_,
    ficoll = "no",
    rbc_lysis = "yes",
    sample_storage = "Fresh whole blood",
    frozen_or_fresh = "Fresh",
    platform = "10x Genomics",
    sequencer_model = "Illumina NovaSeq 6000",
    library_type = "10x 5prime GEX",
    chemistry_version = "Chromium Single Cell 5' v5.1",
    feature_barcoding = "none",
    vdj_capture = "no",
    read_len = "R1=28; R2=98",
    aligner = "Cell Ranger",
    cellranger_or_equivalent_version = "3.0.2",
    reference_genome = "GRCh38",
    gene_annotation = "Ensembl v85",
    seq_depth = NA_character_,
    n_loaded_cells = 15000L,
    paper_loaded_cells_per_sample = 15000L,
    paper_qc_min_genes_per_cell = 100L,
    paper_qc_max_percent_mt = 20,
    paper_qc_max_percent_ribo = 50,
    paper_qc_max_percent_hb = NA_real_,
    enrollment_window = "within 3 days of hospitalization; sampled on first day after admission",
    batch_unit_recommended = "library",
    pool_id = NA_character_,
    pool_members = NA_character_,
    pooled_flag = FALSE,
    demux_status = "not_needed",
    source_note = paste(
      "Paper-side experimental metadata embedded from GEO title + article methods + supplementary tables.",
      "Patient-level COMET/SRA linkage uses external CSV tables."
    )
  ) %>%
  ensure_columns(common_meta_columns)

unpooled_samples <- atlas_standardize_metadata(unpooled_samples, project_id = project_id)
write_csv(unpooled_samples, file.path(meta_dir, "GSE163668_sample_metadata_unpooled.csv"))

# ------------------------------------------------------------------------------
# Read BioSample-level matrices and build objects
# ------------------------------------------------------------------------------
input_dirs_all <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
input_dirs <- input_dirs_all[basename(input_dirs_all) %in% unpooled_samples$biosample_accession]
input_dirs <- input_dirs[order(match(basename(input_dirs), unpooled_samples$biosample_accession))]

if (nrow(unpooled_samples) == 0) {
  stop(
    "No unpooled sample metadata rows were generated.\n",
    "Check joins among SRA/GEO/clinical tables."
  )
}

if (length(input_dirs) == 0) {
  stop(
    "No BioSample directories under input_root matched unpooled metadata.\n",
    "Expected examples: ", paste(head(unpooled_samples$biosample_accession, 10), collapse = ", "), "\n",
    "Available under input_root: ", paste(head(basename(input_dirs_all), 20), collapse = ", ")
  )
}

message("Unpooled metadata rows: ", nrow(unpooled_samples))
message("Matched BioSample dirs : ", length(input_dirs))

unpooled_samples <- atlas_standardize_metadata(unpooled_samples, project_id = project_id)

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
  working_rds_path = character(),
  preprocess_status = character(),
  exclusion_reason = character(),
  n_called_after_droplet_qc = integer(),
  minimum_called_cells = integer(),
  validrops_strategy = character(),
  validrops_fallback_source = character(),
  droplet_qc_summary_path = character()
)

for (samn_dir in input_dirs) {
  biosample_id <- basename(samn_dir)
  raw_dir <- file.path(samn_dir, "outs", "raw_feature_bc_matrix")
  if (!dir.exists(raw_dir)) {
    message("Skipping missing raw_feature_bc_matrix: ", raw_dir)
    next
  }

  meta_row <- unpooled_samples %>% filter(biosample_accession == biosample_id)
  if (nrow(meta_row) != 1) {
    stop("Sample metadata row is not unique for ", biosample_id)
  }

  message("Reading ", biosample_id, " -> ", meta_row$geo_accession)
  # ATLAS_LOW_CELL_SKIP_V5: low-recovery libraries are recorded and skipped
  # before any Seurat RDS is created. Other preprocessing errors still stop.
  res <- tryCatch(
    read_gene_expression(raw_dir, sample_id = biosample_id),
    atlas_sample_exclusion = function(e) e
  )

  if (inherits(res, "atlas_sample_exclusion")) {
    message(
      "[EXCLUDED] ", biosample_id, " -> ", meta_row$geo_accession[[1]],
      ": ", res$reason,
      " (called=", res$n_called_after_droplet_qc,
      ", minimum=", res$minimum_called_cells, ")"
    )

    qc_summary <- bind_rows(
      qc_summary,
      tibble(
        biosample_accession = biosample_id,
        geo_accession = meta_row$geo_accession[[1]],
        paper_sample = meta_row$paper_sample[[1]],
        n_raw = as.integer(res$n_called_after_droplet_qc),
        n_after_qc = 0L,
        saved = FALSE,
        preqc_rds_path = NA_character_,
        qc_rds_path = NA_character_,
        working_rds_path = NA_character_,
        preprocess_status = "excluded_before_rds",
        exclusion_reason = res$reason,
        n_called_after_droplet_qc = as.integer(res$n_called_after_droplet_qc),
        minimum_called_cells = as.integer(res$minimum_called_cells),
        validrops_strategy = res$validrops_strategy,
        validrops_fallback_source = res$validrops_fallback_source,
        droplet_qc_summary_path = res$qc_summary_path
      )
    )
    next
  }

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
    percent.ribo <= paper_qc_max_percent_ribo
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

  obj <- run_optional_helpers(obj)

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

  # Working object
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

write_csv(qc_summary, file.path(qc_dir, "GSE163668_qc_summary.csv"))

excluded_samples <- qc_summary %>%
  dplyr::filter(preprocess_status %in% "excluded_before_rds")
write_csv(
  excluded_samples,
  file.path(qc_dir, "GSE163668_excluded_before_rds.csv")
)

sample_registry <- qc_summary %>% dplyr::filter(saved %in% TRUE)
write_csv(sample_registry, file.path(output_root, "GSE163668_preintegration_registry.csv"))
saveRDS(obj_list, file = file.path(output_root, "GSE163668_preintegration_objlist_unpooled.rds"))

if (nrow(sample_registry) == 0) {
  warning("No unpooled samples produced saved objects. Check qc/GSE163668_qc_summary.csv and metadata joins.")
}

message("Done. Saved pre-QC raw objects to   : ", rds_preqc_dir)
message("Done. Saved QC-filtered raw objects to: ", rds_qc_dir)
message("Done. Saved working objects to      : ", rds_work_dir)
message("Registry: ", file.path(output_root, "GSE163668_preintegration_registry.csv"))
