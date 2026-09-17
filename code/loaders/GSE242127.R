#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(Matrix)
  library(tidyr)
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

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) PROJECT_ROOT <- normalizePath(args[[1]], winslash = "/", mustWork = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}


project_id <- "GSE242127"
bioproject <- "PRJNA1011802"
organism <- "Homo sapiens"
expected_n <- 6L

normalize_names <- function(x) {
  x <- make.unique(x, sep = "_dup")
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "") %>%
    tolower()
}

pick_column <- function(df, exact = NULL, regex = NULL, required = FALSE, label = "table") {
  nms <- names(df)
  if (!is.null(exact)) {
    hit <- nms[nms %in% exact]
    if (length(hit) > 0) return(hit[[1]])
  }
  if (!is.null(regex)) {
    hit <- nms[stringr::str_detect(nms, regex)]
    if (length(hit) > 0) return(hit[[1]])
  }
  if (required) stop(label, " is missing required column; exact=", paste(exact, collapse=","), "; regex=", regex)
  NA_character_
}

col_or_na <- function(df, col, type = "character") {
  n <- nrow(df)
  if (is.na(col) || !col %in% names(df)) {
    if (type == "numeric") return(rep(NA_real_, n))
    if (type == "integer") return(rep(NA_integer_, n))
    return(rep(NA_character_, n))
  }
  x <- df[[col]]
  if (type == "numeric") return(suppressWarnings(as.numeric(x)))
  if (type == "integer") return(suppressWarnings(as.integer(x)))
  as.character(x)
}

read_optional_sra <- function(project_root, strict = TRUE) {
  candidates <- c(
    file.path(project_root, "configs", project_id, "GSE242127_SraRunTable.csv"),
    file.path(project_root, "configs", project_id, "GSE242127__SraRunTable.csv"),
    file.path(project_root, "configs", project_id, "SraRunTable.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NULL)

  sra_file <- normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
  sra_raw <- readr::read_csv(sra_file, show_col_types = FALSE, name_repair = "unique")
  names(sra_raw) <- normalize_names(names(sra_raw))

  bioproject_col <- pick_column(sra_raw, exact = c("bioproject"), regex = "^bioproject$", required = TRUE, label = "SRA table")
  biosample_col <- pick_column(sra_raw, exact = c("biosample"), regex = "^biosample$", required = TRUE, label = "SRA table")
  exp_col <- pick_column(sra_raw, exact = c("experiment"), regex = "^experiment$", required = FALSE)
  run_col <- pick_column(sra_raw, exact = c("run"), regex = "^run$", required = FALSE)
  sra_study_col <- pick_column(sra_raw, exact = c("sra_study"), regex = "sra.*study", required = FALSE)
  geo_col <- pick_column(sra_raw, exact = c("geo_accession_exp", "geo_accession"), regex = "geo.*accession", required = FALSE)
  sample_name_col <- pick_column(sra_raw, exact = c("sample_name"), regex = "^sample.*name$", required = FALSE)
  library_name_col <- pick_column(sra_raw, exact = c("libraryname", "library_name"), regex = "^library.*name$", required = FALSE)
  disease_state_col <- pick_column(sra_raw, exact = c("disease_state"), regex = "disease.*state|treatment", required = FALSE)
  assay_col <- pick_column(sra_raw, exact = c("assay_type"), regex = "assay.*type|library.*strategy", required = FALSE)
  instrument_col <- pick_column(sra_raw, exact = c("instrument", "model"), regex = "instrument|^model$", required = FALSE)
  layout_col <- pick_column(sra_raw, exact = c("librarylayout","library_layout"), regex = "library.*layout", required = FALSE)
  selection_col <- pick_column(sra_raw, exact = c("libraryselection","library_selection"), regex = "library.*selection", required = FALSE)
  source_col <- pick_column(sra_raw, exact = c("librarysource","library_source"), regex = "library.*source", required = FALSE)
  source_name_col <- pick_column(sra_raw, exact = c("source_name"), regex = "source.*name", required = FALSE)
  tissue_col <- pick_column(sra_raw, exact = c("tissue"), regex = "^tissue$", required = FALSE)
  cell_type_col <- pick_column(sra_raw, exact = c("cell_type"), regex = "cell.*type", required = FALSE)
  center_col <- pick_column(sra_raw, exact = c("centername", "center_name"), regex = "center.*name", required = FALSE)
  scientific_name_col <- pick_column(sra_raw, exact = c("scientificname", "scientific_name"), regex = "scientific.*name", required = FALSE)

  sra_norm <- tibble::tibble(
    run_accession = col_or_na(sra_raw, run_col),
    experiment_accession_sra = col_or_na(sra_raw, exp_col),
    biosample_accession = col_or_na(sra_raw, biosample_col),
    geo_accession_raw = col_or_na(sra_raw, geo_col),
    sample_name_sra = col_or_na(sra_raw, sample_name_col),
    library_name_sra = col_or_na(sra_raw, library_name_col),
    disease_state_sra = col_or_na(sra_raw, disease_state_col),
    bioproject_sra = col_or_na(sra_raw, bioproject_col),
    sra_study = col_or_na(sra_raw, sra_study_col),
    assay_type_sra = col_or_na(sra_raw, assay_col),
    instrument_sra = col_or_na(sra_raw, instrument_col),
    library_layout_sra = col_or_na(sra_raw, layout_col),
    library_selection_sra = col_or_na(sra_raw, selection_col),
    library_source_sra = col_or_na(sra_raw, source_col),
    source_name_sra = col_or_na(sra_raw, source_name_col),
    tissue_sra = col_or_na(sra_raw, tissue_col),
    cell_type_sra = col_or_na(sra_raw, cell_type_col),
    center_name_sra = col_or_na(sra_raw, center_col),
    scientific_name_sra = col_or_na(sra_raw, scientific_name_col)
  ) %>% dplyr::mutate(dplyr::across(dplyr::where(is.character), ~dplyr::na_if(.x, "")))

  if (!any(sra_norm$bioproject_sra == bioproject, na.rm = TRUE)) {
    msg <- paste0(
      "SRA RunTable does not contain expected BioProject ", bioproject,
      ". Observed BioProject values: ", paste(unique(stats::na.omit(sra_norm$bioproject_sra)), collapse = ";"),
      ". Replace configs/GSE242127/GSE242127_SraRunTable.csv with the correct GSE242127 table."
    )
    if (strict) stop(msg) else warning(msg)
  }

  sra_sample_agg <- sra_norm %>%
    dplyr::group_by(biosample_accession) %>%
    dplyr::summarise(
      run_accession = paste(unique(stats::na.omit(run_accession)), collapse = ";"),
      experiment_accession_sra = paste(unique(stats::na.omit(experiment_accession_sra)), collapse = ";"),
      geo_accession_raw = paste(unique(stats::na.omit(geo_accession_raw)), collapse = ";"),
      sample_name_sra = paste(unique(stats::na.omit(sample_name_sra)), collapse = ";"),
      library_name_sra = paste(unique(stats::na.omit(library_name_sra)), collapse = ";"),
      disease_state_sra = paste(unique(stats::na.omit(disease_state_sra)), collapse = ";"),
      bioproject_sra = paste(unique(stats::na.omit(bioproject_sra)), collapse = ";"),
      sra_study = paste(unique(stats::na.omit(sra_study)), collapse = ";"),
      assay_type_sra = paste(unique(stats::na.omit(assay_type_sra)), collapse = ";"),
      instrument_sra = paste(unique(stats::na.omit(instrument_sra)), collapse = ";"),
      library_layout_sra = paste(unique(stats::na.omit(library_layout_sra)), collapse = ";"),
      library_selection_sra = paste(unique(stats::na.omit(library_selection_sra)), collapse = ";"),
      library_source_sra = paste(unique(stats::na.omit(library_source_sra)), collapse = ";"),
      source_name_sra = paste(unique(stats::na.omit(source_name_sra)), collapse = ";"),
      tissue_sra = paste(unique(stats::na.omit(tissue_sra)), collapse = ";"),
      cell_type_sra = paste(unique(stats::na.omit(cell_type_sra)), collapse = ";"),
      center_name_sra = paste(unique(stats::na.omit(center_name_sra)), collapse = ";"),
      scientific_name_sra = paste(unique(stats::na.omit(scientific_name_sra)), collapse = ";"),
      n_sra_runs = dplyr::n_distinct(run_accession),
      .groups = "drop"
    ) %>% dplyr::mutate(dplyr::across(dplyr::where(is.character), ~dplyr::na_if(.x, "")))

  list(file = sra_file, raw = sra_norm, sample_agg = sra_sample_agg)
}

find_matrix_dir <- function(input_root, meta_row) {
  split_ids <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & x != ""]
    unique(unlist(strsplit(x, ";", fixed = TRUE)))
  }

  ids <- unique(c(
    split_ids(meta_row$biosample_accession),
    split_ids(meta_row$run_accession),
    split_ids(meta_row$experiment_accession),
    split_ids(meta_row$experiment_accession_sra),
    split_ids(meta_row$geo_accession),
    split_ids(meta_row$geo_title),
    split_ids(meta_row$geo_processed_sample_prefix),
    paste0(as.character(meta_row$geo_accession), "_", as.character(meta_row$geo_processed_sample_prefix))
  ))
  ids <- ids[!is.na(ids) & ids != ""]
  ids <- unique(trimws(ids))

  # Direct candidate search for common Cell Ranger layouts.
  for (id in ids) {
    candidates <- c(
      file.path(input_root, id, "outs", "raw_feature_bc_matrix"),
      file.path(input_root, id, "raw_feature_bc_matrix"),
      file.path(input_root, id)
    )
    for (p in candidates) {
      if (file.exists(file.path(p, "matrix.mtx.gz")) &&
          file.exists(file.path(p, "barcodes.tsv.gz")) &&
          file.exists(file.path(p, "features.tsv.gz"))) {
        return(normalizePath(p, winslash = "/", mustWork = TRUE))
      }
    }
  }

  # Recursive fallback: find all raw_feature_bc_matrix directories and map by
  # nearby directory names. This covers SRR-level Cell Ranger output roots such as
  # SRP457965/SRR25867779/outs/raw_feature_bc_matrix.
  all_dirs <- list.dirs(input_root, recursive = TRUE, full.names = TRUE)
  raw_dirs <- all_dirs[basename(all_dirs) == "raw_feature_bc_matrix"]
  if (length(raw_dirs) > 0) {
    for (p in raw_dirs) {
      if (!(file.exists(file.path(p, "matrix.mtx.gz")) &&
            file.exists(file.path(p, "barcodes.tsv.gz")) &&
            file.exists(file.path(p, "features.tsv.gz")))) next
      parts <- strsplit(normalizePath(p, winslash = "/", mustWork = FALSE), "/", fixed = TRUE)[[1]]
      if (any(parts %in% ids)) {
        return(normalizePath(p, winslash = "/", mustWork = TRUE))
      }
    }
  }

  NA_character_
}

resolve_input_root <- function(project_root, sra_obj = NULL) {
  sra_studies <- character()
  if (!is.null(sra_obj) && "sra_study" %in% names(sra_obj$sample_agg)) {
    sra_studies <- unique(stats::na.omit(sra_obj$sample_agg$sra_study))
    sra_studies <- sra_studies[nzchar(sra_studies)]
  }
  candidates <- unique(c(
    Sys.getenv("CELLRANGER_COUNT_ROOT", unset = ""),
    file.path(project_root, "cellranger_count", "output", project_id),
    file.path(project_root, "cellranger_count", "output", bioproject),
    file.path(project_root, "cellranger_count", "output", sra_studies),
    file.path(Sys.getenv("ATLAS_CELLRANGER_OUTPUT_ROOT", unset = ""), project_id),
    file.path(Sys.getenv("ATLAS_CELLRANGER_OUTPUT_ROOT", unset = ""), bioproject),
    file.path(Sys.getenv("ATLAS_CELLRANGER_OUTPUT_ROOT", unset = ""), project_id),
    file.path(Sys.getenv("ATLAS_CELLRANGER_OUTPUT_ROOT", unset = ""), bioproject)
  ))
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[dir.exists(candidates)]
  if (length(hit) == 0) {
    stop("Cannot find Cell Ranger output root. Tried:\n", paste(candidates, collapse = "\n"),
         "\nSet CELLRANGER_COUNT_ROOT to the directory containing the six sample subdirectories.")
  }
  normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

geo_sample_tbl <- tibble::tribble(
  ~geo_accession, ~geo_title, ~patient_id, ~days, ~timepoint_label, ~biosample_accession, ~experiment_accession, ~geo_processed_sample_prefix,
  "GSM7749422", "01-day1", "001", 1L, "D1_onset", "SAMN37223935", "SRX21589181", "S-001-1d_S4",
  "GSM7749423", "01-day7", "001", 7L, "D7_recovery", "SAMN37223934", "SRX21589182", "S-001-7d_S3",
  "GSM7749424", "02-day1", "002", 1L, "D1_onset", "SAMN37223933", "SRX21589183", "S-005-1d_S1",
  "GSM7749425", "02-day7", "002", 7L, "D7_recovery", "SAMN37223932", "SRX21589184", "S-005-7d_S2",
  "GSM7749426", "03-day1", "003", 1L, "D1_onset", "SAMN37223931", "SRX21589185", "A--S-014-1d_S13",
  "GSM7749427", "03-day7", "003", 7L, "D7_recovery", "SAMN37223930", "SRX21589186", "A--S-014-7d_S12"
)

clinical_patient_tbl <- tibble::tribble(
  ~patient_id, ~age, ~sex, ~is_sepsis, ~hypertensive_disease, ~diabetes_mellitus, ~cad_chf_mi, ~chronic_lung_disease, ~cerebral_infarction_or_mental_disorder, ~chronic_kidney_disease, ~hours_hospital_arrival_to_enrollment, ~hours_initial_antibiotic_to_enrollment, ~documented_temp_ge37, ~documented_sbp_lt90, ~vasopressor_therapy_within_48h, ~duration_mechanical_ventilation_hours, ~icu_stay_days, ~death_index_illness_or_hospitalization, ~second_hospital_admission_within_30d,
  "001", 76L, "female", "yes", "N", "Y", "N", "N", "Y", "N", 10L, 3.5, "Y", "Y", "N", NA_real_, 16L, "N", "N",
  "002", 68L, "male", "yes", "Y", "N", "N", "N", "Y", "N", 26L, 1.5, "Y", "Y", "Y", 156, 15L, "N", "N",
  "003", 76L, "male", "yes", "Y", "Y", "N", "Y", "N", "N", 12L, 3, "Y", "Y", "N", 72, 15L, "N", "Y"
)

clinical_day_tbl <- tibble::tribble(
  ~patient_id, ~days, ~sofa, ~oxygenation_index, ~ards_classification_table_s2_code, ~wbc_10e9_l, ~neutrophil_percent, ~lymphocyte_percent, ~monocyte_percent, ~crp_mg_l, ~procalcitonin_ng_ml, ~endotracheal_intubation, ~body_temperature_c, ~heart_rate_per_min,
  "001", 1L, 4L, 218L, 1L, 26.38, 88.2, 1.7, 8.9, 169.77, 19.02, "N", 37.2, 104L,
  "001", 7L, 0L, 336.36, 0L, 7.53, 82.6, 9.3, 7.0, 36.9, 0.15, "N", 36.5, 65L,
  "002", 1L, 11L, 127L, 2L, 11.02, 86.7, 5.4, 7.7, 136.16, 3.58, "Y", 36.2, 71L,
  "002", 7L, 3L, 185.45, 2L, 8.01, 79.3, 9.9, 8.7, 17.53, 0.07, "N", 36.0, 68L,
  "003", 1L, 3L, 110L, 2L, 17.79, 82.4, 8.9, 8.9, 14.85, 0.26, "Y", 37.4, 92L,
  "003", 7L, 1L, 265.41, 1L, 6.56, 72.6, 13.6, 13.6, 28.43, 0.24, "N", 36.9, 76L
)

ensure_columns <- function(df, wanted) {
  miss <- setdiff(wanted, names(df))
  for (m in miss) df[[m]] <- NA
  df[, unique(c(wanted, names(df))), drop = FALSE]
}

get_counts_matrix <- function(obj, assay = "RNA") {
  tryCatch(
    SeuratObject::GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) SeuratObject::GetAssayData(obj, assay = assay, slot = "counts")
  )
}

compute_percent_ribo <- function(obj) {
  ribo_genes <- grep("^RP[SL]", rownames(obj), value = TRUE)
  if (length(ribo_genes) == 0) return(rep(0, ncol(obj)))
  counts <- get_counts_matrix(obj)
  totals <- Matrix::colSums(counts)
  as.numeric(100 * Matrix::colSums(counts[ribo_genes, , drop = FALSE]) / pmax(totals, 1))
}

compute_percent_hb <- function(obj) {
  hb_genes <- intersect(rownames(obj), c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ"))
  if (length(hb_genes) == 0) return(rep(0, ncol(obj)))
  counts <- get_counts_matrix(obj)
  totals <- Matrix::colSums(counts)
  as.numeric(100 * Matrix::colSums(counts[hb_genes, , drop = FALSE]) / pmax(totals, 1))
}

read_gene_expression <- function(raw_dir, sample_id) {
  atlas_preprocess_10x_raw(
    raw_dir = raw_dir,
    sample_id = sample_id,
    qc_output_dir = file.path(qc_dir, "droplet_qc"),
    project_id = project_id
  )
}

attach_sample_metadata <- function(obj, meta_row) {
  atlas_attach_sample_metadata(obj, meta_row)
}

sra_obj <- read_optional_sra(PROJECT_ROOT, strict = TRUE)
input_root <- resolve_input_root(PROJECT_ROOT, sra_obj)

output_root <- file.path(PROJECT_ROOT, "pre_integration", project_id)
rds_preqc_dir <- file.path(output_root, "rds_preqc_raw")
rds_qc_dir <- file.path(output_root, "rds_qcfiltered_raw")
rds_work_dir <- file.path(output_root, "rds_working")
meta_dir <- file.path(output_root, "metadata")
qc_dir <- file.path(output_root, "qc")
manifest_dir <- file.path(output_root, "manifest")
for (d in c(output_root, rds_preqc_dir, rds_qc_dir, rds_work_dir, meta_dir, qc_dir, manifest_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

message("PROJECT_ROOT: ", PROJECT_ROOT)
message("input_root  : ", input_root)
message("SRA table   : ", sra_obj$file)
message("output_root : ", output_root)

apply_paper_qc_filter <- TRUE
min_cells_after_qc <- 200L

paper_methods_tbl <- tibble::tribble(
  ~field, ~value, ~note,
  "series_accession", "GSE242127", "Geriatric sepsis-induced ARDS PBMC scRNA-seq",
  "bioproject", "PRJNA1011802", "GEO relation",
  "cohort_scRNA_subjects", "3", "Table S2 and GEO design",
  "samples", "6", "Day 1 and Day 7 from three patients",
  "covid_status", "no SARS-CoV-2 infection", "Methods",
  "source_material", "PBMC", "GEO and methods",
  "scRNA_protocol", "10x Genomics Single Cell 3prime Library and Gel Bead Kit V3.1", "Methods/GEO",
  "sequencer", "Illumina NovaSeq 6000", "Methods/GEO",
  "read_mode", "150 bp paired-end", "Methods/GEO",
  "paper_original_cellranger_version", "5.0.0", "GEO/methods",
  "paper_original_seurat_version", "3.1.1", "GEO/methods",
  "paper_qc_min_genes_per_cell", "200", "GEO/methods",
  "paper_qc_max_genes_per_cell", "5000", "GEO/methods",
  "paper_qc_max_percent_mt", "15", "GEO/methods",
  "timepoint_definition", "Day 1 onset phase; Day 7 recovery phase", "Methods"
)

write_csv(geo_sample_tbl, file.path(manifest_dir, "GSE242127_geo_sample_embedded.csv"))
write_csv(clinical_patient_tbl, file.path(manifest_dir, "GSE242127_clinical_patient_tableS2_embedded.csv"))
write_csv(clinical_day_tbl, file.path(manifest_dir, "GSE242127_clinical_day_tableS2_embedded.csv"))
write_csv(paper_methods_tbl, file.path(manifest_dir, "GSE242127_paper_methods_and_results_embedded.csv"))
write_csv(sra_obj$raw, file.path(manifest_dir, "GSE242127_sra_table_normalized.csv"))
write_csv(sra_obj$sample_agg, file.path(manifest_dir, "GSE242127_sra_sample_aggregated.csv"))

sample_meta <- geo_sample_tbl %>%
  left_join(clinical_patient_tbl, by = "patient_id") %>%
  left_join(clinical_day_tbl, by = c("patient_id", "days")) %>%
  left_join(sra_obj$sample_agg, by = "biosample_accession") %>%
  mutate(
    database_accession = project_id,
    sra_study = ifelse(is.na(sra_study) | sra_study == "", "from_valid_GSE242127_SraRunTable", sra_study),
    bioproject = bioproject,
    organism = organism,
    experiment_accession = coalesce(experiment_accession_sra, experiment_accession),
    run_accession = coalesce(run_accession, NA_character_),
    sample_id = geo_accession,
    sample_alias = geo_title,
    paper_sample = geo_title,
    subject_id = patient_id,
    participant_id = patient_id,
    cellranger_dir_id = biosample_accession,

    study_group = "sepsis_ards",
    factor = "sepsis",
    cov19 = "no",
    is_sepsis = "yes",
    ards = "yes",
    pneumonia = NA_character_,
    infection_site = "not_reported",
    pathogen_etiology = "not_reported",
    mortality = ifelse(death_index_illness_or_hospitalization == "Y", "yes", "no"),
    final_outcome = ifelse(death_index_illness_or_hospitalization == "Y", "death_index_hospitalization", "survived_index_hospitalization"),
    mortality_90d = NA_character_,
    icu_admission = "yes",
    severity = paste0(study_group, "_", timepoint_label),
    severityatday = severity,

    ards_classification = paste0("tableS2_code_", ards_classification_table_s2_code),
    clinical_mapping_status = "exact_patient_and_day_from_GEO_title_and_TableS2",
    clinical_match_confidence = "exact_subject_day",

    qsofa = NA_real_,
    apache_ii = NA_real_,
    race = NA_character_,

    source_name = coalesce(source_name_sra, "blood"),
    tissue_label = coalesce(tissue_sra, "blood"),
    source_blood_fraction = "PBMC",
    cell_type = coalesce(cell_type_sra, "PBMC"),
    extracted_molecule = "total RNA",
    frozen_or_fresh = "flash frozen on dry ice / PBMC processed for 10x 3prime scRNA-seq",
    fresh_or_frozen = frozen_or_fresh,
    pbmc_isolation_method = "PBMCs isolated from 5 mL EDTA blood within two hours; scRNA-seq samples prepared with 10x Chromium",
    sample_storage = "flash frozen on dry ice according to GEO extraction protocol",
    cell_suspension_buffer = NA_character_,
    viability_reported = NA_character_,
    sorted_population = "PBMC",

    platform = "10x Genomics",
    sequencer_model = coalesce(instrument_sra, "Illumina NovaSeq 6000"),
    library_type = "10x Genomics single-cell gene expression",
    chemistry_version = "10x Genomics Single Cell 3prime V3.1",
    feature_barcoding = "none for RNA object",
    vdj_capture = "not_reported",
    read_len = "150 bp paired-end",
    assay_type = coalesce(assay_type_sra, "RNA-Seq"),
    library_layout = coalesce(library_layout_sra, "PAIRED"),
    library_selection = coalesce(library_selection_sra, "cDNA"),
    library_source = coalesce(library_source_sra, "TRANSCRIPTOMIC SINGLE CELL"),

    aligner = "Cell Ranger / user common atlas pipeline for this object; paper original used Cell Ranger v5.0.0",
    cellranger_or_equivalent_version = "atlas common reprocessing; paper provenance: Cell Ranger v5.0.0",
    paper_original_aligner = "10x Genomics Cell Ranger",
    paper_original_cellranger_version = "5.0.0",
    paper_original_reference_genome = "not provided/requested in GEO",
    reference_genome = "not_specified_in_GEO",
    gene_annotation = NA_character_,
    seq_depth = NA_character_,
    n_loaded_cells = NA_integer_,
    paper_qc_min_genes_per_cell = 200L,
    paper_qc_max_genes_per_cell = 5000L,
    paper_qc_max_percent_mt = 15,
    paper_seurat_version = "3.1.1",

    matrix_dir = vapply(seq_len(n()), function(i) find_matrix_dir(input_root, .[i, , drop = FALSE]), character(1)),
    enrollment_window = "Blood collected on day 1 and day 7 of admission from geriatric sepsis-induced ARDS patients",
    batch_unit_recommended = "sample",
    source_note = paste(
      "GSE242127 contains PBMC scRNA-seq from three geriatric sepsis-induced ARDS patients at day 1 and day 7.",
      "GEO title maps directly to Table S2 patient/day identifiers: 01/02/03 to 001/002/003 and day1/day7 to D1/D7.",
      "SARS-CoV-2 infection was absent in the enrolled sepsis-induced ARDS patients.",
      "Clinical values are embedded from Supplementary Table S2 at patient-level and day-specific levels."
    ),
    read_mode = NA_character_,
    predicted.cluster = NA_character_
  )

if (nrow(sample_meta) != expected_n) stop("Expected six sample metadata rows, got ", nrow(sample_meta))
if (any(is.na(sample_meta$matrix_dir))) {
  print(sample_meta %>% filter(is.na(matrix_dir)) %>% select(geo_accession, geo_title, biosample_accession, matrix_dir))
  stop("Some samples do not have a resolvable matrix directory.")
}
if (any(is.na(sample_meta$age) | is.na(sample_meta$sex) | is.na(sample_meta$sofa) | is.na(sample_meta$crp_mg_l))) {
  print(sample_meta %>% filter(is.na(age) | is.na(sex) | is.na(sofa) | is.na(crp_mg_l)))
  stop("Some samples did not receive Table S2 clinical metadata.")
}

common_meta_columns <- c(
  "database_accession", "sra_study", "bioproject", "organism",
  "run_accession", "experiment_accession", "biosample_accession", "geo_accession",
  "cellranger_dir_id", "sample_id", "sample_alias", "paper_sample", "geo_title",
  "patient_id", "subject_id", "participant_id", "clinical_mapping_status", "clinical_match_confidence",
  "study_group", "factor", "cov19", "is_sepsis", "ards", "pneumonia",
  "infection_site", "pathogen_etiology", "mortality", "final_outcome", "mortality_90d",
  "icu_admission", "days", "timepoint_label", "severity", "severityatday",
  "sex", "age", "race", "qsofa", "apache_ii", "sofa", "oxygenation_index",
  "ards_classification_table_s2_code", "ards_classification",
  "wbc_10e9_l", "neutrophil_percent", "lymphocyte_percent", "monocyte_percent",
  "crp_mg_l", "procalcitonin_ng_ml", "endotracheal_intubation", "body_temperature_c", "heart_rate_per_min",
  "hypertensive_disease", "diabetes_mellitus", "cad_chf_mi", "chronic_lung_disease",
  "cerebral_infarction_or_mental_disorder", "chronic_kidney_disease",
  "hours_hospital_arrival_to_enrollment", "hours_initial_antibiotic_to_enrollment",
  "documented_temp_ge37", "documented_sbp_lt90", "vasopressor_therapy_within_48h",
  "duration_mechanical_ventilation_hours", "icu_stay_days",
  "death_index_illness_or_hospitalization", "second_hospital_admission_within_30d",
  "source_name", "tissue_label", "source_blood_fraction", "cell_type", "extracted_molecule",
  "frozen_or_fresh", "fresh_or_frozen", "pbmc_isolation_method", "sample_storage",
  "cell_suspension_buffer", "viability_reported", "sorted_population",
  "platform", "sequencer_model", "library_type", "chemistry_version",
  "feature_barcoding", "vdj_capture", "read_len", "assay_type",
  "library_layout", "library_selection", "library_source",
  "aligner", "cellranger_or_equivalent_version", "paper_original_aligner",
  "paper_original_cellranger_version", "paper_original_reference_genome", "reference_genome",
  "gene_annotation", "seq_depth", "n_loaded_cells",
  "paper_qc_min_genes_per_cell", "paper_qc_max_genes_per_cell", "paper_qc_max_percent_mt",
  "paper_seurat_version", "n_sra_runs", "sample_name_sra", "disease_state_sra",
  "geo_processed_sample_prefix", "matrix_dir",
  "enrollment_window", "batch_unit_recommended", "source_note", "read_mode", "predicted.cluster"
)
sample_meta <- ensure_columns(sample_meta, common_meta_columns)

sample_meta <- atlas_standardize_metadata(sample_meta, project_id = project_id)

obj_list <- list()
qc_summary <- tibble(
  cellranger_dir_id = character(),
  run_accession = character(),
  biosample_accession = character(),
  geo_accession = character(),
  paper_sample = character(),
  patient_id = character(),
  days = numeric(),
  timepoint_label = character(),
  study_group = character(),
  factor = character(),
  n_raw = integer(),
  n_after_qc = integer(),
  saved = logical(),
  preqc_rds_path = character(),
  qc_rds_path = character(),
  working_rds_path = character()
)

for (i in seq_len(nrow(sample_meta))) {
  meta_row <- sample_meta[i, , drop = FALSE]
  raw_dir <- meta_row$matrix_dir
  dir_id <- meta_row$cellranger_dir_id

  message("Reading ", dir_id, " -> ", meta_row$geo_accession, " (", meta_row$geo_title, ")")
  res <- read_gene_expression(raw_dir, sample_id = dir_id)

  obj <- CreateSeuratObject(
    counts = res$counts,
    project = meta_row$geo_accession,
    meta.data = res$meta,
    min.cells = 0,
    min.features = 0
  )
  obj <- atlas_add_droplet_qc_assays(obj, res)

  n_raw <- ncol(obj)
  sample_meta$read_mode[i] <- res$read_mode
  sample_meta$n_loaded_cells[i] <- n_raw
  meta_row$read_mode <- res$read_mode
  meta_row$n_loaded_cells <- n_raw

  obj <- attach_sample_metadata(obj, meta_row)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- compute_percent_ribo(obj)
  obj[["percent.hb"]] <- compute_percent_hb(obj)
  obj[["paper_qc_pass"]] <- with(obj@meta.data,
    nFeature_RNA >= paper_qc_min_genes_per_cell &
      nFeature_RNA <= paper_qc_max_genes_per_cell &
      percent.mt <= paper_qc_max_percent_mt
  )
  obj[["atlas_basic_qc_pass"]] <- with(obj@meta.data, nFeature_RNA >= 200)
  obj <- atlas_standardize_seurat_metadata(obj, project_id = project_id)

  preqc_rds <- file.path(rds_preqc_dir, paste0(meta_row$geo_accession, "__", dir_id, "__preqc_raw.rds"))
  saveRDS(obj, preqc_rds, compress = TRUE)

  qc_subset <- rep(TRUE, ncol(obj))
  if (apply_paper_qc_filter) qc_subset <- qc_subset & obj$paper_qc_pass
  qc_subset <- qc_subset & obj$atlas_basic_qc_pass
  cells_keep <- colnames(obj)[which(qc_subset)]
  n_after_qc <- length(cells_keep)

  if (n_after_qc < min_cells_after_qc) {
    warning(meta_row$geo_accession, " dropped because only ", n_after_qc, " cells remained after QC.")
    qc_summary <- bind_rows(qc_summary, tibble(
      cellranger_dir_id = dir_id, run_accession = meta_row$run_accession,
      biosample_accession = meta_row$biosample_accession, geo_accession = meta_row$geo_accession,
      paper_sample = meta_row$paper_sample, patient_id = meta_row$patient_id,
      days = meta_row$days, timepoint_label = meta_row$timepoint_label,
      study_group = meta_row$study_group, factor = meta_row$factor,
      n_raw = n_raw, n_after_qc = n_after_qc, saved = FALSE,
      preqc_rds_path = preqc_rds, qc_rds_path = NA_character_, working_rds_path = NA_character_
    ))
    next
  }

  obj_qc <- subset(obj, cells = cells_keep)
  qc_rds <- file.path(rds_qc_dir, paste0(meta_row$geo_accession, "__", dir_id, "__qcfiltered_raw.rds"))
  saveRDS(obj_qc, qc_rds, compress = TRUE)

  obj_work <- NormalizeData(obj_qc, verbose = FALSE)
  obj_work <- FindVariableFeatures(obj_work, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  work_rds <- file.path(rds_work_dir, paste0(meta_row$geo_accession, "__", dir_id, "__working.rds"))
  saveRDS(obj_work, work_rds, compress = TRUE)

  obj_list[[as.character(meta_row$geo_accession)]] <- obj_work

  qc_summary <- bind_rows(qc_summary, tibble(
    cellranger_dir_id = dir_id, run_accession = meta_row$run_accession,
    biosample_accession = meta_row$biosample_accession, geo_accession = meta_row$geo_accession,
    paper_sample = meta_row$paper_sample, patient_id = meta_row$patient_id,
    days = meta_row$days, timepoint_label = meta_row$timepoint_label,
    study_group = meta_row$study_group, factor = meta_row$factor,
    n_raw = n_raw, n_after_qc = n_after_qc, saved = TRUE,
    preqc_rds_path = preqc_rds, qc_rds_path = qc_rds, working_rds_path = work_rds
  ))
}

sample_meta <- atlas_standardize_metadata(sample_meta, project_id = project_id)
write_csv(sample_meta, file.path(meta_dir, "GSE242127_sample_metadata.csv"))
write_csv(sample_meta, file.path(meta_dir, "GSE242127_run_metadata.csv"))
write_csv(qc_summary, file.path(qc_dir, "GSE242127_qc_summary.csv"))

registry <- qc_summary %>%
  filter(saved) %>%
  left_join(sample_meta %>%
              select(cellranger_dir_id, run_accession, biosample_accession, geo_accession,
                     paper_sample, patient_id, days, timepoint_label, study_group, factor, is_sepsis,
                     ards, cov19, sex, age, sofa, oxygenation_index, crp_mg_l,
                     clinical_mapping_status, batch_unit_recommended),
            by = c("cellranger_dir_id", "run_accession", "biosample_accession", "geo_accession",
                   "paper_sample", "patient_id", "days", "timepoint_label", "study_group", "factor")) %>%
  mutate(project_id = project_id, preintegration_tier = "working", object_grain = "sample")

write_csv(registry, file.path(output_root, "GSE242127_preintegration_registry.csv"))
saveRDS(obj_list, file.path(output_root, "GSE242127_preintegration_objlist.rds"))

message("Done. Saved pre-QC raw objects to    : ", rds_preqc_dir)
message("Done. Saved QC-filtered raw objects to: ", rds_qc_dir)
message("Done. Saved working objects to       : ", rds_work_dir)
message("Registry: ", file.path(output_root, "GSE242127_preintegration_registry.csv"))
