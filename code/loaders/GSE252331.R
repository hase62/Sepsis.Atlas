#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
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

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) PROJECT_ROOT <- normalizePath(args[[1]], winslash = "/", mustWork = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

geo_gex_tbl <- tibble::tribble(
  ~geo_accession, ~geo_title, ~raw_internal_sample_id, ~participant_id_base, ~subject_label, ~factor, ~study_group, ~timepoint_label, ~days, ~rap_cci_classification, ~severity, ~paired_adt_geo_accession,
  "GSM7999550", "Sepsis 5, CCI, Day 14, GEX", "P0657", "P0657", "Sepsis 5", "sepsis", "post_sepsis_CCI_D14", "D14_post_sepsis_CCI", 14L, "CCI", "chronic_critical_illness", "GSM7999567",
  "GSM7999551", "Sepsis 6, RAP, Day 14, GEX", "P0658", "P0658", "Sepsis 6", "sepsis", "post_sepsis_RAP_D14", "D14_post_sepsis_RAP", 14L, "RAP", "rapid_recovery", "GSM7999568",
  "GSM7999552", "Sepsis 7, CCI, Day 14, GEX", "P0660", "P0660", "Sepsis 7", "sepsis", "post_sepsis_CCI_D14", "D14_post_sepsis_CCI", 14L, "CCI", "chronic_critical_illness", "GSM7999569",
  "GSM7999553", "Sepsis 8, RAP, Day 14, GEX", "P0661", "P0661", "Sepsis 8", "sepsis", "post_sepsis_RAP_D14", "D14_post_sepsis_RAP", 14L, "RAP", "rapid_recovery", "GSM7999570",
  "GSM7999554", "Sepsis 9, RAP, Day 14, GEX", "P0663", "P0663", "Sepsis 9", "sepsis", "post_sepsis_RAP_D14", "D14_post_sepsis_RAP", 14L, "RAP", "rapid_recovery", "GSM7999571",
  "GSM7999555", "Sepsis 10, Day 4, GEX", "P0664", "P0664", "Sepsis 10", "sepsis", "acute_sepsis_D4", "D4_acute_sepsis", 4L, "acute_sepsis", "acute_sepsis", "GSM7999572",
  "GSM7999556", "Sepsis 11, Day 4, GEX", "P0665", "P0665", "Sepsis 11", "sepsis", "acute_sepsis_D4", "D4_acute_sepsis", 4L, "acute_sepsis", "acute_sepsis", "GSM7999573",
  "GSM7999557", "Sepsis 12, RAP, Day 14, GEX", "P0665-T5", "P0665", "Sepsis 12", "sepsis", "post_sepsis_RAP_D14", "D14_post_sepsis_RAP", 14L, "RAP", "rapid_recovery", "GSM7999574",
  "GSM7999558", "Sepsis 13, Day 4, GEX", "P0667", "P0667", "Sepsis 13", "sepsis", "acute_sepsis_D4", "D4_acute_sepsis", 4L, "acute_sepsis", "acute_sepsis", "GSM7999575",
  "GSM7999559", "Sepsis 14, Day 4, GEX", "P0668", "P0668", "Sepsis 14", "sepsis", "acute_sepsis_D4", "D4_acute_sepsis", 4L, "acute_sepsis", "acute_sepsis", "GSM7999576",
  "GSM7999560", "Healthy Control 6, GEX", "PC003", "PC003", "Healthy Control 6", "healthy_control", "healthy_control", "healthy_control", 0L, "healthy_control", "healthy_control", "GSM7999577",
  "GSM7999561", "Healthy Control 7, GEX", "PC009", "PC009", "Healthy Control 7", "healthy_control", "healthy_control", "healthy_control", 0L, "healthy_control", "healthy_control", "GSM7999578",
  "GSM7999562", "Healthy Control 8, GEX", "PC011", "PC011", "Healthy Control 8", "healthy_control", "healthy_control", "healthy_control", 0L, "healthy_control", "healthy_control", "GSM7999579",
  "GSM7999563", "Healthy Control 9, GEX", "PC012", "PC012", "Healthy Control 9", "healthy_control", "healthy_control", "healthy_control", 0L, "healthy_control", "healthy_control", "GSM7999580",
  "GSM7999564", "Healthy Control 10, GEX", "PC032", "PC032", "Healthy Control 10", "healthy_control", "healthy_control", "healthy_control", 0L, "healthy_control", "healthy_control", "GSM7999581",
  "GSM7999565", "Healthy Control 11, GEX", "PC034", "PC034", "Healthy Control 11", "healthy_control", "healthy_control", "healthy_control", 0L, "healthy_control", "healthy_control", "GSM7999582",
  "GSM7999566", "Healthy Control 12, GEX", "PC057", "PC057", "Healthy Control 12", "healthy_control", "healthy_control", "healthy_control", 0L, "healthy_control", "healthy_control", "GSM7999583"
)

geo_adt_tbl <- tibble::tribble(
  ~adt_geo_accession, ~adt_geo_title, ~adt_internal_sample_id, ~adt_participant_id_base, ~paired_gex_geo_accession,
  "GSM7999567", "Sepsis 5, CCI, Day 14, ADT", "P0657", "P0657", "GSM7999550",
  "GSM7999568", "Sepsis 6, RAP, Day 14, ADT", "P0658", "P0658", "GSM7999551",
  "GSM7999569", "Sepsis 7, CCI, Day 14, ADT", "P0660", "P0660", "GSM7999552",
  "GSM7999570", "Sepsis 8, RAP, Day 14, ADT", "P0661", "P0661", "GSM7999553",
  "GSM7999571", "Sepsis 9, RAP, Day 14, ADT", "P0663", "P0663", "GSM7999554",
  "GSM7999572", "Sepsis 10, Day 4, ADT", "P0664", "P0664", "GSM7999555",
  "GSM7999573", "Sepsis 11, Day 4, ADT", "P0665", "P0665", "GSM7999556",
  "GSM7999574", "Sepsis 12, RAP, Day 14, ADT", "P0665-T5", "P0665", "GSM7999557",
  "GSM7999575", "Sepsis 13, Day 4, ADT", "P0667", "P0667", "GSM7999558",
  "GSM7999576", "Sepsis 14, Day 4, ADT", "P0668", "P0668", "GSM7999559",
  "GSM7999577", "Healthy Control 6, ADT", "PC003", "PC003", "GSM7999560",
  "GSM7999578", "Healthy Control 7, ADT", "PC009", "PC009", "GSM7999561",
  "GSM7999579", "Healthy Control 8, ADT", "PC011", "PC011", "GSM7999562",
  "GSM7999580", "Healthy Control 9, ADT", "PC012", "PC012", "GSM7999563",
  "GSM7999581", "Healthy Control 10, ADT", "PC032", "PC032", "GSM7999564",
  "GSM7999582", "Healthy Control 11, ADT", "PC034", "PC034", "GSM7999565",
  "GSM7999583", "Healthy Control 12, ADT", "PC057", "PC057", "GSM7999566"
)

adt_panel_tbl <- tibble::tribble(
  ~adt_reagent_id, ~adt_marker, ~adt_oligo_sequence,
  "Biolegend_307663_C0159", "HLADR", "AATAGCGAGCAAGTA",
  "Biolegend_300479_C0034", "CD3", "CTCATTGTAACTCCT",
  "Biolegend_323053_C0392", "CD15", "TCACCAGTACCTAGT",
  "Biolegend_351356_C0390", "CD127", "GTGTGTTGTCCTATG",
  "Biolegend_353747_C0140", "CD183", "GCGATGGTAGATTAT",
  "Biolegend_392909_C0166", "CD66B", "AGCTGTAAGTTTCGG",
  "Biolegend_300567_C0072", "CD4", "TGTTCCCGCTCAACT",
  "Biolegend_301859_C0081", "CD14", "TCTCAGACCTCCGTA",
  "Biolegend_353440_C0143", "CD196", "GATCCCTTTGTCACT",
  "Biolegend_302649_C0085", "CD25", "TTTGTCCTGTACGCC",
  "Biolegend_301359_C0161", "CD11B", "GACAAGTGATCTGCA"
)


project_id <- "GSE252331"
sra_study_expected <- "SRP481033"
bioproject <- "PRJNA1060229"
organism <- "Homo sapiens"
expected_n_gex <- 17L
expected_n_adt <- 17L

fixed_cellranger_root <- file.path(Sys.getenv("ATLAS_CELLRANGER_OUTPUT_ROOT", unset = ""), "SRP481033_R1R2_swapfix_SC5P-R2")
fixed_cellranger_reprocessing_id <- "SRP481033_R1R2_swapfix_SC5P-R2"
fixed_cellranger_chemistry_parameter <- "SC5P-R2"

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

list_existing_dirs <- function(paths) {
  paths <- unique(paths[nzchar(paths)])
  paths[dir.exists(paths)]
}

read_sra_table <- function(project_root) {
  candidates <- c(
    file.path(project_root, "configs", project_id, "GSE252331_SraRunTable.csv"),
    file.path(project_root, "configs", project_id, "SraRunTable.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("Cannot find GSE252331 SraRunTable. Put it at configs/GSE252331/GSE252331_SraRunTable.csv")
  }

  sra_file <- normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
  sra_raw <- readr::read_csv(sra_file, show_col_types = FALSE, name_repair = "unique")
  names(sra_raw) <- normalize_names(names(sra_raw))

  run_col <- pick_column(sra_raw, exact = c("run"), regex = "^run$", required = TRUE, label = "SRA table")
  bioproject_col <- pick_column(sra_raw, exact = c("bioproject"), regex = "^bioproject$", required = TRUE, label = "SRA table")
  biosample_col <- pick_column(sra_raw, exact = c("biosample"), regex = "^biosample$", required = TRUE, label = "SRA table")
  exp_col <- pick_column(sra_raw, exact = c("experiment"), regex = "^experiment$", required = TRUE, label = "SRA table")
  sra_study_col <- pick_column(sra_raw, exact = c("sra_study"), regex = "sra.*study", required = TRUE, label = "SRA table")
  library_name_col <- pick_column(sra_raw, exact = c("library_name", "libraryname"), regex = "^library.*name$", required = TRUE, label = "SRA table")
  sample_name_col <- pick_column(sra_raw, exact = c("sample_name", "samplename"), regex = "^sample.*name$", required = FALSE)
  cci_col <- pick_column(sra_raw, exact = c("cci"), regex = "^cci$", required = FALSE)
  timepoint_col <- pick_column(sra_raw, exact = c("timepoint"), regex = "^timepoint$", required = FALSE)
  cell_type_col <- pick_column(sra_raw, exact = c("cell_type"), regex = "cell.*type", required = FALSE)
  tissue_col <- pick_column(sra_raw, exact = c("tissue"), regex = "^tissue$", required = FALSE)
  source_name_col <- pick_column(sra_raw, exact = c("source_name"), regex = "source.*name", required = FALSE)
  instrument_col <- pick_column(sra_raw, exact = c("instrument"), regex = "instrument|model", required = FALSE)
  assay_col <- pick_column(sra_raw, exact = c("assay_type"), regex = "assay.*type|library.*strategy", required = FALSE)
  layout_col <- pick_column(sra_raw, exact = c("librarylayout","library_layout"), regex = "library.*layout", required = FALSE)
  selection_col <- pick_column(sra_raw, exact = c("libraryselection","library_selection"), regex = "library.*selection", required = FALSE)
  source_col <- pick_column(sra_raw, exact = c("librarysource","library_source"), regex = "library.*source", required = FALSE)
  center_col <- pick_column(sra_raw, exact = c("center_name"), regex = "center.*name", required = FALSE)
  datastore_col <- pick_column(sra_raw, exact = c("datastore_filetype"), regex = "datastore.*filetype", required = FALSE)
  release_col <- pick_column(sra_raw, exact = c("releasedate"), regex = "release", required = FALSE)
  avgspot_col <- pick_column(sra_raw, exact = c("avgspotlen"), regex = "avg.*spot", required = FALSE)
  bases_col <- pick_column(sra_raw, exact = c("bases"), regex = "^bases$", required = FALSE)
  bytes_col <- pick_column(sra_raw, exact = c("bytes"), regex = "^bytes$", required = FALSE)

  sra_norm <- tibble::tibble(
    run_accession = col_or_na(sra_raw, run_col),
    experiment_accession = col_or_na(sra_raw, exp_col),
    biosample_accession = col_or_na(sra_raw, biosample_col),
    bioproject_sra = col_or_na(sra_raw, bioproject_col),
    sra_study = col_or_na(sra_raw, sra_study_col),
    library_name_sra = col_or_na(sra_raw, library_name_col),
    sample_name_sra = col_or_na(sra_raw, sample_name_col),
    cci_sra = col_or_na(sra_raw, cci_col),
    timepoint_sra = col_or_na(sra_raw, timepoint_col),
    cell_type_sra = col_or_na(sra_raw, cell_type_col),
    tissue_sra = col_or_na(sra_raw, tissue_col),
    source_name_sra = col_or_na(sra_raw, source_name_col),
    instrument_sra = col_or_na(sra_raw, instrument_col),
    assay_type_sra = col_or_na(sra_raw, assay_col),
    library_layout_sra = col_or_na(sra_raw, layout_col),
    library_selection_sra = col_or_na(sra_raw, selection_col),
    library_source_sra = col_or_na(sra_raw, source_col),
    center_name_sra = col_or_na(sra_raw, center_col),
    datastore_filetype_sra = col_or_na(sra_raw, datastore_col),
    release_date_sra = col_or_na(sra_raw, release_col),
    avg_spot_len_sra = suppressWarnings(as.numeric(col_or_na(sra_raw, avgspot_col))),
    bases_sra = suppressWarnings(as.numeric(col_or_na(sra_raw, bases_col))),
    bytes_sra = suppressWarnings(as.numeric(col_or_na(sra_raw, bytes_col)))
  ) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.character), ~dplyr::na_if(.x, "")))

  if (!any(sra_norm$bioproject_sra == bioproject, na.rm = TRUE)) {
    stop("SRA RunTable does not contain expected BioProject ", bioproject,
         ". Observed: ", paste(unique(stats::na.omit(sra_norm$bioproject_sra)), collapse = ";"))
  }
  if (!any(sra_norm$sra_study == sra_study_expected, na.rm = TRUE)) {
    stop("SRA RunTable does not contain expected SRA Study ", sra_study_expected,
         ". Observed: ", paste(unique(stats::na.omit(sra_norm$sra_study)), collapse = ";"))
  }

  list(file = sra_file, raw = sra_norm)
}

cellranger_roots <- function(project_root) {
  # Fixed validated input: R1/R2-corrected Cell Ranger reprocessing using SC5P-R2.
  # Do not fall back to the original SRP481033 output because that output was generated
  # with the SRA reads in the wrong Cell Ranger roles and contained no usable cells.
  root <- fixed_cellranger_root
  if (!dir.exists(root)) {
    stop("Validated GSE252331 Cell Ranger root is not accessible: ", root)
  }
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

has_10x_triplet <- function(p) {
  file.exists(file.path(p, "matrix.mtx.gz")) &&
    file.exists(file.path(p, "barcodes.tsv.gz")) &&
    file.exists(file.path(p, "features.tsv.gz"))
}

find_matrix_dir <- function(roots, meta_row) {
  ids <- unique(c(
    meta_row$run_accession,
    meta_row$biosample_accession,
    meta_row$experiment_accession,
    meta_row$geo_accession,
    meta_row$sample_id,
    meta_row$raw_internal_sample_id
  ))
  ids <- ids[!is.na(ids) & ids != ""]
  ids <- unique(unlist(strsplit(paste(ids, collapse = ";"), ";", fixed = TRUE)))
  ids <- trimws(ids)

  for (root in roots) {
    for (id in ids) {
      candidates <- c(
        file.path(root, id, "outs", "raw_feature_bc_matrix"),
        file.path(root, id, "raw_feature_bc_matrix"),
        file.path(root, id)
      )
      for (p in candidates) {
        if (has_10x_triplet(p)) return(normalizePath(p, winslash = "/", mustWork = TRUE))
      }
    }

    raw_dirs <- list.dirs(root, recursive = TRUE, full.names = TRUE)
    raw_dirs <- raw_dirs[basename(raw_dirs) == "raw_feature_bc_matrix"]
    for (p in raw_dirs) {
      if (!has_10x_triplet(p)) next
      parts <- strsplit(normalizePath(p, winslash = "/", mustWork = FALSE), "/", fixed = TRUE)[[1]]
      if (any(parts %in% ids)) return(normalizePath(p, winslash = "/", mustWork = TRUE))
    }
  }

  NA_character_
}

build_sample_metadata <- function(project_root) {
  sra_obj <- read_sra_table(project_root)
  roots <- cellranger_roots(project_root)
  if (length(roots) == 0) {
    stop("Validated fixed Cell Ranger root for GSE252331 is unavailable.")
  }

  adt_sra <- sra_obj$raw %>%
    dplyr::rename(
      adt_run_accession = run_accession,
      adt_experiment_accession = experiment_accession,
      adt_biosample_accession = biosample_accession,
      adt_bioproject_sra = bioproject_sra,
      adt_sra_study = sra_study,
      adt_library_name_sra = library_name_sra,
      adt_sample_name_sra = sample_name_sra,
      adt_cci_sra = cci_sra,
      adt_timepoint_sra = timepoint_sra,
      adt_cell_type_sra = cell_type_sra,
      adt_tissue_sra = tissue_sra,
      adt_source_name_sra = source_name_sra,
      adt_instrument_sra = instrument_sra,
      adt_assay_type_sra = assay_type_sra,
      adt_library_layout_sra = library_layout_sra,
      adt_library_selection_sra = library_selection_sra,
      adt_library_source_sra = library_source_sra,
      adt_center_name_sra = center_name_sra,
      adt_datastore_filetype_sra = datastore_filetype_sra,
      adt_release_date_sra = release_date_sra,
      adt_avg_spot_len_sra = avg_spot_len_sra,
      adt_bases_sra = bases_sra,
      adt_bytes_sra = bytes_sra
    )

  sample_meta <- geo_gex_tbl %>%
    dplyr::left_join(geo_adt_tbl, by = c("paired_adt_geo_accession" = "adt_geo_accession")) %>%
    dplyr::left_join(sra_obj$raw, by = c("geo_accession" = "library_name_sra")) %>%
    dplyr::left_join(adt_sra, by = c("paired_adt_geo_accession" = "adt_library_name_sra")) %>%
    dplyr::mutate(
      database_accession = project_id,
      sra_study = sra_study_expected,
      bioproject = bioproject,
      organism = organism,
      sample_id = geo_accession,
      paper_sample = geo_accession,
      cellranger_dir_id = run_accession,
      cellranger_input_root = fixed_cellranger_root,
      cellranger_reprocessing_id = fixed_cellranger_reprocessing_id,
      cellranger_chemistry_parameter = fixed_cellranger_chemistry_parameter,
      sra_read1_original_role = "91bp_cDNA_read_used_as_CellRanger_R2",
      sra_read2_original_role = "26bp_cell_barcode_plus_UMI_used_as_CellRanger_R1",
      sra_r1_r2_swap_applied = "yes",
      cellranger_reprocessing_validation = "validated_all_17_GEX_runs_raw_and_filtered_matrices_issue_ok",
      cellranger_output_access = "read_only_external_reprocessing",
      parent_series = project_id,

      internal_sample_id = raw_internal_sample_id,
      patient_id = raw_internal_sample_id,
      patient_id_source = "GEO_RAW_filename_internal_sample_id_exact",
      participant_id = participant_id_base,
      participant_id_source = "GEO_RAW_filename_embedded_mapping",
      participant_id_base = participant_id_base,
      participant_id_base_inferred_from_filename = ifelse(raw_internal_sample_id != participant_id_base, "yes", "no"),
      longitudinal_pairing_status = ifelse(raw_internal_sample_id != participant_id_base, "possible_same_participant_from_filename_suffix_not_clinically_asserted", "single_public_sample_id"),
      subject_label_source = "GEO_title",
      paired_adt_internal_sample_id = adt_internal_sample_id,
      paired_adt_participant_id_base = adt_participant_id_base,
      adt_pairing_confidence = "filename_internal_sample_id_and_GEO_order_exact",
      adt_pairing_source = "GEO_RAW_filename_list_user_provided_and_embedded",
      adt_available = "yes",
      adt_loaded = "no",
      adt_excluded_reason = "RNA atlas uses GEX only; ADT retained as provenance/pairing manifest only",

      sex = NA_character_,
      sex_source = "unavailable_public_metadata",
      age = NA_real_,
      age_source = "unavailable_public_metadata",
      age_group = NA_character_,
      bmi = NA_real_,
      bmi_source = "unavailable_public_metadata",
      individual_clinical_covariates_available = "no_public_individual_demographic_table_found",
      demographic_note = "Age/sex are available only as cohort-level summaries in publications, not as sample-level public metadata; not assigned to individual samples.",
      sex_inferred_from_expression = NA_character_,
      sex_inference_score = NA_real_,
      sex_inference_method = NA_character_,

      source_name = dplyr::coalesce(source_name_sra, "whole blood"),
      tissue_label = dplyr::coalesce(tissue_sra, "whole blood"),
      source_blood_fraction = "Ficoll PBMC plus RosetteSep HLA myeloid-enriched whole-blood fraction mixture",
      sample_preparation = "Ficoll PBMC and RosetteSep HLA Myeloid Cell Enrichment Kit fractions mixed for CITE-seq",
      pbmc_isolation_method = "Ficoll-Paque PLUS density gradient centrifugation",
      myeloid_enrichment_method = "RosetteSep HLA Myeloid Cell Enrichment Kit",
      pbmc_myeloid_mixing_ratio = "enriched_PBMC:myeloid_cells = 1:3",
      fresh_or_frozen = "fresh",
      rbc_lysis = NA_character_,
      ficoll = "yes",
      frozen_or_fresh = "fresh",
      cell_type = dplyr::coalesce(cell_type_sra, "circulating immune cells"),
      extracted_molecule = "polyA RNA / total RNA listed in GEO sample metadata",

      platform = "10x Genomics",
      sequencer_model = dplyr::coalesce(instrument_sra, "Illumina NovaSeq 6000"),
      library_type = "10x Genomics 5prime single-cell gene expression; CITE-seq GEX component",
      chemistry_version = "10x Genomics v1.1 5prime chemistry",
      cite_seq = "yes",
      cite_seq_panel_n = 11L,
      cite_seq_panel_markers = "HLADR;CD3;CD15;CD127;CD183;CD66B;CD4;CD14;CD196;CD25;CD11B",
      feature_barcoding = "yes; ADT feature barcode library available as paired GEO/SRA samples",
      gex_component_loaded = "yes",
      adt_component_available = "yes",
      adt_component_loaded = "no",
      loaded_assay_component = "GEX_only",
      adt_included_in_atlas_object = "no",
      assay_type = dplyr::coalesce(assay_type_sra, "RNA-Seq"),
      library_layout = dplyr::coalesce(library_layout_sra, "PAIRED"),
      library_selection = dplyr::coalesce(library_selection_sra, "cDNA"),
      library_source = dplyr::coalesce(library_source_sra, "TRANSCRIPTOMIC SINGLE CELL"),

      aligner = "Cell Ranger common atlas reprocessing",
      cellranger_version = "atlas common Cell Ranger reprocessing from SRA FASTQ after validated R1/R2 role correction; chemistry=SC5P-R2",
      cellranger_or_equivalent_version = "atlas common Cell Ranger reprocessing from SRA FASTQ after validated R1/R2 role correction; chemistry=SC5P-R2",
      paper_original_mkfastq = "Cell Ranger v7.0.1 mkfastq according to GEO sample processing",
      paper_original_aligner = "alevin-fry v0.7.0 according to GEO sample processing; paper text may report alevin-fry v0.8.1 depending on analysis branch",
      paper_original_cellranger_version = "Cell Ranger v7.0.1 mkfastq according to GEO sample processing",
      paper_original_reference_genome = "GRCh38/hg38 spliced+intronic/splici reference in paper analysis",
      reference_genome = "GRCh38/hg38",
      gene_annotation = "atlas common Cell Ranger reference; paper used Ensembl transcript IDs aggregated to gene IDs",
      paper_qc_empty_droplet_method = "DropletUtils FDR < 0.01",
      paper_qc_max_percent_mt = 5,
      paper_qc_cells_reported_myeloid_paper = 119062L,
      paper_qc_genes_reported_myeloid_paper = 28952L,
      paper_qc_cells_reported_lymphocyte_paper = 66225L,
      paper_qc_genes_reported_lymphocyte_paper = 36601L,
      paper_cell_subset_note = "GSE252331 GEO contains 17 GEX and 17 ADT records for the GEX subset used here; ADT is excluded from RNA atlas object.",
      paper_study_design = "whole blood from healthy controls and sepsis patients at day 4 or post-sepsis day 14-21; PBMC and myeloid-enriched fractions mixed for CITE-seq",
      sepsis_definition = "Sepsis-3 sepsis or septic shock; electronic MEWS-SRS recognition and clinical protocol",
      cci_definition = "ICU length of stay >=14 days with persistent organ dysfunction or transfer/discharge with persistent organ dysfunction",
      cov19 = "no",
      is_sepsis = ifelse(factor == "sepsis", "yes", "no"),
      ards = "not_reported",
      infection_site = NA_character_,
      pathogen_etiology = NA_character_,
      septic_shock = ifelse(timepoint_label == "D4_acute_sepsis", "yes_group_level_all_day4_in_paper", NA_character_),
      mortality = NA_character_,
      disease_phase = dplyr::case_when(
        factor == "healthy_control" ~ "healthy_control",
        timepoint_label == "D4_acute_sepsis" ~ "acute_sepsis",
        TRUE ~ "post_sepsis_day14"
      ),
      clinical_mapping_status = "GEO_title_SRA_and_RAW_filename_internal_ID_mapping_no_individual_demographic_table",
      clinical_match_confidence = "sample-level condition/timepoint/outcome exact; individual demographics unavailable",
      outcome_group = rap_cci_classification,
      batch_unit_recommended = "sample/run",
      inclusion_in_primary_atlas = "yes_sepsis_core",
      processed_geo_matrix_available = "yes_alevin_fry_quants_mat_not_used_for_atlas_expression",
      processed_geo_matrix_policy = "not_used_for_primary_atlas_expression_common_cellranger_raw_feature_bc_matrix_required",
      read_mode = NA_character_
    )

  matrix_dirs <- vapply(seq_len(nrow(sample_meta)), function(i) find_matrix_dir(roots, sample_meta[i, , drop = FALSE]), character(1))
  sample_meta <- sample_meta %>%
    dplyr::mutate(
      matrix_dir = dplyr::na_if(matrix_dirs, ""),
      matrix_found = !is.na(matrix_dir)
    )

  if (nrow(sample_meta) != expected_n_gex) stop("Expected 17 GEX sample rows, got ", nrow(sample_meta))
  if (any(is.na(sample_meta$run_accession) | sample_meta$run_accession == "")) stop("SRA join failed for some GEX GEO accessions.")
  if (any(is.na(sample_meta$adt_run_accession) | sample_meta$adt_run_accession == "")) stop("SRA join failed for some paired ADT GEO accessions.")
  if (any(!sample_meta$matrix_found)) {
    print(sample_meta %>% dplyr::filter(!matrix_found) %>% dplyr::select(geo_accession, geo_title, run_accession, biosample_accession))
    stop("Some GEX samples do not have common Cell Ranger raw_feature_bc_matrix directories.")
  }

  list(sra = sra_obj, roots = roots, sample_meta = sample_meta)
}

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

read_gene_expression <- function(meta_row) {
  atlas_preprocess_10x_raw(
    raw_dir = meta_row$matrix_dir,
    sample_id = as.character(meta_row$cellranger_dir_id),
    qc_output_dir = file.path(qc_dir, "droplet_qc"),
    project_id = project_id
  )
}

attach_sample_metadata <- function(obj, meta_row) {
  atlas_attach_sample_metadata(obj, meta_row)
}

bundle <- build_sample_metadata(PROJECT_ROOT)
sample_meta <- bundle$sample_meta

output_root <- file.path(PROJECT_ROOT, "pre_integration", project_id)
rds_preqc_dir <- file.path(output_root, "rds_preqc_raw")
rds_qc_dir <- file.path(output_root, "rds_qcfiltered_raw")
rds_work_dir <- file.path(output_root, "rds_working")
meta_dir <- file.path(output_root, "metadata")
qc_dir <- file.path(output_root, "qc")
manifest_dir <- file.path(output_root, "manifest")
for (d in c(output_root, rds_preqc_dir, rds_qc_dir, rds_work_dir, meta_dir, qc_dir, manifest_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

message("PROJECT_ROOT: ", PROJECT_ROOT)
message("Cell Ranger roots: ", paste(bundle$roots, collapse = " | "))
message("SRA table   : ", bundle$sra$file)
message("output_root : ", output_root)

adt_pairing_manifest <- sample_meta %>%
  transmute(
    internal_sample_id,
    participant_id,
    participant_id_base,
    participant_id_base_inferred_from_filename,
    longitudinal_pairing_status,
    gex_geo_accession = geo_accession,
    gex_geo_title = geo_title,
    gex_run_accession = run_accession,
    gex_experiment_accession = experiment_accession,
    gex_biosample_accession = biosample_accession,
    adt_geo_accession = paired_adt_geo_accession,
    adt_geo_title,
    adt_internal_sample_id,
    adt_run_accession,
    adt_experiment_accession,
    adt_biosample_accession,
    adt_pairing_confidence,
    adt_loaded,
    adt_excluded_reason
  )

paper_methods_tbl <- tibble::tribble(
  ~field, ~value, ~note,
  "series_accession", "GSE252331", "post-septic myeloid and linked lymphocyte CITE-seq GEX/ADT series",
  "bioproject", "PRJNA1060229", "GEO relation",
  "sra_study", "SRP481033", "SRA RunTable",
  "gex_samples_loaded", "17", "GSM7999550-GSM7999566; ADT samples excluded from RNA object",
  "adt_samples_excluded", "17", "GSM7999567-GSM7999583; pairing manifest saved",
  "study_groups_in_loaded_gex", "7 healthy control, 4 acute sepsis day 4, 4 RAP day 14, 2 CCI day 14", "GEO GEX subset",
  "overall_design", "whole blood myeloid-enriched and Ficoll-enriched PBMC mixture from sepsis and healthy subjects underwent CITE-seq", "GEO/paper",
  "scRNA_protocol", "10x Genomics v1.1 5prime chemistry with feature barcoding", "paper methods",
  "sample_preparation", "Ficoll PBMC plus RosetteSep HLA myeloid-enriched cells mixed at enriched_PBMC:myeloid_cells = 1:3", "paper methods",
  "fresh_or_frozen", "fresh", "paper methods",
  "adt_panel", "HLADR;CD3;CD15;CD127;CD183;CD66B;CD4;CD14;CD196;CD25;CD11B", "GEO/paper ADT oligo panel",
  "common_atlas_matrix_source", "SRA FASTQ reprocessed after swapping SRA read_2 to Cell Ranger R1 and SRA read_1 to Cell Ranger R2; SC5P-R2; raw_feature_bc_matrix only", "atlas policy",
  "processed_geo_matrix_policy", "GEO alevin-fry quants_mat processed matrices not used for primary atlas expression", "atlas policy",
  "paper_original_quantification", "alevin-fry on GRCh38/hg38 splici reference after Cell Ranger mkfastq", "paper/GEO methods",
  "paper_qc", "DropletUtils FDR <0.01 and mitochondrial percent <5%", "paper methods",
  "individual_age_sex", "unavailable_public_metadata", "do not assign cohort-level summaries to individual samples"
)

write_csv(geo_gex_tbl, file.path(manifest_dir, "GSE252331_geo_gex_samples_embedded.v1_1.csv"))
write_csv(geo_adt_tbl, file.path(manifest_dir, "GSE252331_geo_adt_samples_excluded_embedded.v1_1.csv"))
write_csv(adt_panel_tbl, file.path(manifest_dir, "GSE252331_citeseq_adt_panel_embedded.v1_1.csv"))
write_csv(adt_pairing_manifest, file.path(manifest_dir, "GSE252331_gex_adt_pairing_manifest.v1_1.csv"))
write_csv(bundle$sra$raw, file.path(manifest_dir, "GSE252331_sra_table_normalized_all_34.v1_1.csv"))
sample_meta <- atlas_standardize_metadata(sample_meta, project_id = project_id)
write_csv(sample_meta, file.path(manifest_dir, "GSE252331_pool_sample_metadata_preloop.v1_1.csv"))
write_csv(paper_methods_tbl, file.path(manifest_dir, "GSE252331_paper_methods_and_design_embedded.v1_1.csv"))

common_meta_columns <- c(
  "database_accession", "sra_study", "bioproject", "organism",
  "run_accession", "experiment_accession", "biosample_accession", "geo_accession",
  "cellranger_dir_id", "cellranger_input_root", "cellranger_reprocessing_id",
  "cellranger_chemistry_parameter", "sra_read1_original_role", "sra_read2_original_role",
  "sra_r1_r2_swap_applied", "cellranger_reprocessing_validation", "cellranger_output_access",
  "sample_id", "paper_sample", "geo_title",
  "internal_sample_id", "raw_internal_sample_id", "patient_id", "patient_id_source",
  "participant_id", "participant_id_base",
  "participant_id_source", "participant_id_base_inferred_from_filename", "longitudinal_pairing_status",
  "subject_label", "subject_label_source", "factor", "study_group",
  "timepoint_label", "days", "rap_cci_classification", "outcome_group", "severity",
  "cov19", "is_sepsis", "ards", "disease_phase",
  "clinical_mapping_status", "clinical_match_confidence",
  "sex", "sex_source", "age", "age_source", "age_group", "bmi", "bmi_source",
  "individual_clinical_covariates_available", "demographic_note",
  "sex_inferred_from_expression", "sex_inference_score", "sex_inference_method",
  "infection_site", "pathogen_etiology", "septic_shock", "mortality",
  "source_name", "tissue_label", "source_blood_fraction", "cell_type", "extracted_molecule",
  "sample_preparation", "pbmc_isolation_method", "myeloid_enrichment_method", "pbmc_myeloid_mixing_ratio",
  "fresh_or_frozen", "frozen_or_fresh", "ficoll", "rbc_lysis",
  "platform", "sequencer_model", "library_type", "chemistry_version",
  "cite_seq", "cite_seq_panel_n", "cite_seq_panel_markers", "feature_barcoding",
  "gex_component_loaded", "adt_component_available", "adt_component_loaded",
  "loaded_assay_component", "adt_included_in_atlas_object",
  "paired_adt_geo_accession", "paired_adt_internal_sample_id", "paired_adt_participant_id_base",
  "adt_run_accession", "adt_experiment_accession", "adt_biosample_accession",
  "adt_pairing_confidence", "adt_pairing_source", "adt_available", "adt_loaded", "adt_excluded_reason",
  "assay_type", "library_layout", "library_selection", "library_source",
  "aligner", "cellranger_version", "cellranger_or_equivalent_version",
  "paper_original_mkfastq", "paper_original_aligner", "paper_original_cellranger_version",
  "paper_original_reference_genome", "reference_genome", "gene_annotation",
  "paper_qc_empty_droplet_method", "paper_qc_max_percent_mt",
  "paper_cell_subset_note", "paper_study_design", "sepsis_definition", "cci_definition",
  "batch_unit_recommended", "inclusion_in_primary_atlas",
  "processed_geo_matrix_available", "processed_geo_matrix_policy",
  "matrix_dir", "matrix_found", "read_mode"
)
sample_meta <- ensure_columns(sample_meta, common_meta_columns)

sample_meta <- atlas_standardize_metadata(sample_meta, project_id = project_id)

obj_list <- list()
qc_summary <- tibble(
  cellranger_dir_id = character(),
  run_accession = character(),
  biosample_accession = character(),
  geo_accession = character(),
  internal_sample_id = character(),
  patient_id = character(),
  participant_id = character(),
  factor = character(),
  study_group = character(),
  timepoint_label = character(),
  days = integer(),
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

for (i in seq_len(nrow(sample_meta))) {
  meta_row <- sample_meta[i, , drop = FALSE]
  dir_id <- meta_row$cellranger_dir_id

  message("Reading ", dir_id, " -> ", meta_row$geo_accession, " / ", meta_row$internal_sample_id, " (", meta_row$geo_title, ")")
  # ATLAS_RESOURCE_SAFE_EXCLUSION_GSE252331
  res <- tryCatch(
    read_gene_expression(meta_row),
    atlas_sample_exclusion = function(e) e
  )
  if (inherits(res, "atlas_sample_exclusion")) {
    message(
      "[EXCLUDED] ", dir_id, " -> ", meta_row$geo_accession,
      ": ", res$reason, " (called=", res$n_called_after_droplet_qc,
      ", minimum=", res$minimum_called_cells, ")"
    )
    qc_summary <- bind_rows(qc_summary, tibble(
        cellranger_dir_id = dir_id,
        run_accession = meta_row$run_accession,
        biosample_accession = meta_row$biosample_accession,
        geo_accession = meta_row$geo_accession,
        internal_sample_id = meta_row$internal_sample_id,
        patient_id = meta_row$patient_id,
        participant_id = meta_row$participant_id,
        factor = meta_row$factor,
        study_group = meta_row$study_group,
        timepoint_label = meta_row$timepoint_label,
        days = meta_row$days,
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
      ))
    invisible(gc())
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

  n_raw <- ncol(obj)
  sample_meta$read_mode[i] <- res$read_mode
  meta_row$read_mode <- res$read_mode

  obj <- attach_sample_metadata(obj, meta_row)

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- compute_percent_ribo(obj)
  obj[["percent.hb"]] <- compute_percent_hb(obj)
  obj[["paper_qc_pass"]] <- with(obj@meta.data, percent.mt <= paper_qc_max_percent_mt)
  obj[["atlas_basic_qc_pass"]] <- with(obj@meta.data, nFeature_RNA >= 200)
  obj[["atlas_qc_pass"]] <- with(obj@meta.data, paper_qc_pass & atlas_basic_qc_pass)
  obj <- atlas_standardize_seurat_metadata(obj, project_id = project_id)

  preqc_rds <- file.path(rds_preqc_dir, paste0(meta_row$geo_accession, "__", meta_row$internal_sample_id, "__", meta_row$cellranger_dir_id, "__preqc_raw.rds"))
  saveRDS(obj, preqc_rds, compress = TRUE)

  cells_keep <- colnames(obj)[which(obj$atlas_qc_pass)]
  n_after_qc <- length(cells_keep)

  if (n_after_qc < 200L) {
    warning(meta_row$geo_accession, " dropped because only ", n_after_qc, " cells remained after QC.")
    qc_summary <- bind_rows(qc_summary, tibble(
      cellranger_dir_id = dir_id, run_accession = meta_row$run_accession,
      biosample_accession = meta_row$biosample_accession, geo_accession = meta_row$geo_accession,
      internal_sample_id = meta_row$internal_sample_id, patient_id = meta_row$patient_id,
      participant_id = meta_row$participant_id,
      factor = meta_row$factor, study_group = meta_row$study_group,
      timepoint_label = meta_row$timepoint_label, days = meta_row$days,
      n_raw = n_raw, n_after_qc = n_after_qc, saved = FALSE,
      preqc_rds_path = preqc_rds, qc_rds_path = NA_character_, working_rds_path = NA_character_
    ))
    next
  }

  obj_qc <- subset(obj, cells = cells_keep)
  qc_rds <- file.path(rds_qc_dir, paste0(meta_row$geo_accession, "__", meta_row$internal_sample_id, "__", meta_row$cellranger_dir_id, "__qcfiltered_raw.rds"))
  saveRDS(obj_qc, qc_rds, compress = TRUE)

  obj_work <- NormalizeData(obj_qc, verbose = FALSE)
  obj_work <- FindVariableFeatures(obj_work, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  work_rds <- file.path(rds_work_dir, paste0(meta_row$geo_accession, "__", meta_row$internal_sample_id, "__", meta_row$cellranger_dir_id, "__working.rds"))
  saveRDS(obj_work, work_rds, compress = TRUE)

  obj_list[[as.character(meta_row$geo_accession)]] <- work_rds

  qc_summary <- bind_rows(qc_summary, tibble(
    cellranger_dir_id = dir_id, run_accession = meta_row$run_accession,
    biosample_accession = meta_row$biosample_accession, geo_accession = meta_row$geo_accession,
    internal_sample_id = meta_row$internal_sample_id, patient_id = meta_row$patient_id,
    participant_id = meta_row$participant_id,
    factor = meta_row$factor, study_group = meta_row$study_group,
    timepoint_label = meta_row$timepoint_label, days = meta_row$days,
    n_raw = n_raw, n_after_qc = n_after_qc, saved = TRUE,
    preqc_rds_path = preqc_rds, qc_rds_path = qc_rds, working_rds_path = work_rds
  ))
  # ATLAS_RESOURCE_SAFE_GC_GSE252331
  rm(res, obj, obj_qc, obj_work)
  invisible(gc())
}

write_csv(sample_meta, file.path(meta_dir, "GSE252331_sample_metadata.csv"))
write_csv(sample_meta, file.path(meta_dir, "GSE252331_run_metadata.csv"))
write_csv(qc_summary, file.path(qc_dir, "GSE252331_qc_summary.csv"))
write_csv(
  qc_summary %>% dplyr::filter(preprocess_status %in% "excluded_before_rds"),
  file.path(qc_dir, "GSE252331_excluded_before_rds.csv")
)

registry <- qc_summary %>%
  filter(saved) %>%
  left_join(sample_meta %>%
              select(cellranger_dir_id, run_accession, biosample_accession, geo_accession,
                     internal_sample_id, patient_id, patient_id_source,
                     participant_id, participant_id_base,
                     factor, study_group, timepoint_label, days,
                     rap_cci_classification, severity, cov19, is_sepsis,
                     sex, sex_source, age, age_source,
                     cite_seq, loaded_assay_component, adt_component_available,
                     paired_adt_geo_accession, adt_run_accession,
                     clinical_mapping_status, batch_unit_recommended, inclusion_in_primary_atlas),
            by = c("cellranger_dir_id", "run_accession", "biosample_accession", "geo_accession",
                   "internal_sample_id", "patient_id", "participant_id",
                   "factor", "study_group", "timepoint_label", "days")) %>%
  mutate(project_id = project_id, preintegration_tier = "working")

write_csv(registry, file.path(output_root, "GSE252331_preintegration_registry.csv"))
saveRDS(obj_list, file.path(output_root, "GSE252331_preintegration_objlist.rds"))

message("Done. Saved pre-QC raw objects to    : ", rds_preqc_dir)
message("Done. Saved QC-filtered raw objects to: ", rds_qc_dir)
message("Done. Saved working objects to       : ", rds_work_dir)
message("Registry: ", file.path(output_root, "GSE252331_preintegration_registry.csv"))
