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

# ==============================================================================
# GSE216020 pre-integration loader
#
# Dataset:
#   GSE216020 / SRP403283 / PRJNA891717
#   "Heterogeneity of neutrophils and inflammatory responses in patients with
#    COVID-19 and healthy controls"
#
# Atlas-side decision:
#   Public raw reads were reprocessed through the user's common atlas pipeline.
#   Paper-side Cell Ranger / Seurat / DoubletFinder / QC details are stored as
#   provenance metadata only; the raw matrices used here are the user's common
#   Cell Ranger outputs.
#
# Expected external input:
#   configs/GSE216020/GSE216020_SraRunTable.csv
#   or configs/GSE216020/SraRunTable.csv
#
# Expected Cell Ranger output root:
#   cellranger_count/output/SRP403283/<BioSample>/outs/raw_feature_bc_matrix
#
# Output follows the same pre-integration contract as prior atlas loaders:
#   rds_preqc_raw/
#   rds_qcfiltered_raw/
#   rds_working/
#   metadata/
#   manifest/
#   qc/
# ==============================================================================

project_id <- "GSE216020"
sra_study  <- "SRP403283"
bioproject <- "PRJNA891717"
organism   <- "Homo sapiens"

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

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
  if (required) {
    stop(label, " is missing a required column. Tried exact = ",
         paste(exact, collapse = ", "), "; regex = ", regex)
  }
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

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

parse_day <- function(x) {
  x <- as.character(x)
  out <- suppressWarnings(as.numeric(stringr::str_extract(x, "[0-9]+")))
  out
}

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
      paste0(" - ", candidates, collapse = "\n")
    )
  }
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

resolve_sra_table <- function(project_root, project_id) {
  candidates <- c(
    file.path(project_root, "configs", project_id, paste0(project_id, "_SraRunTable.csv")),
    file.path(project_root, "configs", project_id, "SraRunTable.csv")
  )
  hits <- candidates[file.exists(candidates)]
  if (length(hits) == 0) {
    stop(
      "Missing SRA run table. Tried:\n",
      paste0(" - ", candidates, collapse = "\n")
    )
  }
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

ensure_columns <- function(df, wanted) {
  miss <- setdiff(wanted, names(df))
  if (length(miss) > 0) {
    for (m in miss) df[[m]] <- NA
  }
  df[, unique(c(wanted, names(df))), drop = FALSE]
}

safe_save_rds <- function(obj, path) {
  saveRDS(obj, file = path, compress = TRUE)
  invisible(path)
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

run_optional_helpers <- function(obj) {
  # Reserved hook for future common helper insertion.
  obj
}

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
input_root    <- resolve_input_root(PROJECT_ROOT, sra_study)
sra_run_table <- resolve_sra_table(PROJECT_ROOT, project_id)

output_root   <- file.path(PROJECT_ROOT, "pre_integration", project_id)
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

message("PROJECT_ROOT   : ", PROJECT_ROOT)
message("input_root     : ", input_root)
message("SRA run table  : ", sra_run_table)
message("output_root    : ", output_root)

# ------------------------------------------------------------------------------
# Atlas-side QC/output policy
# ------------------------------------------------------------------------------
apply_paper_qc_filter <- TRUE
min_cells_after_qc    <- 200L

# ------------------------------------------------------------------------------
# Embedded GEO sample mapping
# ------------------------------------------------------------------------------
geo_sample_tbl <- tibble::tribble(
  ~geo_accession, ~paper_sample, ~participant_id, ~library_batch_id, ~library_unit_id,
  "GSM6656081", "VCM01 1577-JX-1",  "VCM01", "1577", "1577-JX-1",
  "GSM6656082", "UCM01 1577-JX-2",  "UCM01", "1577", "1577-JX-2",
  "GSM6656083", "UCM01 1577-JX-3",  "UCM01", "1577", "1577-JX-3",
  "GSM6656084", "VCM02 1577-JX-4",  "VCM02", "1577", "1577-JX-4",
  "GSM6656085", "VCM03 1577-JX-5",  "VCM03", "1577", "1577-JX-5",
  "GSM6656086", "VCM02 1577-JX-6",  "VCM02", "1577", "1577-JX-6",
  "GSM6656087", "VCM03 1577-JX-7",  "VCM03", "1577", "1577-JX-7",
  "GSM6656088", "VCM04 2654-JX-1",  "VCM04", "2654", "2654-JX-1",
  "GSM6656089", "VCM04 2654-JX-2",  "VCM04", "2654", "2654-JX-2",
  "GSM6656090", "VCM05 2654-JX-3",  "VCM05", "2654", "2654-JX-3",
  "GSM6656091", "VCM06 2654-JX-4",  "VCM06", "2654", "2654-JX-4",
  "GSM6656092", "VCM05 2654-JX-5",  "VCM05", "2654", "2654-JX-5",
  "GSM6656093", "VCM06 2654-JX-6",  "VCM06", "2654", "2654-JX-6",
  "GSM6656094", "VSC01 2654-JX-7",  "VSC01", "2654", "2654-JX-7",
  "GSM6656095", "VCH01 2654-JX-8",  "VCH01", "2654", "2654-JX-8",
  "GSM6656096", "VCH02 2654-JX-9",  "VCH02", "2654", "2654-JX-9",
  "GSM6656097", "VCH03 2654-JX-10", "VCH03", "2654", "2654-JX-10",
  "GSM6656098", "VCH04 2654-JX-11", "VCH04", "2654", "2654-JX-11",
  "GSM6656099", "VCH05 2654-JX-12", "VCH05", "2654", "2654-JX-12",
  "GSM6656100", "VCM08 2972-JX-1",  "VCM08", "2972", "2972-JX-1",
  "GSM6656101", "VCM09 2972-JX-2",  "VCM09", "2972", "2972-JX-2",
  "GSM6656102", "VCM08 2972-JX-3",  "VCM08", "2972", "2972-JX-3",
  "GSM6656103", "VCM09 2972-JX-4",  "VCM09", "2972", "2972-JX-4",
  "GSM6656104", "VCM10 2972-JX-5",  "VCM10", "2972", "2972-JX-5"
)

# Subject metadata extracted from Supplementary Fig. S1B.
# outcome_status is exact category from the figure (discharge/death/healthy);
# outcome days are not encoded as exact numeric values because the figure is a timeline.
subject_tbl <- tibble::tribble(
  ~participant_id, ~study_group_from_supp_fig1b, ~sex,   ~age, ~race,  ~final_outcome, ~mortality, ~covid19_severity_numeric_paper,
  "VCH01", "healthy_control",      "male", 37, "White", "healthy",   "no", 1,
  "VCH02", "healthy_control",      "male", 74, "White", "healthy",   "no", 1,
  "VCH03", "healthy_control",      "male", 62, "White", "healthy",   "no", 1,
  "VCH04", "healthy_control",      "male", 65, "White", "healthy",   "no", 1,
  "VCH05", "healthy_control",      "male", 43, "Black", "healthy",   "no", 1,
  "VCM01", "mild_covid19",         "male", 70, "White", "recovered", "no", 2,
  "VCM03", "mild_covid19",         "male", 73, "White", "recovered", "no", 2,
  "VCM09", "mild_covid19",         "male", 74, "White", "recovered", "no", 2,
  "VCM10", "mild_covid19",         "male", 49, "White", "recovered", "no", 2,
  "UCM01", "severe_covid19",       "male", 42, "Black", "recovered", "no", 3,
  "VCM02", "severe_covid19",       "male", 78, "White", "death",     "yes", 4,
  "VCM04", "severe_covid19",       "male", 83, "White", "death",     "yes", 4,
  "VCM05", "severe_covid19",       "male", 72, "White", "death",     "yes", 4,
  "VCM06", "severe_covid19",       "male", 72, "White", "death",     "yes", 4,
  "VCM08", "severe_covid19",       "male", 64, "White", "recovered", "no", 3,
  "VSC01", "severe_covid19",       "male", 65, "White", "recovered", "no", 3
)

paper_methods_tbl <- tibble::tribble(
  ~field, ~value, ~note,
  "paper_title", "Heterogeneity of neutrophils and inflammatory responses in patients with COVID-19 and healthy controls", "Frontiers in Immunology 2022",
  "cohort_total_covid19_patients", "11", "adult hospitalized SARS-CoV-2 infection patients",
  "cohort_healthy_controls", "5", "healthy controls, all outpatients",
  "mild_definition", "hospitalized but needing <=50% oxygen", "paper Methods",
  "severe_definition", "needing >50% oxygen or ICU", "paper Methods",
  "severe_deaths", "4 of 7 severe patients", "paper Methods",
  "sampling_t1", "within 72 hours of hospital admission", "paper Methods",
  "sampling_t2", "5-7 days later", "paper Methods",
  "sampling_exception", "one patient sampled 15 days later", "paper Methods; in SRA table this appears as day 20",
  "blood_collection_tube", "lavender top EDTA tube", "paper Methods",
  "isolation_medium", "Lymphocyte-poly isolation media", "paper Methods",
  "centrifugation_1", "500g, 35 minutes, room temperature", "paper Methods",
  "centrifugation_2", "350g, 10 minutes, room temperature", "paper Methods",
  "rbc_lysis", "1x ACK lysis buffer for 2 minutes", "paper Methods",
  "cell_suspension_buffer", "PBS containing 0.04% BSA", "paper Methods",
  "cell_mixture", "mononuclear and granulocyte populations combined at 1:1 ratio", "paper Methods",
  "fresh_processing", "immediate processing for single-cell RNA library", "paper Methods",
  "minimum_viability", ">90%", "paper Methods",
  "chromium_kit", "Chromium Next GEM Single Cell 3prime Reagent Kit v3.1", "paper Methods",
  "cell_suspension_concentration", "700-1,200 cells per ul", "paper Methods",
  "sequencer", "Illumina NovaSeq", "paper Methods and GEO GPL24676",
  "paper_original_cellranger", "Cell Ranger 6.0.0", "paper provenance only",
  "paper_original_reference", "GRCh38", "paper provenance",
  "paper_qc_min_umi", "500", "paper provenance",
  "paper_qc_min_genes", "200", "paper provenance",
  "paper_qc_max_percent_mito", "20", "paper provenance",
  "paper_doubletfinder_expected_rate", "0.075", "paper provenance; not applied here",
  "paper_cells_after_qc", "108597", "paper Results",
  "paper_hd_cells_after_qc", "30429", "paper Results",
  "paper_mild_cells_after_qc", "22188", "paper Results",
  "paper_severe_cells_after_qc", "55980", "paper Results",
  "paper_neutrophil_cells", "45463", "paper Results",
  "paper_major_cell_types", "15", "paper Results",
  "paper_neutrophil_clusters", "9", "7 mature and 2 immature neutrophil clusters"
)

write_csv(geo_sample_tbl, file.path(manifest_dir, "GSE216020_geo_sample_embedded.csv"))
write_csv(subject_tbl, file.path(manifest_dir, "GSE216020_subject_metadata_from_supp_fig1b_embedded.csv"))
write_csv(paper_methods_tbl, file.path(manifest_dir, "GSE216020_paper_methods_and_results_embedded.csv"))

# ------------------------------------------------------------------------------
# Read and normalize SRA run table
# ------------------------------------------------------------------------------
sra_raw <- readr::read_csv(sra_run_table, show_col_types = FALSE)
names(sra_raw) <- normalize_names(names(sra_raw))

run_col       <- pick_column(sra_raw, exact = c("run", "run_accession"), regex = "(^|_)run(_|$)", required = TRUE, label = "SRA run table")
exp_col       <- pick_column(sra_raw, exact = c("experiment", "experiment_accession"), regex = "(^|_)experiment(_|$)", required = TRUE, label = "SRA run table")
biosample_col <- pick_column(sra_raw, exact = c("biosample", "biosample_accession", "bio_sample"), regex = "(^|_)biosample(_|$)|(^|_)bio_sample(_|$)", required = TRUE, label = "SRA run table")

sample_name_col <- pick_column(sra_raw, exact = c("sample_name"), regex = "^sample_name$", required = FALSE)
library_name_col <- pick_column(sra_raw, exact = c("library_name"), regex = "^library_name$", required = FALSE)
geo_col <- pick_column(sra_raw, exact = c("geo_accession_exp", "geo_accession"), regex = "geo.*accession", required = FALSE)
bioproject_col <- pick_column(sra_raw, exact = c("bioproject"), regex = "^bioproject$", required = FALSE)
sra_study_col <- pick_column(sra_raw, exact = c("sra_study"), regex = "^sra_study$", required = FALSE)
assay_col <- pick_column(sra_raw, exact = c("assay_type"), regex = "assay.*type", required = FALSE)
avgspot_col <- pick_column(sra_raw, exact = c("avgspotlen"), regex = "avgspotlen", required = FALSE)
bases_col <- pick_column(sra_raw, exact = c("bases"), regex = "^bases$", required = FALSE)
bytes_col <- pick_column(sra_raw, exact = c("bytes"), regex = "^bytes$", required = FALSE)
cell_line_col <- pick_column(sra_raw, exact = c("cell_line"), regex = "^cell_line$", required = FALSE)
cell_type_col <- pick_column(sra_raw, exact = c("cell_type"), regex = "^cell_type$", required = FALSE)
center_col <- pick_column(sra_raw, exact = c("center_name"), regex = "center.*name", required = FALSE)
consent_col <- pick_column(sra_raw, exact = c("consent"), regex = "^consent$", required = FALSE)
instrument_col <- pick_column(sra_raw, exact = c("instrument"), regex = "^instrument$", required = FALSE)
library_layout_col <- pick_column(sra_raw, exact = c("librarylayout"), regex = "librarylayout", required = FALSE)
library_selection_col <- pick_column(sra_raw, exact = c("libraryselection"), regex = "libraryselection", required = FALSE)
library_source_col <- pick_column(sra_raw, exact = c("librarysource"), regex = "librarysource", required = FALSE)
organism_col <- pick_column(sra_raw, exact = c("organism"), regex = "^organism$", required = FALSE)
platform_col <- pick_column(sra_raw, exact = c("platform"), regex = "^platform$", required = FALSE)
release_col <- pick_column(sra_raw, exact = c("releasedate"), regex = "releasedate", required = FALSE)
create_col <- pick_column(sra_raw, exact = c("create_date"), regex = "create_date", required = FALSE)
version_col <- pick_column(sra_raw, exact = c("version"), regex = "^version$", required = FALSE)
source_name_col <- pick_column(sra_raw, exact = c("source_name"), regex = "^source_name$", required = FALSE)
time_col <- pick_column(sra_raw, exact = c("time"), regex = "^time$", required = FALSE)
tissue_col <- pick_column(sra_raw, exact = c("tissue"), regex = "^tissue$", required = FALSE)
genotype_col <- pick_column(sra_raw, exact = c("genotype"), regex = "^genotype$", required = FALSE)
filetype_col <- pick_column(sra_raw, exact = c("datastore_filetype"), regex = "datastore.*filetype", required = FALSE)
provider_col <- pick_column(sra_raw, exact = c("datastore_provider"), regex = "datastore.*provider", required = FALSE)
region_col <- pick_column(sra_raw, exact = c("datastore_region"), regex = "datastore.*region", required = FALSE)

geo_accession_from_sra <- dplyr::coalesce(
  col_or_na(sra_raw, geo_col),
  col_or_na(sra_raw, sample_name_col),
  col_or_na(sra_raw, library_name_col)
)

sra_norm <- tibble(
  run_accession = col_or_na(sra_raw, run_col),
  experiment_accession = col_or_na(sra_raw, exp_col),
  biosample_accession = col_or_na(sra_raw, biosample_col),
  geo_accession = geo_accession_from_sra,
  sample_name_sra = col_or_na(sra_raw, sample_name_col),
  library_name_sra = col_or_na(sra_raw, library_name_col),
  bioproject_sra = col_or_na(sra_raw, bioproject_col),
  sra_study_sra = col_or_na(sra_raw, sra_study_col),
  assay_type_sra = col_or_na(sra_raw, assay_col),
  avg_spot_len_sra = col_or_na(sra_raw, avgspot_col, "numeric"),
  bases_sra = col_or_na(sra_raw, bases_col, "numeric"),
  bytes_sra = col_or_na(sra_raw, bytes_col, "numeric"),
  cell_line_sra = col_or_na(sra_raw, cell_line_col),
  cell_type_sra = col_or_na(sra_raw, cell_type_col),
  center_name_sra = col_or_na(sra_raw, center_col),
  consent_sra = col_or_na(sra_raw, consent_col),
  instrument_sra = col_or_na(sra_raw, instrument_col),
  library_layout_sra = col_or_na(sra_raw, library_layout_col),
  library_selection_sra = col_or_na(sra_raw, library_selection_col),
  library_source_sra = col_or_na(sra_raw, library_source_col),
  organism_sra = col_or_na(sra_raw, organism_col),
  platform_sra = col_or_na(sra_raw, platform_col),
  release_date_sra = col_or_na(sra_raw, release_col),
  create_date_sra = col_or_na(sra_raw, create_col),
  version_sra = col_or_na(sra_raw, version_col),
  source_name_sra = col_or_na(sra_raw, source_name_col),
  genotype_sra = col_or_na(sra_raw, genotype_col),
  time_sra = col_or_na(sra_raw, time_col),
  tissue_sra = col_or_na(sra_raw, tissue_col),
  datastore_filetype_sra = col_or_na(sra_raw, filetype_col),
  datastore_provider_sra = col_or_na(sra_raw, provider_col),
  datastore_region_sra = col_or_na(sra_raw, region_col)
) %>%
  mutate(across(where(is.character), ~na_if(.x, "")))

write_csv(sra_norm, file.path(manifest_dir, "GSE216020_sra_table_normalized.csv"))

# Strict validation of SRA table identity.
if (!any(sra_norm$bioproject_sra == bioproject, na.rm = TRUE)) {
  stop("SRA run table does not contain expected BioProject ", bioproject)
}
if (!any(sra_norm$sra_study_sra == sra_study, na.rm = TRUE)) {
  stop("SRA run table does not contain expected SRA Study ", sra_study)
}

sra_sample_agg <- sra_norm %>%
  filter(geo_accession %in% geo_sample_tbl$geo_accession) %>%
  group_by(geo_accession, biosample_accession, experiment_accession) %>%
  summarise(
    run_accession = paste(unique(na.omit(run_accession)), collapse = ";"),
    n_sra_runs = n_distinct(run_accession),
    sample_name_sra = paste(unique(na.omit(sample_name_sra)), collapse = ";"),
    library_name_sra = paste(unique(na.omit(library_name_sra)), collapse = ";"),
    bioproject_sra = paste(unique(na.omit(bioproject_sra)), collapse = ";"),
    sra_study_sra = paste(unique(na.omit(sra_study_sra)), collapse = ";"),
    assay_type_sra = paste(unique(na.omit(assay_type_sra)), collapse = ";"),
    avg_spot_len_sra = suppressWarnings(mean(avg_spot_len_sra, na.rm = TRUE)),
    bases_sra = suppressWarnings(sum(bases_sra, na.rm = TRUE)),
    bytes_sra = suppressWarnings(sum(bytes_sra, na.rm = TRUE)),
    cell_line_sra = paste(unique(na.omit(cell_line_sra)), collapse = ";"),
    cell_type_sra = paste(unique(na.omit(cell_type_sra)), collapse = ";"),
    center_name_sra = paste(unique(na.omit(center_name_sra)), collapse = ";"),
    consent_sra = paste(unique(na.omit(consent_sra)), collapse = ";"),
    instrument_sra = paste(unique(na.omit(instrument_sra)), collapse = ";"),
    library_layout_sra = paste(unique(na.omit(library_layout_sra)), collapse = ";"),
    library_selection_sra = paste(unique(na.omit(library_selection_sra)), collapse = ";"),
    library_source_sra = paste(unique(na.omit(library_source_sra)), collapse = ";"),
    organism_sra = paste(unique(na.omit(organism_sra)), collapse = ";"),
    platform_sra = paste(unique(na.omit(platform_sra)), collapse = ";"),
    release_date_sra = paste(unique(na.omit(release_date_sra)), collapse = ";"),
    create_date_sra = paste(unique(na.omit(create_date_sra)), collapse = ";"),
    version_sra = paste(unique(na.omit(version_sra)), collapse = ";"),
    source_name_sra = paste(unique(na.omit(source_name_sra)), collapse = ";"),
    genotype_sra = paste(unique(na.omit(genotype_sra)), collapse = ";"),
    time_sra = paste(unique(na.omit(time_sra)), collapse = ";"),
    tissue_sra = paste(unique(na.omit(tissue_sra)), collapse = ";"),
    datastore_filetype_sra = paste(unique(na.omit(datastore_filetype_sra)), collapse = ";"),
    datastore_provider_sra = paste(unique(na.omit(datastore_provider_sra)), collapse = ";"),
    datastore_region_sra = paste(unique(na.omit(datastore_region_sra)), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(
    across(where(is.character), ~na_if(.x, "")),
    avg_spot_len_sra = ifelse(is.nan(avg_spot_len_sra) | is.infinite(avg_spot_len_sra), NA_real_, avg_spot_len_sra),
    bases_sra = ifelse(is.nan(bases_sra) | is.infinite(bases_sra), NA_real_, bases_sra),
    bytes_sra = ifelse(is.nan(bytes_sra) | is.infinite(bytes_sra), NA_real_, bytes_sra)
  )

missing_geo <- setdiff(geo_sample_tbl$geo_accession, sra_sample_agg$geo_accession)
if (length(missing_geo) > 0 || nrow(sra_sample_agg) != nrow(geo_sample_tbl)) {
  stop(
    "SRA-to-GEO matching is incomplete for GSE216020.\n",
    "Matched GEO: ", paste(sra_sample_agg$geo_accession, collapse = ", "), "\n",
    "Missing GEO: ", paste(missing_geo, collapse = ", "), "\n",
    "Expected exactly 24 GEO samples."
  )
}
if (any(is.na(sra_sample_agg$run_accession) | sra_sample_agg$run_accession == "")) {
  stop("One or more GSE216020 sample rows lack run_accession after SRA aggregation.")
}

write_csv(sra_sample_agg, file.path(manifest_dir, "GSE216020_sra_sample_aggregated.csv"))

# ------------------------------------------------------------------------------
# Build sample metadata
# ------------------------------------------------------------------------------
sample_meta <- geo_sample_tbl %>%
  left_join(sra_sample_agg, by = "geo_accession") %>%
  left_join(subject_tbl, by = "participant_id") %>%
  mutate(
    database_accession = project_id,
    sra_study = sra_study,
    bioproject = bioproject,
    organism = organism,

    cellranger_dir_id = biosample_accession,
    sample_id = geo_accession,
    sample_alias = paper_sample,
    patient_id = participant_id,
    subject_id = participant_id,

    days_after_hospital_admission = parse_day(time_sra),
    days = days_after_hospital_admission,
    timepoint_label = dplyr::case_when(
      days_after_hospital_admission == 0 ~ "T1_day0",
      days_after_hospital_admission == 5 ~ "T2_day5",
      days_after_hospital_admission > 5 ~ "T_late",
      TRUE ~ NA_character_
    ),
    timepoint_note = dplyr::case_when(
      days_after_hospital_admission == 0 ~ "within 72 h of hospital admission in paper design",
      days_after_hospital_admission == 5 ~ "5-7 days after first timepoint in paper design",
      days_after_hospital_admission > 5 ~ "late sample; paper notes one patient sampled substantially later",
      TRUE ~ NA_character_
    ),

    study_group = dplyr::case_when(
      genotype_sra == "Healthy" ~ "healthy_control",
      genotype_sra == "Mild" ~ "covid19_mild",
      genotype_sra == "Severe" & mortality == "yes" ~ "covid19_severe_deceased",
      genotype_sra == "Severe" & mortality == "no" ~ "covid19_severe_recovered",
      TRUE ~ NA_character_
    ),
    outcome_label = dplyr::case_when(
      final_outcome == "healthy" ~ "healthy_control",
      final_outcome == "recovered" ~ "recovered",
      final_outcome == "death" ~ "death",
      TRUE ~ NA_character_
    ),
    factor = dplyr::case_when(
      genotype_sra == "Healthy" ~ "healthy_control",
      genotype_sra == "Mild" ~ "covid19_mild",
      genotype_sra == "Severe" ~ "covid19_severe",
      TRUE ~ NA_character_
    ),
    cov19 = dplyr::case_when(
      genotype_sra == "Healthy" ~ "NEG",
      genotype_sra %in% c("Mild", "Severe") ~ "POS",
      TRUE ~ NA_character_
    ),
    is_sepsis = "no",
    ards = dplyr::case_when(
      genotype_sra == "Severe" ~ "possible_or_not_reported",
      TRUE ~ NA_character_
    ),
    pneumonia = dplyr::case_when(
      genotype_sra %in% c("Mild", "Severe") ~ "covid19_respiratory_infection",
      TRUE ~ NA_character_
    ),
    infection_site = dplyr::case_when(
      genotype_sra %in% c("Mild", "Severe") ~ "respiratory/systemic blood immune response",
      TRUE ~ NA_character_
    ),
    pathogen_etiology = dplyr::case_when(
      genotype_sra %in% c("Mild", "Severe") ~ "SARS-CoV-2",
      TRUE ~ NA_character_
    ),
    icu_admission = dplyr::case_when(
      genotype_sra == "Severe" ~ "yes_or_high_oxygen_requirement",
      genotype_sra == "Mild" ~ "no_or_not_reported",
      TRUE ~ NA_character_
    ),
    severity = dplyr::case_when(
      genotype_sra == "Healthy" ~ "healthy_control",
      genotype_sra == "Mild" ~ "mild_hospitalized_covid19",
      genotype_sra == "Severe" & mortality == "yes" ~ "severe_covid19_deceased",
      genotype_sra == "Severe" & mortality == "no" ~ "severe_covid19_recovered",
      TRUE ~ NA_character_
    ),
    severityatday = paste0(severity, "_", timepoint_label),
    ventilation = NA_character_,
    qsofa = NA_real_,
    apache_ii = NA_real_,
    sofa = NA_real_,

    # GEO/SRA and sample-processing metadata.
    source_name = dplyr::coalesce(source_name_sra, "blood"),
    tissue_label = dplyr::coalesce(tissue_sra, "blood"),
    source_blood_fraction = "fresh blood leukocyte PMN/non-PMN mixture",
    cell_type = dplyr::coalesce(cell_type_sra, "immune cells/neutrophils from fresh blood"),
    cell_line = dplyr::coalesce(cell_line_sra, "Not a cell line"),
    extracted_molecule = "total RNA",
    fresh_or_frozen = "fresh, immediately processed",
    blood_collection_tube = "lavender top EDTA tube",
    pbmc_isolation_method = "Lymphocyte-poly isolation medium; mononuclear and granulocyte bands collected",
    granulocyte_isolation_method = "Lymphocyte-poly isolation medium; granulocyte-containing band collected",
    ficoll = "Lymphocyte-poly isolation media",
    rbc_lysis = "1x ACK lysis buffer for 2 minutes",
    sample_storage = "fresh, not cryopreserved",
    cell_suspension_buffer = "PBS with 0.04% BSA",
    cell_mixture = "mononuclear and granulocyte populations combined at 1:1 ratio",
    viability_reported = ">90%",
    sorted_population = "fresh peripheral blood immune cells enriched for neutrophil analysis; PMN and non-PMN mixed 1:1",

    # Technology/provenance.
    platform = "10x Genomics",
    sequencer_model = dplyr::coalesce(instrument_sra, "Illumina NovaSeq 6000"),
    library_type = "10x Genomics single-cell 3prime gene expression",
    chemistry_version = "Chromium Next GEM Single Cell 3prime Reagent Kit v3.1",
    feature_barcoding = "none_reported",
    vdj_capture = "no",
    read_len = NA_character_,
    assay_type = dplyr::coalesce(assay_type_sra, "RNA-Seq"),
    library_layout = dplyr::coalesce(library_layout_sra, "PAIRED"),
    library_selection = dplyr::coalesce(library_selection_sra, "cDNA"),
    library_source = dplyr::coalesce(library_source_sra, "TRANSCRIPTOMIC"),

    aligner = "Cell Ranger / user common atlas pipeline for this object; paper original used Cell Ranger 6.0.0",
    cellranger_or_equivalent_version = "atlas common reprocessing; paper provenance: Cell Ranger 6.0.0",
    paper_original_cellranger_version = "6.0.0",
    paper_original_reference_genome = "GRCh38",
    reference_genome = "GRCh38",
    gene_annotation = NA_character_,
    seq_depth = NA_character_,

    # Paper QC provenance.
    n_loaded_cells = NA_integer_,
    paper_loaded_cells_per_sample = NA_integer_,
    paper_qc_min_umi_per_cell = 500L,
    paper_qc_min_genes_per_cell = 200L,
    paper_qc_max_percent_mt = 20,
    paper_doubletfinder_expected_rate = 0.075,
    paper_doubletfinder_applied_here = FALSE,
    paper_cells_after_qc_total = 108597L,
    paper_hd_cells_after_qc = 30429L,
    paper_mild_cells_after_qc = 22188L,
    paper_severe_cells_after_qc = 55980L,
    paper_neutrophil_cells = 45463L,

    enrollment_window = "T1 within 72 h of hospital admission; T2 5-7 days later; healthy controls one sample; one late severe sample",
    batch_unit_recommended = "sample",
    source_note = paste(
      "GSE216020 is a COVID-19 neutrophil/whole blood scRNA-seq dataset.",
      "Metadata embedded from GEO, SRA Run Table, manuscript Methods/Results, and Supplementary Fig. S1B.",
      "Paper-side Cell Ranger/QC/DoubletFinder settings are provenance only; atlas matrices are from the user's common raw-data reprocessing."
    ),
    read_mode = NA_character_,
    predicted.cluster = NA_character_
  )

if (nrow(sample_meta) != 24) {
  stop("Expected 24 sample metadata rows for GSE216020, got ", nrow(sample_meta))
}

# Ensure all samples have subject metadata.
if (any(is.na(sample_meta$sex) | is.na(sample_meta$age) | is.na(sample_meta$race))) {
  bad <- sample_meta %>% filter(is.na(sex) | is.na(age) | is.na(race)) %>% pull(geo_accession)
  stop("Some samples did not receive subject-level sex/age/race metadata: ", paste(bad, collapse = ", "))
}

run_meta <- sample_meta

common_meta_columns <- c(
  "database_accession", "sra_study", "bioproject", "organism",
  "biosample_accession", "experiment_accession", "run_accession", "geo_accession",
  "cellranger_dir_id", "sample_id", "sample_alias", "paper_sample",
  "patient_id", "subject_id", "participant_id", "library_batch_id", "library_unit_id",
  "study_group", "outcome_label", "factor", "cov19", "is_sepsis", "ards", "pneumonia",
  "infection_site", "pathogen_etiology", "mortality", "final_outcome", "icu_admission",
  "days", "days_after_hospital_admission", "timepoint_label", "timepoint_note",
  "severity", "severityatday", "covid19_severity_numeric_paper", "ventilation",
  "sex", "age", "race", "qsofa", "apache_ii", "sofa",
  "source_name", "tissue_label", "source_blood_fraction", "cell_type", "cell_line",
  "extracted_molecule", "fresh_or_frozen", "blood_collection_tube",
  "pbmc_isolation_method", "granulocyte_isolation_method", "ficoll", "rbc_lysis",
  "sample_storage", "cell_suspension_buffer", "cell_mixture", "viability_reported",
  "sorted_population",
  "platform", "sequencer_model", "library_type", "chemistry_version",
  "feature_barcoding", "vdj_capture", "read_len", "assay_type",
  "library_layout", "library_selection", "library_source",
  "aligner", "cellranger_or_equivalent_version", "paper_original_cellranger_version",
  "paper_original_reference_genome", "reference_genome", "gene_annotation", "seq_depth",
  "n_loaded_cells", "paper_loaded_cells_per_sample",
  "paper_qc_min_umi_per_cell", "paper_qc_min_genes_per_cell",
  "paper_qc_max_percent_mt", "paper_doubletfinder_expected_rate",
  "paper_doubletfinder_applied_here",
  "paper_cells_after_qc_total", "paper_hd_cells_after_qc",
  "paper_mild_cells_after_qc", "paper_severe_cells_after_qc", "paper_neutrophil_cells",
  "n_sra_runs", "sample_name_sra", "library_name_sra",
  "bioproject_sra", "sra_study_sra", "avg_spot_len_sra", "bases_sra", "bytes_sra",
  "instrument_sra", "platform_sra", "consent_sra", "center_name_sra",
  "release_date_sra", "create_date_sra", "version_sra",
  "genotype_sra", "time_sra", "source_name_sra", "tissue_sra", "cell_type_sra",
  "datastore_filetype_sra", "datastore_provider_sra", "datastore_region_sra",
  "enrollment_window", "batch_unit_recommended", "source_note",
  "read_mode", "predicted.cluster"
)

sample_meta <- ensure_columns(sample_meta, common_meta_columns)
run_meta    <- ensure_columns(run_meta, common_meta_columns)

sample_meta <- atlas_standardize_metadata(sample_meta, project_id = project_id)
run_meta <- atlas_standardize_metadata(run_meta, project_id = project_id)
write_csv(sample_meta, file.path(meta_dir, "GSE216020_sample_metadata.csv"))
write_csv(run_meta,    file.path(meta_dir, "GSE216020_run_metadata.csv"))

# ------------------------------------------------------------------------------
# Load each BioSample-level matrix and build Seurat objects
# ------------------------------------------------------------------------------
input_dirs_all <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
input_dirs <- input_dirs_all[basename(input_dirs_all) %in% sample_meta$cellranger_dir_id]
input_dirs <- input_dirs[order(match(basename(input_dirs), sample_meta$cellranger_dir_id))]

if (length(input_dirs) == 0) {
  stop(
    "No Cell Ranger directories under input_root matched metadata.\n",
    "Expected BioSample dirs: ", paste(sample_meta$cellranger_dir_id, collapse = ", "), "\n",
    "Available under input_root: ", paste(head(basename(input_dirs_all), 50), collapse = ", ")
  )
}
if (length(input_dirs) != nrow(sample_meta)) {
  stop(
    "Matched ", length(input_dirs), " Cell Ranger directories but metadata has ",
    nrow(sample_meta), " rows. Missing: ",
    paste(setdiff(sample_meta$cellranger_dir_id, basename(input_dirs)), collapse = ", ")
  )
}

message("Sample metadata rows : ", nrow(sample_meta))
message("Matched BioSample dirs: ", length(input_dirs))

sample_meta <- atlas_standardize_metadata(sample_meta, project_id = project_id)

obj_list <- list()
qc_summary <- tibble(
  cellranger_dir_id = character(),
  run_accession = character(),
  biosample_accession = character(),
  geo_accession = character(),
  paper_sample = character(),
  participant_id = character(),
  study_group = character(),
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

for (run_dir in input_dirs) {
  dir_id <- basename(run_dir)
  raw_dir <- file.path(run_dir, "outs", "raw_feature_bc_matrix")
  if (!dir.exists(raw_dir)) {
    stop("Missing raw_feature_bc_matrix: ", raw_dir)
  }

  meta_row <- sample_meta %>% filter(cellranger_dir_id == dir_id)
  if (nrow(meta_row) != 1) {
    stop("Sample metadata row is not unique for cellranger_dir_id ", dir_id)
  }

  message("Reading ", dir_id, " -> ", meta_row$geo_accession, " (", meta_row$paper_sample, ")")
  # ATLAS_RESOURCE_SAFE_EXCLUSION_GSE216020
  res <- tryCatch(
    read_gene_expression(raw_dir, sample_id = dir_id),
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
        paper_sample = meta_row$paper_sample,
        participant_id = meta_row$participant_id,
        study_group = meta_row$study_group,
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

  meta_row$read_mode <- res$read_mode
  meta_row$predicted.cluster <- NA_character_
  obj <- attach_sample_metadata(obj, meta_row)

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- compute_percent_ribo(obj)
  obj[["percent.hb"]] <- compute_percent_hb(obj)

  obj[["paper_qc_pass"]] <- with(obj@meta.data,
    nCount_RNA >= paper_qc_min_umi_per_cell &
      nFeature_RNA >= paper_qc_min_genes_per_cell &
      percent.mt <= paper_qc_max_percent_mt
  )

  obj[["atlas_basic_qc_pass"]] <- with(obj@meta.data,
    nFeature_RNA >= 200
  )

  n_raw <- ncol(obj)
  obj <- atlas_standardize_seurat_metadata(obj, project_id = project_id)

  preqc_rds <- file.path(rds_preqc_dir, paste0(meta_row$geo_accession, "__", dir_id, "__preqc_raw.rds"))
  safe_save_rds(obj, preqc_rds)

  obj <- run_optional_helpers(obj)

  qc_subset <- rep(TRUE, ncol(obj))
  if (apply_paper_qc_filter) qc_subset <- qc_subset & obj$paper_qc_pass
  qc_subset <- qc_subset & obj$atlas_basic_qc_pass
  cells_keep <- colnames(obj)[which(qc_subset)]
  n_after_qc <- length(cells_keep)

  if (n_after_qc < min_cells_after_qc) {
    warning(meta_row$geo_accession, " dropped because only ", n_after_qc, " cells remained after QC.")
    qc_summary <- bind_rows(
      qc_summary,
      tibble(
        cellranger_dir_id = dir_id,
        run_accession = meta_row$run_accession,
        biosample_accession = meta_row$biosample_accession,
        geo_accession = meta_row$geo_accession,
        paper_sample = meta_row$paper_sample,
        participant_id = meta_row$participant_id,
        study_group = meta_row$study_group,
        n_raw = n_raw,
        n_after_qc = n_after_qc,
        saved = FALSE,
        preqc_rds_path = preqc_rds,
        qc_rds_path = NA_character_,
        working_rds_path = NA_character_
      )
    )
    next
  }

  obj_qc <- subset(obj, cells = cells_keep)

  qc_rds <- file.path(rds_qc_dir, paste0(meta_row$geo_accession, "__", dir_id, "__qcfiltered_raw.rds"))
  safe_save_rds(obj_qc, qc_rds)

  obj_work <- NormalizeData(obj_qc, verbose = FALSE)
  obj_work <- FindVariableFeatures(obj_work, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

  work_rds <- file.path(rds_work_dir, paste0(meta_row$geo_accession, "__", dir_id, "__working.rds"))
  safe_save_rds(obj_work, work_rds)

  obj_list[[as.character(meta_row$geo_accession)]] <- work_rds

  qc_summary <- bind_rows(
    qc_summary,
    tibble(
      cellranger_dir_id = dir_id,
      run_accession = meta_row$run_accession,
      biosample_accession = meta_row$biosample_accession,
      geo_accession = meta_row$geo_accession,
      paper_sample = meta_row$paper_sample,
      participant_id = meta_row$participant_id,
      study_group = meta_row$study_group,
      n_raw = n_raw,
      n_after_qc = n_after_qc,
      saved = TRUE,
      preqc_rds_path = preqc_rds,
      qc_rds_path = qc_rds,
      working_rds_path = work_rds
    )
  )
  # ATLAS_RESOURCE_SAFE_GC_GSE216020
  rm(res, obj, obj_qc, obj_work)
  invisible(gc())
}

write_csv(qc_summary, file.path(qc_dir, "GSE216020_qc_summary.csv"))
write_csv(
  qc_summary %>% dplyr::filter(preprocess_status %in% "excluded_before_rds"),
  file.path(qc_dir, "GSE216020_excluded_before_rds.csv")
)

registry <- qc_summary %>%
  filter(saved) %>%
  left_join(
    sample_meta %>%
      select(
        cellranger_dir_id, run_accession, biosample_accession, geo_accession,
        paper_sample, participant_id, study_group, factor, cov19, severity,
        timepoint_label, days_after_hospital_admission, final_outcome,
        mortality, batch_unit_recommended
      ),
    by = c("cellranger_dir_id", "run_accession", "biosample_accession", "geo_accession", "paper_sample", "participant_id", "study_group")
  ) %>%
  mutate(
    project_id = project_id,
    preintegration_tier = "working",
    object_grain = "sample"
  )

write_csv(registry, file.path(output_root, "GSE216020_preintegration_registry.csv"))
saveRDS(obj_list, file.path(output_root, "GSE216020_preintegration_objlist.rds"))

message("Done. Saved pre-QC raw objects to    : ", rds_preqc_dir)
message("Done. Saved QC-filtered raw objects to: ", rds_qc_dir)
message("Done. Saved working objects to       : ", rds_work_dir)
message("Registry: ", file.path(output_root, "GSE216020_preintegration_registry.csv"))
