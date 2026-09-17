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
# GSE217906 pre-integration loader
#
# Dataset:
#   GSE217906 / SRP407793 / PRJNA901378
#   Immune-cell signatures of persistent inflammation, immunosuppression,
#   and catabolism syndrome after sepsis
#
# Atlas-side decision:
#   Public raw reads are assumed to have been reprocessed by the common atlas
#   Cell Ranger pipeline. Paper-side Cell Ranger/Seurat/QC settings are kept
#   as provenance metadata only.
#
# Expected external input:
#   configs/GSE217906/GSE217906_SraRunTable.csv
#   or configs/GSE217906/SraRunTable.csv
#
# Expected Cell Ranger output root:
#   cellranger_count/output/SRP407793/<BioSample>/outs/raw_feature_bc_matrix
#
# Output follows the same pre-integration contract as previous atlas loaders:
#   rds_preqc_raw/
#   rds_qcfiltered_raw/
#   rds_working/
#   metadata/
#   manifest/
#   qc/
# ==============================================================================

project_id <- "GSE217906"
sra_study  <- "SRP407793"
bioproject <- "PRJNA901378"
organism   <- "Homo sapiens"

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

normalize_names <- function(x) {
  x <- make.unique(x, sep = "_dup")
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "") %>%
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

resolve_input_root <- function(project_root, sra_study) {
  candidates <- c(
    Sys.getenv("CELLRANGER_COUNT_ROOT", unset = ""),
    file.path(project_root, "cellranger_count", "output", sra_study),
    file.path(Sys.getenv("ATLAS_CELLRANGER_OUTPUT_ROOT", unset = ""), sra_study),
    file.path(Sys.getenv("ATLAS_CELLRANGER_OUTPUT_ROOT", unset = ""), sra_study)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hits <- candidates[dir.exists(candidates)]
  if (length(hits) == 0) {
    stop("Could not find Cell Ranger output directory. Tried:\n",
         paste0(" - ", candidates, collapse = "\n"))
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
    stop("Missing SRA run table. Tried:\n", paste0(" - ", candidates, collapse = "\n"))
  }
  normalizePath(hits[1], winslash = "/", mustWork = TRUE)
}

ensure_columns <- function(df, wanted) {
  miss <- setdiff(wanted, names(df))
  for (m in miss) df[[m]] <- NA
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

for (d in c(output_root, rds_preqc_dir, rds_qc_dir, rds_work_dir, meta_dir, qc_dir, manifest_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

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
# Embedded GEO sample metadata
# GEO gives 16 samples. First 9 use paper-style IDs; 7 later GEO titles are
# group labels and are conservatively crosswalked to the remaining Table S1
# rows with a "needs_manual_verification" flag.
# ------------------------------------------------------------------------------
geo_sample_tbl <- tibble::tribble(
  ~geo_accession, ~geo_title, ~paper_sample, ~geo_disease_state, ~clinical_id_table_s1, ~clinical_mapping_status,
  "GSM6729707", "PD-01",   "PD-01",   "PICS-death", "PD-01", "exact_geo_title_to_table_s1",
  "GSM6729708", "PD-02",   "PD-02",   "PICS-death", "PD-02", "exact_geo_title_to_table_s1",
  "GSM6729709", "PA-01",   "PA-01",   "PICS-alive", "PA-01", "exact_geo_title_to_table_s1",
  "GSM6729710", "PA-02",   "PA-02",   "PICS-alive", "PA-02", "exact_geo_title_to_table_s1",
  "GSM6729711", "AS-01",   "AS-01",   "Sepsis",     "AS-01", "exact_geo_title_to_table_s1",
  "GSM6729712", "AS-02",   "AS-02",   "Sepsis",     "AS-02", "exact_geo_title_to_table_s1",
  "GSM6729713", "HC-01",   "HC-01",   "HC",         "HC-01", "exact_geo_title_to_table_s1",
  "GSM6729714", "HC-02",   "HC-02",   "HC",         "HC-02", "exact_geo_title_to_table_s1",
  "GSM6729715", "HC-03",   "HC-03",   "HC",         "HC-03", "exact_geo_title_to_table_s1",
  "GSM8217322", "PICS5",   "PICS5",   "PICS-alive", NA_character_, "group_only_sra_disease_state__clinical_id_unresolved",
  "GSM8217323", "Sepsis4", "Sepsis4", "Sepsis",     NA_character_, "group_only_sra_disease_state__clinical_id_unresolved",
  "GSM8217324", "Sepsis2", "Sepsis2", "Sepsis",     NA_character_, "group_only_sra_disease_state__clinical_id_unresolved",
  "GSM8217325", "Sepsis3", "Sepsis3", "Sepsis",     NA_character_, "group_only_sra_disease_state__clinical_id_unresolved",
  "GSM8217326", "HC3",     "HC3",     "HC",         NA_character_, "group_only_sra_disease_state__clinical_id_unresolved",
  "GSM8217327", "HC5",     "HC5",     "HC",         NA_character_, "group_only_sra_disease_state__clinical_id_unresolved",
  "GSM8217328", "PICS3",   "PICS3",   "PICS-alive", NA_character_, "group_only_sra_disease_state__clinical_id_unresolved"
)

# Embedded Supplementary Table S1 clinical characteristics.
clinical_tbl <- tibble::tribble(
  ~clinical_id, ~age, ~sex, ~diagnosis, ~past_medical_history_comorbidities, ~icu_stay_day, ~crp_mg_l, ~albumin_g_l, ~prealbumin_g_l, ~lymphocyte_10e9_l, ~creatinine_height_index, ~retinol_binding_protein_mg_dl, ~weight_loss_or_bmi_note, ~prognosis_28d, ~sofa,
  "PD-01", 82, "male",   "Pneumonia, PICS secondary to sepsis", "Four months after surgery for intracranial tumor, and hypertension for 15 years", 17, 65.9, 29, 0.08, 0.5, 0.71, 0.6, "Weight loss=11.5% BMI=23.1", "death", 3,
  "PD-02", 74, "male",   "Pneumonia, Acute respiratory distress syndrome (ARDS), PICS secondary to sepsis", "One month after surgery for gastric cancer, hypertension for 10 years and diabetes for one year", 21, 208, 28, 0.08, 0.7, 0.72, 0.7, "Weight loss=12.4% BMI=17.4", "death", 4,
  "PA-01", 81, "female", "Pneumonia, PICS secondary to sepsis", "Three months after appendectomy, diabetes for 10 years, hypertension and atrial fibrillation for 15 years", 28, 85.7, 22, 0.08, 0.7, 0.76, 0.78, "Weight loss=13.1% BMI=17.9", "alive", 2,
  "PA-02", 65, "male",   "Septic shock, PICS secondary to sepsis", "Half a month after surgery for gastric cancer", 18, 220, 24.5, 0.06, 0.7, 0.75, 0.89, "Weight loss=13.8% BMI=17.7", "alive", 3,
  "PA-03", 77, "male",   "Pneumonia, PICS secondary to sepsis", "Two months after hernia repair surgery, hypertension and atrial fibrillation for 30 years", 19, 67.7, 25, 0.08, 0.1, 0.75, 0.86, "Weight loss=11.4% BMI=15.5", "alive", 2,
  "PA-04", 71, "male",   "Septic shock, PICS secondary to sepsis", "Half a month after surgery for postoperative esophageal cancer, hypertension 12 years", 19, 192, 26.1, 0.05, 0.3, 0.8, 0.78, "Weight loss=21.0% BMI=14.7", "alive", 2,
  "AS-01", 80, "female", "Acute cholecystitis, laparotomy, sepsis", "15 days after exploratory laparotomy, hypertension 15 years", 3, 92.4, 30, 0.08, 1.4, 0.83, 0.75, "BMI=25.1", "alive", 3,
  "AS-02", 67, "male",   "Postoperative anastomotic leakage in colon cancer, sepsis", "10 days after colon cancer surgery, hypertension for five years and atrial fibrillation for three years", 3, 234.7, 31.7, 0.11, 0.3, 0.83, 23, "BMI=25.4", "alive", 3,
  "AS-03", 75, "male",   "Gastric perforation, sepsis", "Postoperative history of lung cancer, hypertension for 30 years, coronary heart disease for 10 years, and cerebral infarction for 20 years", 2, 759.4, 30, 0.08, 1.1, 0.75, 12, "BMI=23.0", "alive", 3,
  "AS-04", 55, "male",   "Abdominal infection, pneumonia, sepsis", "10 days after appendectomy, pulmonary tuberculosis for 40 years", 2, 396, 29, 0.08, 1.3, 0.82, 6, "BMI=24.8", "alive", 4,
  "AS-05", 68, "male",   "Diffuse peritonitis, septic shock, sepsis", "Anastomotic leakage after rectal cancer surgery, hypertension 5 years and hepatitis B", 3, 240.4, 26.5, 0.08, 0.8, 0.85, 0.95, "BMI=17.1", "alive", 4,
  "HC-01", 65, "female", "HC", "Hypertension 1 year", 0, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, "BMI=20.2", "alive", 0,
  "HC-02", 72, "male",   "HC", "Hypertension 12 years", 0, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, "BMI=23.3", "alive", 0,
  "HC-03", 88, "female", "HC", "Hypertension 10 years", 0, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, "BMI=22.1", "alive", 0,
  "HC-04", 69, "male",   "HC", "Hypertension 3 years", 0, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, "BMI=21.7", "alive", 0,
  "HC-05", 70, "female", "HC", "None", 0, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, "BMI=20.8", "alive", 0
)

paper_methods_tbl <- tibble::tribble(
  ~field, ~value, ~note,
  "paper_title", "Immune-cell signatures of persistent inflammation, immunosuppression, and catabolism syndrome after sepsis", "Med 2025",
  "cohort_pics_n", "6", "six PICS secondary to sepsis participants",
  "cohort_acute_sepsis_n", "5", "five acute sepsis participants",
  "cohort_healthy_control_n", "5", "five healthy controls",
  "pics_death_n", "2", "PD-01 and PD-02",
  "pics_alive_n", "4", "PA-01 to PA-04",
  "high_quality_cells_reported", "91180", "paper Results",
  "participants_age_range_requirement", "all participants age >=65", "paper Methods",
  "sepsis_definition", "Sepsis 3.0; diagnosis within 24 h; infection + SOFA organ dysfunction", "paper Methods",
  "pics_definition", "ICU >14 days plus inflammation, immunosuppression, and catabolism criteria", "paper Methods",
  "sample_timing", "blood samples within 24 h after diagnosis of sepsis or PICS", "paper Methods",
  "source_material", "peripheral blood mononuclear cells", "paper/GEO",
  "pbmc_isolation", "whole blood diluted 1:1 with PBS, Ficoll-Paque Plus density gradient at 400g for 30 min, DPBS washes, red blood cell lysis", "paper Methods",
  "cell_viability", "approximately 90% for each sample", "paper Methods",
  "geo_library_protocol", "10x Genomics Single Cell 3prime V3.1 according to GEO sample pages", "GEO provenance",
  "paper_library_protocol", "10x Chromium Single Cell 5prime library and V(D)J enrichment kit", "paper Methods; keep as provenance",
  "sequencer", "Illumina NovaSeq 6000", "GEO/Paper",
  "paper_original_cellranger", "CellRanger v3.1.0", "paper Methods/GEO",
  "paper_original_reference", "GRCh38 / GRCh38 Ensembl v91", "paper Methods/GEO",
  "paper_original_seurat", "Seurat v3.1.4", "paper Methods/GEO",
  "paper_qc_min_genes_per_cell", "200", "GEO data processing",
  "paper_qc_max_percent_mt", "20", "GEO data processing",
  "paper_note_common_reprocessing", "Atlas objects are generated from user common Cell Ranger outputs; paper processing settings are provenance only.", "Atlas policy"
)

write_csv(geo_sample_tbl, file.path(manifest_dir, "GSE217906_geo_sample_embedded.csv"))
write_csv(clinical_tbl, file.path(manifest_dir, "GSE217906_clinical_table_s1_embedded.csv"))
write_csv(paper_methods_tbl, file.path(manifest_dir, "GSE217906_paper_methods_and_results_embedded.csv"))

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
disease_state_col <- pick_column(sra_raw, exact = c("disease_state"), regex = "^disease_state$", required = TRUE, label = "SRA run table")
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
  disease_state_sra = col_or_na(sra_raw, disease_state_col),
  time_sra = col_or_na(sra_raw, time_col),
  tissue_sra = col_or_na(sra_raw, tissue_col),
  datastore_filetype_sra = col_or_na(sra_raw, filetype_col),
  datastore_provider_sra = col_or_na(sra_raw, provider_col),
  datastore_region_sra = col_or_na(sra_raw, region_col)
) %>%
  mutate(across(where(is.character), ~na_if(.x, "")))

write_csv(sra_norm, file.path(manifest_dir, "GSE217906_sra_table_normalized.csv"))

if (!any(sra_norm$bioproject_sra == bioproject, na.rm = TRUE)) {
  stop("SRA run table does not contain expected BioProject ", bioproject,
       ". You may have supplied a SRA table for another dataset.")
}
if (!any(sra_norm$sra_study_sra == sra_study, na.rm = TRUE)) {
  stop("SRA run table does not contain expected SRA Study ", sra_study,
       ". You may have supplied a SRA table for another dataset.")
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
    disease_state_sra = paste(unique(na.omit(disease_state_sra)), collapse = ";"),
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
  stop("SRA-to-GEO matching is incomplete for GSE217906.\n",
       "Matched GEO: ", paste(sra_sample_agg$geo_accession, collapse = ", "), "\n",
       "Missing GEO: ", paste(missing_geo, collapse = ", "), "\n",
       "Expected exactly 16 GEO samples.")
}

write_csv(sra_sample_agg, file.path(manifest_dir, "GSE217906_sra_sample_aggregated.csv"))

# ------------------------------------------------------------------------------
# Build sample metadata
# ------------------------------------------------------------------------------
sample_meta <- geo_sample_tbl %>%
  left_join(sra_sample_agg, by = "geo_accession") %>%
  left_join(clinical_tbl, by = c("clinical_id_table_s1" = "clinical_id")) %>%
  mutate(
    database_accession = project_id,
    sra_study = sra_study,
    bioproject = bioproject,
    organism = organism,

    disease_state_sra_primary = disease_state_sra,
    disease_state_norm = case_when(
      str_detect(str_to_lower(disease_state_sra_primary), "pics[-_ ]?death") ~ "PICS_death",
      str_detect(str_to_lower(disease_state_sra_primary), "pics[-_ ]?alive") ~ "PICS_alive",
      str_to_lower(disease_state_sra_primary) %in% c("sepsis", "acute sepsis", "acute_sepsis") ~ "acute_sepsis",
      str_to_lower(disease_state_sra_primary) %in% c("hc", "healthy control", "healthy_control") ~ "healthy_control",
      TRUE ~ NA_character_
    ),

    cellranger_dir_id = biosample_accession,
    sample_id = geo_accession,
    sample_alias = geo_title,

    # clinical_id is only populated when GEO title exactly matches Table S1 ID.
    # For later GEO update samples (PICS5/PICS3/Sepsis2/3/4/HC3/HC5), group labels
    # are taken from SRA disease_state but individual clinical covariates remain NA.
    clinical_id = clinical_id_table_s1,
    clinical_id_source = ifelse(
      clinical_mapping_status == "exact_geo_title_to_table_s1",
      "GEO title exactly matches Table S1 participant ID",
      "unresolved; SRA disease_state gives group only"
    ),
    clinical_match_confidence = ifelse(
      clinical_mapping_status == "exact_geo_title_to_table_s1",
      "exact_individual_clinical_match",
      "group_only_individual_clinical_covariates_missing"
    ),
    patient_id = ifelse(!is.na(clinical_id), clinical_id, paste0("unresolved_", geo_accession)),
    subject_id = patient_id,
    participant_id = patient_id,

    study_group = disease_state_norm,
    factor = case_when(
      disease_state_norm %in% c("PICS_death", "PICS_alive") ~ "PICS",
      disease_state_norm == "acute_sepsis" ~ "sepsis",
      disease_state_norm == "healthy_control" ~ "healthy_control",
      TRUE ~ NA_character_
    ),
    cov19 = "not_applicable",
    is_sepsis = case_when(
      factor %in% c("PICS", "sepsis") ~ "yes",
      factor == "healthy_control" ~ "no",
      TRUE ~ NA_character_
    ),
    ards = case_when(
      !is.na(diagnosis) & str_detect(diagnosis, "ARDS") ~ "yes",
      !is.na(diagnosis) ~ "not_reported_or_no",
      TRUE ~ NA_character_
    ),
    pneumonia = case_when(
      !is.na(diagnosis) & str_detect(str_to_lower(diagnosis), "pneumonia") ~ "yes",
      !is.na(diagnosis) ~ "not_reported_or_no",
      TRUE ~ NA_character_
    ),
    infection_site = case_when(
      !is.na(diagnosis) & str_detect(str_to_lower(diagnosis), "pneumonia") ~ "respiratory",
      !is.na(diagnosis) & str_detect(str_to_lower(diagnosis), "abdominal|peritonitis|perforation|cholecystitis|anastomotic") ~ "abdominal",
      factor == "healthy_control" ~ "none",
      factor %in% c("PICS", "sepsis") ~ "not_resolvable_sample_level",
      TRUE ~ NA_character_
    ),
    pathogen_etiology = case_when(
      factor %in% c("PICS", "sepsis") ~ "gram-negative bacterial infection",
      factor == "healthy_control" ~ "none",
      TRUE ~ NA_character_
    ),
    mortality = case_when(
      disease_state_norm == "PICS_death" ~ "yes",
      disease_state_norm %in% c("PICS_alive", "acute_sepsis", "healthy_control") ~ "no",
      TRUE ~ NA_character_
    ),
    final_outcome = case_when(
      disease_state_norm == "PICS_death" ~ "death",
      disease_state_norm %in% c("PICS_alive", "acute_sepsis") ~ "alive",
      disease_state_norm == "healthy_control" ~ "healthy_control",
      TRUE ~ NA_character_
    ),
    icu_admission = case_when(
      factor %in% c("PICS", "sepsis") ~ "yes",
      factor == "healthy_control" ~ "no",
      TRUE ~ NA_character_
    ),
    days = NA_real_,
    timepoint_label = "diagnosis_window_24h",
    severity = study_group,
    severityatday = paste0(severity, "_", timepoint_label),
    ventilation = NA_character_,
    qsofa = NA_real_,
    apache_ii = NA_real_,

    pics_criteria_icu_gt14d = case_when(
      factor == "PICS" & !is.na(clinical_id) ~ "yes",
      factor == "PICS" & is.na(clinical_id) ~ "unknown_individual_clinical_mapping",
      TRUE ~ "not_applicable"
    ),
    pics_criteria_crp_gt150 = case_when(
      factor == "PICS" & !is.na(crp_mg_l) & crp_mg_l > 150 ~ "yes",
      factor == "PICS" & !is.na(crp_mg_l) ~ "no",
      factor == "PICS" ~ "unknown_individual_clinical_mapping",
      TRUE ~ "not_applicable"
    ),
    pics_criteria_lym_lt0_8 = case_when(
      factor == "PICS" & !is.na(lymphocyte_10e9_l) & lymphocyte_10e9_l < 0.8 ~ "yes",
      factor == "PICS" & !is.na(lymphocyte_10e9_l) ~ "no",
      factor == "PICS" ~ "unknown_individual_clinical_mapping",
      TRUE ~ "not_applicable"
    ),
    pics_criteria_albumin_lt30 = case_when(
      factor == "PICS" & !is.na(albumin_g_l) & albumin_g_l < 30 ~ "yes",
      factor == "PICS" & !is.na(albumin_g_l) ~ "no",
      factor == "PICS" ~ "unknown_individual_clinical_mapping",
      TRUE ~ "not_applicable"
    ),
    pics_criteria_prealbumin_lt0_1 = case_when(
      factor == "PICS" & !is.na(prealbumin_g_l) & prealbumin_g_l < 0.1 ~ "yes",
      factor == "PICS" & !is.na(prealbumin_g_l) ~ "no",
      factor == "PICS" ~ "unknown_individual_clinical_mapping",
      TRUE ~ "not_applicable"
    ),
    pics_criteria_rbp_lt1 = case_when(
      factor == "PICS" & !is.na(retinol_binding_protein_mg_dl) & retinol_binding_protein_mg_dl < 1 ~ "yes",
      factor == "PICS" & !is.na(retinol_binding_protein_mg_dl) ~ "no",
      factor == "PICS" ~ "unknown_individual_clinical_mapping",
      TRUE ~ "not_applicable"
    ),
    pics_criteria_chi_lt80pct = case_when(
      factor == "PICS" & !is.na(creatinine_height_index) & creatinine_height_index < 0.8 ~ "yes",
      factor == "PICS" & !is.na(creatinine_height_index) ~ "no",
      factor == "PICS" ~ "unknown_individual_clinical_mapping",
      TRUE ~ "not_applicable"
    ),
    pics_criteria_weightloss_or_bmi = case_when(
      factor == "PICS" & !is.na(weight_loss_or_bmi_note) & str_detect(str_to_lower(weight_loss_or_bmi_note), "weight loss|bmi=1[0-7]|bmi = 1[0-7]") ~ "yes_or_reported",
      factor == "PICS" & !is.na(weight_loss_or_bmi_note) ~ "not_reported_or_no",
      factor == "PICS" ~ "unknown_individual_clinical_mapping",
      TRUE ~ "not_applicable"
    ),

    source_name = coalesce(source_name_sra, "Blood"),
    tissue_label = coalesce(tissue_sra, "Blood"),
    source_blood_fraction = "PBMC",
    cell_type = coalesce(cell_type_sra, "peripheral blood mononuclear cells"),
    cell_line = coalesce(cell_line_sra, "Not a cell line"),
    extracted_molecule = "total RNA",
    frozen_or_fresh = "fresh PBMC according to GEO/paper processing",
    fresh_or_frozen = frozen_or_fresh,
    blood_collection_tube = NA_character_,
    pbmc_isolation_method = "Ficoll-Paque Plus density-gradient centrifugation; RBC lysis; PBS/0.04% BSA wash/filter",
    ficoll = "Ficoll-Paque Plus",
    rbc_lysis = "yes; red blood cell lysis buffer",
    sample_storage = "fresh processing for single-cell library",
    cell_suspension_buffer = "PBS containing 0.04% BSA",
    viability_reported = "approximately 90%",
    sorted_population = "PBMC",

    platform = "10x Genomics",
    sequencer_model = coalesce(instrument_sra, "Illumina NovaSeq 6000"),
    library_type = "10x Genomics single-cell gene expression",
    chemistry_version = "GEO: Single Cell 3prime V3.1; paper: Single Cell 5prime plus V(D)J provenance",
    feature_barcoding = "V(D)J libraries reported in paper; RNA object only here",
    vdj_capture = "reported in paper; not loaded here",
    read_len = "150 bp paired-end in paper",
    assay_type = coalesce(assay_type_sra, "RNA-Seq"),
    library_layout = coalesce(library_layout_sra, "PAIRED"),
    library_selection = coalesce(library_selection_sra, "cDNA"),
    library_source = coalesce(library_source_sra, "TRANSCRIPTOMIC SINGLE CELL"),

    aligner = "Cell Ranger / user common atlas pipeline for this object; paper original used CellRanger v3.1.0",
    cellranger_or_equivalent_version = "atlas common reprocessing; paper provenance: CellRanger v3.1.0",
    paper_original_cellranger_version = "3.1.0",
    paper_original_reference_genome = "GRCh38 Ensembl v91",
    reference_genome = "GRCh38",
    gene_annotation = "GRCh38 Ensembl v91 in paper provenance",
    seq_depth = NA_character_,

    n_loaded_cells = NA_integer_,
    paper_cells_after_qc_total = 91180L,
    paper_qc_min_genes_per_cell = 200L,
    paper_qc_max_percent_mt = 20,
    paper_mitochondria_genes_removed_from_expression_table = TRUE,

    enrollment_window = "blood samples within 24 h after diagnosis of sepsis or PICS; acute sepsis diagnosed within 24 h",
    batch_unit_recommended = "sample",
    source_note = paste(
      "GSE217906 PBMC scRNA-seq dataset with PICS, acute sepsis, and healthy controls.",
      "Group labels are assigned primarily from SRA disease_state.",
      "Only samples whose GEO title exactly matches Supplementary Table S1 IDs receive individual clinical covariates.",
      "Later GEO update samples with titles PICS5/PICS3/Sepsis2/3/4/HC3/HC5 retain group labels but individual clinical covariates are intentionally NA."
    ),
    read_mode = NA_character_,
    predicted.cluster = NA_character_
  )

if (nrow(sample_meta) != 16) {
  stop("Expected 16 sample metadata rows for GSE217906, got ", nrow(sample_meta))
}

if (any(is.na(sample_meta$disease_state_sra_primary) | is.na(sample_meta$study_group) | is.na(sample_meta$factor))) {
  bad <- sample_meta %>% filter(is.na(disease_state_sra_primary) | is.na(study_group) | is.na(factor)) %>% pull(geo_accession)
  stop("Some samples did not receive SRA-derived group metadata: ", paste(bad, collapse = ", "))
}

exact_rows <- sample_meta %>% filter(clinical_mapping_status == "exact_geo_title_to_table_s1")
if (nrow(exact_rows) != 9) {
  stop("Expected 9 exact Table S1 clinical matches, got ", nrow(exact_rows))
}
if (any(is.na(exact_rows$sex) | is.na(exact_rows$age) | is.na(exact_rows$sofa))) {
  bad <- exact_rows %>% filter(is.na(sex) | is.na(age) | is.na(sofa)) %>% pull(geo_accession)
  stop("Some exact clinical matches did not receive sex/age/SOFA metadata: ", paste(bad, collapse = ", "))
}

unresolved_rows <- sample_meta %>% filter(clinical_mapping_status == "group_only_sra_disease_state__clinical_id_unresolved")
if (nrow(unresolved_rows) != 7) {
  stop("Expected 7 group-only unresolved clinical matches, got ", nrow(unresolved_rows))
}
if (any(!is.na(unresolved_rows$clinical_id) | !is.na(unresolved_rows$age) | !is.na(unresolved_rows$sex) | !is.na(unresolved_rows$sofa))) {
  bad <- unresolved_rows %>% filter(!is.na(clinical_id) | !is.na(age) | !is.na(sex) | !is.na(sofa)) %>% pull(geo_accession)
  stop("Unresolved samples unexpectedly received individual clinical metadata: ", paste(bad, collapse = ", "))
}

common_meta_columns <- c(
  "database_accession", "sra_study", "bioproject", "organism",
  "biosample_accession", "experiment_accession", "run_accession", "geo_accession",
  "cellranger_dir_id", "sample_id", "sample_alias", "paper_sample", "geo_title",
  "patient_id", "subject_id", "participant_id", "clinical_id", "clinical_id_source", "clinical_match_confidence", "clinical_mapping_status",
  "disease_state_sra_primary", "disease_state_norm",
  "study_group", "factor", "cov19", "is_sepsis", "ards", "pneumonia",
  "infection_site", "pathogen_etiology", "mortality", "final_outcome", "icu_admission",
  "days", "timepoint_label", "severity", "severityatday", "ventilation",
  "sex", "age", "qsofa", "apache_ii", "sofa",
  "diagnosis", "past_medical_history_comorbidities", "icu_stay_day",
  "crp_mg_l", "albumin_g_l", "prealbumin_g_l", "lymphocyte_10e9_l",
  "creatinine_height_index", "retinol_binding_protein_mg_dl", "weight_loss_or_bmi_note",
  "pics_criteria_icu_gt14d", "pics_criteria_crp_gt150", "pics_criteria_lym_lt0_8",
  "pics_criteria_albumin_lt30", "pics_criteria_prealbumin_lt0_1",
  "pics_criteria_rbp_lt1", "pics_criteria_chi_lt80pct", "pics_criteria_weightloss_or_bmi",
  "source_name", "tissue_label", "source_blood_fraction", "cell_type", "cell_line",
  "extracted_molecule", "frozen_or_fresh", "fresh_or_frozen", "blood_collection_tube",
  "pbmc_isolation_method", "ficoll", "rbc_lysis", "sample_storage",
  "cell_suspension_buffer", "viability_reported", "sorted_population",
  "platform", "sequencer_model", "library_type", "chemistry_version",
  "feature_barcoding", "vdj_capture", "read_len", "assay_type",
  "library_layout", "library_selection", "library_source",
  "aligner", "cellranger_or_equivalent_version", "paper_original_cellranger_version",
  "paper_original_reference_genome", "reference_genome", "gene_annotation", "seq_depth",
  "n_loaded_cells", "paper_cells_after_qc_total",
  "paper_qc_min_genes_per_cell", "paper_qc_max_percent_mt",
  "paper_mitochondria_genes_removed_from_expression_table",
  "n_sra_runs", "sample_name_sra", "library_name_sra",
  "bioproject_sra", "sra_study_sra", "avg_spot_len_sra", "bases_sra", "bytes_sra",
  "instrument_sra", "platform_sra", "consent_sra", "center_name_sra",
  "release_date_sra", "create_date_sra", "version_sra",
  "genotype_sra", "disease_state_sra", "time_sra", "source_name_sra", "tissue_sra", "cell_type_sra",
  "datastore_filetype_sra", "datastore_provider_sra", "datastore_region_sra",
  "enrollment_window", "batch_unit_recommended", "source_note",
  "read_mode", "predicted.cluster"
)

sample_meta <- ensure_columns(sample_meta, common_meta_columns)
run_meta <- sample_meta

sample_meta <- atlas_standardize_metadata(sample_meta, project_id = project_id)
run_meta <- atlas_standardize_metadata(run_meta, project_id = project_id)
write_csv(sample_meta, file.path(meta_dir, "GSE217906_sample_metadata.csv"))
write_csv(run_meta,    file.path(meta_dir, "GSE217906_run_metadata.csv"))

# ------------------------------------------------------------------------------
# Load each BioSample-level matrix and build Seurat objects
# ------------------------------------------------------------------------------
input_dirs_all <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
input_dirs <- input_dirs_all[basename(input_dirs_all) %in% sample_meta$cellranger_dir_id]
input_dirs <- input_dirs[order(match(basename(input_dirs), sample_meta$cellranger_dir_id))]

if (length(input_dirs) == 0) {
  stop("No Cell Ranger directories under input_root matched metadata.\n",
       "Expected BioSample dirs: ", paste(sample_meta$cellranger_dir_id, collapse = ", "), "\n",
       "Available under input_root: ", paste(head(basename(input_dirs_all), 50), collapse = ", "))
}
if (length(input_dirs) != nrow(sample_meta)) {
  stop("Matched ", length(input_dirs), " Cell Ranger directories but metadata has ",
       nrow(sample_meta), " rows. Missing: ",
       paste(setdiff(sample_meta$cellranger_dir_id, basename(input_dirs)), collapse = ", "))
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
  clinical_id = character(),
  study_group = character(),
  factor = character(),
  n_raw = integer(),
  n_after_qc = integer(),
  saved = logical(),
  preqc_rds_path = character(),
  qc_rds_path = character(),
  working_rds_path = character()
)

for (run_dir in input_dirs) {
  dir_id <- basename(run_dir)
  raw_dir <- file.path(run_dir, "outs", "raw_feature_bc_matrix")
  if (!dir.exists(raw_dir)) stop("Missing raw_feature_bc_matrix: ", raw_dir)

  meta_row <- sample_meta %>% filter(cellranger_dir_id == dir_id)
  if (nrow(meta_row) != 1) stop("Sample metadata row is not unique for cellranger_dir_id ", dir_id)

  message("Reading ", dir_id, " -> ", meta_row$geo_accession, " (", meta_row$paper_sample, ")")
  res <- read_gene_expression(raw_dir, sample_id = dir_id)

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
    nFeature_RNA >= paper_qc_min_genes_per_cell &
      percent.mt <= paper_qc_max_percent_mt
  )

  obj[["atlas_basic_qc_pass"]] <- with(obj@meta.data, nFeature_RNA >= 200)

  n_raw <- ncol(obj)
  obj <- atlas_standardize_seurat_metadata(obj, project_id = project_id)

  preqc_rds <- file.path(rds_preqc_dir, paste0(meta_row$geo_accession, "__", dir_id, "__preqc_raw.rds"))
  safe_save_rds(obj, preqc_rds)

  qc_subset <- rep(TRUE, ncol(obj))
  if (apply_paper_qc_filter) qc_subset <- qc_subset & obj$paper_qc_pass
  qc_subset <- qc_subset & obj$atlas_basic_qc_pass
  cells_keep <- colnames(obj)[which(qc_subset)]
  n_after_qc <- length(cells_keep)

  if (n_after_qc < min_cells_after_qc) {
    warning(meta_row$geo_accession, " dropped because only ", n_after_qc, " cells remained after QC.")
    qc_summary <- bind_rows(qc_summary, tibble(
      cellranger_dir_id = dir_id,
      run_accession = meta_row$run_accession,
      biosample_accession = meta_row$biosample_accession,
      geo_accession = meta_row$geo_accession,
      paper_sample = meta_row$paper_sample,
      clinical_id = meta_row$clinical_id,
      study_group = meta_row$study_group,
      factor = meta_row$factor,
      n_raw = n_raw,
      n_after_qc = n_after_qc,
      saved = FALSE,
      preqc_rds_path = preqc_rds,
      qc_rds_path = NA_character_,
      working_rds_path = NA_character_
    ))
    next
  }

  obj_qc <- subset(obj, cells = cells_keep)

  qc_rds <- file.path(rds_qc_dir, paste0(meta_row$geo_accession, "__", dir_id, "__qcfiltered_raw.rds"))
  safe_save_rds(obj_qc, qc_rds)

  obj_work <- NormalizeData(obj_qc, verbose = FALSE)
  obj_work <- FindVariableFeatures(obj_work, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

  work_rds <- file.path(rds_work_dir, paste0(meta_row$geo_accession, "__", dir_id, "__working.rds"))
  safe_save_rds(obj_work, work_rds)

  obj_list[[as.character(meta_row$geo_accession)]] <- obj_work

  qc_summary <- bind_rows(qc_summary, tibble(
    cellranger_dir_id = dir_id,
    run_accession = meta_row$run_accession,
    biosample_accession = meta_row$biosample_accession,
    geo_accession = meta_row$geo_accession,
    paper_sample = meta_row$paper_sample,
    clinical_id = meta_row$clinical_id,
    study_group = meta_row$study_group,
    factor = meta_row$factor,
    n_raw = n_raw,
    n_after_qc = n_after_qc,
    saved = TRUE,
    preqc_rds_path = preqc_rds,
    qc_rds_path = qc_rds,
    working_rds_path = work_rds
  ))
}

write_csv(qc_summary, file.path(qc_dir, "GSE217906_qc_summary.csv"))

registry <- qc_summary %>%
  filter(saved) %>%
  left_join(
    sample_meta %>%
      select(
        cellranger_dir_id, run_accession, biosample_accession, geo_accession,
        paper_sample, clinical_id, study_group, factor, is_sepsis,
        severity, final_outcome, mortality, clinical_mapping_status,
        batch_unit_recommended
      ),
    by = c("cellranger_dir_id", "run_accession", "biosample_accession", "geo_accession", "paper_sample", "clinical_id", "study_group", "factor")
  ) %>%
  mutate(
    project_id = project_id,
    preintegration_tier = "working",
    object_grain = "sample"
  )

write_csv(registry, file.path(output_root, "GSE217906_preintegration_registry.csv"))
saveRDS(obj_list, file.path(output_root, "GSE217906_preintegration_objlist.rds"))

message("Done. Saved pre-QC raw objects to    : ", rds_preqc_dir)
message("Done. Saved QC-filtered raw objects to: ", rds_qc_dir)
message("Done. Saved working objects to       : ", rds_work_dir)
message("Registry: ", file.path(output_root, "GSE217906_preintegration_registry.csv"))
