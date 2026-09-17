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

project_id <- "GSE220189"
parent_series <- "GSE220190"
sra_study <- "SRP411630"
bioproject <- "PRJNA909218"
organism <- "Homo sapiens"
expected_n <- 44L

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1 && nzchar(args[[1]])) {
  PROJECT_ROOT <- normalizePath(args[[1]], winslash = "/", mustWork = FALSE)
}

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
    hit <- nms[str_detect(nms, regex)]
    if (length(hit) > 0) return(hit[[1]])
  }
  if (required) stop(label, " is missing required column; exact=", paste(exact, collapse = ","), "; regex=", regex)
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

ensure_columns <- function(df, wanted) {
  miss <- setdiff(wanted, names(df))
  for (m in miss) df[[m]] <- NA
  df[, unique(c(wanted, names(df))), drop = FALSE]
}

resolve_sra_table <- function(project_root) {
  candidates <- c(
    file.path(project_root, "configs", project_id, "GSE220189_SraRunTable.csv"),
    file.path(project_root, "configs", project_id, "SraRunTable.csv")
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) {
    stop("No SRA RunTable found. Expected one of:\n", paste(candidates, collapse = "\n"))
  }
  normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
}

resolve_input_root <- function(project_root) {
  candidates <- c(
    Sys.getenv("CELLRANGER_COUNT_ROOT", unset = ""),
    file.path(project_root, "cellranger_count", "output", sra_study),
    file.path(Sys.getenv("ATLAS_CELLRANGER_OUTPUT_ROOT", unset = ""), sra_study),
    file.path(Sys.getenv("ATLAS_CELLRANGER_OUTPUT_ROOT", unset = ""), sra_study)
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hit <- candidates[dir.exists(candidates)]
  if (length(hit) == 0) stop("Cannot find Cell Ranger output root. Tried:\n", paste(candidates, collapse = "\n"))
  normalizePath(hit[[1]], winslash = "/", mustWork = TRUE)
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

input_root <- resolve_input_root(PROJECT_ROOT)
sra_file <- resolve_sra_table(PROJECT_ROOT)

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
message("SRA table   : ", sra_file)
message("output_root : ", output_root)

apply_paper_qc_filter <- TRUE
min_cells_after_qc <- 200L

geo_sample_tbl <- tibble::tribble(
  ~geo_accession, ~geo_title, ~paper_sample, ~geo_group,
  "GSM6793461", "AS08-09890, Control, scRNAseq", "AS08-09890", "Control",
  "GSM6793462", "AS09-13278, Control, scRNAseq", "AS09-13278", "Control",
  "GSM6793463", "AS10-21035, Control, scRNAseq", "AS10-21035", "Control",
  "GSM6793464", "AS11-07049, Control, scRNAseq", "AS11-07049", "Control",
  "GSM6793465", "AS11-07881, Control, scRNAseq", "AS11-07881", "Control",
  "GSM6793466", "AS11-12162, Control, scRNAseq", "AS11-12162", "Control",
  "GSM6793467", "AS11-18755, Control, scRNAseq", "AS11-18755", "Control",
  "GSM6793468", "AS13-13951, Control, scRNAseq", "AS13-13951", "Control",
  "GSM6793469", "AS13-08590, Control, scRNAseq", "AS13-08590", "Control",
  "GSM6793470", "AS14-00902, Control, scRNAseq", "AS14-00902", "Control",
  "GSM6793471", "AS14-03700, Control, scRNAseq", "AS14-03700", "Control",
  "GSM6793472", "AS17-00144, Control, scRNAseq", "AS17-00144", "Control",
  "GSM6793473", "AS17-02129, Control, scRNAseq", "AS17-02129", "Control",
  "GSM6793474", "AS18-00669, Control, scRNAseq", "AS18-00669", "Control",
  "GSM6793475", "BMI0037, Control, scRNAseq", "BMI0037", "Control",
  "GSM6793476", "BMI0040, Control, scRNAseq", "BMI0040", "Control",
  "GSM6793477", "BMI0093, Control, scRNAseq", "BMI0093", "Control",
  "GSM6793478", "BMI0094, Control, scRNAseq", "BMI0094", "Control",
  "GSM6793479", "BMI0095, Control, scRNAseq", "BMI0095", "Control",
  "GSM6793480", "BMI0099, Control, scRNAseq", "BMI0099", "Control",
  "GSM6793481", "BMI0101, Control, scRNAseq", "BMI0101", "Control",
  "GSM6793482", "BMI0102, Control, scRNAseq", "BMI0102", "Control",
  "GSM6793483", "BWJ0023, Control, scRNAseq", "BWJ0023", "Control",
  "GSM6793484", "DU19-01S0003453, MSSA, scRNAseq", "DU19-01S0003453", "MSSA",
  "GSM6793485", "DU19-01S0003462, MSSA, scRNAseq", "DU19-01S0003462", "MSSA",
  "GSM6793486", "DU19-01S0003464, MSSA, scRNAseq", "DU19-01S0003464", "MSSA",
  "GSM6793487", "DU19-01S0003466, MSSA, scRNAseq", "DU19-01S0003466", "MSSA",
  "GSM6793488", "DU19-01S0003482, MSSA, scRNAseq", "DU19-01S0003482", "MSSA",
  "GSM6793489", "DU19-01S0003492, MSSA, scRNAseq", "DU19-01S0003492", "MSSA",
  "GSM6793490", "DU19-01S0003507, MSSA, scRNAseq", "DU19-01S0003507", "MSSA",
  "GSM6793491", "DU19-01S0003509, MSSA, scRNAseq", "DU19-01S0003509", "MSSA",
  "GSM6793492", "DU19-01S0003515, MSSA, scRNAseq", "DU19-01S0003515", "MSSA",
  "GSM6793493", "DU19-01S0003527, MSSA, scRNAseq", "DU19-01S0003527", "MSSA",
  "GSM6793494", "DU19-01S0003542, MRSA, scRNAseq", "DU19-01S0003542", "MRSA",
  "GSM6793495", "DU19-01S0003549, MRSA, scRNAseq", "DU19-01S0003549", "MRSA",
  "GSM6793496", "DU19-01S0003978, MRSA, scRNAseq", "DU19-01S0003978", "MRSA",
  "GSM6793497", "DU19-01S0003987, MRSA, scRNAseq", "DU19-01S0003987", "MRSA",
  "GSM6793498", "DU19-01S0003989, MRSA, scRNAseq", "DU19-01S0003989", "MRSA",
  "GSM6793499", "DU19-01S0003992, MRSA, scRNAseq", "DU19-01S0003992", "MRSA",
  "GSM6793500", "DU19-01S0003994, MRSA, scRNAseq", "DU19-01S0003994", "MRSA",
  "GSM6793501", "DU19-01S0004011, MRSA, scRNAseq", "DU19-01S0004011", "MRSA",
  "GSM6793502", "DU19-01S0004013, MRSA, scRNAseq", "DU19-01S0004013", "MRSA",
  "GSM6793503", "DU19-01S0004015, MRSA, scRNAseq", "DU19-01S0004015", "MRSA",
  "GSM6793504", "DU19-01S0004017, MSSA, scRNAseq", "DU19-01S0004017", "MSSA"
) %>%
  mutate(
    sample_title_short = paste(paper_sample, geo_group, "scRNAseq", sep = ", ")
  )

table_s6_tbl <- tibble::tribble(
  ~table_s6_aliquot_id, ~sex_table_s6, ~age_group, ~condition_table_s6, ~paper_scRNAseq_qc_status, ~paper_scATACseq_qc_status,
  "AS08-09890", "M", "<35", "Control", "PASS", "FAIL",
  "AS09-13278", "M", "<35", "Control", "PASS", "PASS",
  "AS10-21035", "M", "<35", "Control", "PASS", "FAIL",
  "AS11-07049", "M", "<35", "Control", "PASS", "FAIL",
  "AS11-07881", "M", "<35", "Control", "PASS", "FAIL",
  "AS11-12162", "M", "35-65", "Control", "PASS", "FAIL",
  "AS11-18755", "M", "35-65", "Control", "PASS", "PASS",
  "AS13-08590", "M", "35-65", "Control", "FAIL", "PASS",
  "AS13-13951", "M", "<35", "Control", "PASS", "FAIL",
  "AS14-00902", "M", "35-65", "Control", "PASS", "FAIL",
  "AS14-03700", "M", "35-65", "Control", "PASS", "PASS",
  "AS17-00144", "M", "35-65", "Control", "PASS", "PASS",
  "AS17-02129", "M", "35-65", "Control", "PASS", "FAIL",
  "AS18-00669", "M", "<35", "Control", "PASS", "FAIL",
  "BMI0037-M03_2", "F", "<35", "Control", "PASS", "FAIL",
  "BMI0040-M03_2", "M", "<35", "Control", "PASS", "FAIL",
  "BMI0093-M03_2", "F", "35-65", "Control", "PASS", "PASS",
  "BMI0094-M03_2", "M", "<35", "Control", "PASS", "PASS",
  "BMI0095-M03_2", "M", "35-65", "Control", "PASS", "PASS",
  "BMI0099-M03_2", "F", "35-65", "Control", "PASS", "PASS",
  "BMI0101-M03_2", "M", "<35", "Control", "PASS", "PASS",
  "BMI0102-M03_2", "M", "<35", "Control", "PASS", "FAIL",
  "BWJ0023-M03_2", "F", "<35", "Control", "PASS", "PASS",
  "DU19-01S0003453", "F", "<35", "MSSA", "PASS", "PASS",
  "DU19-01S0003462", "F", "<35", "MSSA", "PASS", "PASS",
  "DU19-01S0003464", "F", "35-65", "MSSA", "PASS", "PASS",
  "DU19-01S0003466", "M", "<35", "MSSA", "PASS", "PASS",
  "DU19-01S0003482", "M", "35-65", "MSSA", "PASS", "PASS",
  "DU19-01S0003492", "M", "35-65", "MSSA", "PASS", "PASS",
  "DU19-01S0003507", "F", "35-65", "MSSA", "PASS", "PASS",
  "DU19-01S0003509", "F", "35-65", "MSSA", "PASS", "PASS",
  "DU19-01S0003515", "M", "35-65", "MSSA", "PASS", "PASS",
  "DU19-01S0003527", "F", "35-65", "MSSA", "PASS", "PASS",
  "DU19-01S0003542", "M", "35-65", "MRSA", "PASS", "PASS",
  "DU19-01S0003549", "F", ">65", "MRSA", "PASS", "PASS",
  "DU19-01S0003978", "F", "<35", "MRSA", "PASS", "PASS",
  "DU19-01S0003987", "F", "35-65", "MRSA", "PASS", "PASS",
  "DU19-01S0003989", "F", "<35", "MRSA", "PASS", "PASS",
  "DU19-01S0003992", "M", ">65", "MRSA", "PASS", "PASS",
  "DU19-01S0003994", "M", "35-65", "MRSA", "PASS", "PASS",
  "DU19-01S0004011", "F", "35-65", "MRSA", "PASS", "PASS",
  "DU19-01S0004013", "F", "35-65", "MRSA", "PASS", "PASS",
  "DU19-01S0004015", "M", ">65", "MRSA", "PASS", "PASS",
  "DU19-01S0004017", "M", ">65", "MSSA", "PASS", "PASS"
) %>%
  mutate(
    table_s6_sample_key = stringr::str_remove(table_s6_aliquot_id, "-M03_2$"),
    sex = dplyr::case_when(
      sex_table_s6 == "M" ~ "male",
      sex_table_s6 == "F" ~ "female",
      TRUE ~ NA_character_
    ),
    age_range = age_group,
    age_numeric_resolution = "range_only_not_exact_age",
    paper_scRNAseq_pass = ifelse(paper_scRNAseq_qc_status == "PASS", "yes", "no"),
    paper_scATACseq_pass = ifelse(paper_scATACseq_qc_status == "PASS", "yes", "no")
  )


paper_methods_tbl <- tibble::tribble(
  ~field, ~value, ~note,
  "paper_title", "Mapping disease regulatory circuits at cell-type resolution from single-cell multiomics data", "Nature Computational Science 2023",
  "series_accession", "GSE220189", "S. aureus scRNA-seq SubSeries",
  "parent_series", "GSE220190", "SuperSeries containing S. aureus scATAC-seq, S. aureus scRNA-seq, and COVID-19 scATAC-seq",
  "cohort_mrsa_n", "10", "S. aureus bloodstream infection, MRSA",
  "cohort_mssa_n", "11", "S. aureus bloodstream infection, MSSA",
  "cohort_control_n", "23", "uninfected controls",
  "high_quality_scRNA_cells_reported", "276200", "paper Results",
  "source_material", "PBMC", "paper Methods and GEO",
  "sample_state", "frozen PBMC vials thawed before assay", "paper Methods",
  "scRNA_protocol", "10x Genomics Single Cell 3prime Reagents Kits V3.1; Chromium Single Cell 3prime Chip G", "paper Methods/GEO",
  "target_cells_loaded", "5000-10000 final cells", "paper Methods",
  "paper_original_aligner", "10x Genomics Cell Ranger", "paper Methods",
  "paper_original_cellranger_version", "1.2", "paper Methods",
  "paper_original_reference_genome", "hg38", "paper Methods",
  "paper_original_seurat_version", "Seurat V4 for DEG; Seurat-based integration", "paper Methods",
  "paper_qc_min_genes_per_cell", "400", "paper scRNA QC",
  "paper_qc_max_genes_per_cell", "5000", "paper scRNA QC",
  "paper_qc_max_percent_mt", "10", "paper scRNA QC",
  "paper_normalization", "SCTransform; mitochondrial reads and cell cycle heterogeneity regressed out", "paper Methods",
  "table_s6_note", "Table S6 provides aliquot ID, gender, age range, condition, scRNAseq status, scATACseq status", "embedded"
)

write_csv(geo_sample_tbl, file.path(manifest_dir, "GSE220189_geo_sample_embedded.csv"))
write_csv(table_s6_tbl, file.path(manifest_dir, "GSE220189_supplementary_table_s6_embedded.csv"))
write_csv(paper_methods_tbl, file.path(manifest_dir, "GSE220189_paper_methods_and_results_embedded.csv"))

sra_raw <- read_csv(sra_file, show_col_types = FALSE, name_repair = "unique")
names(sra_raw) <- normalize_names(names(sra_raw))

run_col <- pick_column(sra_raw, exact = c("run"), regex = "^run$", required = FALSE)
assay_col <- pick_column(sra_raw, exact = c("assay_type"), regex = "assay.*type", required = FALSE)
avgspot_col <- pick_column(sra_raw, exact = c("avgspotlen"), regex = "avgspotlen", required = FALSE)
bases_col <- pick_column(sra_raw, exact = c("bases"), regex = "^bases$", required = FALSE)
bytes_col <- pick_column(sra_raw, exact = c("bytes"), regex = "^bytes$", required = FALSE)
biosample_col <- pick_column(sra_raw, exact = c("biosample"), regex = "^biosample$", required = TRUE, label = "SRA table")
bioproject_col <- pick_column(sra_raw, exact = c("bioproject"), regex = "^bioproject$", required = TRUE, label = "SRA table")
exp_col <- pick_column(sra_raw, exact = c("experiment"), regex = "^experiment$", required = FALSE)
geo_col <- pick_column(sra_raw, exact = c("geo_accession_exp", "geo_accession"), regex = "geo.*accession", required = FALSE)
sample_name_col <- pick_column(sra_raw, exact = c("sample_name"), regex = "^sample.*name$", required = FALSE)
library_name_col <- pick_column(sra_raw, exact = c("library_name"), regex = "^library.*name$", required = FALSE)
disease_state_col <- pick_column(sra_raw, exact = c("disease_state"), regex = "disease.*state", required = FALSE)
sra_study_col <- pick_column(sra_raw, exact = c("sra_study"), regex = "sra.*study", required = TRUE, label = "SRA table")
cell_type_col <- pick_column(sra_raw, exact = c("cell_type"), regex = "cell.*type", required = FALSE)
source_name_col <- pick_column(sra_raw, exact = c("source_name"), regex = "source.*name", required = FALSE)
tissue_col <- pick_column(sra_raw, exact = c("tissue"), regex = "^tissue$", required = FALSE)
instrument_col <- pick_column(sra_raw, exact = c("instrument"), regex = "instrument", required = FALSE)
library_layout_col <- pick_column(sra_raw, exact = c("librarylayout", "library_layout"), regex = "library.*layout", required = FALSE)
library_selection_col <- pick_column(sra_raw, exact = c("libraryselection", "library_selection"), regex = "library.*selection", required = FALSE)
library_source_col <- pick_column(sra_raw, exact = c("librarysource", "library_source"), regex = "library.*source", required = FALSE)
organism_col <- pick_column(sra_raw, exact = c("organism"), regex = "^organism$", required = FALSE)
platform_col <- pick_column(sra_raw, exact = c("platform"), regex = "^platform$", required = FALSE)
center_col <- pick_column(sra_raw, exact = c("center_name"), regex = "center.*name", required = FALSE)
consent_col <- pick_column(sra_raw, exact = c("consent"), regex = "^consent$", required = FALSE)
release_col <- pick_column(sra_raw, exact = c("releasedate"), regex = "releasedate", required = FALSE)
create_col <- pick_column(sra_raw, exact = c("create_date"), regex = "create_date", required = FALSE)
version_col <- pick_column(sra_raw, exact = c("version"), regex = "^version$", required = FALSE)
filetype_col <- pick_column(sra_raw, exact = c("datastore_filetype"), regex = "datastore.*filetype", required = FALSE)
provider_col <- pick_column(sra_raw, exact = c("datastore_provider"), regex = "datastore.*provider", required = FALSE)
region_col <- pick_column(sra_raw, exact = c("datastore_region"), regex = "datastore.*region", required = FALSE)

sra_norm <- tibble(
  run_accession = col_or_na(sra_raw, run_col),
  experiment_accession = col_or_na(sra_raw, exp_col),
  biosample_accession = col_or_na(sra_raw, biosample_col),
  geo_accession_raw = col_or_na(sra_raw, geo_col),
  sample_name_sra = col_or_na(sra_raw, sample_name_col),
  library_name_sra = col_or_na(sra_raw, library_name_col),
  disease_state_sra = col_or_na(sra_raw, disease_state_col),
  bioproject_sra = col_or_na(sra_raw, bioproject_col),
  sra_study_sra = col_or_na(sra_raw, sra_study_col),
  assay_type_sra = col_or_na(sra_raw, assay_col),
  avg_spot_len_sra = col_or_na(sra_raw, avgspot_col, "numeric"),
  bases_sra = col_or_na(sra_raw, bases_col, "numeric"),
  bytes_sra = col_or_na(sra_raw, bytes_col, "numeric"),
  cell_type_sra = col_or_na(sra_raw, cell_type_col),
  source_name_sra = col_or_na(sra_raw, source_name_col),
  tissue_sra = col_or_na(sra_raw, tissue_col),
  instrument_sra = col_or_na(sra_raw, instrument_col),
  library_layout_sra = col_or_na(sra_raw, library_layout_col),
  library_selection_sra = col_or_na(sra_raw, library_selection_col),
  library_source_sra = col_or_na(sra_raw, library_source_col),
  organism_sra = col_or_na(sra_raw, organism_col),
  platform_sra = col_or_na(sra_raw, platform_col),
  center_name_sra = col_or_na(sra_raw, center_col),
  consent_sra = col_or_na(sra_raw, consent_col),
  release_date_sra = col_or_na(sra_raw, release_col),
  create_date_sra = col_or_na(sra_raw, create_col),
  version_sra = col_or_na(sra_raw, version_col),
  datastore_filetype_sra = col_or_na(sra_raw, filetype_col),
  datastore_provider_sra = col_or_na(sra_raw, provider_col),
  datastore_region_sra = col_or_na(sra_raw, region_col)
) %>% mutate(across(where(is.character), ~na_if(.x, "")))

write_csv(sra_norm, file.path(manifest_dir, "GSE220189_sra_table_normalized.csv"))

if (!any(sra_norm$bioproject_sra == bioproject, na.rm = TRUE)) stop("SRA RunTable does not contain expected BioProject ", bioproject)
if (!any(sra_norm$sra_study_sra == sra_study, na.rm = TRUE)) stop("SRA RunTable does not contain expected SRA Study ", sra_study)

sra_sample_agg0 <- sra_norm %>%
  group_by(biosample_accession) %>%
  summarise(
    run_accession = paste(unique(na.omit(run_accession)), collapse = ";"),
    experiment_accession = paste(unique(na.omit(experiment_accession)), collapse = ";"),
    n_sra_runs = n_distinct(run_accession),
    geo_accession_raw = paste(unique(na.omit(geo_accession_raw)), collapse = ";"),
    sample_name_sra = paste(unique(na.omit(sample_name_sra)), collapse = ";"),
    library_name_sra = paste(unique(na.omit(library_name_sra)), collapse = ";"),
    disease_state_sra = paste(unique(na.omit(disease_state_sra)), collapse = ";"),
    bioproject_sra = paste(unique(na.omit(bioproject_sra)), collapse = ";"),
    sra_study_sra = paste(unique(na.omit(sra_study_sra)), collapse = ";"),
    assay_type_sra = paste(unique(na.omit(assay_type_sra)), collapse = ";"),
    avg_spot_len_sra = suppressWarnings(mean(avg_spot_len_sra, na.rm = TRUE)),
    bases_sra = suppressWarnings(sum(bases_sra, na.rm = TRUE)),
    bytes_sra = suppressWarnings(sum(bytes_sra, na.rm = TRUE)),
    cell_type_sra = paste(unique(na.omit(cell_type_sra)), collapse = ";"),
    source_name_sra = paste(unique(na.omit(source_name_sra)), collapse = ";"),
    tissue_sra = paste(unique(na.omit(tissue_sra)), collapse = ";"),
    instrument_sra = paste(unique(na.omit(instrument_sra)), collapse = ";"),
    library_layout_sra = paste(unique(na.omit(library_layout_sra)), collapse = ";"),
    library_selection_sra = paste(unique(na.omit(library_selection_sra)), collapse = ";"),
    library_source_sra = paste(unique(na.omit(library_source_sra)), collapse = ";"),
    organism_sra = paste(unique(na.omit(organism_sra)), collapse = ";"),
    platform_sra = paste(unique(na.omit(platform_sra)), collapse = ";"),
    center_name_sra = paste(unique(na.omit(center_name_sra)), collapse = ";"),
    consent_sra = paste(unique(na.omit(consent_sra)), collapse = ";"),
    release_date_sra = paste(unique(na.omit(release_date_sra)), collapse = ";"),
    create_date_sra = paste(unique(na.omit(create_date_sra)), collapse = ";"),
    version_sra = paste(unique(na.omit(version_sra)), collapse = ";"),
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

if (nrow(sra_sample_agg0) != expected_n) stop("Expected 44 unique BioSamples, got ", nrow(sra_sample_agg0))

key_tbl <- bind_rows(
  geo_sample_tbl %>% transmute(geo_accession, key = geo_accession, key_type = "geo_accession"),
  geo_sample_tbl %>% transmute(geo_accession, key = paper_sample, key_type = "paper_sample"),
  geo_sample_tbl %>% transmute(geo_accession, key = geo_title, key_type = "geo_title"),
  geo_sample_tbl %>% transmute(geo_accession, key = sample_title_short, key_type = "sample_title_short")
) %>%
  mutate(key = na_if(key, "")) %>%
  filter(!is.na(key))

sra_key_tbl <- sra_sample_agg0 %>%
  mutate(row_id = row_number()) %>%
  select(row_id, biosample_accession, geo_accession_raw, sample_name_sra, library_name_sra) %>%
  pivot_longer(cols = c(geo_accession_raw, sample_name_sra, library_name_sra), names_to = "sra_key_source", values_to = "key") %>%
  mutate(key = na_if(key, "")) %>%
  filter(!is.na(key)) %>%
  separate_rows(key, sep = ";") %>%
  mutate(key = str_trim(key)) %>%
  filter(key != "")

match_tbl_raw <- sra_key_tbl %>%
  inner_join(key_tbl, by = "key", relationship = "many-to-many") %>%
  distinct(row_id, biosample_accession, geo_accession, key, key_type, sra_key_source)

match_summary <- match_tbl_raw %>%
  group_by(row_id, biosample_accession) %>%
  summarise(
    n_geo = n_distinct(geo_accession),
    geo_accession = paste(unique(geo_accession), collapse = ";"),
    matched_key = paste(unique(key), collapse = ";"),
    matched_key_type = paste(unique(key_type), collapse = ";"),
    matched_sra_key_source = paste(unique(sra_key_source), collapse = ";"),
    .groups = "drop"
  )

ambig_biosample <- match_summary %>% filter(n_geo > 1)
ambig_geo <- match_summary %>% separate_rows(geo_accession, sep = ";") %>% count(geo_accession) %>% filter(!is.na(geo_accession), n > 1)

if (nrow(ambig_biosample) > 0 || nrow(ambig_geo) > 0) {
  write_csv(match_tbl_raw, file.path(manifest_dir, "GSE220189_explicit_key_matches_ambiguous_raw.csv"))
  write_csv(match_summary, file.path(manifest_dir, "GSE220189_explicit_key_matches_ambiguous_summary.csv"))
  stop("Explicit SRA-to-GEO matching is biologically ambiguous. Inspect manifest files.")
}

write_csv(match_tbl_raw, file.path(manifest_dir, "GSE220189_explicit_key_matches_raw.csv"))
write_csv(match_summary, file.path(manifest_dir, "GSE220189_explicit_key_matches_summary.csv"))

sra_sample_agg <- sra_sample_agg0 %>%
  mutate(row_id = row_number()) %>%
  left_join(match_summary %>% select(row_id, geo_accession, matched_key, matched_key_type, matched_sra_key_source), by = "row_id") %>%
  select(-row_id)

write_csv(sra_sample_agg, file.path(manifest_dir, "GSE220189_sra_sample_aggregated.csv"))

missing_geo <- setdiff(geo_sample_tbl$geo_accession, sra_sample_agg$geo_accession)
if (length(missing_geo) > 0 || any(is.na(sra_sample_agg$geo_accession))) {
  stop("Could not fully establish explicit GSM-to-BioSample mapping. Missing GEO: ", paste(missing_geo, collapse = ", "))
}

sample_meta <- geo_sample_tbl %>%
  left_join(sra_sample_agg, by = "geo_accession") %>%
  left_join(table_s6_tbl, by = c("paper_sample" = "table_s6_sample_key")) %>%
  mutate(
    database_accession = project_id,
    parent_series = parent_series,
    sra_study = sra_study,
    bioproject = bioproject,
    organism = organism,
    cellranger_dir_id = biosample_accession,
    sample_id = geo_accession,
    sample_alias = paper_sample,
    patient_id = paper_sample,
    subject_id = paper_sample,
    participant_id = paper_sample,
    disease_state_sra_primary = coalesce(disease_state_sra, geo_group),
    disease_state_norm = case_when(
      str_to_lower(disease_state_sra_primary) == "control" ~ "healthy_control",
      str_to_lower(disease_state_sra_primary) == "mssa" ~ "MSSA_sepsis",
      str_to_lower(disease_state_sra_primary) == "mrsa" ~ "MRSA_sepsis",
      TRUE ~ NA_character_
    ),
    study_group = disease_state_norm,
    factor = ifelse(study_group == "healthy_control", "healthy_control", "sepsis"),
    cov19 = "not_applicable",
    is_sepsis = ifelse(factor == "sepsis", "yes", "no"),
    ards = NA_character_,
    pneumonia = NA_character_,
    infection_site = ifelse(factor == "sepsis", "bloodstream", "none"),
    pathogen_etiology = ifelse(factor == "sepsis", "Staphylococcus aureus", "none"),
    antibiotic_resistance = case_when(
      study_group == "MRSA_sepsis" ~ "MRSA",
      study_group == "MSSA_sepsis" ~ "MSSA",
      TRUE ~ NA_character_
    ),
    pathogen_species = ifelse(factor == "sepsis", "Staphylococcus aureus", "none"),
    culture_confirmed_bloodstream_infection = ifelse(factor == "sepsis", "yes", "no"),
    mortality = NA_character_,
    final_outcome = NA_character_,
    icu_admission = NA_character_,
    days = NA_real_,
    timepoint_label = "bloodstream_infection_or_control_baseline",
    severity = study_group,
    severityatday = paste0(severity, "_", timepoint_label),
    ventilation = NA_character_,
    age = NA_real_,
    race = NA_character_,
    qsofa = NA_real_,
    apache_ii = NA_real_,
    sofa = NA_real_,
    diagnosis = ifelse(factor == "sepsis", "culture-confirmed S. aureus bloodstream infection", "uninfected healthy control"),
    source_name = coalesce(source_name_sra, "PBMC"),
    tissue_label = coalesce(tissue_sra, "peripheral blood"),
    source_blood_fraction = "PBMC",
    cell_type = coalesce(cell_type_sra, "peripheral blood mononuclear cells"),
    cell_line = "Not a cell line",
    extracted_molecule = "total RNA",
    frozen_or_fresh = "frozen PBMC thawed before 10x scRNA-seq",
    fresh_or_frozen = frozen_or_fresh,
    pbmc_isolation_method = "PBMC vials thawed; RPMI/FBS dilution and wash; Countess viability count",
    ficoll = NA_character_,
    rbc_lysis = NA_character_,
    sample_storage = "frozen PBMC vial",
    cell_suspension_buffer = "RPMI/FBS during thawing; 10x loading after filtering/counting",
    viability_reported = "assessed by Trypan Blue on Countess II",
    sorted_population = "PBMC",
    platform = "10x Genomics",
    sequencer_model = coalesce(instrument_sra, "NextSeq 2000"),
    library_type = "10x Genomics single-cell gene expression",
    chemistry_version = "10x Genomics Single Cell 3prime Reagents Kits V3.1",
    feature_barcoding = "none for RNA object",
    vdj_capture = "not loaded here",
    read_len = NA_character_,
    assay_type = coalesce(assay_type_sra, "RNA-Seq"),
    library_layout = coalesce(library_layout_sra, "PAIRED"),
    library_selection = coalesce(library_selection_sra, "cDNA"),
    library_source = coalesce(library_source_sra, "TRANSCRIPTOMIC SINGLE CELL"),
    aligner = "Cell Ranger / user common atlas pipeline for this object; paper original used Cell Ranger v1.2",
    cellranger_or_equivalent_version = "atlas common reprocessing; paper provenance: Cell Ranger v1.2",
    paper_original_aligner = "10x Genomics Cell Ranger",
    paper_original_cellranger_version = "1.2",
    paper_original_reference_genome = "hg38",
    reference_genome = "hg38",
    gene_annotation = "hg38 in paper provenance",
    seq_depth = NA_character_,
    n_loaded_cells = NA_integer_,
    paper_cells_after_qc_total = 276200L,
    paper_qc_min_genes_per_cell = 400L,
    paper_qc_max_genes_per_cell = 5000L,
    paper_qc_max_percent_mt = 10,
    paper_sctransform_regressors = "mitochondrial reads and cell cycle heterogeneity",
    clinical_mapping_status = "group_sex_age_range_from_GEO_SRA_TableS6_no_exact_numeric_age_or_outcome",
    enrollment_window = "S. aureus bloodstream infection subjects and uninfected controls; exact sampling time not encoded in RunTable",
    batch_unit_recommended = "sample",
    source_note = paste(
      "GSE220189 is the S. aureus scRNA-seq SubSeries of GSE220190.",
      "GSM-to-BioSample mapping is established from explicit SRA Sample Name/Library Name/GSM keys.",
      "Disease group is assigned primarily from SRA disease_state and cross-checked against GEO title and Supplementary Table S6.",
      "Sex and age-range metadata are embedded from Nat Comput Sci Supplementary Table 6.",
      "Exact numeric age, mortality, and outcome are not embedded because Table S6 only reports age bins and no outcome fields."
    ),
    read_mode = NA_character_,
    predicted.cluster = NA_character_
  )

if (nrow(sample_meta) != expected_n) stop("Expected 44 sample metadata rows, got ", nrow(sample_meta))
if (any(is.na(sample_meta$geo_accession) | is.na(sample_meta$biosample_accession) | is.na(sample_meta$factor) | is.na(sample_meta$study_group))) {
  bad <- sample_meta %>% filter(is.na(geo_accession) | is.na(biosample_accession) | is.na(factor) | is.na(study_group))
  print(bad)
  stop("Missing key sample metadata.")
}
if (any(is.na(sample_meta$sex) | is.na(sample_meta$age_group) | is.na(sample_meta$table_s6_aliquot_id))) {
  bad <- sample_meta %>% filter(is.na(sex) | is.na(age_group) | is.na(table_s6_aliquot_id))
  print(bad %>% select(geo_accession, paper_sample, geo_group, table_s6_aliquot_id, sex, age_group))
  stop("Some samples did not receive Supplementary Table S6 sex/age-range metadata.")
}
table_s6_conflict <- sample_meta %>% filter(condition_table_s6 != geo_group)
if (nrow(table_s6_conflict) > 0) {
  write_csv(table_s6_conflict, file.path(manifest_dir, "GSE220189_table_s6_geo_condition_conflicts.csv"))
  stop("Supplementary Table S6 condition conflicts with GEO group.")
}
conflict <- sample_meta %>%
  mutate(
    geo_expected_norm = case_when(
      geo_group == "Control" ~ "healthy_control",
      geo_group == "MSSA" ~ "MSSA_sepsis",
      geo_group == "MRSA" ~ "MRSA_sepsis",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(study_group != geo_expected_norm)
if (nrow(conflict) > 0) {
  write_csv(conflict, file.path(manifest_dir, "GSE220189_sra_geo_group_conflicts.csv"))
  stop("SRA disease_state conflicts with embedded GEO group.")
}

common_meta_columns <- c(
  "database_accession", "parent_series", "sra_study", "bioproject", "organism",
  "biosample_accession", "experiment_accession", "run_accession", "geo_accession",
  "cellranger_dir_id", "sample_id", "sample_alias", "paper_sample", "geo_title",
  "patient_id", "subject_id", "participant_id", "clinical_mapping_status",
  "disease_state_sra", "disease_state_sra_primary", "disease_state_norm",
  "geo_group", "study_group", "factor", "cov19", "is_sepsis", "ards", "pneumonia",
  "infection_site", "pathogen_etiology", "pathogen_species", "antibiotic_resistance",
  "culture_confirmed_bloodstream_infection", "mortality", "final_outcome", "icu_admission",
  "days", "timepoint_label", "severity", "severityatday", "ventilation",
  "sex", "age", "age_group", "age_range", "age_numeric_resolution", "race", "qsofa", "apache_ii", "sofa", "diagnosis",
  "source_name", "tissue_label", "source_blood_fraction", "cell_type", "cell_line",
  "extracted_molecule", "frozen_or_fresh", "fresh_or_frozen", "pbmc_isolation_method",
  "ficoll", "rbc_lysis", "sample_storage", "cell_suspension_buffer", "viability_reported",
  "sorted_population", "platform", "sequencer_model", "library_type", "chemistry_version",
  "feature_barcoding", "vdj_capture", "read_len", "assay_type",
  "library_layout", "library_selection", "library_source",
  "aligner", "cellranger_or_equivalent_version", "paper_original_aligner",
  "paper_original_cellranger_version", "paper_original_reference_genome",
  "reference_genome", "gene_annotation", "seq_depth", "n_loaded_cells",
  "paper_cells_after_qc_total", "paper_qc_min_genes_per_cell",
  "paper_qc_max_genes_per_cell", "paper_qc_max_percent_mt", "paper_sctransform_regressors",
  "n_sra_runs", "sample_name_sra", "library_name_sra",
  "table_s6_aliquot_id", "sex_table_s6", "condition_table_s6",
  "paper_scRNAseq_qc_status", "paper_scATACseq_qc_status", "paper_scRNAseq_pass", "paper_scATACseq_pass",
  "bioproject_sra", "sra_study_sra", "avg_spot_len_sra", "bases_sra", "bytes_sra",
  "instrument_sra", "platform_sra", "consent_sra", "center_name_sra",
  "release_date_sra", "create_date_sra", "version_sra",
  "source_name_sra", "tissue_sra", "cell_type_sra",
  "datastore_filetype_sra", "datastore_provider_sra", "datastore_region_sra",
  "matched_key", "matched_key_type", "matched_sra_key_source",
  "enrollment_window", "batch_unit_recommended", "source_note",
  "read_mode", "predicted.cluster"
)
sample_meta <- ensure_columns(sample_meta, common_meta_columns)

sample_meta <- atlas_standardize_metadata(sample_meta, project_id = project_id)
write_csv(sample_meta %>% count(factor, study_group, geo_group, condition_table_s6, disease_state_sra, antibiotic_resistance, sex, age_group, paper_scRNAseq_qc_status, paper_scATACseq_qc_status),
          file.path(manifest_dir, "GSE220189_group_tableS6_distribution_crosscheck.csv"))

input_dirs_all <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
input_dirs <- input_dirs_all[basename(input_dirs_all) %in% sample_meta$cellranger_dir_id]
input_dirs <- input_dirs[order(match(basename(input_dirs), sample_meta$cellranger_dir_id))]
if (length(input_dirs) != expected_n) {
  stop("Matched ", length(input_dirs), " Cell Ranger directories but metadata has 44 rows. Missing: ",
       paste(setdiff(sample_meta$cellranger_dir_id, basename(input_dirs)), collapse = ", "))
}

sample_meta <- atlas_standardize_metadata(sample_meta, project_id = project_id)

obj_list <- list()
qc_summary <- tibble(
  cellranger_dir_id = character(),
  run_accession = character(),
  biosample_accession = character(),
  geo_accession = character(),
  paper_sample = character(),
  study_group = character(),
  factor = character(),
  sex = character(),
  age_group = character(),
  paper_scRNAseq_qc_status = character(),
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
  if (!dir.exists(raw_dir)) stop("Missing raw_feature_bc_matrix: ", raw_dir)

  idx <- match(dir_id, sample_meta$cellranger_dir_id)
  if (is.na(idx)) stop("No sample metadata for ", dir_id)
  meta_row <- sample_meta[idx, , drop = FALSE]
  if (nrow(meta_row) != 1) stop("Sample metadata row is not unique for ", dir_id)

  message("Reading ", dir_id, " -> ", meta_row$geo_accession, " (", meta_row$paper_sample, ", ", meta_row$study_group, ")")
  # ATLAS_RESOURCE_SAFE_EXCLUSION_GSE220189
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
        study_group = meta_row$study_group,
        factor = meta_row$factor,
        sex = meta_row$sex,
        age_group = meta_row$age_group,
        paper_scRNAseq_qc_status = meta_row$paper_scRNAseq_qc_status,
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
  sample_meta$read_mode[idx] <- res$read_mode
  sample_meta$n_loaded_cells[idx] <- n_raw
  meta_row$read_mode <- res$read_mode
  meta_row$n_loaded_cells <- n_raw

  obj <- attach_sample_metadata(obj, meta_row)

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.ribo"]] <- compute_percent_ribo(obj)
  obj[["percent.hb"]] <- compute_percent_hb(obj)

  obj[["paper_qc_pass"]] <- with(obj@meta.data,
    nFeature_RNA > paper_qc_min_genes_per_cell &
      nFeature_RNA < paper_qc_max_genes_per_cell &
      percent.mt < paper_qc_max_percent_mt
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
      cellranger_dir_id = dir_id,
      run_accession = meta_row$run_accession,
      biosample_accession = meta_row$biosample_accession,
      geo_accession = meta_row$geo_accession,
      paper_sample = meta_row$paper_sample,
      study_group = meta_row$study_group,
      factor = meta_row$factor,
      sex = meta_row$sex,
      age_group = meta_row$age_group,
      paper_scRNAseq_qc_status = meta_row$paper_scRNAseq_qc_status,
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
  saveRDS(obj_qc, qc_rds, compress = TRUE)

  obj_work <- NormalizeData(obj_qc, verbose = FALSE)
  obj_work <- FindVariableFeatures(obj_work, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  work_rds <- file.path(rds_work_dir, paste0(meta_row$geo_accession, "__", dir_id, "__working.rds"))
  saveRDS(obj_work, work_rds, compress = TRUE)

  obj_list[[as.character(meta_row$geo_accession)]] <- work_rds

  qc_summary <- bind_rows(qc_summary, tibble(
    cellranger_dir_id = dir_id,
    run_accession = meta_row$run_accession,
    biosample_accession = meta_row$biosample_accession,
    geo_accession = meta_row$geo_accession,
    paper_sample = meta_row$paper_sample,
    study_group = meta_row$study_group,
    factor = meta_row$factor,
    sex = meta_row$sex,
    age_group = meta_row$age_group,
    paper_scRNAseq_qc_status = meta_row$paper_scRNAseq_qc_status,
    n_raw = n_raw,
    n_after_qc = n_after_qc,
    saved = TRUE,
    preqc_rds_path = preqc_rds,
    qc_rds_path = qc_rds,
    working_rds_path = work_rds
  ))
  # ATLAS_RESOURCE_SAFE_GC_GSE220189
  rm(res, obj, obj_qc, obj_work)
  invisible(gc())
}

write_csv(sample_meta, file.path(meta_dir, "GSE220189_sample_metadata.csv"))
write_csv(sample_meta, file.path(meta_dir, "GSE220189_run_metadata.csv"))
write_csv(qc_summary, file.path(qc_dir, "GSE220189_qc_summary.csv"))
write_csv(
  qc_summary %>% dplyr::filter(preprocess_status %in% "excluded_before_rds"),
  file.path(qc_dir, "GSE220189_excluded_before_rds.csv")
)

registry <- qc_summary %>%
  filter(saved) %>%
  left_join(sample_meta %>%
              select(cellranger_dir_id, run_accession, biosample_accession, geo_accession,
                     paper_sample, study_group, factor, is_sepsis, disease_state_sra,
                     antibiotic_resistance, pathogen_etiology, clinical_mapping_status,
                     sex, age_group, paper_scRNAseq_qc_status, batch_unit_recommended),
            by = c("cellranger_dir_id", "run_accession", "biosample_accession", "geo_accession", "paper_sample", "study_group", "factor", "sex", "age_group", "paper_scRNAseq_qc_status")) %>%
  mutate(project_id = project_id, preintegration_tier = "working", object_grain = "sample")

write_csv(registry, file.path(output_root, "GSE220189_preintegration_registry.csv"))
saveRDS(obj_list, file.path(output_root, "GSE220189_preintegration_objlist.rds"))

message("Done. Saved pre-QC raw objects to    : ", rds_preqc_dir)
message("Done. Saved QC-filtered raw objects to: ", rds_qc_dir)
message("Done. Saved working objects to       : ", rds_work_dir)
message("Registry: ", file.path(output_root, "GSE220189_preintegration_registry.csv"))
