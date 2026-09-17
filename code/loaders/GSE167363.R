#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(purrr)
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
# GSE167363 pre-integration loader (run-level Cell Ranger outputs)
#
# Design principles
# - Atlas-side raw data were reprocessed by a common pipeline. Original paper
#   processing settings are stored as provenance metadata only.
# - Output structure is aligned with prior loaders:
#     rds_preqc_raw/       metadata + QC metrics, before filtering
#     rds_qcfiltered_raw/  post-filter raw counts
#     rds_working/         normalized/HVG working object
# - External input file required:
#     configs/GSE167363/GSE167363_SraRunTable.csv
#   (fallback also allows configs/GSE167363/SraRunTable.csv)
# - Metadata are maximally embedded from:
#     GEO sample labels
#     SRA Run Table
#     manuscript main text / Table 1 / methods
#     Supplementary Table S1
#
# Important study-specific note:
#   GEO/SRA sample labels:
#     GSM5102902/03 = NS LS T0/T6 (male, age 90-95, donor 3)
#     GSM5511351/52 = NS ES_T0/T6 (female, age 65-70, donor 5)
#   Manuscript text states:
#     female nonsurvivor P50 = NS LS
#     male nonsurvivor P34   = NS ES
#   Therefore this loader stores BOTH:
#     - stage_label_geo_raw     (exact GEO/SRA sample-side label)
#     - patient_id_inferred     (P34/P50 inferred from SRA age/sex + paper text)
#     - stage_label_manuscript_inferred
#     - label_conflict_flag
#   so downstream analyses can choose which layer to use explicitly.
# ==============================================================================

project_id <- "GSE167363"
sra_study  <- "SRP307774"
bioproject <- "PRJNA704394"
organism   <- "Homo sapiens"

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

normalize_names <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "") %>%
    tolower()
}

pick_column <- function(df, exact = NULL, regex = NULL, required = TRUE, label = "table") {
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

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

midpoint_from_range <- function(x) {
  x <- as.character(x)
  out <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    xi <- x[[i]]
    if (is.na(xi) || xi == "") next
    if (stringr::str_detect(xi, "^[0-9.]+\\s*-\\s*[0-9.]+$")) {
      vals <- as.numeric(unlist(strsplit(gsub("\\s+", "", xi), "-")))
      out[[i]] <- mean(vals)
    } else {
      out[[i]] <- suppressWarnings(as.numeric(xi))
    }
  }
  out
}

canonical_label <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_trim() %>%
    stringr::str_replace_all("\\s+", "_") %>%
    stringr::str_replace_all("__+", "_")
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
  # Hook reserved for future common helper insertions.
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
# Atlas-side QC / output policy
# ------------------------------------------------------------------------------
apply_paper_qc_filter   <- TRUE
apply_atlas_mt10_filter <- FALSE
min_cells_after_qc      <- 200L

# ------------------------------------------------------------------------------
# Embedded metadata from paper + supplement + GEO
# ------------------------------------------------------------------------------
# Table 1 (paper):
# controls: male 35-40, female 45-50
# nonsurvivors: male 90-95, female 65-70
# survivors: male 45-50, female 65-70, female 70-75
#
# Paper text:
# female nonsurvivor P50 = NS LS
# male nonsurvivor   P34 = NS ES
#
# SRA sample labels:
# male 90-95 -> NS LS T0/T6
# female 65-70 -> NS ES_T0/T6
#
# Therefore patient_id_inferred is defensible using SRA age/sex + paper text,
# but stage labels differ between GEO/SRA and manuscript interpretation. We store
# both and explicitly flag the conflict.

paper_subject_tbl <- tibble::tribble(
  ~subject_id_standard, ~patient_id_inferred, ~subject_type,         ~sex_paper, ~age_range_paper, ~age_midpoint_paper,
  ~outcome_label,      ~manuscript_stage_inferred, ~manuscript_identity_note,
  ~sepsis_etiology,    ~apache_ii, ~sofa, ~qsofa, ~time_of_death_days_post_enrollment_reported,
  ~plasma_resistin_reported, ~plasma_il6_reported, ~plasma_il8_reported, ~plasma_il10_reported, ~lps_induced_tnf_alpha_reported,
  "HC1",               "HC1",      "healthy_control", "male",   "35-40", 37.5,
  "control",           NA_character_,
  "Healthy control donor from Table 1 and GEO label HC1.",
  NA_character_,       NA_real_,   NA_real_, NA_real_, NA_character_,
  "22.5",              "N.D",      "0.03",             "N.D.",              "0.656",
  "HC2",               "HC2",      "healthy_control", "female", "45-50", 47.5,
  "control",           NA_character_,
  "Healthy control donor from Table 1 and GEO label HC2.",
  NA_character_,       NA_real_,   NA_real_, NA_real_, NA_character_,
  "36.7",              "0.002",    "0.026",            "N.D.",              "0.979",
  "NS_male_90_95",     "P34",      "sepsis_patient",   "male",   "90-95", 92.5,
  "non_survivor",      "NS_ES",
  "P34 inferred from manuscript text: male nonsurvivor designated NS ES.",
  "E.coli bacteremia", 18,         11,       NA_real_, "<30",
  "202",               "30.2",     "6.65",             "0.13",              "0.45",
  "NS_female_65_70",   "P50",      "sepsis_patient",   "female", "65-70", 67.5,
  "non_survivor",      "NS_LS",
  "P50 inferred from manuscript text: female nonsurvivor designated NS LS.",
  "E.coli bacteremia", 38,         16,       NA_real_, "1",
  "135.9",             "142.3",    "27.2",             "9.71",              "0.0046",
  "S1",                "S1",       "sepsis_patient",   "female", "65-70", 67.5,
  "survivor",          "survivor",
  "Survivor subject S1 from GEO labels and Table 1 sex/age combination.",
  "E.coli bacteremia", 41,         15,       NA_real_, NA_character_,
  "281",               "2.48",     "0.61",             "0.15",              "0.047",
  "S2",                "S2",       "sepsis_patient",   "male",   "45-50", 47.5,
  "survivor",          "survivor",
  "Survivor subject S2 from GEO labels and Table 1 sex/age combination.",
  "E.coli bacteremia", 31,         11,       NA_real_, NA_character_,
  "147",               "133",      "41.7",             "0.52",              NA_character_,
  "S3",                "S3",       "sepsis_patient",   "female", "70-75", 72.5,
  "survivor",          "survivor",
  "Survivor subject S3 from GEO labels and Table 1 sex/age combination.",
  "E.coli bacteremia", 19,         7,        NA_real_, NA_character_,
  "92",                "0.31",     "0.4",              "0.39",              "0.1"
) %>%
  mutate(
    plasma_resistin_ng_ml = safe_numeric(plasma_resistin_reported),
    plasma_il6_ng_ml = safe_numeric(plasma_il6_reported),
    plasma_il8_ng_ml = safe_numeric(plasma_il8_reported),
    plasma_il10_ng_ml = safe_numeric(plasma_il10_reported),
    lps_induced_tnf_alpha_ng_ml = safe_numeric(lps_induced_tnf_alpha_reported)
  )

supp_table_s1_tbl <- tibble::tribble(
  ~paper_sample, ~paper_b_cells, ~paper_cd14_mono, ~paper_cd4_t, ~paper_fcgr3a_mono,
  ~paper_cd8_t, ~paper_nk, ~paper_dc, ~paper_platelet, ~paper_erythroid_precursors,
  ~paper_neutrophil, ~paper_cmp, ~paper_cells_after_qc_total,
  "HC1",       1474L, 1091L, 2549L, 621L, 879L, 662L,  35L,   98L,  27L,   7L, 10L, 7453L,
  "HC2",       2462L, 1201L, 2010L, 830L, 520L, 656L, 206L,  126L,  64L,   3L, 15L, 8093L,
  "NS LS T0",   156L,  201L,  214L,  77L, 328L, 410L,  35L,  779L, 232L,   5L,  3L, 2440L,
  "NS LS T6",   194L,  441L,  258L, 259L, 205L, 260L, 131L, 1195L, 364L,   8L,  3L, 3318L,
  "NS ES_T0",   492L,  753L,  904L, 559L, 215L, 257L,  21L,  286L,  46L, 152L, 27L, 3712L,
  "NS ES_T6",   998L,  730L, 2049L, 549L, 272L, 382L,  54L,  272L,  59L, 145L, 36L, 5546L,
  "S1 T0",      147L,  929L,  337L, 444L, 195L, 229L,  23L,   97L,  42L,  25L,  0L, 2468L,
  "S1 T6",      514L, 1095L,  477L, 635L, 363L, 396L, 135L,  198L,  37L,  59L,  4L, 3913L,
  "S2_T0",     2744L,  278L, 1726L,  84L, 278L, 430L, 948L,   81L,  66L,  39L,  1L, 6675L,
  "S2_T6",     1061L,  755L, 1309L, 835L, 491L, 625L,  47L,   32L,  39L,  90L,  9L, 5293L,
  "S3_T0",      262L,  929L,  409L, 781L, 172L, 163L,  24L,   27L,   3L,  16L,  2L, 2788L,
  "S3_T6",      760L,  934L, 1577L, 887L, 638L, 455L,  80L,   69L,   7L,  26L,  1L, 5434L
)

# GEO-facing labels
geo_sample_tbl <- tibble::tribble(
  ~geo_accession, ~paper_sample, ~paper_sample_canonical, ~timepoint_label, ~hours_from_sepsis_recognition,
  "GSM5102900", "HC1",       "HC1",       "T0", NA_real_,
  "GSM5102901", "HC2",       "HC2",       "T0", NA_real_,
  "GSM5102902", "NS LS T0",  "NS_LS_T0",  "T0", 0,
  "GSM5102903", "NS LS T6",  "NS_LS_T6",  "T6", 6,
  "GSM5102904", "S1 T0",     "S1_T0",     "T0", 0,
  "GSM5102905", "S1 T6",     "S1_T6",     "T6", 6,
  "GSM5511351", "NS ES_T0",  "NS_ES_T0",  "T0", 0,
  "GSM5511352", "NS ES_T6",  "NS_ES_T6",  "T6", 6,
  "GSM5511353", "S2_T0",     "S2_T0",     "T0", 0,
  "GSM5511354", "S2_T6",     "S2_T6",     "T6", 6,
  "GSM5511355", "S3_T0",     "S3_T0",     "T0", 0,
  "GSM5511356", "S3_T6",     "S3_T6",     "T6", 6
) %>%
  mutate(
    stage_label_geo_raw = dplyr::case_when(
      stringr::str_detect(paper_sample_canonical, "^NS_LS_") ~ "NS_LS",
      stringr::str_detect(paper_sample_canonical, "^NS_ES_") ~ "NS_ES",
      stringr::str_detect(paper_sample_canonical, "^S[123]_") ~ "S",
      stringr::str_detect(paper_sample_canonical, "^HC") ~ "HC",
      TRUE ~ NA_character_
    )
  )

write_csv(paper_subject_tbl, file.path(manifest_dir, "GSE167363_paper_subject_embedded.csv"))
write_csv(supp_table_s1_tbl, file.path(manifest_dir, "GSE167363_supp_table_s1_embedded.csv"))
write_csv(geo_sample_tbl,    file.path(manifest_dir, "GSE167363_geo_sample_embedded.csv"))

# ------------------------------------------------------------------------------
# SRA normalization
# ------------------------------------------------------------------------------
sra_raw <- readr::read_csv(sra_run_table, show_col_types = FALSE)
names(sra_raw) <- normalize_names(names(sra_raw))

run_col        <- pick_column(sra_raw, exact = c("run", "run_accession"), regex = "(^|_)run(_|$)", label = "SRA run table")
exp_col        <- pick_column(sra_raw, exact = c("experiment", "experiment_accession"), regex = "(^|_)experiment(_|$)", label = "SRA run table")
biosample_col  <- pick_column(sra_raw, exact = c("biosample", "biosample_accession", "bio_sample"), regex = "(^|_)biosample(_|$)|(^|_)bio_sample(_|$)", label = "SRA run table")
geo_col        <- pick_column(sra_raw, exact = c("geo_accession_exp", "geo_accession"), regex = "geo.*accession", label = "SRA run table")

sample_name_col   <- pick_column(sra_raw, exact = c("sample_name"), regex = "^sample_name$", required = FALSE)
sample_name2_col  <- pick_column(sra_raw, exact = c("sample_name_1"), regex = "^sample_name_1$", required = FALSE)
sample_name_alt   <- pick_column(sra_raw, exact = c("sample_name_2"), regex = "^sample_name_2$", required = FALSE)
source_name_col   <- pick_column(sra_raw, exact = c("source_name"), regex = "^source_name$", required = FALSE)
age_range_col     <- pick_column(sra_raw, exact = c("age_range"), regex = "age.*range", required = FALSE)
sex_col           <- pick_column(sra_raw, exact = c("sex"), regex = "^sex$", required = FALSE)
disease_state_col <- pick_column(sra_raw, exact = c("disease_state"), regex = "disease.*state", required = FALSE)
donor_id_col      <- pick_column(sra_raw, exact = c("donor_id"), regex = "^donor_id$", required = FALSE)
donor_col         <- pick_column(sra_raw, exact = c("donor"), regex = "^donor$", required = FALSE)
avgspot_col       <- pick_column(sra_raw, exact = c("avgspotlen"), regex = "avgspotlen", required = FALSE)
bases_col         <- pick_column(sra_raw, exact = c("bases"), regex = "^bases$", required = FALSE)
bytes_col         <- pick_column(sra_raw, exact = c("bytes"), regex = "^bytes$", required = FALSE)
bioproject_col    <- pick_column(sra_raw, exact = c("bioproject"), regex = "^bioproject$", required = FALSE)
assay_type_col    <- pick_column(sra_raw, exact = c("assay_type"), regex = "assay.*type", required = FALSE)
instrument_col    <- pick_column(sra_raw, exact = c("instrument"), regex = "^instrument$", required = FALSE)
platform_col      <- pick_column(sra_raw, exact = c("platform"), regex = "^platform$", required = FALSE)
library_layout_col<- pick_column(sra_raw, exact = c("librarylayout"), regex = "librarylayout", required = FALSE)
library_sel_col   <- pick_column(sra_raw, exact = c("libraryselection"), regex = "libraryselection", required = FALSE)
library_src_col   <- pick_column(sra_raw, exact = c("librarysource"), regex = "librarysource", required = FALSE)
organism_col      <- pick_column(sra_raw, exact = c("organism"), regex = "^organism$", required = FALSE)
consent_col       <- pick_column(sra_raw, exact = c("consent"), regex = "^consent$", required = FALSE)
center_col        <- pick_column(sra_raw, exact = c("center_name"), regex = "center.*name", required = FALSE)
release_col       <- pick_column(sra_raw, exact = c("releasedate"), regex = "releasedate", required = FALSE)
create_col        <- pick_column(sra_raw, exact = c("create_date"), regex = "create_date", required = FALSE)
version_col       <- pick_column(sra_raw, exact = c("version"), regex = "^version$", required = FALSE)
tissue_col        <- pick_column(sra_raw, exact = c("tissue"), regex = "^tissue$", required = FALSE)
filetype_col      <- pick_column(sra_raw, exact = c("datastore_filetype"), regex = "datastore.*filetype", required = FALSE)
provider_col      <- pick_column(sra_raw, exact = c("datastore_provider"), regex = "datastore.*provider", required = FALSE)
region_col        <- pick_column(sra_raw, exact = c("datastore_region"), regex = "datastore.*region", required = FALSE)

sample_label_from_sra <- dplyr::coalesce(
  if (!is.na(sample_name_col))  as.character(sra_raw[[sample_name_col]])  else NA_character_,
  if (!is.na(sample_name2_col)) as.character(sra_raw[[sample_name2_col]]) else NA_character_,
  if (!is.na(sample_name_alt))  as.character(sra_raw[[sample_name_alt]])  else NA_character_
)

sra_tbl <- tibble(
  run_accession              = as.character(sra_raw[[run_col]]),
  experiment_accession       = as.character(sra_raw[[exp_col]]),
  biosample_accession        = as.character(sra_raw[[biosample_col]]),
  geo_accession              = as.character(sra_raw[[geo_col]]),
  sample_name_sra_primary    = sample_label_from_sra,
  source_name_sra            = if (!is.na(source_name_col)) as.character(sra_raw[[source_name_col]]) else NA_character_,
  age_range_sra              = if (!is.na(age_range_col)) as.character(sra_raw[[age_range_col]]) else NA_character_,
  sex_sra                    = if (!is.na(sex_col)) as.character(sra_raw[[sex_col]]) else NA_character_,
  disease_state_sra          = if (!is.na(disease_state_col)) as.character(sra_raw[[disease_state_col]]) else NA_character_,
  donor_id_sra               = if (!is.na(donor_id_col)) as.character(sra_raw[[donor_id_col]]) else NA_character_,
  donor_sra                  = if (!is.na(donor_col)) as.character(sra_raw[[donor_col]]) else NA_character_,
  assay_type_sra             = if (!is.na(assay_type_col)) as.character(sra_raw[[assay_type_col]]) else NA_character_,
  avg_spot_len_sra           = if (!is.na(avgspot_col)) safe_numeric(sra_raw[[avgspot_col]]) else NA_real_,
  bases_sra                  = if (!is.na(bases_col)) safe_numeric(sra_raw[[bases_col]]) else NA_real_,
  bytes_sra                  = if (!is.na(bytes_col)) safe_numeric(sra_raw[[bytes_col]]) else NA_real_,
  bioproject_sra             = if (!is.na(bioproject_col)) as.character(sra_raw[[bioproject_col]]) else NA_character_,
  instrument_sra             = if (!is.na(instrument_col)) as.character(sra_raw[[instrument_col]]) else NA_character_,
  platform_sra               = if (!is.na(platform_col)) as.character(sra_raw[[platform_col]]) else NA_character_,
  library_layout_sra         = if (!is.na(library_layout_col)) as.character(sra_raw[[library_layout_col]]) else NA_character_,
  library_selection_sra      = if (!is.na(library_sel_col)) as.character(sra_raw[[library_sel_col]]) else NA_character_,
  library_source_sra         = if (!is.na(library_src_col)) as.character(sra_raw[[library_src_col]]) else NA_character_,
  organism_sra               = if (!is.na(organism_col)) as.character(sra_raw[[organism_col]]) else NA_character_,
  consent_sra                = if (!is.na(consent_col)) as.character(sra_raw[[consent_col]]) else NA_character_,
  center_name_sra            = if (!is.na(center_col)) as.character(sra_raw[[center_col]]) else NA_character_,
  release_date_sra           = if (!is.na(release_col)) as.character(sra_raw[[release_col]]) else NA_character_,
  create_date_sra            = if (!is.na(create_col)) as.character(sra_raw[[create_col]]) else NA_character_,
  version_sra                = if (!is.na(version_col)) as.character(sra_raw[[version_col]]) else NA_character_,
  tissue_sra                 = if (!is.na(tissue_col)) as.character(sra_raw[[tissue_col]]) else NA_character_,
  datastore_filetype_sra     = if (!is.na(filetype_col)) as.character(sra_raw[[filetype_col]]) else NA_character_,
  datastore_provider_sra     = if (!is.na(provider_col)) as.character(sra_raw[[provider_col]]) else NA_character_,
  datastore_region_sra       = if (!is.na(region_col)) as.character(sra_raw[[region_col]]) else NA_character_
) %>%
  mutate(
    sample_name_sra_primary = stringr::str_trim(sample_name_sra_primary),
    sample_name_sra_canonical = canonical_label(sample_name_sra_primary),
    age_sra_midpoint = midpoint_from_range(age_range_sra),
    sex_sra = stringr::str_to_lower(sex_sra)
  )

write_csv(sra_tbl, file.path(manifest_dir, "GSE167363_sra_table_normalized.csv"))

# ------------------------------------------------------------------------------
# Merge SRA <-> GEO labels <-> paper identities
# ------------------------------------------------------------------------------
sample_meta <- sra_tbl %>%
  left_join(geo_sample_tbl, by = "geo_accession") %>%
  mutate(
    sample_name_matches_geo_label = canonical_label(sample_name_sra_primary) == canonical_label(paper_sample),
    subject_mapping_key = dplyr::case_when(
      paper_sample %in% c("HC1", "HC2") ~ paper_sample,
      stringr::str_detect(paper_sample, "^S1") ~ "S1",
      stringr::str_detect(paper_sample, "^S2") ~ "S2",
      stringr::str_detect(paper_sample, "^S3") ~ "S3",
      stringr::str_detect(paper_sample, "^NS LS") & sex_sra == "male"   & age_range_sra == "90-95" ~ "NS_male_90_95",
      stringr::str_detect(paper_sample, "^NS ES") & sex_sra == "female" & age_range_sra == "65-70" ~ "NS_female_65_70",
      TRUE ~ NA_character_
    )
  ) %>%
  left_join(paper_subject_tbl, by = c("subject_mapping_key" = "subject_id_standard")) %>%
  left_join(supp_table_s1_tbl, by = "paper_sample") %>%
  mutate(
    database_accession = project_id,
    sra_study = sra_study,
    bioproject = coalesce(bioproject_sra, bioproject),
    organism = coalesce(organism_sra, organism),
    sample_id = geo_accession,
    sample_alias = coalesce(paper_sample, sample_name_sra_primary, geo_accession),

    # Atlas grouping fields
    study_group = dplyr::case_when(
      subject_type == "healthy_control" ~ "healthy_control",
      outcome_label == "survivor" ~ "sepsis_survivor",
      outcome_label == "non_survivor" ~ "sepsis_nonsurvivor",
      TRUE ~ NA_character_
    ),
    factor = dplyr::case_when(
      subject_type == "healthy_control" ~ "healthy_control",
      TRUE ~ "gram_negative_bacterial_sepsis"
    ),
    cov19 = "NEG",
    is_sepsis = dplyr::case_when(
      subject_type == "healthy_control" ~ "no",
      TRUE ~ "yes"
    ),
    ards = NA_character_,
    pneumonia = NA_character_,
    infection_site = dplyr::case_when(
      subject_type == "healthy_control" ~ NA_character_,
      TRUE ~ "blood / systemic"
    ),
    pathogen_etiology = sepsis_etiology,
    mortality = dplyr::case_when(
      outcome_label == "non_survivor" ~ "yes",
      outcome_label == "survivor" ~ "no",
      TRUE ~ NA_character_
    ),
    icu_admission = dplyr::case_when(
      subject_type == "healthy_control" ~ "no",
      TRUE ~ "yes"
    ),
    days = NA_real_,
    severity = dplyr::case_when(
      subject_type == "healthy_control" ~ "healthy_control",
      outcome_label == "survivor" ~ "survivor_outcome",
      outcome_label == "non_survivor" ~ "fatal_outcome",
      TRUE ~ NA_character_
    ),
    severityatday = coalesce(stage_label_geo_raw, manuscript_stage_inferred),
    sex = coalesce(sex_sra, sex_paper),
    age = coalesce(age_sra_midpoint, age_midpoint_paper),
    ventilation = NA_character_,

    # Provenance from SRA/GEO/raw labels vs inferred manuscript identity
    sample_name_geo_or_sra = sample_name_sra_primary,
    donor_id = coalesce(donor_id_sra, donor_sra),
    donor_label_sra = donor_sra,
    patient_id = patient_id_inferred,
    patient_id_inference_basis = dplyr::case_when(
      subject_type == "healthy_control" ~ "direct HC label from GEO/SRA",
      stringr::str_detect(patient_id_inferred, "^S[123]$") ~ "direct survivor label from GEO/SRA",
      patient_id_inferred %in% c("P34", "P50") ~ "inferred from SRA sex+age_range plus manuscript text",
      TRUE ~ NA_character_
    ),
    patient_id_inference_confidence = dplyr::case_when(
      patient_id_inferred %in% c("P34", "P50") ~ "high_but_indirect",
      TRUE ~ "direct"
    ),
    stage_label_manuscript_inferred = manuscript_stage_inferred,
    label_conflict_flag = dplyr::case_when(
      !is.na(stage_label_geo_raw) & !is.na(stage_label_manuscript_inferred) &
        stage_label_geo_raw != stage_label_manuscript_inferred ~ TRUE,
      TRUE ~ FALSE
    ),
    label_conflict_note = dplyr::case_when(
      label_conflict_flag ~ paste0(
        "GEO/SRA sample label indicates ", stage_label_geo_raw,
        " whereas manuscript identity inference indicates ", stage_label_manuscript_inferred,
        " for patient_id=", patient_id_inferred, ". Both are retained."
      ),
      TRUE ~ NA_character_
    ),

    # Sample handling / technology
    source_name = coalesce(source_name_sra, "PBMC"),
    tissue_label = "PBMC",
    source_blood_fraction = "PBMC",
    blood_collection_tube = "Heparin vacutainer glass",
    pbmc_isolation_method = "Histopaque-1077 density gradient centrifugation",
    ficoll = "no (Histopaque-1077 used)",
    rbc_lysis = "no",
    sample_storage = "PBMC isolated within 24 h, plasma recovered, aliquots frozen in liquid nitrogen, thawed before 10x",
    frozen_or_fresh = "cryopreserved_PBMC_then_thawed",
    platform = "10x Genomics",
    sequencer_model = coalesce(instrument_sra, "Illumina NovaSeq 6000"),
    library_type = "10x Genomics 5prime/3prime GEX not explicitly stated in GEO; manuscript says Chromium Next GEM Single Cell v3.1 cDNA libraries",
    chemistry_version = "Chromium Next GEM Single Cell v3.1",
    feature_barcoding = "none_reported",
    vdj_capture = "no",
    read_len = NA_character_,
    aligner = "Cell Ranger/STAR (paper original pipeline provenance)",
    cellranger_or_equivalent_version = "3.1.0 (paper provenance; atlas raw data were reprocessed independently)",
    reference_genome = NA_character_,
    gene_annotation = NA_character_,
    seq_depth = NA_character_,
    assay_type = assay_type_sra,
    library_layout = library_layout_sra,
    library_selection = library_selection_sra,
    library_source = library_source_sra,
    n_loaded_cells = 15000L,
    paper_loaded_cells_per_sample = 15000L,

    # Paper-side QC provenance
    paper_qc_min_genes_per_cell = 200L,
    paper_qc_max_genes_per_cell = 6000L,
    paper_qc_min_umi_per_cell = 1000L,
    paper_qc_max_percent_mt = 20,
    paper_qc_max_percent_ribo = NA_real_,
    paper_qc_max_percent_hb = NA_real_,

    enrollment_window = "Sepsis recognition (0 h) and 6 h after diagnosis/resuscitation window; controls single time point",
    batch_unit_recommended = "sample",
    source_note = paste(
      "Atlas raw matrices were generated by the user's common Cell Ranger pipeline from public raw reads.",
      "Paper-side Cell Ranger/10x/QC metadata are stored as provenance only, not as the atlas processing definition.",
      "Metadata embedded from manuscript main text/Table 1/Methods, GEO labels, SRA Run Table, and Supplementary Table S1."
    )
  )

# sanity
if (nrow(sample_meta) != 12) {
  stop("Expected 12 sample rows for GSE167363, got ", nrow(sample_meta))
}
if (any(is.na(sample_meta$paper_sample))) {
  stop("Some rows failed to map GEO accession -> paper sample.")
}

# run metadata here is effectively run-level metadata; one run per sample
run_meta <- sample_meta %>%
  mutate(run_title = sample_name_sra_primary)

# explicit output column set
common_meta_columns <- c(
  "database_accession", "sra_study", "bioproject", "organism",
  "biosample_accession", "experiment_accession", "run_accession", "geo_accession",
  "sample_id", "sample_alias", "paper_sample", "sample_name_geo_or_sra",
  "patient_id", "patient_id_inference_basis", "patient_id_inference_confidence",
  "donor_id", "donor_label_sra", "subject_type",
  "study_group", "outcome_label", "factor", "cov19", "is_sepsis", "ards", "pneumonia",
  "infection_site", "pathogen_etiology", "mortality", "icu_admission",
  "days", "timepoint_label", "hours_from_sepsis_recognition",
  "severity", "severityatday", "ventilation",
  "sex", "age_range_paper", "age_range_sra", "age", "age_midpoint_paper", "age_sra_midpoint",
  "qsofa", "apache_ii", "sofa", "time_of_death_days_post_enrollment_reported",
  "plasma_resistin_reported", "plasma_resistin_ng_ml",
  "plasma_il6_reported", "plasma_il6_ng_ml",
  "plasma_il8_reported", "plasma_il8_ng_ml",
  "plasma_il10_reported", "plasma_il10_ng_ml",
  "lps_induced_tnf_alpha_reported", "lps_induced_tnf_alpha_ng_ml",
  "paper_b_cells", "paper_cd14_mono", "paper_cd4_t", "paper_fcgr3a_mono",
  "paper_cd8_t", "paper_nk", "paper_dc", "paper_platelet",
  "paper_erythroid_precursors", "paper_neutrophil", "paper_cmp",
  "paper_cells_after_qc_total",
  "stage_label_geo_raw", "stage_label_manuscript_inferred", "label_conflict_flag", "label_conflict_note",
  "sample_name_matches_geo_label",
  "source_name", "tissue_label", "source_blood_fraction",
  "blood_collection_tube", "pbmc_isolation_method", "ficoll", "rbc_lysis",
  "sample_storage", "frozen_or_fresh",
  "platform", "sequencer_model", "library_type", "chemistry_version",
  "feature_barcoding", "vdj_capture", "read_len",
  "aligner", "cellranger_or_equivalent_version",
  "reference_genome", "gene_annotation", "seq_depth",
  "assay_type", "library_layout", "library_selection", "library_source",
  "n_loaded_cells", "paper_loaded_cells_per_sample",
  "paper_qc_min_genes_per_cell", "paper_qc_max_genes_per_cell",
  "paper_qc_min_umi_per_cell", "paper_qc_max_percent_mt",
  "paper_qc_max_percent_ribo", "paper_qc_max_percent_hb",
  "avg_spot_len_sra", "bases_sra", "bytes_sra",
  "bioproject_sra", "instrument_sra", "platform_sra", "organism_sra",
  "consent_sra", "center_name_sra", "release_date_sra", "create_date_sra", "version_sra",
  "disease_state_sra", "source_name_sra", "tissue_sra",
  "datastore_filetype_sra", "datastore_provider_sra", "datastore_region_sra",
  "enrollment_window", "batch_unit_recommended", "source_note",
  "read_mode", "predicted.cluster"
)

sample_meta <- ensure_columns(sample_meta, common_meta_columns)
run_meta    <- ensure_columns(run_meta, common_meta_columns)

sample_meta <- atlas_standardize_metadata(sample_meta, project_id = project_id)
run_meta <- atlas_standardize_metadata(run_meta, project_id = project_id)
write_csv(sample_meta, file.path(meta_dir, "GSE167363_sample_metadata.csv"))
write_csv(run_meta,    file.path(meta_dir, "GSE167363_run_metadata.csv"))

# ------------------------------------------------------------------------------
# Read run-level matrices and build Seurat objects
# ------------------------------------------------------------------------------
input_dirs_all <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
input_dirs <- input_dirs_all[basename(input_dirs_all) %in% run_meta$run_accession]
input_dirs <- input_dirs[order(match(basename(input_dirs), run_meta$run_accession))]

if (length(input_dirs) == 0) {
  stop(
    "No run directories under input_root matched run metadata.\n",
    "Expected: ", paste(run_meta$run_accession, collapse = ", "), "\n",
    "Found: ", paste(head(basename(input_dirs_all), 50), collapse = ", ")
  )
}

message("Sample metadata rows : ", nrow(sample_meta))
message("Matched run dirs     : ", length(input_dirs))

sample_meta <- atlas_standardize_metadata(sample_meta, project_id = project_id)

obj_list <- list()
qc_summary <- tibble(
  run_accession = character(),
  biosample_accession = character(),
  geo_accession = character(),
  paper_sample = character(),
  patient_id = character(),
  n_raw = integer(),
  n_after_qc = integer(),
  saved = logical(),
  preqc_rds_path = character(),
  qc_rds_path = character(),
  working_rds_path = character()
)

for (run_dir in input_dirs) {
  run_id <- basename(run_dir)
  raw_dir <- file.path(run_dir, "outs", "raw_feature_bc_matrix")
  if (!dir.exists(raw_dir)) {
    message("Skipping missing raw_feature_bc_matrix: ", raw_dir)
    next
  }

  meta_row <- sample_meta %>% filter(run_accession == run_id)
  if (nrow(meta_row) != 1) {
    stop("Sample metadata row is not unique for ", run_id)
  }

  message("Reading ", run_id, " -> ", meta_row$geo_accession, " (", meta_row$paper_sample, ")")
  res <- read_gene_expression(raw_dir, sample_id = run_id)

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
    nFeature_RNA <= paper_qc_max_genes_per_cell &
    nCount_RNA > paper_qc_min_umi_per_cell &
    percent.mt <= paper_qc_max_percent_mt
  )

  obj[["atlas_qc_pass"]] <- with(obj@meta.data,
    nFeature_RNA >= 200 &
    percent.mt < if (apply_atlas_mt10_filter) 10 else 20
  )

  n_raw <- ncol(obj)
  obj <- atlas_standardize_seurat_metadata(obj, project_id = project_id)

  preqc_rds <- file.path(
    rds_preqc_dir,
    paste0(meta_row$geo_accession, "__", run_id, "__preqc_raw.rds")
  )
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
        run_accession = run_id,
        biosample_accession = meta_row$biosample_accession,
        geo_accession = meta_row$geo_accession,
        paper_sample = meta_row$paper_sample,
        patient_id = meta_row$patient_id,
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

  qc_rds <- file.path(
    rds_qc_dir,
    paste0(meta_row$geo_accession, "__", run_id, "__qcfiltered_raw.rds")
  )
  safe_save_rds(obj_qc, qc_rds)

  obj_work <- NormalizeData(obj_qc, verbose = FALSE)
  obj_work <- FindVariableFeatures(obj_work, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

  work_rds <- file.path(
    rds_work_dir,
    paste0(meta_row$geo_accession, "__", run_id, "__working.rds")
  )
  safe_save_rds(obj_work, work_rds)

  obj_list[[as.character(meta_row$geo_accession)]] <- obj_work

  qc_summary <- bind_rows(
    qc_summary,
    tibble(
      run_accession = run_id,
      biosample_accession = meta_row$biosample_accession,
      geo_accession = meta_row$geo_accession,
      paper_sample = meta_row$paper_sample,
      patient_id = meta_row$patient_id,
      n_raw = n_raw,
      n_after_qc = n_after_qc,
      saved = TRUE,
      preqc_rds_path = preqc_rds,
      qc_rds_path = qc_rds,
      working_rds_path = work_rds
    )
  )
}

write_csv(qc_summary, file.path(qc_dir, "GSE167363_qc_summary.csv"))

registry <- qc_summary %>%
  filter(saved) %>%
  left_join(
    sample_meta %>%
      select(
        run_accession, biosample_accession, geo_accession, paper_sample, patient_id,
        study_group, outcome_label, timepoint_label, hours_from_sepsis_recognition,
        stage_label_geo_raw, stage_label_manuscript_inferred,
        label_conflict_flag, label_conflict_note
      ),
    by = c("run_accession", "biosample_accession", "geo_accession", "paper_sample", "patient_id")
  ) %>%
  mutate(
    project_id = project_id,
    preintegration_tier = "working"
  )

write_csv(registry, file.path(output_root, "GSE167363_preintegration_registry.csv"))
saveRDS(obj_list, file.path(output_root, "GSE167363_preintegration_objlist.rds"))

message("Done. Saved pre-QC raw objects to    : ", rds_preqc_dir)
message("Done. Saved QC-filtered raw objects to: ", rds_qc_dir)
message("Done. Saved working objects to       : ", rds_work_dir)
message("Registry: ", file.path(output_root, "GSE167363_preintegration_registry.csv"))
