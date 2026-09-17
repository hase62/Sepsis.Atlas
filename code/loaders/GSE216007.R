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
# GSE216007 pre-integration loader
#
# Dataset:
#   GSE216007 / SRP576936 / PRJNA891700
#   HSPC_ECG multiome subseries from:
#   "Neutrophils and emergency granulopoiesis drive immune suppression and an
#    extreme response endotype during sepsis"
#
# Important atlas-side decision:
#   The user has already reprocessed public raw reads through the common atlas
#   pipeline. Paper-side Cell Ranger ARC / ArchR / vireo / QC settings are stored
#   here as provenance metadata, not as the definition of the atlas pipeline.
#
# Important study-specific decision:
#   GEO RNA samples are gPlexA_RNA ... gPlexE_RNA. Each represents a multiplexed
#   pool containing HSPCs from sepsis patients and healthy volunteers, not a single
#   clinical subject. Therefore this script creates one Seurat object per gPlex
#   RNA pool and explicitly marks the objects as pooled / not demultiplexed.
#   Individual patient-level labels must be added later only if a valid cell-level
#   demultiplexing table is available.
#
# Expected Cell Ranger output directory names:
#   cellranger_count/output/SRP576936/SAMN47814517/outs/raw_feature_bc_matrix
#   cellranger_count/output/SRP576936/SAMN47814519/outs/raw_feature_bc_matrix
#   cellranger_count/output/SRP576936/SAMN47814583/outs/raw_feature_bc_matrix
#   cellranger_count/output/SRP576936/SAMN47814584/outs/raw_feature_bc_matrix
#   cellranger_count/output/SRP576936/SAMN47815104/outs/raw_feature_bc_matrix
#
# External required input:
#   configs/GSE216007/GSE216007_SraRunTable.csv
#   or configs/GSE216007/SraRunTable.csv
# The table is validated to match SRP576936 / PRJNA891700. GEO-derived fields
# are embedded in this script, but SRA/run-level fields are intentionally read
# from the external run table, following the previous atlas loaders.
# ==============================================================================

project_id <- "GSE216007"
sra_study  <- "SRP576936"
bioproject <- "PRJNA891700"
organism   <- "Homo sapiens"

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------
`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

normalize_names <- function(x) {
  out <- x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "") %>%
    tolower()
  make.unique(out, sep = "_")
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
  if (is.na(col) || !col %in% names(df)) {
    n <- nrow(df)
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

bool_str <- function(x) {
  ifelse(is.na(x), NA_character_, ifelse(x, "yes", "no"))
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
      "Missing required SRA run table. Tried:\n",
      paste0(" - ", candidates, collapse = "\n"),
      "\nPlace the correct GSE216007 / SRP576936 SraRunTable.csv under configs/GSE216007/."
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
apply_paper_rna_qc_filter <- TRUE
min_cells_after_qc        <- 200L

# ------------------------------------------------------------------------------
# Embedded GEO / paper / supplement mapping
# ------------------------------------------------------------------------------
# RNA samples are the five GEX components actually loaded from raw_feature_bc_matrix.
geo_rna_tbl <- tibble::tribble(
  ~geo_accession, ~paper_sample, ~gplex_id, ~biosample_accession, ~experiment_accession_from_geo, ~modality, ~geo_title,    ~paired_atac_geo_accession,
  "GSM6654340",   "gPlexA_RNA",  "gPlexA",  "SAMN47814583",      "SRX28270570",                 "RNA",     "gPlexA_RNA",  "GSM6654345",
  "GSM6654341",   "gPlexB_RNA",  "gPlexB",  "SAMN47814584",      "SRX28270571",                 "RNA",     "gPlexB_RNA",  "GSM6654346",
  "GSM6654342",   "gPlexC_RNA",  "gPlexC",  "SAMN47814517",      "SRX28270572",                 "RNA",     "gPlexC_RNA",  "GSM6654347",
  "GSM6654343",   "gPlexD_RNA",  "gPlexD",  "SAMN47814519",      "SRX28270573",                 "RNA",     "gPlexD_RNA",  "GSM6654348",
  "GSM6654344",   "gPlexE_RNA",  "gPlexE",  "SAMN47815104",      "SRX28270574",                 "RNA",     "gPlexE_RNA",  "GSM6654349"
) %>%
  mutate(
    geo_status = "Public on Oct 16, 2024",
    geo_source_name = "PBMC",
    geo_characteristics_tissue = "PBMC",
    geo_characteristics_cell_type = "Hematopoietic stem and progenitor cells",
    geo_characteristics_disease_state = "sepsis patients and healthy volunteers",
    geo_extracted_molecule = "nuclear RNA",
    geo_extraction_protocol = "PBMCs were isolated from patients and healthy volunteers by Ficoll density gradient centrifugation. Cells were cryopreserved in 10% DMSO + 90% FCS and thawed for use on the day of the experiment. Cells were enriched for CD34+ via MACS twice, before being FACS sorted for CD34+CD45+ live singlets, and then nuclei were isolated for 10X multiome.",
    geo_protocol = "10X multiome protocol; scRNA-seq + scATAC-seq",
    geo_library_strategy = "RNA-Seq",
    geo_library_source = "transcriptomic single cell",
    geo_library_selection = "cDNA",
    geo_instrument_model = "Illumina NovaSeq 6000",
    geo_description = "Sample represents both sepsis patients and healthy controls.",
    geo_processed_object = "gex_obj.rds",
    geo_data_processing = paste(
      "Raw FASTQ files of scRNA-seq and scATAC-seq were aligned to the GRCh38 reference genome using 10X-arc (v2).",
      "Genetic demultiplexing and doublet removal was performed for each batch with vireo.",
      "HSPC multi-omic data were input into ArchR (v1.0.2) with mintss=4 and minFrags=1000.",
      "Homotypic doublets were removed with removeDoublet.",
      "Cells expressing <100 or >6000 genes, >25000 UMIs or with log10(UMI/gene) <0.8 were removed.",
      "Cells with TSS enrichment <7 and <1000 unique fragments were filtered out.",
      "Genes expressed in <10 cells or <3 total counts were removed.",
      "HSPC identity was assigned by mapping scRNA-seq data to Hao et al. 2021 and Granja et al. 2019 bone-marrow reference datasets.",
      "Non-HSPC labels from either mapping were filtered out.",
      "Multimodal dimensionality reduction used iterativeLSI for scRNA and scATAC, Harmony batch correction, addCombinedDims, and cluster majority RNA mapping.",
      "Each individual sample was pseudobulked for MACS2 peak calling and iterative peak overlap removal within ArchR.",
      "Dimensionality reduction and batch correction were repeated on gene expression and peak matrices before HSC reclustering.",
      sep = " "
    ),
    geo_assembly = "GRCh38",
    geo_supplementary_files_format_content_rna = "RNA count matrix (single cell experiment object)",
    geo_supplementary_files_format_content_atac = "ATAC peak matrix (single cell experiment object)",
    geo_series = "GSE216007; GSE216011"
  )

# ATAC samples are tracked in manifest only. They are not read by this RNA loader.
# Note that ATAC BioSamples are NOT the same as the paired RNA BioSamples.
geo_atac_tbl <- tibble::tribble(
  ~geo_accession, ~paper_sample, ~gplex_id, ~biosample_accession, ~experiment_accession_from_geo, ~modality, ~geo_title,      ~paired_rna_geo_accession,
  "GSM6654345",   "gPlexA_ATAC", "gPlexA",  "SAMN47815022",      "SRX28270575",                 "ATAC",    "gPlexA_ATAC",  "GSM6654340",
  "GSM6654346",   "gPlexB_ATAC", "gPlexB",  "SAMN47814693",      "SRX28270576",                 "ATAC",    "gPlexB_ATAC",  "GSM6654341",
  "GSM6654347",   "gPlexC_ATAC", "gPlexC",  "SAMN47815105",      "SRX28270577",                 "ATAC",    "gPlexC_ATAC",  "GSM6654342",
  "GSM6654348",   "gPlexD_ATAC", "gPlexD",  "SAMN47815106",      "SRX28270578",                 "ATAC",    "gPlexD_ATAC",  "GSM6654343",
  "GSM6654349",   "gPlexE_ATAC", "gPlexE",  "SAMN47815107",      "SRX28270579",                 "ATAC",    "gPlexE_ATAC",  "GSM6654344"
) %>%
  mutate(
    geo_status = "Public on Oct 16, 2024",
    geo_source_name = "PBMC",
    geo_characteristics_tissue = "PBMC",
    geo_characteristics_cell_type = "Hematopoietic stem and progenitor cells",
    geo_characteristics_disease_state = "sepsis patients and healthy volunteers",
    geo_extracted_molecule = "genomic DNA",
    geo_extraction_protocol = "PBMCs were isolated from patients and healthy volunteers by Ficoll density gradient centrifugation. Cells were cryopreserved in 10% DMSO + 90% FCS and thawed for use on the day of the experiment. Cells were enriched for CD34+ via MACS twice, before being FACS sorted for CD34+CD45+ live singlets, and then nuclei were isolated for 10X multiome.",
    geo_protocol = "10X multiome protocol; scRNA-seq + scATAC-seq",
    geo_library_strategy = "ATAC-seq",
    geo_library_source = "genomic single cell",
    geo_library_selection = "other",
    geo_instrument_model = "Illumina NextSeq 500",
    geo_description = "Sample represents both sepsis patients and healthy controls.",
    geo_processed_object = "peak_obj.rds",
    geo_data_processing = geo_rna_tbl$geo_data_processing[[1]],
    geo_assembly = "GRCh38",
    geo_supplementary_files_format_content_rna = "RNA count matrix (single cell experiment object)",
    geo_supplementary_files_format_content_atac = "ATAC peak matrix (single cell experiment object)",
    geo_series = "GSE216007; GSE216011"
  )

paper_scHSPC_tbl <- tibble::tribble(
  ~field, ~value, ~note,
  "paper_title", "Neutrophils and emergency granulopoiesis drive immune suppression and an extreme response endotype during sepsis", "Nature Immunology 2023",
  "series_title", "Immature neutrophil subsets and emergency granulopoiesis drive sepsis immune suppression and a specific extreme response to infection [HSPC_ECG]", "GSE216007",
  "super_series", "GSE216011", "GSE216007 is a SubSeries",
  "scWB_atlas", "272,993 cells; n=39 individuals", "whole-blood scRNA-seq + cell-surface protein profiling",
  "scWB_groups", "HC n=6; cardiac surgery control n=7; sepsis n=26", "Fig. 1 and Methods",
  "scHSPC_acute_sepsis_n", "15", "CD34+ HSPC isolation cohort described in Methods",
  "scHSPC_healthy_control_n", "7", "age- and sex-matched healthy controls",
  "scHSPC_convalescent_input_n", "8", "included before post-processing exclusions",
  "scHSPC_convalescent_removed_n", "3", "removed after preprocessing/demultiplexing because samples were >6 months after discharge by retrospective clinical evaluation",
  "scHSPC_final_individuals_n", "27", "15 acute sepsis + 7 healthy controls + 5 retained convalescents",
  "scHSPC_batches_reported", "6", "each batch contained at least one comparator group",
  "cell_source", "PBMC", "PBMCs isolated from sepsis acute/convalescent and healthy-control whole blood",
  "pbmc_isolation", "density gradient centrifugation with Leucosep tubes and lymphoprep", "Methods",
  "cryopreservation", "10% dimethylsulfoxide", "Methods",
  "enrichment", "CD34+ magnetic activated cell sorting twice; FACS live singlet CD34+CD45+ HSPCs", "Methods",
  "nuclei_protocol", "10x Genomics Demonstrated Protocol CG000365 Rev B low-input workflow", "Methods",
  "library_protocol", "10x Genomics Multiome RNA + ATAC library preparation, 1000285", "Methods",
  "rna_sequencer_paper", "Illumina NovaSeq 6000", "GEO RNA sample metadata",
  "atac_sequencer_paper", "Illumina NextSeq 500", "GEO ATAC sample metadata",
  "paper_original_aligner", "10x Cell Ranger ARC v2", "Paper/GEO provenance only; atlas objects are from user common raw-data processing",
  "paper_original_reference", "GRCh38", "Paper/GEO provenance",
  "paper_demux", "vireo genetic demultiplexing", "Paper/GEO provenance",
  "paper_doublet_removal", "vireo doublet removal + ArchR removeDoublet for homotypic doublets", "Paper/GEO provenance",
  "paper_qc_min_genes_per_cell", "100", "RNA-side paper/GEO QC provenance",
  "paper_qc_max_genes_per_cell", "6000", "RNA-side paper/GEO QC provenance",
  "paper_qc_max_umi_per_cell", "25000", "RNA-side paper/GEO QC provenance",
  "paper_qc_min_log10_umi_per_gene", "0.8", "RNA-side paper/GEO QC provenance",
  "paper_qc_min_tss_enrichment", "7", "ATAC-side paper/GEO QC, not applied to RNA-only atlas objects",
  "paper_qc_min_unique_fragments", "1000", "ATAC-side paper/GEO QC, not applied to RNA-only atlas objects",
  "paper_genes_after_filter", "26660", "paper processed dataset after filtering",
  "paper_hspcs_after_non_hspc_filter", "46156", "paper processed dataset after reference mapping",
  "paper_hscs_for_downstream", "29336", "paper processed dataset after progenitor exclusion",
  "supplementary_table_1_upload_note", "Uploaded Supplementary Table 1 contains cohort-level scWB discovery and mmV validation summaries, not gPlex-level or demultiplexed HSPC metadata.", "Stored below as cohort-level manifest tables"
)

supp_table1_discovery_scWB_tbl <- tibble::tribble(
  ~variable, ~sepsis_patients_n26, ~cardiac_surgery_n7, ~healthy_controls_n6,
  "Age", "63 (17)", "60 (16)", "63 (22)",
  "Male sex", "16 (62%)", "5 (71%)", "3 (50%)",
  "SOFA score", "4.8 (3.6)", "-", "-",
  "Mortality (14d)", "3 (12%)", "-", "-",
  "Mortality (28d)", "3 (12%)", "-", "-",
  "Infection source: Community acquired pneumonia", "9 (35%)", "-", "-",
  "Infection source: Urosepsis", "5 (19%)", "-", "-",
  "Infection source: Intra-abdominal sepsis", "8 (31%)", "-", "-",
  "Infection source: Meningitis", "1 (4%)", "-", "-",
  "Infection source: Necrotising fasciitis", "1 (4%)", "-", "-",
  "Infection source: Infective endocarditis", "1 (4%)", "-", "-",
  "Infection source: Bacteraemia", "1 (4%)", "-", "-",
  "Microbiology: E. coli", "4 (15%)", "-", "-",
  "Microbiology: E. coli + Clostridium perfringens", "1 (4%)", "-", "-",
  "Microbiology: S. pneumoniae", "1 (4%)", "-", "-",
  "Microbiology: S. aureus", "2 (8%)", "-", "-",
  "Microbiology: Proteus mirabilis", "2 (8%)", "-", "-",
  "Microbiology: Coagulase negative staphylococcus", "1 (4%)", "-", "-",
  "Microbiology: Pseudomonas aeruginosa", "1 (4%)", "-", "-",
  "Microbiology: Citrobacter", "1 (4%)", "-", "-",
  "Microbiology: Unknown", "13 (50%)", "-", "-",
  "Intensive care at point of sampling", "14 (54%)", "-", "-",
  "Mechanical ventilation", "5 (19%)", "-", "-",
  "Vasopressors", "8 (31%)", "-", "-",
  "Renal replacement therapy", "1 (4%)", "-", "-",
  "White cell count", "15.5 (6.4)", "-", "-",
  "Proportion neutrophils", "0.82", "-", "-",
  "Proportion lymphocytes", "0.09", "-", "-",
  "Proportion monocytes", "0.06", "-", "-",
  "Proportion eosinophils", "0.01", "-", "-",
  "Days from hospital admission until sampling", "2 (1.5)", "-", "-",
  "Footnote", "Data are n (%) or mean (SD) unless otherwise specified. SOFA=Sequential Organ Failure Assessment on day of sampling.", NA_character_, NA_character_
)

supp_table1_validation_mmV_tbl <- tibble::tribble(
  ~variable, ~sepsis_patients_n36, ~healthy_controls_n11,
  "Age", "67 (19)", "62 (12)",
  "Male sex", "21 (58%)", "6 (55%)",
  "SOFA score", "4 (3)", "-",
  "Mortality (14d)", "9 (25%)", "-",
  "Mortality (28d)", "9 (25%)", "-",
  "Infection source: Community acquired pneumonia", "15 (42%)", "-",
  "Infection source: Urosepsis", "9 (25%)", "-",
  "Infection source: Intra-abdominal sepsis", "8 (22%)", "-",
  "Infection source: Meningitis", "2 (6%)", "-",
  "Infection source: Necrotising fasciitis", "2 (6%)", "-",
  "Intensive care at point of sampling", "23 (64%)", NA_character_,
  "Mechanical ventilation", "14 (39%)", "-",
  "Vasopressors", "15 (42%)", "-",
  "Renal replacement therapy", "5 (14%)", "-",
  "White cell count", "15.8 (9.6)", "-",
  "Proportion neutrophils", "0.86", NA_character_,
  "Proportion lymphocytes", "0.08", NA_character_,
  "Proportion monocytes", "0.05", NA_character_,
  "Proportion eosinophils", "0.01", NA_character_,
  "Footnote", "Data are n (%) or mean (SD) unless otherwise specified. SOFA=Sequential Organ Failure Assessment on day of sampling.", NA_character_
)

write_csv(geo_rna_tbl,     file.path(manifest_dir, "GSE216007_geo_rna_embedded.csv"))
write_csv(geo_atac_tbl,    file.path(manifest_dir, "GSE216007_geo_atac_embedded.csv"))
write_csv(paper_scHSPC_tbl,file.path(manifest_dir, "GSE216007_paper_scHSPC_embedded.csv"))
write_csv(supp_table1_discovery_scWB_tbl, file.path(manifest_dir, "GSE216007_supp_table1_discovery_scWB_embedded.csv"))
write_csv(supp_table1_validation_mmV_tbl, file.path(manifest_dir, "GSE216007_supp_table1_validation_mmV_embedded.csv"))

# ------------------------------------------------------------------------------
# Required SRA table normalization and validation
# ------------------------------------------------------------------------------
sra_norm <- tibble(
  run_accession = character(),
  experiment_accession = character(),
  biosample_accession = character(),
  geo_accession = character(),
  sra_sample_name = character(),
  sra_sample_alias = character(),
  assay_type_sra = character(),
  avg_spot_len_sra = numeric(),
  bases_sra = numeric(),
  bytes_sra = numeric(),
  bioproject_sra = character(),
  instrument_sra = character(),
  platform_sra = character(),
  library_layout_sra = character(),
  library_selection_sra = character(),
  library_source_sra = character(),
  organism_sra = character(),
  consent_sra = character(),
  center_name_sra = character(),
  release_date_sra = character(),
  create_date_sra = character(),
  version_sra = character(),
  source_name_sra = character(),
  tissue_sra = character(),
  cell_type_sra = character(),
  disease_state_sra = character(),
  datastore_filetype_sra = character(),
  datastore_provider_sra = character(),
  datastore_region_sra = character(),
  sra_study_sra = character()
)

if (!is.na(sra_run_table) && file.exists(sra_run_table)) {
  sra_raw <- read_csv(sra_run_table, show_col_types = FALSE)
  names(sra_raw) <- normalize_names(names(sra_raw))

  run_col        <- pick_column(sra_raw, exact = c("run", "run_accession"), regex = "(^|_)run(_|$)", required = TRUE, label = "SRA run table")
  exp_col        <- pick_column(sra_raw, exact = c("experiment", "experiment_accession"), regex = "(^|_)experiment(_|$)", required = TRUE, label = "SRA run table")
  biosample_col  <- pick_column(sra_raw, exact = c("biosample", "biosample_accession", "bio_sample"), regex = "(^|_)biosample(_|$)|(^|_)bio_sample(_|$)", required = TRUE, label = "SRA run table")
  geo_col        <- pick_column(sra_raw, exact = c("geo_accession_exp", "geo_accession"), regex = "geo.*accession")
  sample_name_col<- pick_column(sra_raw, exact = c("sample_name"), regex = "^sample_name$")
  sample_alias_col <- pick_column(sra_raw, exact = c("sample_name_1", "sample_alias"), regex = "sample.*name|sample.*alias")
  assay_col      <- pick_column(sra_raw, exact = c("assay_type"), regex = "assay.*type")
  avgspot_col    <- pick_column(sra_raw, exact = c("avgspotlen"), regex = "avgspotlen")
  bases_col      <- pick_column(sra_raw, exact = c("bases"), regex = "^bases$")
  bytes_col      <- pick_column(sra_raw, exact = c("bytes"), regex = "^bytes$")
  bioproject_col <- pick_column(sra_raw, exact = c("bioproject"), regex = "^bioproject$")
  instr_col      <- pick_column(sra_raw, exact = c("instrument"), regex = "^instrument$")
  platform_col   <- pick_column(sra_raw, exact = c("platform"), regex = "^platform$")
  liblayout_col  <- pick_column(sra_raw, exact = c("librarylayout"), regex = "librarylayout")
  libsel_col     <- pick_column(sra_raw, exact = c("libraryselection"), regex = "libraryselection")
  libsrc_col     <- pick_column(sra_raw, exact = c("librarysource"), regex = "librarysource")
  organism_col   <- pick_column(sra_raw, exact = c("organism"), regex = "^organism$")
  consent_col    <- pick_column(sra_raw, exact = c("consent"), regex = "^consent$")
  center_col     <- pick_column(sra_raw, exact = c("center_name"), regex = "center.*name")
  release_col    <- pick_column(sra_raw, exact = c("releasedate"), regex = "releasedate")
  create_col     <- pick_column(sra_raw, exact = c("create_date"), regex = "create_date")
  version_col    <- pick_column(sra_raw, exact = c("version"), regex = "^version$")
  source_col     <- pick_column(sra_raw, exact = c("source_name"), regex = "^source_name$")
  tissue_col     <- pick_column(sra_raw, exact = c("tissue"), regex = "^tissue$")
  cell_type_col  <- pick_column(sra_raw, exact = c("cell_type"), regex = "cell.*type")
  disease_col    <- pick_column(sra_raw, exact = c("disease_state"), regex = "disease.*state")
  filetype_col   <- pick_column(sra_raw, exact = c("datastore_filetype"), regex = "datastore.*filetype")
  provider_col   <- pick_column(sra_raw, exact = c("datastore_provider"), regex = "datastore.*provider")
  region_col     <- pick_column(sra_raw, exact = c("datastore_region"), regex = "datastore.*region")
  sra_study_col  <- pick_column(sra_raw, exact = c("sra_study"), regex = "^sra_study$")

  sra_norm <- tibble(
    run_accession = col_or_na(sra_raw, run_col),
    experiment_accession = col_or_na(sra_raw, exp_col),
    biosample_accession = col_or_na(sra_raw, biosample_col),
    geo_accession = col_or_na(sra_raw, geo_col),
    sra_sample_name = col_or_na(sra_raw, sample_name_col),
    sra_sample_alias = col_or_na(sra_raw, sample_alias_col),
    assay_type_sra = col_or_na(sra_raw, assay_col),
    avg_spot_len_sra = col_or_na(sra_raw, avgspot_col, "numeric"),
    bases_sra = col_or_na(sra_raw, bases_col, "numeric"),
    bytes_sra = col_or_na(sra_raw, bytes_col, "numeric"),
    bioproject_sra = col_or_na(sra_raw, bioproject_col),
    instrument_sra = col_or_na(sra_raw, instr_col),
    platform_sra = col_or_na(sra_raw, platform_col),
    library_layout_sra = col_or_na(sra_raw, liblayout_col),
    library_selection_sra = col_or_na(sra_raw, libsel_col),
    library_source_sra = col_or_na(sra_raw, libsrc_col),
    organism_sra = col_or_na(sra_raw, organism_col),
    consent_sra = col_or_na(sra_raw, consent_col),
    center_name_sra = col_or_na(sra_raw, center_col),
    release_date_sra = col_or_na(sra_raw, release_col),
    create_date_sra = col_or_na(sra_raw, create_col),
    version_sra = col_or_na(sra_raw, version_col),
    source_name_sra = col_or_na(sra_raw, source_col),
    tissue_sra = col_or_na(sra_raw, tissue_col),
    cell_type_sra = col_or_na(sra_raw, cell_type_col),
    disease_state_sra = col_or_na(sra_raw, disease_col),
    datastore_filetype_sra = col_or_na(sra_raw, filetype_col),
    datastore_provider_sra = col_or_na(sra_raw, provider_col),
    datastore_region_sra = col_or_na(sra_raw, region_col),
    sra_study_sra = col_or_na(sra_raw, sra_study_col)
  ) %>%
    mutate(
      across(where(is.character), ~na_if(.x, ""))
    )
}

if (nrow(sra_norm) == 0) {
  stop("SRA run table was read but no rows were parsed. Check the CSV content.")
}

valid_project <- any(sra_norm$bioproject_sra %in% bioproject, na.rm = TRUE)
valid_study <- any(sra_norm$sra_study_sra %in% sra_study, na.rm = TRUE)
valid_biosamples <- any(sra_norm$biosample_accession %in% geo_rna_tbl$biosample_accession, na.rm = TRUE)

if (!valid_project && !valid_study && !valid_biosamples) {
  stop(
    "The supplied SRA run table does not appear to be GSE216007/SRP576936.\n",
    "Expected BioProject=", bioproject, " or SRA Study=", sra_study,
    " or at least one RNA BioSample among: ", paste(geo_rna_tbl$biosample_accession, collapse = ", "), "\n",
    "Observed BioProject values: ", paste(unique(na.omit(sra_norm$bioproject_sra)), collapse = ", "), "\n",
    "Observed SRA Study values: ", paste(unique(na.omit(sra_norm$sra_study_sra)), collapse = ", "), "\n",
    "Observed first BioSamples: ", paste(head(unique(na.omit(sra_norm$biosample_accession)), 10), collapse = ", ")
  )
}

write_csv(sra_norm, file.path(manifest_dir, "GSE216007_sra_table_normalized.csv"))

# Prefer matched SRA rows, but do not require them because embedded GEO mapping
# includes the BioSample-level directory IDs.
sra_rna_match <- sra_norm %>%
  filter(
    geo_accession %in% geo_rna_tbl$geo_accession |
      biosample_accession %in% geo_rna_tbl$biosample_accession |
      sra_sample_name %in% geo_rna_tbl$paper_sample |
      sra_sample_alias %in% geo_rna_tbl$paper_sample
  ) %>%
  group_by(biosample_accession) %>%
  summarise(
    run_accession_sra = paste(unique(na.omit(run_accession)), collapse = ";"),
    experiment_accession_sra = paste(unique(na.omit(experiment_accession)), collapse = ";"),
    geo_accession_sra = paste(unique(na.omit(geo_accession)), collapse = ";"),
    assay_type_sra = paste(unique(na.omit(assay_type_sra)), collapse = ";"),
    avg_spot_len_sra = suppressWarnings(mean(avg_spot_len_sra, na.rm = TRUE)),
    bases_sra = suppressWarnings(sum(bases_sra, na.rm = TRUE)),
    bytes_sra = suppressWarnings(sum(bytes_sra, na.rm = TRUE)),
    bioproject_sra = paste(unique(na.omit(bioproject_sra)), collapse = ";"),
    instrument_sra = paste(unique(na.omit(instrument_sra)), collapse = ";"),
    platform_sra = paste(unique(na.omit(platform_sra)), collapse = ";"),
    library_layout_sra = paste(unique(na.omit(library_layout_sra)), collapse = ";"),
    library_selection_sra = paste(unique(na.omit(library_selection_sra)), collapse = ";"),
    library_source_sra = paste(unique(na.omit(library_source_sra)), collapse = ";"),
    organism_sra = paste(unique(na.omit(organism_sra)), collapse = ";"),
    consent_sra = paste(unique(na.omit(consent_sra)), collapse = ";"),
    center_name_sra = paste(unique(na.omit(center_name_sra)), collapse = ";"),
    source_name_sra = paste(unique(na.omit(source_name_sra)), collapse = ";"),
    tissue_sra = paste(unique(na.omit(tissue_sra)), collapse = ";"),
    cell_type_sra = paste(unique(na.omit(cell_type_sra)), collapse = ";"),
    disease_state_sra = paste(unique(na.omit(disease_state_sra)), collapse = ";"),
    sra_study_sra = paste(unique(na.omit(sra_study_sra)), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(
    across(where(is.character), ~na_if(.x, "")),
    avg_spot_len_sra = ifelse(is.infinite(avg_spot_len_sra) | is.nan(avg_spot_len_sra), NA_real_, avg_spot_len_sra),
    bases_sra = ifelse(is.infinite(bases_sra) | is.nan(bases_sra), NA_real_, bases_sra),
    bytes_sra = ifelse(is.infinite(bytes_sra) | is.nan(bytes_sra), NA_real_, bytes_sra)
  )

if (nrow(sra_rna_match) == 0) {
  stop(
    "SRA run table was read and passed broad study validation, but no rows matched the five GSE216007 RNA BioSamples/GSMs.\n",
    "Expected RNA BioSamples: ", paste(geo_rna_tbl$biosample_accession, collapse = ", "), "\n",
    "Expected RNA GSMs: ", paste(geo_rna_tbl$geo_accession, collapse = ", ")
  )
}

missing_sra_rna_biosamples <- setdiff(geo_rna_tbl$biosample_accession, sra_rna_match$biosample_accession)
extra_sra_rna_biosamples <- setdiff(sra_rna_match$biosample_accession, geo_rna_tbl$biosample_accession)

if (length(missing_sra_rna_biosamples) > 0 || length(extra_sra_rna_biosamples) > 0 || nrow(sra_rna_match) != nrow(geo_rna_tbl)) {
  stop(
    "SRA-to-GEO RNA matching is incomplete for GSE216007.\n",
    "Matched RNA BioSamples: ", paste(sra_rna_match$biosample_accession, collapse = ", "), "\n",
    "Missing expected RNA BioSamples: ", paste(missing_sra_rna_biosamples, collapse = ", "), "\n",
    "Unexpected matched BioSamples: ", paste(extra_sra_rna_biosamples, collapse = ", "), "\n",
    "Expected exactly five RNA gPlex BioSamples."
  )
}

if (any(is.na(sra_rna_match$run_accession_sra) | sra_rna_match$run_accession_sra == "")) {
  stop("SRA-to-GEO RNA matching succeeded, but one or more matched rows lack run_accession_sra.")
}

# ------------------------------------------------------------------------------
# Build sample metadata
# ------------------------------------------------------------------------------
sample_meta <- geo_rna_tbl %>%
  left_join(sra_rna_match, by = "biosample_accession") %>%
  mutate(
    database_accession = project_id,
    sra_study = sra_study,
    bioproject = ifelse(!is.na(bioproject_sra), bioproject_sra, bioproject),
    organism = ifelse(!is.na(organism_sra), organism_sra, organism),
    experiment_accession = dplyr::coalesce(experiment_accession_sra, experiment_accession_from_geo),
    # run_accession is the actual SRR or semicolon-joined SRRs from the external SRA table.
    # This is intentionally not inferred from the BioSample directory once strict SRA matching has passed.
    run_accession = run_accession_sra,
    cellranger_dir_id = biosample_accession,

    sample_id = geo_accession,
    sample_alias = paper_sample,

    # Dataset grain.
    pooled_flag = TRUE,
    pool_id = gplex_id,
    pool_members = NA_character_,
    demux_status = "not_demultiplexed",
    demux_required_for_subject_assignment = TRUE,
    demux_method_original = "vireo genetic demultiplexing in paper pipeline",
    atlas_cell_assignment_status = "pool-level only; no individual donor assignment applied",
    sample_represents = "multiplexed pool of HSPCs from sepsis patients and healthy volunteers",
    study_group = "mixed_hc_sepsis_pool",
    factor = "mixed_healthy_control_and_sepsis",
    cov19 = NA_character_,
    is_sepsis = "mixed",
    ards = NA_character_,
    pneumonia = NA_character_,
    infection_site = NA_character_,
    pathogen_etiology = NA_character_,
    mortality = NA_character_,
    icu_admission = NA_character_,
    days = NA_real_,
    timepoint_label = "mixed_or_not_applicable",
    severity = NA_character_,
    severityatday = NA_character_,
    ventilation = NA_character_,
    sex = NA_character_,
    age = NA_real_,
    age_range = NA_character_,
    qsofa = NA_real_,
    sofa = NA_real_,

    # Cohort-level information from paper.
    paper_scHSPC_acute_sepsis_n = 15L,
    paper_scHSPC_healthy_control_n = 7L,
    paper_scHSPC_convalescent_input_n = 8L,
    paper_scHSPC_convalescent_removed_n = 3L,
    paper_scHSPC_final_individuals_n = 27L,
    paper_scHSPC_hspcs_after_non_hspc_filter = 46156L,
    paper_scHSPC_hscs_downstream_n = 29336L,
    paper_scHSPC_batches_reported = 6L,
    supp_table1_discovery_scWB_available = TRUE,
    supp_table1_validation_mmV_available = TRUE,
    supp_table1_note = "Supplementary Table 1 is cohort-level for scWB/mmV and is stored as manifest tables; it is not gPlex-level HSPC demultiplexed metadata.",

    # Explicit GEO sample fields.
    geo_sample_status = geo_status,
    geo_source_name = geo_source_name,
    geo_characteristics_tissue = geo_characteristics_tissue,
    geo_characteristics_cell_type = geo_characteristics_cell_type,
    geo_characteristics_disease_state = geo_characteristics_disease_state,
    geo_extracted_molecule = geo_extracted_molecule,
    geo_extraction_protocol = geo_extraction_protocol,
    geo_protocol = geo_protocol,
    geo_library_strategy = geo_library_strategy,
    geo_library_source = geo_library_source,
    geo_library_selection = geo_library_selection,
    geo_instrument_model = geo_instrument_model,
    geo_description = geo_description,
    geo_processed_object = geo_processed_object,
    geo_data_processing = geo_data_processing,
    geo_assembly = geo_assembly,
    geo_supplementary_files_format_content_rna = geo_supplementary_files_format_content_rna,
    geo_supplementary_files_format_content_atac = geo_supplementary_files_format_content_atac,

    # Sample handling and cell source.
    source_name = ifelse(!is.na(source_name_sra), source_name_sra, geo_source_name),
    tissue_label = "PBMC",
    source_blood_fraction = "PBMC",
    cell_type = "hematopoietic stem and progenitor cells",
    extracted_molecule = "nuclear RNA",
    blood_collection_tube = NA_character_,
    pbmc_isolation_method = "density gradient centrifugation with Leucosep tubes and lymphoprep",
    ficoll = "lymphoprep density gradient",
    rbc_lysis = "no",
    sample_storage = "PBMC cryopreserved in 10% DMSO; thawed for CD34+ HSPC enrichment and nuclei isolation",
    frozen_or_fresh = "cryopreserved PBMC, thawed",
    enrichment_method = "CD34+ MACS enrichment twice, followed by FACS sorting for live singlet CD34+CD45+ HSPCs",
    nuclei_isolation = "10x Genomics low-input nuclei workflow CG000365 Rev B",
    sorted_population = "live singlet CD34+CD45+ HSPCs",

    # Technology / provenance.
    platform = "10x Genomics Multiome",
    sequencer_model = "Illumina NovaSeq 6000",
    paired_atac_sequencer_model = "Illumina NextSeq 500",
    library_type = "10x Genomics single-cell multiome GEX",
    chemistry_version = "10x Genomics Multiome RNA+ATAC",
    feature_barcoding = "ATAC paired in same multiome assay; RNA object only here",
    vdj_capture = "no",
    read_len = NA_character_,
    assay_type = dplyr::coalesce(assay_type_sra, "RNA-Seq"),
    library_layout = dplyr::coalesce(library_layout_sra, NA_character_),
    library_selection = dplyr::coalesce(library_selection_sra, "cDNA"),
    library_source = dplyr::coalesce(library_source_sra, "TRANSCRIPTOMIC"),

    # Original paper pipeline provenance only.
    aligner = "Cell Ranger / user common atlas pipeline for this object; paper original used Cell Ranger ARC v2",
    cellranger_or_equivalent_version = "atlas common reprocessing; paper provenance: 10x Cell Ranger ARC v2",
    paper_original_aligner = "10x Cell Ranger ARC v2",
    paper_original_reference_genome = "GRCh38",
    reference_genome = "GRCh38",
    gene_annotation = NA_character_,
    seq_depth = NA_character_,

    # RNA-side paper QC provenance. ATAC-side QC is recorded but not applied here.
    n_loaded_cells = NA_integer_,
    paper_loaded_cells_per_sample = NA_integer_,
    paper_qc_min_genes_per_cell = 100L,
    paper_qc_max_genes_per_cell = 6000L,
    paper_qc_max_umi_per_cell = 25000L,
    paper_qc_min_log10_umi_per_gene = 0.8,
    paper_qc_max_percent_mt = NA_real_,
    paper_qc_max_percent_ribo = NA_real_,
    paper_qc_max_percent_hb = NA_real_,
    paper_qc_min_tss_enrichment = 7,
    paper_qc_min_unique_fragments = 1000,

    avg_spot_len_sra = avg_spot_len_sra,
    bases_sra = bases_sra,
    bytes_sra = bytes_sra,
    instrument_sra = instrument_sra,
    platform_sra = platform_sra,
    consent_sra = consent_sra,
    center_name_sra = center_name_sra,
    source_name_sra = source_name_sra,
    tissue_sra = tissue_sra,
    cell_type_sra = cell_type_sra,
    disease_state_sra = disease_state_sra,

    enrollment_window = "mixed pool; individual acute/convalescent/control assignment requires demultiplexing",
    batch_unit_recommended = "gPlex pool before demultiplexing; individual donor after demultiplexing",
    source_note = paste(
      "GSE216007 RNA samples are five gPlex RNA pools from CD34+CD45+ HSPC multiome.",
      "Each GEO RNA sample represents both sepsis patients and healthy controls; no patient-level metadata are assigned at this stage.",
      "ATAC samples are tracked in manifest only. Atlas RNA matrices were generated by user common reprocessing."
    ),
    read_mode = NA_character_,
    predicted.cluster = NA_character_
  )

pool_manifest <- sample_meta %>%
  select(
    database_accession, gplex_id, pool_id, geo_accession, paper_sample,
    biosample_accession, cellranger_dir_id, paired_atac_geo_accession,
    sample_represents, pooled_flag, demux_status,
    demux_required_for_subject_assignment, demux_method_original,
    paper_scHSPC_acute_sepsis_n, paper_scHSPC_healthy_control_n,
    paper_scHSPC_convalescent_input_n, paper_scHSPC_final_individuals_n,
    source_note
  )

write_csv(pool_manifest, file.path(manifest_dir, "GSE216007_pool_manifest.csv"))

# run metadata has same grain here because each Cell Ranger output directory is one BioSample/gPlex RNA pool.
run_meta <- sample_meta

common_meta_columns <- c(
  "database_accession", "sra_study", "bioproject", "organism",
  "biosample_accession", "experiment_accession", "run_accession", "geo_accession",
  "cellranger_dir_id", "sample_id", "sample_alias", "paper_sample", "geo_title",
  "gplex_id", "pool_id", "pool_members", "modality", "paired_atac_geo_accession",
  "pooled_flag", "demux_status", "demux_required_for_subject_assignment",
  "demux_method_original", "atlas_cell_assignment_status", "sample_represents",
  "study_group", "factor", "cov19", "is_sepsis", "ards", "pneumonia",
  "infection_site", "pathogen_etiology", "mortality", "icu_admission",
  "days", "timepoint_label", "severity", "severityatday", "ventilation",
  "sex", "age", "age_range", "qsofa", "sofa",
  "paper_scHSPC_acute_sepsis_n", "paper_scHSPC_healthy_control_n",
  "paper_scHSPC_convalescent_input_n", "paper_scHSPC_convalescent_removed_n",
  "paper_scHSPC_final_individuals_n", "paper_scHSPC_hspcs_after_non_hspc_filter",
  "paper_scHSPC_hscs_downstream_n", "paper_scHSPC_batches_reported",
  "supp_table1_discovery_scWB_available", "supp_table1_validation_mmV_available", "supp_table1_note",
  "geo_sample_status", "geo_source_name", "geo_characteristics_tissue", "geo_characteristics_cell_type",
  "geo_characteristics_disease_state", "geo_extracted_molecule", "geo_extraction_protocol", "geo_protocol",
  "geo_library_strategy", "geo_library_source", "geo_library_selection", "geo_instrument_model",
  "geo_description", "geo_processed_object", "geo_data_processing", "geo_assembly",
  "geo_supplementary_files_format_content_rna", "geo_supplementary_files_format_content_atac",
  "source_name", "tissue_label", "source_blood_fraction", "cell_type",
  "extracted_molecule", "blood_collection_tube", "pbmc_isolation_method",
  "ficoll", "rbc_lysis", "sample_storage", "frozen_or_fresh",
  "enrichment_method", "nuclei_isolation", "sorted_population",
  "platform", "sequencer_model", "paired_atac_sequencer_model",
  "library_type", "chemistry_version", "feature_barcoding", "vdj_capture",
  "read_len", "assay_type", "library_layout", "library_selection", "library_source",
  "aligner", "cellranger_or_equivalent_version", "paper_original_aligner",
  "paper_original_reference_genome", "reference_genome", "gene_annotation", "seq_depth",
  "n_loaded_cells", "paper_loaded_cells_per_sample",
  "paper_qc_min_genes_per_cell", "paper_qc_max_genes_per_cell",
  "paper_qc_max_umi_per_cell", "paper_qc_min_log10_umi_per_gene",
  "paper_qc_max_percent_mt", "paper_qc_max_percent_ribo", "paper_qc_max_percent_hb",
  "paper_qc_min_tss_enrichment", "paper_qc_min_unique_fragments",
  "avg_spot_len_sra", "bases_sra", "bytes_sra",
  "instrument_sra", "platform_sra", "consent_sra", "center_name_sra",
  "source_name_sra", "tissue_sra", "cell_type_sra", "disease_state_sra", "sra_study_sra",
  "enrollment_window", "batch_unit_recommended", "source_note",
  "read_mode", "predicted.cluster"
)

sample_meta <- ensure_columns(sample_meta, common_meta_columns)
run_meta    <- ensure_columns(run_meta, common_meta_columns)

sample_meta <- atlas_standardize_metadata(sample_meta, project_id = project_id)
run_meta <- atlas_standardize_metadata(run_meta, project_id = project_id)
write_csv(sample_meta, file.path(meta_dir, "GSE216007_sample_metadata_gplexRNA.csv"))
write_csv(run_meta,    file.path(meta_dir, "GSE216007_run_metadata_gplexRNA.csv"))

# ------------------------------------------------------------------------------
# Load each BioSample-level matrix and build objects
# ------------------------------------------------------------------------------
input_dirs_all <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
input_dirs <- input_dirs_all[basename(input_dirs_all) %in% sample_meta$cellranger_dir_id]
input_dirs <- input_dirs[order(match(basename(input_dirs), sample_meta$cellranger_dir_id))]

if (nrow(sample_meta) != 5) {
  stop("Expected 5 GSE216007 RNA gPlex rows, got ", nrow(sample_meta))
}
if (length(input_dirs) == 0) {
  stop(
    "No Cell Ranger directories under input_root matched metadata.\n",
    "Expected BioSample dirs: ", paste(sample_meta$cellranger_dir_id, collapse = ", "), "\n",
    "Available under input_root: ", paste(head(basename(input_dirs_all), 50), collapse = ", ")
  )
}
if (length(input_dirs) != nrow(sample_meta)) {
  warning(
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
  gplex_id = character(),
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
  if (!dir.exists(raw_dir)) {
    message("Skipping missing raw_feature_bc_matrix: ", raw_dir)
    next
  }

  meta_row <- sample_meta %>% filter(cellranger_dir_id == dir_id)
  if (nrow(meta_row) != 1) {
    stop("Sample metadata row is not unique for cellranger_dir_id ", dir_id)
  }

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
  # Paper provenance describes a log10 genes-per-UMI / UMI-per-gene style QC metric.
  # Implement this as log10(nUMI) / log10(nGenes), not log10(nUMI / nGenes).
  # The latter is too stringent for sparse single-cell data and can remove all cells.
  obj[["log10_umi_per_gene_ratio"]] <- log10(pmax(obj$nCount_RNA, 1)) / log10(pmax(obj$nFeature_RNA, 2))

  # RNA-only subset of the paper's original HSPC QC.
  obj[["paper_rna_qc_pass"]] <- with(obj@meta.data,
    nFeature_RNA >= paper_qc_min_genes_per_cell &
      nFeature_RNA <= paper_qc_max_genes_per_cell &
      nCount_RNA <= paper_qc_max_umi_per_cell &
      log10_umi_per_gene_ratio >= paper_qc_min_log10_umi_per_gene
  )

  # Minimal atlas sanity flag; actual atlas-wide validrop/doublet handling can be
  # done downstream from rds_preqc_raw if desired.
  obj[["atlas_basic_qc_pass"]] <- with(obj@meta.data,
    nFeature_RNA >= 100
  )

  n_raw <- ncol(obj)
  obj <- atlas_standardize_seurat_metadata(obj, project_id = project_id)

  preqc_rds <- file.path(rds_preqc_dir, paste0(meta_row$geo_accession, "__", dir_id, "__preqc_raw.rds"))
  safe_save_rds(obj, preqc_rds)

  obj <- run_optional_helpers(obj)

  qc_subset <- rep(TRUE, ncol(obj))
  if (apply_paper_rna_qc_filter) qc_subset <- qc_subset & obj$paper_rna_qc_pass
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
        gplex_id = meta_row$gplex_id,
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

  obj_list[[as.character(meta_row$geo_accession)]] <- obj_work

  qc_summary <- bind_rows(
    qc_summary,
    tibble(
      cellranger_dir_id = dir_id,
      run_accession = meta_row$run_accession,
      biosample_accession = meta_row$biosample_accession,
      geo_accession = meta_row$geo_accession,
      paper_sample = meta_row$paper_sample,
      gplex_id = meta_row$gplex_id,
      n_raw = n_raw,
      n_after_qc = n_after_qc,
      saved = TRUE,
      preqc_rds_path = preqc_rds,
      qc_rds_path = qc_rds,
      working_rds_path = work_rds
    )
  )
}

write_csv(qc_summary, file.path(qc_dir, "GSE216007_qc_summary_gplexRNA.csv"))

registry <- qc_summary %>%
  filter(saved) %>%
  left_join(
    sample_meta %>%
      select(
        cellranger_dir_id, run_accession, biosample_accession, geo_accession,
        paper_sample, gplex_id, pool_id, study_group, sample_represents,
        pooled_flag, demux_status, demux_required_for_subject_assignment,
        batch_unit_recommended
      ),
    by = c("cellranger_dir_id", "run_accession", "biosample_accession", "geo_accession", "paper_sample", "gplex_id")
  ) %>%
  mutate(
    project_id = project_id,
    preintegration_tier = "working",
    object_grain = "gPlex RNA pool"
  )

write_csv(registry, file.path(output_root, "GSE216007_preintegration_registry_gplexRNA.csv"))
saveRDS(obj_list, file.path(output_root, "GSE216007_preintegration_objlist_gplexRNA.rds"))

message("Done. Saved pre-QC raw objects to    : ", rds_preqc_dir)
message("Done. Saved QC-filtered raw objects to: ", rds_qc_dir)
message("Done. Saved working objects to       : ", rds_work_dir)
message("Registry: ", file.path(output_root, "GSE216007_preintegration_registry_gplexRNA.csv"))
