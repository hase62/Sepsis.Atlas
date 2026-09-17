#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 2L){
  stop("Usage: script <root> <tag>")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]

atlas_root <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1"
)

out_dir <- file.path(
  atlas_root,
  paste0(
    "global_deferred_reconciliation_inventory_v1__",
    tag
  )
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

done_file <- file.path(
  out_dir,
  "GLOBAL_DEFERRED_RECONCILIATION_INVENTORY_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop("Already completed: ", done_file)
}

# ============================================================
# I/O
# ============================================================

read_tsv <- function(path){

  if(!file.exists(path)){
    stop("Missing file: ", path)
  }

  con <- if(grepl("\\.gz$", path)){
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }

  tryCatch(
    read.delim(
      con,
      sep="\t",
      quote="\"",
      comment.char="",
      stringsAsFactors=FALSE,
      check.names=FALSE
    ),
    finally=close(con)
  )
}

write_gz_tsv <- function(x, path){

  con <- gzfile(path, "wt")

  tryCatch(
    write.table(
      x,
      con,
      sep="\t",
      quote=TRUE,
      qmethod="double",
      row.names=FALSE
    ),
    finally=close(con)
  )

  status <- system2(
    "gzip",
    c("-t", path),
    stdout=FALSE,
    stderr=FALSE
  )

  if(status != 0L){
    stop("gzip integrity failure: ", path)
  }
}

latest_dir <- function(marker){

  files <- list.files(
    atlas_root,
    pattern=paste0("^", marker, "$"),
    recursive=TRUE,
    full.names=TRUE
  )

  if(!length(files)){
    stop("Completion marker not found: ", marker)
  }

  info <- file.info(files)

  files <- files[
    order(
      info$mtime,
      decreasing=TRUE
    )
  ]

  dirname(files[[1]])
}

single_file <- function(dir, pattern){

  f <- list.files(
    dir,
    pattern=pattern,
    full.names=TRUE
  )

  if(length(f) != 1L){
    stop(
      "Expected exactly one file in ",
      dir,
      " pattern=",
      pattern,
      " found=",
      length(f)
    )
  }

  f[[1]]
}

# ============================================================
# Resolve canonical freeze directories
# ============================================================

b_dir <- latest_dir(
  "B_PLASMA_FULL_ANNOTATION_FREEZE_COMPLETE.ok"
)

tnk_dir <- latest_dir(
  "ANNOTATION_FREEZE_COMPLETE.ok"
)

# The generic marker exists for both TNK and MonoDC, so
# explicitly resolve their known canonical directories.
tnk_dir <- file.path(
  atlas_root,
  "tnk_annotation_freeze_v1"
)

monodc_dir <- file.path(
  atlas_root,
  "monodc_annotation_freeze_v1"
)

neut_dir <- latest_dir(
  "NEUTROPHIL_ANNOTATION_FREEZE_COMPLETE.ok"
)

platelet_dir <- latest_dir(
  "PLATELET_FULL_ANNOTATION_FREEZE_COMPLETE.ok"
)

ery_dir <- latest_dir(
  "ERYTHROID_FULL_ANNOTATION_FREEZE_COMPLETE.ok"
)

prog_dir <- latest_dir(
  "PROGENITOR_FULL_ANNOTATION_FREEZE_COMPLETE.ok"
)

stopifnot(
  dir.exists(tnk_dir),
  dir.exists(monodc_dir)
)

# ============================================================
# Read canonical frozen annotations
# ============================================================

b <- readRDS(
  file.path(
    b_dir,
    "b_plasma_full_annotation_freeze_v1.rds"
  )
)

stopifnot(
  is.data.frame(b),
  nrow(b) == 45955L
)

tnk <- read_tsv(
  file.path(
    tnk_dir,
    "tnk_cell_annotation_v1.tsv.gz"
  )
)

monodc <- read_tsv(
  file.path(
    monodc_dir,
    "monodc_cell_annotation_v1.tsv.gz"
  )
)

neut <- read_tsv(
  single_file(
    neut_dir,
    "^neutrophil_cell_annotation_freeze_v1__.*\\.tsv\\.gz$"
  )
)

platelet <- read_tsv(
  single_file(
    platelet_dir,
    "^platelet_full_cell_annotation_v1__.*\\.tsv\\.gz$"
  )
)

ery <- read_tsv(
  single_file(
    ery_dir,
    "^erythroid_full_cell_annotation_v1__.*\\.tsv\\.gz$"
  )
)

prog <- read_tsv(
  single_file(
    prog_dir,
    "^progenitor_full_cell_annotation_v1__.*\\.tsv\\.gz$"
  )
)

# ============================================================
# Standardization helper
# ============================================================

make_std <- function(
  d,
  source_compartment,
  core_col,
  axis_col=NULL,
  state_col=NULL,
  confidence_col=NULL,
  scope_col=NULL
){

  req <- c(
    "global_cell",
    "project_id",
    "library_key",
    core_col
  )

  stopifnot(
    all(req %in% names(d))
  )

  out <- data.frame(
    global_cell=
      as.character(d$global_cell),

    project_id=
      as.character(d$project_id),

    library_key=
      as.character(d$library_key),

    source_compartment=
      source_compartment,

    frozen_core=
      as.character(d[[core_col]]),

    frozen_axis=
      if(!is.null(axis_col) &&
         axis_col %in% names(d))
        as.character(d[[axis_col]])
      else
        NA_character_,

    frozen_state=
      if(!is.null(state_col) &&
         state_col %in% names(d))
        as.character(d[[state_col]])
      else
        NA_character_,

    frozen_confidence=
      if(!is.null(confidence_col) &&
         confidence_col %in% names(d))
        as.character(d[[confidence_col]])
      else
        NA_character_,

    frozen_scope=
      if(!is.null(scope_col) &&
         scope_col %in% names(d))
        as.character(d[[scope_col]])
      else
        NA_character_,

    stringsAsFactors=FALSE
  )

  stopifnot(
    !anyNA(out$global_cell),
    !anyDuplicated(out$global_cell)
  )

  out
}

# ============================================================
# Standardize without changing taxonomy granularity
# ============================================================

std_b <- make_std(
  b,
  source_compartment="B_plasma",
  core_col="b_plasma_atlas_core_identity_v1",
  state_col="b_plasma_atlas_state_v1",
  confidence_col=
    "b_plasma_atlas_biological_evidence_confidence_v1",
  scope_col=
    "b_plasma_atlas_annotation_scope_v1"
)

std_tnk <- make_std(
  tnk,
  source_compartment="T_NK",
  core_col="tnk_core_v1",
  state_col="tnk_state_v1",
  confidence_col=
    "tnk_annotation_confidence_v1"
)

std_monodc <- make_std(
  monodc,
  source_compartment="Monocyte_DC",
  core_col="monodc_core_v1",
  state_col="monodc_state_v1",
  confidence_col=
    "monodc_core_confidence_v1"
)

std_neut <- make_std(
  neut,
  source_compartment="Neutrophil",
  core_col="neutrophil_core_v1",
  axis_col="neutrophil_maturation_v1",
  state_col="neutrophil_state_v1",
  confidence_col=
    "neutrophil_core_evidence_confidence_v1",
  scope_col=
    "neutrophil_annotation_scope_v1"
)

std_platelet <- make_std(
  platelet,
  source_compartment="Platelet_megakaryocyte",
  core_col="platelet_core_v1",
  axis_col=
    "platelet_transcriptional_axis_v1",
  state_col="platelet_state_v1",
  confidence_col=
    "platelet_core_evidence_confidence_v1",
  scope_col=
    "platelet_annotation_scope_v1"
)

std_ery <- make_std(
  ery,
  source_compartment="Erythroid",
  core_col="erythroid_core_v1",
  axis_col=
    "erythroid_maturation_axis_v1",
  state_col="erythroid_state_v1",
  confidence_col=
    "erythroid_evidence_confidence_v1"
)

std_prog <- make_std(
  prog,
  source_compartment="Progenitor",
  core_col="core_identity",
  axis_col="lineage_priming_axis",
  state_col="state",
  confidence_col=
    "annotation_confidence"
)

all_std <- rbind(
  std_b,
  std_tnk,
  std_monodc,
  std_neut,
  std_platelet,
  std_ery,
  std_prog
)

rownames(all_std) <- NULL

# ============================================================
# Global invariants
# ============================================================

within_dup <- do.call(
  rbind,
  lapply(
    split(
      all_std,
      all_std$source_compartment
    ),
    function(d){
      data.frame(
        source_compartment=
          d$source_compartment[[1]],
        n_cells=nrow(d),
        n_unique_global_cell=
          length(unique(d$global_cell)),
        n_duplicate_global_cell=
          sum(duplicated(d$global_cell)),
        stringsAsFactors=FALSE
      )
    }
  )
)

global_dup_ids <- unique(
  all_std$global_cell[
    duplicated(all_std$global_cell) |
      duplicated(
        all_std$global_cell,
        fromLast=TRUE
      )
  ]
)

global_dup <- if(length(global_dup_ids)){

  all_std[
    all_std$global_cell %in%
      global_dup_ids,
    ,
    drop=FALSE
  ]

} else {

  all_std[0, , drop=FALSE]
}

# ============================================================
# Frozen label census
# ============================================================

core_counts <- as.data.frame(
  table(
    source_compartment=
      all_std$source_compartment,
    frozen_core=
      all_std$frozen_core,
    useNA="ifany"
  ),
  stringsAsFactors=FALSE
)

core_counts <- core_counts[
  core_counts$Freq > 0,
  ,
  drop=FALSE
]

axis_counts <- as.data.frame(
  table(
    source_compartment=
      all_std$source_compartment,
    frozen_axis=
      all_std$frozen_axis,
    useNA="ifany"
  ),
  stringsAsFactors=FALSE
)

axis_counts <- axis_counts[
  axis_counts$Freq > 0,
  ,
  drop=FALSE
]

# ============================================================
# Deferred / out-of-compartment inventory
#
# This is deliberately broad for INVENTORY only.
# No cell is reassigned here.
# ============================================================

deferred_pattern <- paste(
  c(
    "^Deferred",
    "deferred",
    "non[_ -]?",
    "not_applicable"
  ),
  collapse="|"
)

core_text <- ifelse(
  is.na(all_std$frozen_core),
  "",
  all_std$frozen_core
)

scope_text <- ifelse(
  is.na(all_std$frozen_scope),
  "",
  all_std$frozen_scope
)

all_std$inventory_deferred_candidate <-
  grepl(
    deferred_pattern,
    core_text,
    ignore.case=TRUE
  ) |
  grepl(
    "deferred|non[_ -]?",
    scope_text,
    ignore.case=TRUE
  )

deferred <- all_std[
  all_std$inventory_deferred_candidate,
  ,
  drop=FALSE
]

# ============================================================
# Exact deferred-label counts
# ============================================================

deferred_counts <- as.data.frame(
  table(
    source_compartment=
      deferred$source_compartment,
    frozen_core=
      deferred$frozen_core,
    frozen_scope=
      deferred$frozen_scope,
    useNA="ifany"
  ),
  stringsAsFactors=FALSE
)

deferred_counts <- deferred_counts[
  deferred_counts$Freq > 0,
  ,
  drop=FALSE
]

# ============================================================
# Project distribution of deferred cells
# ============================================================

deferred_project_counts <- as.data.frame(
  table(
    source_compartment=
      deferred$source_compartment,
    frozen_core=
      deferred$frozen_core,
    project_id=
      deferred$project_id,
    useNA="ifany"
  ),
  stringsAsFactors=FALSE
)

deferred_project_counts <-
  deferred_project_counts[
    deferred_project_counts$Freq > 0,
    ,
    drop=FALSE
  ]

# ============================================================
# Source manifest
# ============================================================

source_manifest <- data.frame(
  source_compartment=c(
    "B_plasma",
    "T_NK",
    "Monocyte_DC",
    "Neutrophil",
    "Platelet_megakaryocyte",
    "Erythroid",
    "Progenitor"
  ),
  source_path=c(
    file.path(
      b_dir,
      "b_plasma_full_annotation_freeze_v1.rds"
    ),
    file.path(
      tnk_dir,
      "tnk_cell_annotation_v1.tsv.gz"
    ),
    file.path(
      monodc_dir,
      "monodc_cell_annotation_v1.tsv.gz"
    ),
    single_file(
      neut_dir,
      "^neutrophil_cell_annotation_freeze_v1__.*\\.tsv\\.gz$"
    ),
    single_file(
      platelet_dir,
      "^platelet_full_cell_annotation_v1__.*\\.tsv\\.gz$"
    ),
    single_file(
      ery_dir,
      "^erythroid_full_cell_annotation_v1__.*\\.tsv\\.gz$"
    ),
    single_file(
      prog_dir,
      "^progenitor_full_cell_annotation_v1__.*\\.tsv\\.gz$"
    )
  ),
  stringsAsFactors=FALSE
)

# ============================================================
# Write
# ============================================================

inventory_file <- file.path(
  out_dir,
  paste0(
    "atlas_frozen_annotation_inventory_v1__",
    tag,
    ".tsv.gz"
  )
)

write_gz_tsv(
  all_std,
  inventory_file
)

deferred_file <- file.path(
  out_dir,
  paste0(
    "atlas_deferred_candidate_cells_v1__",
    tag,
    ".tsv.gz"
  )
)

write_gz_tsv(
  deferred,
  deferred_file
)

write.table(
  core_counts,
  file.path(
    out_dir,
    "atlas_frozen_core_label_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  axis_counts,
  file.path(
    out_dir,
    "atlas_frozen_axis_label_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  deferred_counts,
  file.path(
    out_dir,
    "atlas_deferred_candidate_label_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  deferred_project_counts,
  file.path(
    out_dir,
    "atlas_deferred_candidate_project_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  within_dup,
  file.path(
    out_dir,
    "atlas_compartment_cell_uniqueness_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  global_dup,
  file.path(
    out_dir,
    "atlas_cross_compartment_duplicate_cells_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  source_manifest,
  file.path(
    out_dir,
    "atlas_deferred_reconciliation_source_manifest_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Re-read gzip outputs
# ============================================================

chk_inventory <- read_tsv(
  inventory_file
)

chk_deferred <- read_tsv(
  deferred_file
)

stopifnot(
  nrow(chk_inventory) ==
    nrow(all_std),
  identical(
    as.character(
      chk_inventory$global_cell
    ),
    as.character(
      all_std$global_cell
    )
  ),
  nrow(chk_deferred) ==
    nrow(deferred)
)

# ============================================================
# SHA256
# ============================================================

sha_files <- c(
  inventory_file,
  deferred_file,
  file.path(
    out_dir,
    "atlas_frozen_core_label_counts_v1.tsv"
  ),
  file.path(
    out_dir,
    "atlas_frozen_axis_label_counts_v1.tsv"
  ),
  file.path(
    out_dir,
    "atlas_deferred_candidate_label_counts_v1.tsv"
  ),
  file.path(
    out_dir,
    "atlas_deferred_candidate_project_counts_v1.tsv"
  ),
  file.path(
    out_dir,
    "atlas_compartment_cell_uniqueness_v1.tsv"
  ),
  file.path(
    out_dir,
    "atlas_cross_compartment_duplicate_cells_v1.tsv"
  ),
  file.path(
    out_dir,
    "atlas_deferred_reconciliation_source_manifest_v1.tsv"
  )
)

sha_file <- file.path(
  out_dir,
  paste0(
    "SHA256SUMS_v1__",
    tag,
    ".txt"
  )
)

status <- system2(
  "sha256sum",
  sha_files,
  stdout=sha_file
)

if(status != 0L){
  stop("SHA256 creation failed")
}

status <- system2(
  "sha256sum",
  c("-c", sha_file),
  stdout=FALSE,
  stderr=FALSE
)

if(status != 0L){
  stop("SHA256 verification failed")
}

# ============================================================
# Completion LAST
# ============================================================

writeLines(
  c(
    "PASS",
    "Global deferred reconciliation inventory v1",
    "",
    paste0(
      "n_frozen_annotation_rows=",
      nrow(all_std)
    ),
    paste0(
      "n_unique_global_cells=",
      length(unique(all_std$global_cell))
    ),
    paste0(
      "n_cross_compartment_duplicate_cell_ids=",
      length(global_dup_ids)
    ),
    paste0(
      "n_inventory_deferred_candidate_cells=",
      nrow(deferred)
    ),
    "",
    "annotation granularity unchanged",
    "no cell reassignment performed",
    "no new clustering",
    "no marker discovery",
    "frozen compartment annotations used as input",
    "cellranger_count not accessed",
    "",
    "gzip_integrity=PASS",
    "reread_count_order=PASS",
    "sha256_manifest=PASS"
  ),
  done_file
)

cat("\n===== COMPARTMENT CELL COUNTS =====\n")
print(within_dup, row.names=FALSE)

cat("\n===== CORE LABEL COUNTS =====\n")
print(core_counts, row.names=FALSE)

cat("\n===== DEFERRED CANDIDATE LABEL COUNTS =====\n")
print(deferred_counts, row.names=FALSE)

cat(
  "\nCross-compartment duplicate global_cell IDs = ",
  length(global_dup_ids),
  "\n",
  sep=""
)

cat(
  "Deferred candidate cells = ",
  nrow(deferred),
  "\n",
  sep=""
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: global deferred reconciliation inventory completed\n"
)
