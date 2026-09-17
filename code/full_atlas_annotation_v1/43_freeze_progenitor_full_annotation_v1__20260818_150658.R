#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 6L){
  stop(
    paste(
      "Usage: script <root> <tag>",
      "<discovery_dir> <loo_dir>",
      "<cluster_dir> <replication_dir>"
    )
  )
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]

discovery_dir <- normalizePath(args[[3]], mustWork=TRUE)
loo_dir <- normalizePath(args[[4]], mustWork=TRUE)
cluster_dir <- normalizePath(args[[5]], mustWork=TRUE)
replication_dir <- normalizePath(args[[6]], mustWork=TRUE)

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

source(
  file.path(
    root,
    "full_atlas_primary_integration_v1",
    "00_common.R"
  )
)

set.seed(20260818)

DISCOVERY_PROJECT <- "GSE216007"

N_DISCOVERY <- 26962L
N_QUERY <- 813L
N_FULL <- 27775L

N_HVG <- 3000L
N_PCS <- 30L

TARGET_AXIS_PRECISION <- 0.90
MIN_CALIBRATION_CALLS <- 100L

POSITIVE_AXES <- c(
  "lymphoid_priming",
  "granulocytic_priming",
  "megakaryocytic_priming",
  "megakaryocytic_erythroid_priming",
  "eos_baso_mast_priming"
)

# Intentionally NOT query-transferable.
NONTRANSFERABLE_AXIS <- "primitive_HSC_like"

# ============================================================
# Output
# ============================================================

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "progenitor_full_annotation_freeze_v1__",
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
  "PROGENITOR_FULL_ANNOTATION_FREEZE_COMPLETE.ok"
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

technical_gene <- function(g){

  grepl("^MT-", g) |
    grepl("^RPL[0-9]", g) |
    grepl("^RPS[0-9]", g) |
    g %in% c(
      "MALAT1",
      "NEAT1",
      "HBA1",
      "HBA2",
      "HBB",
      "HBD",
      "HBM"
    )
}

align_to_features <- function(m, features){

  out <- Matrix(
    0,
    nrow=length(features),
    ncol=ncol(m),
    sparse=TRUE,
    dimnames=list(
      features,
      colnames(m)
    )
  )

  present <- intersect(
    features,
    rownames(m)
  )

  if(length(present)){
    out[
      match(present, features),
      ,
      drop=FALSE
    ] <- m[
      present,
      ,
      drop=FALSE
    ]
  }

  attr(out, "n_features_present") <-
    length(present)

  out
}

# ============================================================
# Validate prerequisite stages
# ============================================================

prereq <- c(
  file.path(
    discovery_dir,
    "PROGENITOR_DISCOVERY_ANNOTATION_FREEZE_COMPLETE.ok"
  ),
  file.path(
    loo_dir,
    "PROGENITOR_LEAVE_ONE_LIBRARY_OUT_VALIDATION_COMPLETE.ok"
  ),
  file.path(
    cluster_dir,
    "PROGENITOR_DOMINANT_PROJECT_NATIVE_CLUSTERING_COMPLETE.ok"
  ),
  file.path(
    replication_dir,
    "PROGENITOR_CROSS_PROJECT_PROGRAM_REPLICATION_COMPLETE.ok"
  )
)

stopifnot(
  all(file.exists(prereq))
)

# ============================================================
# Discovery annotation
# ============================================================

disc_files <- list.files(
  discovery_dir,
  pattern=
    "^progenitor_discovery_cell_annotation_v1__.*\\.tsv\\.gz$",
  full.names=TRUE
)

stopifnot(
  length(disc_files) == 1L
)

disc <- read_tsv(
  disc_files[[1]]
)

stopifnot(
  nrow(disc) == N_DISCOVERY,
  !anyDuplicated(disc$global_cell)
)

for(nm in names(disc)){
  if(is.factor(disc[[nm]])){
    disc[[nm]] <- as.character(disc[[nm]])
  }
}

disc$global_cell <- as.character(disc$global_cell)
disc$project_id <- as.character(disc$project_id)
disc$library_key <- as.character(disc$library_key)

stopifnot(
  all(
    disc$project_id ==
      DISCOVERY_PROJECT
  )
)

disc$transfer_class <- paste(
  disc$core_identity,
  disc$lineage_priming_axis,
  sep="|||AXIS|||"
)

class_lookup <- unique(
  disc[
    ,
    c(
      "transfer_class",
      "core_identity",
      "lineage_priming_axis"
    ),
    drop=FALSE
  ]
)

stopifnot(
  !anyDuplicated(
    class_lookup$transfer_class
  )
)

# ============================================================
# Step42 LOO predictions: operational threshold calibration
#
# This is calibration only.
# The original unthresholded LOO metrics remain the
# Technical Validation result.
# ============================================================

loo_files <- list.files(
  loo_dir,
  pattern=
    "^progenitor_loo_predictions_v1__.*\\.tsv\\.gz$",
  full.names=TRUE
)

stopifnot(
  length(loo_files) == 1L
)

loo <- read_tsv(
  loo_files[[1]]
)

stopifnot(
  nrow(loo) == N_DISCOVERY,
  !anyDuplicated(loo$global_cell)
)

threshold_grid <- seq(
  0,
  0.99,
  by=0.01
)

cal_rows <- list()

for(axis in POSITIVE_AXES){

  selected <- NULL

  all_rows <- list()

  for(th in threshold_grid){

    idx <- which(
      loo$predicted_axis == axis &
        loo$class_prediction_score >= th
    )

    n_called <- length(idx)

    n_correct <- if(n_called){
      sum(
        loo$actual_axis[idx] == axis
      )
    } else {
      0L
    }

    precision <- if(n_called){
      n_correct / n_called
    } else {
      NA_real_
    }

    actual_n <- sum(
      loo$actual_axis == axis
    )

    recall <- if(actual_n){
      n_correct / actual_n
    } else {
      NA_real_
    }

    z <- data.frame(
      axis=axis,
      threshold=th,
      n_called=n_called,
      precision=precision,
      recall=recall,
      stringsAsFactors=FALSE
    )

    all_rows[[length(all_rows)+1L]] <- z

    if(
      is.null(selected) &&
      n_called >= MIN_CALIBRATION_CALLS &&
      is.finite(precision) &&
      precision >= TARGET_AXIS_PRECISION
    ){
      selected <- z
    }
  }

  audit_axis <- do.call(
    rbind,
    all_rows
  )

  cal_rows[[
    length(cal_rows)+1L
  ]] <- audit_axis

  if(is.null(selected)){

    selected <- data.frame(
      axis=axis,
      threshold=NA_real_,
      n_called=NA_integer_,
      precision=NA_real_,
      recall=NA_real_,
      stringsAsFactors=FALSE
    )
  }

  assign(
    paste0(
      "selected_",
      axis
    ),
    selected
  )
}

threshold_audit <- do.call(
  rbind,
  cal_rows
)

selected_thresholds <- do.call(
  rbind,
  lapply(
    POSITIVE_AXES,
    function(axis){
      get(
        paste0(
          "selected_",
          axis
        )
      )
    }
  )
)

rownames(selected_thresholds) <- NULL

selected_thresholds$policy <- ifelse(
  is.finite(
    selected_thresholds$threshold
  ),
  "score_threshold",
  ifelse(
    selected_thresholds$axis ==
      "granulocytic_priming",
    "native_program_gate_fallback",
    "not_directly_transferable"
  )
)

# Primitive HSC is deliberately not calibrated into a
# hard-transfer class because LOO showed poor specificity
# versus HSC-like library-confounded unresolved populations.

primitive_policy <- data.frame(
  axis=NONTRANSFERABLE_AXIS,
  threshold=NA_real_,
  n_called=NA_integer_,
  precision=NA_real_,
  recall=NA_real_,
  policy=
    "prediction_retained_as_candidate_only",
  stringsAsFactors=FALSE
)

policy <- rbind(
  selected_thresholds,
  primitive_policy
)

# ============================================================
# Query/support cells
# ============================================================

support_file <- file.path(
  cluster_dir,
  "progenitor_non_discovery_projects_replication_support_v1.tsv"
)

query_meta <- read_tsv(
  support_file
)

required_query <- c(
  "global_cell",
  "project_id",
  "library_key",
  "original_cell"
)

stopifnot(
  all(required_query %in% names(query_meta)),
  nrow(query_meta) == N_QUERY,
  !anyDuplicated(query_meta$global_cell)
)

for(nm in required_query){
  query_meta[[nm]] <- as.character(
    query_meta[[nm]]
  )
}

stopifnot(
  !any(
    query_meta$project_id ==
      DISCOVERY_PROJECT
  )
)

# Step40 native query program evidence.
query_program_files <- list.files(
  replication_dir,
  pattern=
    "^progenitor_replication_support_cell_programs_v1__.*\\.tsv\\.gz$",
  full.names=TRUE
)

stopifnot(
  length(query_program_files) == 1L
)

qprog <- read_tsv(
  query_program_files[[1]]
)

qprog$global_cell <-
  as.character(qprog$global_cell)

stopifnot(
  nrow(qprog) == N_QUERY,
  !anyDuplicated(qprog$global_cell)
)

qi <- match(
  query_meta$global_cell,
  qprog$global_cell
)

stopifnot(
  !anyNA(qi)
)

qprog <- qprog[
  qi,
  ,
  drop=FALSE
]

stopifnot(
  identical(
    query_meta$global_cell,
    qprog$global_cell
  )
)

# ============================================================
# Recover discovery original barcodes
# ============================================================

transfer_dir <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "annotation_transfer_by_library"
)

transfer_files <- sort(
  list.files(
    transfer_dir,
    pattern="\\.tsv\\.gz$",
    full.names=TRUE
  )
)

meta_list <- list()

for(f in transfer_files){

  d <- read_tsv(f)

  req <- c(
    "global_cell",
    "original_cell",
    "project_id",
    "library_key",
    "integration_compartment_primary_v1"
  )

  stopifnot(
    all(req %in% names(d))
  )

  d <- d[
    as.character(d$project_id) ==
      DISCOVERY_PROJECT &
      as.character(
        d$integration_compartment_primary_v1
      ) == "Progenitor",
    req,
    drop=FALSE
  ]

  if(nrow(d)){
    meta_list[[
      length(meta_list)+1L
    ]] <- d
  }
}

disc_meta <- do.call(
  rbind,
  meta_list
)

rm(meta_list)

for(nm in c(
  "global_cell",
  "original_cell",
  "project_id",
  "library_key"
)){
  disc_meta[[nm]] <- as.character(
    disc_meta[[nm]]
  )
}

stopifnot(
  nrow(disc_meta) == N_DISCOVERY,
  !anyDuplicated(disc_meta$global_cell)
)

di <- match(
  disc$global_cell,
  disc_meta$global_cell
)

stopifnot(
  !anyNA(di)
)

disc$original_cell <-
  disc_meta$original_cell[di]

stopifnot(
  disc$library_key ==
    disc_meta$library_key[di]
)

rm(disc_meta, di)
invisible(gc())

# ============================================================
# Source library lookup
# ============================================================

lib <- read_tsv(
  file.path(
    root,
    "pre_integration",
    "pilot_unintegrated_large_v1",
    "resolved_library_table.tsv"
  )
)

stopifnot(
  all(
    c(
      "project_id",
      "library_key",
      "final_rds_resolved"
    ) %in% names(lib)
  )
)

lib$project_id <- as.character(
  lib$project_id
)

lib$library_key <- as.character(
  lib$library_key
)

lib$key <- paste(
  lib$project_id,
  lib$library_key,
  sep="|||"
)

# ============================================================
# Build native discovery reference counts
# ============================================================

disc_libraries <- sort(
  unique(
    disc$library_key
  )
)

ref_list <- list()

for(ii in seq_along(disc_libraries)){

  lk <- disc_libraries[[ii]]

  key <- paste(
    DISCOVERY_PROJECT,
    lk,
    sep="|||"
  )

  li <- match(
    key,
    lib$key
  )

  if(is.na(li)){
    stop(
      "Discovery library lookup failed: ",
      key
    )
  }

  source_rds <- as.character(
    lib$final_rds_resolved[[li]]
  )

  d <- disc[
    disc$library_key == lk,
    ,
    drop=FALSE
  ]

  obj <- readRDS(source_rds)
  counts0 <- get_rna_counts(obj)

  stopifnot(
    all(
      d$original_cell %in%
        colnames(counts0)
    )
  )

  m <- counts0[
    ,
    d$original_cell,
    drop=FALSE
  ]

  colnames(m) <- d$global_cell

  ref_list[[lk]] <- m

  cat(
    sprintf(
      "REF %d/%d %s cells=%d features=%d\n",
      ii,
      length(disc_libraries),
      lk,
      ncol(m),
      nrow(m)
    )
  )

  rm(obj, counts0, m)
  invisible(gc())
}

ref_common <- Reduce(
  intersect,
  lapply(
    ref_list,
    rownames
  )
)

first_features <- rownames(
  ref_list[[1]]
)

ref_common <- first_features[
  first_features %in%
    ref_common
]

if(length(ref_common) < 10000L){
  stop(
    "Too few discovery common features: ",
    length(ref_common)
  )
}

ref_list <- lapply(
  ref_list,
  function(m){
    m[
      ref_common,
      ,
      drop=FALSE
    ]
  }
)

ref_counts <- do.call(
  cbind,
  ref_list
)

rm(ref_list)
invisible(gc())

ri <- match(
  disc$global_cell,
  colnames(ref_counts)
)

stopifnot(
  !anyNA(ri)
)

ref_counts <- ref_counts[
  ,
  ri,
  drop=FALSE
]

stopifnot(
  identical(
    colnames(ref_counts),
    disc$global_cell
  )
)

# ============================================================
# Reference-only HVG selection
# ============================================================

ref_hvg <- CreateSeuratObject(
  counts=ref_counts,
  min.cells=0,
  min.features=0,
  project="progenitor_reference_HVG"
)

ref_hvg <- NormalizeData(
  ref_hvg,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

ref_hvg <- FindVariableFeatures(
  ref_hvg,
  selection.method="vst",
  nfeatures=N_HVG,
  verbose=FALSE
)

hvg_raw <- VariableFeatures(
  ref_hvg
)

transfer_features <- hvg_raw[
  !technical_gene(hvg_raw)
]

if(length(transfer_features) < 1500L){
  stop(
    "Too few transfer features: ",
    length(transfer_features)
  )
}

rm(ref_hvg)
invisible(gc())

ref_counts <- ref_counts[
  transfer_features,
  ,
  drop=FALSE
]

# ============================================================
# Build query counts only for reference-selected features
# ============================================================

query_meta$key <- paste(
  query_meta$project_id,
  query_meta$library_key,
  sep="|||"
)

query_keys <- sort(
  unique(
    query_meta$key
  )
)

query_list <- list()
coverage_rows <- list()

for(ii in seq_along(query_keys)){

  key <- query_keys[[ii]]

  d <- query_meta[
    query_meta$key == key,
    ,
    drop=FALSE
  ]

  li <- match(
    key,
    lib$key
  )

  if(is.na(li)){
    stop(
      "Query library lookup failed: ",
      key
    )
  }

  source_rds <- as.character(
    lib$final_rds_resolved[[li]]
  )

  if(!file.exists(source_rds)){
    stop(
      "Missing query source RDS: ",
      source_rds
    )
  }

  obj <- readRDS(source_rds)
  counts0 <- get_rna_counts(obj)

  stopifnot(
    all(
      d$original_cell %in%
        colnames(counts0)
    )
  )

  raw <- counts0[
    ,
    d$original_cell,
    drop=FALSE
  ]

  colnames(raw) <-
    d$global_cell

  aligned <- align_to_features(
    raw,
    transfer_features
  )

  n_present <- attr(
    aligned,
    "n_features_present"
  )

  coverage_rows[[
    length(coverage_rows)+1L
  ]] <- data.frame(
    project_id=d$project_id[[1]],
    library_key=d$library_key[[1]],
    n_cells=nrow(d),
    n_transfer_features=
      length(transfer_features),
    n_features_present=
      n_present,
    feature_fraction_present=
      n_present /
      length(transfer_features),
    stringsAsFactors=FALSE
  )

  if(
    n_present /
      length(transfer_features) <
      0.80
  ){
    stop(
      "Query library has <80% transfer features: ",
      key
    )
  }

  query_list[[
    length(query_list)+1L
  ]] <- aligned

  cat(
    sprintf(
      "QUERY %d/%d %s cells=%d feature_coverage=%.4f\n",
      ii,
      length(query_keys),
      key,
      ncol(aligned),
      n_present /
        length(transfer_features)
    )
  )

  rm(
    obj,
    counts0,
    raw,
    aligned
  )

  invisible(gc())
}

query_feature_coverage <- do.call(
  rbind,
  coverage_rows
)

query_counts <- do.call(
  cbind,
  query_list
)

rm(
  query_list,
  coverage_rows
)

invisible(gc())

qci <- match(
  query_meta$global_cell,
  colnames(query_counts)
)

stopifnot(
  !anyNA(qci)
)

query_counts <- query_counts[
  ,
  qci,
  drop=FALSE
]

stopifnot(
  identical(
    colnames(query_counts),
    query_meta$global_cell
  )
)

# ============================================================
# Full discovery reference -> 813-cell query mapping
# ============================================================

ref <- CreateSeuratObject(
  counts=ref_counts,
  min.cells=0,
  min.features=0,
  project="progenitor_discovery_reference"
)

query <- CreateSeuratObject(
  counts=query_counts,
  min.cells=0,
  min.features=0,
  project="progenitor_external_query"
)

rm(
  ref_counts,
  query_counts
)

invisible(gc())

ref <- NormalizeData(
  ref,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

query <- NormalizeData(
  query,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

ref <- ScaleData(
  ref,
  features=transfer_features,
  verbose=FALSE
)

npcs <- min(
  N_PCS,
  length(transfer_features)-1L,
  ncol(ref)-1L
)

stopifnot(
  npcs >= 10L
)

ref <- RunPCA(
  ref,
  features=transfer_features,
  npcs=npcs,
  verbose=FALSE
)

anchors <- FindTransferAnchors(
  reference=ref,
  query=query,
  normalization.method="LogNormalize",
  reduction="pcaproject",
  reference.reduction="pca",
  features=transfer_features,
  dims=seq_len(npcs),
  k.anchor=5,
  k.filter=50,
  verbose=TRUE
)

ref_labels <- disc$transfer_class

names(ref_labels) <-
  disc$global_cell

pred <- TransferData(
  anchorset=anchors,
  refdata=ref_labels,
  dims=seq_len(npcs),
  k.weight=50,
  verbose=TRUE
)

pred <- as.data.frame(pred)

pi <- match(
  query_meta$global_cell,
  rownames(pred)
)

stopifnot(
  !anyNA(pi)
)

pred <- pred[
  pi,
  ,
  drop=FALSE
]

predicted_class <- as.character(
  pred$predicted.id
)

pli <- match(
  predicted_class,
  class_lookup$transfer_class
)

stopifnot(
  !anyNA(pli)
)

query_result <- data.frame(
  global_cell=
    query_meta$global_cell,
  project_id=
    query_meta$project_id,
  library_key=
    query_meta$library_key,
  predicted_class=
    predicted_class,
  prediction_score=
    as.numeric(
      pred$prediction.score.max
    ),
  predicted_core=
    class_lookup$core_identity[pli],
  predicted_axis=
    class_lookup$lineage_priming_axis[pli],
  stringsAsFactors=FALSE
)

rm(
  ref,
  query,
  anchors,
  pred
)

invisible(gc())

# ============================================================
# Apply frozen operational policy
# ============================================================

threshold_map <- setNames(
  policy$threshold,
  policy$axis
)

policy_map <- setNames(
  policy$policy,
  policy$axis
)

query_result$core_identity <-
  query_result$predicted_core

query_result$lineage_priming_axis <-
  "unresolved"

query_result$axis_policy_reason <-
  "prediction_not_accepted"

for(i in seq_len(nrow(query_result))){

  core <- query_result$core_identity[[i]]
  ax <- query_result$predicted_axis[[i]]
  score <- query_result$prediction_score[[i]]

  if(
    core ==
      "Deferred_non_progenitor_lymphoid_like"
  ){
    query_result$lineage_priming_axis[[i]] <-
      "not_applicable"

    query_result$axis_policy_reason[[i]] <-
      "deferred_core"

    next
  }

  if(ax == "unresolved"){
    query_result$axis_policy_reason[[i]] <-
      "reference_unresolved"
    next
  }

  if(ax == NONTRANSFERABLE_AXIS){

    query_result$axis_policy_reason[[i]] <-
      "primitive_HSC_prediction_candidate_only"

    next
  }

  if(!(ax %in% POSITIVE_AXES)){

    query_result$axis_policy_reason[[i]] <-
      "noncanonical_prediction"

    next
  }

  pol <- policy_map[[ax]]
  th <- threshold_map[[ax]]

  if(
    identical(
      pol,
      "score_threshold"
    )
  ){

    if(
      is.finite(th) &&
      score >= th
    ){
      query_result$lineage_priming_axis[[i]] <-
        ax

      query_result$axis_policy_reason[[i]] <-
        paste0(
          "accepted_score_threshold_",
          sprintf("%.2f", th)
        )
    } else {
      query_result$axis_policy_reason[[i]] <-
        paste0(
          "below_score_threshold_",
          sprintf("%.2f", th)
        )
    }

    next
  }

  if(
    ax == "granulocytic_priming" &&
    identical(
      pol,
      "native_program_gate_fallback"
    )
  ){

    if(
      isTRUE(
        as.logical(
          qprog$granulocytic_flag[[i]]
        )
      )
    ){
      query_result$lineage_priming_axis[[i]] <-
        ax

      query_result$axis_policy_reason[[i]] <-
        "accepted_native_granulocytic_gate"
    } else {
      query_result$axis_policy_reason[[i]] <-
        "granulocytic_gate_negative"
    }

    next
  }

  query_result$axis_policy_reason[[i]] <-
    "axis_not_directly_transferable"
}

# ============================================================
# State: use native query cell-level program, not transferred
# classifier state.
# ============================================================

cycling_flag <- as.logical(
  qprog$cycling_flag
)

query_result$state <- ifelse(
  query_result$core_identity ==
    "Deferred_non_progenitor_lymphoid_like",
  "not_applicable",
  ifelse(
    cycling_flag,
    "cycling",
    "none"
  )
)

# ============================================================
# Evidence flags
# ============================================================

query_result$evidence_flag <- ""

for(i in seq_len(nrow(query_result))){

  ev <- character(0)

  if(
    query_result$predicted_axis[[i]] ==
      NONTRANSFERABLE_AXIS
  ){
    ev <- c(
      ev,
      "primitive_HSC_like_prediction_not_hard_transferred"
    )
  }

  if(
    "primitive_HSC_low_commitment_flag" %in%
      names(qprog) &&
    isTRUE(
      as.logical(
        qprog$primitive_HSC_low_commitment_flag[[i]]
      )
    )
  ){
    ev <- c(
      ev,
      "native_primitive_HSC_low_commitment_candidate"
    )
  }

  if(
    "B_commitment_flag" %in%
      names(qprog) &&
    isTRUE(
      as.logical(
        qprog$B_commitment_flag[[i]]
      )
    )
  ){
    ev <- c(
      ev,
      "B_commitment_candidate_not_canonical"
    )
  }

  if(
    isTRUE(
      cycling_flag[[i]]
    )
  ){
    ev <- c(
      ev,
      "native_cycling_program"
    )
  }

  ev <- c(
    ev,
    query_result$axis_policy_reason[[i]]
  )

  query_result$evidence_flag[[i]] <-
    paste(
      unique(ev),
      collapse=";"
    )
}

# ============================================================
# Confidence
# ============================================================

query_result$core_transfer_confidence <- ifelse(
  query_result$prediction_score >= 0.80,
  "high",
  ifelse(
    query_result$prediction_score >= 0.50,
    "medium",
    "low"
  )
)

query_result$axis_transfer_confidence <- ifelse(
  query_result$lineage_priming_axis %in%
    c(
      "unresolved",
      "not_applicable"
    ),
  ifelse(
    query_result$lineage_priming_axis ==
      "not_applicable",
    query_result$core_transfer_confidence,
    "unresolved"
  ),
  ifelse(
    query_result$prediction_score >= 0.90,
    "high",
    "medium"
  )
)

query_result$annotation_confidence <- ifelse(
  query_result$core_identity ==
    "Deferred_non_progenitor_lymphoid_like",
  query_result$core_transfer_confidence,
  ifelse(
    query_result$lineage_priming_axis ==
      "unresolved",
    "low_medium",
    query_result$axis_transfer_confidence
  )
)

query_result$annotation_source <-
  "full_discovery_reference_transfer_with_native_policy"

# ============================================================
# Standardize discovery + query output
# ============================================================

disc_final <- data.frame(
  global_cell=disc$global_cell,
  project_id=disc$project_id,
  library_key=disc$library_key,
  annotation_origin=
    "discovery_reference",
  discovery_cluster=
    disc$discovery_cluster,
  core_identity=
    disc$core_identity,
  lineage_priming_axis=
    disc$lineage_priming_axis,
  state=
    disc$state,
  annotation_confidence=
    disc$annotation_confidence,
  core_transfer_confidence=
    "reference",
  axis_transfer_confidence=
    "reference",
  prediction_score=
    NA_real_,
  evidence_flag=
    disc$evidence_flag,
  annotation_source=
    disc$annotation_source,
  stringsAsFactors=FALSE
)

query_final <- data.frame(
  global_cell=query_result$global_cell,
  project_id=query_result$project_id,
  library_key=query_result$library_key,
  annotation_origin=
    "external_query",
  discovery_cluster=
    NA_character_,
  core_identity=
    query_result$core_identity,
  lineage_priming_axis=
    query_result$lineage_priming_axis,
  state=
    query_result$state,
  annotation_confidence=
    query_result$annotation_confidence,
  core_transfer_confidence=
    query_result$core_transfer_confidence,
  axis_transfer_confidence=
    query_result$axis_transfer_confidence,
  prediction_score=
    query_result$prediction_score,
  evidence_flag=
    query_result$evidence_flag,
  annotation_source=
    query_result$annotation_source,
  stringsAsFactors=FALSE
)

full <- rbind(
  disc_final,
  query_final
)

stopifnot(
  nrow(full) == N_FULL,
  !anyDuplicated(full$global_cell),
  sum(
    full$annotation_origin ==
      "discovery_reference"
  ) == N_DISCOVERY,
  sum(
    full$annotation_origin ==
      "external_query"
  ) == N_QUERY
)

# ============================================================
# Summaries
# ============================================================

core_summary <- as.data.frame(
  table(
    core_identity=
      full$core_identity
  ),
  stringsAsFactors=FALSE
)

axis_summary <- as.data.frame(
  table(
    lineage_priming_axis=
      full$lineage_priming_axis
  ),
  stringsAsFactors=FALSE
)

state_summary <- as.data.frame(
  table(
    state=
      full$state
  ),
  stringsAsFactors=FALSE
)

origin_summary <- as.data.frame(
  table(
    annotation_origin=
      full$annotation_origin
  ),
  stringsAsFactors=FALSE
)

query_axis_summary <- as.data.frame(
  table(
    predicted_axis=
      query_result$predicted_axis,
    final_axis=
      query_result$lineage_priming_axis
  ),
  stringsAsFactors=FALSE
)

query_axis_summary <- query_axis_summary[
  query_axis_summary$Freq > 0,
  ,
  drop=FALSE
]

query_project_summary <- as.data.frame(
  table(
    project_id=
      query_final$project_id,
    core_identity=
      query_final$core_identity,
    lineage_priming_axis=
      query_final$lineage_priming_axis
  ),
  stringsAsFactors=FALSE
)

query_project_summary <- query_project_summary[
  query_project_summary$Freq > 0,
  ,
  drop=FALSE
]

# ============================================================
# Write
# ============================================================

policy_file <- file.path(
  out_dir,
  "progenitor_full_transfer_policy_v1.tsv"
)

write.table(
  policy,
  policy_file,
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  threshold_audit,
  file.path(
    out_dir,
    "progenitor_full_transfer_threshold_calibration_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  query_feature_coverage,
  file.path(
    out_dir,
    "progenitor_query_transfer_feature_coverage_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

query_file <- file.path(
  out_dir,
  paste0(
    "progenitor_query_cell_annotation_v1__",
    tag,
    ".tsv.gz"
  )
)

write_gz_tsv(
  query_final,
  query_file
)

full_file <- file.path(
  out_dir,
  paste0(
    "progenitor_full_cell_annotation_v1__",
    tag,
    ".tsv.gz"
  )
)

write_gz_tsv(
  full,
  full_file
)

write.table(
  core_summary,
  file.path(
    out_dir,
    "progenitor_full_core_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  axis_summary,
  file.path(
    out_dir,
    "progenitor_full_axis_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  state_summary,
  file.path(
    out_dir,
    "progenitor_full_state_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  origin_summary,
  file.path(
    out_dir,
    "progenitor_full_origin_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  query_axis_summary,
  file.path(
    out_dir,
    "progenitor_query_predicted_to_final_axis_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  query_project_summary,
  file.path(
    out_dir,
    "progenitor_query_project_annotation_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Provenance
# ============================================================

source_manifest <- data.frame(
  stage=c(
    "Step41_discovery_freeze",
    "Step42_leave_one_library_out",
    "Step38_query_support_definition",
    "Step40_cross_project_native_programs"
  ),
  path=c(
    discovery_dir,
    loo_dir,
    cluster_dir,
    replication_dir
  ),
  stringsAsFactors=FALSE
)

write.table(
  source_manifest,
  file.path(
    out_dir,
    "progenitor_full_annotation_source_manifest_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Canonical freeze RDS
# ============================================================

freeze <- list(
  annotation_version=
    "Progenitor_full_annotation_v1",

  n_full_cells=N_FULL,
  n_discovery_cells=N_DISCOVERY,
  n_query_cells=N_QUERY,

  discovery_project=
    DISCOVERY_PROJECT,

  reference_mapping=
    "Seurat pcaproject FindTransferAnchors TransferData",

  transfer_features=
    transfer_features,

  transfer_policy=
    policy,

  discovery_annotation=
    disc_final,

  query_annotation=
    query_final,

  full_annotation=
    full,

  core_summary=
    core_summary,

  axis_summary=
    axis_summary,

  state_summary=
    state_summary,

  query_feature_coverage=
    query_feature_coverage,

  source_manifest=
    source_manifest
)

rds_file <- file.path(
  out_dir,
  paste0(
    "progenitor_full_annotation_freeze_v1__",
    tag,
    ".rds"
  )
)

tmp_rds <- paste0(
  rds_file,
  ".tmp"
)

saveRDS(
  freeze,
  tmp_rds,
  compress=TRUE
)

if(!file.rename(
  tmp_rds,
  rds_file
)){
  stop("Atomic RDS rename failed")
}

reread <- readRDS(
  rds_file
)

stopifnot(
  reread$n_full_cells == N_FULL,
  identical(
    reread$full_annotation$global_cell,
    full$global_cell
  )
)

# ============================================================
# Re-read gzip outputs
# ============================================================

query_check <- read_tsv(
  query_file
)

full_check <- read_tsv(
  full_file
)

stopifnot(
  nrow(query_check) == N_QUERY,
  identical(
    as.character(query_check$global_cell),
    as.character(query_final$global_cell)
  ),
  nrow(full_check) == N_FULL,
  identical(
    as.character(full_check$global_cell),
    as.character(full$global_cell)
  ),
  identical(
    as.character(full_check$core_identity),
    as.character(full$core_identity)
  ),
  identical(
    as.character(
      full_check$lineage_priming_axis
    ),
    as.character(
      full$lineage_priming_axis
    )
  ),
  identical(
    as.character(full_check$state),
    as.character(full$state)
  )
)

# ============================================================
# SHA256
# ============================================================

sha_files <- c(
  policy_file,
  file.path(
    out_dir,
    "progenitor_full_transfer_threshold_calibration_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_query_transfer_feature_coverage_v1.tsv"
  ),
  query_file,
  full_file,
  file.path(
    out_dir,
    "progenitor_full_core_summary_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_full_axis_summary_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_full_state_summary_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_full_origin_summary_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_query_predicted_to_final_axis_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_query_project_annotation_summary_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_full_annotation_source_manifest_v1.tsv"
  ),
  rds_file
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
    "Progenitor full annotation freeze v1",
    "annotation_version=Progenitor_full_annotation_v1",
    "n_full_cells=27775",
    "n_discovery_cells=26962",
    "n_query_cells=813",
    "discovery_project=GSE216007",
    "",
    "discovery annotation=frozen Step41",
    "query mapping=Seurat reference pcaproject",
    "query transfer features selected from discovery reference only",
    "query libraries not used for HVG selection or PCA",
    "",
    "LOO calibration target precision=0.90",
    "minimum calibration calls=100",
    "primitive_HSC_like not directly hard-transferred to query",
    "primitive HSC predictions retained as candidate evidence",
    "granulocytic native-program fallback permitted if required",
    "B commitment retained as evidence only, not canonical axis",
    "query cycling state assigned from native cell-level program",
    "",
    "condition not used",
    "native SoupX RNA",
    "no Harmony",
    "no CCA",
    "no joint cross-project integration",
    "cellranger_count not accessed",
    "",
    "gzip_integrity=PASS",
    "query_reread_count_order=PASS",
    "full_reread_count_order_annotation=PASS",
    "freeze_rds_reread=PASS",
    "sha256_manifest=PASS"
  ),
  done_file
)

cat("\n===== TRANSFER POLICY =====\n")
print(policy, row.names=FALSE)

cat("\n===== CORE =====\n")
print(core_summary, row.names=FALSE)

cat("\n===== AXIS =====\n")
print(axis_summary, row.names=FALSE)

cat("\n===== STATE =====\n")
print(state_summary, row.names=FALSE)

cat("\n===== QUERY PREDICTED -> FINAL AXIS =====\n")
print(query_axis_summary, row.names=FALSE)

cat("\n===== QUERY FEATURE COVERAGE =====\n")
print(query_feature_coverage, row.names=FALSE)

cat("\nOUT_DIR=", out_dir, "\n", sep="")
cat("\nPASS: Progenitor full annotation freeze completed\n")
