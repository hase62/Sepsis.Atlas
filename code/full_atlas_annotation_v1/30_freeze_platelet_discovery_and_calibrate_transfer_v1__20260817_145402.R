#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 5L){
  stop(
    "Usage: script <root> <tag> <discovery_dir> <marker_dir> <screen_dir>"
  )
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]
discovery_dir <- normalizePath(args[[3]], mustWork=TRUE)
marker_dir <- normalizePath(args[[4]], mustWork=TRUE)
screen_dir <- normalizePath(args[[5]], mustWork=TRUE)

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

set.seed(20260817)

N_DISCOVERY <- 26668L
N_FULL <- 49585L
N_PCS <- 20L

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "platelet_discovery_freeze_transfer_calibration_v1__",
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
  "PLATELET_DISCOVERY_FREEZE_TRANSFER_CALIBRATION_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop("Already completed: ", done_file)
}

read_tsv <- function(path){

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

  con <- gzfile(
    path,
    "wt"
  )

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

to_logical <- function(x){

  if(is.logical(x)){
    return(x)
  }

  if(is.numeric(x)){
    return(x != 0)
  }

  tolower(
    as.character(x)
  ) %in% c(
    "true","t","1","yes"
  )
}

# ============================================================
# Inputs
# ============================================================

assignment_file <- file.path(
  discovery_dir,
  "platelet_balanced_discovery_cluster_assignments_v1.tsv.gz"
)

discovery_rds <- file.path(
  discovery_dir,
  "platelet_balanced_discovery_clustering_v1.rds"
)

marker_file <- file.path(
  marker_dir,
  "platelet_r0p4_project_balanced_marker_stats_v1.tsv.gz"
)

screen_file <- file.path(
  screen_dir,
  "platelet_r0p4_corrected_consensus_marker_screen_v1.tsv"
)

transfer_dir <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "annotation_transfer_by_library"
)

lib_file <- file.path(
  root,
  "pre_integration",
  "pilot_unintegrated_large_v1",
  "resolved_library_table.tsv"
)

for(f in c(
  assignment_file,
  discovery_rds,
  marker_file,
  screen_file,
  lib_file
)){
  if(!file.exists(f)){
    stop("Missing input: ", f)
  }
}

a <- read_tsv(
  assignment_file
)

screen <- read_tsv(
  screen_file
)

stopifnot(
  nrow(a) == N_DISCOVERY,
  !anyDuplicated(a$global_cell)
)

a$global_cell <- as.character(a$global_cell)
a$project_id <- as.character(a$project_id)
a$library_key <- as.character(a$library_key)
a$cluster_r0p4 <- as.character(a$cluster_r0p4)

screen$cluster <- as.character(
  screen$cluster
)

# ============================================================
# Verify expected 29b evidence before freezing decision
# ============================================================

get_screen <- function(k, column){

  i <- match(
    as.character(k),
    screen$cluster
  )

  if(is.na(i)){
    stop("Missing screen cluster: ", k)
  }

  screen[[column]][[i]]
}

stopifnot(

  get_screen(
    "5",
    "corrected_scope_screen"
  ) ==
    "Deferred_non_platelet_myeloid_like",

  get_screen(
    "9",
    "corrected_scope_screen"
  ) ==
    "Deferred_non_platelet_T_NK_like",

  get_screen(
    "2",
    "megakaryocytic_transcription_screen"
  ) ==
    "megakaryocytic_transcription_high",

  get_screen(
    "6",
    "megakaryocytic_transcription_screen"
  ) ==
    "megakaryocytic_transcription_high",

  get_screen(
    "8",
    "megakaryocytic_transcription_screen"
  ) ==
    "megakaryocytic_transcription_candidate",

  get_screen(
    "3",
    "corrected_state_screen"
  ) ==
    "activation_degranulation_candidate",

  get_screen(
    "4",
    "corrected_state_screen"
  ) ==
    "activation_degranulation_candidate",

  as.integer(
    get_screen(
      "7",
      "n_interferon_specific_hits"
    )
  ) == 3L
)

# ============================================================
# Frozen discovery decision
#
# Three independent dimensions:
#   core
#   megakaryocytic transcriptional axis
#   orthogonal state
# ============================================================

decision <- data.frame(
  cluster=as.character(0:9),

  platelet_core_v1=c(
    "Platelet_like",                     # 0
    "Platelet_like",                     # 1
    "Platelet_like",                     # 2
    "Platelet_like",                     # 3
    "Platelet_like",                     # 4
    "Deferred_non_platelet_myeloid_like",# 5
    "Platelet_like",                     # 6
    "Platelet_like",                     # 7
    "Platelet_like",                     # 8
    "Deferred_non_platelet_T_NK_like"    # 9
  ),

  platelet_transcriptional_axis_v1=c(
    "none",
    "none",
    "megakaryocytic_transcription_high",
    "none",
    "none",
    "not_applicable",
    "megakaryocytic_transcription_high",
    "none",
    "megakaryocytic_transcription_candidate",
    "not_applicable"
  ),

  platelet_state_v1=c(
    "none",
    "none",
    "none",
    "activation_degranulation_candidate",
    "activation_degranulation_candidate",
    "not_applicable",
    "none",
    "interferon_stimulated_candidate",
    "none",
    "not_applicable"
  ),

  core_evidence_confidence_v1=c(
    "high","high","high","high","high",
    "high",
    "high","high","high",
    "high"
  ),

  transcriptional_axis_evidence_confidence_v1=c(
    "not_applicable",
    "not_applicable",
    "high",
    "not_applicable",
    "not_applicable",
    "not_applicable",
    "high",
    "not_applicable",
    "medium",
    "not_applicable"
  ),

  state_evidence_confidence_v1=c(
    "not_applicable",
    "not_applicable",
    "not_applicable",
    "medium",
    "medium",
    "not_applicable",
    "not_applicable",
    "low_medium",
    "not_applicable",
    "not_applicable"
  ),

  decision_basis_v1=c(
    "baseline platelet cluster; no specific contaminant consensus program",
    "baseline platelet cluster; no specific contaminant consensus program",
    "8 formal megakaryocytic-transcription consensus markers",
    "platelet identity plus 4 specific activation/degranulation consensus markers",
    "platelet identity plus 4 specific activation/degranulation consensus markers",
    "11 specific project-balanced myeloid consensus markers",
    "platelet identity plus 13 megakaryocytic-transcription consensus markers",
    paste0(
      "platelet identity plus IFI6/IFI27/LY6E formal IFN markers; ",
      "IFITM1/IFITM3 additionally supported in native top markers"
    ),
    "4 formal megakaryocytic-transcription consensus markers",
    "7 specific project-balanced T/NK consensus markers"
  ),

  stringsAsFactors=FALSE
)

decision$transfer_class_v1 <- paste(
  decision$platelet_core_v1,
  decision$platelet_transcriptional_axis_v1,
  decision$platelet_state_v1,
  sep="|||"
)

stopifnot(
  !anyDuplicated(decision$cluster),
  setequal(
    decision$cluster,
    unique(a$cluster_r0p4)
  )
)

# ============================================================
# Discovery-cell annotation
# ============================================================

di <- match(
  a$cluster_r0p4,
  decision$cluster
)

stopifnot(
  !anyNA(di)
)

ann <- data.frame(
  global_cell=a$global_cell,
  project_id=a$project_id,
  library_key=a$library_key,
  condition_binary=
    as.character(
      a$condition_binary
    ),
  cluster_r0p4=
    a$cluster_r0p4,

  platelet_core_v1=
    decision$platelet_core_v1[di],

  platelet_transcriptional_axis_v1=
    decision$platelet_transcriptional_axis_v1[di],

  platelet_state_v1=
    decision$platelet_state_v1[di],

  platelet_core_evidence_confidence_v1=
    decision$core_evidence_confidence_v1[di],

  platelet_transcriptional_axis_evidence_confidence_v1=
    decision$transcriptional_axis_evidence_confidence_v1[di],

  platelet_state_evidence_confidence_v1=
    decision$state_evidence_confidence_v1[di],

  platelet_transfer_class_v1=
    decision$transfer_class_v1[di],

  stringsAsFactors=FALSE
)

# ============================================================
# Write discovery freeze BEFORE calibration
# ============================================================

decision_file <- file.path(
  out_dir,
  "platelet_r0p4_annotation_decision_v1.tsv"
)

write.table(
  decision,
  decision_file,
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

ann_file <- file.path(
  out_dir,
  "platelet_discovery_cell_annotation_v1.tsv.gz"
)

write_gz_tsv(
  ann,
  ann_file
)

chk <- read_tsv(
  ann_file
)

stopifnot(
  nrow(chk) == N_DISCOVERY,
  identical(
    as.character(chk$global_cell),
    ann$global_cell
  )
)

rm(chk)

# ============================================================
# Transfer features
# ============================================================

disc <- readRDS(
  discovery_rds
)

stopifnot(
  identical(
    disc$analysis_version,
    "Platelet_balanced_discovery_clustering_v1"
  )
)

integration_features <- as.character(
  disc$integration_features
)

cons <- read_tsv(
  marker_file
)

cons$marker_candidate <- to_logical(
  cons$marker_candidate
)

formal_marker_genes <- unique(
  as.character(
    cons$gene[
      cons$marker_candidate
    ]
  )
)

feature_candidates <- unique(
  c(
    integration_features,
    formal_marker_genes
  )
)

cat(
  "candidate transfer features=",
  length(feature_candidates),
  "\n",
  sep=""
)

# ============================================================
# Map pilot global cells back to native source cells
# ============================================================

transfer_files <- list.files(
  transfer_dir,
  pattern="\\.tsv\\.gz$",
  full.names=TRUE
)

meta_list <- list()

for(f in transfer_files){

  z <- read_tsv(f)

  required <- c(
    "global_cell",
    "original_cell",
    "project_id",
    "library_key",
    "integration_compartment_primary_v1"
  )

  stopifnot(
    all(required %in% names(z))
  )

  z <- z[
    as.character(
      z$integration_compartment_primary_v1
    ) ==
      "Platelet_megakaryocyte",
    required,
    drop=FALSE
  ]

  if(nrow(z)){
    meta_list[[length(meta_list)+1L]] <- z
  }
}

full_meta <- do.call(
  rbind,
  meta_list
)

full_meta$global_cell <-
  as.character(
    full_meta$global_cell
  )

full_meta$original_cell <-
  as.character(
    full_meta$original_cell
  )

full_meta$project_id <-
  as.character(
    full_meta$project_id
  )

full_meta$library_key <-
  as.character(
    full_meta$library_key
  )

stopifnot(
  nrow(full_meta) == N_FULL,
  !anyDuplicated(
    full_meta$global_cell
  )
)

mi <- match(
  ann$global_cell,
  full_meta$global_cell
)

stopifnot(
  !anyNA(mi)
)

ann$original_cell <-
  full_meta$original_cell[mi]

stopifnot(
  ann$project_id ==
    full_meta$project_id[mi],
  ann$library_key ==
    full_meta$library_key[mi]
)

# ============================================================
# Library lookup
# ============================================================

lib <- read_tsv(
  lib_file
)

lib$project_id <-
  as.character(
    lib$project_id
  )

lib$library_key <-
  as.character(
    lib$library_key
  )

lib$key <- paste(
  lib$project_id,
  lib$library_key,
  sep="|||"
)

# ============================================================
# Build native discovery count matrix using ONLY transfer genes
# ============================================================

library_keys <- unique(
  paste(
    ann$project_id,
    ann$library_key,
    sep="|||"
  )
)

matrices <- list()

transfer_features <- NULL

for(ii in seq_along(library_keys)){

  key <- library_keys[[ii]]

  sp <- strsplit(
    key,
    "\\|\\|\\|"
  )[[1]]

  p <- sp[[1]]
  lk <- sp[[2]]

  d <- ann[
    ann$project_id == p &
      ann$library_key == lk,
    ,
    drop=FALSE
  ]

  li <- match(
    key,
    lib$key
  )

  if(is.na(li)){
    stop(
      "Missing library lookup: ",
      key
    )
  }

  source_rds <- as.character(
    lib$final_rds_resolved[[li]]
  )

  if(!file.exists(source_rds)){
    stop(
      "Missing source RDS: ",
      source_rds
    )
  }

  obj0 <- readRDS(
    source_rds
  )

  counts0 <- get_rna_counts(
    obj0
  )

  if(is.null(transfer_features)){

    transfer_features <- feature_candidates[
      feature_candidates %in%
        rownames(counts0)
    ]

    if(length(transfer_features) < 1200L){
      stop(
        "Too few transfer features: ",
        length(transfer_features)
      )
    }

  } else {

    if(!all(
      transfer_features %in%
        rownames(counts0)
    )){
      stop(
        "Transfer feature mismatch: ",
        key
      )
    }
  }

  cells <- d$original_cell

  if(!all(
    cells %in%
      colnames(counts0)
  )){
    stop(
      "Missing pilot cells: ",
      key
    )
  }

  m <- counts0[
    transfer_features,
    cells,
    drop=FALSE
  ]

  colnames(m) <-
    d$global_cell

  matrices[[length(matrices)+1L]] <- m

  rm(
    obj0,
    counts0,
    m
  )

  invisible(gc())

  cat(
    sprintf(
      "READ %3d/%3d %s\n",
      ii,
      length(library_keys),
      key
    )
  )
}

counts <- do.call(
  cbind,
  matrices
)

rm(matrices)

counts <- counts[
  ,
  ann$global_cell,
  drop=FALSE
]

stopifnot(
  ncol(counts) == N_DISCOVERY,
  identical(
    colnames(counts),
    ann$global_cell
  )
)

# ============================================================
# Native Seurat object
# ============================================================

md <- ann
rownames(md) <- md$global_cell

obj <- CreateSeuratObject(
  counts=counts,
  meta.data=md,
  min.cells=0,
  min.features=0
)

rm(counts)
invisible(gc())

obj <- NormalizeData(
  obj,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

# ============================================================
# Deterministic stratified 80/20 calibration split
#
# Stratification = project x transfer class.
# Small strata remain in training.
# ============================================================

stratum <- paste(
  obj$project_id,
  obj$platelet_transfer_class_v1,
  sep="|||"
)

groups <- split(
  colnames(obj),
  stratum
)

validation_cells <- character()

set.seed(20260817)

for(g in groups){

  n <- length(g)

  if(n < 6L){
    next
  }

  n_val <- max(
    1L,
    floor(
      0.20 * n
    )
  )

  n_val <- min(
    n_val,
    n - 4L
  )

  validation_cells <- c(
    validation_cells,
    sample(
      g,
      n_val,
      replace=FALSE
    )
  )
}

validation_cells <- unique(
  validation_cells
)

training_cells <- setdiff(
  colnames(obj),
  validation_cells
)

stopifnot(
  length(training_cells) +
    length(validation_cells) ==
      N_DISCOVERY,
  length(validation_cells) > 3000L
)

ref <- subset(
  obj,
  cells=training_cells
)

query <- subset(
  obj,
  cells=validation_cells
)

rm(obj)
invisible(gc())

# ============================================================
# Reference PCA on native RNA
# ============================================================

pca_features <- transfer_features[
  transfer_features %in%
    rownames(ref)
]

ref <- ScaleData(
  ref,
  features=pca_features,
  verbose=FALSE
)

ref <- RunPCA(
  ref,
  features=pca_features,
  npcs=N_PCS,
  verbose=FALSE
)

# Use only features actually retained by PCA.
pca_features <- rownames(
  Loadings(
    ref[["pca"]]
  )
)

stopifnot(
  length(pca_features) >= 1000L
)

# ============================================================
# Transfer calibration
# ============================================================

anchors <- FindTransferAnchors(
  reference=ref,
  query=query,
  normalization.method="LogNormalize",
  reduction="pcaproject",
  reference.reduction="pca",
  features=pca_features,
  dims=seq_len(N_PCS),
  k.anchor=5,
  k.filter=50,
  verbose=TRUE
)

ref_labels <- as.character(
  ref$platelet_transfer_class_v1
)

names(ref_labels) <- colnames(ref)

pred <- TransferData(
  anchorset=anchors,
  refdata=ref_labels,
  dims=seq_len(N_PCS),
  k.weight=50,
  verbose=TRUE
)

pred$global_cell <- rownames(pred)

qi <- match(
  pred$global_cell,
  colnames(query)
)

stopifnot(
  !anyNA(qi)
)

pred$actual_class <-
  as.character(
    query$platelet_transfer_class_v1[
      qi
    ]
  )

pred$project_id <-
  as.character(
    query$project_id[
      qi
    ]
  )

pred$exact_correct <-
  pred$predicted.id ==
    pred$actual_class

# ============================================================
# Decode classes into independent dimensions
# ============================================================

class_map <- unique(
  decision[
    ,
    c(
      "transfer_class_v1",
      "platelet_core_v1",
      "platelet_transcriptional_axis_v1",
      "platelet_state_v1"
    ),
    drop=FALSE
  ]
)

actual_i <- match(
  pred$actual_class,
  class_map$transfer_class_v1
)

pred_i <- match(
  pred$predicted.id,
  class_map$transfer_class_v1
)

stopifnot(
  !anyNA(actual_i),
  !anyNA(pred_i)
)

pred$actual_core <-
  class_map$platelet_core_v1[
    actual_i
  ]

pred$predicted_core <-
  class_map$platelet_core_v1[
    pred_i
  ]

pred$actual_axis <-
  class_map$platelet_transcriptional_axis_v1[
    actual_i
  ]

pred$predicted_axis <-
  class_map$platelet_transcriptional_axis_v1[
    pred_i
  ]

pred$actual_state <-
  class_map$platelet_state_v1[
    actual_i
  ]

pred$predicted_state <-
  class_map$platelet_state_v1[
    pred_i
  ]

pred$core_correct <-
  pred$actual_core ==
    pred$predicted_core

pred$axis_correct <-
  pred$actual_axis ==
    pred$predicted_axis

pred$state_correct <-
  pred$actual_state ==
    pred$predicted_state

# ============================================================
# Metrics
# ============================================================

overall <- data.frame(
  metric=c(
    "n_reference",
    "n_validation",
    "n_transfer_features",
    "n_pca_features",
    "exact_class_accuracy",
    "core_accuracy",
    "transcriptional_axis_accuracy",
    "state_accuracy",
    "median_prediction_score"
  ),
  value=c(
    ncol(ref),
    nrow(pred),
    length(transfer_features),
    length(pca_features),
    mean(pred$exact_correct),
    mean(pred$core_correct),
    mean(pred$axis_correct),
    mean(pred$state_correct),
    median(pred$prediction.score.max)
  ),
  stringsAsFactors=FALSE
)

# Confidence-bin accuracy.
pred$score_bin <- cut(
  pred$prediction.score.max,
  breaks=c(
    -Inf,
    0.70,
    0.80,
    Inf
  ),
  labels=c(
    "<0.70",
    "0.70-0.80",
    ">=0.80"
  ),
  right=FALSE
)

score_metrics <- do.call(
  rbind,
  lapply(
    levels(pred$score_bin),
    function(b){

      d <- pred[
        pred$score_bin == b,
        ,
        drop=FALSE
      ]

      data.frame(
        score_bin=b,
        n_cells=nrow(d),
        fraction=
          nrow(d) /
          nrow(pred),
        exact_accuracy=
          if(nrow(d))
            mean(d$exact_correct)
          else
            NA_real_,
        core_accuracy=
          if(nrow(d))
            mean(d$core_correct)
          else
            NA_real_,
        axis_accuracy=
          if(nrow(d))
            mean(d$axis_correct)
          else
            NA_real_,
        state_accuracy=
          if(nrow(d))
            mean(d$state_correct)
          else
            NA_real_,
        stringsAsFactors=FALSE
      )
    }
  )
)

# Per-class recall.
class_metrics <- do.call(
  rbind,
  lapply(
    sort(
      unique(
        pred$actual_class
      )
    ),
    function(cl){

      d <- pred[
        pred$actual_class == cl,
        ,
        drop=FALSE
      ]

      data.frame(
        transfer_class_v1=cl,
        n_validation=nrow(d),
        recall_exact=
          mean(
            d$predicted.id == cl
          ),
        median_prediction_score=
          median(
            d$prediction.score.max
          ),
        stringsAsFactors=FALSE
      )
    }
  )
)

confusion <- as.data.frame.matrix(
  table(
    actual=pred$actual_class,
    predicted=pred$predicted.id
  )
)

confusion$actual_class <-
  rownames(confusion)

confusion <- confusion[
  ,
  c(
    "actual_class",
    setdiff(
      names(confusion),
      "actual_class"
    )
  ),
  drop=FALSE
]

# Project-level technical validation.
project_metrics <- do.call(
  rbind,
  lapply(
    sort(
      unique(
        pred$project_id
      )
    ),
    function(p){

      d <- pred[
        pred$project_id == p,
        ,
        drop=FALSE
      ]

      data.frame(
        project_id=p,
        n_validation=nrow(d),
        exact_accuracy=
          mean(d$exact_correct),
        core_accuracy=
          mean(d$core_correct),
        axis_accuracy=
          mean(d$axis_correct),
        state_accuracy=
          mean(d$state_correct),
        median_prediction_score=
          median(
            d$prediction.score.max
          ),
        stringsAsFactors=FALSE
      )
    }
  )
)

# ============================================================
# Save
# ============================================================

write.table(
  decision,
  file=file.path(
    out_dir,
    "platelet_r0p4_annotation_decision_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  overall,
  file=file.path(
    out_dir,
    "platelet_transfer_calibration_overall_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  score_metrics,
  file=file.path(
    out_dir,
    "platelet_transfer_calibration_by_score_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  class_metrics,
  file=file.path(
    out_dir,
    "platelet_transfer_calibration_by_class_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  project_metrics,
  file=file.path(
    out_dir,
    "platelet_transfer_calibration_by_project_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  confusion,
  file=file.path(
    out_dir,
    "platelet_transfer_calibration_confusion_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write_gz_tsv(
  pred,
  file.path(
    out_dir,
    "platelet_transfer_calibration_predictions_v1.tsv.gz"
  )
)

saveRDS(
  list(
    annotation_version=
      "Platelet_discovery_annotation_v1",
    n_discovery_cells=
      N_DISCOVERY,
    decision=
      decision,
    discovery_annotation=
      ann,
    transfer_features=
      transfer_features,
    pca_features=
      pca_features,
    training_cells=
      training_cells,
    validation_cells=
      validation_cells,
    calibration_overall=
      overall,
    calibration_by_score=
      score_metrics,
    calibration_by_class=
      class_metrics,
    calibration_by_project=
      project_metrics
  ),
  file.path(
    out_dir,
    "platelet_discovery_freeze_transfer_calibration_v1.rds"
  ),
  compress=TRUE
)

# ============================================================
# Completion
# ============================================================

writeLines(
  c(
    "PASS",
    "Platelet discovery freeze + transfer calibration v1",
    "discovery_cells=26668",
    "full_cells=49585",
    "backbone_resolution=r0.4",
    "",
    "core / megakaryocytic transcription / state kept separate",
    "c5 deferred myeloid-like",
    "c9 deferred T/NK-like",
    "c7 IFN candidate retained with low_medium biological evidence",
    "",
    "transfer calibration uses native SoupX-corrected RNA",
    "reference-query mapping=Seurat pcaproject",
    "split=project x transfer-class stratified held-out validation",
    "no full-query transfer performed yet",
    "cellranger_count not accessed",
    "gzip_integrity=PASS",
    "prediction re-read ready"
  ),
  done_file
)

cat(
  "\n===== DECISION =====\n"
)

print(
  decision,
  row.names=FALSE
)

cat(
  "\n===== CALIBRATION OVERALL =====\n"
)

print(
  overall,
  row.names=FALSE
)

cat(
  "\n===== CALIBRATION BY SCORE =====\n"
)

print(
  score_metrics,
  row.names=FALSE
)

cat(
  "\n===== CALIBRATION BY CLASS =====\n"
)

print(
  class_metrics,
  row.names=FALSE
)

cat(
  "\n===== CALIBRATION BY PROJECT =====\n"
)

print(
  project_metrics,
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat("\nPASS\n")
