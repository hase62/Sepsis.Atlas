#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3L){
  stop("Usage: script <root> <tag> <freeze_dir>")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]
freeze_dir <- normalizePath(args[[3]], mustWork=TRUE)

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

N_DISCOVERY <- 8692L
N_FULL <- 12177L
VALIDATION_FRACTION <- 0.20
N_HVG_PER_PROJECT <- 2500L
N_INTEGRATION_FEATURES <- 4000L
MAX_PCS <- 20L

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "erythroid_discovery_transfer_calibration_v1__",
    tag
  )
)

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

done_file <- file.path(
  out_dir,
  "ERYTHROID_DISCOVERY_TRANSFER_CALIBRATION_COMPLETE.ok"
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

# ============================================================
# Load frozen discovery annotation
# ============================================================

freeze_files <- list.files(
  freeze_dir,
  pattern="^erythroid_discovery_annotation_freeze_v1__.*\\.rds$",
  full.names=TRUE
)

if(length(freeze_files) != 1L){
  stop(
    "Expected exactly one discovery freeze RDS; found ",
    length(freeze_files)
  )
}

freeze <- readRDS(freeze_files[[1]])

ann <- freeze$annotation

stopifnot(
  nrow(ann) == N_DISCOVERY,
  !anyDuplicated(ann$global_cell)
)

required_ann <- c(
  "global_cell",
  "project_id",
  "library_key",
  "cluster_r0p4",
  "erythroid_core_v1",
  "erythroid_maturation_axis_v1",
  "erythroid_state_v1",
  "erythroid_transfer_class_v1"
)

stopifnot(
  all(required_ann %in% names(ann))
)

for(nm in required_ann){
  ann[[nm]] <- as.character(ann[[nm]])
}

# ============================================================
# Stratified train/validation split
#
# Stratum = project x frozen transfer class
# No annotation information from validation is used for training.
# ============================================================

ann$stratum <- paste(
  ann$project_id,
  ann$erythroid_transfer_class_v1,
  sep="|||STRATUM|||"
)

strata <- split(
  seq_len(nrow(ann)),
  ann$stratum
)

validation_idx <- integer(0)

for(s in sort(names(strata))){

  idx <- strata[[s]]
  n <- length(idx)

  n_val <- if(n >= 10L){
    max(1L, floor(n * VALIDATION_FRACTION))
  } else if(n >= 5L){
    1L
  } else {
    0L
  }

  # Always retain at least 3 reference cells in a stratum.
  n_val <- min(
    n_val,
    max(0L, n - 3L)
  )

  if(n_val > 0L){
    validation_idx <- c(
      validation_idx,
      sample(idx, n_val)
    )
  }
}

validation_idx <- sort(unique(validation_idx))
reference_idx <- setdiff(seq_len(nrow(ann)), validation_idx)

ref_ann <- ann[reference_idx, , drop=FALSE]
val_ann <- ann[validation_idx, , drop=FALSE]

stopifnot(
  nrow(ref_ann) + nrow(val_ann) == N_DISCOVERY,
  !length(intersect(ref_ann$global_cell, val_ann$global_cell)),
  nrow(val_ann) > 500L
)

cat(
  "Reference train cells = ",
  nrow(ref_ann),
  "\n",
  sep=""
)

cat(
  "Validation cells = ",
  nrow(val_ann),
  "\n",
  sep=""
)

# ============================================================
# Recover original cell IDs
# ============================================================

transfer_dir <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "annotation_transfer_by_library"
)

transfer_files <- list.files(
  transfer_dir,
  pattern="\\.tsv\\.gz$",
  full.names=TRUE
)

if(!length(transfer_files)){
  stop("No annotation transfer TSV files found")
}

meta_list <- list()

for(f in transfer_files){

  z <- read_tsv(f)

  req <- c(
    "global_cell",
    "original_cell",
    "project_id",
    "library_key",
    "integration_compartment_primary_v1"
  )

  stopifnot(all(req %in% names(z)))

  z <- z[
    as.character(
      z$integration_compartment_primary_v1
    ) == "Erythroid",
    req,
    drop=FALSE
  ]

  if(nrow(z)){
    meta_list[[length(meta_list)+1L]] <- z
  }
}

full_meta <- do.call(rbind, meta_list)
rm(meta_list)

for(nm in c(
  "global_cell",
  "original_cell",
  "project_id",
  "library_key"
)){
  full_meta[[nm]] <- as.character(full_meta[[nm]])
}

stopifnot(
  nrow(full_meta) == N_FULL,
  !anyDuplicated(full_meta$global_cell)
)

mi <- match(
  ann$global_cell,
  full_meta$global_cell
)

stopifnot(!anyNA(mi))

ann$original_cell <- full_meta$original_cell[mi]

stopifnot(
  ann$project_id == full_meta$project_id[mi],
  ann$library_key == full_meta$library_key[mi]
)

# Update split tables after adding original_cell.
ref_ann <- ann[reference_idx, , drop=FALSE]
val_ann <- ann[validation_idx, , drop=FALSE]

# ============================================================
# Library lookup
# ============================================================

lib_file <- file.path(
  root,
  "pre_integration",
  "pilot_unintegrated_large_v1",
  "resolved_library_table.tsv"
)

lib <- read_tsv(lib_file)

stopifnot(
  all(
    c(
      "project_id",
      "library_key",
      "final_rds_resolved"
    ) %in% names(lib)
  )
)

lib$project_id <- as.character(lib$project_id)
lib$library_key <- as.character(lib$library_key)

lib$key <- paste(
  lib$project_id,
  lib$library_key,
  sep="|||"
)

# ============================================================
# Read native discovery RNA once per library
# ============================================================

library_keys <- unique(
  paste(
    ann$project_id,
    ann$library_key,
    sep="|||"
  )
)

mat_list <- list()

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

  li <- match(key, lib$key)

  if(is.na(li)){
    stop("Library lookup failed: ", key)
  }

  source_rds <- as.character(
    lib$final_rds_resolved[[li]]
  )

  if(!file.exists(source_rds)){
    stop("Missing source RDS: ", source_rds)
  }

  obj0 <- readRDS(source_rds)
  counts0 <- get_rna_counts(obj0)

  if(!all(d$original_cell %in% colnames(counts0))){
    stop("Missing cells in native RNA: ", key)
  }

  m <- counts0[
    ,
    d$original_cell,
    drop=FALSE
  ]

  # Original barcodes can repeat across libraries.
  # Global cell IDs are unique.
  colnames(m) <- d$global_cell

  mat_list[[key]] <- m

  rm(obj0, counts0, m)
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

# ============================================================
# Harmonize features
# ============================================================

common_features <- Reduce(
  intersect,
  lapply(mat_list, rownames)
)

if(length(common_features) < 10000L){
  stop(
    "Too few common features: ",
    length(common_features)
  )
}

mat_list <- lapply(
  mat_list,
  function(m){
    m[
      common_features,
      ,
      drop=FALSE
    ]
  }
)

counts <- do.call(
  cbind,
  mat_list
)

rm(mat_list)
invisible(gc())

stopifnot(
  ncol(counts) == N_DISCOVERY,
  !anyDuplicated(colnames(counts))
)

oi <- match(
  ann$global_cell,
  colnames(counts)
)

stopifnot(!anyNA(oi))

counts <- counts[
  ,
  oi,
  drop=FALSE
]

stopifnot(
  identical(
    colnames(counts),
    ann$global_cell
  )
)

# ============================================================
# Training-only project-wise HVG discovery
# ============================================================

technical_gene <- function(g){

  grepl("^MT-", g) |
    grepl("^RPL[0-9]", g) |
    grepl("^RPS[0-9]", g) |
    g %in% c(
      "MALAT1",
      "NEAT1"
    )
}

# IMPORTANT:
# HBA/HBB/HBD/HBM are NOT technical genes here.

project_objects <- list()

for(p in sort(unique(ref_ann$project_id))){

  cells <- ref_ann$global_cell[
    ref_ann$project_id == p
  ]

  if(length(cells) < 50L){
    next
  }

  obj <- CreateSeuratObject(
    counts=counts[, cells, drop=FALSE],
    min.cells=0,
    min.features=0,
    project=p
  )

  obj <- NormalizeData(
    obj,
    normalization.method="LogNormalize",
    scale.factor=10000,
    verbose=FALSE
  )

  nfeatures <- min(
    N_HVG_PER_PROJECT,
    nrow(obj) - 1L
  )

  obj <- FindVariableFeatures(
    obj,
    selection.method="vst",
    nfeatures=nfeatures,
    verbose=FALSE
  )

  project_objects[[p]] <- obj
}

if(length(project_objects) < 4L){
  stop(
    "Too few projects for feature discovery: ",
    length(project_objects)
  )
}

transfer_features <- SelectIntegrationFeatures(
  object.list=project_objects,
  nfeatures=N_INTEGRATION_FEATURES
)

transfer_features <- transfer_features[
  transfer_features %in% rownames(counts) &
    !technical_gene(transfer_features)
]

transfer_features <- unique(transfer_features)

if(length(transfer_features) < 1000L){
  stop(
    "Too few retained transfer features: ",
    length(transfer_features)
  )
}

# Explicit safety check: globins must not be excluded simply
# because they are globins.
globins_available <- intersect(
  c("HBA1","HBA2","HBB","HBD","HBM"),
  rownames(counts)
)

cat(
  "Transfer features retained = ",
  length(transfer_features),
  "\n",
  sep=""
)

cat(
  "Globins available in native matrix = ",
  paste(globins_available, collapse=","),
  "\n",
  sep=""
)

rm(project_objects)
invisible(gc())

# ============================================================
# Build reference and validation Seurat objects
# ============================================================

ref <- CreateSeuratObject(
  counts=counts[
    transfer_features,
    ref_ann$global_cell,
    drop=FALSE
  ],
  min.cells=0,
  min.features=0,
  project="erythroid_reference_train"
)

val <- CreateSeuratObject(
  counts=counts[
    transfer_features,
    val_ann$global_cell,
    drop=FALSE
  ],
  min.cells=0,
  min.features=0,
  project="erythroid_validation"
)

rm(counts)
invisible(gc())

ref <- NormalizeData(
  ref,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

val <- NormalizeData(
  val,
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
  MAX_PCS,
  length(transfer_features) - 1L,
  ncol(ref) - 1L
)

if(npcs < 10L){
  stop("Too few PCs possible: ", npcs)
}

ref <- RunPCA(
  ref,
  features=transfer_features,
  npcs=npcs,
  verbose=FALSE
)

# ============================================================
# Reference mapping
# ============================================================

anchors <- FindTransferAnchors(
  reference=ref,
  query=val,
  normalization.method="LogNormalize",
  reduction="pcaproject",
  reference.reduction="pca",
  features=transfer_features,
  dims=seq_len(npcs),
  k.anchor=5,
  k.filter=50,
  verbose=TRUE
)

ref_labels <- ref_ann$erythroid_transfer_class_v1
names(ref_labels) <- ref_ann$global_cell

pred <- TransferData(
  anchorset=anchors,
  refdata=ref_labels,
  dims=seq_len(npcs),
  k.weight=50,
  verbose=TRUE
)

pred <- as.data.frame(pred)

stopifnot(
  nrow(pred) == nrow(val_ann),
  all(rownames(pred) %in% val_ann$global_cell)
)

# Restore validation annotation order.
pi <- match(
  val_ann$global_cell,
  rownames(pred)
)

stopifnot(!anyNA(pi))

pred <- pred[pi, , drop=FALSE]

stopifnot(
  identical(
    rownames(pred),
    val_ann$global_cell
  )
)

# ============================================================
# Decode transfer classes
# ============================================================

decode_class <- function(x){

  z <- strsplit(
    as.character(x),
    "\\|\\|\\|"
  )

  if(any(lengths(z) != 3L)){
    stop("Malformed transfer class")
  }

  data.frame(
    core=vapply(z, `[[`, character(1), 1L),
    maturation=vapply(z, `[[`, character(1), 2L),
    state=vapply(z, `[[`, character(1), 3L),
    stringsAsFactors=FALSE
  )
}

actual_dec <- decode_class(
  val_ann$erythroid_transfer_class_v1
)

pred_dec <- decode_class(
  pred$predicted.id
)

result <- data.frame(
  global_cell=val_ann$global_cell,
  project_id=val_ann$project_id,
  library_key=val_ann$library_key,
  cluster_r0p4=val_ann$cluster_r0p4,

  actual_class=
    val_ann$erythroid_transfer_class_v1,

  predicted_class=
    as.character(pred$predicted.id),

  prediction_score=
    as.numeric(pred$prediction.score.max),

  actual_core=
    actual_dec$core,

  predicted_core=
    pred_dec$core,

  actual_maturation=
    actual_dec$maturation,

  predicted_maturation=
    pred_dec$maturation,

  actual_state=
    actual_dec$state,

  predicted_state=
    pred_dec$state,

  stringsAsFactors=FALSE
)

result$correct_exact <-
  result$actual_class ==
  result$predicted_class

result$correct_core <-
  result$actual_core ==
  result$predicted_core

result$correct_maturation <-
  result$actual_maturation ==
  result$predicted_maturation

result$correct_state <-
  result$actual_state ==
  result$predicted_state

# ============================================================
# Overall
# ============================================================

overall <- data.frame(
  metric=c(
    "exact_accuracy",
    "core_accuracy",
    "maturation_accuracy",
    "state_accuracy",
    "median_prediction_score"
  ),
  value=c(
    mean(result$correct_exact),
    mean(result$correct_core),
    mean(result$correct_maturation),
    mean(result$correct_state),
    median(result$prediction_score)
  ),
  stringsAsFactors=FALSE
)

# ============================================================
# Threshold audit
# ============================================================

thresholds <- c(
  0,
  0.60,
  0.70,
  0.75,
  0.80,
  0.85,
  0.90,
  0.95
)

threshold_summary <- do.call(
  rbind,
  lapply(
    thresholds,
    function(th){

      keep <- result$prediction_score >= th

      data.frame(
        threshold=th,
        n_accepted=sum(keep),
        coverage=mean(keep),

        exact_accuracy=
          if(any(keep))
            mean(result$correct_exact[keep])
          else NA_real_,

        core_accuracy=
          if(any(keep))
            mean(result$correct_core[keep])
          else NA_real_,

        maturation_accuracy=
          if(any(keep))
            mean(result$correct_maturation[keep])
          else NA_real_,

        state_accuracy=
          if(any(keep))
            mean(result$correct_state[keep])
          else NA_real_,

        stringsAsFactors=FALSE
      )
    }
  )
)

# ============================================================
# Actual-class recall by threshold
# ============================================================

classes <- sort(
  unique(result$actual_class)
)

class_recall <- do.call(
  rbind,
  lapply(
    thresholds,
    function(th){

      do.call(
        rbind,
        lapply(
          classes,
          function(cl){

            actual <- result$actual_class == cl
            accepted <- result$prediction_score >= th
            correct <- result$predicted_class == cl

            data.frame(
              threshold=th,
              actual_class=cl,
              n_actual=sum(actual),

              n_correct_accepted=
                sum(
                  actual &
                  accepted &
                  correct
                ),

              accepted_recall=
                sum(
                  actual &
                  accepted &
                  correct
                ) /
                sum(actual),

              stringsAsFactors=FALSE
            )
          }
        )
      )
    }
  )
)

# ============================================================
# Predicted-class precision by threshold
# ============================================================

predicted_classes <- sort(
  unique(result$predicted_class)
)

class_precision <- do.call(
  rbind,
  lapply(
    thresholds,
    function(th){

      do.call(
        rbind,
        lapply(
          predicted_classes,
          function(cl){

            accepted_prediction <-
              result$predicted_class == cl &
              result$prediction_score >= th

            n <- sum(accepted_prediction)

            data.frame(
              threshold=th,
              predicted_class=cl,
              n_accepted=n,

              precision=
                if(n > 0L){
                  mean(
                    result$actual_class[
                      accepted_prediction
                    ] == cl
                  )
                } else {
                  NA_real_
                },

              stringsAsFactors=FALSE
            )
          }
        )
      )
    }
  )
)

# ============================================================
# Dimension precision
# ============================================================

dimension_precision_one <- function(
  actual,
  predicted,
  score,
  dimension
){

  labs <- sort(unique(predicted))

  do.call(
    rbind,
    lapply(
      thresholds,
      function(th){

        do.call(
          rbind,
          lapply(
            labs,
            function(lb){

              keep <-
                predicted == lb &
                score >= th

              n <- sum(keep)

              data.frame(
                dimension=dimension,
                threshold=th,
                predicted_label=lb,
                n_accepted=n,

                precision=
                  if(n > 0L){
                    mean(actual[keep] == lb)
                  } else {
                    NA_real_
                  },

                stringsAsFactors=FALSE
              )
            }
          )
        )
      }
    )
  )
}

dimension_precision <- rbind(
  dimension_precision_one(
    result$actual_core,
    result$predicted_core,
    result$prediction_score,
    "core"
  ),

  dimension_precision_one(
    result$actual_maturation,
    result$predicted_maturation,
    result$prediction_score,
    "maturation"
  ),

  dimension_precision_one(
    result$actual_state,
    result$predicted_state,
    result$prediction_score,
    "state"
  )
)

# ============================================================
# Project performance
# ============================================================

project_split <- split(
  result,
  result$project_id
)

project_performance <- do.call(
  rbind,
  lapply(
    project_split,
    function(d){

      data.frame(
        project_id=
          as.character(
            d$project_id[[1]]
          ),

        n_validation=nrow(d),

        exact_accuracy=
          mean(d$correct_exact),

        core_accuracy=
          mean(d$correct_core),

        maturation_accuracy=
          mean(d$correct_maturation),

        state_accuracy=
          mean(d$correct_state),

        median_prediction_score=
          median(d$prediction_score),

        stringsAsFactors=FALSE
      )
    }
  )
)

rownames(project_performance) <- NULL

# ============================================================
# Correct vs incorrect score by actual class
# ============================================================

score_by_actual_class <- do.call(
  rbind,
  lapply(
    classes,
    function(cl){

      d <- result[
        result$actual_class == cl,
        ,
        drop=FALSE
      ]

      ok <- d$correct_exact

      data.frame(
        actual_class=cl,
        n=nrow(d),
        correct_fraction=mean(ok),

        median_score_correct=
          if(any(ok))
            median(d$prediction_score[ok])
          else NA_real_,

        median_score_incorrect=
          if(any(!ok))
            median(d$prediction_score[!ok])
          else NA_real_,

        q25_score_correct=
          if(any(ok))
            unname(
              quantile(
                d$prediction_score[ok],
                0.25
              )
            )
          else NA_real_,

        q75_score_incorrect=
          if(any(!ok))
            unname(
              quantile(
                d$prediction_score[!ok],
                0.75
              )
            )
          else NA_real_,

        stringsAsFactors=FALSE
      )
    }
  )
)

# ============================================================
# Split manifest
# ============================================================

split_manifest <- aggregate(
  rep(1L, nrow(ann)),
  by=list(
    project_id=ann$project_id,
    transfer_class=ann$erythroid_transfer_class_v1,
    split=ifelse(
      seq_len(nrow(ann)) %in% validation_idx,
      "validation",
      "reference_train"
    )
  ),
  FUN=sum
)

names(split_manifest)[4] <- "n_cells"

# ============================================================
# Save
# ============================================================

write.table(
  overall,
  file.path(
    out_dir,
    "erythroid_transfer_calibration_overall_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  threshold_summary,
  file.path(
    out_dir,
    "erythroid_transfer_threshold_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  class_recall,
  file.path(
    out_dir,
    "erythroid_transfer_class_recall_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  class_precision,
  file.path(
    out_dir,
    "erythroid_transfer_class_precision_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  dimension_precision,
  file.path(
    out_dir,
    "erythroid_transfer_dimension_precision_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  project_performance,
  file.path(
    out_dir,
    "erythroid_transfer_project_performance_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  score_by_actual_class,
  file.path(
    out_dir,
    "erythroid_transfer_score_by_actual_class_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  split_manifest,
  file.path(
    out_dir,
    "erythroid_transfer_split_manifest_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  data.frame(
    feature=transfer_features,
    stringsAsFactors=FALSE
  ),
  file.path(
    out_dir,
    "erythroid_transfer_features_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write_gz_tsv(
  result,
  file.path(
    out_dir,
    paste0(
      "erythroid_transfer_validation_predictions_v1__",
      tag,
      ".tsv.gz"
    )
  )
)

calibration <- list(
  calibration_version=
    "Erythroid_discovery_transfer_calibration_v1",

  n_discovery=N_DISCOVERY,
  n_reference_train=nrow(ref_ann),
  n_validation=nrow(val_ann),

  transfer_features=transfer_features,
  pca_features=transfer_features,
  npcs=npcs,

  reference_cells=ref_ann$global_cell,
  validation_cells=val_ann$global_cell,

  overall=overall,
  threshold_summary=threshold_summary,
  class_recall=class_recall,
  class_precision=class_precision,
  dimension_precision=dimension_precision,
  project_performance=project_performance,
  score_by_actual_class=score_by_actual_class,

  source_freeze_dir=freeze_dir,
  source_freeze_rds=freeze_files[[1]]
)

rds_file <- file.path(
  out_dir,
  paste0(
    "erythroid_discovery_transfer_calibration_v1__",
    tag,
    ".rds"
  )
)

tmp_rds <- paste0(rds_file, ".tmp")

saveRDS(
  calibration,
  tmp_rds,
  compress=TRUE
)

if(!file.rename(tmp_rds, rds_file)){
  stop("Atomic RDS rename failed")
}

chk <- readRDS(rds_file)

stopifnot(
  chk$n_discovery == N_DISCOVERY,
  chk$n_reference_train == nrow(ref_ann),
  chk$n_validation == nrow(val_ann),
  identical(
    chk$transfer_features,
    transfer_features
  )
)

pred_file <- file.path(
  out_dir,
  paste0(
    "erythroid_transfer_validation_predictions_v1__",
    tag,
    ".tsv.gz"
  )
)

pred_chk <- read_tsv(pred_file)

stopifnot(
  nrow(pred_chk) == nrow(result),
  identical(
    as.character(pred_chk$global_cell),
    as.character(result$global_cell)
  )
)

sha_files <- c(
  file.path(out_dir, "erythroid_transfer_calibration_overall_v1.tsv"),
  file.path(out_dir, "erythroid_transfer_threshold_summary_v1.tsv"),
  file.path(out_dir, "erythroid_transfer_class_recall_v1.tsv"),
  file.path(out_dir, "erythroid_transfer_class_precision_v1.tsv"),
  file.path(out_dir, "erythroid_transfer_dimension_precision_v1.tsv"),
  file.path(out_dir, "erythroid_transfer_project_performance_v1.tsv"),
  file.path(out_dir, "erythroid_transfer_score_by_actual_class_v1.tsv"),
  file.path(out_dir, "erythroid_transfer_split_manifest_v1.tsv"),
  file.path(out_dir, "erythroid_transfer_features_v1.tsv"),
  pred_file,
  rds_file
)

sha_out <- file.path(
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
  stdout=sha_out
)

if(status != 0L){
  stop("sha256sum failed")
}

writeLines(
  c(
    "PASS",
    "Erythroid discovery held-out transfer calibration v1",
    paste0("n_discovery=", N_DISCOVERY),
    paste0("n_reference_train=", nrow(ref_ann)),
    paste0("n_validation=", nrow(val_ann)),
    paste0("n_transfer_features=", length(transfer_features)),
    paste0("n_pcs=", npcs),
    "",
    "split=project x frozen transfer class stratified",
    "validation annotations withheld from training",
    "native SoupX-corrected RNA",
    "normalization=LogNormalize",
    "mapping=Seurat pcaproject + TransferData",
    "",
    "core / maturation / state evaluated separately",
    "class precision and recall audited across score thresholds",
    "full-query transfer NOT performed",
    "",
    "globin genes not removed as technical features by rule",
    "integrated assay NOT used as transfer expression input",
    "cellranger_count not accessed",
    "gzip_integrity=PASS",
    "re_read_check=PASS",
    "freeze_rds_re_read=PASS",
    "sha256_manifest=PASS"
  ),
  done_file
)

cat("\n===== OVERALL =====\n")
print(overall, row.names=FALSE)

cat("\n===== THRESHOLDS =====\n")
print(threshold_summary, row.names=FALSE)

cat("\n===== CLASS PRECISION: >=0.70 / >=0.80 / >=0.90 =====\n")
print(
  class_precision[
    class_precision$threshold %in% c(0.70,0.80,0.90),
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat("\n===== CLASS RECALL: >=0.70 / >=0.80 / >=0.90 =====\n")
print(
  class_recall[
    class_recall$threshold %in% c(0.70,0.80,0.90),
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat("\n===== DIMENSION PRECISION: >=0.70 / >=0.80 / >=0.90 =====\n")
print(
  dimension_precision[
    dimension_precision$threshold %in% c(0.70,0.80,0.90),
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat("\n===== PROJECT PERFORMANCE =====\n")
print(project_performance, row.names=FALSE)

cat("\n===== SCORE BY ACTUAL CLASS =====\n")
print(score_by_actual_class, row.names=FALSE)

cat("\nOUT_DIR=", out_dir, "\n", sep="")
cat("\nPASS: Erythroid transfer calibration completed\n")
