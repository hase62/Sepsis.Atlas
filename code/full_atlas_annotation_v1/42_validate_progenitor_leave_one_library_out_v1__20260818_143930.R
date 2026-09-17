#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3L){
  stop("Usage: script <root> <tag> <freeze_dir>")
}

root <- normalizePath(
  args[[1]],
  mustWork=TRUE
)

tag <- args[[2]]

freeze_dir <- normalizePath(
  args[[3]],
  mustWork=TRUE
)

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

# ============================================================
# Frozen validation design
# ============================================================

DISCOVERY_PROJECT <- "GSE216007"
N_EXPECTED <- 26962L
N_LIBRARIES_EXPECTED <- 5L

N_HVG <- 3000L
N_PCS <- 30L

SCORE_THRESHOLDS <- c(
  0,
  0.50,
  0.60,
  0.70,
  0.80,
  0.90
)

# ============================================================
# Paths
# ============================================================

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "progenitor_leave_one_library_out_validation_v1__",
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
  "PROGENITOR_LEAVE_ONE_LIBRARY_OUT_VALIDATION_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop("Already completed: ", done_file)
}

# ============================================================
# I/O helpers
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

safe_mean <- function(x){

  if(!length(x)){
    return(NA_real_)
  }

  mean(x)
}

safe_median <- function(x){

  if(!length(x)){
    return(NA_real_)
  }

  median(x)
}

# ============================================================
# Classification metrics
# ============================================================

class_metrics <- function(actual, predicted){

  lev <- sort(
    unique(
      c(
        as.character(actual),
        as.character(predicted)
      )
    )
  )

  a <- factor(
    actual,
    levels=lev
  )

  p <- factor(
    predicted,
    levels=lev
  )

  cm <- table(
    actual=a,
    predicted=p
  )

  tp <- diag(cm)
  support <- rowSums(cm)
  predicted_n <- colSums(cm)

  recall <- ifelse(
    support > 0,
    tp / support,
    NA_real_
  )

  precision <- ifelse(
    predicted_n > 0,
    tp / predicted_n,
    NA_real_
  )

  f1 <- ifelse(
    is.finite(precision) &
      is.finite(recall) &
      (precision + recall) > 0,
    2 * precision * recall /
      (precision + recall),
    NA_real_
  )

  data.frame(
    class=lev,
    support=as.integer(support),
    predicted_n=as.integer(predicted_n),
    precision=as.numeric(precision),
    recall=as.numeric(recall),
    F1=as.numeric(f1),
    stringsAsFactors=FALSE
  )
}

macro_f1 <- function(actual, predicted){

  d <- class_metrics(
    actual,
    predicted
  )

  mean(
    d$F1,
    na.rm=TRUE
  )
}

confusion_long <- function(actual, predicted){

  lev <- sort(
    unique(
      c(
        as.character(actual),
        as.character(predicted)
      )
    )
  )

  cm <- as.data.frame(
    table(
      actual=factor(
        actual,
        levels=lev
      ),
      predicted=factor(
        predicted,
        levels=lev
      )
    ),
    stringsAsFactors=FALSE
  )

  cm[
    cm$Freq > 0,
    ,
    drop=FALSE
  ]
}

# ============================================================
# Read frozen Step41 discovery annotation
# ============================================================

completion_file <- file.path(
  freeze_dir,
  "PROGENITOR_DISCOVERY_ANNOTATION_FREEZE_COMPLETE.ok"
)

stopifnot(
  file.exists(completion_file)
)

ann_files <- list.files(
  freeze_dir,
  pattern=
    "^progenitor_discovery_cell_annotation_v1__.*\\.tsv\\.gz$",
  full.names=TRUE
)

stopifnot(
  length(ann_files) == 1L
)

ann <- read_tsv(
  ann_files[[1]]
)

required_ann <- c(
  "global_cell",
  "project_id",
  "library_key",
  "discovery_cluster",
  "core_identity",
  "lineage_priming_axis",
  "state",
  "annotation_confidence",
  "evidence_flag"
)

stopifnot(
  all(required_ann %in% names(ann)),
  nrow(ann) == N_EXPECTED,
  !anyDuplicated(ann$global_cell)
)

for(nm in required_ann){
  ann[[nm]] <- as.character(
    ann[[nm]]
  )
}

stopifnot(
  all(
    ann$project_id ==
      DISCOVERY_PROJECT
  )
)

libraries <- sort(
  unique(
    ann$library_key
  )
)

stopifnot(
  length(libraries) ==
    N_LIBRARIES_EXPECTED
)

# Transfer target deliberately excludes orthogonal state.
ann$transfer_class <- paste(
  ann$core_identity,
  ann$lineage_priming_axis,
  sep="|||AXIS|||"
)

class_lookup <- unique(
  ann[
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
# Recover original_cell
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

if(!length(transfer_files)){
  stop("No annotation transfer files")
}

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

meta <- do.call(
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
  meta[[nm]] <- as.character(
    meta[[nm]]
  )
}

stopifnot(
  nrow(meta) == N_EXPECTED,
  !anyDuplicated(meta$global_cell)
)

mi <- match(
  ann$global_cell,
  meta$global_cell
)

stopifnot(
  !anyNA(mi)
)

ann$original_cell <-
  meta$original_cell[mi]

stopifnot(
  ann$project_id ==
    meta$project_id[mi],
  ann$library_key ==
    meta$library_key[mi]
)

rm(meta, mi)
invisible(gc())

# ============================================================
# Source library table
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
# Load all five GSE216007 libraries in native SoupX space
# ============================================================

counts_by_library <- list()

for(ii in seq_along(libraries)){

  lk <- libraries[[ii]]

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
      "Library lookup failed: ",
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

  d <- ann[
    ann$library_key == lk,
    ,
    drop=FALSE
  ]

  obj <- readRDS(
    source_rds
  )

  counts0 <- get_rna_counts(
    obj
  )

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

  colnames(m) <-
    d$global_cell

  counts_by_library[[lk]] <- m

  cat(
    sprintf(
      "READ %d/%d %s cells=%d features=%d\n",
      ii,
      length(libraries),
      lk,
      ncol(m),
      nrow(m)
    )
  )

  rm(
    obj,
    counts0,
    m
  )

  invisible(gc())
}

# ============================================================
# Common feature space
# ============================================================

feature_sets <- lapply(
  counts_by_library,
  rownames
)

common_features <- Reduce(
  intersect,
  feature_sets
)

first_features <- rownames(
  counts_by_library[[1]]
)

common_features <- first_features[
  first_features %in%
    common_features
]

if(length(common_features) < 10000L){
  stop(
    "Too few common features: ",
    length(common_features)
  )
}

counts_by_library <- lapply(
  counts_by_library,
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
  counts_by_library
)

rm(
  counts_by_library,
  feature_sets
)

invisible(gc())

ci <- match(
  ann$global_cell,
  colnames(counts)
)

stopifnot(
  !anyNA(ci)
)

counts <- counts[
  ,
  ci,
  drop=FALSE
]

stopifnot(
  identical(
    colnames(counts),
    ann$global_cell
  )
)

cat(
  "Common features=",
  nrow(counts),
  "\n",
  sep=""
)

# ============================================================
# Leave-one-library-out folds
#
# IMPORTANT:
# HVGs are selected from the four training libraries ONLY.
# Held-out library is never used for feature selection/PCA.
# ============================================================

prediction_list <- list()
fold_metric_list <- list()
fold_manifest_list <- list()
fold_confusion_list <- list()

for(fi in seq_along(libraries)){

  heldout <- libraries[[fi]]

  cat(
    "\n============================================\n",
    "FOLD ",
    fi,
    "/",
    length(libraries),
    " HELDOUT=",
    heldout,
    "\n",
    "============================================\n",
    sep=""
  )

  test_idx <- which(
    ann$library_key == heldout
  )

  train_idx <- which(
    ann$library_key != heldout
  )

  ref_ann <- ann[
    train_idx,
    ,
    drop=FALSE
  ]

  val_ann <- ann[
    test_idx,
    ,
    drop=FALSE
  ]

  stopifnot(
    nrow(ref_ann) +
      nrow(val_ann) ==
      N_EXPECTED,
    !length(
      intersect(
        ref_ann$global_cell,
        val_ann$global_cell
      )
    ),
    all(
      val_ann$library_key ==
        heldout
    ),
    !any(
      ref_ann$library_key ==
        heldout
    )
  )

  # Every canonical transfer class must exist in training.
  train_classes <- unique(
    ref_ann$transfer_class
  )

  test_classes <- unique(
    val_ann$transfer_class
  )

  missing_train_classes <- setdiff(
    test_classes,
    train_classes
  )

  if(length(missing_train_classes)){
    stop(
      "Held-out fold has classes absent from training: ",
      heldout,
      " :: ",
      paste(
        missing_train_classes,
        collapse=","
      )
    )
  }

  # ----------------------------------------------------------
  # Training object for HVG discovery.
  # ----------------------------------------------------------

  ref_all <- CreateSeuratObject(
    counts=counts[
      ,
      ref_ann$global_cell,
      drop=FALSE
    ],
    min.cells=0,
    min.features=0,
    project="progenitor_LOO_reference"
  )

  ref_all <- NormalizeData(
    ref_all,
    normalization.method="LogNormalize",
    scale.factor=10000,
    verbose=FALSE
  )

  ref_all <- FindVariableFeatures(
    ref_all,
    selection.method="vst",
    nfeatures=N_HVG,
    verbose=FALSE
  )

  hvg_raw <- VariableFeatures(
    ref_all
  )

  transfer_features <- hvg_raw[
    !technical_gene(
      hvg_raw
    )
  ]

  if(length(transfer_features) < 1500L){
    stop(
      "Too few transfer features in fold ",
      heldout,
      ": ",
      length(transfer_features)
    )
  }

  # ----------------------------------------------------------
  # Rebuild compact reference/query using only frozen
  # training-derived transfer features.
  # ----------------------------------------------------------

  ref <- CreateSeuratObject(
    counts=counts[
      transfer_features,
      ref_ann$global_cell,
      drop=FALSE
    ],
    min.cells=0,
    min.features=0,
    project="progenitor_LOO_reference"
  )

  val <- CreateSeuratObject(
    counts=counts[
      transfer_features,
      val_ann$global_cell,
      drop=FALSE
    ],
    min.cells=0,
    min.features=0,
    project="progenitor_LOO_validation"
  )

  rm(ref_all)
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
    N_PCS,
    length(transfer_features) - 1L,
    ncol(ref) - 1L
  )

  if(npcs < 10L){
    stop(
      "Too few PCs in fold: ",
      heldout
    )
  }

  ref <- RunPCA(
    ref,
    features=transfer_features,
    npcs=npcs,
    verbose=FALSE
  )

  # ----------------------------------------------------------
  # pcaproject mapping
  # ----------------------------------------------------------

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

  # ----------------------------------------------------------
  # Core + lineage axis transfer
  # ----------------------------------------------------------

  ref_class <- ref_ann$transfer_class
  names(ref_class) <-
    ref_ann$global_cell

  pred_class <- TransferData(
    anchorset=anchors,
    refdata=ref_class,
    dims=seq_len(npcs),
    k.weight=50,
    verbose=TRUE
  )

  pred_class <- as.data.frame(
    pred_class
  )

  pi <- match(
    val_ann$global_cell,
    rownames(pred_class)
  )

  stopifnot(
    !anyNA(pi)
  )

  pred_class <- pred_class[
    pi,
    ,
    drop=FALSE
  ]

  # ----------------------------------------------------------
  # Orthogonal state transfer
  # ----------------------------------------------------------

  ref_state <- ref_ann$state
  names(ref_state) <-
    ref_ann$global_cell

  pred_state <- TransferData(
    anchorset=anchors,
    refdata=ref_state,
    dims=seq_len(npcs),
    k.weight=50,
    verbose=FALSE
  )

  pred_state <- as.data.frame(
    pred_state
  )

  si <- match(
    val_ann$global_cell,
    rownames(pred_state)
  )

  stopifnot(
    !anyNA(si)
  )

  pred_state <- pred_state[
    si,
    ,
    drop=FALSE
  ]

  # ----------------------------------------------------------
  # Decode predicted canonical class
  # ----------------------------------------------------------

  predicted_class <- as.character(
    pred_class$predicted.id
  )

  li <- match(
    predicted_class,
    class_lookup$transfer_class
  )

  stopifnot(
    !anyNA(li)
  )

  predicted_core <-
    class_lookup$core_identity[li]

  predicted_axis <-
    class_lookup$lineage_priming_axis[li]

  result <- data.frame(

    global_cell=
      val_ann$global_cell,

    heldout_library=
      heldout,

    discovery_cluster=
      val_ann$discovery_cluster,

    actual_class=
      val_ann$transfer_class,

    predicted_class=
      predicted_class,

    class_prediction_score=
      as.numeric(
        pred_class$prediction.score.max
      ),

    actual_core=
      val_ann$core_identity,

    predicted_core=
      predicted_core,

    actual_axis=
      val_ann$lineage_priming_axis,

    predicted_axis=
      predicted_axis,

    actual_state=
      val_ann$state,

    predicted_state=
      as.character(
        pred_state$predicted.id
      ),

    state_prediction_score=
      as.numeric(
        pred_state$prediction.score.max
      ),

    stringsAsFactors=FALSE
  )

  result$correct_exact <-
    result$actual_class ==
      result$predicted_class

  result$correct_core <-
    result$actual_core ==
      result$predicted_core

  result$correct_axis <-
    result$actual_axis ==
      result$predicted_axis

  result$correct_state <-
    result$actual_state ==
      result$predicted_state

  prediction_list[[
    length(prediction_list)+1L
  ]] <- result

  exact_class_metrics <- class_metrics(
    result$actual_class,
    result$predicted_class
  )

  fold_metric_list[[
    length(fold_metric_list)+1L
  ]] <- data.frame(
    heldout_library=heldout,
    n_train=nrow(ref_ann),
    n_test=nrow(val_ann),
    n_transfer_features=
      length(transfer_features),
    n_pcs=npcs,
    exact_accuracy=
      mean(result$correct_exact),
    core_accuracy=
      mean(result$correct_core),
    axis_accuracy=
      mean(result$correct_axis),
    state_accuracy=
      mean(result$correct_state),
    macro_F1_exact=
      mean(
        exact_class_metrics$F1,
        na.rm=TRUE
      ),
    median_class_prediction_score=
      median(
        result$class_prediction_score
      ),
    median_state_prediction_score=
      median(
        result$state_prediction_score
      ),
    stringsAsFactors=FALSE
  )

  fold_manifest_list[[
    length(fold_manifest_list)+1L
  ]] <- data.frame(
    heldout_library=heldout,
    n_train=nrow(ref_ann),
    n_test=nrow(val_ann),
    n_train_libraries=
      length(
        unique(
          ref_ann$library_key
        )
      ),
    n_test_libraries=
      length(
        unique(
          val_ann$library_key
        )
      ),
    n_train_classes=
      length(
        unique(
          ref_ann$transfer_class
        )
      ),
    n_test_classes=
      length(
        unique(
          val_ann$transfer_class
        )
      ),
    n_transfer_features=
      length(transfer_features),
    stringsAsFactors=FALSE
  )

  cm <- confusion_long(
    result$actual_class,
    result$predicted_class
  )

  cm$heldout_library <- heldout

  fold_confusion_list[[
    length(fold_confusion_list)+1L
  ]] <- cm[
    ,
    c(
      "heldout_library",
      "actual",
      "predicted",
      "Freq"
    ),
    drop=FALSE
  ]

  cat(
    sprintf(
      paste0(
        "FOLD_RESULT heldout=%s ",
        "n=%d exact=%.4f core=%.4f ",
        "axis=%.4f state=%.4f ",
        "macroF1=%.4f score=%.4f\n"
      ),
      heldout,
      nrow(result),
      mean(result$correct_exact),
      mean(result$correct_core),
      mean(result$correct_axis),
      mean(result$correct_state),
      mean(
        exact_class_metrics$F1,
        na.rm=TRUE
      ),
      median(
        result$class_prediction_score
      )
    )
  )

  rm(
    ref,
    val,
    anchors,
    pred_class,
    pred_state,
    result
  )

  invisible(gc())
}

# ============================================================
# All cells must appear exactly once as held-out validation
# ============================================================

predictions <- do.call(
  rbind,
  prediction_list
)

rownames(predictions) <- NULL

stopifnot(
  nrow(predictions) ==
    N_EXPECTED,
  !anyDuplicated(
    predictions$global_cell
  ),
  setequal(
    predictions$global_cell,
    ann$global_cell
  )
)

fold_metrics <- do.call(
  rbind,
  fold_metric_list
)

fold_manifest <- do.call(
  rbind,
  fold_manifest_list
)

fold_confusion <- do.call(
  rbind,
  fold_confusion_list
)

# Restore frozen discovery order.
oi <- match(
  ann$global_cell,
  predictions$global_cell
)

stopifnot(
  !anyNA(oi)
)

predictions <- predictions[
  oi,
  ,
  drop=FALSE
]

stopifnot(
  identical(
    predictions$global_cell,
    ann$global_cell
  )
)

# ============================================================
# Overall metrics
# ============================================================

exact_metrics <- class_metrics(
  predictions$actual_class,
  predictions$predicted_class
)

overall <- data.frame(
  metric=c(
    "exact_accuracy",
    "core_accuracy",
    "axis_accuracy",
    "state_accuracy",
    "macro_F1_exact",
    "median_class_prediction_score",
    "median_state_prediction_score"
  ),
  value=c(
    mean(predictions$correct_exact),
    mean(predictions$correct_core),
    mean(predictions$correct_axis),
    mean(predictions$correct_state),
    mean(
      exact_metrics$F1,
      na.rm=TRUE
    ),
    median(
      predictions$class_prediction_score
    ),
    median(
      predictions$state_prediction_score
    )
  ),
  stringsAsFactors=FALSE
)

# ============================================================
# Per-class metrics with canonical decoding
# ============================================================

exact_metrics$transfer_class <-
  exact_metrics$class

ci <- match(
  exact_metrics$transfer_class,
  class_lookup$transfer_class
)

exact_metrics$core_identity <-
  class_lookup$core_identity[ci]

exact_metrics$lineage_priming_axis <-
  class_lookup$lineage_priming_axis[ci]

exact_metrics <- exact_metrics[
  ,
  c(
    "transfer_class",
    "core_identity",
    "lineage_priming_axis",
    "support",
    "predicted_n",
    "precision",
    "recall",
    "F1"
  ),
  drop=FALSE
]

# ============================================================
# Overall confusion
# ============================================================

overall_confusion <- confusion_long(
  predictions$actual_class,
  predictions$predicted_class
)

# ============================================================
# Prediction-score threshold audit
# ============================================================

threshold_rows <- list()

for(th in SCORE_THRESHOLDS){

  keep <-
    predictions$class_prediction_score >= th

  n_keep <- sum(keep)

  threshold_rows[[
    length(threshold_rows)+1L
  ]] <- data.frame(
    threshold=th,
    n_retained=n_keep,
    coverage=
      n_keep /
      nrow(predictions),
    exact_accuracy=
      if(n_keep)
        mean(
          predictions$correct_exact[keep]
        )
      else
        NA_real_,
    core_accuracy=
      if(n_keep)
        mean(
          predictions$correct_core[keep]
        )
      else
        NA_real_,
    axis_accuracy=
      if(n_keep)
        mean(
          predictions$correct_axis[keep]
        )
      else
        NA_real_,
    macro_F1_exact=
      if(n_keep)
        macro_f1(
          predictions$actual_class[keep],
          predictions$predicted_class[keep]
        )
      else
        NA_real_,
    stringsAsFactors=FALSE
  )
}

threshold_audit <- do.call(
  rbind,
  threshold_rows
)

# ============================================================
# State score threshold audit
# ============================================================

state_threshold_rows <- list()

for(th in SCORE_THRESHOLDS){

  keep <-
    predictions$state_prediction_score >= th

  n_keep <- sum(keep)

  state_threshold_rows[[
    length(state_threshold_rows)+1L
  ]] <- data.frame(
    threshold=th,
    n_retained=n_keep,
    coverage=
      n_keep /
      nrow(predictions),
    state_accuracy=
      if(n_keep)
        mean(
          predictions$correct_state[keep]
        )
      else
        NA_real_,
    stringsAsFactors=FALSE
  )
}

state_threshold_audit <- do.call(
  rbind,
  state_threshold_rows
)

# ============================================================
# Per-library class metrics
# ============================================================

per_library_class_list <- list()

for(lk in libraries){

  d <- predictions[
    predictions$heldout_library == lk,
    ,
    drop=FALSE
  ]

  m <- class_metrics(
    d$actual_class,
    d$predicted_class
  )

  m$heldout_library <- lk

  per_library_class_list[[
    length(per_library_class_list)+1L
  ]] <- m[
    ,
    c(
      "heldout_library",
      "class",
      "support",
      "predicted_n",
      "precision",
      "recall",
      "F1"
    ),
    drop=FALSE
  ]
}

per_library_class <- do.call(
  rbind,
  per_library_class_list
)

# ============================================================
# Write outputs
# ============================================================

prediction_file <- file.path(
  out_dir,
  paste0(
    "progenitor_loo_predictions_v1__",
    tag,
    ".tsv.gz"
  )
)

write_gz_tsv(
  predictions,
  prediction_file
)

write.table(
  fold_manifest,
  file.path(
    out_dir,
    "progenitor_loo_fold_manifest_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  fold_metrics,
  file.path(
    out_dir,
    "progenitor_loo_fold_metrics_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  overall,
  file.path(
    out_dir,
    "progenitor_loo_overall_metrics_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  exact_metrics,
  file.path(
    out_dir,
    "progenitor_loo_per_class_metrics_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  per_library_class,
  file.path(
    out_dir,
    "progenitor_loo_per_library_class_metrics_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  threshold_audit,
  file.path(
    out_dir,
    "progenitor_loo_class_score_threshold_audit_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  state_threshold_audit,
  file.path(
    out_dir,
    "progenitor_loo_state_score_threshold_audit_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  overall_confusion,
  file.path(
    out_dir,
    "progenitor_loo_overall_confusion_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  fold_confusion,
  file.path(
    out_dir,
    "progenitor_loo_fold_confusion_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Small validation freeze object
# ============================================================

validation_freeze <- list(
  validation_version=
    "Progenitor_leave_one_library_out_validation_v1",

  discovery_project=
    DISCOVERY_PROJECT,

  n_cells=
    N_EXPECTED,

  n_libraries=
    N_LIBRARIES_EXPECTED,

  design=
    "five-fold leave-one-library-out",

  transfer_target=
    "core_identity + lineage_priming_axis",

  state_target=
    "orthogonal state",

  feature_selection=
    "training libraries only within each fold",

  fold_manifest=
    fold_manifest,

  fold_metrics=
    fold_metrics,

  overall_metrics=
    overall,

  per_class_metrics=
    exact_metrics,

  class_threshold_audit=
    threshold_audit,

  state_threshold_audit=
    state_threshold_audit,

  source_freeze_dir=
    freeze_dir
)

rds_file <- file.path(
  out_dir,
  paste0(
    "progenitor_leave_one_library_out_validation_v1__",
    tag,
    ".rds"
  )
)

tmp_rds <- paste0(
  rds_file,
  ".tmp"
)

saveRDS(
  validation_freeze,
  tmp_rds,
  compress=TRUE
)

if(!file.rename(
  tmp_rds,
  rds_file
)){
  stop("Atomic RDS rename failed")
}

chk_rds <- readRDS(
  rds_file
)

stopifnot(
  chk_rds$n_cells ==
    N_EXPECTED,
  chk_rds$n_libraries ==
    N_LIBRARIES_EXPECTED
)

# ============================================================
# gzip reread
# ============================================================

chk_pred <- read_tsv(
  prediction_file
)

stopifnot(
  nrow(chk_pred) ==
    N_EXPECTED,
  identical(
    as.character(
      chk_pred$global_cell
    ),
    as.character(
      predictions$global_cell
    )
  )
)

# ============================================================
# SHA256
# ============================================================

sha_files <- c(
  prediction_file,
  file.path(
    out_dir,
    "progenitor_loo_fold_manifest_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_loo_fold_metrics_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_loo_overall_metrics_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_loo_per_class_metrics_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_loo_per_library_class_metrics_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_loo_class_score_threshold_audit_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_loo_state_score_threshold_audit_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_loo_overall_confusion_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_loo_fold_confusion_v1.tsv"
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
  stop("sha256sum creation failed")
}

status <- system2(
  "sha256sum",
  c("-c", sha_file),
  stdout=FALSE,
  stderr=FALSE
)

if(status != 0L){
  stop("sha256 verification failed")
}

# ============================================================
# Completion marker LAST
# ============================================================

writeLines(
  c(
    "PASS",
    "Progenitor leave-one-library-out validation v1",
    "validation_version=Progenitor_leave_one_library_out_validation_v1",
    "discovery_project=GSE216007",
    "n_cells=26962",
    "n_libraries=5",
    "validation=five-fold leave-one-library-out",
    "",
    "held-out library excluded from training",
    "held-out library excluded from HVG selection",
    "held-out library excluded from PCA construction",
    "training-only HVG selection per fold",
    "normalization=Seurat LogNormalize",
    "reduction=pcaproject",
    "transfer=FindTransferAnchors + TransferData",
    "PCA_max=30",
    "",
    "primary transfer target=core_identity + lineage_priming_axis",
    "state transferred independently",
    "globins excluded from transfer HVGs",
    "no cross-library integration",
    "no Harmony",
    "no CCA",
    "condition not used",
    "native SoupX-corrected RNA",
    "cellranger_count not accessed",
    "",
    "gzip_integrity=PASS",
    "prediction_reread_count_order=PASS",
    "validation_rds_reread=PASS",
    "sha256_manifest=PASS"
  ),
  done_file
)

cat(
  "\n===== FOLD METRICS =====\n"
)

print(
  fold_metrics,
  row.names=FALSE
)

cat(
  "\n===== OVERALL =====\n"
)

print(
  overall,
  row.names=FALSE
)

cat(
  "\n===== CLASS THRESHOLD AUDIT =====\n"
)

print(
  threshold_audit,
  row.names=FALSE
)

cat(
  "\n===== PER CLASS =====\n"
)

print(
  exact_metrics,
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Progenitor leave-one-library-out validation completed\n"
)
