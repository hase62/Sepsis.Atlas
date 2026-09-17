#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 4L){
  stop(
    "Usage: script <root> <tag> <calibration_dir> <precision_dir>"
  )
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]
cal_dir <- normalizePath(args[[3]], mustWork=TRUE)
precision_dir <- normalizePath(args[[4]], mustWork=TRUE)

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

N_FULL <- 49585L
N_REFERENCE <- 26668L
N_QUERY <- N_FULL - N_REFERENCE
N_PCS <- 20L

CORE_THRESHOLD <- 0.70
POSITIVE_AXIS_THRESHOLD <- 0.70
NONE_AXIS_THRESHOLD <- 0.80
ACTIVATION_THRESHOLD <- 0.80
NONE_STATE_THRESHOLD <- 0.80

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "platelet_full_annotation_freeze_v1__",
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
  "PLATELET_FULL_ANNOTATION_FREEZE_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop("Already completed: ", done_file)
}

# ============================================================
# Helpers
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
    open="wt"
  )

  tryCatch(
    write.table(
      x,
      file=con,
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
    stop(
      "gzip integrity failure: ",
      path
    )
  }

  invisible(TRUE)
}

count_table <- function(x, name){

  z <- as.data.frame(
    table(x),
    stringsAsFactors=FALSE
  )

  names(z) <- c(
    name,
    "n_cells"
  )

  z
}

# ============================================================
# Required inputs
# ============================================================

cal_rds_file <- file.path(
  cal_dir,
  "platelet_discovery_freeze_transfer_calibration_v1.rds"
)

cal_done <- file.path(
  cal_dir,
  "PLATELET_DISCOVERY_FREEZE_TRANSFER_CALIBRATION_COMPLETE.ok"
)

precision_done <- file.path(
  precision_dir,
  "PLATELET_TRANSFER_PRECISION_AUDIT_COMPLETE.ok"
)

precision_file <- file.path(
  precision_dir,
  "platelet_transfer_class_precision_by_threshold_v1.tsv"
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
  cal_rds_file,
  cal_done,
  precision_done,
  precision_file,
  lib_file
)){
  if(!file.exists(f)){
    stop(
      "Missing required input: ",
      f
    )
  }
}

# ============================================================
# Calibration freeze
# ============================================================

cal <- readRDS(
  cal_rds_file
)

stopifnot(
  identical(
    cal$annotation_version,
    "Platelet_discovery_annotation_v1"
  ),
  cal$n_discovery_cells == N_REFERENCE
)

decision <- cal$decision
ref_ann <- cal$discovery_annotation

transfer_features <- as.character(
  cal$transfer_features
)

pca_features <- as.character(
  cal$pca_features
)

stopifnot(
  nrow(ref_ann) == N_REFERENCE,
  !anyDuplicated(
    ref_ann$global_cell
  ),
  length(transfer_features) == 3697L,
  length(pca_features) == 3697L
)

ref_ann$global_cell <- as.character(
  ref_ann$global_cell
)

ref_ann$project_id <- as.character(
  ref_ann$project_id
)

ref_ann$library_key <- as.character(
  ref_ann$library_key
)

# ============================================================
# Verify threshold evidence from step 30b
# ============================================================

prec <- read_tsv(
  precision_file
)

get_precision <- function(
  class_pattern,
  threshold
){

  d <- prec[
    grepl(
      class_pattern,
      prec$predicted_class,
      fixed=TRUE
    ) &
      abs(
        prec$threshold -
          threshold
      ) < 1e-10,
    ,
    drop=FALSE
  ]

  if(nrow(d) != 1L){
    stop(
      "Precision lookup failed: ",
      class_pattern,
      " @ ",
      threshold
    )
  }

  as.numeric(
    d$exact_precision[[1]]
  )
}

stopifnot(
  get_precision(
    "Deferred_non_platelet_myeloid_like",
    0.70
  ) >= 0.99,

  get_precision(
    "Deferred_non_platelet_T_NK_like",
    0.70
  ) >= 0.99,

  get_precision(
    "megakaryocytic_transcription_candidate",
    0.70
  ) >= 0.96,

  get_precision(
    "megakaryocytic_transcription_high",
    0.70
  ) >= 0.98,

  get_precision(
    "activation_degranulation_candidate",
    0.80
  ) >= 0.95
)

# IFN transfer is intentionally NOT promoted.
ifn80 <- prec[
  grepl(
    "interferon_stimulated_candidate",
    prec$predicted_class,
    fixed=TRUE
  ) &
    abs(
      prec$threshold -
        0.80
    ) < 1e-10,
  ,
  drop=FALSE
]

stopifnot(
  nrow(ifn80) == 1L,
  as.integer(
    ifn80$n_predicted_accepted[[1]]
  ) == 0L
)

# ============================================================
# Full Platelet/megakaryocyte metadata
# ============================================================

transfer_files <- list.files(
  transfer_dir,
  pattern="\\.tsv\\.gz$",
  full.names=TRUE
)

stopifnot(
  length(transfer_files) == 158L
)

meta_list <- list()

for(f in transfer_files){

  z <- read_tsv(f)

  req <- c(
    "global_cell",
    "original_cell",
    "project_id",
    "library_key",
    "condition_binary",
    "integration_compartment_primary_v1"
  )

  stopifnot(
    all(
      req %in%
        names(z)
    )
  )

  z <- z[
    as.character(
      z$integration_compartment_primary_v1
    ) ==
      "Platelet_megakaryocyte",
    req,
    drop=FALSE
  ]

  if(nrow(z)){
    meta_list[[
      length(meta_list)+1L
    ]] <- z
  }
}

full_meta <- do.call(
  rbind,
  meta_list
)

rm(meta_list)

for(nm in c(
  "global_cell",
  "original_cell",
  "project_id",
  "library_key",
  "condition_binary"
)){
  full_meta[[nm]] <- as.character(
    full_meta[[nm]]
  )
}

stopifnot(
  nrow(full_meta) == N_FULL,
  !anyDuplicated(
    full_meta$global_cell
  ),
  all(
    ref_ann$global_cell %in%
      full_meta$global_cell
  )
)

query_meta <- full_meta[
  !full_meta$global_cell %in%
    ref_ann$global_cell,
  ,
  drop=FALSE
]

stopifnot(
  nrow(query_meta) == N_QUERY,
  !anyDuplicated(
    query_meta$global_cell
  )
)

cat(
  "\nreference=",
  nrow(ref_ann),
  " query=",
  nrow(query_meta),
  " full=",
  nrow(full_meta),
  "\n",
  sep=""
)

# ============================================================
# Library lookup
# ============================================================

lib <- read_tsv(
  lib_file
)

stopifnot(
  all(
    c(
      "project_id",
      "library_key",
      "final_rds_resolved"
    ) %in%
      names(lib)
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
# Attach original cell IDs to reference annotation
# ============================================================

ri <- match(
  ref_ann$global_cell,
  full_meta$global_cell
)

stopifnot(
  !anyNA(ri)
)

ref_ann$original_cell <-
  full_meta$original_cell[ri]

stopifnot(
  ref_ann$project_id ==
    full_meta$project_id[ri],
  ref_ann$library_key ==
    full_meta$library_key[ri]
)

# ============================================================
# Read each native library ONCE.
#
# Only 3697 frozen transfer features are retained.
# No Cell Ranger output is accessed.
# ============================================================

ref_matrices <- list()
query_matrices <- list()

library_keys <- unique(
  paste(
    full_meta$project_id,
    full_meta$library_key,
    sep="|||"
  )
)

for(ii in seq_along(
  library_keys
)){

  key <- library_keys[[ii]]

  sp <- strsplit(
    key,
    "\\|\\|\\|"
  )[[1]]

  p <- sp[[1]]
  lk <- sp[[2]]

  r <- ref_ann[
    ref_ann$project_id == p &
      ref_ann$library_key == lk,
    ,
    drop=FALSE
  ]

  q <- query_meta[
    query_meta$project_id == p &
      query_meta$library_key == lk,
    ,
    drop=FALSE
  ]

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

  obj0 <- readRDS(
    source_rds
  )

  counts0 <- get_rna_counts(
    obj0
  )

  if(!all(
    transfer_features %in%
      rownames(counts0)
  )){
    stop(
      "Frozen transfer feature mismatch: ",
      key
    )
  }

  if(nrow(r)){

    cells <- r$original_cell

    if(!all(
      cells %in%
        colnames(counts0)
    )){
      stop(
        "Missing reference cells: ",
        key
      )
    }

    m <- counts0[
      transfer_features,
      cells,
      drop=FALSE
    ]

    colnames(m) <- r$global_cell

    ref_matrices[[
      length(ref_matrices)+1L
    ]] <- m
  }

  if(nrow(q)){

    cells <- q$original_cell

    if(!all(
      cells %in%
        colnames(counts0)
    )){
      stop(
        "Missing query cells: ",
        key
      )
    }

    m <- counts0[
      transfer_features,
      cells,
      drop=FALSE
    ]

    colnames(m) <- q$global_cell

    query_matrices[[
      length(query_matrices)+1L
    ]] <- m
  }

  rm(
    obj0,
    counts0
  )

  if(exists("m")){
    rm(m)
  }

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

ref_counts <- do.call(
  cbind,
  ref_matrices
)

query_counts <- do.call(
  cbind,
  query_matrices
)

rm(
  ref_matrices,
  query_matrices
)

invisible(gc())

# ============================================================
# Exact order / count checks
# ============================================================

ref_counts <- ref_counts[
  ,
  ref_ann$global_cell,
  drop=FALSE
]

query_counts <- query_counts[
  ,
  query_meta$global_cell,
  drop=FALSE
]

stopifnot(
  ncol(ref_counts) == N_REFERENCE,
  ncol(query_counts) == N_QUERY,

  identical(
    colnames(ref_counts),
    ref_ann$global_cell
  ),

  identical(
    colnames(query_counts),
    query_meta$global_cell
  )
)

# ============================================================
# Build native reference/query objects
# ============================================================

ref_md <- ref_ann
rownames(ref_md) <- ref_md$global_cell

query_md <- query_meta
rownames(query_md) <- query_md$global_cell

ref <- CreateSeuratObject(
  counts=ref_counts,
  meta.data=ref_md,
  min.cells=0,
  min.features=0
)

query <- CreateSeuratObject(
  counts=query_counts,
  meta.data=query_md,
  min.cells=0,
  min.features=0
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

# ============================================================
# Frozen reference PCA
# ============================================================

stopifnot(
  all(
    pca_features %in%
      rownames(ref)
  ),
  all(
    pca_features %in%
      rownames(query)
  )
)

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

# ============================================================
# Full query transfer
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

pi <- match(
  query_meta$global_cell,
  pred$global_cell
)

stopifnot(
  !anyNA(pi)
)

pred <- pred[
  pi,
  ,
  drop=FALSE
]

stopifnot(
  identical(
    pred$global_cell,
    query_meta$global_cell
  )
)

# ============================================================
# Decode raw transfer class
# ============================================================

class_map <- unique(
  decision[
    ,
    c(
      "transfer_class_v1",
      "platelet_core_v1",
      "platelet_transcriptional_axis_v1",
      "platelet_state_v1",
      "core_evidence_confidence_v1",
      "transcriptional_axis_evidence_confidence_v1",
      "state_evidence_confidence_v1"
    ),
    drop=FALSE
  ]
)

ci <- match(
  as.character(
    pred$predicted.id
  ),
  class_map$transfer_class_v1
)

stopifnot(
  !anyNA(ci)
)

score <- as.numeric(
  pred$prediction.score.max
)

raw_core <-
  class_map$platelet_core_v1[ci]

raw_axis <-
  class_map$platelet_transcriptional_axis_v1[ci]

raw_state <-
  class_map$platelet_state_v1[ci]

raw_core_bio_conf <-
  class_map$core_evidence_confidence_v1[ci]

raw_axis_bio_conf <-
  class_map$transcriptional_axis_evidence_confidence_v1[ci]

raw_state_bio_conf <-
  class_map$state_evidence_confidence_v1[ci]

# ============================================================
# Transfer confidence
# ============================================================

transfer_conf <- ifelse(
  score >= 0.80,
  "high",
  ifelse(
    score >= CORE_THRESHOLD,
    "medium",
    "unresolved"
  )
)

# ============================================================
# CORE
#
# Calibration precision at >=0.70:
#   myeloid 0.992
#   T/NK    1.000
#   Platelet 0.985
# ============================================================

final_core <- ifelse(
  score >= CORE_THRESHOLD,
  raw_core,
  "Platelet_transfer_unresolved"
)

final_scope <- ifelse(
  grepl(
    "^Deferred_non_platelet_",
    final_core
  ),
  "deferred_non_platelet",
  ifelse(
    final_core == "Platelet_like",
    "Platelet_megakaryocyte",
    "transfer_unresolved"
  )
)

final_core_bio_conf <- ifelse(
  score >= CORE_THRESHOLD,
  raw_core_bio_conf,
  "unresolved"
)

# ============================================================
# TRANSCRIPTIONAL AXIS
#
# Positive axis:
#   accept >=0.70
#
# "none":
#   only freeze >=0.80
#
# Lower-score absence is unresolved rather than negative.
# ============================================================

final_axis <- rep(
  "unresolved",
  length(score)
)

final_axis_bio_conf <- rep(
  "unresolved",
  length(score)
)

deferred_ok <-
  score >= CORE_THRESHOLD &
  grepl(
    "^Deferred_non_platelet_",
    raw_core
  )

final_axis[
  deferred_ok
] <- "not_applicable"

final_axis_bio_conf[
  deferred_ok
] <- "not_applicable"

platelet_ok <-
  score >= CORE_THRESHOLD &
  raw_core == "Platelet_like"

positive_axis <-
  raw_axis %in%
    c(
      "megakaryocytic_transcription_high",
      "megakaryocytic_transcription_candidate"
    )

ii <-
  platelet_ok &
  positive_axis &
  score >= POSITIVE_AXIS_THRESHOLD

final_axis[ii] <-
  raw_axis[ii]

final_axis_bio_conf[ii] <-
  raw_axis_bio_conf[ii]

ii <-
  platelet_ok &
  raw_axis == "none" &
  score >= NONE_AXIS_THRESHOLD

final_axis[ii] <- "none"

final_axis_bio_conf[ii] <-
  "not_applicable"

# ============================================================
# STATE
#
# activation:
#   accept >=0.80
#
# IFN:
#   NEVER promote in query.
#   Raw prediction is retained separately.
#
# none:
#   freeze only >=0.80.
# ============================================================

final_state <- rep(
  "unresolved",
  length(score)
)

final_state_bio_conf <- rep(
  "unresolved",
  length(score)
)

final_state[
  deferred_ok
] <- "not_applicable"

final_state_bio_conf[
  deferred_ok
] <- "not_applicable"

ii <-
  platelet_ok &
  raw_state ==
    "activation_degranulation_candidate" &
  score >= ACTIVATION_THRESHOLD

final_state[ii] <-
  "activation_degranulation_candidate"

final_state_bio_conf[ii] <-
  raw_state_bio_conf[ii]

ii <-
  platelet_ok &
  raw_state == "none" &
  score >= NONE_STATE_THRESHOLD

final_state[ii] <- "none"

final_state_bio_conf[ii] <-
  "not_applicable"

# Raw IFN prediction remains unresolved by design.
if(any(
  final_state[
    raw_state ==
      "interferon_stimulated_candidate"
  ] ==
    "interferon_stimulated_candidate"
)){
  stop(
    "Query IFN state was unexpectedly promoted"
  )
}

# ============================================================
# Query final annotation
# ============================================================

query_ann <- data.frame(
  global_cell=
    query_meta$global_cell,

  project_id=
    query_meta$project_id,

  library_key=
    query_meta$library_key,

  condition_binary=
    query_meta$condition_binary,

  platelet_core_v1=
    final_core,

  platelet_transcriptional_axis_v1=
    final_axis,

  platelet_state_v1=
    final_state,

  platelet_annotation_scope_v1=
    final_scope,

  platelet_core_evidence_confidence_v1=
    final_core_bio_conf,

  platelet_transcriptional_axis_evidence_confidence_v1=
    final_axis_bio_conf,

  platelet_state_evidence_confidence_v1=
    final_state_bio_conf,

  platelet_transfer_score_v1=
    score,

  platelet_transfer_confidence_v1=
    transfer_conf,

  platelet_raw_predicted_class_v1=
    as.character(
      pred$predicted.id
    ),

  platelet_raw_predicted_core_v1=
    raw_core,

  platelet_raw_predicted_transcriptional_axis_v1=
    raw_axis,

  platelet_raw_predicted_state_v1=
    raw_state,

  platelet_annotation_basis_v1=
    ifelse(
      score >= CORE_THRESHOLD,
      "full_query_supervised_transfer",
      "full_query_transfer_unresolved"
    ),

  stringsAsFactors=FALSE,
  check.names=FALSE
)

# ============================================================
# Reference annotation
#
# Discovery calls are NOT altered by transfer calibration.
# ============================================================

ref_final <- data.frame(
  global_cell=
    as.character(
      ref_ann$global_cell
    ),

  project_id=
    as.character(
      ref_ann$project_id
    ),

  library_key=
    as.character(
      ref_ann$library_key
    ),

  condition_binary=
    as.character(
      ref_ann$condition_binary
    ),

  platelet_core_v1=
    as.character(
      ref_ann$platelet_core_v1
    ),

  platelet_transcriptional_axis_v1=
    as.character(
      ref_ann$platelet_transcriptional_axis_v1
    ),

  platelet_state_v1=
    as.character(
      ref_ann$platelet_state_v1
    ),

  platelet_annotation_scope_v1=
    ifelse(
      grepl(
        "^Deferred_non_platelet_",
        ref_ann$platelet_core_v1
      ),
      "deferred_non_platelet",
      "Platelet_megakaryocyte"
    ),

  platelet_core_evidence_confidence_v1=
    as.character(
      ref_ann$platelet_core_evidence_confidence_v1
    ),

  platelet_transcriptional_axis_evidence_confidence_v1=
    as.character(
      ref_ann$platelet_transcriptional_axis_evidence_confidence_v1
    ),

  platelet_state_evidence_confidence_v1=
    as.character(
      ref_ann$platelet_state_evidence_confidence_v1
    ),

  platelet_transfer_score_v1=
    NA_real_,

  platelet_transfer_confidence_v1=
    "reference_discovery",

  platelet_raw_predicted_class_v1=
    NA_character_,

  platelet_raw_predicted_core_v1=
    NA_character_,

  platelet_raw_predicted_transcriptional_axis_v1=
    NA_character_,

  platelet_raw_predicted_state_v1=
    NA_character_,

  platelet_annotation_basis_v1=
    "project_balanced_native_discovery",

  stringsAsFactors=FALSE,
  check.names=FALSE
)

# ============================================================
# Full 49,585 annotation
# ============================================================

full_ann <- rbind(
  ref_final,
  query_ann
)

fi <- match(
  full_meta$global_cell,
  full_ann$global_cell
)

stopifnot(
  !anyNA(fi)
)

full_ann <- full_ann[
  fi,
  ,
  drop=FALSE
]

rownames(full_ann) <- NULL

stopifnot(
  nrow(full_ann) == N_FULL,
  !anyDuplicated(
    full_ann$global_cell
  ),
  identical(
    full_ann$global_cell,
    full_meta$global_cell
  )
)

# ============================================================
# Semantic invariants
# ============================================================

deferred <-
  full_ann$platelet_annotation_scope_v1 ==
    "deferred_non_platelet"

stopifnot(
  all(
    full_ann$platelet_transcriptional_axis_v1[
      deferred
    ] ==
      "not_applicable"
  ),

  all(
    full_ann$platelet_state_v1[
      deferred
    ] ==
      "not_applicable"
  )
)

# Discovery IFN must remain exactly c7-derived annotation.
n_ref_ifn <- sum(
  ref_final$platelet_state_v1 ==
    "interferon_stimulated_candidate"
)

stopifnot(
  n_ref_ifn == 850L
)

# No query IFN promotion.
stopifnot(
  !any(
    query_ann$platelet_state_v1 ==
      "interferon_stimulated_candidate"
  )
)

# ============================================================
# Frozen policy table
# ============================================================

policy <- data.frame(
  dimension=c(
    "core",
    "positive_megakaryocytic_axis",
    "axis_none",
    "activation_state",
    "state_none",
    "IFN_state"
  ),

  threshold=c(
    CORE_THRESHOLD,
    POSITIVE_AXIS_THRESHOLD,
    NONE_AXIS_THRESHOLD,
    ACTIVATION_THRESHOLD,
    NONE_STATE_THRESHOLD,
    NA_real_
  ),

  query_policy=c(
    "accept predicted core",
    "accept positive axis",
    "freeze none only above threshold",
    "accept activation candidate",
    "freeze none only above threshold",
    "do not promote from transfer"
  ),

  calibration_basis=c(
    "core precision >=0.985 at score >=0.70",
    "positive-axis precision >=0.969 at score >=0.70",
    "axis-none precision 0.971 at score >=0.80",
    "activation precision 0.962 at score >=0.80",
    "state-none precision 0.949 at score >=0.80",
    "no IFN predictions accepted at score >=0.80 in validation"
  ),

  stringsAsFactors=FALSE
)

# ============================================================
# Summaries
# ============================================================

core_counts <- count_table(
  full_ann$platelet_core_v1,
  "platelet_core_v1"
)

axis_counts <- count_table(
  full_ann$platelet_transcriptional_axis_v1,
  "platelet_transcriptional_axis_v1"
)

state_counts <- count_table(
  full_ann$platelet_state_v1,
  "platelet_state_v1"
)

transfer_counts <- count_table(
  full_ann$platelet_transfer_confidence_v1,
  "platelet_transfer_confidence_v1"
)

joint_counts <- aggregate(
  list(
    n_cells=
      rep(
        1L,
        nrow(full_ann)
      )
  ),
  by=list(
    platelet_core_v1=
      full_ann$platelet_core_v1,

    platelet_transcriptional_axis_v1=
      full_ann$platelet_transcriptional_axis_v1,

    platelet_state_v1=
      full_ann$platelet_state_v1,

    annotation_basis=
      full_ann$platelet_annotation_basis_v1
  ),
  FUN=sum
)

stopifnot(
  sum(core_counts$n_cells) == N_FULL,
  sum(axis_counts$n_cells) == N_FULL,
  sum(state_counts$n_cells) == N_FULL,
  sum(transfer_counts$n_cells) == N_FULL,
  sum(joint_counts$n_cells) == N_FULL
)

# Query raw-class distribution.
raw_query_counts <- count_table(
  query_ann$platelet_raw_predicted_class_v1,
  "platelet_raw_predicted_class_v1"
)

# ============================================================
# Output files
# ============================================================

annotation_file <- file.path(
  out_dir,
  paste0(
    "platelet_full_cell_annotation_v1__",
    tag,
    ".tsv.gz"
  )
)

raw_pred_file <- file.path(
  out_dir,
  paste0(
    "platelet_full_query_raw_transfer_predictions_v1__",
    tag,
    ".tsv.gz"
  )
)

freeze_rds <- file.path(
  out_dir,
  paste0(
    "platelet_full_annotation_freeze_v1__",
    tag,
    ".rds"
  )
)

policy_file <- file.path(
  out_dir,
  "platelet_full_transfer_policy_v1.tsv"
)

decision_file <- file.path(
  out_dir,
  "platelet_r0p4_annotation_decision_v1.tsv"
)

core_file <- file.path(
  out_dir,
  "platelet_full_core_counts_v1.tsv"
)

axis_file <- file.path(
  out_dir,
  "platelet_full_transcriptional_axis_counts_v1.tsv"
)

state_file <- file.path(
  out_dir,
  "platelet_full_state_counts_v1.tsv"
)

transfer_file <- file.path(
  out_dir,
  "platelet_full_transfer_confidence_counts_v1.tsv"
)

joint_file <- file.path(
  out_dir,
  "platelet_full_core_axis_state_counts_v1.tsv"
)

raw_query_count_file <- file.path(
  out_dir,
  "platelet_full_query_raw_class_counts_v1.tsv"
)

# ============================================================
# Write TSVs
# ============================================================

write_gz_tsv(
  full_ann,
  annotation_file
)

pred_out <- pred[
  ,
  c(
    setdiff(
      names(pred),
      "global_cell"
    ),
    "global_cell"
  ),
  drop=FALSE
]

write_gz_tsv(
  pred_out,
  raw_pred_file
)

write.table(
  policy,
  policy_file,
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  decision,
  decision_file,
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  core_counts,
  core_file,
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  axis_counts,
  axis_file,
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  state_counts,
  state_file,
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  transfer_counts,
  transfer_file,
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  joint_counts,
  joint_file,
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  raw_query_counts,
  raw_query_count_file,
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Freeze RDS
# ============================================================

freeze <- list(

  annotation_version=
    "Platelet_megakaryocyte_annotation_v1",

  n_cells=
    N_FULL,

  n_discovery_reference_cells=
    N_REFERENCE,

  n_full_query_cells=
    N_QUERY,

  annotation=
    full_ann,

  discovery_decision=
    decision,

  transfer_policy=
    policy,

  transfer_features=
    transfer_features,

  pca_features=
    pca_features,

  calibration_directory=
    cal_dir,

  precision_audit_directory=
    precision_dir,

  rules=list(

    discovery_backbone=
      "balanced RPCA r0.4",

    biological_marker_evidence=
      "native library -> project-balanced consensus markers",

    core_axis_state_separate=
      TRUE,

    query_mapping=
      "Seurat FindTransferAnchors pcaproject + TransferData",

    core_acceptance_threshold=
      CORE_THRESHOLD,

    positive_axis_acceptance_threshold=
      POSITIVE_AXIS_THRESHOLD,

    axis_none_threshold=
      NONE_AXIS_THRESHOLD,

    activation_state_threshold=
      ACTIVATION_THRESHOLD,

    state_none_threshold=
      NONE_STATE_THRESHOLD,

    IFN_query_transfer=
      "not promoted",

    IFN_discovery_reference=
      paste0(
        "retained as interferon_stimulated_candidate ",
        "from native discovery evidence"
      ),

    low_score_policy=
      paste0(
        "score <0.70: transfer unresolved; ",
        "raw prediction retained"
      ),

    none_definition=
      paste0(
        "no frozen positive program/state claim; ",
        "not proof of biological absence"
      )
  )
)

tmp_rds <- paste0(
  freeze_rds,
  ".tmp.",
  Sys.getpid()
)

saveRDS(
  freeze,
  tmp_rds,
  compress=TRUE
)

if(!file.rename(
  tmp_rds,
  freeze_rds
)){
  stop(
    "Atomic freeze RDS rename failed"
  )
}

# ============================================================
# Strong re-read checks
# ============================================================

check_tsv <- read_tsv(
  annotation_file
)

stopifnot(
  nrow(check_tsv) == N_FULL,
  !anyDuplicated(
    check_tsv$global_cell
  ),
  identical(
    as.character(
      check_tsv$global_cell
    ),
    as.character(
      full_ann$global_cell
    )
  ),
  identical(
    as.character(
      check_tsv$platelet_core_v1
    ),
    as.character(
      full_ann$platelet_core_v1
    )
  ),
  identical(
    as.character(
      check_tsv$platelet_transcriptional_axis_v1
    ),
    as.character(
      full_ann$platelet_transcriptional_axis_v1
    )
  ),
  identical(
    as.character(
      check_tsv$platelet_state_v1
    ),
    as.character(
      full_ann$platelet_state_v1
    )
  )
)

check_pred <- read_tsv(
  raw_pred_file
)

stopifnot(
  nrow(check_pred) == N_QUERY,
  !anyDuplicated(
    check_pred$global_cell
  )
)

check_rds <- readRDS(
  freeze_rds
)

stopifnot(
  identical(
    check_rds$annotation_version,
    "Platelet_megakaryocyte_annotation_v1"
  ),
  check_rds$n_cells == N_FULL,
  nrow(
    check_rds$annotation
  ) == N_FULL,
  identical(
    check_rds$annotation$global_cell,
    full_ann$global_cell
  )
)

rm(
  check_tsv,
  check_pred,
  check_rds
)

invisible(gc())

# ============================================================
# SHA256 provenance
# ============================================================

manifest_files <- c(
  annotation_file,
  raw_pred_file,
  freeze_rds,
  policy_file,
  decision_file,
  core_file,
  axis_file,
  state_file,
  transfer_file,
  joint_file,
  raw_query_count_file
)

sha_lines <- unlist(
  lapply(
    manifest_files,
    function(f){

      out <- system2(
        "sha256sum",
        f,
        stdout=TRUE,
        stderr=TRUE
      )

      if(length(out) != 1L){
        stop(
          "sha256sum failed: ",
          f
        )
      }

      out
    }
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

writeLines(
  sha_lines,
  sha_file
)

# ============================================================
# Completion marker LAST
# ============================================================

writeLines(
  c(
    "PASS",
    "Platelet/megakaryocyte full annotation freeze v1",
    "annotation_version=Platelet_megakaryocyte_annotation_v1",
    "n_cells=49585",
    "n_discovery_reference_cells=26668",
    "n_full_query_cells=22917",
    "",
    "discovery_backbone=balanced RPCA r0.4",
    "biological evidence=native project-balanced markers",
    "full query mapping=Seurat pcaproject + TransferData",
    "",
    "core acceptance threshold=0.70",
    "positive megakaryocytic-axis threshold=0.70",
    "axis none threshold=0.80",
    "activation state threshold=0.80",
    "state none threshold=0.80",
    "IFN state NOT promoted by full-query transfer",
    "",
    "discovery IFN candidate retained from native evidence",
    "raw query predictions retained",
    "core / transcriptional axis / state remain separate",
    "",
    "gzip_integrity=PASS",
    "full_annotation_re_read_count_order=PASS",
    "raw_prediction_re_read=PASS",
    "freeze_rds_re_read=PASS",
    "sha256_manifest=PASS",
    "cellranger_count not accessed",
    paste0(
      "freeze_rds=",
      freeze_rds
    ),
    paste0(
      "freeze_tsv=",
      annotation_file
    )
  ),
  done_file
)

# ============================================================
# Console
# ============================================================

cat(
  "\n===== CORE =====\n"
)

print(
  core_counts,
  row.names=FALSE
)

cat(
  "\n===== TRANSCRIPTIONAL AXIS =====\n"
)

print(
  axis_counts,
  row.names=FALSE
)

cat(
  "\n===== STATE =====\n"
)

print(
  state_counts,
  row.names=FALSE
)

cat(
  "\n===== TRANSFER CONFIDENCE =====\n"
)

print(
  transfer_counts,
  row.names=FALSE
)

cat(
  "\n===== RAW QUERY CLASS =====\n"
)

print(
  raw_query_counts,
  row.names=FALSE
)

cat(
  "\n===== CORE x AXIS x STATE x BASIS =====\n"
)

print(
  joint_counts[
    order(
      joint_counts$platelet_core_v1,
      joint_counts$platelet_transcriptional_axis_v1,
      joint_counts$platelet_state_v1,
      joint_counts$annotation_basis
    ),
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Platelet full annotation freeze completed\n"
)
