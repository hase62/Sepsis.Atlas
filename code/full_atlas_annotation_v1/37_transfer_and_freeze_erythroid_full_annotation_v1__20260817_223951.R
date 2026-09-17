#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 4L){
  stop("Usage: script <root> <tag> <freeze_dir> <calibration_dir>")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]
freeze_dir <- normalizePath(args[[3]], mustWork=TRUE)
cal_dir <- normalizePath(args[[4]], mustWork=TRUE)

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

N_DISCOVERY <- 8692L
N_FULL <- 12177L
N_QUERY <- N_FULL - N_DISCOVERY

# Frozen transfer policy.
CORE_ERY_THRESHOLD <- 0.70
CORE_MIXED_THRESHOLD <- 0.70

# Myeloid query predictions are deliberately NOT promoted.
# Held-out precision remained 0.913 / 0.926 / 0.939
# at score thresholds 0.70 / 0.80 / 0.90.

MAT_REGULATORY_THRESHOLD <- 0.70
MAT_LATE_THRESHOLD <- 0.90
MAT_TERMINAL_THRESHOLD <- 0.90
STATE_NONE_THRESHOLD <- 0.80

TRANSFER_HIGH_THRESHOLD <- 0.90
TRANSFER_MEDIUM_THRESHOLD <- 0.70

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "erythroid_full_annotation_freeze_v1__",
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
  "ERYTHROID_FULL_ANNOTATION_FREEZE_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop("Already completed: ", done_file)
}

read_tsv <- function(path){

  if(!file.exists(path)){
    stop("Missing: ", path)
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

count1 <- function(x, name){

  z <- as.data.frame(
    table(x),
    stringsAsFactors=FALSE
  )

  names(z) <- c(name, "n_cells")
  z
}

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

# ============================================================
# Frozen discovery annotation
# ============================================================

freeze_files <- list.files(
  freeze_dir,
  pattern="^erythroid_discovery_annotation_freeze_v1__.*\\.rds$",
  full.names=TRUE
)

if(length(freeze_files) != 1L){
  stop("Expected exactly one discovery freeze RDS")
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
  "erythroid_evidence_confidence_v1",
  "erythroid_evidence_flag_v1",
  "erythroid_transfer_class_v1",
  "rationale"
)

stopifnot(all(required_ann %in% names(ann)))

for(nm in required_ann){
  ann[[nm]] <- as.character(ann[[nm]])
}

# ============================================================
# Frozen calibration
# ============================================================

cal_files <- list.files(
  cal_dir,
  pattern="^erythroid_discovery_transfer_calibration_v1__.*\\.rds$",
  full.names=TRUE
)

if(length(cal_files) != 1L){
  stop("Expected exactly one calibration RDS")
}

cal <- readRDS(cal_files[[1]])

stopifnot(
  cal$n_discovery == N_DISCOVERY,
  cal$n_validation > 0L
)

transfer_features <- as.character(
  cal$transfer_features
)

pca_features <- as.character(
  cal$pca_features
)

npcs <- as.integer(
  cal$npcs
)

stopifnot(
  length(transfer_features) > 1000L,
  identical(transfer_features, pca_features),
  npcs >= 10L
)

cat(
  "Frozen transfer features = ",
  length(transfer_features),
  "\n",
  sep=""
)

cat(
  "Frozen PCs = ",
  npcs,
  "\n",
  sep=""
)

# ============================================================
# Full Erythroid metadata
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
  stop("No annotation-transfer files found")
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

full_meta <- do.call(
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
  full_meta[[nm]] <- as.character(full_meta[[nm]])
}

stopifnot(
  nrow(full_meta) == N_FULL,
  !anyDuplicated(full_meta$global_cell)
)

discovery_match <- match(
  ann$global_cell,
  full_meta$global_cell
)

stopifnot(!anyNA(discovery_match))

stopifnot(
  ann$project_id ==
    full_meta$project_id[discovery_match],

  ann$library_key ==
    full_meta$library_key[discovery_match]
)

query_mask <- !(
  full_meta$global_cell %in%
    ann$global_cell
)

query_meta <- full_meta[
  query_mask,
  ,
  drop=FALSE
]

stopifnot(
  nrow(query_meta) == N_QUERY,
  !anyDuplicated(query_meta$global_cell)
)

cat(
  "Full cells = ",
  N_FULL,
  "\nDiscovery reference = ",
  N_DISCOVERY,
  "\nFull query = ",
  N_QUERY,
  "\n",
  sep=""
)

# ============================================================
# Resolved library lookup
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
    ) %in%
      names(lib)
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
# Read each native SoupX RNA library ONCE
# Only frozen transfer features are retained.
# ============================================================

library_keys <- unique(
  paste(
    full_meta$project_id,
    full_meta$library_key,
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

  d <- full_meta[
    full_meta$project_id == p &
      full_meta$library_key == lk,
    ,
    drop=FALSE
  ]

  li <- match(
    key,
    lib$key
  )

  if(is.na(li)){
    stop("Library lookup failed: ", key)
  }

  source_rds <- as.character(
    lib$final_rds_resolved[[li]]
  )

  if(!file.exists(source_rds)){
    stop("Missing native RDS: ", source_rds)
  }

  obj0 <- readRDS(source_rds)

  counts0 <- get_rna_counts(
    obj0
  )

  missing_features <- setdiff(
    transfer_features,
    rownames(counts0)
  )

  if(length(missing_features)){
    stop(
      "Frozen transfer features missing in ",
      key,
      ": ",
      paste(head(missing_features, 20L), collapse=",")
    )
  }

  if(!all(
    d$original_cell %in%
      colnames(counts0)
  )){
    stop("Missing full Erythroid cells in ", key)
  }

  m <- counts0[
    transfer_features,
    d$original_cell,
    drop=FALSE
  ]

  colnames(m) <- d$global_cell

  mat_list[[key]] <- m

  rm(
    obj0,
    counts0,
    m
  )

  invisible(gc())

  cat(
    sprintf(
      "READ %3d/%3d %s cells=%d\n",
      ii,
      length(library_keys),
      key,
      nrow(d)
    )
  )
}

counts <- do.call(
  cbind,
  mat_list
)

rm(mat_list)
invisible(gc())

stopifnot(
  ncol(counts) == N_FULL,
  !anyDuplicated(colnames(counts))
)

oi <- match(
  full_meta$global_cell,
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
    full_meta$global_cell
  )
)

# ============================================================
# Reference = ALL frozen discovery cells
# Query = remaining 3,485 cells
# ============================================================

ref <- CreateSeuratObject(
  counts=counts[
    ,
    ann$global_cell,
    drop=FALSE
  ],
  min.cells=0,
  min.features=0,
  project="erythroid_full_reference"
)

query <- CreateSeuratObject(
  counts=counts[
    ,
    query_meta$global_cell,
    drop=FALSE
  ],
  min.cells=0,
  min.features=0,
  project="erythroid_full_query"
)

rm(counts)
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
  features=pca_features,
  verbose=FALSE
)

ref <- RunPCA(
  ref,
  features=pca_features,
  npcs=npcs,
  verbose=FALSE
)

anchors <- FindTransferAnchors(
  reference=ref,
  query=query,
  normalization.method="LogNormalize",
  reduction="pcaproject",
  reference.reduction="pca",
  features=pca_features,
  dims=seq_len(npcs),
  k.anchor=5,
  k.filter=50,
  verbose=TRUE
)

ref_labels <- ann$erythroid_transfer_class_v1

names(ref_labels) <- ann$global_cell

pred <- TransferData(
  anchorset=anchors,
  refdata=ref_labels,
  dims=seq_len(npcs),
  k.weight=50,
  verbose=TRUE
)

pred <- as.data.frame(pred)

stopifnot(
  nrow(pred) == N_QUERY,
  all(rownames(pred) %in% query_meta$global_cell)
)

pi <- match(
  query_meta$global_cell,
  rownames(pred)
)

stopifnot(!anyNA(pi))

pred <- pred[
  pi,
  ,
  drop=FALSE
]

stopifnot(
  identical(
    rownames(pred),
    query_meta$global_cell
  )
)

raw_class <- as.character(
  pred$predicted.id
)

raw_score <- as.numeric(
  pred$prediction.score.max
)

raw_dec <- decode_class(
  raw_class
)

# ============================================================
# Frozen full-query policy
# ============================================================

canonical_core <- rep(
  "Erythroid_transfer_unresolved",
  N_QUERY
)

# Erythroid core is very high precision at >=0.70.
take_ery <-
  raw_dec$core == "Erythroid_like" &
  raw_score >= CORE_ERY_THRESHOLD

canonical_core[take_ery] <-
  "Erythroid_like"

# Mixed lineage transferred extremely precisely.
take_mixed <-
  raw_dec$core == "Deferred_mixed_lineage_like" &
  raw_score >= CORE_MIXED_THRESHOLD

canonical_core[take_mixed] <-
  "Deferred_mixed_lineage_like"

# IMPORTANT:
# Deferred_non_erythroid_myeloid_like is NEVER promoted in query.
# Raw prediction is retained for later global reconciliation.

canonical_maturation <- rep(
  "unresolved",
  N_QUERY
)

# Deferred mixed-lineage cells have no Erythroid maturation claim.
canonical_maturation[
  canonical_core == "Deferred_mixed_lineage_like"
] <- "not_applicable"

ery_ok <- canonical_core == "Erythroid_like"

take_reg <-
  ery_ok &
  raw_dec$maturation == "regulatory_rich_erythroid" &
  raw_score >= MAT_REGULATORY_THRESHOLD

canonical_maturation[take_reg] <-
  "regulatory_rich_erythroid"

take_late <-
  ery_ok &
  raw_dec$maturation ==
    "late_maturation_globin_dominant_candidate" &
  raw_score >= MAT_LATE_THRESHOLD

canonical_maturation[take_late] <-
  "late_maturation_globin_dominant_candidate"

take_terminal <-
  ery_ok &
  raw_dec$maturation ==
    "terminal_hemoglobinization_high" &
  raw_score >= MAT_TERMINAL_THRESHOLD

canonical_maturation[take_terminal] <-
  "terminal_hemoglobinization_high"

# Raw unresolved remains unresolved by construction.

canonical_state <- rep(
  "unresolved",
  N_QUERY
)

canonical_state[
  canonical_core == "Deferred_mixed_lineage_like"
] <- "not_applicable"

take_state_none <-
  canonical_core == "Erythroid_like" &
  raw_score >= STATE_NONE_THRESHOLD

canonical_state[take_state_none] <- "none"

transfer_confidence <- ifelse(
  raw_score >= TRANSFER_HIGH_THRESHOLD,
  "high",
  ifelse(
    raw_score >= TRANSFER_MEDIUM_THRESHOLD,
    "medium",
    "unresolved"
  )
)

evidence_flag <- rep(
  "query_supervised_transfer",
  N_QUERY
)

evidence_flag[
  raw_dec$core ==
    "Deferred_non_erythroid_myeloid_like"
] <- "raw_myeloid_prediction_not_promoted"

annotation_basis <- ifelse(
  canonical_core == "Erythroid_transfer_unresolved",
  "full_query_transfer_unresolved",
  "full_query_supervised_transfer"
)

query_out <- data.frame(
  global_cell=query_meta$global_cell,
  original_cell=query_meta$original_cell,
  project_id=query_meta$project_id,
  library_key=query_meta$library_key,

  erythroid_discovery_cluster_r0p4=NA_character_,

  erythroid_core_v1=canonical_core,
  erythroid_maturation_axis_v1=canonical_maturation,
  erythroid_state_v1=canonical_state,

  erythroid_evidence_confidence_v1=
    "not_assessed_native_query",

  erythroid_evidence_flag_v1=
    evidence_flag,

  rationale=
    ifelse(
      evidence_flag ==
        "raw_myeloid_prediction_not_promoted",
      "raw myeloid transfer prediction retained but not promoted by frozen calibration policy",
      "full-query supervised transfer under frozen class-specific thresholds"
    ),

  annotation_basis=
    annotation_basis,

  erythroid_transfer_confidence_v1=
    transfer_confidence,

  erythroid_raw_predicted_class_v1=
    raw_class,

  erythroid_raw_prediction_score_v1=
    raw_score,

  stringsAsFactors=FALSE
)

query_out$erythroid_annotation_class_v1 <- paste(
  query_out$erythroid_core_v1,
  query_out$erythroid_maturation_axis_v1,
  query_out$erythroid_state_v1,
  sep="|||"
)

# ============================================================
# Discovery annotation, unchanged
# ============================================================

disc_meta_idx <- match(
  ann$global_cell,
  full_meta$global_cell
)

disc_out <- data.frame(
  global_cell=ann$global_cell,
  original_cell=full_meta$original_cell[disc_meta_idx],
  project_id=ann$project_id,
  library_key=ann$library_key,

  erythroid_discovery_cluster_r0p4=
    ann$cluster_r0p4,

  erythroid_core_v1=
    ann$erythroid_core_v1,

  erythroid_maturation_axis_v1=
    ann$erythroid_maturation_axis_v1,

  erythroid_state_v1=
    ann$erythroid_state_v1,

  erythroid_evidence_confidence_v1=
    ann$erythroid_evidence_confidence_v1,

  erythroid_evidence_flag_v1=
    ann$erythroid_evidence_flag_v1,

  rationale=
    ann$rationale,

  annotation_basis=
    "project_balanced_native_discovery",

  erythroid_transfer_confidence_v1=
    "reference_discovery",

  erythroid_raw_predicted_class_v1=
    NA_character_,

  erythroid_raw_prediction_score_v1=
    NA_real_,

  stringsAsFactors=FALSE
)

disc_out$erythroid_annotation_class_v1 <- paste(
  disc_out$erythroid_core_v1,
  disc_out$erythroid_maturation_axis_v1,
  disc_out$erythroid_state_v1,
  sep="|||"
)

# ============================================================
# Combine in full metadata order
# ============================================================

combined <- rbind(
  disc_out,
  query_out
)

ci <- match(
  full_meta$global_cell,
  combined$global_cell
)

stopifnot(
  !anyNA(ci),
  !anyDuplicated(combined$global_cell)
)

combined <- combined[
  ci,
  ,
  drop=FALSE
]

rownames(combined) <- NULL

stopifnot(
  nrow(combined) == N_FULL,
  identical(
    combined$global_cell,
    full_meta$global_cell
  )
)

# ============================================================
# Semantic invariants
# ============================================================

# Discovery must be bit-for-bit identical in semantic labels.
di <- match(
  ann$global_cell,
  combined$global_cell
)

stopifnot(
  identical(
    combined$erythroid_core_v1[di],
    ann$erythroid_core_v1
  ),

  identical(
    combined$erythroid_maturation_axis_v1[di],
    ann$erythroid_maturation_axis_v1
  ),

  identical(
    combined$erythroid_state_v1[di],
    ann$erythroid_state_v1
  ),

  identical(
    combined$erythroid_evidence_flag_v1[di],
    ann$erythroid_evidence_flag_v1
  )
)

# Discovery myeloid count must remain exactly 352.
stopifnot(
  sum(
    combined$erythroid_core_v1[di] ==
      "Deferred_non_erythroid_myeloid_like"
  ) == 352L
)

# Query myeloid raw predictions are NEVER promoted.
qi <- match(
  query_out$global_cell,
  combined$global_cell
)

qraw <- raw_dec$core

stopifnot(
  all(
    combined$erythroid_core_v1[qi][
      qraw == "Deferred_non_erythroid_myeloid_like"
    ] == "Erythroid_transfer_unresolved"
  )
)

# Deferred canonical populations carry no maturation/state claim.
deferred <- grepl(
  "^Deferred_",
  combined$erythroid_core_v1
)

stopifnot(
  all(
    combined$erythroid_maturation_axis_v1[deferred] ==
      "not_applicable"
  ),

  all(
    combined$erythroid_state_v1[deferred] ==
      "not_applicable"
  )
)

# Transfer unresolved must remain unresolved in both dimensions.
tu <-
  combined$erythroid_core_v1 ==
  "Erythroid_transfer_unresolved"

stopifnot(
  all(
    combined$erythroid_maturation_axis_v1[tu] ==
      "unresolved"
  ),

  all(
    combined$erythroid_state_v1[tu] ==
      "unresolved"
  )
)

# Erythroid cells cannot have not_applicable maturation/state.
ery <-
  combined$erythroid_core_v1 ==
  "Erythroid_like"

stopifnot(
  all(
    combined$erythroid_maturation_axis_v1[ery] !=
      "not_applicable"
  ),

  all(
    combined$erythroid_state_v1[ery] !=
      "not_applicable"
  )
)

# ============================================================
# Policy table
# ============================================================

policy <- data.frame(
  dimension=c(
    "core",
    "core",
    "core",
    "maturation",
    "maturation",
    "maturation",
    "maturation",
    "state"
  ),

  raw_label=c(
    "Erythroid_like",
    "Deferred_mixed_lineage_like",
    "Deferred_non_erythroid_myeloid_like",
    "regulatory_rich_erythroid",
    "late_maturation_globin_dominant_candidate",
    "terminal_hemoglobinization_high",
    "unresolved",
    "none"
  ),

  threshold=c(
    0.70,
    0.70,
    NA,
    0.70,
    0.90,
    0.90,
    NA,
    0.80
  ),

  action=c(
    "accept",
    "accept",
    "never_promote_query; retain raw prediction",
    "accept_if_threshold_met",
    "accept_if_threshold_met",
    "accept_if_threshold_met",
    "remain_unresolved",
    "accept_if_threshold_met"
  ),

  heldout_precision_reference=c(
    "0.9985 at 0.70",
    "1.0000 at 0.70",
    "0.9394 even at 0.90",
    "0.9753 at 0.70",
    "0.9442 at 0.90",
    "0.9176 at 0.90",
    "0.9583 predicted-class precision at 0.70; actual-class recall poor",
    "1.0000 at 0.80"
  ),

  stringsAsFactors=FALSE
)

# ============================================================
# Summaries
# ============================================================

core_counts <- count1(
  combined$erythroid_core_v1,
  "erythroid_core_v1"
)

maturation_counts <- count1(
  combined$erythroid_maturation_axis_v1,
  "erythroid_maturation_axis_v1"
)

state_counts <- count1(
  combined$erythroid_state_v1,
  "erythroid_state_v1"
)

transfer_confidence_counts <- count1(
  combined$erythroid_transfer_confidence_v1,
  "erythroid_transfer_confidence_v1"
)

raw_query_class_counts <- count1(
  query_out$erythroid_raw_predicted_class_v1,
  "erythroid_raw_predicted_class_v1"
)

joint <- as.data.frame(
  table(
    combined$erythroid_core_v1,
    combined$erythroid_maturation_axis_v1,
    combined$erythroid_state_v1,
    combined$annotation_basis
  ),
  stringsAsFactors=FALSE
)

names(joint) <- c(
  "erythroid_core_v1",
  "erythroid_maturation_axis_v1",
  "erythroid_state_v1",
  "annotation_basis",
  "n_cells"
)

joint <- joint[
  joint$n_cells > 0L,
  ,
  drop=FALSE
]

stopifnot(
  sum(core_counts$n_cells) == N_FULL,
  sum(maturation_counts$n_cells) == N_FULL,
  sum(state_counts$n_cells) == N_FULL,
  sum(transfer_confidence_counts$n_cells) == N_FULL,
  sum(raw_query_class_counts$n_cells) == N_QUERY,
  sum(joint$n_cells) == N_FULL
)

# ============================================================
# Outputs
# ============================================================

annotation_file <- file.path(
  out_dir,
  paste0(
    "erythroid_full_cell_annotation_v1__",
    tag,
    ".tsv.gz"
  )
)

raw_file <- file.path(
  out_dir,
  paste0(
    "erythroid_full_query_raw_transfer_predictions_v1__",
    tag,
    ".tsv.gz"
  )
)

write_gz_tsv(
  combined,
  annotation_file
)

write_gz_tsv(
  query_out,
  raw_file
)

write.table(
  policy,
  file.path(
    out_dir,
    "erythroid_full_transfer_policy_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  freeze$decision,
  file.path(
    out_dir,
    "erythroid_r0p4_annotation_decision_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  core_counts,
  file.path(
    out_dir,
    "erythroid_full_core_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  maturation_counts,
  file.path(
    out_dir,
    "erythroid_full_maturation_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  state_counts,
  file.path(
    out_dir,
    "erythroid_full_state_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  transfer_confidence_counts,
  file.path(
    out_dir,
    "erythroid_full_transfer_confidence_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  raw_query_class_counts,
  file.path(
    out_dir,
    "erythroid_full_query_raw_class_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  joint,
  file.path(
    out_dir,
    "erythroid_full_core_maturation_state_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Freeze RDS
# ============================================================

freeze_full <- list(
  annotation_version=
    "Erythroid_annotation_v1",

  n_cells=N_FULL,
  n_discovery_reference_cells=N_DISCOVERY,
  n_full_query_cells=N_QUERY,

  annotation=combined,
  query_raw_transfer=query_out,
  transfer_policy=policy,

  transfer_features=transfer_features,
  pca_features=pca_features,
  npcs=npcs,

  source_discovery_freeze_dir=freeze_dir,
  source_discovery_freeze_rds=freeze_files[[1]],
  source_calibration_dir=cal_dir,
  source_calibration_rds=cal_files[[1]]
)

rds_file <- file.path(
  out_dir,
  paste0(
    "erythroid_full_annotation_freeze_v1__",
    tag,
    ".rds"
  )
)

tmp_rds <- paste0(
  rds_file,
  ".tmp"
)

saveRDS(
  freeze_full,
  tmp_rds,
  compress=TRUE
)

if(!file.rename(
  tmp_rds,
  rds_file
)){
  stop("Atomic RDS rename failed")
}

# ============================================================
# Re-read integrity
# ============================================================

chk_ann <- read_tsv(
  annotation_file
)

stopifnot(
  nrow(chk_ann) == N_FULL,
  identical(
    as.character(chk_ann$global_cell),
    as.character(combined$global_cell)
  ),
  identical(
    as.character(chk_ann$erythroid_core_v1),
    as.character(combined$erythroid_core_v1)
  ),
  identical(
    as.character(chk_ann$erythroid_maturation_axis_v1),
    as.character(combined$erythroid_maturation_axis_v1)
  ),
  identical(
    as.character(chk_ann$erythroid_state_v1),
    as.character(combined$erythroid_state_v1)
  )
)

chk_raw <- read_tsv(
  raw_file
)

stopifnot(
  nrow(chk_raw) == N_QUERY,
  identical(
    as.character(chk_raw$global_cell),
    as.character(query_out$global_cell)
  )
)

chk_rds <- readRDS(
  rds_file
)

stopifnot(
  chk_rds$n_cells == N_FULL,
  chk_rds$n_discovery_reference_cells == N_DISCOVERY,
  chk_rds$n_full_query_cells == N_QUERY,
  nrow(chk_rds$annotation) == N_FULL,
  identical(
    as.character(chk_rds$annotation$global_cell),
    as.character(combined$global_cell)
  )
)

# ============================================================
# SHA manifest
# ============================================================

sha_files <- c(
  annotation_file,
  raw_file,
  rds_file,
  file.path(out_dir, "erythroid_full_transfer_policy_v1.tsv"),
  file.path(out_dir, "erythroid_r0p4_annotation_decision_v1.tsv"),
  file.path(out_dir, "erythroid_full_core_counts_v1.tsv"),
  file.path(out_dir, "erythroid_full_maturation_counts_v1.tsv"),
  file.path(out_dir, "erythroid_full_state_counts_v1.tsv"),
  file.path(out_dir, "erythroid_full_transfer_confidence_counts_v1.tsv"),
  file.path(out_dir, "erythroid_full_query_raw_class_counts_v1.tsv"),
  file.path(out_dir, "erythroid_full_core_maturation_state_counts_v1.tsv")
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
  stop("sha256sum creation failed")
}

status <- system2(
  "sha256sum",
  c("-c", sha_out),
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
    "Erythroid full annotation freeze v1",
    "annotation_version=Erythroid_annotation_v1",
    paste0("n_cells=", N_FULL),
    paste0("n_discovery_reference_cells=", N_DISCOVERY),
    paste0("n_full_query_cells=", N_QUERY),
    "",
    "discovery_backbone=balanced RPCA r0.4",
    "biological evidence=native project-balanced markers + absolute audit + cell-level co-detection audit",
    "full query mapping=Seurat pcaproject + TransferData",
    "",
    "core Erythroid_like threshold=0.70",
    "core mixed-lineage threshold=0.70",
    "query myeloid prediction NOT promoted",
    "regulatory-rich maturation threshold=0.70",
    "late-globin maturation threshold=0.90",
    "terminal-hemoglobinization threshold=0.90",
    "state none threshold=0.80",
    "",
    "discovery annotation retained unchanged",
    "raw query predictions retained",
    "core / maturation / state / evidence flag remain separate",
    "",
    "gzip_integrity=PASS",
    "full_annotation_re_read_count_order=PASS",
    "raw_prediction_re_read=PASS",
    "freeze_rds_re_read=PASS",
    "sha256_manifest=PASS",
    "cellranger_count not accessed",
    paste0("freeze_rds=", rds_file),
    paste0("freeze_tsv=", annotation_file)
  ),
  done_file
)

cat("\n===== CORE =====\n")
print(core_counts, row.names=FALSE)

cat("\n===== MATURATION =====\n")
print(maturation_counts, row.names=FALSE)

cat("\n===== STATE =====\n")
print(state_counts, row.names=FALSE)

cat("\n===== TRANSFER CONFIDENCE =====\n")
print(transfer_confidence_counts, row.names=FALSE)

cat("\n===== RAW QUERY CLASSES =====\n")
print(raw_query_class_counts, row.names=FALSE)

cat("\n===== JOINT =====\n")
print(joint, row.names=FALSE)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Erythroid full annotation freeze completed\n"
)
