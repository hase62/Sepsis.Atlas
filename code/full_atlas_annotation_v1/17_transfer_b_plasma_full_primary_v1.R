#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

set.seed(20260814)

source(
  file.path(
    root,
    "full_atlas_primary_integration_v1",
    "00_common.R"
  )
)

# ============================================================
# Frozen thresholds from pilot holdout calibration
#
# score >= 0.80       high
# 0.70 <= score <0.80 medium
# score < 0.70        unresolved
#
# Calibration:
#   overall exact taxonomy = 0.975709
#   overall core identity  = 0.981275
#   overall state          = 0.984312
#
# score 0.70-0.80:
#   exact taxonomy = 0.965517
#   core identity  = 0.988506
#
# score 0.60-0.70:
#   exact taxonomy = 0.828125
#
# Therefore <0.70 is NOT accepted as a frozen annotation.
# ============================================================

HIGH_THRESHOLD <- 0.80
ACCEPT_THRESHOLD <- 0.70

# ============================================================
# Paths
# ============================================================

pilot_rds <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_pilot_annotation_freeze_v1",
  "b_plasma_pilot_native_subclustering__annotation_freeze_v1.rds"
)

taxonomy_file <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_pilot_annotation_freeze_v1",
  "b_plasma_pilot_cluster_taxonomy_freeze_v1.tsv"
)

cal_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_full_transfer_calibration_v1"
)

feature_file <- file.path(
  cal_dir,
  "b_plasma_transfer_features_v1.tsv"
)

cal_done <- file.path(
  cal_dir,
  "B_PLASMA_FULL_TRANSFER_CALIBRATION_COMPLETE.ok"
)

lib_file <- file.path(
  root,
  "pre_integration",
  "pilot_unintegrated_large_v1",
  "resolved_library_table.tsv"
)

transfer_dir <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "annotation_transfer_by_library"
)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_full_primary_transfer_v1"
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

done_file <- file.path(
  out_dir,
  "B_PLASMA_FULL_PRIMARY_TRANSFER_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop(
    "Output already completed; refusing to overwrite: ",
    done_file
  )
}

for(f in c(
  pilot_rds,
  taxonomy_file,
  feature_file,
  cal_done,
  lib_file
)){
  if(!file.exists(f)){
    stop("Missing required input: ", f)
  }
}

# ============================================================
# Load frozen pilot
# ============================================================

pilot <- readRDS(pilot_rds)

stopifnot(
  nrow(pilot) == 38606L,
  ncol(pilot) == 9997L
)

pilot_md <- pilot[[]]
pilot_cells <- colnames(pilot)

required_pilot_md <- c(
  "project_id",
  "cluster_res_0p4",
  "b_plasma_core_identity_v1",
  "b_plasma_state_v1",
  "b_plasma_annotation_confidence_v1",
  "b_plasma_annotation_scope_v1"
)

stopifnot(
  all(required_pilot_md %in% names(pilot_md))
)

pilot_counts <- tryCatch(
  LayerData(
    pilot[["RNA"]],
    layer="counts"
  ),
  error=function(e) NULL
)

if(is.null(pilot_counts)){
  pilot_counts <- GetAssayData(
    pilot,
    assay="RNA",
    layer="counts"
  )
}

stopifnot(
  nrow(pilot_counts) == 38606L,
  ncol(pilot_counts) == 9997L
)

# ============================================================
# Frozen taxonomy
# ============================================================

tax <- read.delim(
  taxonomy_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

stopifnot(
  nrow(tax) == 15L
)

tax$cluster <- as.character(
  tax$cluster
)

tax$taxonomy_key <- paste(
  tax$core_identity,
  tax$state,
  sep="|||"
)

# Collapse identical cluster-level taxonomies.
combo <- unique(
  tax[,c(
    "taxonomy_key",
    "core_identity",
    "state",
    "annotation_scope"
  )]
)

# Keep deterministic order based on first cluster occurrence.
combo <- combo[
  match(
    unique(tax$taxonomy_key),
    combo$taxonomy_key
  ),
  ,
  drop=FALSE
]

combo$taxonomy_id <- sprintf(
  "T%02d",
  seq_len(nrow(combo))
)

combo <- combo[,c(
  "taxonomy_id",
  "taxonomy_key",
  "core_identity",
  "state",
  "annotation_scope"
)]

stopifnot(
  nrow(combo) == 12L
)

tax$taxonomy_id <- combo$taxonomy_id[
  match(
    tax$taxonomy_key,
    combo$taxonomy_key
  )
]

# Each transferable taxonomy should have one biological
# confidence level in the frozen taxonomy.
confidence_by_taxonomy <- lapply(
  combo$taxonomy_id,
  function(k){

    clusters <- tax$cluster[
      tax$taxonomy_id == k
    ]

    z <- unique(
      as.character(
        tax$confidence[
          tax$cluster %in% clusters
        ]
      )
    )

    if(length(z) != 1L){
      stop(
        "Non-unique frozen biological confidence for taxonomy ",
        k,
        ": ",
        paste(z, collapse=",")
      )
    }

    z
  }
)

combo$frozen_taxonomy_confidence <-
  unlist(
    confidence_by_taxonomy,
    use.names=FALSE
  )

# ============================================================
# Verify pilot against taxonomy
# ============================================================

pilot_cluster <- as.character(
  pilot_md$cluster_res_0p4
)

ii <- match(
  pilot_cluster,
  tax$cluster
)

stopifnot(
  !anyNA(ii)
)

pilot_taxonomy_id <- tax$taxonomy_id[ii]

stopifnot(
  identical(
    as.character(
      pilot_md$b_plasma_core_identity_v1
    ),
    as.character(
      tax$core_identity[ii]
    )
  ),
  identical(
    as.character(
      pilot_md$b_plasma_state_v1
    ),
    as.character(
      tax$state[ii]
    )
  )
)

names(pilot_taxonomy_id) <-
  rownames(pilot_md)

# ============================================================
# Frozen transfer features from calibration
# ============================================================

feature_table <- read.delim(
  feature_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

stopifnot(
  "gene" %in% names(feature_table)
)

transfer_features <- unique(
  as.character(
    feature_table$gene
  )
)

stopifnot(
  length(transfer_features) == 2849L,
  all(
    transfer_features %in%
      rownames(pilot_counts)
  )
)

# ============================================================
# Reconstruct FULL B/plasma metadata and QUERY native counts
#
# Pilot cells are excluded from query counts.
# We do not duplicate the 9,997-cell reference in query.
# ============================================================

lib <- read_tsv(lib_file)

stopifnot(
  all(c(
    "project_id",
    "library_key",
    "final_rds_resolved"
  ) %in% names(lib))
)

query_count_list <- list()
full_md_list <- list()

n_full <- 0L
n_query <- 0L
n_bp_libraries <- 0L

for(i in seq_len(nrow(lib))){

  project_id <- as.character(
    lib$project_id[[i]]
  )

  library_key <- as.character(
    lib$library_key[[i]]
  )

  tag <- safe_name(
    paste(
      project_id,
      library_key,
      sep="__"
    )
  )

  tr_file <- file.path(
    transfer_dir,
    paste0(tag, ".tsv.gz")
  )

  if(!file.exists(tr_file)){
    stop("Missing transfer file: ", tr_file)
  }

  tr <- read_tsv(tr_file)

  tr <- tr[
    as.character(
      tr$integration_compartment_primary_v1
    ) == "B_plasma",
    ,
    drop=FALSE
  ]

  if(!nrow(tr)){
    next
  }

  n_bp_libraries <- n_bp_libraries + 1L
  n_full <- n_full + nrow(tr)

  global_cells <- as.character(
    tr$global_cell
  )

  if(anyDuplicated(global_cells)){
    stop(
      "Duplicated global_cell within library: ",
      tag
    )
  }

  # Keep full metadata, including pilot/reference cells.
  md_block <- tr
  rownames(md_block) <- global_cells

  full_md_list[[
    length(full_md_list)+1L
  ]] <- md_block

  # Only non-pilot cells become transfer query.
  keep_query <- !global_cells %in%
    pilot_cells

  if(!any(keep_query)){
    next
  }

  obj <- readRDS(
    as.character(
      lib$final_rds_resolved[[i]]
    )
  )

  counts <- get_rna_counts(obj)

  stopifnot(
    identical(
      rownames(counts),
      rownames(pilot_counts)
    )
  )

  original_query <- as.character(
    tr$original_cell[keep_query]
  )

  global_query <- global_cells[
    keep_query
  ]

  missing_cells <- setdiff(
    original_query,
    colnames(counts)
  )

  if(length(missing_cells)){
    stop(
      "Missing query cells in source object ",
      tag,
      ": ",
      length(missing_cells)
    )
  }

  # Only transfer features are needed for reference mapping.
  # Avoid materializing the full 38,606-gene x 35,958-cell query matrix.
  q <- counts[
    transfer_features,
    original_query,
    drop=FALSE
  ]

  colnames(q) <- global_query

  query_count_list[[
    length(query_count_list)+1L
  ]] <- q

  n_query <- n_query + ncol(q)

  cat(
    sprintf(
      "LOAD %3d/%3d  %-12s %-40s full=%4d query=%4d\n",
      i,
      nrow(lib),
      project_id,
      library_key,
      nrow(tr),
      ncol(q)
    )
  )

  rm(
    obj,
    counts,
    q,
    tr
  )

  gc(verbose=FALSE)
}

full_md <- do.call(
  rbind,
  full_md_list
)

query_counts <- do.call(
  cbind,
  query_count_list
)

rm(
  query_count_list,
  full_md_list
)

gc(verbose=FALSE)

stopifnot(
  n_full == 45955L,
  n_bp_libraries == 134L,
  nrow(full_md) == 45955L,
  !anyDuplicated(
    rownames(full_md)
  ),
  n_query == 35958L,
  ncol(query_counts) == 35958L,
  nrow(query_counts) == length(transfer_features),
  !anyDuplicated(
    colnames(query_counts)
  )
)

stopifnot(
  setequal(
    pilot_cells,
    intersect(
      pilot_cells,
      rownames(full_md)
    )
  )
)

query_cells <- rownames(full_md)[
  !rownames(full_md) %in%
    pilot_cells
]

stopifnot(
  length(query_cells) == 35958L,
  setequal(
    query_cells,
    colnames(query_counts)
  )
)

# Reorder query matrix to full metadata order.
query_counts <- query_counts[
  ,
  query_cells,
  drop=FALSE
]

# ============================================================
# Build reference
# ============================================================

ref_md <- pilot_md

ref_md$taxonomy_id <-
  pilot_taxonomy_id[
    rownames(ref_md)
  ]

ref <- CreateSeuratObject(
  counts=pilot_counts[
    transfer_features,
    ,
    drop=FALSE
  ],
  meta.data=ref_md,
  assay="RNA",
  project="B_plasma_frozen_reference"
)

ref <- NormalizeData(
  ref,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

ref <- ScaleData(
  ref,
  features=transfer_features,
  verbose=FALSE
)

ref <- RunPCA(
  ref,
  features=transfer_features,
  npcs=30,
  verbose=FALSE
)

# ============================================================
# Build query
# ============================================================

query_md <- full_md[
  query_cells,
  ,
  drop=FALSE
]

qry <- CreateSeuratObject(
  counts=query_counts,
  meta.data=query_md,
  assay="RNA",
  project="B_plasma_full_primary_query"
)

rm(
  query_counts,
  pilot_counts
)

gc(verbose=FALSE)

qry <- NormalizeData(
  qry,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

qry <- ScaleData(
  qry,
  features=transfer_features,
  verbose=FALSE
)

# ============================================================
# Find anchors and transfer frozen taxonomy
# ============================================================

anchors <- FindTransferAnchors(
  reference=ref,
  query=qry,
  normalization.method="LogNormalize",
  reference.assay="RNA",
  query.assay="RNA",
  reduction="pcaproject",
  reference.reduction="pca",
  features=transfer_features,
  npcs=30,
  dims=1:30,
  k.anchor=5,
  k.score=30,
  nn.method="annoy",
  verbose=TRUE
)

ref_labels <- setNames(
  as.character(
    ref$taxonomy_id
  ),
  colnames(ref)
)

pred <- TransferData(
  anchorset=anchors,
  refdata=ref_labels,
  dims=1:30,
  k.weight=50,
  verbose=TRUE
)

stopifnot(
  nrow(pred) == 35958L,
  identical(
    rownames(pred),
    colnames(qry)
  )
)

# Checkpoint raw Seurat prediction output immediately.
saveRDS(
  pred,
  file=file.path(
    out_dir,
    "b_plasma_query_transferdata_raw_v1.rds"
  ),
  compress=FALSE
)

# ============================================================
# Score matrix and margin
# ============================================================

expected_score_cols <- paste0(
  "prediction.score.",
  combo$taxonomy_id
)

missing_score_cols <- setdiff(
  expected_score_cols,
  names(pred)
)

if(length(missing_score_cols)){
  stop(
    "Missing prediction score columns: ",
    paste(
      missing_score_cols,
      collapse=", "
    )
  )
}

score_mat <- as.matrix(
  pred[
    ,
    expected_score_cols,
    drop=FALSE
  ]
)

colnames(score_mat) <- sub(
  "^prediction\\.score\\.",
  "",
  colnames(score_mat)
)

top_two <- t(
  apply(
    score_mat,
    1,
    function(z){

      s <- sort(
        as.numeric(z),
        decreasing=TRUE
      )

      if(length(s) < 2L){
        stop(
          "Fewer than two taxonomy scores"
        )
      }

      c(
        top1=s[[1]],
        top2=s[[2]]
      )
    }
  )
)

colnames(top_two) <- c(
  "top1",
  "top2"
)

prediction_margin <-
  top_two[,"top1"] -
  top_two[,"top2"]

# ============================================================
# Decode raw query predictions
# ============================================================

pred_id <- as.character(
  pred$predicted.id
)

pred_idx <- match(
  pred_id,
  combo$taxonomy_id
)

stopifnot(
  !anyNA(pred_idx)
)

prediction_score <- as.numeric(
  pred$prediction.score.max
)

transfer_confidence <- ifelse(
  prediction_score >= HIGH_THRESHOLD,
  "high",
  ifelse(
    prediction_score >= ACCEPT_THRESHOLD,
    "medium",
    "unresolved"
  )
)

accepted <- prediction_score >=
  ACCEPT_THRESHOLD

query_result <- data.frame(
  cell_id=colnames(qry),

  project_id=as.character(
    qry$project_id
  ),

  library_key=as.character(
    qry$library_key
  ),

  raw_predicted_taxonomy_id=
    pred_id,

  raw_predicted_core_identity=
    combo$core_identity[pred_idx],

  raw_predicted_state=
    combo$state[pred_idx],

  raw_predicted_annotation_scope=
    combo$annotation_scope[pred_idx],

  raw_predicted_frozen_taxonomy_confidence=
    combo$frozen_taxonomy_confidence[pred_idx],

  prediction_score=
    prediction_score,

  prediction_margin=
    as.numeric(
      prediction_margin
    ),

  transfer_confidence=
    transfer_confidence,

  transfer_accepted=
    accepted,

  stringsAsFactors=FALSE
)

# ============================================================
# Final query calls
#
# Low-score raw prediction is retained but NOT promoted
# into the final frozen taxonomy.
# ============================================================

query_result$final_taxonomy_id <-
  ifelse(
    accepted,
    query_result$raw_predicted_taxonomy_id,
    "TRANSFER_UNRESOLVED"
  )

query_result$final_core_identity <-
  ifelse(
    accepted,
    query_result$raw_predicted_core_identity,
    "B_plasma_transfer_unresolved"
  )

query_result$final_state <-
  ifelse(
    accepted,
    query_result$raw_predicted_state,
    "unresolved"
  )

query_result$final_annotation_scope <-
  ifelse(
    accepted,
    query_result$raw_predicted_annotation_scope,
    "unresolved"
  )

query_result$final_frozen_taxonomy_confidence <-
  ifelse(
    accepted,
    query_result$raw_predicted_frozen_taxonomy_confidence,
    "unresolved"
  )

query_result$annotation_source <-
  "full_transfer_v1"

rownames(query_result) <-
  query_result$cell_id

# ============================================================
# Reference/pilot final rows
# ============================================================

pilot_idx <- match(
  pilot_taxonomy_id,
  combo$taxonomy_id
)

stopifnot(
  !anyNA(pilot_idx)
)

pilot_result <- data.frame(
  cell_id=pilot_cells,

  project_id=as.character(
    pilot_md$project_id
  ),

  library_key=if(
    "library_key" %in%
      names(pilot_md)
  ){
    as.character(
      pilot_md$library_key
    )
  } else {
    NA_character_
  },

  raw_predicted_taxonomy_id=
    NA_character_,

  raw_predicted_core_identity=
    NA_character_,

  raw_predicted_state=
    NA_character_,

  raw_predicted_annotation_scope=
    NA_character_,

  raw_predicted_frozen_taxonomy_confidence=
    NA_character_,

  prediction_score=
    NA_real_,

  prediction_margin=
    NA_real_,

  transfer_confidence=
    "reference_frozen",

  transfer_accepted=
    TRUE,

  final_taxonomy_id=
    pilot_taxonomy_id,

  final_core_identity=
    as.character(
      pilot_md$b_plasma_core_identity_v1
    ),

  final_state=
    as.character(
      pilot_md$b_plasma_state_v1
    ),

  final_annotation_scope=
    as.character(
      pilot_md$b_plasma_annotation_scope_v1
    ),

  final_frozen_taxonomy_confidence=
    as.character(
      pilot_md$b_plasma_annotation_confidence_v1
    ),

  annotation_source=
    "pilot_freeze_v1",

  stringsAsFactors=FALSE
)

rownames(pilot_result) <-
  pilot_result$cell_id

# ============================================================
# Combine in canonical full-cell order
# ============================================================

result_all <- rbind(
  pilot_result,
  query_result
)

stopifnot(
  nrow(result_all) == 45955L,
  !anyDuplicated(
    result_all$cell_id
  ),
  setequal(
    result_all$cell_id,
    rownames(full_md)
  )
)

result_all <- result_all[
  rownames(full_md),
  ,
  drop=FALSE
]

stopifnot(
  identical(
    result_all$cell_id,
    rownames(full_md)
  )
)

# ============================================================
# Threshold provenance
# ============================================================

thresholds <- data.frame(
  parameter=c(
    "high_threshold",
    "accept_threshold"
  ),
  value=c(
    HIGH_THRESHOLD,
    ACCEPT_THRESHOLD
  ),
  interpretation=c(
    "score >= 0.80 -> high confidence accepted transfer",
    "score >= 0.70 -> accepted; below 0.70 -> unresolved"
  ),
  calibration_basis=c(
    "0.80-0.90 exact taxonomy accuracy 0.988095; >=0.90 accuracy 1.0",
    "0.70-0.80 exact taxonomy accuracy 0.965517; 0.60-0.70 falls to 0.828125"
  ),
  stringsAsFactors=FALSE
)

write.table(
  thresholds,
  file=file.path(
    out_dir,
    "b_plasma_transfer_thresholds_freeze_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

write.table(
  combo,
  file=file.path(
    out_dir,
    "b_plasma_full_transfer_taxonomy_map_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# Write query predictions
# ============================================================

write.table(
  query_result,
  file=gzfile(
    file.path(
      out_dir,
      "b_plasma_full_query_predictions_v1.tsv.gz"
    ),
    open="wt"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# Full 45,955-cell annotation table
# ============================================================

full_annotation <- cbind(
  full_md,
  result_all[
    ,
    setdiff(
      names(result_all),
      c(
        "cell_id",
        "project_id",
        "library_key"
      )
    ),
    drop=FALSE
  ]
)

full_annotation$cell_id <-
  rownames(full_annotation)

# Put cell_id first.
full_annotation <- full_annotation[
  ,
  c(
    "cell_id",
    setdiff(
      names(full_annotation),
      "cell_id"
    )
  ),
  drop=FALSE
]

write.table(
  full_annotation,
  file=gzfile(
    file.path(
      out_dir,
      "b_plasma_full_primary_annotation_v1.tsv.gz"
    ),
    open="wt"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# Summary
# ============================================================

query_summary <- as.data.frame(
  table(
    transfer_confidence=
      query_result$transfer_confidence
  ),
  stringsAsFactors=FALSE
)

taxonomy_summary <- as.data.frame(
  table(
    annotation_source=
      result_all$annotation_source,
    final_taxonomy_id=
      result_all$final_taxonomy_id
  ),
  stringsAsFactors=FALSE
)

taxonomy_summary <- taxonomy_summary[
  taxonomy_summary$Freq > 0,
  ,
  drop=FALSE
]

project_summary <- do.call(
  rbind,
  lapply(
    sort(
      unique(
        query_result$project_id
      )
    ),
    function(p){

      d <- query_result[
        query_result$project_id == p,
        ,
        drop=FALSE
      ]

      data.frame(
        project_id=p,
        n_query=nrow(d),
        n_high=sum(
          d$transfer_confidence ==
            "high"
        ),
        n_medium=sum(
          d$transfer_confidence ==
            "medium"
        ),
        n_unresolved=sum(
          d$transfer_confidence ==
            "unresolved"
        ),
        accepted_fraction=
          mean(
            d$transfer_accepted
          ),
        median_prediction_score=
          median(
            d$prediction_score
          ),
        median_prediction_margin=
          median(
            d$prediction_margin
          ),
        stringsAsFactors=FALSE
      )
    }
  )
)

write.table(
  query_summary,
  file=file.path(
    out_dir,
    "b_plasma_full_query_confidence_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

write.table(
  taxonomy_summary,
  file=file.path(
    out_dir,
    "b_plasma_full_taxonomy_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

write.table(
  project_summary,
  file=file.path(
    out_dir,
    "b_plasma_full_transfer_by_project_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# Final integrity
# ============================================================

n_high <- sum(
  query_result$transfer_confidence ==
    "high"
)

n_medium <- sum(
  query_result$transfer_confidence ==
    "medium"
)

n_unresolved <- sum(
  query_result$transfer_confidence ==
    "unresolved"
)

stopifnot(
  n_high +
    n_medium +
    n_unresolved ==
    35958L
)

# ============================================================
# Completion marker
# ============================================================

writeLines(
  c(
    "PASS",
    "B/plasma full primary transfer v1",
    "n_full_cells=45955",
    "n_reference_pilot_cells=9997",
    "n_query_cells=35958",
    "n_features=38606",
    "n_transfer_features=2849",
    "n_projects=9",
    "n_B_plasma_libraries=134",
    "",
    "mapping=Seurat pcaproject + TransferData",
    "reference=all 9997 frozen pilot cells",
    "query=remaining 35958 full-primary B/plasma cells",
    "",
    "frozen transfer thresholds:",
    "score>=0.80 high accepted",
    "0.70<=score<0.80 medium accepted",
    "score<0.70 unresolved; raw prediction retained but not promoted",
    "",
    paste0(
      "query_high=",
      n_high
    ),
    paste0(
      "query_medium=",
      n_medium
    ),
    paste0(
      "query_unresolved=",
      n_unresolved
    ),
    paste0(
      "query_accepted_fraction=",
      sprintf(
        "%.6f",
        mean(
          query_result$transfer_accepted
        )
      )
    ),
    "",
    "pilot labels are never overwritten by prediction",
    "core identity and orthogonal state remain separate",
    "deferred non-B classes remain explicit",
    "native RNA expression is not corrected or overwritten"
  ),
  done_file
)

cat(
  "\n===== QUERY TRANSFER CONFIDENCE =====\n"
)

print(
  query_summary,
  row.names=FALSE
)

cat(
  "\n===== FULL TAXONOMY COUNTS =====\n"
)

print(
  taxonomy_summary,
  row.names=FALSE
)

cat(
  "\n===== TRANSFER BY PROJECT =====\n"
)

print(
  project_summary,
  row.names=FALSE
)

cat(
  "\nPASS: B/plasma full primary transfer v1 completed\n"
)

