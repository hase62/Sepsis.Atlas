#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

set.seed(20260814)

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

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_full_transfer_calibration_v1"
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

done_file <- file.path(
  out_dir,
  "B_PLASMA_FULL_TRANSFER_CALIBRATION_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop(
    "Calibration already completed; refusing to overwrite: ",
    done_file
  )
}

for(f in c(
  pilot_rds,
  taxonomy_file
)){
  if(!file.exists(f)){
    stop("Missing input: ", f)
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

md <- pilot[[]]

required_md <- c(
  "project_id",
  "cluster_res_0p4",
  "b_plasma_core_identity_v1",
  "b_plasma_state_v1",
  "b_plasma_annotation_confidence_v1",
  "b_plasma_annotation_scope_v1"
)

stopifnot(
  all(required_md %in% names(md))
)

counts <- tryCatch(
  LayerData(
    pilot[["RNA"]],
    layer="counts"
  ),
  error=function(e) NULL
)

if(is.null(counts)){
  counts <- GetAssayData(
    pilot,
    assay="RNA",
    layer="counts"
  )
}

stopifnot(
  nrow(counts) == 38606L,
  ncol(counts) == 9997L,
  identical(
    colnames(counts),
    rownames(md)
  )
)

# ============================================================
# Frozen taxonomy -> compact transferable taxonomy IDs
#
# Multiple pilot clusters with the same final
# core_identity + state are intentionally collapsed.
# ============================================================

tax <- read.delim(
  taxonomy_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

required_tax <- c(
  "cluster",
  "core_identity",
  "state",
  "confidence",
  "annotation_scope"
)

stopifnot(
  all(required_tax %in% names(tax)),
  nrow(tax) == 15L,
  setequal(
    as.character(tax$cluster),
    as.character(1:15)
  )
)

tax$cluster <- as.character(
  tax$cluster
)

tax <- tax[
  order(as.integer(tax$cluster)),
  ,
  drop=FALSE
]

tax$taxonomy_key <- paste(
  tax$core_identity,
  tax$state,
  sep="|||"
)

combo <- tax[
  !duplicated(tax$taxonomy_key),
  c(
    "taxonomy_key",
    "core_identity",
    "state",
    "annotation_scope"
  ),
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

tax$taxonomy_id <- combo$taxonomy_id[
  match(
    tax$taxonomy_key,
    combo$taxonomy_key
  )
]

stopifnot(
  !anyNA(tax$taxonomy_id)
)

# ============================================================
# Verify the frozen object still agrees exactly with taxonomy
# ============================================================

cluster <- as.character(
  md$cluster_res_0p4
)

ii <- match(
  cluster,
  tax$cluster
)

stopifnot(
  !anyNA(ii),
  identical(
    as.character(md$b_plasma_core_identity_v1),
    as.character(tax$core_identity[ii])
  ),
  identical(
    as.character(md$b_plasma_state_v1),
    as.character(tax$state[ii])
  )
)

md$taxonomy_id <- tax$taxonomy_id[ii]
md$taxonomy_key <- tax$taxonomy_key[ii]

stopifnot(
  !anyNA(md$taxonomy_id)
)

cat(
  "\n===== FROZEN TAXONOMY CLASSES =====\n"
)

print(
  combo,
  row.names=FALSE
)

cat(
  "\n===== PILOT TAXONOMY COUNTS =====\n"
)

print(
  table(md$taxonomy_id)
)

# ============================================================
# Deterministic stratified holdout
#
# Stratify by project + frozen taxonomy.
#
# n >= 10 : ~20% holdout
# n 4..9  : 1 holdout
# n < 4   : reference only
#
# This preserves tiny project/taxonomy strata in reference.
# ============================================================

stratum <- paste(
  as.character(md$project_id),
  as.character(md$taxonomy_id),
  sep="|||"
)

validation_cells <- character()

for(s in sort(unique(stratum))){

  cells <- rownames(md)[
    stratum == s
  ]

  n <- length(cells)

  n_val <- if(n >= 10L){
    max(
      1L,
      as.integer(
        floor(0.20 * n)
      )
    )
  } else if(n >= 4L){
    1L
  } else {
    0L
  }

  n_val <- min(
    n_val,
    n - 1L
  )

  if(n_val > 0L){

    validation_cells <- c(
      validation_cells,
      sample(
        cells,
        size=n_val,
        replace=FALSE
      )
    )
  }
}

validation_cells <- unique(
  validation_cells
)

reference_cells <- setdiff(
  colnames(pilot),
  validation_cells
)

stopifnot(
  length(reference_cells) +
    length(validation_cells) ==
    9997L,
  !length(
    intersect(
      reference_cells,
      validation_cells
    )
  )
)

# Every taxonomy class must remain in reference.
stopifnot(
  setequal(
    unique(md[reference_cells,"taxonomy_id"]),
    unique(md$taxonomy_id)
  )
)

# Ensure every sufficiently represented taxonomy has validation.
tax_counts <- table(
  md$taxonomy_id
)

for(k in names(tax_counts)[
  tax_counts >= 20
]){

  if(!any(
    md[validation_cells,"taxonomy_id"] == k
  )){
    stop(
      "No validation cells for sufficiently represented taxonomy: ",
      k
    )
  }
}

cat(
  "\nreference cells  =",
  length(reference_cells),
  "\n"
)

cat(
  "validation cells =",
  length(validation_cells),
  "\n"
)

# ============================================================
# Create clean reference/query objects from native counts
#
# Do NOT use integrated/adjusted expression.
# ============================================================

ref_md <- md[
  reference_cells,
  ,
  drop=FALSE
]

val_md <- md[
  validation_cells,
  ,
  drop=FALSE
]

ref <- CreateSeuratObject(
  counts=counts[
    ,
    reference_cells,
    drop=FALSE
  ],
  meta.data=ref_md,
  assay="RNA",
  project="B_plasma_pilot_reference"
)

qry <- CreateSeuratObject(
  counts=counts[
    ,
    validation_cells,
    drop=FALSE
  ],
  meta.data=val_md,
  assay="RNA",
  project="B_plasma_pilot_validation"
)

rm(pilot)
gc(verbose=FALSE)

# ============================================================
# Normalize
# ============================================================

ref <- NormalizeData(
  ref,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

qry <- NormalizeData(
  qry,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

ref <- FindVariableFeatures(
  ref,
  selection.method="vst",
  nfeatures=3000,
  verbose=FALSE
)

# ============================================================
# Frozen transfer feature set
#
# Start with reference-variable genes.
# Add canonical B identity/state/contamination axes.
#
# Constant immunoglobulin genes are retained.
# Rearranged V/J genes are excluded.
# ============================================================

canonical_transfer <- c(

  # naive B
  "TCL1A",
  "IGHD",
  "IGHM",
  "FCER2",
  "IL4R",
  "BACH2",
  "FCRL1",
  "VPREB3",

  # memory B
  "CD27",
  "TNFRSF13B",
  "AIM2",
  "GPR183",
  "CD82",
  "ITGB1",
  "POU2AF1",

  # class-switched memory
  "IGHA1",
  "IGHA2",
  "IGHG1",
  "IGHG2",
  "IGHG3",
  "IGHG4",

  # interferon
  "IFI6",
  "IFI27",
  "IFI35",
  "IFI44",
  "IFI44L",
  "IFIT1",
  "IFIT2",
  "IFIT3",
  "IFITM1",
  "IFITM2",
  "IFITM3",
  "ISG15",
  "MX1",
  "MX2",
  "OAS1",
  "OAS2",
  "OAS3",
  "OASL",
  "XAF1",
  "BST2",
  "LY6E",
  "EPSTI1",

  # early activation / stress
  "FOS",
  "FOSB",
  "JUN",
  "JUNB",
  "EGR1",
  "EGR2",
  "ATF3",
  "NFKBIA",
  "DUSP1",
  "DUSP2",
  "CD69",
  "NR4A1",
  "NR4A2",
  "PPP1R15A",

  # activated / atypical candidate
  "ITGAX",
  "CD86",
  "TBX21",
  "FCRL5",

  # secretory differentiation
  "XBP1",
  "PRDM1",
  "JCHAIN",
  "MZB1",
  "FKBP11",
  "DERL3",
  "SEC11C",
  "ELL2",
  "IRF4",
  "TNFRSF17",
  "SDC1",

  # platelet-like deferred non-B
  "PPBP",
  "PF4",
  "PF4V1",
  "GP1BB",
  "GP9",
  "ITGA2B",
  "TUBB1",
  "TREML1",
  "MPIG6B",

  # T-cell-like deferred non-B
  "CD3D",
  "CD3E",
  "CD3G",
  "TRAC",
  "TRBC1",
  "TRBC2",
  "CD247",
  "CD2",
  "LCK",
  "NKG7",
  "GNLY",
  "CCL5"
)

candidate_features <- unique(c(
  VariableFeatures(ref),
  canonical_transfer
))

candidate_features <- intersect(
  candidate_features,
  rownames(ref)
)

technical <- grepl(
  paste0(
    "^MT-|",
    "^RPL[0-9]|",
    "^RPS[0-9]|",
    "^HBA[12]$|",
    "^HBB$|",
    "^HBD$|",
    "^HBG[12]$|",
    "^MALAT1$|",
    "^NEAT1$"
  ),
  candidate_features
)

rearranged_ig <- grepl(
  paste0(
    "^IGHV|",
    "^IGKV|",
    "^IGLV|",
    "^IGHJ[0-9]|",
    "^IGKJ[0-9]|",
    "^IGLJ[0-9]"
  ),
  candidate_features
)

transfer_features <- candidate_features[
  !technical &
  !rearranged_ig
]

transfer_features <- unique(
  transfer_features
)

stopifnot(
  length(transfer_features) >= 1000L
)

cat(
  "transfer features =",
  length(transfer_features),
  "\n"
)

write.table(
  data.frame(
    gene=transfer_features,
    variable_feature=
      transfer_features %in%
        VariableFeatures(ref),
    canonical_transfer=
      transfer_features %in%
        canonical_transfer,
    stringsAsFactors=FALSE
  ),
  file=file.path(
    out_dir,
    "b_plasma_transfer_features_v1.tsv"
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
    "b_plasma_transfer_taxonomy_map_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# PCA projection reference mapping
# ============================================================

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

qry <- ScaleData(
  qry,
  features=transfer_features,
  verbose=FALSE
)

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
  as.character(ref$taxonomy_id),
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
  nrow(pred) == ncol(qry),
  identical(
    rownames(pred),
    colnames(qry)
  )
)

# ============================================================
# Prediction score margin
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
  pred[,expected_score_cols,drop=FALSE]
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

      # Drop class names before taking the two largest scores.
      # Otherwise c(top1=s[1]) produces names such as
      # "top1.T01", which breaks downstream column lookup.
      s <- sort(
        as.numeric(z),
        decreasing=TRUE
      )

      if(length(s) < 2L){
        stop("Fewer than two taxonomy prediction scores")
      }

      c(
        top1=s[[1]],
        top2=s[[2]]
      )
    }
  )
)

stopifnot(
  is.matrix(top_two),
  nrow(top_two) == nrow(score_mat),
  ncol(top_two) == 2L
)

colnames(top_two) <- c(
  "top1",
  "top2"
)

prediction_margin <-
  top_two[,"top1"] -
  top_two[,"top2"]

# ============================================================
# Decode truth / prediction
# ============================================================

true_id <- as.character(
  qry$taxonomy_id
)

pred_id <- as.character(
  pred$predicted.id
)

true_idx <- match(
  true_id,
  combo$taxonomy_id
)

pred_idx <- match(
  pred_id,
  combo$taxonomy_id
)

stopifnot(
  !anyNA(true_idx),
  !anyNA(pred_idx)
)

result <- data.frame(
  cell_id=colnames(qry),
  project_id=as.character(qry$project_id),

  true_taxonomy_id=true_id,
  predicted_taxonomy_id=pred_id,

  true_core_identity=
    combo$core_identity[true_idx],

  predicted_core_identity=
    combo$core_identity[pred_idx],

  true_state=
    combo$state[true_idx],

  predicted_state=
    combo$state[pred_idx],

  true_scope=
    combo$annotation_scope[true_idx],

  predicted_scope=
    combo$annotation_scope[pred_idx],

  prediction_score=
    as.numeric(
      pred$prediction.score.max
    ),

  prediction_margin=
    as.numeric(
      prediction_margin
    ),

  exact_taxonomy_correct=
    true_id == pred_id,

  core_identity_correct=
    combo$core_identity[true_idx] ==
      combo$core_identity[pred_idx],

  state_correct=
    combo$state[true_idx] ==
      combo$state[pred_idx],

  scope_correct=
    combo$annotation_scope[true_idx] ==
      combo$annotation_scope[pred_idx],

  stringsAsFactors=FALSE
)

# ============================================================
# Main accuracy metrics
# ============================================================

summary_df <- data.frame(
  metric=c(
    "n_pilot_cells",
    "n_reference_cells",
    "n_validation_cells",
    "n_taxonomy_classes",
    "n_transfer_features",
    "exact_taxonomy_accuracy",
    "core_identity_accuracy",
    "state_accuracy",
    "annotation_scope_accuracy",
    "median_prediction_score",
    "median_prediction_margin"
  ),
  value=c(
    nrow(md),
    length(reference_cells),
    length(validation_cells),
    nrow(combo),
    length(transfer_features),
    mean(result$exact_taxonomy_correct),
    mean(result$core_identity_correct),
    mean(result$state_correct),
    mean(result$scope_correct),
    median(result$prediction_score),
    median(result$prediction_margin)
  ),
  stringsAsFactors=FALSE
)

# ============================================================
# Per-taxonomy metrics
# ============================================================

per_taxonomy <- do.call(
  rbind,
  lapply(
    combo$taxonomy_id,
    function(k){

      d <- result[
        result$true_taxonomy_id == k,
        ,
        drop=FALSE
      ]

      map <- combo[
        combo$taxonomy_id == k,
        ,
        drop=FALSE
      ]

      data.frame(
        taxonomy_id=k,
        core_identity=
          map$core_identity,
        state=
          map$state,
        annotation_scope=
          map$annotation_scope,
        n_validation=
          nrow(d),
        exact_accuracy=
          if(nrow(d))
            mean(d$exact_taxonomy_correct)
          else NA_real_,
        core_accuracy=
          if(nrow(d))
            mean(d$core_identity_correct)
          else NA_real_,
        state_accuracy=
          if(nrow(d))
            mean(d$state_correct)
          else NA_real_,
        median_prediction_score=
          if(nrow(d))
            median(d$prediction_score)
          else NA_real_,
        median_prediction_margin=
          if(nrow(d))
            median(d$prediction_margin)
          else NA_real_,
        stringsAsFactors=FALSE
      )
    }
  )
)

# ============================================================
# Per-project metrics
# ============================================================

per_project <- do.call(
  rbind,
  lapply(
    sort(unique(result$project_id)),
    function(p){

      d <- result[
        result$project_id == p,
        ,
        drop=FALSE
      ]

      data.frame(
        project_id=p,
        n_validation=nrow(d),
        exact_taxonomy_accuracy=
          mean(d$exact_taxonomy_correct),
        core_identity_accuracy=
          mean(d$core_identity_correct),
        state_accuracy=
          mean(d$state_correct),
        annotation_scope_accuracy=
          mean(d$scope_correct),
        median_prediction_score=
          median(d$prediction_score),
        median_prediction_margin=
          median(d$prediction_margin),
        stringsAsFactors=FALSE
      )
    }
  )
)

# ============================================================
# Confusion matrix
# ============================================================

conf <- as.data.frame(
  table(
    true_taxonomy_id=
      result$true_taxonomy_id,
    predicted_taxonomy_id=
      result$predicted_taxonomy_id
  ),
  stringsAsFactors=FALSE
)

conf <- conf[
  conf$Freq > 0,
  ,
  drop=FALSE
]

# ============================================================
# Score calibration bins
#
# Thresholds are NOT frozen here.
# These tables are diagnostic evidence for choosing them.
# ============================================================

result$score_bin <- cut(
  result$prediction_score,
  breaks=c(
    -Inf,
    0.50,
    0.60,
    0.70,
    0.80,
    0.90,
    0.95,
    Inf
  ),
  right=FALSE
)

score_bins <- do.call(
  rbind,
  lapply(
    levels(result$score_bin),
    function(b){

      d <- result[
        result$score_bin == b,
        ,
        drop=FALSE
      ]

      data.frame(
        score_bin=b,
        n_cells=nrow(d),
        exact_taxonomy_accuracy=
          if(nrow(d))
            mean(d$exact_taxonomy_correct)
          else NA_real_,
        core_identity_accuracy=
          if(nrow(d))
            mean(d$core_identity_correct)
          else NA_real_,
        state_accuracy=
          if(nrow(d))
            mean(d$state_correct)
          else NA_real_,
        median_margin=
          if(nrow(d))
            median(d$prediction_margin)
          else NA_real_,
        stringsAsFactors=FALSE
      )
    }
  )
)

# ============================================================
# Score/margin quantiles by correctness
# ============================================================

quantile_rows <- list()

for(correct in c(TRUE,FALSE)){

  d <- result[
    result$exact_taxonomy_correct == correct,
    ,
    drop=FALSE
  ]

  if(!nrow(d)){
    next
  }

  probs <- c(
    0,
    0.05,
    0.10,
    0.25,
    0.50,
    0.75,
    0.90,
    0.95,
    1
  )

  qs <- quantile(
    d$prediction_score,
    probs=probs,
    na.rm=TRUE
  )

  qm <- quantile(
    d$prediction_margin,
    probs=probs,
    na.rm=TRUE
  )

  quantile_rows[[
    length(quantile_rows)+1L
  ]] <- data.frame(
    exact_correct=correct,
    quantile=probs,
    prediction_score=
      as.numeric(qs),
    prediction_margin=
      as.numeric(qm),
    stringsAsFactors=FALSE
  )
}

quantiles_df <- do.call(
  rbind,
  quantile_rows
)

# ============================================================
# Write outputs
# ============================================================

write.table(
  result,
  file=file.path(
    out_dir,
    "b_plasma_transfer_calibration_predictions_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

write.table(
  summary_df,
  file=file.path(
    out_dir,
    "b_plasma_transfer_calibration_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

write.table(
  per_taxonomy,
  file=file.path(
    out_dir,
    "b_plasma_transfer_calibration_per_taxonomy_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

write.table(
  per_project,
  file=file.path(
    out_dir,
    "b_plasma_transfer_calibration_per_project_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

write.table(
  conf,
  file=file.path(
    out_dir,
    "b_plasma_transfer_calibration_confusion_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

write.table(
  score_bins,
  file=file.path(
    out_dir,
    "b_plasma_transfer_calibration_score_bins_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

write.table(
  quantiles_df,
  file=file.path(
    out_dir,
    "b_plasma_transfer_calibration_score_quantiles_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# Print
# ============================================================

cat(
  "\n===== TRANSFER CALIBRATION SUMMARY =====\n"
)

print(
  summary_df,
  row.names=FALSE
)

cat(
  "\n===== PER TAXONOMY =====\n"
)

print(
  per_taxonomy,
  row.names=FALSE
)

cat(
  "\n===== PER PROJECT =====\n"
)

print(
  per_project,
  row.names=FALSE
)

cat(
  "\n===== SCORE BINS =====\n"
)

print(
  score_bins,
  row.names=FALSE
)

cat(
  "\n===== SCORE / MARGIN QUANTILES =====\n"
)

print(
  quantiles_df,
  row.names=FALSE
)

# ============================================================
# Completion marker
#
# IMPORTANT:
# No high/medium/low transfer threshold is frozen here.
# ============================================================

writeLines(
  c(
    "PASS",
    "B/plasma full-transfer calibration v1",
    paste0(
      "n_pilot_cells=",
      nrow(md)
    ),
    paste0(
      "n_reference_cells=",
      length(reference_cells)
    ),
    paste0(
      "n_validation_cells=",
      length(validation_cells)
    ),
    paste0(
      "n_taxonomy_classes=",
      nrow(combo)
    ),
    paste0(
      "n_transfer_features=",
      length(transfer_features)
    ),
    paste0(
      "exact_taxonomy_accuracy=",
      sprintf(
        "%.6f",
        mean(result$exact_taxonomy_correct)
      )
    ),
    paste0(
      "core_identity_accuracy=",
      sprintf(
        "%.6f",
        mean(result$core_identity_correct)
      )
    ),
    paste0(
      "state_accuracy=",
      sprintf(
        "%.6f",
        mean(result$state_correct)
      )
    ),
    paste0(
      "annotation_scope_accuracy=",
      sprintf(
        "%.6f",
        mean(result$scope_correct)
      )
    ),
    "",
    "reference=query split is deterministic and stratified by project + frozen taxonomy",
    "native RNA counts only",
    "mapping method=Seurat FindTransferAnchors reduction=pcaproject + TransferData",
    "taxonomy transfer preserves frozen core_identity + state combinations",
    "transfer thresholds are intentionally NOT frozen in this step"
  ),
  done_file
)

cat(
  "\nPASS: B/plasma full-transfer calibration v1 completed\n"
)

