#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

HIGH_THRESHOLD <- 0.80
ACCEPT_THRESHOLD <- 0.70

source(
  file.path(
    root,
    "full_atlas_primary_integration_v1",
    "00_common.R"
  )
)

# ============================================================
# Paths
# ============================================================

old_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_full_primary_transfer_v1"
)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_full_primary_transfer_v1__repaired_20260814_1430"
)

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

raw_pred_file <- file.path(
  old_dir,
  "b_plasma_query_transferdata_raw_v1.rds"
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

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

done_file <- file.path(
  out_dir,
  "B_PLASMA_FULL_PRIMARY_TRANSFER_REPAIR_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop(
    "Repair already completed; refusing to overwrite: ",
    done_file
  )
}

for(f in c(
  pilot_rds,
  taxonomy_file,
  raw_pred_file,
  lib_file
)){
  if(!file.exists(f)){
    stop("Missing input: ", f)
  }
}

# ============================================================
# Explicit gzip writer
#
# Critical:
# gzfile(open="wt") passed inline to write.table() can remain
# open because write.table did not create the connection.
# Always explicitly close it before proceeding.
# ============================================================

write_gz_tsv <- function(x, path){

  con <- gzfile(
    path,
    open="wt"
  )

  tryCatch(
    {
      write.table(
        x,
        file=con,
        sep="\t",
        quote=TRUE,
        qmethod="double",
        row.names=FALSE
      )
    },
    finally={
      close(con)
    }
  )

  status <- system2(
    "gzip",
    c(
      "-t",
      path
    ),
    stdout=FALSE,
    stderr=FALSE
  )

  if(status != 0L){
    stop(
      "gzip integrity check failed: ",
      path
    )
  }

  invisible(path)
}

# ============================================================
# Load frozen pilot
# ============================================================

pilot <- readRDS(
  pilot_rds
)

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
  all(
    required_pilot_md %in%
      names(pilot_md)
  )
)

# ============================================================
# Frozen taxonomy map
# ============================================================

tax <- read.delim(
  taxonomy_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

tax$cluster <- as.character(
  tax$cluster
)

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

stopifnot(
  nrow(combo) == 12L
)

tax$taxonomy_id <- combo$taxonomy_id[
  match(
    tax$taxonomy_key,
    combo$taxonomy_key
  )
]

# Biological confidence for each frozen taxonomy.
combo$frozen_taxonomy_confidence <- vapply(
  combo$taxonomy_id,
  function(k){

    z <- unique(
      as.character(
        tax$confidence[
          tax$taxonomy_id == k
        ]
      )
    )

    if(length(z) != 1L){
      stop(
        "Non-unique biological confidence for ",
        k,
        ": ",
        paste(z, collapse=",")
      )
    }

    z
  },
  character(1)
)

# ============================================================
# Frozen pilot taxonomy IDs
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

pilot_taxonomy_id <-
  tax$taxonomy_id[ii]

names(pilot_taxonomy_id) <-
  rownames(pilot_md)

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

# ============================================================
# Reconstruct canonical full B/plasma metadata only
#
# No expression matrices are needed for this repair.
# ============================================================

lib <- read_tsv(
  lib_file
)

full_md_list <- list()

n_full <- 0L
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
    stop(
      "Missing transfer file: ",
      tr_file
    )
  }

  tr <- read_tsv(
    tr_file
  )

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

  global_cells <- as.character(
    tr$global_cell
  )

  if(anyDuplicated(global_cells)){
    stop(
      "Duplicated global_cell within ",
      tag
    )
  }

  rownames(tr) <-
    global_cells

  full_md_list[[
    length(full_md_list)+1L
  ]] <- tr

  n_full <- n_full + nrow(tr)
  n_bp_libraries <-
    n_bp_libraries + 1L
}

full_md <- do.call(
  rbind,
  full_md_list
)

rm(full_md_list)

stopifnot(
  n_full == 45955L,
  n_bp_libraries == 134L,
  nrow(full_md) == 45955L,
  !anyDuplicated(
    rownames(full_md)
  ),
  all(
    pilot_cells %in%
      rownames(full_md)
  )
)

query_cells <- rownames(full_md)[
  !rownames(full_md) %in%
    pilot_cells
]

stopifnot(
  length(query_cells) == 35958L
)

# ============================================================
# Load RAW TransferData result
#
# This is the intact checkpoint written before the corrupted
# text exports.
# ============================================================

pred <- readRDS(
  raw_pred_file
)

cat(
  "raw prediction rows =",
  nrow(pred),
  "\n"
)

stopifnot(
  nrow(pred) == 35958L,
  !anyDuplicated(
    rownames(pred)
  ),
  setequal(
    rownames(pred),
    query_cells
  )
)

# Canonical query order.
pred <- pred[
  query_cells,
  ,
  drop=FALSE
]

stopifnot(
  identical(
    rownames(pred),
    query_cells
  )
)

# ============================================================
# Recompute prediction margin
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

top_two <- t(
  apply(
    score_mat,
    1,
    function(z){

      s <- sort(
        as.numeric(z),
        decreasing=TRUE
      )

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
# Decode query
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

accepted <-
  prediction_score >=
  ACCEPT_THRESHOLD

query_result <- data.frame(
  cell_id=query_cells,

  project_id=as.character(
    full_md[
      query_cells,
      "project_id"
    ]
  ),

  library_key=as.character(
    full_md[
      query_cells,
      "library_key"
    ]
  ),

  b_plasma_raw_predicted_taxonomy_id_v1=
    pred_id,

  b_plasma_raw_predicted_core_identity_v1=
    combo$core_identity[pred_idx],

  b_plasma_raw_predicted_state_v1=
    combo$state[pred_idx],

  b_plasma_raw_predicted_annotation_scope_v1=
    combo$annotation_scope[pred_idx],

  b_plasma_raw_predicted_frozen_taxonomy_confidence_v1=
    combo$frozen_taxonomy_confidence[pred_idx],

  b_plasma_prediction_score_v1=
    prediction_score,

  b_plasma_prediction_margin_v1=
    as.numeric(
      prediction_margin
    ),

  b_plasma_transfer_confidence_v1=
    transfer_confidence,

  b_plasma_transfer_accepted_v1=
    accepted,

  stringsAsFactors=FALSE
)

query_result$b_plasma_final_taxonomy_id_v1 <-
  ifelse(
    accepted,
    query_result$b_plasma_raw_predicted_taxonomy_id_v1,
    "TRANSFER_UNRESOLVED"
  )

query_result$b_plasma_final_core_identity_v1 <-
  ifelse(
    accepted,
    query_result$b_plasma_raw_predicted_core_identity_v1,
    "B_plasma_transfer_unresolved"
  )

query_result$b_plasma_final_state_v1 <-
  ifelse(
    accepted,
    query_result$b_plasma_raw_predicted_state_v1,
    "unresolved"
  )

query_result$b_plasma_final_annotation_scope_v1 <-
  ifelse(
    accepted,
    query_result$b_plasma_raw_predicted_annotation_scope_v1,
    "unresolved"
  )

query_result$b_plasma_final_frozen_taxonomy_confidence_v1 <-
  ifelse(
    accepted,
    query_result$b_plasma_raw_predicted_frozen_taxonomy_confidence_v1,
    "unresolved"
  )

query_result$b_plasma_annotation_source_v1 <-
  "full_transfer_v1"

rownames(query_result) <-
  query_result$cell_id

# ============================================================
# Frozen reference rows
# ============================================================

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
    as.character(
      full_md[
        pilot_cells,
        "library_key"
      ]
    )
  },

  b_plasma_raw_predicted_taxonomy_id_v1=
    NA_character_,

  b_plasma_raw_predicted_core_identity_v1=
    NA_character_,

  b_plasma_raw_predicted_state_v1=
    NA_character_,

  b_plasma_raw_predicted_annotation_scope_v1=
    NA_character_,

  b_plasma_raw_predicted_frozen_taxonomy_confidence_v1=
    NA_character_,

  b_plasma_prediction_score_v1=
    NA_real_,

  b_plasma_prediction_margin_v1=
    NA_real_,

  b_plasma_transfer_confidence_v1=
    "reference_frozen",

  b_plasma_transfer_accepted_v1=
    TRUE,

  b_plasma_final_taxonomy_id_v1=
    pilot_taxonomy_id,

  b_plasma_final_core_identity_v1=
    as.character(
      pilot_md$b_plasma_core_identity_v1
    ),

  b_plasma_final_state_v1=
    as.character(
      pilot_md$b_plasma_state_v1
    ),

  b_plasma_final_annotation_scope_v1=
    as.character(
      pilot_md$b_plasma_annotation_scope_v1
    ),

  b_plasma_final_frozen_taxonomy_confidence_v1=
    as.character(
      pilot_md$b_plasma_annotation_confidence_v1
    ),

  b_plasma_annotation_source_v1=
    "pilot_freeze_v1",

  stringsAsFactors=FALSE
)

rownames(pilot_result) <-
  pilot_result$cell_id

# ============================================================
# Combine
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

# ============================================================
# Build full annotation
#
# IMPORTANT:
# Existing upstream column `transfer_confidence` is preserved.
# New B/plasma fields all use `b_plasma_*_v1`, so there are
# no duplicate column names.
# ============================================================

bp_cols <- setdiff(
  names(result_all),
  c(
    "cell_id",
    "project_id",
    "library_key"
  )
)

stopifnot(
  !any(
    bp_cols %in%
      names(full_md)
  )
)

full_annotation <- cbind(
  full_md,
  result_all[
    ,
    bp_cols,
    drop=FALSE
  ]
)

full_annotation$cell_id <-
  rownames(full_annotation)

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

stopifnot(
  nrow(full_annotation) == 45955L,
  !anyDuplicated(
    full_annotation$cell_id
  ),
  !anyDuplicated(
    names(full_annotation)
  )
)

# ============================================================
# Save lossless checkpoints BEFORE text exports
# ============================================================

saveRDS(
  query_result,
  file=file.path(
    out_dir,
    "b_plasma_full_query_predictions_repaired_v1.rds"
  ),
  compress=FALSE
)

saveRDS(
  full_annotation,
  file=file.path(
    out_dir,
    "b_plasma_full_primary_annotation_repaired_v1.rds"
  ),
  compress=FALSE
)

# ============================================================
# Write gzip exports with EXPLICIT CLOSE
# ============================================================

query_tsv <- file.path(
  out_dir,
  "b_plasma_full_query_predictions_repaired_v1.tsv.gz"
)

full_tsv <- file.path(
  out_dir,
  "b_plasma_full_primary_annotation_repaired_v1.tsv.gz"
)

write_gz_tsv(
  query_result,
  query_tsv
)

write_gz_tsv(
  full_annotation,
  full_tsv
)

# ============================================================
# Re-read exported files: hard integrity check
# ============================================================

q_check <- read.delim(
  gzfile(query_tsv),
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

a_check <- read.delim(
  gzfile(full_tsv),
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

stopifnot(
  nrow(q_check) == 35958L,
  nrow(a_check) == 45955L,
  !anyDuplicated(q_check$cell_id),
  !anyDuplicated(a_check$cell_id),
  identical(
    q_check$cell_id,
    query_result$cell_id
  ),
  identical(
    a_check$cell_id,
    full_annotation$cell_id
  )
)

# Critical confidence consistency.
score_conf_check <- ifelse(
  q_check$b_plasma_prediction_score_v1 >= HIGH_THRESHOLD,
  "high",
  ifelse(
    q_check$b_plasma_prediction_score_v1 >= ACCEPT_THRESHOLD,
    "medium",
    "unresolved"
  )
)

stopifnot(
  identical(
    as.character(
      q_check$b_plasma_transfer_confidence_v1
    ),
    as.character(
      score_conf_check
    )
  )
)

# ============================================================
# Summaries
# ============================================================

query_summary <- as.data.frame(
  table(
    b_plasma_transfer_confidence_v1=
      query_result$b_plasma_transfer_confidence_v1
  ),
  stringsAsFactors=FALSE
)

taxonomy_summary <- as.data.frame(
  table(
    b_plasma_annotation_source_v1=
      result_all$b_plasma_annotation_source_v1,
    b_plasma_final_taxonomy_id_v1=
      result_all$b_plasma_final_taxonomy_id_v1
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
          d$b_plasma_transfer_confidence_v1 ==
            "high"
        ),
        n_medium=sum(
          d$b_plasma_transfer_confidence_v1 ==
            "medium"
        ),
        n_unresolved=sum(
          d$b_plasma_transfer_confidence_v1 ==
            "unresolved"
        ),
        accepted_fraction=mean(
          d$b_plasma_transfer_accepted_v1
        ),
        median_prediction_score=median(
          d$b_plasma_prediction_score_v1
        ),
        median_prediction_margin=median(
          d$b_plasma_prediction_margin_v1
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
    "b_plasma_full_query_confidence_summary_repaired_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  taxonomy_summary,
  file=file.path(
    out_dir,
    "b_plasma_full_taxonomy_counts_repaired_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  project_summary,
  file=file.path(
    out_dir,
    "b_plasma_full_transfer_by_project_repaired_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

n_high <- sum(
  query_result$b_plasma_transfer_confidence_v1 ==
    "high"
)

n_medium <- sum(
  query_result$b_plasma_transfer_confidence_v1 ==
    "medium"
)

n_unresolved <- sum(
  query_result$b_plasma_transfer_confidence_v1 ==
    "unresolved"
)

stopifnot(
  n_high + n_medium + n_unresolved ==
    35958L
)

# ============================================================
# Completion marker only AFTER gzip verification
# ============================================================

writeLines(
  c(
    "PASS",
    "B/plasma full primary transfer output repair v1",
    "source_mapping=b_plasma_query_transferdata_raw_v1.rds",
    "mapping_recomputed=NO",
    "n_full_cells=45955",
    "n_reference_pilot_cells=9997",
    "n_query_cells=35958",
    "n_projects=9",
    "n_B_plasma_libraries=134",
    "",
    "repair reason:",
    "original gzip TSV outputs were truncated because gzip connections were not explicitly closed",
    "original full annotation also contained duplicate transfer_confidence column names",
    "",
    "repair actions:",
    "raw TransferData checkpoint reused",
    "all B/plasma-specific columns renamed with b_plasma_*_v1 prefix",
    "gzip connections explicitly closed",
    "gzip -t passed for both exported TSV.gz files",
    "both exported files re-read and exact row counts verified",
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
          query_result$b_plasma_transfer_accepted_v1
        )
      )
    )
  ),
  done_file
)

cat(
  "\n===== REPAIRED QUERY CONFIDENCE =====\n"
)

print(
  query_summary,
  row.names=FALSE
)

cat(
  "\n===== REPAIRED FULL TAXONOMY COUNTS =====\n"
)

print(
  taxonomy_summary,
  row.names=FALSE
)

cat(
  "\n===== REPAIRED PROJECT SUMMARY =====\n"
)

print(
  project_summary,
  row.names=FALSE
)

cat(
  "\nPASS: B/plasma output repair completed\n"
)

