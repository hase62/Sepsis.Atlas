#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

suppressPackageStartupMessages({
  library(SeuratObject)
})

# ============================================================
# Paths
# ============================================================

in_rds <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_native_subclustering_v1_primary9997",
  "b_plasma_pilot_native_subclustering_v1.rds"
)

identity_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_underlying_identity_assignment_v1"
)

identity_summary_file <- file.path(
  identity_dir,
  "b_plasma_target_identity_summary_v1.tsv"
)

anchor_validation_file <- file.path(
  identity_dir,
  "b_plasma_anchor_loocv_validation_v1.tsv"
)

project_similarity_file <- file.path(
  identity_dir,
  "b_plasma_target_project_similarity_v1.tsv"
)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_pilot_annotation_freeze_v1"
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

done_file <- file.path(
  out_dir,
  "B_PLASMA_PILOT_ANNOTATION_FREEZE_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop(
    "Freeze already completed; refusing to overwrite: ",
    done_file
  )
}

for(f in c(
  in_rds,
  identity_summary_file,
  anchor_validation_file,
  project_similarity_file
)){
  if(!file.exists(f)){
    stop("Missing required input: ", f)
  }
}

# ============================================================
# Load pilot object
# ============================================================

bp <- readRDS(in_rds)

stopifnot(
  nrow(bp) == 38606,
  ncol(bp) == 9997
)

md <- bp[[]]

stopifnot(
  all(c(
    "project_id",
    "cluster_res_0p4"
  ) %in% names(md))
)

cluster <- as.character(
  md$cluster_res_0p4
)

stopifnot(
  setequal(
    unique(cluster),
    as.character(1:15)
  )
)

# ============================================================
# Frozen taxonomy
#
# Important decisions:
#
# cluster 5:
#   NOT frozen as Class_switched_memory_B_like.
#   Project-level identity is discordant:
#       GSE167363 -> Naive_B_like (n=691)
#       GSE216020 -> Class_switched_memory_B_like (n=19)
#       GSE220189 -> Class_switched_memory_B_like (n=83)
#   Therefore core identity remains unresolved.
#
# cluster 14:
#   Automated state-excluded similarity is essentially tied:
#       Memory_B_like                0.5176297
#       Class_switched_memory_B_like 0.5042927
#       margin                       0.0133370
#   Native RNA markers show clear IGHG1/IGHG2/IGHG3
#   class-switch evidence.
#   Therefore manually freeze as
#   Class_switched_memory_B_like.
# ============================================================

taxonomy <- data.frame(
  cluster = as.character(1:15),

  core_identity = c(
    "Naive_B_like",                       # 1
    "Memory_B_like",                      # 2
    "Class_switched_memory_B_like",       # 3
    "Naive_B_like",                       # 4
    "B_cell_unresolved",                  # 5
    "Naive_B_like",                       # 6
    "Naive_B_like",                       # 7
    "Naive_B_like",                       # 8
    "Naive_B_like",                       # 9
    "Naive_B_like",                       # 10
    "Naive_B_like",                       # 11
    "Naive_B_like",                       # 12
    "Deferred_non_B_platelet_like",       # 13
    "Class_switched_memory_B_like",       # 14
    "Deferred_non_B_T_cell_like"          # 15
  ),

  state = c(
    "none",                               # 1
    "none",                               # 2
    "none",                               # 3
    "none",                               # 4
    "activated_atypical_B_candidate",     # 5
    "interferon_stimulated",              # 6
    "none",                               # 7
    "early_activation",                   # 8
    "interferon_stimulated",              # 9
    "secretory_differentiation_candidate",# 10
    "activation_stress",                  # 11
    "unresolved_state",                   # 12
    "not_applicable",                     # 13
    "interferon_stimulated",              # 14
    "not_applicable"                      # 15
  ),

  confidence = c(
    "high",       # 1
    "high",       # 2
    "high",       # 3
    "high",       # 4
    "low_medium", # 5
    "medium",     # 6
    "high",       # 7
    "medium",     # 8
    "medium",     # 9
    "low_medium", # 10
    "medium",     # 11
    "low",        # 12
    "high",       # 13
    "low_medium", # 14
    "high"        # 15
  ),

  annotation_scope = c(
    "B",              # 1
    "B",              # 2
    "B",              # 3
    "B",              # 4
    "B_unresolved",   # 5
    "B",              # 6
    "B",              # 7
    "B",              # 8
    "B",              # 9
    "B",              # 10
    "B",              # 11
    "B",              # 12
    "deferred_non_B", # 13
    "B",              # 14
    "deferred_non_B"  # 15
  ),

  decision_basis = c(
    "anchor_high_confidence",
    "anchor_high_confidence",
    "anchor_high_confidence",
    "state_excluded_identity_strong_positive_control",
    "project_discordant_identity_keep_core_unresolved",
    "state_excluded_identity_plus_IFN_state",
    "state_excluded_identity_strong_positive_control",
    "state_excluded_identity_plus_activation_state",
    "state_excluded_identity_plus_IFN_state",
    "state_excluded_naive_identity_plus_secretory_candidate",
    "state_excluded_identity_plus_activation_stress",
    "state_excluded_naive_identity_state_unresolved",
    "non_B_platelet_marker_program",
    "native_IGHG_class_switch_evidence_overrides_similarity_near_tie",
    "non_B_T_cell_marker_program"
  ),

  stringsAsFactors=FALSE
)

stopifnot(
  nrow(taxonomy) == 15,
  !anyDuplicated(taxonomy$cluster)
)

# ============================================================
# Verify frozen LOPO diagnostic
# ============================================================

v <- read.delim(
  anchor_validation_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

calls <- unique(
  v[,c(
    "true_identity",
    "project_id",
    "n_cells",
    "winner",
    "correct"
  )]
)

stopifnot(
  nrow(calls) == 12,
  sum(calls$correct) == 11
)

lopo_accuracy <- mean(
  calls$correct
)

failed <- calls[
  !calls$correct,
  ,
  drop=FALSE
]

stopifnot(
  nrow(failed) == 1,
  failed$true_identity == "Memory_B_like",
  failed$project_id == "GSE220189",
  failed$winner == "Class_switched_memory_B_like"
)

# ============================================================
# Verify critical project-level calls
# ============================================================

s <- read.delim(
  project_similarity_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

get_project_calls <- function(k){

  d <- s[
    as.character(s$cluster) == k,
    ,
    drop=FALSE
  ]

  unique(
    d[,c(
      "project_id",
      "n_cells",
      "winner",
      "winner_combined",
      "winner_margin"
    )]
  )
}

c5 <- get_project_calls("5")
c10 <- get_project_calls("10")
c14 <- get_project_calls("14")

stopifnot(
  nrow(c5) == 3,
  nrow(c10) == 1,
  nrow(c14) == 1
)

# cluster 5 must remain project-discordant.
stopifnot(
  c5$winner[
    c5$project_id == "GSE167363"
  ] == "Naive_B_like",

  c5$winner[
    c5$project_id == "GSE216020"
  ] == "Class_switched_memory_B_like",

  c5$winner[
    c5$project_id == "GSE220189"
  ] == "Class_switched_memory_B_like"
)

# cluster 10 state-excluded profile must be naive-like.
stopifnot(
  c10$project_id == "GSE151263",
  c10$winner == "Naive_B_like"
)

# cluster 14 similarity must remain weak/near-tied.
stopifnot(
  c14$project_id == "GSE151263",
  c14$winner == "Memory_B_like",
  c14$winner_margin < 0.05
)

# ============================================================
# Map taxonomy to every pilot cell
# ============================================================

ii <- match(
  cluster,
  taxonomy$cluster
)

if(anyNA(ii)){
  stop(
    "Unmapped r0.4 cluster detected"
  )
}

ann <- data.frame(
  b_plasma_core_identity_v1 =
    taxonomy$core_identity[ii],

  b_plasma_state_v1 =
    taxonomy$state[ii],

  b_plasma_annotation_confidence_v1 =
    taxonomy$confidence[ii],

  b_plasma_annotation_scope_v1 =
    taxonomy$annotation_scope[ii],

  b_plasma_annotation_decision_basis_v1 =
    taxonomy$decision_basis[ii],

  stringsAsFactors=FALSE,
  row.names=rownames(md)
)

stopifnot(
  nrow(ann) == 9997,
  identical(
    rownames(ann),
    colnames(bp)
  ),
  !anyNA(ann)
)

bp <- AddMetaData(
  bp,
  metadata=ann
)

# ============================================================
# Freeze outputs
# ============================================================

taxonomy_file <- file.path(
  out_dir,
  "b_plasma_pilot_cluster_taxonomy_freeze_v1.tsv"
)

cell_metadata_file <- file.path(
  out_dir,
  "b_plasma_pilot_cell_metadata_freeze_v1.tsv"
)

summary_file <- file.path(
  out_dir,
  "b_plasma_pilot_annotation_counts_freeze_v1.tsv"
)

out_rds <- file.path(
  out_dir,
  "b_plasma_pilot_native_subclustering__annotation_freeze_v1.rds"
)

write.table(
  taxonomy,
  file=taxonomy_file,
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

cell_out <- data.frame(
  cell_id=rownames(bp[[]]),
  project_id=as.character(
    bp$project_id
  ),
  cluster_res_0p4=as.character(
    bp$cluster_res_0p4
  ),
  bp[[]][,c(
    "b_plasma_core_identity_v1",
    "b_plasma_state_v1",
    "b_plasma_annotation_confidence_v1",
    "b_plasma_annotation_scope_v1",
    "b_plasma_annotation_decision_basis_v1"
  )],
  stringsAsFactors=FALSE,
  check.names=FALSE
)

write.table(
  cell_out,
  file=cell_metadata_file,
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

summary_df <- as.data.frame(
  table(
    core_identity=
      bp$b_plasma_core_identity_v1,
    state=
      bp$b_plasma_state_v1
  ),
  stringsAsFactors=FALSE
)

summary_df <- summary_df[
  summary_df$Freq > 0,
  ,
  drop=FALSE
]

summary_df <- summary_df[
  order(
    summary_df$core_identity,
    summary_df$state
  ),
  ,
  drop=FALSE
]

write.table(
  summary_df,
  file=summary_file,
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

saveRDS(
  bp,
  file=out_rds,
  compress=FALSE
)

# ============================================================
# Final integrity checks
# ============================================================

stopifnot(
  file.exists(taxonomy_file),
  file.exists(cell_metadata_file),
  file.exists(summary_file),
  file.exists(out_rds)
)

bp_check <- readRDS(out_rds)

stopifnot(
  nrow(bp_check) == 38606,
  ncol(bp_check) == 9997,
  all(c(
    "b_plasma_core_identity_v1",
    "b_plasma_state_v1",
    "b_plasma_annotation_confidence_v1",
    "b_plasma_annotation_scope_v1",
    "b_plasma_annotation_decision_basis_v1"
  ) %in% names(bp_check[[]]))
)

# ============================================================
# Provenance / completion marker
# ============================================================

writeLines(
  c(
    "PASS",
    "B/plasma pilot annotation freeze v1",
    "n_cells=9997",
    "n_features=38606",
    paste0(
      "anchor_LOPO_accuracy=",
      sprintf("%.6f", lopo_accuracy)
    ),
    "anchor_LOPO_correct=11/12",
    "anchor_LOPO_only_failure=GSE220189 Memory_B_like -> Class_switched_memory_B_like",
    "",
    "critical_manual_decisions:",
    "cluster5=B_cell_unresolved + activated_atypical_B_candidate",
    "cluster5_reason=project-level underlying identity discordant; no single core identity frozen",
    "cluster10=Naive_B_like + secretory_differentiation_candidate",
    "cluster14=Class_switched_memory_B_like + interferon_stimulated",
    "cluster14_reason=native IGHG1/IGHG2/IGHG3 evidence overrides near-tied Memory/Class-switched similarity",
    "cluster13=Deferred_non_B_platelet_like",
    "cluster15=Deferred_non_B_T_cell_like",
    "",
    paste0(
      "source_rds=",
      in_rds
    ),
    paste0(
      "output_rds=",
      out_rds
    ),
    "",
    "This taxonomy is frozen for B/plasma pilot v1.",
    "Future full-atlas transfer must preserve core identity and orthogonal state as separate fields."
  ),
  done_file
)

cat(
  "\n===== FROZEN B/PLASMA TAXONOMY =====\n"
)

print(
  taxonomy[,c(
    "cluster",
    "core_identity",
    "state",
    "confidence"
  )],
  row.names=FALSE
)

cat(
  "\n===== CELL COUNTS =====\n"
)

print(
  summary_df,
  row.names=FALSE
)

cat(
  "\nPASS: B/plasma pilot annotation freeze v1 completed\n"
)

