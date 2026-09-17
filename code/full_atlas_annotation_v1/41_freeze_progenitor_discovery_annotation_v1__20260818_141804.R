#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 7L){
  stop(
    paste(
      "Usage:",
      "script <root> <tag> <cluster_dir>",
      "<marker_dir> <absolute_dir>",
      "<codetect_dir> <replication_dir>"
    )
  )
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]

cluster_dir <- normalizePath(args[[3]], mustWork=TRUE)
marker_dir <- normalizePath(args[[4]], mustWork=TRUE)
absolute_dir <- normalizePath(args[[5]], mustWork=TRUE)
codetect_dir <- normalizePath(args[[6]], mustWork=TRUE)
replication_dir <- normalizePath(args[[7]], mustWork=TRUE)

N_EXPECTED <- 26962L
DISCOVERY_PROJECT <- "GSE216007"
RESOLUTION <- "native_cluster_r0p4"

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "progenitor_discovery_annotation_freeze_v1__",
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
  "PROGENITOR_DISCOVERY_ANNOTATION_FREEZE_COMPLETE.ok"
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

# ============================================================
# Validate all evidence stages
# ============================================================

required_completion <- c(
  file.path(
    cluster_dir,
    "PROGENITOR_DOMINANT_PROJECT_NATIVE_CLUSTERING_COMPLETE.ok"
  ),
  file.path(
    marker_dir,
    "PROGENITOR_NATIVE_MARKER_AGGREGATION_COMPLETE.ok"
  ),
  file.path(
    absolute_dir,
    "PROGENITOR_NATIVE_ABSOLUTE_PROGRAM_AUDIT_COMPLETE.ok"
  ),
  file.path(
    codetect_dir,
    "PROGENITOR_CELL_LEVEL_PROGRAM_CODETECTION_COMPLETE.ok"
  ),
  file.path(
    replication_dir,
    "PROGENITOR_CROSS_PROJECT_PROGRAM_REPLICATION_COMPLETE.ok"
  )
)

stopifnot(
  all(file.exists(required_completion))
)

# ============================================================
# Inputs
# ============================================================

assignment_file <- file.path(
  cluster_dir,
  "progenitor_dominant_project_native_cluster_assignments_v1.tsv.gz"
)

marker_summary_file <- file.path(
  marker_dir,
  "progenitor_r0p4_cluster_marker_summary_v1.tsv"
)

codetect_files <- list.files(
  codetect_dir,
  pattern="^progenitor_cell_level_program_flags_v1__.*\\.tsv\\.gz$",
  full.names=TRUE
)

stopifnot(
  length(codetect_files) == 1L
)

codetect_file <- codetect_files[[1]]

replication_file <- file.path(
  replication_dir,
  "progenitor_program_cross_project_replication_summary_v1.tsv"
)

assign <- read_tsv(assignment_file)
marker_summary <- read_tsv(marker_summary_file)
flags <- read_tsv(codetect_file)
replication <- read_tsv(replication_file)

# ============================================================
# Basic invariants
# ============================================================

req_assign <- c(
  "global_cell",
  "project_id",
  "library_key",
  RESOLUTION
)

stopifnot(
  all(req_assign %in% names(assign)),
  nrow(assign) == N_EXPECTED,
  !anyDuplicated(assign$global_cell)
)

assign$global_cell <- as.character(assign$global_cell)
assign$project_id <- as.character(assign$project_id)
assign$library_key <- as.character(assign$library_key)
assign$cluster <- as.character(assign[[RESOLUTION]])

stopifnot(
  all(assign$project_id == DISCOVERY_PROJECT),
  setequal(
    unique(assign$cluster),
    as.character(0:13)
  )
)

expected_counts <- c(
  "0"=5883L,
  "1"=5627L,
  "2"=2704L,
  "3"=2533L,
  "4"=2027L,
  "5"=1862L,
  "6"=1836L,
  "7"=1222L,
  "8"=956L,
  "9"=679L,
  "10"=576L,
  "11"=526L,
  "12"=490L,
  "13"=41L
)

observed_counts <- table(assign$cluster)

stopifnot(
  identical(
    as.integer(
      observed_counts[names(expected_counts)]
    ),
    as.integer(expected_counts)
  )
)

stopifnot(
  nrow(flags) == N_EXPECTED,
  !anyDuplicated(flags$global_cell),
  "cycling_flag" %in% names(flags)
)

flags$global_cell <- as.character(
  flags$global_cell
)

# ============================================================
# Verify cross-project replication evidence
# ============================================================

stopifnot(
  all(
    c(
      "signal",
      "n_projects_supported",
      "replicated_ge2_projects"
    ) %in% names(replication)
  )
)

replication$signal <- as.character(
  replication$signal
)

rep_bool <- as.logical(
  replication$replicated_ge2_projects
)

names(rep_bool) <- replication$signal

rep_n <- as.integer(
  replication$n_projects_supported
)

names(rep_n) <- replication$signal

required_replication <- c(
  "strict_HSC",
  "primitive_HSC_low_commitment",
  "lymphoid_priming",
  "lymphoid_primed_nonB",
  "granulocytic",
  "megakaryocytic",
  "MegE",
  "meg_only",
  "eos_baso_mast",
  "cycling",
  "B_commitment"
)

stopifnot(
  all(required_replication %in% names(rep_bool))
)

# Expected evidence from frozen Step40.
stopifnot(
  rep_bool[["strict_HSC"]],
  rep_bool[["primitive_HSC_low_commitment"]],
  rep_bool[["lymphoid_priming"]],
  rep_bool[["lymphoid_primed_nonB"]],
  rep_bool[["granulocytic"]],
  rep_bool[["megakaryocytic"]],
  rep_bool[["MegE"]],
  rep_bool[["meg_only"]],
  rep_bool[["eos_baso_mast"]],
  rep_bool[["cycling"]],
  !rep_bool[["B_commitment"]]
)

stopifnot(
  rep_n[["primitive_HSC_low_commitment"]] == 7L,
  rep_n[["lymphoid_primed_nonB"]] == 6L,
  rep_n[["granulocytic"]] == 2L,
  rep_n[["megakaryocytic"]] == 7L,
  rep_n[["MegE"]] == 6L,
  rep_n[["eos_baso_mast"]] == 6L,
  rep_n[["cycling"]] == 3L,
  rep_n[["B_commitment"]] == 0L
)

# ============================================================
# Frozen discovery decision
#
# Core identity is separated from lineage-priming axis.
# Library-confounded clusters are not promoted to
# independent biological subtypes.
# ============================================================

decision <- data.frame(

  cluster=as.character(0:13),

  core_identity=c(
    "Progenitor_like",                      # c0
    "Progenitor_like",                      # c1
    "Progenitor_like",                      # c2
    "Progenitor_like",                      # c3
    "Progenitor_like",                      # c4
    "Progenitor_like",                      # c5
    "Progenitor_like",                      # c6
    "Progenitor_like",                      # c7
    "Progenitor_like",                      # c8
    "Deferred_non_progenitor_lymphoid_like",# c9
    "Progenitor_like",                      # c10
    "Progenitor_like",                      # c11
    "Progenitor_like",                      # c12
    "Progenitor_like"                       # c13
  ),

  lineage_priming_axis=c(
    "megakaryocytic_priming",            # c0
    "lymphoid_priming",                   # c1
    "megakaryocytic_erythroid_priming",  # c2
    "lymphoid_priming",                   # c3
    "granulocytic_priming",               # c4
    "unresolved",                         # c5
    "unresolved",                         # c6
    "primitive_HSC_like",                 # c7
    "eos_baso_mast_priming",              # c8
    "not_applicable",                     # c9
    "unresolved",                         # c10
    "unresolved",                         # c11
    "primitive_HSC_like",                 # c12
    "unresolved"                          # c13
  ),

  annotation_confidence=c(
    "medium_high", # c0
    "high",        # c1
    "high",        # c2
    "high",        # c3
    "medium_high", # c4
    "low_medium",  # c5
    "low_medium",  # c6
    "high",        # c7
    "high",        # c8
    "medium",      # c9
    "low_medium",  # c10
    "low_medium",  # c11
    "high",        # c12
    "low_medium"   # c13
  ),

  evidence_flag=c(
    "megakaryocytic_program_replicated;within_cluster_eos_baso_overlap",
    "lymphoid_priming_replicated",
    "MegE_program_replicated;cycling_enriched_subset",
    "lymphoid_priming_replicated",
    "granulocytic_program_replicated_two_projects",
    "library_confounded;primitive_HSC_like_signal;axis_not_promoted",
    "library_confounded;primitive_HSC_like_signal;axis_not_promoted",
    "primitive_HSC_low_commitment_replicated",
    "eos_baso_mast_program_replicated",
    "low_HSC;lymphoid_TNK_like_signal;core_deferred",
    "library_confounded;primitive_HSC_like_signal;axis_not_promoted",
    "library_confounded;primitive_HSC_like_signal;axis_not_promoted",
    "primitive_HSC_low_commitment_replicated",
    "B_commitment_strong_within_cluster;tiny;not_cross_project_replicated;axis_not_promoted"
  ),

  stringsAsFactors=FALSE
)

# ============================================================
# Attach discovery technical evidence
# ============================================================

marker_summary$cluster <- as.character(
  marker_summary$cluster
)

mi <- match(
  decision$cluster,
  marker_summary$cluster
)

stopifnot(
  !anyNA(mi)
)

decision$n_cells <-
  marker_summary$n_cells[mi]

decision$n_libraries_present <-
  marker_summary$n_libraries_present[mi]

decision$max_library_fraction <-
  marker_summary$max_library_fraction[mi]

decision$effective_n_libraries <-
  marker_summary$effective_n_libraries[mi]

decision$discovery_evidence_scope <-
  marker_summary$evidence_scope[mi]

stopifnot(
  all(
    decision$n_cells ==
      as.integer(
        expected_counts[decision$cluster]
      )
  )
)

# ============================================================
# Merge cell-level state
# ============================================================

fi <- match(
  assign$global_cell,
  flags$global_cell
)

stopifnot(
  !anyNA(fi)
)

cycling_flag <- as.logical(
  flags$cycling_flag[fi]
)

di <- match(
  assign$cluster,
  decision$cluster
)

stopifnot(
  !anyNA(di)
)

cell_annotation <- data.frame(
  global_cell=assign$global_cell,
  project_id=assign$project_id,
  library_key=assign$library_key,
  discovery_cluster=assign$cluster,

  core_identity=
    decision$core_identity[di],

  lineage_priming_axis=
    decision$lineage_priming_axis[di],

  state=ifelse(
    decision$core_identity[di] ==
      "Deferred_non_progenitor_lymphoid_like",
    "not_applicable",
    ifelse(
      cycling_flag,
      "cycling",
      "none"
    )
  ),

  annotation_confidence=
    decision$annotation_confidence[di],

  evidence_flag=
    decision$evidence_flag[di],

  annotation_source=
    "GSE216007_native_r0p4_discovery_reference",

  stringsAsFactors=FALSE
)

stopifnot(
  nrow(cell_annotation) == N_EXPECTED,
  !anyDuplicated(
    cell_annotation$global_cell
  )
)

# c9 must be completely deferred.
stopifnot(
  all(
    cell_annotation$core_identity[
      cell_annotation$discovery_cluster == "9"
    ] ==
      "Deferred_non_progenitor_lymphoid_like"
  ),
  all(
    cell_annotation$lineage_priming_axis[
      cell_annotation$discovery_cluster == "9"
    ] ==
      "not_applicable"
  ),
  all(
    cell_annotation$state[
      cell_annotation$discovery_cluster == "9"
    ] ==
      "not_applicable"
  )
)

# B commitment is intentionally NOT a canonical axis.
stopifnot(
  !"B_commitment" %in%
    unique(
      cell_annotation$lineage_priming_axis
    )
)

# ============================================================
# Summary tables
# ============================================================

core_summary <- as.data.frame(
  table(
    core_identity=
      cell_annotation$core_identity
  ),
  stringsAsFactors=FALSE
)

axis_summary <- as.data.frame(
  table(
    lineage_priming_axis=
      cell_annotation$lineage_priming_axis
  ),
  stringsAsFactors=FALSE
)

state_summary <- as.data.frame(
  table(
    state=
      cell_annotation$state
  ),
  stringsAsFactors=FALSE
)

joint_summary <- as.data.frame(
  table(
    cluster=
      cell_annotation$discovery_cluster,
    core_identity=
      cell_annotation$core_identity,
    lineage_priming_axis=
      cell_annotation$lineage_priming_axis,
    state=
      cell_annotation$state
  ),
  stringsAsFactors=FALSE
)

joint_summary <- joint_summary[
  joint_summary$Freq > 0,
  ,
  drop=FALSE
]

# ============================================================
# Write
# ============================================================

decision_file <- file.path(
  out_dir,
  "progenitor_discovery_cluster_decision_v1.tsv"
)

write.table(
  decision,
  decision_file,
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

cell_file <- file.path(
  out_dir,
  paste0(
    "progenitor_discovery_cell_annotation_v1__",
    tag,
    ".tsv.gz"
  )
)

write_gz_tsv(
  cell_annotation,
  cell_file
)

write.table(
  core_summary,
  file.path(
    out_dir,
    "progenitor_discovery_core_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  axis_summary,
  file.path(
    out_dir,
    "progenitor_discovery_axis_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  state_summary,
  file.path(
    out_dir,
    "progenitor_discovery_state_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  joint_summary,
  file.path(
    out_dir,
    "progenitor_discovery_joint_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Provenance manifest
# ============================================================

source_manifest <- data.frame(
  stage=c(
    "Step38_clustering",
    "Step39_native_markers",
    "Step39b_absolute_programs",
    "Step39c_cell_codetection",
    "Step40_cross_project_replication"
  ),
  path=c(
    cluster_dir,
    marker_dir,
    absolute_dir,
    codetect_dir,
    replication_dir
  ),
  stringsAsFactors=FALSE
)

write.table(
  source_manifest,
  file.path(
    out_dir,
    "progenitor_discovery_annotation_source_manifest_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Freeze RDS
# ============================================================

freeze <- list(
  annotation_version=
    "Progenitor_discovery_annotation_v1",

  discovery_project=
    DISCOVERY_PROJECT,

  n_cells=
    N_EXPECTED,

  backbone_resolution=
    RESOLUTION,

  decision=
    decision,

  cell_annotation=
    cell_annotation,

  core_summary=
    core_summary,

  axis_summary=
    axis_summary,

  state_summary=
    state_summary,

  source_manifest=
    source_manifest,

  cross_project_replication=
    replication
)

rds_file <- file.path(
  out_dir,
  paste0(
    "progenitor_discovery_annotation_freeze_v1__",
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
  reread$n_cells == N_EXPECTED,
  identical(
    reread$cell_annotation$global_cell,
    cell_annotation$global_cell
  ),
  identical(
    reread$decision$cluster,
    decision$cluster
  )
)

# ============================================================
# Gzip re-read
# ============================================================

cell_check <- read_tsv(
  cell_file
)

stopifnot(
  nrow(cell_check) == N_EXPECTED,
  identical(
    as.character(cell_check$global_cell),
    as.character(cell_annotation$global_cell)
  ),
  identical(
    as.character(cell_check$core_identity),
    as.character(cell_annotation$core_identity)
  ),
  identical(
    as.character(cell_check$lineage_priming_axis),
    as.character(cell_annotation$lineage_priming_axis)
  ),
  identical(
    as.character(cell_check$state),
    as.character(cell_annotation$state)
  )
)

# ============================================================
# SHA256
# ============================================================

sha_files <- c(
  decision_file,
  cell_file,
  file.path(
    out_dir,
    "progenitor_discovery_core_summary_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_discovery_axis_summary_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_discovery_state_summary_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_discovery_joint_summary_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_discovery_annotation_source_manifest_v1.tsv"
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
  stop("sha256 creation failed")
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
# Completion LAST
# ============================================================

writeLines(
  c(
    "PASS",
    "Progenitor discovery annotation freeze v1",
    "annotation_version=Progenitor_discovery_annotation_v1",
    "discovery_project=GSE216007",
    "n_discovery_cells=26962",
    "backbone=native_cluster_r0p4",
    "",
    "core identity separated from lineage-priming axis",
    "cycling assigned as orthogonal cell-level state",
    "",
    "accepted replicated axes:",
    "primitive_HSC_like",
    "lymphoid_priming",
    "granulocytic_priming",
    "megakaryocytic_priming",
    "megakaryocytic_erythroid_priming",
    "eos_baso_mast_priming",
    "",
    "B commitment not promoted because cross-project replication failed",
    "c5,c6,c10,c11 lineage axis unresolved because library-confounded",
    "c13 lineage axis unresolved because tiny and unreplicated",
    "c9 deferred as non-progenitor lymphoid-like",
    "",
    "no classical HSC/MPP/CMP/GMP taxonomy forced",
    "condition not used",
    "native SoupX evidence",
    "cellranger_count not accessed",
    "",
    "gzip_integrity=PASS",
    "cell_annotation_reread_count_order=PASS",
    "freeze_rds_reread=PASS",
    "sha256_manifest=PASS"
  ),
  done_file
)

cat("\n===== DECISION =====\n")
print(decision, row.names=FALSE)

cat("\n===== CORE =====\n")
print(core_summary, row.names=FALSE)

cat("\n===== AXIS =====\n")
print(axis_summary, row.names=FALSE)

cat("\n===== STATE =====\n")
print(state_summary, row.names=FALSE)

cat("\nOUT_DIR=", out_dir, "\n", sep="")
cat("\nPASS: Progenitor discovery annotation freeze completed\n")
