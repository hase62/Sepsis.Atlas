#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

base <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1"
)

cl_rds <- file.path(
  base,
  "integrated_compartment_clustering",
  "Monocyte_DC_combined",
  "integrated_clustering_v1.rds"
)

meta_file <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "scmerge2_primary_runs",
  "Monocyte_DC_combined__full_v1_cells.tsv.gz"
)

audit_file <- file.path(
  base,
  "monodc_project_balanced_marker_audit_v1",
  "monodc_r0p6_project_balanced_summary_v1.tsv"
)

out_dir <- file.path(
  base,
  "monodc_annotation_freeze_v1"
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

for(f in c(cl_rds, meta_file, audit_file)){
  if(!file.exists(f)){
    stop("Missing input: ", f)
  }
}

# ============================================================
# Load authoritative clustering
# ============================================================

x <- readRDS(cl_rds)

cells <- as.character(x$cells)

stopifnot(
  "cluster_res_0p6" %in%
    names(x$clusters)
)

r06 <- as.character(
  x$clusters$cluster_res_0p6
)

stopifnot(
  length(cells) == 192040,
  length(r06) == length(cells)
)

# ============================================================
# Load authoritative cell metadata
# ============================================================

meta <- read.delim(
  gzfile(meta_file),
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

required <- c(
  "global_cell",
  "project_id",
  "library_key",
  "clinical_domain_std",
  "condition_binary",
  "integration_celltype_full_v1",
  "transfer_confidence"
)

miss <- setdiff(
  required,
  names(meta)
)

if(length(miss)){
  stop(
    "Missing metadata columns: ",
    paste(miss, collapse=", ")
  )
}

m <- match(
  cells,
  as.character(meta$global_cell)
)

if(anyNA(m)){
  stop(
    "Cell metadata mismatch: ",
    sum(is.na(m))
  )
}

meta <- meta[m,,drop=FALSE]

stopifnot(
  identical(
    as.character(meta$global_cell),
    cells
  )
)

backbone <- as.character(
  meta$integration_celltype_full_v1
)

allowed <- c(
  "Classical_monocyte_like",
  "Nonclassical_monocyte_like",
  "DC2_monocyte_boundary"
)

if(any(!backbone %in% allowed)){
  stop(
    "Unexpected backbone label(s): ",
    paste(
      sort(unique(
        backbone[
          !backbone %in% allowed
        ]
      )),
      collapse=", "
    )
  )
}

# ============================================================
# Start from transferred biological backbone
# ============================================================

core <- backbone

core_basis <- rep(
  "transferred_biological_backbone",
  length(cells)
)

core_confidence <- as.character(
  meta$transfer_confidence
)

# Normalize unexpected missing confidence.
core_confidence[
  is.na(core_confidence) |
  !nzchar(core_confidence)
] <- "medium"

# ============================================================
# Cluster 14:
# reproducible FCER1A/CD1C cDC2 identity
#
# Only refine cells whose transferred backbone already lies
# on the DC2/monocyte boundary.
# ============================================================

ii_cdc2 <- (
  r06 == "14" &
  backbone == "DC2_monocyte_boundary"
)

core[ii_cdc2] <-
  "CD1C_FCER1A_cDC2_like"

core_basis[ii_cdc2] <-
  "r0.6_cluster14_project_balanced_native_markers"

core_confidence[ii_cdc2] <-
  "high"

# ============================================================
# Orthogonal transcriptional states
# ============================================================

state <- rep(
  "none",
  length(cells)
)

state_confidence <- rep(
  "not_applicable",
  length(cells)
)

state_basis <- rep(
  "none",
  length(cells)
)

# Cluster 21:
# strong canonical interferon-stimulated program,
# reproducible in two eligible projects.
ii <- r06 == "21"

state[ii] <-
  "interferon_stimulated"

state_confidence[ii] <-
  "medium"

state_basis[ii] <-
  "r0.6_cluster21_two_project_ISG_program"

# Cluster 30:
# reproducible across four eligible projects.
# Keep descriptive marker-defined name.
ii <- r06 == "30"

state[ii] <-
  "LPAR1_RFX3_high"

state_confidence[ii] <-
  "high"

state_basis[ii] <-
  "r0.6_cluster30_four_project_marker_program"

# ============================================================
# Candidate state:
# recorded for provenance, NOT promoted to final state taxonomy
# ============================================================

candidate_state <- rep(
  "none",
  length(cells)
)

candidate_confidence <- rep(
  "not_applicable",
  length(cells)
)

ii <- r06 == "3"

candidate_state[ii] <-
  "RETN_PADI4_inflammatory_candidate"

candidate_confidence[ii] <-
  "low"

# ============================================================
# Cell-level output
# ============================================================

out <- data.frame(
  global_cell=cells,
  project_id=meta$project_id,
  library_key=meta$library_key,
  clinical_domain_std=
    meta$clinical_domain_std,
  condition_binary=
    meta$condition_binary,
  integration_celltype_full_v1=
    backbone,
  transfer_confidence=
    meta$transfer_confidence,
  cluster_r0p6=r06,
  monodc_core_v1=core,
  monodc_core_confidence_v1=
    core_confidence,
  monodc_core_basis_v1=
    core_basis,
  monodc_state_v1=
    state,
  monodc_state_confidence_v1=
    state_confidence,
  monodc_state_basis_v1=
    state_basis,
  monodc_candidate_state_v1=
    candidate_state,
  monodc_candidate_confidence_v1=
    candidate_confidence,
  stringsAsFactors=FALSE
)

# ============================================================
# Validation
# ============================================================

stopifnot(
  nrow(out) == 192040,
  length(unique(out$global_cell)) ==
    192040,
  !anyNA(out$monodc_core_v1),
  !anyNA(out$monodc_state_v1),
  sum(out$cluster_r0p6 == "14") ==
    5348,
  sum(out$cluster_r0p6 == "21") ==
    3521,
  sum(out$cluster_r0p6 == "30") ==
    1241
)

# cDC2 refinement should affect only the dominant
# DC2/monocyte-boundary cells in cluster 14.
stopifnot(
  all(
    out$integration_celltype_full_v1[
      out$monodc_core_v1 ==
        "CD1C_FCER1A_cDC2_like"
    ] ==
      "DC2_monocyte_boundary"
  )
)

# ============================================================
# Summaries
# ============================================================

cat("\n===== CORE =====\n")

core_tab <- sort(
  table(out$monodc_core_v1),
  decreasing=TRUE
)

print(core_tab)

cat("\n===== STATE =====\n")

state_tab <- sort(
  table(out$monodc_state_v1),
  decreasing=TRUE
)

print(state_tab)

cat("\n===== CANDIDATE STATE =====\n")

candidate_tab <- sort(
  table(out$monodc_candidate_state_v1),
  decreasing=TRUE
)

print(candidate_tab)

cat("\n===== CORE x STATE =====\n")

print(
  addmargins(
    table(
      out$monodc_core_v1,
      out$monodc_state_v1
    )
  )
)

cat("\n===== STATE x CONDITION =====\n")

print(
  addmargins(
    table(
      out$monodc_state_v1,
      out$condition_binary
    )
  )
)

# ============================================================
# Write cell-level annotation
# ============================================================

cell_file <- file.path(
  out_dir,
  "monodc_cell_annotation_v1.tsv.gz"
)

con <- gzfile(
  cell_file,
  open="wt"
)

tryCatch(
  write.table(
    out,
    con,
    sep="\t",
    quote=1,
    qmethod="double",
    row.names=FALSE
  ),
  finally=close(con)
)

# ============================================================
# Summary tables
# ============================================================

core_summary <- as.data.frame(
  table(out$monodc_core_v1),
  stringsAsFactors=FALSE
)

names(core_summary) <- c(
  "monodc_core_v1",
  "n_cells"
)

core_summary$fraction <-
  core_summary$n_cells /
  nrow(out)

write.table(
  core_summary,
  file=file.path(
    out_dir,
    "monodc_core_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

state_summary <- as.data.frame(
  table(out$monodc_state_v1),
  stringsAsFactors=FALSE
)

names(state_summary) <- c(
  "monodc_state_v1",
  "n_cells"
)

state_summary$fraction <-
  state_summary$n_cells /
  nrow(out)

write.table(
  state_summary,
  file=file.path(
    out_dir,
    "monodc_state_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

candidate_summary <- as.data.frame(
  table(out$monodc_candidate_state_v1),
  stringsAsFactors=FALSE
)

names(candidate_summary) <- c(
  "monodc_candidate_state_v1",
  "n_cells"
)

candidate_summary$fraction <-
  candidate_summary$n_cells /
  nrow(out)

write.table(
  candidate_summary,
  file=file.path(
    out_dir,
    "monodc_candidate_state_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# Preserve audit snapshot
# ============================================================

audit <- read.delim(
  audit_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

selected_audit <- audit[
  as.character(audit$cluster) %in%
    c("3","14","21","30"),
  ,
  drop=FALSE
]

write.table(
  selected_audit,
  file=file.path(
    out_dir,
    "selected_cluster_audit_snapshot_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# RDS freeze
# ============================================================

saveRDS(
  list(
    annotation_version=
      "MonoDC_annotation_v1",
    n_cells=nrow(out),
    annotation=out,
    source_clustering_rds=
      cl_rds,
    source_metadata=
      meta_file,
    source_project_balanced_audit=
      audit_file,
    rules=list(
      clustering_resolution_scaffold=
        0.6,
      clustering_role=
        "state-discovery scaffold, not final taxonomy",
      core_backbone=
        allowed,
      core_refinement=
        "cluster14 DC2-boundary cells -> CD1C_FCER1A_cDC2_like",
      final_states=c(
        "cluster21=interferon_stimulated; medium confidence",
        "cluster30=LPAR1_RFX3_high; high confidence"
      ),
      candidate_only=
        "cluster3=RETN_PADI4_inflammatory_candidate; low confidence"
    )
  ),
  file.path(
    out_dir,
    "monodc_annotation_freeze_v1.rds"
  ),
  compress=FALSE
)

writeLines(
  c(
    "MONO/DC ANNOTATION FREEZE v1",
    "n_cells=192040",
    "Leiden r0.6 used as state-discovery scaffold only",
    "core backbone preserved from integration_celltype_full_v1",
    "cluster14 DC2-boundary cells refined to CD1C_FCER1A_cDC2_like",
    "cluster21 state=interferon_stimulated; confidence=medium",
    "cluster30 state=LPAR1_RFX3_high; confidence=high",
    "cluster3 RETN_PADI4 signature retained as candidate only",
    "project-specific Leiden clusters are not promoted to taxonomy"
  ),
  file.path(
    out_dir,
    "ANNOTATION_FREEZE_v1.txt"
  )
)

writeLines(
  "PASS",
  file.path(
    out_dir,
    "ANNOTATION_FREEZE_COMPLETE.ok"
  )
)

cat(
  "\nPASS: Mono/DC annotation freeze v1 completed\n"
)

