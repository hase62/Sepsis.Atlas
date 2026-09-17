#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

base <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1"
)

rds_file <- file.path(
  base,
  "integrated_compartment_clustering",
  "T_NK_combined",
  "integrated_clustering_v1.rds"
)

meta_file <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "scmerge2_primary_runs",
  "T_NK_combined__full_v1_cells.tsv.gz"
)

out_dir <- file.path(
  base,
  "tnk_annotation_freeze_v1"
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

x <- readRDS(rds_file)

cells <- as.character(x$cells)
r04 <- as.character(x$clusters$cluster_res_0p4)
r06 <- as.character(x$clusters$cluster_res_0p6)

meta <- read.delim(
  gzfile(meta_file),
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

m <- match(cells, meta$global_cell)

stopifnot(
  !anyNA(m),
  length(cells) == length(r04),
  length(cells) == length(r06)
)

meta <- meta[m,,drop=FALSE]

stopifnot(
  identical(
    as.character(meta$global_cell),
    cells
  )
)

n <- length(cells)

core <- rep(
  "Deferred_T_NK_ambiguous",
  n
)

state <- rep(
  "none",
  n
)

confidence <- rep(
  "medium",
  n
)

basis <- rep(
  "r0.4_native_marker",
  n
)

# ============================================================
# r0.4 core taxonomy
# ============================================================

set_core <- function(cluster, label, conf="high") {
  ii <- r04 == as.character(cluster)
  core[ii] <<- label
  confidence[ii] <<- conf
}

set_core(1,  "CD4_T_memory_like")
set_core(2,  "Cytotoxic_T_like")
# cluster 3 handled below with r0.6
set_core(4,  "CD8_T_GZMK_like")
set_core(5,  "CD4_T_naive_like")
set_core(6,  "Cytotoxic_NK_like")
set_core(7,  "NK_like")
set_core(8,  "CD8_T_naive_like")
set_core(9,  "MAIT_like")
# cluster 10 handled below
set_core(11, "Regulatory_T")
set_core(12, "Cytotoxic_NK_like")
set_core(13, "CD4_T_like", "medium")
set_core(14, "CD4_T_like", "medium")
set_core(15, "CD4_T_like", "medium")
set_core(16, "NK_like", "medium")
set_core(17, "Cycling_T_NK", "high")
set_core(18, "CD8_T_naive_like")
set_core(19, "CD4_T_like")
set_core(20, "CD4_T_like", "medium")
set_core(21, "Deferred_T_NK_ambiguous", "low")
set_core(22, "CD4_T_like")
set_core(23, "NK_like")
set_core(24, "CD4_T_memory_like", "medium")
set_core(25, "NK_like", "medium")
set_core(26, "Deferred_non_TNK", "high")
set_core(27, "Rare_lymphoid_ambiguous", "low")

# ============================================================
# Genuine biological r0.6 split:
# r0.4 cluster 3
# ============================================================

ii <- r04 == "3" & r06 == "7"
core[ii] <- "CD4_T_memory_like"
confidence[ii] <- "high"
basis[ii] <- "r0.4_plus_r0.6_subsplit"

ii <- r04 == "3" & r06 == "8"
core[ii] <- "CD4_T_naive_like"
confidence[ii] <- "high"
basis[ii] <- "r0.4_plus_r0.6_subsplit"

# Any tiny residual cells from parent 3 remain broad CD4.
ii <- r04 == "3" &
      !r06 %in% c("7","8")

core[ii] <- "CD4_T_like"
confidence[ii] <- "medium"
basis[ii] <- "r0.4_residual_after_subsplit"

# ============================================================
# r0.4 cluster 10:
# mixed / project-sensitive structure
# ============================================================

ii <- r04 == "10" & r06 == "13"

core[ii] <- "CD4_T_like"
state[ii] <- "project_specific_low_marker"
confidence[ii] <- "low"
basis[ii] <- "r0.6_project_specific_low_marker"

ii <- r04 == "10" & r06 == "26"

core[ii] <- "Deferred_T_NK_ambiguous"
state[ii] <- "cytotoxic_innate_skewed"
confidence[ii] <- "low"
basis[ii] <- "r0.6_mixed_NK_cytotoxic_T"

# tiny parent-10 residuals
ii <- r04 == "10" &
      !r06 %in% c("13","26")

core[ii] <- "Deferred_T_NK_ambiguous"
confidence[ii] <- "low"
basis[ii] <- "r0.4_unresolved"

# ============================================================
# State annotations
# ============================================================

# Early-response / activation-like CD4
state[r04 == "13"] <- "early_response"

# Disease/project-enriched CCR7/quiescent CD4 states
state[r04 == "14"] <- "project_specific_low_marker"
state[r04 == "20"] <- "project_specific_low_marker"

# CXCR4-high CD4 state
state[r04 == "15"] <- "CXCR4_high"

# Cycling
state[r04 == "17"] <- "cycling"

# IFN-stimulated T
state[r04 %in% c("19","22")] <-
  "interferon_stimulated"

# Activated / IFN-like NK
state[r04 == "23"] <-
  "interferon_stimulated"

# r0.6 child 32:
# clear IFN NK program, but strongly GSE216020/disease enriched.
ii <- r04 == "7" & r06 == "32"

state[ii] <- "IFN_NK_candidate"
basis[ii] <- "r0.6_state_project_confounded"
confidence[ii] <- "medium"

# ============================================================
# Final output
# ============================================================

out <- data.frame(
  global_cell=cells,
  project_id=meta$project_id,
  library_key=meta$library_key,
  condition_binary=meta$condition_binary,
  integration_celltype_full_v1=
    meta$integration_celltype_full_v1,
  transfer_confidence=
    meta$transfer_confidence,
  cluster_r0p4=r04,
  cluster_r0p6=r06,
  tnk_core_v1=core,
  tnk_state_v1=state,
  tnk_annotation_confidence_v1=
    confidence,
  tnk_annotation_basis_v1=
    basis,
  stringsAsFactors=FALSE
)

# ------------------------------------------------------------
# Validation
# ------------------------------------------------------------

stopifnot(
  nrow(out) == 281508,
  !anyNA(out$tnk_core_v1),
  !anyNA(out$tnk_state_v1),
  !anyNA(out$tnk_annotation_confidence_v1),
  length(unique(out$global_cell)) ==
    nrow(out)
)

cat("\n===== CORE =====\n")
print(
  sort(
    table(out$tnk_core_v1),
    decreasing=TRUE
  )
)

cat("\n===== STATE =====\n")
print(
  sort(
    table(out$tnk_state_v1),
    decreasing=TRUE
  )
)

cat("\n===== CONFIDENCE =====\n")
print(
  table(out$tnk_annotation_confidence_v1)
)

cat("\n===== CORE x STATE =====\n")
print(
  addmargins(
    table(
      out$tnk_core_v1,
      out$tnk_state_v1
    )
  )
)

out_file <- file.path(
  out_dir,
  "tnk_cell_annotation_v1.tsv.gz"
)

con <- gzfile(
  out_file,
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

summary_core <- as.data.frame(
  table(out$tnk_core_v1),
  stringsAsFactors=FALSE
)
names(summary_core) <- c(
  "tnk_core_v1",
  "n_cells"
)

summary_core$fraction <-
  summary_core$n_cells /
  nrow(out)

write.table(
  summary_core,
  file=file.path(
    out_dir,
    "tnk_core_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

summary_state <- as.data.frame(
  table(out$tnk_state_v1),
  stringsAsFactors=FALSE
)
names(summary_state) <- c(
  "tnk_state_v1",
  "n_cells"
)

summary_state$fraction <-
  summary_state$n_cells /
  nrow(out)

write.table(
  summary_state,
  file=file.path(
    out_dir,
    "tnk_state_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

saveRDS(
  list(
    annotation_version="T_NK_annotation_v1",
    n_cells=nrow(out),
    annotation=out,
    source_clustering_rds=rds_file,
    source_metadata=meta_file,
    rules=list(
      backbone_resolution=0.4,
      accepted_subsplit=
        "r0.4 cluster 3 -> r0.6 clusters 7/8",
      project_confounded_state=
        "r0.4 cluster 7 / r0.6 cluster 32",
      unresolved_parent=
        "r0.4 cluster 10",
      non_TNK_cluster=
        "r0.4 cluster 26"
    )
  ),
  file.path(
    out_dir,
    "tnk_annotation_freeze_v1.rds"
  ),
  compress=FALSE
)

writeLines(
  c(
    "T/NK ANNOTATION FREEZE v1",
    paste0(
      "n_cells=",
      nrow(out)
    ),
    "backbone_resolution=0.4",
    "accepted_subsplit=parent3:r0.6_7_vs_8",
    "cluster17=cycling_T_NK_not_lineage_subsplit",
    "cluster7_child32=IFN_NK_state_not_core_subtype",
    "cluster26=Deferred_non_TNK",
    "cluster10=low-confidence mixed structure"
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
  "\nPASS: T/NK annotation freeze v1 completed\n"
)

