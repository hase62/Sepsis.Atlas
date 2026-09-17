#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".", mustWork=TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(SeuratObject)
})

tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

out <- file.path(
  root, "publication", "scientific_data",
  paste0("reuse_reference_scmerge2_plan_v1__", tag)
)
dir.create(out, recursive=TRUE)

# ============================================================
# Canonical inputs
# ============================================================

p05b <- file.path(
  root, "publication", "scientific_data",
  "processed_data_release_candidate_v1__20260819_194718"
)

cell_file <- file.path(
  p05b, "metadata",
  "SepsisAtlas_cell_metadata_v1__processed_data_release_candidate_v1__20260819_194718.tsv.gz"
)

lib_file <- file.path(
  p05b, "metadata",
  "SepsisAtlas_library_metadata_with_source_accessions_v1__processed_data_release_candidate_v1__20260819_194718.tsv"
)

stopifnot(file.exists(cell_file), file.exists(lib_file))

cells <- fread(cell_file)
libs  <- fread(lib_file)

stopifnot(
  nrow(cells) == 665816L,
  nrow(libs) == 158L,
  !anyDuplicated(cells$global_cell)
)

# ============================================================
# P05E1 policy
# ============================================================

p1dirs <- list.dirs(
  file.path(root, "publication", "scientific_data"),
  recursive=FALSE,
  full.names=TRUE
)

p1dirs <- sort(
  p1dirs[
    grepl(
      "reuse_reference_design_v1__[0-9]{8}_[0-9]{6}$",
      basename(p1dirs)
    )
  ]
)

stopifnot(length(p1dirs) >= 1L)
p1 <- tail(p1dirs, 1L)

policy <- fread(
  file.path(p1, "compartment_reuse_reference_policy.tsv")
)

dc_support <- fread(
  file.path(p1, "deconvolution_identity_support.tsv")
)

# ============================================================
# Condition mapping
#
# Frozen to match the existing scMerge2 primary integration.
# ============================================================

condition_map <- data.table(
  clinical_domain_std=c(
    "healthy_control",
    "sepsis",
    "covid19",
    "non_covid_respiratory",
    "mixed"
  ),
  condition_binary=c(
    "healthy",
    "disease",
    "disease",
    "disease",
    "disease"
  )
)

fwrite(
  condition_map,
  file.path(out, "condition_mapping_v1.tsv"),
  sep="\t"
)

# ============================================================
# Library metadata used for scMerge2
# ============================================================

req_lib <- c(
  "project_id",
  "library_key",
  "analysis_subject_id",
  "clinical_domain_std"
)

stopifnot(all(req_lib %in% names(libs)))

lib_meta <- libs[, ..req_lib]

lib_meta <- merge(
  lib_meta,
  condition_map,
  by="clinical_domain_std",
  all.x=TRUE
)

if(any(is.na(lib_meta$condition_binary))) {
  bad <- unique(
    lib_meta[
      is.na(condition_binary),
      clinical_domain_std
    ]
  )
  stop(
    "Unmapped clinical_domain_std: ",
    paste(bad, collapse=", ")
  )
}

# ============================================================
# Cell-level final reference universe
# ============================================================

ref_cells <- merge(
  cells,
  lib_meta,
  by=c("project_id","library_key"),
  all.x=TRUE,
  sort=FALSE
)

stopifnot(
  nrow(ref_cells) == 665816L,
  !any(is.na(ref_cells$condition_binary))
)

# ============================================================
# Exact original barcode reconstruction
#
# global_cell is:
#   <library_key>___<original_cell>
# ============================================================

prefix <- paste0(
  ref_cells$library_key,
  "___"
)

good_prefix <- startsWith(
  ref_cells$global_cell,
  prefix
)

if(!all(good_prefix)) {
  stop(
    "global_cell/library_key prefix mismatch: ",
    sum(!good_prefix)
  )
}

ref_cells[, original_cell :=
  substring(
    global_cell,
    nchar(library_key) + 4L
  )
]

stopifnot(
  all(nzchar(ref_cells$original_cell))
)

# ============================================================
# Project x condition support
# ============================================================

pc <- ref_cells[
  final_compartment_v1 != "Deferred_unresolved",
  .N,
  by=.(
    final_compartment_v1,
    project_id,
    condition_binary
  )
]

fwrite(
  pc,
  file.path(out, "compartment_project_condition_counts.tsv"),
  sep="\t"
)

pc_wide <- dcast(
  pc,
  final_compartment_v1 + project_id ~ condition_binary,
  value.var="N",
  fill=0
)

if(!"healthy" %in% names(pc_wide)) pc_wide[, healthy := 0L]
if(!"disease" %in% names(pc_wide)) pc_wide[, disease := 0L]

pc_wide[, has_both_conditions :=
  healthy > 0 & disease > 0
]

fwrite(
  pc_wide,
  file.path(out, "compartment_project_condition_balance.tsv"),
  sep="\t"
)

confound <- pc_wide[, .(
  n_projects=.N,
  n_projects_with_healthy=sum(healthy > 0),
  n_projects_with_disease=sum(disease > 0),
  n_projects_with_both=sum(has_both_conditions),
  healthy_cells=sum(healthy),
  disease_cells=sum(disease)
), by=final_compartment_v1]

# ============================================================
# Identity support per compartment/project
# ============================================================

id_support <- ref_cells[
  final_compartment_v1 != "Deferred_unresolved",
  .N,
  by=.(
    final_compartment_v1,
    atlas_core_identity_v1,
    project_id
  )
]

fwrite(
  id_support,
  file.path(out, "compartment_identity_project_support.tsv"),
  sep="\t"
)

id_summary <- id_support[, .(
  n_cells=sum(N),
  n_projects=uniqueN(project_id),
  max_project_fraction=max(N)/sum(N),
  n_projects_ge50=sum(N >= 50),
  n_projects_ge200=sum(N >= 200)
), by=.(
  final_compartment_v1,
  atlas_core_identity_v1
)]

fwrite(
  id_summary,
  file.path(out, "compartment_identity_support_summary.tsv"),
  sep="\t"
)

# ============================================================
# Locate 3000 consensus HVG reference
# ============================================================

hvg_candidates <- c(
  file.path(
    root,
    "pre_integration",
    "pilot_unintegrated_large_v1",
    "large_pilot_consensus_hvg_3000.tsv"
  ),
  file.path(
    root,
    "archive",
    "internal_provenance",
    "cold_archive_v1__20260820_233312",
    "preintegration_derived_detail_v1.tar.gz"
  )
)

hvg_direct <- hvg_candidates[
  file.exists(hvg_candidates) &
  grepl("\\.tsv$", hvg_candidates)
]

hvg_status <- if(length(hvg_direct)) {
  "DIRECT_FILE_AVAILABLE"
} else {
  "DIRECT_FILE_NOT_AVAILABLE"
}

# ============================================================
# Audit final-QC feature space
# ============================================================

rds <- Sys.glob(
  file.path(
    root,
    "pre_integration",
    "GSE*",
    "rds_final_qc",
    "*.rds"
  )
)

stopifnot(length(rds) == 158L)

first <- readRDS(rds[[1]])
feature_ref <- rownames(first[["RNA"]])
rm(first)
invisible(gc())

feature_audit <- data.table(
  n_final_qc_rds=length(rds),
  n_reference_features=length(feature_ref),
  hvg_source_status=hvg_status,
  hvg_direct_path=if(length(hvg_direct))
    hvg_direct[[1]] else NA_character_
)

fwrite(
  feature_audit,
  file.path(out, "feature_space_readiness.tsv"),
  sep="\t"
)

# ============================================================
# Exact final scMerge2 plan
# ============================================================

run_plan <- merge(
  policy,
  confound,
  by="final_compartment_v1",
  all.x=TRUE,
  suffixes=c("", "_condition")
)

run_plan[, final_action :=
  reuse_reference_action
]

# Strong guardrails.
run_plan[
  final_compartment_v1 == "Progenitor",
  final_action := "NATIVE_ONLY_STRUCTURALLY_CONFOUNDED"
]

run_plan[
  final_compartment_v1 == "Deferred_unresolved",
  final_action := "EXCLUDE_FROM_CORRECTED_REFERENCE"
]

# Exact parameter freeze.
run_plan[, selected_ruvK := fifelse(
  final_compartment_v1 == "T_NK", 5L,
  fifelse(
    final_compartment_v1 == "Monocyte_DC", 2L,
    fifelse(
      final_compartment_v1 == "B_plasma", 3L,
      fifelse(
        final_compartment_v1 == "Platelet_megakaryocyte", 3L,
        fifelse(
          final_compartment_v1 == "Neutrophil", 3L,
          fifelse(
            final_compartment_v1 == "Erythroid", 2L,
            NA_integer_
          )
        )
      )
    )
  )
)]

run_plan[, selected_k_pseudoBulk :=
  fifelse(
    grepl("^FIT|^REFIT", final_action),
    5L,
    NA_integer_
  )
]

run_plan[, n_reference_genes :=
  fifelse(
    grepl("^FIT|^REFIT", final_action),
    3000L,
    NA_integer_
  )
]

run_plan[, correction_batch :=
  fifelse(
    grepl("^FIT|^REFIT", final_action),
    "project_id",
    NA_character_
  )
]

run_plan[, correction_condition :=
  fifelse(
    grepl("^FIT|^REFIT", final_action),
    "condition_binary",
    NA_character_
  )
]

run_plan[, correction_celltype :=
  fifelse(
    grepl("^FIT|^REFIT", final_action),
    "atlas_core_identity_v1",
    NA_character_
  )
]

setorder(run_plan, -n_cells)

fwrite(
  run_plan,
  file.path(out, "scmerge2_final_run_plan.tsv"),
  sep="\t"
)

# ============================================================
# Deconvolution CORE vs EXTENDED policy
# ============================================================

dc <- copy(dc_support)

unresolved_pattern <- paste(
  c(
    "unresolved",
    "ambiguous",
    "^Deferred_",
    "transfer_unresolved"
  ),
  collapse="|"
)

dc[, recommended_reference :=
  fifelse(
    support_tier == "HIGH" &
    !grepl(
      unresolved_pattern,
      atlas_core_identity_v1,
      ignore.case=TRUE
    ),
    "CORE",
    "EXTENDED_ONLY"
  )
]

# Progenitor_like remains valid for native deconvolution reference.
dc[
  atlas_core_identity_v1 == "Progenitor_like" &
  support_tier == "HIGH",
  recommended_reference := "CORE"
]

fwrite(
  dc,
  file.path(out, "deconvolution_reference_policy.tsv"),
  sep="\t"
)

# ============================================================
# Summary
# ============================================================

n_fit <- run_plan[
  grepl("^FIT|^REFIT", final_action),
  .N
]

n_native <- run_plan[
  grepl("^NATIVE_ONLY", final_action),
  .N
]

n_excluded <- run_plan[
  grepl("^EXCLUDE", final_action),
  .N
]

summary_lines <- c(
  "===== P05E2 SCMERGE2 RUN PLAN FREEZE =====",
  "status=PASS",
  paste0("n_cells=", nrow(ref_cells)),
  paste0("n_libraries=", nrow(libs)),
  paste0("n_final_qc_rds=", length(rds)),
  paste0("feature_space=", length(feature_ref)),
  paste0("hvg_source_status=", hvg_status),
  paste0("compartments_to_scMerge2=", n_fit),
  paste0("compartments_native_only=", n_native),
  paste0("compartments_excluded=", n_excluded),
  paste0(
    "deconvolution_CORE_identities=",
    dc[recommended_reference=="CORE", .N]
  ),
  paste0(
    "deconvolution_EXTENDED_ONLY_identities=",
    dc[recommended_reference=="EXTENDED_ONLY", .N]
  ),
  "scMerge2_batch=project_id",
  "scMerge2_condition=condition_binary",
  "scMerge2_cellTypes=atlas_core_identity_v1",
  "scMerge2_k_pseudoBulk=5",
  "corrected_reference_DE_use=NO",
  "native_processed_counts_DE_use=YES",
  "cellranger_count=NOT_ACCESSED",
  "expression_modified=NO",
  "annotation_modified=NO"
)

writeLines(
  summary_lines,
  file.path(out, "P05E2_SUMMARY.txt")
)

# ============================================================
# Root handoff copies
# ============================================================

handoff <- c(
  "P05E2_SUMMARY.txt",
  "scmerge2_final_run_plan.tsv",
  "compartment_project_condition_balance.tsv",
  "compartment_identity_support_summary.tsv",
  "deconvolution_reference_policy.tsv",
  "feature_space_readiness.tsv"
)

for(f in handoff) {
  src <- file.path(out, f)
  dst <- file.path(
    root,
    paste0(
      "P05E2_",
      sub("^P05E2_", "", f)
    )
  )
  file.copy(src, dst, overwrite=TRUE)
}

cat(
  readLines(
    file.path(out, "P05E2_SUMMARY.txt")
  ),
  sep="\n"
)

cat(
  "\nOUT_DIR=",
  out,
  "\n",
  sep=""
)

cat(
  "Root handoff files copied: ",
  length(handoff),
  "\n",
  sep=""
)
