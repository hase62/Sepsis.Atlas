#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".", mustWork=TRUE)

suppressPackageStartupMessages(library(data.table))

tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

out <- file.path(
  root, "publication", "scientific_data",
  paste0("reuse_reference_scmerge2_selection_v1__", tag)
)
dir.create(out, recursive=TRUE)

# ------------------------------------------------------------
# Latest P05E3
# ------------------------------------------------------------

dirs <- list.dirs(
  file.path(root, "publication", "scientific_data"),
  recursive=FALSE,
  full.names=TRUE
)

dirs <- sort(dirs[
  grepl(
    "reuse_reference_scmerge2_model_selection_v1__[0-9]{8}_[0-9]{6}$",
    basename(dirs)
  )
])

stopifnot(length(dirs) >= 1L)

p3 <- tail(dirs, 1L)

m <- fread(file.path(
  p3,
  "scmerge2_model_selection_metrics.tsv"
))

stopifnot(
  nrow(m) == 6L,
  setequal(unique(m$ruvK), c(3L,5L))
)

# ------------------------------------------------------------
# Freeze selected candidates
# ------------------------------------------------------------

selection <- data.table(
  compartment=c(
    "B_plasma",
    "Neutrophil",
    "Platelet_megakaryocyte"
  ),
  selected_ruvK=c(
    5L,
    3L,
    3L
  ),
  selection_status=c(
    "SELECTED",
    "SELECTED",
    "SELECTED_CAUTIOUS"
  ),
  rationale=c(
    paste(
      "Higher project-effect reduction and better",
      "cell-type preservation than ruvK=3;",
      "small distance-preservation trade-off."
    ),
    paste(
      "ruvK=3 and 5 nearly equivalent;",
      "ruvK=3 gives slightly greater project-effect",
      "reduction with comparable biological preservation."
    ),
    paste(
      "ruvK=3 preferred to avoid overcorrection;",
      "ruvK=5 improves project removal but more strongly",
      "reduces condition/cell-type and distance preservation."
    )
  )
)

selected <- merge(
  selection,
  m,
  by.x=c("compartment","selected_ruvK"),
  by.y=c("compartment","ruvK"),
  all.x=TRUE
)

stopifnot(nrow(selected)==3L)
stopifnot(!anyNA(selected$n_cells))

fwrite(
  selected,
  file.path(out, "selected_scmerge2_models.tsv"),
  sep="\t",
  quote=TRUE
)

# ------------------------------------------------------------
# Final five-compartment parameter freeze
# ------------------------------------------------------------

final_plan <- data.table(
  compartment=c(
    "T_NK",
    "Monocyte_DC",
    "B_plasma",
    "Neutrophil",
    "Platelet_megakaryocyte"
  ),
  ruvK=c(
    5L,
    2L,
    5L,
    3L,
    3L
  ),
  k_pseudoBulk=5L,
  parameter_source=c(
    "existing_primary_integration_anchor",
    "existing_primary_integration_anchor",
    "P05E3_model_selection",
    "P05E3_model_selection",
    "P05E3_model_selection"
  ),
  release_status=c(
    "FINAL",
    "FINAL",
    "FINAL",
    "FINAL",
    "FINAL_CAUTIOUS"
  ),
  correction_scope="within_compartment",
  intended_use="reference_mapping_annotation",
  DE_use="NO"
)

fwrite(
  final_plan,
  file.path(out, "final_corrected_reference_plan.tsv"),
  sep="\t"
)

# ------------------------------------------------------------
# Explicit non-corrected compartments
# ------------------------------------------------------------

native_only <- data.table(
  compartment=c(
    "Erythroid",
    "Progenitor",
    "Deferred_unresolved"
  ),
  corrected_reference_status=c(
    "NATIVE_ONLY",
    "NATIVE_ONLY",
    "EXCLUDED"
  ),
  reason=c(
    "condition not identifiable for reliable correction",
    "source-study enrichment structurally confounded",
    "unresolved annotation not suitable for corrected reference training"
  )
)

fwrite(
  native_only,
  file.path(out, "noncorrected_compartment_policy.tsv"),
  sep="\t"
)

summary <- c(
  "===== P05E3b SCMERGE2 SELECTION FREEZE =====",
  "status=PASS",
  "T_NK_ruvK=5",
  "Monocyte_DC_ruvK=2",
  "B_plasma_ruvK=5",
  "Neutrophil_ruvK=3",
  "Platelet_megakaryocyte_ruvK=3",
  "Platelet_megakaryocyte_status=CAUTIOUS",
  "Erythroid=NATIVE_ONLY",
  "Progenitor=NATIVE_ONLY",
  "Deferred_unresolved=EXCLUDED",
  "correction_scope=WITHIN_COMPARTMENT",
  "corrected_reference_DE_use=NO",
  "native_processed_counts_DE_use=YES",
  "expression_modified=NO",
  "annotation_modified=NO",
  "cellranger_count=NOT_ACCESSED"
)

writeLines(
  summary,
  file.path(out, "P05E3b_SUMMARY.txt")
)

# root handoff
for(f in c(
  "P05E3b_SUMMARY.txt",
  "selected_scmerge2_models.tsv",
  "final_corrected_reference_plan.tsv",
  "noncorrected_compartment_policy.tsv"
)) {
  file.copy(
    file.path(out,f),
    file.path(
      root,
      paste0(
        "P05E3b_",
        sub("^P05E3b_","",f)
      )
    ),
    overwrite=TRUE
  )
}

cat(readLines(file.path(out,"P05E3b_SUMMARY.txt")), sep="\n")
cat("\nOUT_DIR=", out, "\n", sep="")
