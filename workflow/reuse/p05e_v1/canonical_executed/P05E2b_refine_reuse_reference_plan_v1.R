#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".", mustWork=TRUE)

suppressPackageStartupMessages(library(data.table))

tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

out <- file.path(
  root, "publication", "scientific_data",
  paste0("reuse_reference_plan_refined_v1__", tag)
)
dir.create(out, recursive=TRUE)

# ------------------------------------------------------------
# Latest P05E2
# ------------------------------------------------------------

dirs <- list.dirs(
  file.path(root, "publication", "scientific_data"),
  recursive=FALSE,
  full.names=TRUE
)

dirs <- sort(
  dirs[
    grepl(
      "reuse_reference_scmerge2_plan_v1__[0-9]{8}_[0-9]{6}$",
      basename(dirs)
    )
  ]
)

stopifnot(length(dirs) >= 1L)
p2 <- tail(dirs, 1L)

plan <- fread(
  file.path(p2, "scmerge2_final_run_plan.tsv")
)

ids <- fread(
  file.path(p2, "compartment_identity_support_summary.tsv")
)

dc <- fread(
  file.path(p2, "deconvolution_reference_policy.tsv")
)

# ------------------------------------------------------------
# Condition support
# ------------------------------------------------------------

plan[, healthy_fraction :=
  healthy_cells / (healthy_cells + disease_cells)
]

plan[, disease_fraction :=
  disease_cells / (healthy_cells + disease_cells)
]

plan[, min_condition_cells :=
  pmin(healthy_cells, disease_cells)
]

plan[, min_condition_fraction :=
  pmin(healthy_fraction, disease_fraction)
]

plan[, condition_support_ok :=
  !is.na(min_condition_cells) &
  min_condition_cells >= 100 &
  min_condition_fraction >= 0.01
]

# ------------------------------------------------------------
# Corrected-reference policy
# ------------------------------------------------------------

# T/NK and Mono/DC: existing validated parameters, refit on final cells.
plan[
  final_compartment_v1 == "T_NK",
  `:=`(
    final_action="REFIT_FINAL_CELL_UNIVERSE_FIXED_PARAMETERS",
    selected_ruvK=5L,
    ruvK_candidates="5"
  )
]

plan[
  final_compartment_v1 == "Monocyte_DC",
  `:=`(
    final_action="REFIT_FINAL_CELL_UNIVERSE_FIXED_PARAMETERS",
    selected_ruvK=2L,
    ruvK_candidates="2"
  )
]

# New reference models: select 3 vs 5 empirically.
for(x in c(
  "B_plasma",
  "Platelet_megakaryocyte",
  "Neutrophil"
)) {
  plan[
    final_compartment_v1 == x,
    `:=`(
      final_action=ifelse(
        x == "Neutrophil",
        "MODEL_SELECT_SCMERGE2_CAUTIOUS",
        "MODEL_SELECT_SCMERGE2"
      ),
      selected_ruvK=NA_integer_,
      ruvK_candidates="3,5",
      selected_k_pseudoBulk=5L,
      n_reference_genes=3000L,
      correction_batch="project_id",
      correction_condition="condition_binary",
      correction_celltype="atlas_core_identity_v1"
    )
  ]
}

# Erythroid: disease/healthy condition is effectively unidentifiable.
plan[
  final_compartment_v1 == "Erythroid",
  `:=`(
    final_action="NATIVE_ONLY_CONDITION_NOT_IDENTIFIABLE",
    selected_ruvK=NA_integer_,
    ruvK_candidates=NA_character_,
    selected_k_pseudoBulk=NA_integer_,
    n_reference_genes=NA_integer_,
    correction_batch=NA_character_,
    correction_condition=NA_character_,
    correction_celltype=NA_character_
  )
]

# Progenitor: source-study enrichment dominates.
plan[
  final_compartment_v1 == "Progenitor",
  `:=`(
    final_action="NATIVE_ONLY_STRUCTURALLY_CONFOUNDED",
    selected_ruvK=NA_integer_,
    ruvK_candidates=NA_character_,
    selected_k_pseudoBulk=NA_integer_,
    n_reference_genes=NA_integer_
  )
]

plan[
  final_compartment_v1 == "Deferred_unresolved",
  final_action := "EXCLUDE_FROM_CORRECTED_REFERENCE"
]

setorder(plan, -n_cells)

fwrite(
  plan,
  file.path(out, "scmerge2_refined_run_plan.tsv"),
  sep="\t"
)

# ------------------------------------------------------------
# Strengthen deconvolution CORE criteria
# ------------------------------------------------------------

dc2 <- merge(
  dc,
  ids[, .(
    final_compartment_v1,
    atlas_core_identity_v1,
    identity_max_project_fraction=max_project_fraction,
    identity_n_projects_ge200=n_projects_ge200
  )],
  by=c(
    "final_compartment_v1",
    "atlas_core_identity_v1"
  ),
  all.x=TRUE
)

dc2[, previous_reference_tier :=
  recommended_reference
]

dc2[, recommended_reference :=
  ifelse(
    previous_reference_tier == "CORE" &
    identity_max_project_fraction <= 0.80 &
    identity_n_projects_ge200 >= 3,
    "CORE",
    "EXTENDED_ONLY"
  )
]

dc2[, tier_reason :=
  fifelse(
    recommended_reference == "CORE",
    "resolved_identity_with_cross_project_support",
    fifelse(
      grepl(
        "unresolved|ambiguous|^Deferred_",
        atlas_core_identity_v1,
        ignore.case=TRUE
      ),
      "unresolved_or_deferred_identity",
      fifelse(
        identity_max_project_fraction > 0.80,
        "single_project_dominance_gt_0.80",
        fifelse(
          identity_n_projects_ge200 < 3,
          "fewer_than_3_projects_with_ge200_cells",
          "extended_reference_only"
        )
      )
    )
  )
]

setorder(
  dc2,
  recommended_reference,
  final_compartment_v1,
  -n_cells
)

fwrite(
  dc2,
  file.path(out, "deconvolution_reference_policy_refined.tsv"),
  sep="\t"
)

# ------------------------------------------------------------
# Explicit model-selection plan
# ------------------------------------------------------------

ms <- plan[
  grepl("^MODEL_SELECT", final_action),
  .(
    final_compartment_v1,
    n_cells,
    n_projects,
    max_project_fraction,
    healthy_cells,
    disease_cells,
    healthy_fraction,
    n_projects_with_both,
    ruvK_candidates,
    k_pseudoBulk=5L,
    genes=3000L,
    batch="project_id",
    condition="condition_binary",
    cellTypes="atlas_core_identity_v1"
  )
]

fwrite(
  ms,
  file.path(out, "scmerge2_model_selection_plan.tsv"),
  sep="\t"
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

summary <- c(
  "===== P05E2b REFINED REUSE REFERENCE PLAN =====",
  "status=PASS",
  paste0("n_cells=", sum(plan$n_cells)),
  "corrected_reference_fixed_parameter_compartments=2",
  "corrected_reference_model_selection_compartments=3",
  "corrected_reference_total_eligible_compartments=5",
  "native_only_erythroid=YES",
  "native_only_progenitor=YES",
  "deferred_unresolved_corrected_reference=EXCLUDED",
  paste0(
    "deconvolution_CORE_identities=",
    dc2[recommended_reference=="CORE", .N]
  ),
  paste0(
    "deconvolution_EXTENDED_ONLY_identities=",
    dc2[recommended_reference=="EXTENDED_ONLY", .N]
  ),
  "new_scMerge2_ruvK_candidates=3,5",
  "scMerge2_k_pseudoBulk=5",
  "corrected_reference_DE_use=NO",
  "native_processed_counts_DE_use=YES",
  "expression_modified=NO",
  "annotation_modified=NO",
  "cellranger_count=NOT_ACCESSED"
)

writeLines(
  summary,
  file.path(out, "P05E2b_SUMMARY.txt")
)

# ------------------------------------------------------------
# Root handoff copies
# ------------------------------------------------------------

handoff <- c(
  "P05E2b_SUMMARY.txt",
  "scmerge2_refined_run_plan.tsv",
  "scmerge2_model_selection_plan.tsv",
  "deconvolution_reference_policy_refined.tsv"
)

for(f in handoff) {
  file.copy(
    file.path(out, f),
    file.path(
      root,
      paste0(
        "P05E2b_",
        sub("^P05E2b_", "", f)
      )
    ),
    overwrite=TRUE
  )
}

cat(
  readLines(file.path(out, "P05E2b_SUMMARY.txt")),
  sep="\n"
)

cat("\nOUT_DIR=", out, "\n", sep="")
