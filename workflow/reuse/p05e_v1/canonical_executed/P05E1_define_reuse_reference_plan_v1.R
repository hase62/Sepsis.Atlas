#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".", mustWork=TRUE)

suppressPackageStartupMessages(library(data.table))

tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

out <- file.path(
  root, "publication", "scientific_data",
  paste0("reuse_reference_design_v1__", tag)
)
dir.create(out, recursive=TRUE)

# ------------------------------------------------------------
# Locate canonical P05E0
# ------------------------------------------------------------

p0dirs <- list.dirs(
  file.path(root, "publication", "scientific_data"),
  recursive=FALSE,
  full.names=TRUE
)

p0dirs <- sort(
  p0dirs[
    grepl(
      "reuse_reference_preflight_v1__[0-9]{8}_[0-9]{6}$",
      basename(p0dirs)
    )
  ]
)

stopifnot(length(p0dirs) >= 1L)
p0 <- tail(p0dirs, 1L)

comp_counts <- fread(
  file.path(p0, "compartment_by_project_cell_counts.tsv")
)

inv <- fread(
  file.path(p0, "final_qc_metadata_column_inventory.tsv")
)

existing_models <- fread(
  file.path(p0, "existing_scmerge2_model_audit.tsv")
)

# ------------------------------------------------------------
# P05B metadata
# ------------------------------------------------------------

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

cells <- fread(cell_file)
libs  <- fread(lib_file)

stopifnot(
  nrow(cells) == 665816L,
  nrow(libs) == 158L,
  !anyDuplicated(cells$global_cell)
)

# ------------------------------------------------------------
# Compartment support
# ------------------------------------------------------------

cp <- copy(comp_counts)

tot <- cp[, .(
  n_cells=sum(N),
  n_projects=.N,
  max_project_cells=max(N),
  max_project_fraction=max(N)/sum(N),
  n_projects_ge200=sum(N >= 200),
  n_projects_ge500=sum(N >= 500)
), by=final_compartment_v1]

clib <- cells[, .N, by=.(
  final_compartment_v1,
  project_id,
  library_key
)]

support <- clib[, .(
  n_libraries=.N,
  n_projects=uniqueN(project_id)
), by=final_compartment_v1]

support <- merge(
  tot, support,
  by="final_compartment_v1",
  all.x=TRUE,
  suffixes=c("", "_library")
)

# Subject support from P05B library metadata
subject_lib <- libs[, .(
  project_id,
  library_key,
  analysis_subject_id,
  clinical_domain_std
)]

tmp <- merge(
  unique(clib[, .(
    final_compartment_v1,
    project_id,
    library_key
  )]),
  subject_lib,
  by=c("project_id","library_key"),
  all.x=TRUE
)

subject_support <- tmp[, .(
  n_subjects=uniqueN(
    analysis_subject_id[
      !is.na(analysis_subject_id) &
      nzchar(analysis_subject_id)
    ]
  ),
  n_clinical_domains=uniqueN(
    clinical_domain_std[
      !is.na(clinical_domain_std) &
      nzchar(clinical_domain_std)
    ]
  )
), by=final_compartment_v1]

support <- merge(
  support,
  subject_support,
  by="final_compartment_v1",
  all.x=TRUE
)

# ------------------------------------------------------------
# Existing-model comparison
# ------------------------------------------------------------

existing_models[, final_compartment_v1 :=
  sub("_combined$", "", compartment)
]

existing_small <- existing_models[, .(
  final_compartment_v1,
  existing_model_n_cells=n_cells,
  existing_ruvK=ruvK,
  existing_k_pseudoBulk=k_pseudoBulk,
  existing_n_hvg=n_chosen_hvg,
  existing_n_controls=n_controls
)]

support <- merge(
  support,
  existing_small,
  by="final_compartment_v1",
  all.x=TRUE
)

support[, existing_model_cell_delta :=
  n_cells - existing_model_n_cells
]

# ------------------------------------------------------------
# Freeze reuse policy
# ------------------------------------------------------------

support[, reuse_reference_action := fifelse(
  final_compartment_v1 == "Deferred_unresolved",
  "EXCLUDE_FROM_CORRECTED_REFERENCE",
  fifelse(
    final_compartment_v1 == "Progenitor",
    "NATIVE_ONLY_STRUCTURALLY_CONFOUNDED",
    fifelse(
      final_compartment_v1 %in% c("T_NK","Monocyte_DC"),
      "REFIT_FINAL_CELL_UNIVERSE_FIXED_PARAMETERS",
      fifelse(
        max_project_fraction <= 0.60 &
        n_projects_ge200 >= 3,
        "FIT_NEW_SCMERGE2",
        fifelse(
          max_project_fraction <= 0.85 &
          n_projects_ge200 >= 3,
          "FIT_NEW_SCMERGE2_CAUTIOUS",
          "NATIVE_ONLY_INSUFFICIENT_CROSS_PROJECT_SUPPORT"
        )
      )
    )
  )
)]

support[, ruvK_candidates := fifelse(
  final_compartment_v1 == "T_NK", "5",
  fifelse(
    final_compartment_v1 == "Monocyte_DC", "2",
    fifelse(
      final_compartment_v1 == "B_plasma", "3,5",
      fifelse(
        final_compartment_v1 == "Platelet_megakaryocyte", "3,5",
        fifelse(
          final_compartment_v1 == "Neutrophil", "3,5",
          fifelse(
            final_compartment_v1 == "Erythroid", "1,2,3,5",
            NA_character_
          )
        )
      )
    )
  )
)]

support[, k_pseudoBulk := fifelse(
  grepl("^FIT|^REFIT", reuse_reference_action),
  5L,
  NA_integer_
)]

support[, batch_variable := fifelse(
  grepl("^FIT|^REFIT", reuse_reference_action),
  "project_id",
  NA_character_
)]

support[, celltype_variable := fifelse(
  grepl("^FIT|^REFIT", reuse_reference_action),
  "atlas_core_identity_v1",
  NA_character_
)]

setorder(support, -n_cells)

fwrite(
  support,
  file.path(out, "compartment_reuse_reference_policy.tsv"),
  sep="\t"
)

# ------------------------------------------------------------
# Existing scMerge2 condition mapping audit
# ------------------------------------------------------------

scdir <- file.path(
  root, "pre_integration",
  "full_atlas_primary_v1",
  "scmerge2_primary_runs"
)

cmaps <- list()

for(f in Sys.glob(file.path(scdir, "*__full_v1_cells.tsv.gz"))) {

  d <- fread(f)

  if(all(c(
    "project_id",
    "library_key",
    "clinical_domain_std",
    "condition_binary"
  ) %in% names(d))) {

    cmaps[[length(cmaps)+1L]] <- unique(
      d[, .(
        project_id,
        library_key,
        clinical_domain_std,
        condition_binary
      )]
    )
  }
}

if(length(cmaps)) {

  cmap <- unique(rbindlist(cmaps))

  fwrite(
    cmap,
    file.path(out, "existing_scmerge2_condition_mapping.tsv"),
    sep="\t"
  )

  cmap_summary <- cmap[, .N, by=.(
    clinical_domain_std,
    condition_binary
  )][order(clinical_domain_std, condition_binary)]

  fwrite(
    cmap_summary,
    file.path(out, "existing_scmerge2_condition_mapping_summary.tsv"),
    sep="\t"
  )

} else {

  fwrite(
    data.table(
      clinical_domain_std=character(),
      condition_binary=character(),
      N=integer()
    ),
    file.path(
      out,
      "existing_scmerge2_condition_mapping_summary.tsv"
    ),
    sep="\t"
  )
}

# ------------------------------------------------------------
# Rich metadata schema policy
# ------------------------------------------------------------

internal_pattern <- paste(
  c(
    "(^|_)path($|_)",
    "(^|_)dir($|_)",
    "(^|_)root($|_)",
    "final_rds",
    "input_root",
    "output_root",
    "fallback_filtered_dir",
    "command_hash"
  ),
  collapse="|"
)

schema <- copy(inv)

schema[, source := "final_qc_rds"]

schema[, metadata_location := fifelse(
  total_nonmissing_cells == 0,
  "EXCLUDE_EMPTY",
  fifelse(
    grepl(internal_pattern, column, ignore.case=TRUE),
    "EXCLUDE_INTERNAL_PATH",
    fifelse(
      max_unique_within_library > 1,
      "CELL_RICH",
      "LIBRARY_RICH"
    )
  )
)]

schema[, reason := fifelse(
  metadata_location == "EXCLUDE_EMPTY",
  "no non-missing values in final-QC RDS",
  fifelse(
    metadata_location == "EXCLUDE_INTERNAL_PATH",
    "local/internal filesystem provenance",
    fifelse(
      metadata_location == "CELL_RICH",
      "varies within at least one library",
      "library-constant metadata retained once per library"
    )
  )
)]

# P05B final annotation fields are mandatory cell metadata.
p05b_cell_schema <- data.table(
  column=names(cells),
  n_libraries_present=NA_integer_,
  n_libraries_with_nonmissing=NA_integer_,
  total_nonmissing_cells=NA_integer_,
  max_unique_within_library=NA_integer_,
  classes=vapply(cells, function(x)
    paste(class(x), collapse=";"), character(1)),
  source="P05B_final_annotation",
  metadata_location="CELL_RICH",
  reason="frozen final Atlas annotation metadata"
)

# P05B library metadata retained.
p05b_lib_schema <- data.table(
  column=names(libs),
  n_libraries_present=158L,
  n_libraries_with_nonmissing=NA_integer_,
  total_nonmissing_cells=NA_integer_,
  max_unique_within_library=1L,
  classes=vapply(libs, function(x)
    paste(class(x), collapse=";"), character(1)),
  source="P05B_library_metadata",
  metadata_location="LIBRARY_RICH",
  reason="canonical public library/source metadata"
)

schema_all <- rbindlist(
  list(
    p05b_cell_schema,
    p05b_lib_schema,
    schema
  ),
  fill=TRUE
)

# Selected library-level fields are also convenient in h5ad obs.
propagate <- c(
  "analysis_sample_id",
  "analysis_subject_id",
  "analysis_subject_id_source",
  "clinical_domain_std",
  "study_group",
  "study_group_std",
  "is_sepsis",
  "is_sepsis_std",
  "cov19",
  "cov19_std",
  "ards",
  "ards_std",
  "pneumonia",
  "septic_shock",
  "septic_shock_std",
  "mortality",
  "mortality_90d",
  "severity",
  "severity_label_std",
  "days",
  "sampling_day",
  "timepoint_label",
  "sex",
  "sex_std",
  "age",
  "age_years",
  "sofa",
  "sofa_score",
  "ventilation",
  "tissue_label",
  "tissue_std",
  "source_blood_fraction",
  "blood_fraction_std",
  "frozen_or_fresh",
  "preservation_std",
  "platform",
  "chemistry_version",
  "chemistry_std",
  "cellranger_or_equivalent_version",
  "reference_genome",
  "reference_genome_std"
)

schema_all[, propagate_to_cell_obs :=
  metadata_location == "CELL_RICH" |
  column %in% propagate
]

fwrite(
  schema_all,
  file.path(out, "rich_metadata_schema_candidates.tsv"),
  sep="\t",
  quote=TRUE,
  na="NA"
)

# ------------------------------------------------------------
# Deconvolution identity support
# ------------------------------------------------------------

dc <- cells[
  final_compartment_v1 != "Deferred_unresolved" &
  !is.na(atlas_core_identity_v1) &
  nzchar(atlas_core_identity_v1)
]

dc_lib <- dc[, .N, by=.(
  final_compartment_v1,
  atlas_core_identity_v1,
  project_id,
  library_key
)]

dc_lib <- merge(
  dc_lib,
  libs[, .(
    project_id,
    library_key,
    analysis_subject_id
  )],
  by=c("project_id","library_key"),
  all.x=TRUE
)

dc_support <- dc_lib[, .(
  n_cells=sum(N),
  n_libraries=.N,
  n_projects=uniqueN(project_id),
  n_subjects=uniqueN(
    analysis_subject_id[
      !is.na(analysis_subject_id) &
      nzchar(analysis_subject_id)
    ]
  ),
  max_library_fraction=max(N)/sum(N)
), by=.(
  final_compartment_v1,
  atlas_core_identity_v1
)]

dc_support[, support_tier := fifelse(
  n_projects >= 3 &
  n_libraries >= 5 &
  n_cells >= 200,
  "HIGH",
  fifelse(
    n_projects >= 2 &
    n_libraries >= 3 &
    n_cells >= 100,
    "MODERATE",
    "LIMITED"
  )
)]

setorder(
  dc_support,
  final_compartment_v1,
  -n_cells
)

fwrite(
  dc_support,
  file.path(out, "deconvolution_identity_support.tsv"),
  sep="\t"
)

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------

loc_counts <- schema_all[, .N, by=metadata_location]

summary_lines <- c(
  "===== P05E1 REUSE REFERENCE DESIGN =====",
  "status=PASS",
  paste0("n_cells=", nrow(cells)),
  paste0("n_libraries=", nrow(libs)),
  paste0("metadata_union_final_qc=", nrow(inv)),
  paste0(
    "metadata_cell_rich_candidates=",
    loc_counts[
      metadata_location=="CELL_RICH",
      sum(N)
    ]
  ),
  paste0(
    "metadata_library_rich_candidates=",
    loc_counts[
      metadata_location=="LIBRARY_RICH",
      sum(N)
    ]
  ),
  paste0(
    "deconvolution_identities=",
    nrow(dc_support)
  ),
  "batch_corrected_reference_DE_use=NO",
  "native_processed_counts_DE_use=YES",
  "progenitor_batch_correction=NATIVE_ONLY",
  "deferred_unresolved_corrected_reference=EXCLUDED",
  "cellranger_count=NOT_ACCESSED",
  "expression_modified=NO",
  "annotation_modified=NO"
)

writeLines(
  summary_lines,
  file.path(out, "P05E1_SUMMARY.txt")
)

# ------------------------------------------------------------
# Root handoff copies
# ------------------------------------------------------------

handoff <- c(
  "P05E1_SUMMARY.txt",
  "compartment_reuse_reference_policy.tsv",
  "rich_metadata_schema_candidates.tsv",
  "deconvolution_identity_support.tsv",
  "existing_scmerge2_condition_mapping_summary.tsv"
)

for(f in handoff) {

  src <- file.path(out, f)

  dst <- file.path(
    root,
    paste0(
      "P05E1_",
      sub("^P05E1_", "", f)
    )
  )

  file.copy(
    src,
    dst,
    overwrite=TRUE
  )
}

cat(readLines(file.path(out, "P05E1_SUMMARY.txt")), sep="\n")
cat("\nOUT_DIR=", out, "\n", sep="")
cat("Root handoff files copied: ", length(handoff), "\n", sep="")
