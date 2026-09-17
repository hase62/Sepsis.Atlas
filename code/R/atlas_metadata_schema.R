# Canonical metadata schema and harmonisation helpers for Sepsis.Atlas.
# Study-specific columns are retained after the canonical columns.

ATLAS_CANONICAL_META_COLUMNS <- c(
  # source identifiers
  "database_accession", "sra_study", "bioproject", "organism",
  "geo_accession", "biosample_accession", "experiment_accession", "run_accession",
  "cellranger_dir_id", "sample_id", "sample_alias", "paper_sample",
  "patient_id", "donor_id", "subject_id", "participant_id",
  "library_id", "library_unit_id", "library_batch_id", "pool_id", "pool_members",
  "pooled_flag", "demux_status", "demux_required_for_subject_assignment",
  # biological/clinical design
  "study_group", "factor", "is_sepsis", "cov19", "ards", "pneumonia",
  "septic_shock", "mortality", "icu_admission", "ventilation",
  "severity", "severityatday", "disease_phase", "outcome_group",
  "days", "timepoint_label", "sex", "age", "age_group", "sofa",
  "infection_site", "pathogen_etiology",
  "clinical_mapping_status", "clinical_match_confidence",
  # specimen and assay
  "source_name", "tissue_label", "source_blood_fraction", "cell_type",
  "blood_collection_tube", "pbmc_isolation_method", "ficoll", "rbc_lysis",
  "sample_storage", "frozen_or_fresh", "extracted_molecule",
  "platform", "sequencer_model", "library_type", "chemistry_version",
  "feature_barcoding", "vdj_capture", "assay_type", "library_layout",
  "aligner", "cellranger_or_equivalent_version", "reference_genome", "gene_annotation",
  "batch_unit_recommended", "inclusion_in_primary_atlas", "source_note", "read_mode",
  # atlas droplet-QC provenance
  "cell_calling_method", "cell_call_policy", "emptydrops_fdr_threshold",
  "emptydrops_lower", "validrops_label_dead", "validrops_drop_dead",
  "soupx_method", "soupx_status", "soupx_contamination_fraction",
  "n_raw_barcodes", "n_nonzero_barcodes", "n_emptydrops_pass",
  "n_validrops_pass", "n_called_after_droplet_qc",
  "preqc_stage", "metadata_embedded_in_loader"
)

ATLAS_DEPRECATED_META_ALIASES <- c(
  ARDS = "ards",
  SOFA = "sofa",
  PEEP = "peep",
  fresh_or_frozen = "frozen_or_fresh",
  shock_status = "septic_shock",
  cellranger_version = "cellranger_or_equivalent_version"
)

atlas_missing_like <- function(x) {
  is.na(x) | (is.character(x) & !nzchar(trimws(x)))
}

atlas_coalesce_column <- function(df, canonical, alias) {
  if (!alias %in% names(df)) return(df)
  if (!canonical %in% names(df)) {
    df[[canonical]] <- df[[alias]]
  } else {
    take <- atlas_missing_like(df[[canonical]]) & !atlas_missing_like(df[[alias]])
    if (any(take, na.rm = TRUE)) df[[canonical]][take] <- df[[alias]][take]
  }
  df
}

atlas_standardize_metadata <- function(df, project_id = NA_character_, drop_deprecated = TRUE) {
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)

  for (alias in names(ATLAS_DEPRECATED_META_ALIASES)) {
    df <- atlas_coalesce_column(df, ATLAS_DEPRECATED_META_ALIASES[[alias]], alias)
  }

  for (nm in setdiff(ATLAS_CANONICAL_META_COLUMNS, names(df))) df[[nm]] <- NA

  if (all(atlas_missing_like(df$database_accession)) && !is.na(project_id)) {
    df$database_accession <- project_id
  }
  if (all(atlas_missing_like(df$sample_id)) && "geo_accession" %in% names(df)) {
    df$sample_id <- df$geo_accession
  }
  if (all(atlas_missing_like(df$patient_id))) {
    for (candidate in c("participant_id", "subject_id", "donor_id")) {
      if (candidate %in% names(df) && any(!atlas_missing_like(df[[candidate]]))) {
        df$patient_id <- df[[candidate]]
        break
      }
    }
  }
  if (all(atlas_missing_like(df$metadata_embedded_in_loader))) {
    df$metadata_embedded_in_loader <- TRUE
  }

  if (drop_deprecated) {
    df <- df[, setdiff(names(df), names(ATLAS_DEPRECATED_META_ALIASES)), drop = FALSE]
  }

  canonical_present <- intersect(ATLAS_CANONICAL_META_COLUMNS, names(df))
  extra <- setdiff(names(df), canonical_present)
  df[, c(canonical_present, extra), drop = FALSE]
}

atlas_standardize_seurat_metadata <- function(obj, project_id = NA_character_) {
  md <- atlas_standardize_metadata(obj@meta.data, project_id = project_id, drop_deprecated = TRUE)
  rownames(md) <- colnames(obj)
  obj@meta.data <- md
  obj
}

atlas_metadata_dictionary <- function() {
  data.frame(
    column = ATLAS_CANONICAL_META_COLUMNS,
    category = c(
      rep("identifier", 23),
      rep("clinical_design", 23),
      rep("specimen_assay", 25),
      rep("droplet_qc", length(ATLAS_CANONICAL_META_COLUMNS) - 71)
    ),
    stringsAsFactors = FALSE
  )
}

# Cell-level droplet-QC fields are generated from the count matrix and must not
# be overwritten by sample-level metadata placeholders added by schema padding.
ATLAS_DROPLET_QC_META_COLUMNS <- c(
  "cell_calling_method", "cell_call_policy", "emptydrops_fdr_threshold",
  "emptydrops_lower", "validrops_label_dead", "validrops_drop_dead",
  "soupx_method", "soupx_status", "soupx_contamination_fraction",
  "n_raw_barcodes", "n_nonzero_barcodes", "n_emptydrops_pass",
  "n_validrops_pass", "n_called_after_droplet_qc", "preqc_stage",
  "metadata_embedded_in_loader"
)

atlas_attach_sample_metadata <- function(obj, meta_row) {
  stopifnot(inherits(obj, "Seurat"), nrow(meta_row) == 1L)

  for (nm in colnames(meta_row)) {
    value <- meta_row[[nm]][1]
    exists_in_object <- nm %in% colnames(obj@meta.data)

    # Matrix-derived droplet-QC values always take precedence over the
    # sample-table copy of the canonical schema (usually NA placeholders).
    if (nm %in% ATLAS_DROPLET_QC_META_COLUMNS && exists_in_object) {
      current <- obj@meta.data[[nm]]
      if (any(!atlas_missing_like(current))) next
    }

    # More generally, never replace a populated cell-level field with a
    # missing scalar from sample metadata.
    if (exists_in_object && all(atlas_missing_like(value))) {
      current <- obj@meta.data[[nm]]
      if (any(!atlas_missing_like(current))) next
    }

    obj[[nm]] <- value
  }

  obj
}

