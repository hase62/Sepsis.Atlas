#!/usr/bin/env Rscript

# P05E6 recovery/finalization v1
#
# Recovers the latest failed P05E6 run from the already-written native
# sample pseudobulk MatrixMarket file. It DOES NOT reread the 158 final-QC
# Seurat RDS objects and DOES NOT recompute pseudobulks.
#
# Policy refinement:
#   - CORE identities remain mandatory and must have >=3 eligible projects.
#   - EXTENDED_ONLY identities that are unresolved/deferred/ambiguous are
#     not published as deconvolution targets.
#   - Other resolved EXTENDED_ONLY identities are published if at least
#     one eligible project profile remains after sample/project QC.
#
# Native SoupX-corrected counts remain the expression source.
# scMerge2-adjusted expression is NOT used for deconvolution.
# cellranger_count/ is not accessed.

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".", mustWork=TRUE)

Sys.setenv(
  OMP_NUM_THREADS="1",
  OPENBLAS_NUM_THREADS="1",
  MKL_NUM_THREADS="1",
  NUMEXPR_NUM_THREADS="1"
)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(matrixStats)
})

sd_base <- file.path(
  root,
  "publication",
  "scientific_data"
)

dirs <- sort(
  list.dirs(
    sd_base,
    recursive=FALSE,
    full.names=TRUE
  )
)

dirs <- dirs[
  grepl(
    "reuse_reference_deconvolution_v1__[0-9]{8}_[0-9]{6}$",
    basename(dirs)
  )
]

if(!length(dirs)) {
  stop("No P05E6 deconvolution directory found")
}

out <- tail(dirs, 1L)

cat("RECOVER_DIR=", out, "\n", sep="")

# -------------------------------------------------------------------------
# Required outputs from the completed pseudobulk phase
# -------------------------------------------------------------------------

pb_mtx_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_sample_pseudobulk_counts_v1.mtx.gz"
)

pb_meta_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_sample_pseudobulk_metadata_v1.tsv.gz"
)

feature_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_features_v1.tsv.gz"
)

identity_support_file <- file.path(
  out,
  "deconvolution_identity_support.tsv"
)

project_support_file <- file.path(
  out,
  "deconvolution_project_identity_support.tsv"
)

blood_support_file <- file.path(
  out,
  "deconvolution_identity_blood_fraction_support.tsv"
)

required_files <- c(
  pb_mtx_file,
  pb_meta_file,
  feature_file,
  identity_support_file,
  project_support_file,
  blood_support_file
)

missing <- required_files[
  !file.exists(required_files)
]

if(length(missing)) {
  stop(
    "Recovery inputs missing: ",
    paste(basename(missing), collapse=", ")
  )
}

for(f in c(
  pb_mtx_file,
  pb_meta_file,
  feature_file
)) {
  st <- system2(
    "gzip",
    args=c("-t", f),
    stdout=FALSE,
    stderr=FALSE
  )

  if(st != 0L) {
    stop("gzip integrity failed: ", f)
  }
}

# -------------------------------------------------------------------------
# Read saved pseudobulk objects
# -------------------------------------------------------------------------

cat("Reading saved pseudobulk MatrixMarket...\n")

con <- gzfile(
  pb_mtx_file,
  open="rt"
)

pb_counts <- Matrix::readMM(con)
close(con)

pb_counts <- methods::as(
  pb_counts,
  "dgCMatrix"
)

profile_meta <- fread(
  pb_meta_file
)

feature_map <- fread(
  feature_file
)

identity_support <- fread(
  identity_support_file
)

project_support <- fread(
  project_support_file
)

blood_support <- fread(
  blood_support_file
)

required_profile_cols <- c(
  "project_id",
  "library_key",
  "analysis_sample_id",
  "analysis_subject_id",
  "condition_binary",
  "reference_identity",
  "recommended_reference",
  "n_cells",
  "profile_id"
)

stopifnot(
  all(required_profile_cols %in% names(profile_meta)),
  nrow(feature_map) == 38606L,
  nrow(pb_counts) == 38606L,
  ncol(pb_counts) == nrow(profile_meta),
  !anyDuplicated(profile_meta$profile_id)
)

if(!"gene_symbol" %in% names(feature_map)) {
  stop("gene_symbol missing from feature metadata")
}

rownames(pb_counts) <- as.character(
  feature_map$gene_symbol
)

colnames(pb_counts) <- profile_meta$profile_id

# -------------------------------------------------------------------------
# Refine publishable identity policy
# -------------------------------------------------------------------------

required_identity_cols <- c(
  "reference_identity",
  "final_compartment_v1",
  "atlas_core_identity_v1",
  "recommended_reference",
  "n_eligible_projects"
)

stopifnot(
  all(required_identity_cols %in% names(identity_support))
)

identity_policy <- copy(
  identity_support
)

identity_policy[, unresolved_label :=
  grepl(
    "deferred|unresolved|ambiguous",
    atlas_core_identity_v1,
    ignore.case=TRUE
  )
]

identity_policy[, publish_status :=
  fifelse(
    recommended_reference == "CORE" &
    n_eligible_projects >= 3L,
    "CORE",
    fifelse(
      recommended_reference == "CORE",
      "EXCLUDE_CORE_INSUFFICIENT_SUPPORT",
      fifelse(
        unresolved_label == TRUE,
        "EXCLUDE_UNRESOLVED_EXTENDED",
        fifelse(
          recommended_reference == "EXTENDED_ONLY" &
          n_eligible_projects >= 1L,
          "EXTENDED_ONLY",
          "EXCLUDE_EXTENDED_NO_ELIGIBLE_PROJECT"
        )
      )
    )
  )
]

bad_core <- identity_policy[
  publish_status == "EXCLUDE_CORE_INSUFFICIENT_SUPPORT"
]

if(nrow(bad_core)) {
  fwrite(
    bad_core,
    file.path(
      out,
      "P05E6_recovery_CORE_support_failure.tsv"
    ),
    sep="\t"
  )

  stop(
    "CORE identity support failure: ",
    paste(
      bad_core$reference_identity,
      collapse=", "
    )
  )
}

core_ids <- identity_policy[
  publish_status == "CORE",
  reference_identity
]

extended_only_ids <- identity_policy[
  publish_status == "EXTENDED_ONLY",
  reference_identity
]

published_extended_ids <- c(
  core_ids,
  extended_only_ids
)

if(length(core_ids) != 20L) {
  stop(
    "Expected 20 CORE identities; observed ",
    length(core_ids)
  )
}

if(!length(extended_only_ids)) {
  warning(
    "No resolved EXTENDED_ONLY identities survived support filtering; ",
    "EXTENDED reference will equal CORE."
  )
}

excluded_policy <- identity_policy[
  grepl("^EXCLUDE_", publish_status)
]

fwrite(
  identity_policy,
  file.path(
    out,
    "SepsisAtlas_deconvolution_identity_publication_policy_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  na="NA"
)

fwrite(
  excluded_policy,
  file.path(
    out,
    "SepsisAtlas_deconvolution_excluded_identities_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  na="NA"
)

cat(
  "CORE=", length(core_ids),
  "; EXTENDED_ONLY_PUBLISHABLE=", length(extended_only_ids),
  "; EXTENDED_TOTAL=", length(published_extended_ids),
  "; EXCLUDED=", nrow(excluded_policy),
  "\n",
  sep=""
)

# -------------------------------------------------------------------------
# CPM normalization of sample pseudobulks
# -------------------------------------------------------------------------

library_size <- Matrix::colSums(
  pb_counts
)

if(
  any(!is.finite(library_size)) ||
  any(library_size <= 0)
) {
  stop("Invalid pseudobulk library sizes")
}

pb_cpm <- pb_counts %*%
  Matrix::Diagonal(
    x=1e6 / library_size
  )

# -------------------------------------------------------------------------
# Repeated sample balancing:
# sample -> subject x condition x identity
# -------------------------------------------------------------------------

profile_meta[, subject_condition_id :=
  paste(
    project_id,
    analysis_subject_id,
    condition_binary,
    reference_identity,
    sep="::"
  )
]

subject_levels <- unique(
  profile_meta$subject_condition_id
)

subject_index <- match(
  profile_meta$subject_condition_id,
  subject_levels
)

n_per_subject <- tabulate(
  subject_index,
  nbins=length(subject_levels)
)

A_subject <- Matrix::sparseMatrix(
  i=seq_len(nrow(profile_meta)),
  j=subject_index,
  x=1 / n_per_subject[subject_index],
  dims=c(
    nrow(profile_meta),
    length(subject_levels)
  )
)

subject_cpm <- pb_cpm %*% A_subject

rownames(subject_cpm) <- rownames(pb_counts)
colnames(subject_cpm) <- subject_levels

subject_meta <- profile_meta[
  ,
  .(
    project_id=project_id[[1]],
    analysis_subject_id=analysis_subject_id[[1]],
    condition_binary=condition_binary[[1]],
    reference_identity=reference_identity[[1]],
    recommended_reference=recommended_reference[[1]],
    n_samples=.N,
    n_cells=sum(n_cells)
  ),
  by=subject_condition_id
]

subject_meta <- subject_meta[
  match(
    subject_levels,
    subject_condition_id
  )
]

stopifnot(
  identical(
    subject_meta$subject_condition_id,
    subject_levels
  )
)

# -------------------------------------------------------------------------
# Subject -> project x condition x identity median
# -------------------------------------------------------------------------

pc_groups <- subject_meta[
  ,
  .(
    subject_indices=list(.I),
    n_subjects=.N,
    n_samples=sum(n_samples),
    n_cells=sum(n_cells),
    recommended_reference=recommended_reference[[1]]
  ),
  by=.(
    project_id,
    condition_binary,
    reference_identity
  )
]

pc_profile_names <- paste(
  pc_groups$project_id,
  pc_groups$condition_binary,
  pc_groups$reference_identity,
  sep="::"
)

pc_profiles <- matrix(
  0,
  nrow=nrow(pb_counts),
  ncol=nrow(pc_groups),
  dimnames=list(
    rownames(pb_counts),
    pc_profile_names
  )
)

for(i in seq_len(nrow(pc_groups))) {

  idx <- pc_groups$subject_indices[[i]]

  if(length(idx) == 1L) {
    pc_profiles[, i] <- as.numeric(
      subject_cpm[, idx]
    )
  } else {
    pc_profiles[, i] <- matrixStats::rowMedians(
      as.matrix(
        subject_cpm[
          ,
          idx,
          drop=FALSE
        ]
      )
    )
  }
}

pc_groups[, project_condition_profile_id :=
  pc_profile_names
]

# -------------------------------------------------------------------------
# Equal healthy/disease weighting within project where both exist
# -------------------------------------------------------------------------

project_groups <- pc_groups[
  ,
  .(
    pc_indices=list(.I),
    n_conditions=.N,
    n_subjects=sum(n_subjects),
    n_samples=sum(n_samples),
    n_cells=sum(n_cells),
    conditions=paste(
      sort(
        unique(condition_binary)
      ),
      collapse=","
    ),
    recommended_reference=recommended_reference[[1]]
  ),
  by=.(
    project_id,
    reference_identity
  )
]

project_profile_names <- paste(
  project_groups$project_id,
  project_groups$reference_identity,
  sep="::"
)

project_profiles <- matrix(
  0,
  nrow=nrow(pb_counts),
  ncol=nrow(project_groups),
  dimnames=list(
    rownames(pb_counts),
    project_profile_names
  )
)

for(i in seq_len(nrow(project_groups))) {

  idx <- project_groups$pc_indices[[i]]

  project_profiles[, i] <- rowMeans(
    pc_profiles[
      ,
      idx,
      drop=FALSE
    ]
  )
}

project_groups[, project_profile_id :=
  project_profile_names
]

# Attach frozen project-profile eligibility.
eligible_lookup <- project_support[
  ,
  .(
    reference_identity,
    project_id,
    project_profile_eligible
  )
]

project_groups <- merge(
  project_groups,
  eligible_lookup,
  by=c(
    "reference_identity",
    "project_id"
  ),
  all.x=TRUE,
  sort=FALSE
)

project_groups <- project_groups[
  match(
    project_profile_names,
    project_profile_id
  )
]

if(any(is.na(
  project_groups$project_profile_eligible
))) {
  stop("Missing project eligibility after recovery merge")
}

# -------------------------------------------------------------------------
# Final project-balanced reference
# -------------------------------------------------------------------------

publication_order <- identity_policy[
  publish_status %in% c(
    "CORE",
    "EXTENDED_ONLY"
  )
]

publication_order[, tier_order :=
  fifelse(
    publish_status == "CORE",
    1L,
    2L
  )
]

setorder(
  publication_order,
  tier_order,
  final_compartment_v1,
  atlas_core_identity_v1
)

extended_ids <- publication_order$reference_identity

core_order <- publication_order[
  publish_status == "CORE",
  reference_identity
]

stopifnot(
  setequal(core_order, core_ids),
  length(core_order) == 20L
)

final_profiles <- matrix(
  NA_real_,
  nrow=nrow(pb_counts),
  ncol=length(extended_ids),
  dimnames=list(
    rownames(pb_counts),
    extended_ids
  )
)

final_support <- vector(
  "list",
  length(extended_ids)
)

normalize_to_cpm <- function(x) {
  s <- sum(x)

  if(!is.finite(s) || s <= 0) {
    stop("Non-positive final profile sum")
  }

  x / s * 1e6
}

for(i in seq_along(extended_ids)) {

  ref_id <- extended_ids[[i]]

  idx <- which(
    project_groups$reference_identity == ref_id &
    project_groups$project_profile_eligible == TRUE
  )

  if(!length(idx)) {
    stop(
      "Internal error: publishable identity has no eligible project: ",
      ref_id
    )
  }

  x <- if(length(idx) == 1L) {
    project_profiles[, idx]
  } else {
    matrixStats::rowMedians(
      project_profiles[
        ,
        idx,
        drop=FALSE
      ]
    )
  }

  x <- normalize_to_cpm(
    as.numeric(x)
  )

  final_profiles[, i] <- x

  pol <- identity_policy[
    reference_identity == ref_id
  ]

  final_support[[i]] <- data.table(
    reference_identity=ref_id,
    final_compartment_v1=
      pol$final_compartment_v1[[1]],
    atlas_core_identity_v1=
      pol$atlas_core_identity_v1[[1]],
    published_tier=
      pol$publish_status[[1]],
    original_recommended_reference=
      pol$recommended_reference[[1]],
    n_eligible_projects=length(idx),
    eligible_projects=paste(
      project_groups$project_id[idx],
      collapse=","
    ),
    n_subject_project_pairs=sum(
      project_groups$n_subjects[idx]
    ),
    n_samples=sum(
      project_groups$n_samples[idx]
    ),
    n_cells=sum(
      project_groups$n_cells[idx]
    ),
    conditions_available=paste(
      sort(
        unique(
          project_groups$conditions[idx]
        )
      ),
      collapse=";"
    )
  )
}

final_support <- rbindlist(
  final_support
)

if(
  any(!is.finite(final_profiles)) ||
  any(final_profiles < 0)
) {
  stop("Invalid values in final deconvolution reference")
}

profile_sums <- colSums(
  final_profiles
)

if(max(abs(profile_sums - 1e6)) > 1) {
  stop("Final reference columns do not sum to ~1e6 CPM")
}

core_cpm <- final_profiles[
  ,
  core_order,
  drop=FALSE
]

extended_cpm <- final_profiles[
  ,
  extended_ids,
  drop=FALSE
]

core_log2cpm <- log2(
  core_cpm + 1
)

extended_log2cpm <- log2(
  extended_cpm + 1
)

# -------------------------------------------------------------------------
# Build public tables
# -------------------------------------------------------------------------

feature_cols <- c(
  "feature_index",
  "ensembl_gene_id",
  "source_ensembl_id",
  "gene_symbol",
  "feature_type",
  "mapping_source",
  "mapping_rule"
)

if(!all(feature_cols %in% names(feature_map))) {
  stop("Public feature columns missing")
}

feature_prefix <- feature_map[
  ,
  ..feature_cols
]

make_reference_table <- function(mat) {

  x <- as.data.table(
    mat
  )

  setnames(
    x,
    colnames(mat)
  )

  cbind(
    feature_prefix,
    x
  )
}

write_gzip_tsv <- function(dt, gzfile_path) {

  plain <- sub(
    "\\.gz$",
    "",
    gzfile_path
  )

  fwrite(
    dt,
    plain,
    sep="\t",
    quote=TRUE,
    na="NA"
  )

  st <- system2(
    "gzip",
    args=c("-f", plain),
    stdout=FALSE,
    stderr=FALSE
  )

  if(st != 0L) {
    stop("gzip failed: ", plain)
  }

  st <- system2(
    "gzip",
    args=c("-t", gzfile_path),
    stdout=FALSE,
    stderr=FALSE
  )

  if(st != 0L) {
    stop("gzip -t failed: ", gzfile_path)
  }
}

core_cpm_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_reference_CORE_CPM_v1.tsv.gz"
)

core_log_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_reference_CORE_log2CPM_v1.tsv.gz"
)

ext_cpm_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_reference_EXTENDED_CPM_v1.tsv.gz"
)

ext_log_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_reference_EXTENDED_log2CPM_v1.tsv.gz"
)

write_gzip_tsv(
  make_reference_table(core_cpm),
  core_cpm_file
)

write_gzip_tsv(
  make_reference_table(core_log2cpm),
  core_log_file
)

write_gzip_tsv(
  make_reference_table(extended_cpm),
  ext_cpm_file
)

write_gzip_tsv(
  make_reference_table(extended_log2cpm),
  ext_log_file
)

# Project-level profiles and metadata are useful for sensitivity analysis.
project_profile_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_project_condition_balanced_CPM_v1.tsv.gz"
)

write_gzip_tsv(
  make_reference_table(project_profiles),
  project_profile_file
)

project_meta_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_project_profile_metadata_v1.tsv.gz"
)

write_gzip_tsv(
  project_groups,
  project_meta_file
)

support_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_reference_support_v1.tsv"
)

fwrite(
  final_support,
  support_file,
  sep="\t",
  quote=TRUE,
  na="NA"
)

identity_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_reference_identities_v1.tsv"
)

identity_table <- merge(
  identity_policy,
  final_support[
    ,
    .(
      reference_identity,
      published_tier,
      n_eligible_projects,
      n_subject_project_pairs,
      n_samples,
      n_cells
    )
  ],
  by="reference_identity",
  all.x=TRUE,
  sort=FALSE
)

identity_table[, CORE_column_order :=
  match(
    reference_identity,
    core_order
  )
]

identity_table[, EXTENDED_column_order :=
  match(
    reference_identity,
    extended_ids
  )
]

fwrite(
  identity_table,
  identity_file,
  sep="\t",
  quote=TRUE,
  na="NA"
)

# -------------------------------------------------------------------------
# Re-read verification
# -------------------------------------------------------------------------

check_core <- fread(
  core_cpm_file
)

check_ext <- fread(
  ext_cpm_file
)

stopifnot(
  nrow(check_core) == 38606L,
  nrow(check_ext) == 38606L,
  identical(
    check_core$gene_symbol,
    feature_map$gene_symbol
  ),
  identical(
    check_ext$gene_symbol,
    feature_map$gene_symbol
  ),
  identical(
    names(check_core)[
      (length(feature_cols)+1L):ncol(check_core)
    ],
    core_order
  ),
  identical(
    names(check_ext)[
      (length(feature_cols)+1L):ncol(check_ext)
    ],
    extended_ids
  )
)

# -------------------------------------------------------------------------
# Usage note
# -------------------------------------------------------------------------

readme <- c(
  "Sepsis Atlas deconvolution reference v1",
  "",
  "Expression source:",
  "  Native final-QC SoupX-corrected integer RNA counts.",
  "  scMerge2-adjusted expression is not used.",
  "",
  "Aggregation:",
  "  sample x identity pseudobulk -> CPM ->",
  "  subject/condition balance -> project/condition median ->",
  "  equal healthy/disease weighting within project when available ->",
  "  equal project weighting by final project median.",
  "",
  "Reference tiers:",
  "  CORE contains the 20 pre-frozen, cross-project-supported identities.",
  "  EXTENDED contains CORE plus resolved weaker-support identities that",
  "  retain at least one eligible project after sample/project QC.",
  "  Deferred/unresolved/ambiguous identities are not published as",
  "  deconvolution targets, even if previously marked EXTENDED_ONLY.",
  "",
  "Recommended use:",
  "  CORE_CPM is the default bulk-RNA deconvolution reference.",
  "  EXTENDED is exploratory.",
  "  log2CPM versions are provided for correlation/visualization-oriented use.",
  "",
  "Differential expression:",
  "  Do not use these aggregated matrices for DE inference.",
  "  Use native processed counts with sample/patient-aware pseudobulk models."
)

writeLines(
  readme,
  file.path(
    out,
    "README_SepsisAtlas_deconvolution_reference_v1.txt"
  )
)

# -------------------------------------------------------------------------
# SHA256
# -------------------------------------------------------------------------

sha_file <- file.path(
  out,
  "SHA256SUMS_P05E6.txt"
)

if(file.exists(sha_file)) {
  unlink(sha_file)
}

oldwd <- getwd()
setwd(out)

files <- sort(
  list.files(
    ".",
    recursive=TRUE,
    full.names=FALSE
  )
)

files <- files[
  files != "SHA256SUMS_P05E6.txt"
]

sha <- vapply(
  files,
  function(f) {

    z <- system2(
      "sha256sum",
      args=f,
      stdout=TRUE,
      stderr=TRUE
    )

    st <- attr(
      z,
      "status"
    )

    if(
      !is.null(st) &&
      st != 0L
    ) {
      stop("sha256sum failed: ", f)
    }

    z[[1]]
  },
  character(1)
)

writeLines(
  sha,
  "SHA256SUMS_P05E6.txt"
)

chk <- system2(
  "sha256sum",
  args=c(
    "-c",
    "SHA256SUMS_P05E6.txt"
  ),
  stdout=TRUE,
  stderr=TRUE
)

st <- attr(
  chk,
  "status"
)

if(
  !is.null(st) &&
  st != 0L
) {
  cat(
    paste(
      chk,
      collapse="\n"
    ),
    "\n"
  )
  stop("P05E6 SHA256 verification failed")
}

setwd(oldwd)

# -------------------------------------------------------------------------
# Final summary
# -------------------------------------------------------------------------

core_min_projects <- min(
  final_support[
    published_tier == "CORE",
    n_eligible_projects
  ]
)

summary <- c(
  "===== P05E6 DECONVOLUTION REFERENCE RECOVERY =====",
  "status=PASS",
  "recovered_from_saved_sample_pseudobulk=YES",
  "final_QC_RDS_reread=NO",
  "pseudobulk_recomputed=NO",
  "expression_source=native_final_QC_SoupX_corrected_integer_RNA_counts",
  "scMerge2_expression_used=NO",
  "reference_genes=38606",
  paste0(
    "CORE_identities=",
    length(core_order)
  ),
  paste0(
    "EXTENDED_ONLY_published=",
    length(extended_only_ids)
  ),
  paste0(
    "EXTENDED_total_identities=",
    length(extended_ids)
  ),
  paste0(
    "excluded_identities=",
    nrow(excluded_policy)
  ),
  paste0(
    "CORE_min_eligible_projects=",
    core_min_projects
  ),
  "unresolved_deferred_ambiguous_deconvolution_targets=EXCLUDED",
  "sample_identity_min_cells=20",
  "project_identity_min_cells=50",
  "repeated_samples_balanced_within_subject=YES",
  "healthy_disease_balanced_within_project_when_available=YES",
  "projects_equal_weight_in_final_reference=YES",
  "final_profile_total_CPM=1000000",
  "gzip_integrity=PASS",
  "final_reference_reread_order=PASS",
  "SHA256=PASS",
  "expression_modified=NO",
  "annotation_modified=NO",
  "P05B_modified=NO",
  "cellranger_count=NOT_ACCESSED"
)

summary_file <- file.path(
  out,
  "P05E6_SUMMARY.txt"
)

writeLines(
  summary,
  summary_file
)

# -------------------------------------------------------------------------
# Root handoff files
# -------------------------------------------------------------------------

handoff <- c(
  "P05E6_SUMMARY.txt",
  "SepsisAtlas_deconvolution_reference_identities_v1.tsv",
  "SepsisAtlas_deconvolution_reference_support_v1.tsv",
  "SepsisAtlas_deconvolution_excluded_identities_v1.tsv",
  "deconvolution_identity_blood_fraction_support.tsv"
)

for(f in handoff) {

  src <- file.path(
    out,
    f
  )

  dst <- file.path(
    root,
    paste0(
      "P05E6_",
      sub(
        "^P05E6_",
        "",
        f
      )
    )
  )

  file.copy(
    src,
    dst,
    overwrite=TRUE
  )
}

cat(
  readLines(summary_file),
  sep="\n"
)

cat(
  "\n\n===== PUBLISHED REFERENCE SUPPORT =====\n"
)

print(
  final_support[
    order(
      factor(
        published_tier,
        levels=c(
          "CORE",
          "EXTENDED_ONLY"
        )
      ),
      reference_identity
    )
  ]
)

cat(
  "\n===== EXCLUDED IDENTITIES =====\n"
)

print(
  excluded_policy[
    ,
    .(
      reference_identity,
      recommended_reference,
      n_eligible_projects,
      publish_status
    )
  ]
)

cat(
  "\nOUT_DIR=",
  out,
  "\n",
  sep=""
)
