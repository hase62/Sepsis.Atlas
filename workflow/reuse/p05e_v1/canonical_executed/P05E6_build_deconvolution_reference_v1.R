#!/usr/bin/env Rscript

# P05E6: build study- and condition-balanced deconvolution references
# from native final-QC SoupX-corrected integer RNA counts.
#
# Design:
#   cell -> sample x core-identity pseudobulk counts
#        -> CPM
#        -> subject x condition x identity average across repeated samples
#        -> project x condition x identity median across subjects
#        -> equal condition average within project where both are available
#        -> median across eligible projects
#        -> rescale each final identity profile to 1e6 total CPM
#
# CORE and EXTENDED identity membership is frozen by P05E2b.
# No scMerge2-adjusted expression is used.
# No expression or annotation is modified.
# cellranger_count/ is not accessed.

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".", mustWork=TRUE)

Sys.setenv(
  OMP_NUM_THREADS="1",
  OPENBLAS_NUM_THREADS="1",
  MKL_NUM_THREADS="1",
  NUMEXPR_NUM_THREADS="1"
)

source(file.path(
  root,
  "full_atlas_primary_integration_v1",
  "00_common.R"
))

required <- c(
  "data.table",
  "Matrix",
  "matrixStats"
)

missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly=TRUE)
]

if(length(missing)) {
  stop("Missing packages: ", paste(missing, collapse=", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

out <- file.path(
  root,
  "publication",
  "scientific_data",
  paste0("reuse_reference_deconvolution_v1__", tag)
)

dir.create(out, recursive=TRUE)

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

latest_dir <- function(base, pattern) {
  x <- list.dirs(base, recursive=FALSE, full.names=TRUE)
  x <- sort(x[grepl(pattern, basename(x))])

  if(!length(x)) {
    stop("No directory matching: ", pattern)
  }

  tail(x, 1L)
}

write_gzip_tsv <- function(dt, gzfile_path) {

  plain <- sub("\\.gz$", "", gzfile_path)

  if(identical(plain, gzfile_path)) {
    stop("Output must end in .gz: ", gzfile_path)
  }

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

write_gzip_mtx <- function(mat, gzfile_path) {

  plain <- sub("\\.gz$", "", gzfile_path)

  if(identical(plain, gzfile_path)) {
    stop("MatrixMarket output must end in .gz")
  }

  Matrix::writeMM(
    mat,
    plain
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

row_median_sparse_subset <- function(mat, idx) {

  if(length(idx) == 1L) {
    return(as.numeric(mat[, idx]))
  }

  matrixStats::rowMedians(
    as.matrix(
      mat[, idx, drop=FALSE]
    )
  )
}

normalize_profile_to_cpm <- function(x) {
  s <- sum(x)

  if(!is.finite(s) || s <= 0) {
    stop("Non-positive reference profile sum")
  }

  x / s * 1e6
}

# -------------------------------------------------------------------------
# Canonical inputs
# -------------------------------------------------------------------------

sd_base <- file.path(
  root,
  "publication",
  "scientific_data"
)

p05e5 <- latest_dir(
  sd_base,
  "reuse_reference_rich_metadata_v1__[0-9]{8}_[0-9]{6}$"
)

p05e2b <- latest_dir(
  sd_base,
  "reuse_reference_plan_refined_v1__[0-9]{8}_[0-9]{6}$"
)

p05b <- file.path(
  sd_base,
  "processed_data_release_candidate_v1__20260819_194718"
)

rich_cell_file <- file.path(
  p05e5,
  "SepsisAtlas_rich_cell_metadata_v1.tsv.gz"
)

policy_file <- file.path(
  p05e2b,
  "deconvolution_reference_policy_refined.tsv"
)

feature_map_file <- file.path(
  p05b,
  "metadata",
  "SepsisAtlas_public_feature_mapping_v1__processed_data_release_candidate_v1__20260819_194718.tsv"
)

stopifnot(
  file.exists(rich_cell_file),
  file.exists(policy_file),
  file.exists(feature_map_file)
)

cells <- fread(rich_cell_file)
policy <- fread(policy_file)
feature_map <- fread(feature_map_file)

required_cell <- c(
  "global_cell",
  "project_id",
  "library_key",
  "analysis_sample_id",
  "analysis_subject_id",
  "clinical_domain_std",
  "blood_fraction_std",
  "final_compartment_v1",
  "atlas_core_identity_v1"
)

required_policy <- c(
  "final_compartment_v1",
  "atlas_core_identity_v1",
  "recommended_reference"
)

stopifnot(
  all(required_cell %in% names(cells)),
  all(required_policy %in% names(policy)),
  nrow(cells) == 665816L,
  !anyDuplicated(cells$global_cell),
  nrow(feature_map) == 38606L
)

if(any(
  is.na(cells$analysis_sample_id) |
  !nzchar(cells$analysis_sample_id)
)) {
  stop("analysis_sample_id missing")
}

if(any(
  is.na(cells$analysis_subject_id) |
  !nzchar(cells$analysis_subject_id)
)) {
  stop("analysis_subject_id missing")
}

# -------------------------------------------------------------------------
# Frozen identity policy
# -------------------------------------------------------------------------

policy <- unique(
  policy[
    recommended_reference %in% c(
      "CORE",
      "EXTENDED_ONLY"
    ),
    .(
      final_compartment_v1,
      atlas_core_identity_v1,
      recommended_reference
    )
  ]
)

policy[, reference_identity :=
  paste(
    final_compartment_v1,
    atlas_core_identity_v1,
    sep="::"
  )
]

if(anyDuplicated(policy$reference_identity)) {
  stop("Duplicated reference_identity in policy")
}

n_core <- policy[
  recommended_reference == "CORE",
  .N
]

n_extended_only <- policy[
  recommended_reference == "EXTENDED_ONLY",
  .N
]

stopifnot(
  n_core == 20L,
  n_extended_only == 13L,
  nrow(policy) == 33L
)

# -------------------------------------------------------------------------
# Disease/healthy mapping, identical to integration-side semantic mapping
# -------------------------------------------------------------------------

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

cells <- merge(
  cells,
  condition_map,
  by="clinical_domain_std",
  all.x=TRUE,
  sort=FALSE
)

if(any(is.na(cells$condition_binary))) {
  stop(
    "Unmapped clinical_domain_std: ",
    paste(
      unique(
        cells[
          is.na(condition_binary),
          clinical_domain_std
        ]
      ),
      collapse=", "
    )
  )
}

# -------------------------------------------------------------------------
# Select cells belonging to frozen deconvolution identities
# -------------------------------------------------------------------------

target <- merge(
  cells,
  policy,
  by=c(
    "final_compartment_v1",
    "atlas_core_identity_v1"
  ),
  all=FALSE,
  sort=FALSE
)

stopifnot(nrow(target) > 0L)

# Exact original barcode reconstruction.
prefix <- paste0(
  target$library_key,
  "___"
)

if(!all(startsWith(
  target$global_cell,
  prefix
))) {
  stop("global_cell/library_key prefix mismatch")
}

target[, original_cell :=
  substring(
    global_cell,
    nchar(library_key) + 4L
  )
]

# -------------------------------------------------------------------------
# Sample x identity support and eligibility
# -------------------------------------------------------------------------

# Explicit sample-profile grouping avoids accidental collapsing across samples.
sample_support <- target[
  ,
  .(
    n_cells=.N,
    recommended_reference=
      recommended_reference[[1]],
    clinical_domain_std=
      clinical_domain_std[[1]],
    blood_fraction_std=
      blood_fraction_std[[1]]
  ),
  by=.(
    project_id,
    library_key,
    analysis_sample_id,
    analysis_subject_id,
    condition_binary,
    reference_identity
  )
]

sample_support[, sample_profile_eligible :=
  n_cells >= 20L
]

sample_support[, profile_id :=
  paste(
    project_id,
    analysis_sample_id,
    reference_identity,
    sep="::"
  )
]

if(anyDuplicated(sample_support$profile_id)) {
  stop("Duplicated sample pseudobulk profile_id")
}

# -------------------------------------------------------------------------
# Preflight project-level support after minimum sample-profile cell threshold
# -------------------------------------------------------------------------

eligible_sample_support <- sample_support[
  sample_profile_eligible == TRUE
]

project_support <- eligible_sample_support[
  ,
  .(
    n_cells=sum(n_cells),
    n_samples=uniqueN(analysis_sample_id),
    n_subjects=uniqueN(analysis_subject_id),
    n_healthy_samples=uniqueN(
      analysis_sample_id[
        condition_binary == "healthy"
      ]
    ),
    n_disease_samples=uniqueN(
      analysis_sample_id[
        condition_binary == "disease"
      ]
    ),
    recommended_reference=
      recommended_reference[[1]]
  ),
  by=.(
    reference_identity,
    project_id
  )
]

project_support[, project_profile_eligible :=
  n_cells >= 50L &
  n_subjects >= 1L
]

identity_support <- project_support[
  ,
  .(
    n_projects=.N,
    n_eligible_projects=sum(
      project_profile_eligible == TRUE
    ),
    n_cells=sum(n_cells),
    n_samples=sum(n_samples),
    n_subject_project_pairs=sum(n_subjects),
    max_project_cell_fraction=
      max(n_cells) / sum(n_cells),
    recommended_reference=
      recommended_reference[[1]]
  ),
  by=reference_identity
]

identity_support <- merge(
  policy,
  identity_support,
  by=c(
    "reference_identity",
    "recommended_reference"
  ),
  all.x=TRUE,
  sort=FALSE
)

identity_support[
  is.na(n_eligible_projects),
  n_eligible_projects := 0L
]

fwrite(
  sample_support,
  file.path(
    out,
    "deconvolution_sample_identity_support.tsv"
  ),
  sep="\t"
)

fwrite(
  project_support,
  file.path(
    out,
    "deconvolution_project_identity_support.tsv"
  ),
  sep="\t"
)

fwrite(
  identity_support,
  file.path(
    out,
    "deconvolution_identity_support.tsv"
  ),
  sep="\t"
)

# Blood-fraction support is descriptive, not a weighting variable.
blood_support <- target[
  ,
  .(
    n_cells=.N,
    n_samples=uniqueN(analysis_sample_id),
    n_subjects=uniqueN(analysis_subject_id),
    n_projects=uniqueN(project_id)
  ),
  by=.(
    reference_identity,
    recommended_reference,
    blood_fraction_std
  )
]

fwrite(
  blood_support,
  file.path(
    out,
    "deconvolution_identity_blood_fraction_support.tsv"
  ),
  sep="\t"
)

bad_core <- identity_support[
  recommended_reference == "CORE" &
  n_eligible_projects < 3L
]

if(nrow(bad_core)) {

  fwrite(
    bad_core,
    file.path(
      out,
      "deconvolution_CORE_support_failure.tsv"
    ),
    sep="\t"
  )

  # Root handoff before stopping.
  file.copy(
    file.path(
      out,
      "deconvolution_CORE_support_failure.tsv"
    ),
    file.path(
      root,
      "P05E6_deconvolution_CORE_support_failure.tsv"
    ),
    overwrite=TRUE
  )

  stop(
    "CORE identities with <3 eligible projects after sample-profile QC: ",
    paste(
      bad_core$reference_identity,
      collapse=", "
    )
  )
}

# -------------------------------------------------------------------------
# final-QC RDS map and feature-order audit
# -------------------------------------------------------------------------

rds_files <- Sys.glob(
  file.path(
    root,
    "pre_integration",
    "GSE*",
    "rds_final_qc",
    "*.rds"
  )
)

stopifnot(length(rds_files) == 158L)

rds_map <- rbindlist(
  lapply(rds_files, function(f) {
    data.table(
      project_id=basename(dirname(dirname(f))),
      library_key=sub(
        "__final_qc\\.rds$",
        "",
        basename(f)
      ),
      rds=f
    )
  })
)

stopifnot(
  nrow(rds_map) == 158L,
  !anyDuplicated(
    paste(
      rds_map$project_id,
      rds_map$library_key,
      sep="|||"
    )
  )
)

first_obj <- readRDS(rds_map$rds[[1]])
feature_ref <- rownames(get_rna_counts(first_obj))
rm(first_obj)
invisible(gc())

stopifnot(length(feature_ref) == 38606L)

# P04B states gene_symbol is the exact Atlas feature name/order.
if(!"gene_symbol" %in% names(feature_map)) {
  stop("gene_symbol missing from public feature mapping")
}

if(!identical(
  as.character(feature_map$gene_symbol),
  feature_ref
)) {
  stop("P04B feature map order does not match final-QC RDS feature order")
}

# -------------------------------------------------------------------------
# Build sample x identity native pseudobulk counts
# -------------------------------------------------------------------------

eligible_keys <- eligible_sample_support[
  ,
  .(
    project_id,
    library_key,
    analysis_sample_id,
    analysis_subject_id,
    condition_binary,
    reference_identity,
    recommended_reference,
    n_cells,
    profile_id
  )
]

target_eligible <- merge(
  target,
  eligible_keys[
    ,
    .(
      project_id,
      library_key,
      reference_identity,
      profile_id,
      sample_profile_eligible_n_cells=n_cells
    )
  ],
  by=c(
    "project_id",
    "library_key",
    "reference_identity"
  ),
  all=FALSE,
  sort=FALSE
)

pb_list <- vector(
  "list",
  nrow(eligible_keys)
)

names(pb_list) <- eligible_keys$profile_id

for(i in seq_len(nrow(rds_map))) {

  pid <- rds_map$project_id[[i]]
  key <- rds_map$library_key[[i]]

  mm <- target_eligible[
    project_id == pid &
    library_key == key
  ]

  if(!nrow(mm)) {
    next
  }

  obj <- readRDS(
    rds_map$rds[[i]]
  )

  counts <- get_rna_counts(obj)

  if(!identical(
    rownames(counts),
    feature_ref
  )) {
    stop(
      "Feature order mismatch: ",
      pid, "::", key
    )
  }

  identities <- unique(
    mm$reference_identity
  )

  for(ref_id in identities) {

    mmi <- mm[
      reference_identity == ref_id
    ]

    if(uniqueN(mmi$profile_id) != 1L) {
      stop(
        "Expected one sample profile per library/identity: ",
        pid, "::", key, "::", ref_id
      )
    }

    profile_id <- mmi$profile_id[[1]]

    barcodes <- mmi$original_cell

    if(!all(barcodes %in% colnames(counts))) {
      stop(
        "Missing cells in final-QC RDS: ",
        profile_id
      )
    }

    pb <- Matrix::rowSums(
      counts[
        ,
        barcodes,
        drop=FALSE
      ]
    )

    if(any(!is.finite(pb)) || any(pb < 0)) {
      stop(
        "Invalid pseudobulk counts: ",
        profile_id
      )
    }

    pb_list[[profile_id]] <- Matrix::Matrix(
      pb,
      ncol=1L,
      sparse=TRUE
    )
  }

  rm(obj, counts)
  invisible(gc())

  cat(
    "PSEUDOBULK ",
    i, "/", nrow(rds_map),
    " ",
    pid, "::", key,
    "\n",
    sep=""
  )
}

# Profile metadata follows the frozen eligible_keys order.
profile_meta <- copy(eligible_keys)

# pb_list was preallocated and named; verify every eligible profile was filled.
filled <- vapply(
  pb_list,
  function(x) !is.null(x),
  logical(1)
)

if(!all(filled)) {
  stop(
    "Missing pseudobulk profiles: ",
    paste(
      names(pb_list)[!filled],
      collapse=", "
    )
  )
}

pb_counts <- do.call(
  cbind,
  pb_list
)

colnames(pb_counts) <- names(pb_list)
rownames(pb_counts) <- feature_ref

stopifnot(
  ncol(pb_counts) == nrow(eligible_keys),
  nrow(pb_counts) == 38606L,
  identical(
    colnames(pb_counts),
    eligible_keys$profile_id
  )
)

profile_meta <- copy(eligible_keys)

# -------------------------------------------------------------------------
# Save sample pseudobulk counts in cross-language Matrix Market format
# -------------------------------------------------------------------------

pb_mtx_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_sample_pseudobulk_counts_v1.mtx.gz"
)

write_gzip_mtx(
  pb_counts,
  pb_mtx_file
)

pb_meta_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_sample_pseudobulk_metadata_v1.tsv.gz"
)

write_gzip_tsv(
  profile_meta,
  pb_meta_file
)

feature_out <- file.path(
  out,
  "SepsisAtlas_deconvolution_features_v1.tsv.gz"
)

write_gzip_tsv(
  feature_map,
  feature_out
)

# -------------------------------------------------------------------------
# CPM normalization
# -------------------------------------------------------------------------

library_size <- Matrix::colSums(
  pb_counts
)

if(any(!is.finite(library_size)) ||
   any(library_size <= 0)) {
  stop("Invalid sample pseudobulk library size")
}

scale_factor <- 1e6 / library_size

pb_cpm <- pb_counts %*%
  Matrix::Diagonal(
    x=scale_factor
  )

# -------------------------------------------------------------------------
# Repeated samples -> subject x condition x identity
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
colnames(subject_cpm) <- subject_levels
rownames(subject_cpm) <- feature_ref

subject_meta <- profile_meta[
  ,
  .(
    project_id=project_id[[1]],
    analysis_subject_id=
      analysis_subject_id[[1]],
    condition_binary=
      condition_binary[[1]],
    reference_identity=
      reference_identity[[1]],
    recommended_reference=
      recommended_reference[[1]],
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
# Subject profiles -> project x condition x identity median
# -------------------------------------------------------------------------

pc_groups <- subject_meta[
  ,
  .(
    subject_indices=list(.I),
    n_subjects=.N,
    n_samples=sum(n_samples),
    n_cells=sum(n_cells),
    recommended_reference=
      recommended_reference[[1]]
  ),
  by=.(
    project_id,
    condition_binary,
    reference_identity
  )
]

pc_profiles <- matrix(
  0,
  nrow=length(feature_ref),
  ncol=nrow(pc_groups),
  dimnames=list(
    feature_ref,
    paste(
      pc_groups$project_id,
      pc_groups$condition_binary,
      pc_groups$reference_identity,
      sep="::"
    )
  )
)

for(i in seq_len(nrow(pc_groups))) {

  idx <- pc_groups$subject_indices[[i]]

  pc_profiles[, i] <-
    row_median_sparse_subset(
      subject_cpm,
      idx
    )
}

pc_groups[, project_condition_profile_id :=
  colnames(pc_profiles)
]

# -------------------------------------------------------------------------
# Equalize healthy/disease within project where both are present
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
      sort(unique(condition_binary)),
      collapse=","
    ),
    recommended_reference=
      recommended_reference[[1]]
  ),
  by=.(
    project_id,
    reference_identity
  )
]

project_profiles <- matrix(
  0,
  nrow=length(feature_ref),
  ncol=nrow(project_groups),
  dimnames=list(
    feature_ref,
    paste(
      project_groups$project_id,
      project_groups$reference_identity,
      sep="::"
    )
  )
)

for(i in seq_len(nrow(project_groups))) {

  idx <- project_groups$pc_indices[[i]]

  # Equal condition weighting when healthy + disease are both available.
  project_profiles[, i] <-
    rowMeans(
      pc_profiles[
        ,
        idx,
        drop=FALSE
      ]
    )
}

project_groups[, project_profile_id :=
  colnames(project_profiles)
]

# Project eligibility frozen from metadata-only preflight.
project_groups <- merge(
  project_groups,
  project_support[
    ,
    .(
      reference_identity,
      project_id,
      project_profile_eligible
    )
  ],
  by=c(
    "reference_identity",
    "project_id"
  ),
  all.x=TRUE,
  sort=FALSE
)

project_groups <- project_groups[
  match(
    colnames(project_profiles),
    project_profile_id
  )
]

stopifnot(
  !any(is.na(
    project_groups$project_profile_eligible
  ))
)

# -------------------------------------------------------------------------
# Final study-balanced identity profiles
# -------------------------------------------------------------------------

identity_order_extended <- policy[
  order(
    factor(
      recommended_reference,
      levels=c(
        "CORE",
        "EXTENDED_ONLY"
      )
    ),
    final_compartment_v1,
    atlas_core_identity_v1
  ),
  reference_identity
]

final_profiles <- matrix(
  NA_real_,
  nrow=length(feature_ref),
  ncol=length(identity_order_extended),
  dimnames=list(
    feature_ref,
    identity_order_extended
  )
)

final_support <- vector(
  "list",
  length(identity_order_extended)
)

for(i in seq_along(identity_order_extended)) {

  ref_id <- identity_order_extended[[i]]

  idx <- which(
    project_groups$reference_identity == ref_id &
    project_groups$project_profile_eligible == TRUE
  )

  if(!length(idx)) {
    stop(
      "No eligible project profiles for ",
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

  x <- normalize_profile_to_cpm(x)

  final_profiles[, i] <- x

  pol <- policy[
    reference_identity == ref_id
  ]

  final_support[[i]] <- data.table(
    reference_identity=ref_id,
    final_compartment_v1=
      pol$final_compartment_v1[[1]],
    atlas_core_identity_v1=
      pol$atlas_core_identity_v1[[1]],
    recommended_reference=
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

stopifnot(
  all(is.finite(final_profiles)),
  all(final_profiles >= 0)
)

profile_sums <- colSums(
  final_profiles
)

if(max(abs(profile_sums - 1e6)) > 1) {
  stop("Final CPM columns do not sum to ~1e6")
}

# -------------------------------------------------------------------------
# CORE and EXTENDED matrices
# -------------------------------------------------------------------------

core_ids <- policy[
  recommended_reference == "CORE",
  reference_identity
]

core_ids <- identity_order_extended[
  identity_order_extended %in% core_ids
]

stopifnot(length(core_ids) == 20L)

extended_ids <- identity_order_extended

stopifnot(length(extended_ids) == 33L)

core_cpm <- final_profiles[
  ,
  core_ids,
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
# Attach public feature identifiers
# -------------------------------------------------------------------------

feature_prefix <- feature_map[
  ,
  .(
    feature_index,
    ensembl_gene_id,
    source_ensembl_id,
    gene_symbol,
    feature_type,
    mapping_source,
    mapping_rule
  )
]

make_reference_table <- function(mat) {
  z <- as.data.table(mat)
  setnames(
    z,
    colnames(mat)
  )
  cbind(
    feature_prefix,
    z
  )
}

core_cpm_dt <- make_reference_table(
  core_cpm
)

core_log_dt <- make_reference_table(
  core_log2cpm
)

ext_cpm_dt <- make_reference_table(
  extended_cpm
)

ext_log_dt <- make_reference_table(
  extended_log2cpm
)

# -------------------------------------------------------------------------
# Write final public-facing reference matrices
# -------------------------------------------------------------------------

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
  core_cpm_dt,
  core_cpm_file
)

write_gzip_tsv(
  core_log_dt,
  core_log_file
)

write_gzip_tsv(
  ext_cpm_dt,
  ext_cpm_file
)

write_gzip_tsv(
  ext_log_dt,
  ext_log_file
)

# -------------------------------------------------------------------------
# Project-level normalized profiles and support metadata
# -------------------------------------------------------------------------

project_profile_dt <- make_reference_table(
  project_profiles
)

project_profile_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_project_condition_balanced_CPM_v1.tsv.gz"
)

write_gzip_tsv(
  project_profile_dt,
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

final_support_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_reference_support_v1.tsv"
)

fwrite(
  final_support,
  final_support_file,
  sep="\t",
  quote=TRUE,
  na="NA"
)

identity_file <- file.path(
  out,
  "SepsisAtlas_deconvolution_reference_identities_v1.tsv"
)

identity_table <- merge(
  policy,
  final_support[
    ,
    .(
      reference_identity,
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
    core_ids
  )
]

identity_table[, EXTENDED_column_order :=
  match(
    reference_identity,
    extended_ids
  )
]

setorder(
  identity_table,
  EXTENDED_column_order
)

fwrite(
  identity_table,
  identity_file,
  sep="\t",
  quote=TRUE,
  na="NA"
)

# -------------------------------------------------------------------------
# Re-read final TSV.gz matrices and verify feature/column order
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
    feature_ref
  ),
  identical(
    check_ext$gene_symbol,
    feature_ref
  ),
  identical(
    names(check_core)[
      (ncol(feature_prefix)+1L):ncol(check_core)
    ],
    core_ids
  ),
  identical(
    names(check_ext)[
      (ncol(feature_prefix)+1L):ncol(check_ext)
    ],
    extended_ids
  )
)

# -------------------------------------------------------------------------
# README / usage note
# -------------------------------------------------------------------------

readme <- c(
  "Sepsis Atlas deconvolution reference v1",
  "",
  "Expression source:",
  "  Native final-QC SoupX-corrected integer RNA counts.",
  "  No scMerge2-adjusted expression is used for deconvolution.",
  "",
  "Aggregation:",
  "  1. Sum native counts within sample x Atlas core identity.",
  "  2. Exclude sample-identity profiles with <20 cells.",
  "  3. Convert each sample pseudobulk to CPM.",
  "  4. Average repeated samples within subject x condition x identity.",
  "  5. Median across subjects within project x condition x identity.",
  "  6. Equal-weight healthy/disease profiles within project when both exist.",
  "  7. Median across eligible projects (>=50 cells after sample-profile QC).",
  "  8. Rescale each final identity profile to total CPM = 1e6.",
  "",
  "Recommended use:",
  "  CORE_CPM is the default bulk-RNA deconvolution reference.",
  "  CORE_log2CPM is supplied for correlation/visualization-oriented methods.",
  "  EXTENDED matrices include identities with weaker cross-project support",
  "  and should be treated as exploratory.",
  "",
  "The reference is study-balanced; it is not weighted by raw cell counts.",
  "Source blood-fraction support is provided separately.",
  "",
  "Differential expression:",
  "  Do not use these aggregated reference matrices for DE inference.",
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

    st <- attr(z, "status")

    if(!is.null(st) && st != 0L) {
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

if(!is.null(st) && st != 0L) {
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
# Summary
# -------------------------------------------------------------------------

summary <- c(
  "===== P05E6 DECONVOLUTION REFERENCE =====",
  "status=PASS",
  "expression_source=native_final_QC_SoupX_corrected_integer_RNA_counts",
  "scMerge2_expression_used=NO",
  "reference_genes=38606",
  paste0(
    "CORE_identities=",
    length(core_ids)
  ),
  paste0(
    "EXTENDED_identities=",
    length(extended_ids)
  ),
  paste0(
    "eligible_sample_identity_profiles=",
    nrow(eligible_keys)
  ),
  paste0(
    "subject_condition_profiles=",
    ncol(subject_cpm)
  ),
  paste0(
    "project_condition_profiles=",
    ncol(pc_profiles)
  ),
  paste0(
    "project_balanced_profiles=",
    ncol(project_profiles)
  ),
  paste0(
    "CORE_min_eligible_projects=",
    min(
      final_support[
        recommended_reference == "CORE",
        n_eligible_projects
      ]
    )
  ),
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
# Root handoff copies
# -------------------------------------------------------------------------

handoff <- c(
  "P05E6_SUMMARY.txt",
  "SepsisAtlas_deconvolution_reference_identities_v1.tsv",
  "SepsisAtlas_deconvolution_reference_support_v1.tsv",
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
  "\n\n===== FINAL REFERENCE SUPPORT =====\n"
)

print(
  final_support[
    order(
      factor(
        recommended_reference,
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
