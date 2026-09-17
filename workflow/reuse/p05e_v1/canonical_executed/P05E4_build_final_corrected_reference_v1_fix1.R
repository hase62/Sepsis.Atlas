#!/usr/bin/env Rscript

# P05E4 fix1: explicit data.table logical-column filtering for resolved_core
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
  "scMerge","BiocParallel","BiocSingular",
  "BiocNeighbors","irlba","data.table","Matrix"
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
  library(BiocParallel)
  library(BiocSingular)
  library(BiocNeighbors)
  library(irlba)
})

tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

out <- file.path(
  root,
  "publication",
  "scientific_data",
  paste0("reuse_reference_corrected_expression_v1__", tag)
)

matrix_dir <- file.path(out, "matrices")
model_dir <- file.path(out, "models")
metadata_dir <- file.path(out, "metadata")

dir.create(matrix_dir, recursive=TRUE)
dir.create(model_dir, recursive=TRUE)
dir.create(metadata_dir, recursive=TRUE)

# ============================================================
# Canonical P05B metadata
# ============================================================

p05b <- file.path(
  root,
  "publication",
  "scientific_data",
  "processed_data_release_candidate_v1__20260819_194718"
)

cell_file <- file.path(
  p05b,
  "metadata",
  "SepsisAtlas_cell_metadata_v1__processed_data_release_candidate_v1__20260819_194718.tsv.gz"
)

lib_file <- file.path(
  p05b,
  "metadata",
  "SepsisAtlas_library_metadata_with_source_accessions_v1__processed_data_release_candidate_v1__20260819_194718.tsv"
)

stopifnot(file.exists(cell_file), file.exists(lib_file))

cells <- fread(cell_file)
libs <- fread(lib_file)

required_cell_cols <- c(
  "global_cell",
  "project_id",
  "library_key",
  "final_compartment_v1",
  "atlas_core_identity_v1",
  "atlas_axis_v1",
  "atlas_state_v1"
)

stopifnot(
  all(required_cell_cols %in% names(cells)),
  nrow(cells) == 665816L,
  !anyDuplicated(cells$global_cell)
)

required_lib_cols <- c(
  "project_id",
  "library_key",
  "clinical_domain_std"
)

stopifnot(
  all(required_lib_cols %in% names(libs)),
  nrow(libs) == 158L
)

# ============================================================
# Latest P05E3b parameter freeze
# ============================================================

p3dirs <- list.dirs(
  file.path(root, "publication", "scientific_data"),
  recursive=FALSE,
  full.names=TRUE
)

p3dirs <- sort(
  p3dirs[
    grepl(
      "reuse_reference_scmerge2_selection_v1__[0-9]{8}_[0-9]{6}$",
      basename(p3dirs)
    )
  ]
)

stopifnot(length(p3dirs) >= 1L)

p3 <- tail(p3dirs, 1L)

plan <- fread(
  file.path(
    p3,
    "final_corrected_reference_plan.tsv"
  )
)

expected_compartments <- c(
  "T_NK",
  "Monocyte_DC",
  "B_plasma",
  "Neutrophil",
  "Platelet_megakaryocyte"
)

stopifnot(
  nrow(plan) == 5L,
  setequal(plan$compartment, expected_compartments)
)

plan <- plan[
  match(expected_compartments, compartment)
]

# ============================================================
# condition mapping
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

lib_cond <- merge(
  libs[, .(
    project_id,
    library_key,
    clinical_domain_std
  )],
  condition_map,
  by="clinical_domain_std",
  all.x=TRUE
)

if(any(is.na(lib_cond$condition_binary))) {
  stop(
    "Unmapped clinical_domain_std: ",
    paste(
      unique(
        lib_cond[
          is.na(condition_binary),
          clinical_domain_std
        ]
      ),
      collapse=", "
    )
  )
}

# ============================================================
# final-QC RDS map
# ============================================================

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

# ============================================================
# Feature space: same 3000 consensus HVGs + human scSEG
# ============================================================

hvg_file <- file.path(
  root,
  "pre_integration",
  "pilot_unintegrated_large_v1",
  "large_pilot_consensus_hvg_3000.tsv"
)

stopifnot(file.exists(hvg_file))

hvg <- fread(hvg_file)
stopifnot("feature" %in% names(hvg))

consensus_hvg <- unique(as.character(hvg$feature))
stopifnot(length(consensus_hvg) == 3000L)

first_obj <- readRDS(rds_map$rds[[1]])
feature_reference <- rownames(get_rna_counts(first_obj))
rm(first_obj)
invisible(gc())

data(
  "segList",
  package="scMerge",
  envir=environment()
)

human_seg <- intersect(
  unique(as.character(segList$human$human_scSEG)),
  feature_reference
)

stopifnot(length(human_seg) >= 100L)

input_features <- unique(c(
  consensus_hvg,
  human_seg
))

# ============================================================
# Reference-universe helpers
# ============================================================

is_resolved_core <- function(x) {
  ok <- !is.na(x) & nzchar(trimws(as.character(x)))
  bad <- grepl(
    "unresolved|ambiguous|^Deferred_",
    as.character(x),
    ignore.case=TRUE
  )
  ok & !bad
}

clean_component <- function(x, missing_label) {
  x <- as.character(x)
  bad <- is.na(x) | !nzchar(trimws(x))
  x[bad] <- missing_label
  x
}


matrix_all_finite <- function(mat, chunk_cols=5000L) {
  n <- ncol(mat)
  starts <- seq.int(1L, n, by=chunk_cols)

  for(a in starts) {
    b <- min(n, a + chunk_cols - 1L)

    if(any(
      !is.finite(
        mat[, a:b, drop=FALSE]
      )
    )) {
      return(FALSE)
    }
  }

  TRUE
}

make_integration_celltype <- function(meta, compartment) {

  core_lab <- clean_component(
    meta$atlas_core_identity_v1,
    "core_unresolved"
  )

  axis_lab <- clean_component(
    meta$atlas_axis_v1,
    "axis_not_available"
  )

  state_lab <- clean_component(
    meta$atlas_state_v1,
    "state_not_available"
  )

  if(compartment == "B_plasma") {
    return(core_lab)
  }

  if(compartment == "Platelet_megakaryocyte") {
    return(
      paste(
        core_lab,
        axis_lab,
        state_lab,
        sep="|||"
      )
    )
  }

  if(compartment == "Neutrophil") {
    return(
      paste(
        axis_lab,
        state_lab,
        sep="|||"
      )
    )
  }

  core_lab
}

prepare_reference_meta <- function(compartment) {

  all_meta <- cells[
    final_compartment_v1 == compartment,
    .(
      global_cell,
      project_id,
      library_key,
      final_compartment_v1,
      atlas_core_identity_v1,
      atlas_axis_v1,
      atlas_state_v1
    )
  ]

  if(!nrow(all_meta)) {
    stop("No cells for compartment: ", compartment)
  }

  all_meta[, resolved_core :=
    is_resolved_core(atlas_core_identity_v1)
  ]

  excluded_core <- all_meta[
    resolved_core == FALSE,
    .(
      global_cell,
      project_id,
      library_key,
      final_compartment_v1,
      exclusion_reason="unresolved_or_deferred_core_identity"
    )
  ]

  meta <- all_meta[resolved_core == TRUE]

  meta <- merge(
    meta,
    lib_cond[, .(
      project_id,
      library_key,
      clinical_domain_std,
      condition_binary
    )],
    by=c("project_id","library_key"),
    all.x=TRUE,
    sort=FALSE
  )

  if(any(is.na(meta$condition_binary))) {
    stop(compartment, ": missing condition")
  }

  meta[, integration_celltype :=
    make_integration_celltype(.SD, compartment)
  ]

  meta[, integration_celltype_n := .N,
       by=integration_celltype]

  excluded_rare_group <- meta[
    integration_celltype_n < 10L,
    .(
      global_cell,
      project_id,
      library_key,
      final_compartment_v1,
      exclusion_reason="integration_celltype_total_lt10"
    )
  ]

  meta <- meta[
    integration_celltype_n >= 10L
  ]

  meta[, batch_condition_n := .N,
       by=.(project_id, condition_binary)]

  excluded_tiny_stratum <- meta[
    batch_condition_n < 10L,
    .(
      global_cell,
      project_id,
      library_key,
      final_compartment_v1,
      exclusion_reason="project_condition_stratum_lt10"
    )
  ]

  meta <- meta[
    batch_condition_n >= 10L
  ]

  prefix <- paste0(
    meta$library_key,
    "___"
  )

  if(!all(startsWith(meta$global_cell, prefix))) {
    stop(compartment, ": global_cell/library_key prefix mismatch")
  }

  meta[, original_cell :=
    substring(
      global_cell,
      nchar(library_key) + 4L
    )
  ]

  excluded <- rbindlist(
    list(
      excluded_core,
      excluded_rare_group,
      excluded_tiny_stratum
    ),
    fill=TRUE
  )

  list(
    all_n=nrow(all_meta),
    meta=meta,
    excluded=excluded
  )
}

# ============================================================
# Preflight ALL five compartments before any scMerge2 fit
# ============================================================

prepared <- list()
preflight <- list()

for(compartment in expected_compartments) {

  z <- prepare_reference_meta(compartment)
  prepared[[compartment]] <- z

  m <- z$meta

  pf <- data.table(
    compartment=compartment,
    n_final_atlas_cells=z$all_n,
    n_reference_cells=nrow(m),
    n_excluded_cells=nrow(z$excluded),
    n_projects=uniqueN(m$project_id),
    n_libraries=uniqueN(m$library_key),
    n_conditions=uniqueN(m$condition_binary),
    n_integration_celltypes=uniqueN(m$integration_celltype),
    min_integration_celltype_n=min(
      m[, .N, by=integration_celltype]$N
    ),
    min_project_condition_n=min(
      m[, .N, by=.(project_id, condition_binary)]$N
    )
  )

  if(
    pf$n_reference_cells < 100L ||
    pf$n_projects < 2L ||
    pf$n_conditions < 2L ||
    pf$n_integration_celltypes < 1L ||
    pf$min_integration_celltype_n < 10L ||
    pf$min_project_condition_n < 10L
  ) {
    print(pf)
    stop("P05E4 preflight failed for ", compartment)
  }

  preflight[[length(preflight)+1L]] <- pf
}

preflight <- rbindlist(preflight)

fwrite(
  preflight,
  file.path(out, "corrected_reference_preflight.tsv"),
  sep="\t"
)

cat("\n===== P05E4 PREFLIGHT =====\n")
print(preflight)

# ============================================================
# Fit final selected models and RETAIN adjusted expression
# ============================================================

manifest_rows <- list()
gene_sets <- list()
all_excluded <- list()

stable_seed <- c(
  T_NK=2026082201L,
  Monocyte_DC=2026082202L,
  B_plasma=2026082203L,
  Neutrophil=2026082204L,
  Platelet_megakaryocyte=2026082205L
)

for(ci in seq_len(nrow(plan))) {

  compartment <- as.character(plan$compartment[[ci]])
  ruvK <- as.integer(plan$ruvK[[ci]])
  kPB <- as.integer(plan$k_pseudoBulk[[ci]])

  cat(
    "\n============================================================\n",
    "P05E4 FINAL CORRECTED REFERENCE: ", compartment, "\n",
    "ruvK=", ruvK, "; k_pseudoBulk=", kPB, "\n",
    "============================================================\n",
    sep=""
  )

  z <- prepared[[compartment]]
  meta <- copy(z$meta)

  if(nrow(z$excluded)) {
    all_excluded[[length(all_excluded)+1L]] <- z$excluded
  }

  lib_pairs <- unique(
    meta[, .(
      project_id,
      library_key
    )]
  )

  counts_list <- vector(
    "list",
    nrow(lib_pairs)
  )

  for(i in seq_len(nrow(lib_pairs))) {

    pid <- lib_pairs$project_id[[i]]
    key <- lib_pairs$library_key[[i]]

    mm <- meta[
      project_id == pid &
      library_key == key
    ]

    rr <- rds_map[
      project_id == pid &
      library_key == key
    ]

    if(nrow(rr) != 1L) {
      stop(
        "RDS mapping failure: ",
        pid, " / ", key
      )
    }

    obj <- readRDS(rr$rds[[1]])
    counts <- get_rna_counts(obj)

    miss_f <- setdiff(
      input_features,
      rownames(counts)
    )

    if(length(miss_f)) {
      stop(
        compartment,
        ": missing input features in ",
        key,
        " n=",
        length(miss_f)
      )
    }

    miss_c <- setdiff(
      mm$original_cell,
      colnames(counts)
    )

    if(length(miss_c)) {
      stop(
        compartment,
        ": missing cells in ",
        key,
        " n=",
        length(miss_c)
      )
    }

    x <- counts[
      input_features,
      mm$original_cell,
      drop=FALSE
    ]

    colnames(x) <- mm$global_cell

    counts_list[[i]] <- x

    rm(obj, counts, x)
    invisible(gc())

    cat(
      "READ ",
      i, "/", nrow(lib_pairs),
      " ",
      pid, "::", key,
      "\n",
      sep=""
    )
  }

  counts_sub <- do.call(
    cbind,
    counts_list
  )

  rm(counts_list)
  invisible(gc())

  cells_order <- colnames(counts_sub)

  idx <- match(
    cells_order,
    meta$global_cell
  )

  if(anyNA(idx)) {
    stop(compartment, ": metadata order mismatch")
  }

  meta <- meta[idx]

  expressed <- Matrix::rowSums(counts_sub) > 0

  counts_sub <- counts_sub[
    expressed,
    ,
    drop=FALSE
  ]

  data_sub <- log_normalize_sparse(
    counts_sub
  )

  chosen_hvg <- intersect(
    consensus_hvg,
    rownames(data_sub)
  )

  ctl <- intersect(
    human_seg,
    rownames(data_sub)
  )

  if(length(chosen_hvg) < 2500L) {
    stop(compartment, ": too few retained HVGs")
  }

  if(length(ctl) < 100L) {
    stop(compartment, ": too few retained controls")
  }

  batch <- as.character(meta$project_id)
  condition <- as.character(meta$condition_binary)
  cell_types <- as.character(meta$integration_celltype)

  cat(
    "TRAINING CELLS=", length(cells_order),
    "; projects=", uniqueN(batch),
    "; libraries=", uniqueN(meta$library_key),
    "; cellTypes=", uniqueN(cell_types),
    "; HVG=", length(chosen_hvg),
    "; controls=", length(ctl),
    "\n",
    sep=""
  )

  set.seed(stable_seed[[compartment]])

  start <- proc.time()[["elapsed"]]

  result <- scMerge::scMerge2(
    exprsMat=data_sub,
    batch=batch,
    cellTypes=cell_types,
    condition=condition,
    ctl=ctl,
    chosen.hvg=chosen_hvg,
    ruvK=ruvK,
    use_bpparam=BiocParallel::SerialParam(),
    use_bsparam=BiocSingular::RandomParam(),
    use_bnparam=BiocNeighbors::AnnoyParam(),
    pseudoBulk_fn="create_pseudoBulk",
    k_pseudoBulk=kPB,
    k_celltype=10L,
    exprsMat_counts=counts_sub,
    cosineNorm=TRUE,
    return_subset=TRUE,
    return_subset_genes=chosen_hvg,
    return_matrix=TRUE,
    byChunk=TRUE,
    chunkSize=5000L,
    verbose=TRUE,
    seed=stable_seed[[compartment]]
  )

  elapsed <- proc.time()[["elapsed"]] - start

  adjusted <- result$newY

  if(is.null(adjusted)) {
    stop(compartment, ": scMerge2 returned no newY")
  }

  if(is.null(rownames(adjusted))) {
    rownames(adjusted) <- chosen_hvg
  }

  if(is.null(colnames(adjusted))) {
    colnames(adjusted) <- cells_order
  }

  if(
    !setequal(rownames(adjusted), chosen_hvg) ||
    !setequal(colnames(adjusted), cells_order)
  ) {
    stop(compartment, ": adjusted matrix row/column universe mismatch")
  }

  if(!identical(rownames(adjusted), chosen_hvg)) {
    adjusted <- adjusted[
      match(chosen_hvg, rownames(adjusted)),
      ,
      drop=FALSE
    ]
  }

  if(!identical(colnames(adjusted), cells_order)) {
    adjusted <- adjusted[
      ,
      match(cells_order, colnames(adjusted)),
      drop=FALSE
    ]
  }

  stopifnot(
    nrow(adjusted) == length(chosen_hvg),
    ncol(adjusted) == length(cells_order),
    identical(rownames(adjusted), chosen_hvg),
    identical(colnames(adjusted), cells_order),
    matrix_all_finite(adjusted)
  )

  matrix_file <- file.path(
    matrix_dir,
    paste0(
      compartment,
      "__scMerge2_adjusted_expression_v1.rds"
    )
  )

  saveRDS(
    adjusted,
    matrix_file,
    compress=FALSE
  )

  model_file <- file.path(
    model_dir,
    paste0(
      compartment,
      "__scMerge2_final_model_v1.rds"
    )
  )

  saveRDS(
    list(
      compartment=compartment,
      ruvK=ruvK,
      k_pseudoBulk=kPB,
      cells=cells_order,
      controls=ctl,
      chosen_hvg=chosen_hvg,
      fullalpha=result$fullalpha,
      M=result$M,
      scMerge_version=as.character(
        packageVersion("scMerge")
      ),
      correction_scope="within_compartment",
      batch_variable="project_id",
      condition_variable="condition_binary",
      celltype_variable="integration_celltype",
      seed=stable_seed[[compartment]],
      note=paste(
        "Final P05E4 reuse-reference model.",
        "Adjusted expression is retained separately.",
        "Not intended for differential-expression inference."
      )
    ),
    model_file,
    compress=FALSE
  )

  meta_out <- meta[, .(
    global_cell,
    original_cell,
    project_id,
    library_key,
    clinical_domain_std,
    condition_binary,
    final_compartment_v1,
    atlas_core_identity_v1,
    atlas_axis_v1,
    atlas_state_v1,
    integration_celltype,
    batch_condition_n,
    integration_celltype_n
  )]

  metadata_file <- file.path(
    metadata_dir,
    paste0(
      compartment,
      "__corrected_reference_cells_v1.tsv"
    )
  )

  fwrite(
    meta_out,
    metadata_file,
    sep="\t",
    quote=TRUE,
    na="NA"
  )

  # reread/order verification for plain TSV
  reread_meta <- fread(metadata_file)

  stopifnot(
    nrow(reread_meta) == nrow(meta_out),
    identical(
      reread_meta$global_cell,
      meta_out$global_cell
    ),
    identical(
      colnames(adjusted),
      meta_out$global_cell
    )
  )

  gene_sets[[compartment]] <- chosen_hvg

  manifest_rows[[length(manifest_rows)+1L]] <- data.table(
    compartment=compartment,
    ruvK=ruvK,
    k_pseudoBulk=kPB,
    n_cells=ncol(adjusted),
    n_genes=nrow(adjusted),
    n_projects=uniqueN(meta$project_id),
    n_libraries=uniqueN(meta$library_key),
    n_integration_celltypes=uniqueN(meta$integration_celltype),
    matrix_file=file.path(
      "matrices",
      basename(matrix_file)
    ),
    model_file=file.path(
      "models",
      basename(model_file)
    ),
    metadata_file=file.path(
      "metadata",
      basename(metadata_file)
    ),
    matrix_bytes=file.info(matrix_file)$size,
    elapsed_seconds=elapsed,
    correction_scope="within_compartment",
    intended_use="reference_mapping_annotation",
    DE_use="NO"
  )

  rm(
    result,
    adjusted,
    counts_sub,
    data_sub,
    meta,
    meta_out,
    reread_meta
  )

  invisible(gc())
}

# ============================================================
# Final manifests
# ============================================================

manifest <- rbindlist(manifest_rows)

fwrite(
  manifest,
  file.path(
    out,
    "corrected_reference_compartment_manifest.tsv"
  ),
  sep="\t"
)

common_genes <- consensus_hvg[
  consensus_hvg %in%
    Reduce(intersect, gene_sets)
]

if(length(common_genes) < 2500L) {
  stop(
    "Common corrected-reference gene set unexpectedly small: ",
    length(common_genes)
  )
}

common_gene_table <- data.table(
  common_gene_order=seq_along(common_genes),
  feature=common_genes
)

fwrite(
  common_gene_table,
  file.path(
    out,
    "corrected_reference_common_genes.tsv"
  ),
  sep="\t"
)

feature_manifest <- rbindlist(
  lapply(names(gene_sets), function(compartment) {
    g <- gene_sets[[compartment]]
    data.table(
      compartment=compartment,
      compartment_gene_order=seq_along(g),
      feature=g,
      in_common_gene_set=g %in% common_genes
    )
  })
)

fwrite(
  feature_manifest,
  file.path(
    out,
    "corrected_reference_feature_manifest.tsv"
  ),
  sep="\t"
)

if(length(all_excluded)) {

  excluded <- rbindlist(
    all_excluded,
    fill=TRUE
  )

  fwrite(
    excluded,
    file.path(
      out,
      "corrected_reference_excluded_cells.tsv"
    ),
    sep="\t",
    quote=TRUE
  )

  excluded_summary <- excluded[, .N, by=.(
    final_compartment_v1,
    exclusion_reason
  )]

} else {

  excluded_summary <- data.table(
    final_compartment_v1=character(),
    exclusion_reason=character(),
    N=integer()
  )
}

fwrite(
  excluded_summary,
  file.path(
    out,
    "corrected_reference_exclusion_summary.tsv"
  ),
  sep="\t"
)

# ============================================================
# Summary
# ============================================================

summary_lines <- c(
  "===== P05E4 FINAL CORRECTED REFERENCE EXPRESSION =====",
  "status=PASS",
  paste0(
    "corrected_compartments=",
    nrow(manifest)
  ),
  paste0(
    "corrected_reference_cells=",
    sum(manifest$n_cells)
  ),
  paste0(
    "common_reference_genes=",
    length(common_genes)
  ),
  "correction_scope=WITHIN_COMPARTMENT",
  "T_NK_ruvK=5",
  "Monocyte_DC_ruvK=2",
  "B_plasma_ruvK=5",
  "Neutrophil_ruvK=3",
  "Platelet_megakaryocyte_ruvK=3",
  "Platelet_megakaryocyte_status=CAUTIOUS",
  "Erythroid=NATIVE_ONLY",
  "Progenitor=NATIVE_ONLY",
  "Deferred_unresolved=EXCLUDED",
  "adjusted_expression_retained=YES",
  "corrected_reference_DE_use=NO",
  "native_processed_counts_DE_use=YES",
  "native_full_atlas_cells=665816",
  "annotation_modified=NO",
  "P05B_modified=NO",
  "cellranger_count=NOT_ACCESSED"
)

writeLines(
  summary_lines,
  file.path(
    out,
    "P05E4_SUMMARY.txt"
  )
)

# ============================================================
# SHA256 manifest
# ============================================================

oldwd <- getwd()
setwd(out)

files_for_sha <- sort(
  list.files(
    ".",
    recursive=TRUE,
    full.names=FALSE
  )
)

files_for_sha <- files_for_sha[
  files_for_sha != "SHA256SUMS_P05E4.txt"
]

sha_lines <- vapply(
  files_for_sha,
  function(f) {
    z <- system2(
      "sha256sum",
      args=f,
      stdout=TRUE,
      stderr=TRUE
    )
    status <- attr(z, "status")
    if(!is.null(status) && status != 0L) {
      stop("sha256sum failed for ", f)
    }
    z[[1]]
  },
  character(1)
)

writeLines(
  sha_lines,
  "SHA256SUMS_P05E4.txt"
)

check <- system2(
  "sha256sum",
  args=c("-c", "SHA256SUMS_P05E4.txt"),
  stdout=TRUE,
  stderr=TRUE
)

check_status <- attr(check, "status")

if(!is.null(check_status) && check_status != 0L) {
  cat(paste(check, collapse="\n"), "\n")
  stop("P05E4 SHA256 verification failed")
}

setwd(oldwd)

# ============================================================
# Root handoff copies
# ============================================================

handoff <- c(
  "P05E4_SUMMARY.txt",
  "corrected_reference_compartment_manifest.tsv",
  "corrected_reference_common_genes.tsv",
  "corrected_reference_exclusion_summary.tsv"
)

for(f in handoff) {

  src <- file.path(out, f)

  dst <- file.path(
    root,
    paste0(
      "P05E4_",
      sub("^P05E4_", "", f)
    )
  )

  file.copy(
    src,
    dst,
    overwrite=TRUE
  )
}

cat(
  readLines(
    file.path(out, "P05E4_SUMMARY.txt")
  ),
  sep="\n"
)

cat("\n\n===== COMPARTMENT MANIFEST =====\n")
print(manifest)

cat("\n===== EXCLUSION SUMMARY =====\n")
print(excluded_summary)

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
