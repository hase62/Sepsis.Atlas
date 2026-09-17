#!/usr/bin/env Rscript

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
if(length(missing))
  stop("Missing packages: ", paste(missing, collapse=", "))

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
  root, "publication", "scientific_data",
  paste0("reuse_reference_scmerge2_model_selection_v1__", tag)
)
dir.create(out, recursive=TRUE)

model_dir <- file.path(out, "candidate_models")
dir.create(model_dir)

# ============================================================
# Inputs
# ============================================================

p05b <- file.path(
  root, "publication", "scientific_data",
  "processed_data_release_candidate_v1__20260819_194718"
)

cells <- fread(file.path(
  p05b, "metadata",
  "SepsisAtlas_cell_metadata_v1__processed_data_release_candidate_v1__20260819_194718.tsv.gz"
))

libs <- fread(file.path(
  p05b, "metadata",
  "SepsisAtlas_library_metadata_with_source_accessions_v1__processed_data_release_candidate_v1__20260819_194718.tsv"
))

stopifnot(
  nrow(cells)==665816L,
  nrow(libs)==158L,
  !anyDuplicated(cells$global_cell)
)

# Latest P05E2b
p2dirs <- list.dirs(
  file.path(root, "publication", "scientific_data"),
  recursive=FALSE,
  full.names=TRUE
)

p2dirs <- sort(p2dirs[
  grepl(
    "reuse_reference_plan_refined_v1__[0-9]{8}_[0-9]{6}$",
    basename(p2dirs)
  )
])

stopifnot(length(p2dirs) >= 1L)
p2 <- tail(p2dirs, 1L)

plan <- fread(file.path(
  p2,
  "scmerge2_model_selection_plan.tsv"
))

stopifnot(
  nrow(plan)==3L,
  setequal(
    plan$final_compartment_v1,
    c("B_plasma","Platelet_megakaryocyte","Neutrophil")
  )
)

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

if(any(is.na(lib_cond$condition_binary)))
  stop("Unmapped clinical_domain_std")

# ============================================================
# final-QC RDS map
# ============================================================

rds_files <- Sys.glob(file.path(
  root, "pre_integration", "GSE*",
  "rds_final_qc", "*.rds"
))

stopifnot(length(rds_files)==158L)

rds_map <- rbindlist(lapply(rds_files, function(f) {
  data.table(
    project_id=basename(dirname(dirname(f))),
    library_key=sub(
      "__final_qc\\.rds$",
      "",
      basename(f)
    ),
    rds=f
  )
}))

stopifnot(
  nrow(rds_map)==158L,
  !anyDuplicated(
    paste(rds_map$project_id, rds_map$library_key)
  )
)

# ============================================================
# Feature space: same 3000 consensus HVG + human scSEG
# ============================================================

hvg_file <- file.path(
  root, "pre_integration",
  "pilot_unintegrated_large_v1",
  "large_pilot_consensus_hvg_3000.tsv"
)

stopifnot(file.exists(hvg_file))

hvg <- fread(hvg_file)
stopifnot("feature" %in% names(hvg))

consensus_hvg <- unique(as.character(hvg$feature))
stopifnot(length(consensus_hvg)==3000L)

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
# Metric functions: identical logic to primary integration
# ============================================================

run_pca <- function(mat, n_pcs=30L) {
  n_pcs <- min(
    n_pcs,
    nrow(mat)-1L,
    ncol(mat)-1L
  )
  fit <- irlba::prcomp_irlba(
    t(mat),
    n=n_pcs,
    center=TRUE,
    scale.=FALSE
  )
  v <- fit$sdev^2
  v <- v/sum(v)
  list(
    embeddings=fit$x,
    variance=v
  )
}

eta_squared <- function(score, group) {
  keep <- is.finite(score) & !is.na(group)
  score <- score[keep]
  group <- droplevels(factor(group[keep]))

  if(length(score)<3L || nlevels(group)<2L)
    return(NA_real_)

  total <- sum((score-mean(score))^2)
  if(total<=0) return(0)

  means <- tapply(score, group, mean)
  counts <- table(group)

  sum(
    counts * (means-mean(score))^2
  ) / total
}

weighted_eta <- function(
  pca,
  variance,
  group,
  n_pc=20L
) {
  n <- min(
    n_pc,
    ncol(pca),
    length(variance)
  )

  eta <- vapply(
    seq_len(n),
    function(i)
      eta_squared(pca[,i], group),
    numeric(1)
  )

  w <- variance[seq_len(n)]
  keep <- is.finite(eta) & is.finite(w)

  if(!any(keep))
    return(NA_real_)

  sum(eta[keep]*w[keep]) /
    sum(w[keep])
}

distance_preservation <- function(
  base_pca,
  adj_pca,
  max_cells=1500L
) {
  common <- intersect(
    rownames(base_pca),
    rownames(adj_pca)
  )

  if(length(common)<10L)
    return(NA_real_)

  set.seed(20260821)

  use <- if(length(common)<=max_cells)
    common
  else
    sample(common, max_cells)

  suppressWarnings(
    cor(
      as.numeric(
        dist(base_pca[use,,drop=FALSE])
      ),
      as.numeric(
        dist(adj_pca[use,,drop=FALSE])
      ),
      method="spearman"
    )
  )
}

# ============================================================
# Corrected-reference training universe
#
# Exclude explicitly unresolved/deferred labels.
# Native Atlas remains all 665,816 cells.
# ============================================================

is_reference_identity <- function(x) {
  ok <- !is.na(x) & nzchar(x)

  bad <- grepl(
    "unresolved|ambiguous|^Deferred_",
    x,
    ignore.case=TRUE
  )

  ok & !bad
}

all_metrics <- list()
training_support <- list()

# ============================================================
# Run three compartments
# ============================================================

for(ci in seq_len(nrow(plan))) {

  compartment <- as.character(
    plan$final_compartment_v1[[ci]]
  )

  cat(
    "\n========================================\n",
    "COMPARTMENT: ", compartment, "\n",
    "========================================\n",
    sep=""
  )

  meta <- cells[
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

  meta <- meta[
    is_reference_identity(
      atlas_core_identity_v1
    )
  ]

  meta <- merge(
    meta,
    lib_cond[, .(
      project_id,
      library_key,
      condition_binary
    )],
    by=c("project_id","library_key"),
    all.x=TRUE,
    sort=FALSE
  )

  if(any(is.na(meta$condition_binary)))
    stop(compartment, ": missing condition")

  # ----------------------------------------------------------
  # scMerge2 pseudo-bulk requires supported batch x condition
  # strata.  Very small strata are excluded from the corrected
  # reference fit; they remain present in the native Atlas.
  # ----------------------------------------------------------

  meta[, batch_condition_n := .N,
       by=.(project_id, condition_binary)]

  tiny_strata <- unique(
    meta[
      batch_condition_n < 10L,
      .(
        project_id,
        condition_binary,
        n_cells=batch_condition_n
      )
    ]
  )

  if(nrow(tiny_strata)) {
    cat("\nEXCLUDING TINY PROJECT x CONDITION STRATA (<10 cells):\n")
    print(tiny_strata)
  }

  meta <- meta[
    batch_condition_n >= 10L
  ]

  if(nrow(meta) < 100L)
    stop(compartment, ": too few supported training cells")

  # ----------------------------------------------------------
  # Biological grouping supplied to scMerge2.
  #
  # B/plasma: curated core identity.
  # Platelet: core + transcriptional axis + state.
  # Neutrophil: maturation axis + orthogonal state.
  # ----------------------------------------------------------

  clean_component <- function(x, missing_label) {
    x <- as.character(x)
    bad <- is.na(x) | !nzchar(trimws(x))
    x[bad] <- missing_label
    x
  }

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
    meta[, integration_celltype := core_lab]
  } else if(compartment == "Platelet_megakaryocyte") {
    meta[, integration_celltype :=
      paste(core_lab, axis_lab, state_lab, sep="|||")]
  } else if(compartment == "Neutrophil") {
    meta[, integration_celltype :=
      paste(axis_lab, state_lab, sep="|||")]
  } else {
    meta[, integration_celltype := core_lab]
  }

  cat(
    "BIOLOGICAL GROUPS=",
    uniqueN(meta$integration_celltype),
    "\n",
    sep=""
  )

  prefix <- paste0(
    meta$library_key,
    "___"
  )

  if(!all(startsWith(
    meta$global_cell,
    prefix
  )))
    stop(compartment, ": global_cell prefix mismatch")

  meta[, original_cell :=
    substring(
      global_cell,
      nchar(library_key)+4L
    )
  ]

  training_support[[length(training_support)+1L]] <-
    meta[, .N, by=.(
      final_compartment_v1,
      project_id,
      condition_binary,
      atlas_core_identity_v1
    )]

  # ----------------------------------------------------------
  # reconstruct native counts once
  # ----------------------------------------------------------

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
      project_id==pid &
      library_key==key
    ]

    rr <- rds_map[
      project_id==pid &
      library_key==key
    ]

    if(nrow(rr)!=1L)
      stop(
        "RDS mapping failure: ",
        pid, " / ", key
      )

    obj <- readRDS(rr$rds[[1]])
    counts <- get_rna_counts(obj)

    miss_f <- setdiff(
      input_features,
      rownames(counts)
    )

    if(length(miss_f))
      stop(
        compartment,
        ": missing features in ",
        key
      )

    miss_c <- setdiff(
      mm$original_cell,
      colnames(counts)
    )

    if(length(miss_c))
      stop(
        compartment,
        ": missing cells in ",
        key
      )

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

  if(anyNA(idx))
    stop(compartment, ": metadata order mismatch")

  meta <- meta[idx]

  expressed <- Matrix::rowSums(
    counts_sub
  ) > 0

  counts_sub <- counts_sub[
    expressed,,
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

  if(length(chosen_hvg)<2500L)
    stop(compartment, ": too few HVGs")

  if(length(ctl)<100L)
    stop(compartment, ": too few controls")

  batch <- as.character(
    meta$project_id
  )

  condition <- as.character(
    meta$condition_binary
  )

  cell_types <- as.character(
    meta$integration_celltype
  )

  cat(
    "TRAINING CELLS=", length(cells_order),
    "; projects=", length(unique(batch)),
    "; cellTypes=", length(unique(cell_types)),
    "; HVG=", length(chosen_hvg),
    "; controls=", length(ctl),
    "\n",
    sep=""
  )

  baseline_pca <- run_pca(
    data_sub[
      chosen_hvg,
      cells_order,
      drop=FALSE
    ],
    30L
  )

  rownames(
    baseline_pca$embeddings
  ) <- cells_order

  base_proj <- weighted_eta(
    baseline_pca$embeddings,
    baseline_pca$variance,
    batch
  )

  base_cond <- weighted_eta(
    baseline_pca$embeddings,
    baseline_pca$variance,
    condition
  )

  base_ct <- weighted_eta(
    baseline_pca$embeddings,
    baseline_pca$variance,
    cell_types
  )

  # ----------------------------------------------------------
  # ruvK = 3 and 5
  # ----------------------------------------------------------

  for(ruvK in c(3L,5L)) {

    cat(
      "\n--- ",
      compartment,
      " ruvK=", ruvK,
      " ---\n",
      sep=""
    )

    set.seed(
      20260821L +
      ci*100L +
      ruvK
    )

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
      k_pseudoBulk=5L,
      k_celltype=10L,
      exprsMat_counts=counts_sub,
      cosineNorm=TRUE,
      return_subset=TRUE,
      return_subset_genes=chosen_hvg,
      return_matrix=TRUE,
      byChunk=TRUE,
      chunkSize=5000L,
      verbose=TRUE,
      seed=20260821L + ci*100L + ruvK
    )

    elapsed <-
      proc.time()[["elapsed"]] -
      start

    adjusted <- result$newY

    if(is.null(adjusted))
      stop("scMerge2 returned no newY")

    if(is.null(rownames(adjusted)))
      rownames(adjusted) <- chosen_hvg

    if(is.null(colnames(adjusted)))
      colnames(adjusted) <- cells_order

    if(
      !identical(
        rownames(adjusted),
        chosen_hvg
      ) ||
      !identical(
        colnames(adjusted),
        cells_order
      )
    ) {
      adjusted <- adjusted[
        chosen_hvg,
        cells_order,
        drop=FALSE
      ]
    }

    pca_adj <- run_pca(
      adjusted,
      30L
    )

    rownames(
      pca_adj$embeddings
    ) <- cells_order

    adj_proj <- weighted_eta(
      pca_adj$embeddings,
      pca_adj$variance,
      batch
    )

    adj_cond <- weighted_eta(
      pca_adj$embeddings,
      pca_adj$variance,
      condition
    )

    adj_ct <- weighted_eta(
      pca_adj$embeddings,
      pca_adj$variance,
      cell_types
    )

    dist_cor <- distance_preservation(
      baseline_pca$embeddings,
      pca_adj$embeddings
    )

    met <- data.table(
      compartment=compartment,
      n_cells=length(cells_order),
      n_projects=uniqueN(batch),
      n_celltypes=uniqueN(cell_types),
      ruvK=ruvK,
      k_pseudoBulk=5L,

      baseline_project_eta2=base_proj,
      adjusted_project_eta2=adj_proj,

      project_eta2_reduction_fraction=
        1 - adj_proj/base_proj,

      baseline_condition_eta2=base_cond,
      adjusted_condition_eta2=adj_cond,

      condition_eta2_retention_ratio=
        ifelse(
          base_cond>0,
          adj_cond/base_cond,
          NA_real_
        ),

      baseline_celltype_eta2=base_ct,
      adjusted_celltype_eta2=adj_ct,

      celltype_eta2_retention_ratio=
        ifelse(
          base_ct>0,
          adj_ct/base_ct,
          NA_real_
        ),

      distance_spearman_native_vs_adjusted=
        dist_cor,

      elapsed_seconds=elapsed,
      rss_gb_after_run=rss_gb()
    )

    print(met)

    model_file <- file.path(
      model_dir,
      paste0(
        compartment,
        "__ruvK", ruvK,
        "__kPB5__candidate_v1.rds"
      )
    )

    saveRDS(
      list(
        compartment=compartment,
        ruvK=ruvK,
        k_pseudoBulk=5L,
        cells=cells_order,
        controls=ctl,
        chosen_hvg=chosen_hvg,
        fullalpha=result$fullalpha,
        M=result$M,
        pca_embeddings=pca_adj$embeddings,
        pca_variance=pca_adj$variance,
        metrics=met,
        scMerge_version=
          as.character(
            packageVersion("scMerge")
          ),
        note=paste(
          "P05E3 model-selection candidate.",
          "Adjusted matrix not retained at this stage."
        )
      ),
      model_file,
      compress=FALSE
    )

    all_metrics[[length(all_metrics)+1L]] <- met

    rm(
      adjusted,
      result,
      pca_adj
    )

    invisible(gc())
  }

  rm(
    counts_sub,
    data_sub,
    baseline_pca,
    meta
  )

  invisible(gc())
}

# ============================================================
# Outputs
# ============================================================

metrics <- rbindlist(all_metrics)

setorder(
  metrics,
  compartment,
  ruvK
)

fwrite(
  metrics,
  file.path(
    out,
    "scmerge2_model_selection_metrics.tsv"
  ),
  sep="\t"
)

support <- rbindlist(training_support)

fwrite(
  support,
  file.path(
    out,
    "scmerge2_training_cell_support.tsv"
  ),
  sep="\t"
)

manifest <- data.table(
  file=list.files(
    model_dir,
    pattern="\\.rds$",
    full.names=FALSE
  )
)

manifest[, bytes :=
  file.info(
    file.path(model_dir,file)
  )$size
]

fwrite(
  manifest,
  file.path(
    out,
    "candidate_model_manifest.tsv"
  ),
  sep="\t"
)

summary <- c(
  "===== P05E3 SCMERGE2 MODEL SELECTION =====",
  "status=PASS",
  paste0(
    "candidate_runs=",
    nrow(metrics)
  ),
  paste0(
    "compartments=",
    uniqueN(metrics$compartment)
  ),
  "ruvK_candidates=3,5",
  "k_pseudoBulk=5",
  "training_excludes_explicit_unresolved_deferred=YES",
  "native_full_atlas_cells=665816",
  "selection_frozen=NO",
  "expression_release_modified=NO",
  "annotation_modified=NO",
  "cellranger_count=NOT_ACCESSED"
)

writeLines(
  summary,
  file.path(out,"P05E3_SUMMARY.txt")
)

# Root handoff copies
for(f in c(
  "P05E3_SUMMARY.txt",
  "scmerge2_model_selection_metrics.tsv",
  "scmerge2_training_cell_support.tsv"
)) {

  file.copy(
    file.path(out,f),
    file.path(
      root,
      paste0(
        "P05E3_",
        sub("^P05E3_","",f)
      )
    ),
    overwrite=TRUE
  )
}

cat(
  readLines(
    file.path(out,"P05E3_SUMMARY.txt")
  ),
  sep="\n"
)

cat(
  "\n\n=== MODEL SELECTION METRICS ===\n"
)

print(metrics)

cat(
  "\nOUT_DIR=",
  out,
  "\n",
  sep=""
)
