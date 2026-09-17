#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3L){
  stop("Usage: script <root> <tag> <cluster_dir>")
}

root <- normalizePath(
  args[[1]],
  mustWork=TRUE
)

tag <- args[[2]]

cluster_dir <- normalizePath(
  args[[3]],
  mustWork=TRUE
)

suppressPackageStartupMessages({
  library(Matrix)
})

source(
  file.path(
    root,
    "full_atlas_primary_integration_v1",
    "00_common.R"
  )
)

# ============================================================
# Frozen design
# ============================================================

RESOLUTION <- "native_cluster_r0p4"

MIN_TARGET_CELLS_PER_LIBRARY <- 10L
MIN_BASELINE_CELLS_PER_LIBRARY <- 20L
MIN_ELIGIBLE_LIBRARIES <- 2L

TOP_N <- 40L

MARKER_MIN_MEDIAN_LOG2FC <- 0.50
MARKER_MIN_MEDIAN_DELTA_PCT <- 0.10
MARKER_MIN_POSITIVE_LIBRARY_FRACTION <- 0.60

# ============================================================
# Paths
# ============================================================

assignment_file <- file.path(
  cluster_dir,
  "neutrophil_project_native_cluster_assignments_v1.tsv.gz"
)

transfer_dir <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "annotation_transfer_by_library"
)

lib_file <- file.path(
  root,
  "pre_integration",
  "pilot_unintegrated_large_v1",
  "resolved_library_table.tsv"
)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "neutrophil_native_marker_aggregation_v1__",
    tag
  )
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

done_file <- file.path(
  out_dir,
  "NEUTROPHIL_NATIVE_MARKER_AGGREGATION_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop(
    "Already completed: ",
    done_file
  )
}

for(f in c(
  assignment_file,
  lib_file
)){
  if(!file.exists(f)){
    stop(
      "Missing input: ",
      f
    )
  }
}

# ============================================================
# Helpers
# ============================================================

read_tsv <- function(path){

  read.delim(
    path,
    sep="\t",
    quote="\"",
    comment.char="",
    stringsAsFactors=FALSE,
    check.names=FALSE
  )
}

write_gz_tsv <- function(
  x,
  path
){

  con <- gzfile(
    path,
    open="wt"
  )

  tryCatch(
    {
      write.table(
        x,
        file=con,
        sep="\t",
        quote=TRUE,
        qmethod="double",
        row.names=FALSE
      )
    },
    finally={
      close(con)
    }
  )

  status <- system2(
    "gzip",
    c(
      "-t",
      path
    ),
    stdout=FALSE,
    stderr=FALSE
  )

  if(status != 0L){
    stop(
      "gzip integrity check failed: ",
      path
    )
  }
}

safe_name <- function(x){

  gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    x
  )
}

is_technical_gene <- function(g){

  grepl(
    paste0(
      "^MT-|",
      "^RPL[0-9]|",
      "^RPS[0-9]|",
      "^HBA[12]$|",
      "^HBB$|",
      "^HBM$|",
      "^MALAT1$|",
      "^NEAT1$|",
      "^IG[HKL][VDJ]"
    ),
    g
  )
}

logcpm <- function(
  count,
  total
){

  log2(
    ((count + 0.5) /
       (total + 1)) *
      1e6
  )
}

row_median <- function(x){

  if(
    requireNamespace(
      "matrixStats",
      quietly=TRUE
    )
  ){
    matrixStats::rowMedians(
      x,
      na.rm=TRUE
    )
  } else {
    apply(
      x,
      1L,
      median,
      na.rm=TRUE
    )
  }
}

# ============================================================
# Native r0.4 assignments
# ============================================================

assign <- read_tsv(
  assignment_file
)

required_assignment <- c(
  "global_cell",
  "project_id",
  "library_key",
  "condition_binary",
  RESOLUTION
)

stopifnot(
  all(
    required_assignment %in%
      names(assign)
  ),
  nrow(assign) == 53460L,
  !anyDuplicated(
    assign$global_cell
  )
)

assign$global_cell <-
  as.character(
    assign$global_cell
  )

assign$project_id <-
  as.character(
    assign$project_id
  )

assign$library_key <-
  as.character(
    assign$library_key
  )

assign$native_cluster <-
  as.character(
    assign[[RESOLUTION]]
  )

eligible_projects <- sort(
  unique(
    assign$project_id
  )
)

stopifnot(
  length(eligible_projects) == 6L
)

# ============================================================
# Recover original_cell without touching cellranger_count
# ============================================================

transfer_files <- list.files(
  transfer_dir,
  pattern="\\.tsv\\.gz$",
  full.names=TRUE
)

stopifnot(
  length(transfer_files) == 158L
)

meta_list <- list()

for(i in seq_along(
  transfer_files
)){

  d <- read_tsv(
    transfer_files[[i]]
  )

  req <- c(
    "global_cell",
    "original_cell",
    "project_id",
    "library_key",
    "integration_compartment_primary_v1"
  )

  stopifnot(
    all(req %in% names(d))
  )

  d <- d[
    as.character(
      d$integration_compartment_primary_v1
    ) ==
      "Neutrophil_granulocyte",
    req,
    drop=FALSE
  ]

  if(nrow(d)){
    meta_list[[
      length(meta_list)+1L
    ]] <- d
  }
}

meta <- do.call(
  rbind,
  meta_list
)

rm(meta_list)

meta$global_cell <-
  as.character(meta$global_cell)

meta$original_cell <-
  as.character(meta$original_cell)

stopifnot(
  nrow(meta) == 53501L,
  !anyDuplicated(
    meta$global_cell
  )
)

mi <- match(
  assign$global_cell,
  meta$global_cell
)

stopifnot(
  !anyNA(mi)
)

assign$original_cell <-
  meta$original_cell[
    mi
  ]

stopifnot(
  assign$project_id ==
    as.character(
      meta$project_id[mi]
    ),
  assign$library_key ==
    as.character(
      meta$library_key[mi]
    )
)

rm(meta, mi)

# ============================================================
# Source library table
# ============================================================

lib <- read_tsv(
  lib_file
)

stopifnot(
  all(c(
    "project_id",
    "library_key",
    "final_rds_resolved"
  ) %in% names(lib))
)

lib$project_id <-
  as.character(
    lib$project_id
  )

lib$library_key <-
  as.character(
    lib$library_key
  )

lib$key <- paste(
  lib$project_id,
  lib$library_key,
  sep="|||"
)

assign$key <- paste(
  assign$project_id,
  assign$library_key,
  sep="|||"
)

# ============================================================
# Output accumulators
# ============================================================

all_top <- list()
all_cluster_summary <- list()
all_marker_support <- list()

feature_reference <- NULL

# ============================================================
# Project loop
# ============================================================

for(p in eligible_projects){

  cat(
    "\n============================================\n",
    "PROJECT: ",
    p,
    "\n",
    "============================================\n",
    sep=""
  )

  ap <- assign[
    assign$project_id == p,
    ,
    drop=FALSE
  ]

  project_clusters <- sort(
    unique(
      ap$native_cluster
    )
  )

  project_libraries <- sort(
    unique(
      ap$library_key
    )
  )

  cat(
    "cells=",
    nrow(ap),
    " libraries=",
    length(project_libraries),
    " clusters=",
    length(project_clusters),
    "\n",
    sep=""
  )

  # ----------------------------------------------------------
  # library x cluster native pseudobulk
  # ----------------------------------------------------------

  agg <- list()

  for(li in seq_along(
    project_libraries
  )){

    lk <- project_libraries[[li]]

    key <- paste(
      p,
      lk,
      sep="|||"
    )

    idx <- match(
      key,
      lib$key
    )

    if(is.na(idx)){
      stop(
        "Library missing from resolved table: ",
        key
      )
    }

    source_rds <- as.character(
      lib$final_rds_resolved[[idx]]
    )

    if(!file.exists(source_rds)){
      stop(
        "Missing source RDS: ",
        source_rds
      )
    }

    ml <- ap[
      ap$library_key == lk,
      ,
      drop=FALSE
    ]

    obj <- readRDS(
      source_rds
    )

    counts <- get_rna_counts(
      obj
    )

    if(is.null(
      feature_reference
    )){

      feature_reference <-
        rownames(counts)

    } else {

      if(
        !setequal(
          feature_reference,
          rownames(counts)
        )
      ){
        stop(
          "Feature-set mismatch: ",
          key
        )
      }

      counts <- counts[
        feature_reference,
        ,
        drop=FALSE
      ]
    }

    cells <- as.character(
      ml$original_cell
    )

    if(!all(
      cells %in%
        colnames(counts)
    )){
      stop(
        "Source cells missing for ",
        key
      )
    }

    x <- counts[
      ,
      cells,
      drop=FALSE
    ]

    # Force sparse numeric matrix for predictable aggregation.
    x <- as(
      x,
      "dgCMatrix"
    )

    cluster_factor <- factor(
      ml$native_cluster,
      levels=project_clusters
    )

    design <- Matrix::sparse.model.matrix(
      ~ 0 + cluster_factor
    )

    colnames(design) <-
      project_clusters

    # Gene UMI sum per native cluster.
    count_by_cluster <-
      x %*% design

    # Number of expressing cells per gene / cluster.
    detected <- x
    detected@x[] <- 1

    detect_by_cluster <-
      detected %*% design

    n_by_cluster <- as.integer(
      table(cluster_factor)
    )

    names(n_by_cluster) <-
      project_clusters

    umi_by_cluster <-
      Matrix::colSums(
        count_by_cluster
      )

    names(umi_by_cluster) <-
      project_clusters

    agg[[lk]] <- list(
      project_id=p,
      library_key=lk,
      clusters=project_clusters,
      n_by_cluster=n_by_cluster,
      umi_by_cluster=
        as.numeric(
          umi_by_cluster
        ),
      count_by_cluster=
        count_by_cluster,
      detect_by_cluster=
        detect_by_cluster,
      all_n=sum(
        n_by_cluster
      ),
      all_umi=sum(
        umi_by_cluster
      ),
      all_count=
        Matrix::rowSums(
          count_by_cluster
        ),
      all_detect=
        Matrix::rowSums(
          detect_by_cluster
        )
    )

    names(
      agg[[lk]]$umi_by_cluster
    ) <- project_clusters

    rm(
      obj,
      counts,
      x,
      detected,
      design,
      count_by_cluster,
      detect_by_cluster
    )

    invisible(
      gc()
    )

    cat(
      sprintf(
        "AGG %s %2d/%2d %s n=%d\n",
        p,
        li,
        length(project_libraries),
        lk,
        nrow(ml)
      )
    )
  }

  # ----------------------------------------------------------
  # Preserve project pseudobulk checkpoint
  # ----------------------------------------------------------

  project_agg_file <- file.path(
    out_dir,
    paste0(
      safe_name(p),
      "_native_library_cluster_pseudobulk_r0p4_v1.rds"
    )
  )

  saveRDS(
    list(
      project_id=p,
      resolution="r0.4",
      features=feature_reference,
      aggregates=agg
    ),
    project_agg_file,
    compress=FALSE
  )

  # ----------------------------------------------------------
  # Cluster-level composition diagnostics
  # ----------------------------------------------------------

  for(k in project_clusters){

    ak <- ap[
      ap$native_cluster == k,
      ,
      drop=FALSE
    ]

    lt <- sort(
      table(
        ak$library_key
      ),
      decreasing=TRUE
    )

    ct <- table(
      ak$condition_binary
    )

    all_cluster_summary[[
      length(all_cluster_summary)+1L
    ]] <- data.frame(
      project_id=p,
      cluster=k,
      n_cells=nrow(ak),
      n_libraries=
        length(lt),
      max_library=
        names(lt)[1],
      max_library_fraction=
        as.numeric(
          lt[1]
        ) /
        nrow(ak),
      healthy_fraction=
        if(
          "healthy" %in%
            names(ct)
        ){
          as.numeric(
            ct[["healthy"]]
          ) /
            nrow(ak)
        } else {
          NA_real_
        },
      disease_fraction=
        if(
          "disease" %in%
            names(ct)
        ){
          as.numeric(
            ct[["disease"]]
          ) /
            nrow(ak)
        } else {
          NA_real_
        },
      stringsAsFactors=FALSE
    )
  }

  # ----------------------------------------------------------
  # Native cluster vs rest, matched within library
  # ----------------------------------------------------------

  project_stats <- list()

  for(k in project_clusters){

    eligible_libraries <- character()

    for(lk in names(agg)){

      a <- agg[[lk]]

      target_n <-
        a$n_by_cluster[[k]]

      if(is.na(target_n)){
        target_n <- 0L
      }

      baseline_n <-
        a$all_n -
        target_n

      if(
        target_n >=
          MIN_TARGET_CELLS_PER_LIBRARY &&
        baseline_n >=
          MIN_BASELINE_CELLS_PER_LIBRARY
      ){
        eligible_libraries <- c(
          eligible_libraries,
          lk
        )
      }
    }

    n_eligible <-
      length(
        eligible_libraries
      )

    all_marker_support[[
      length(all_marker_support)+1L
    ]] <- data.frame(
      project_id=p,
      cluster=k,
      n_eligible_libraries=
        n_eligible,
      eligible_libraries=
        paste(
          eligible_libraries,
          collapse=","
        ),
      stringsAsFactors=FALSE
    )

    if(
      n_eligible <
        MIN_ELIGIBLE_LIBRARIES
    ){

      cat(
        "SKIP marker contrast ",
        p,
        " cluster ",
        k,
        ": eligible libraries=",
        n_eligible,
        "\n",
        sep=""
      )

      next
    }

    delta_logcpm <- matrix(
      NA_real_,
      nrow=length(
        feature_reference
      ),
      ncol=n_eligible,
      dimnames=list(
        feature_reference,
        eligible_libraries
      )
    )

    delta_pct <- delta_logcpm

    for(j in seq_along(
      eligible_libraries
    )){

      lk <- eligible_libraries[[j]]
      a <- agg[[lk]]

      target_n <-
        a$n_by_cluster[[k]]

      target_umi <-
        a$umi_by_cluster[[k]]

      target_count <-
        as.numeric(
          a$count_by_cluster[
            ,
            k
          ]
        )

      target_detect <-
        as.numeric(
          a$detect_by_cluster[
            ,
            k
          ]
        )

      baseline_n <-
        a$all_n -
        target_n

      baseline_umi <-
        a$all_umi -
        target_umi

      baseline_count <-
        a$all_count -
        target_count

      baseline_detect <-
        a$all_detect -
        target_detect

      delta_logcpm[,j] <-
        logcpm(
          target_count,
          target_umi
        ) -
        logcpm(
          baseline_count,
          baseline_umi
        )

      delta_pct[,j] <-
        (
          target_detect /
            target_n
        ) -
        (
          baseline_detect /
            baseline_n
        )
    }

    median_log2FC <-
      row_median(
        delta_logcpm
      )

    median_delta_pct <-
      row_median(
        delta_pct
      )

    positive_library_fraction <-
      rowMeans(
        delta_logcpm > 0,
        na.rm=TRUE
      )

    positive_delta_pct_fraction <-
      rowMeans(
        delta_pct > 0,
        na.rm=TRUE
      )

    strong_library_fraction <-
      rowMeans(
        delta_logcpm >=
          MARKER_MIN_MEDIAN_LOG2FC &
        delta_pct >=
          MARKER_MIN_MEDIAN_DELTA_PCT,
        na.rm=TRUE
      )

    stat <- data.frame(
      project_id=p,
      cluster=k,
      gene=feature_reference,
      n_eligible_libraries=
        n_eligible,
      median_log2FC=
        median_log2FC,
      median_delta_pct=
        median_delta_pct,
      positive_library_fraction=
        positive_library_fraction,
      positive_delta_pct_fraction=
        positive_delta_pct_fraction,
      strong_library_fraction=
        strong_library_fraction,
      technical_gene=
        is_technical_gene(
          feature_reference
        ),
      stringsAsFactors=FALSE
    )

    stat$marker_candidate <-
      !stat$technical_gene &
      stat$median_log2FC >=
        MARKER_MIN_MEDIAN_LOG2FC &
      stat$median_delta_pct >=
        MARKER_MIN_MEDIAN_DELTA_PCT &
      stat$positive_library_fraction >=
        MARKER_MIN_POSITIVE_LIBRARY_FRACTION

    stat$evidence_score <-
      pmax(
        stat$median_log2FC,
        0
      ) *
      pmax(
        stat$median_delta_pct,
        0
      ) *
      stat$positive_library_fraction *
      (
        0.5 +
        0.5 *
          stat$strong_library_fraction
      )

    project_stats[[
      length(project_stats)+1L
    ]] <- stat

    top <- stat[
      stat$marker_candidate,
      ,
      drop=FALSE
    ]

    top <- top[
      order(
        top$evidence_score,
        decreasing=TRUE
      ),
      ,
      drop=FALSE
    ]

    top <- head(
      top,
      TOP_N
    )

    if(nrow(top)){

      top$rank_in_cluster <-
        seq_len(
          nrow(top)
        )

      all_top[[
        length(all_top)+1L
      ]] <- top
    }

    cat(
      "MARKERS ",
      p,
      " cluster ",
      k,
      " eligible_libraries=",
      n_eligible,
      " candidates=",
      sum(
        stat$marker_candidate
      ),
      "\n",
      sep=""
    )

    rm(
      delta_logcpm,
      delta_pct,
      stat,
      top
    )

    invisible(
      gc()
    )
  }

  # ----------------------------------------------------------
  # Full gene-level marker statistics, project-specific
  # ----------------------------------------------------------

  if(length(project_stats)){

    ps <- do.call(
      rbind,
      project_stats
    )

    stats_file <- file.path(
      out_dir,
      paste0(
        safe_name(p),
        "_native_marker_stats_r0p4_v1.tsv.gz"
      )
    )

    write_gz_tsv(
      ps,
      stats_file
    )

    rm(ps)
  }

  rm(
    agg,
    project_stats
  )

  invisible(
    gc()
  )
}

# ============================================================
# Combined outputs
# ============================================================

cluster_summary <- do.call(
  rbind,
  all_cluster_summary
)

marker_support <- do.call(
  rbind,
  all_marker_support
)

if(length(all_top)){

  top_markers <- do.call(
    rbind,
    all_top
  )

} else {

  top_markers <- data.frame()
}

cluster_summary$key <- paste(
  cluster_summary$project_id,
  cluster_summary$cluster,
  sep="|||"
)

marker_support$key <- paste(
  marker_support$project_id,
  marker_support$cluster,
  sep="|||"
)

si <- match(
  cluster_summary$key,
  marker_support$key
)

cluster_summary$n_eligible_marker_libraries <-
  marker_support$n_eligible_libraries[
    si
  ]

cluster_summary$marker_support_status <-
  ifelse(
    cluster_summary$n_eligible_marker_libraries >=
      MIN_ELIGIBLE_LIBRARIES,
    "eligible",
    "insufficient_library_support"
  )

cluster_summary$key <- NULL

write.table(
  cluster_summary,
  file=file.path(
    out_dir,
    "neutrophil_native_cluster_summary_r0p4_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  marker_support[
    ,
    setdiff(
      names(marker_support),
      "key"
    ),
    drop=FALSE
  ],
  file=file.path(
    out_dir,
    "neutrophil_native_marker_library_support_r0p4_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  top_markers,
  file=file.path(
    out_dir,
    "neutrophil_native_top_markers_r0p4_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Compact top-marker strings for next-stage matching
# ============================================================

compact <- list()

if(nrow(top_markers)){

  keys <- unique(
    paste(
      top_markers$project_id,
      top_markers$cluster,
      sep="|||"
    )
  )

  for(key in keys){

    sp <- strsplit(
      key,
      "\\|\\|\\|"
    )[[1]]

    p <- sp[[1]]
    k <- sp[[2]]

    d <- top_markers[
      top_markers$project_id == p &
        top_markers$cluster == k,
      ,
      drop=FALSE
    ]

    d <- d[
      order(
        d$rank_in_cluster
      ),
      ,
      drop=FALSE
    ]

    compact[[
      length(compact)+1L
    ]] <- data.frame(
      project_id=p,
      cluster=k,
      n_top_markers=nrow(d),
      top_markers=
        paste(
          d$gene,
          collapse=","
        ),
      top_median_log2FC=
        paste(
          sprintf(
            "%.3f",
            d$median_log2FC
          ),
          collapse=","
        ),
      top_median_delta_pct=
        paste(
          sprintf(
            "%.3f",
            d$median_delta_pct
          ),
          collapse=","
        ),
      stringsAsFactors=FALSE
    )
  }
}

compact_df <- if(
  length(compact)
){
  do.call(
    rbind,
    compact
  )
} else {
  data.frame()
}

write.table(
  compact_df,
  file=file.path(
    out_dir,
    "neutrophil_native_cluster_marker_compact_r0p4_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Completion
# ============================================================

writeLines(
  c(
    "PASS",
    "Neutrophil native pseudobulk marker aggregation v1",
    "n_discovery_cells=53460",
    "n_projects=6",
    "primary_resolution=r0.4",
    "",
    "unit=library x native cluster pseudobulk",
    "contrast=native cluster versus rest within same library",
    "aggregation=median across eligible libraries within project",
    paste0(
      "minimum_target_cells_per_library=",
      MIN_TARGET_CELLS_PER_LIBRARY
    ),
    paste0(
      "minimum_baseline_cells_per_library=",
      MIN_BASELINE_CELLS_PER_LIBRARY
    ),
    paste0(
      "minimum_eligible_libraries=",
      MIN_ELIGIBLE_LIBRARIES
    ),
    "",
    "marker candidate thresholds:",
    paste0(
      "median_log2FC>=",
      MARKER_MIN_MEDIAN_LOG2FC
    ),
    paste0(
      "median_delta_pct>=",
      MARKER_MIN_MEDIAN_DELTA_PCT
    ),
    paste0(
      "positive_library_fraction>=",
      MARKER_MIN_POSITIVE_LIBRARY_FRACTION
    ),
    "",
    "no cross-project integration",
    "no cell-level pseudoreplication for marker evidence",
    "native SoupX-corrected RNA",
    "cellranger_count not accessed"
  ),
  done_file
)

cat(
  "\n===== CLUSTER SUMMARY =====\n"
)

print(
  cluster_summary,
  row.names=FALSE
)

cat(
  "\n===== TOP MARKERS =====\n"
)

if(nrow(top_markers)){

  print(
    top_markers[
      ,
      c(
        "project_id",
        "cluster",
        "gene",
        "n_eligible_libraries",
        "median_log2FC",
        "median_delta_pct",
        "positive_library_fraction",
        "strong_library_fraction",
        "evidence_score",
        "rank_in_cluster"
      ),
      drop=FALSE
    ],
    row.names=FALSE
  )

} else {

  cat(
    "No marker candidates passed thresholds.\n"
  )
}

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Neutrophil native marker aggregation completed\n"
)

