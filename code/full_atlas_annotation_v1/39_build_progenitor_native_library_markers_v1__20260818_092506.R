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

DISCOVERY_PROJECT <- "GSE216007"
N_DISCOVERY_EXPECTED <- 26962L

RESOLUTION <- "native_cluster_r0p4"
N_CLUSTERS_EXPECTED <- 14L
N_LIBRARIES_EXPECTED <- 5L

MIN_TARGET_CELLS_PER_LIBRARY <- 10L
MIN_BASELINE_CELLS_PER_LIBRARY <- 20L
MIN_ELIGIBLE_LIBRARIES_FORMAL <- 3L

TOP_N <- 40L

MARKER_MIN_MEDIAN_LOG2FC <- 0.50
MARKER_MIN_MEDIAN_DELTA_PCT <- 0.10
MARKER_MIN_POSITIVE_LIBRARY_FRACTION <- 0.60

# ============================================================
# Paths
# ============================================================

assignment_file <- file.path(
  cluster_dir,
  "progenitor_dominant_project_native_cluster_assignments_v1.tsv.gz"
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
    "progenitor_native_marker_aggregation_v1__",
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
  "PROGENITOR_NATIVE_MARKER_AGGREGATION_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop("Already completed: ", done_file)
}

# ============================================================
# Helpers
# ============================================================

read_tsv <- function(path){

  if(!file.exists(path)){
    stop("Missing file: ", path)
  }

  con <- if(grepl("\\.gz$", path)){
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }

  tryCatch(
    read.delim(
      con,
      sep="\t",
      quote="\"",
      comment.char="",
      stringsAsFactors=FALSE,
      check.names=FALSE
    ),
    finally=close(con)
  )
}

write_gz_tsv <- function(x, path){

  con <- gzfile(
    path,
    "wt"
  )

  tryCatch(
    write.table(
      x,
      con,
      sep="\t",
      quote=TRUE,
      qmethod="double",
      row.names=FALSE
    ),
    finally=close(con)
  )

  status <- system2(
    "gzip",
    c("-t", path),
    stdout=FALSE,
    stderr=FALSE
  )

  if(status != 0L){
    stop("gzip integrity failure: ", path)
  }
}

technical_gene_for_formal_marker <- function(g){

  grepl("^MT-", g) |
    grepl("^RPL[0-9]", g) |
    grepl("^RPS[0-9]", g) |
    g %in% c(
      "MALAT1",
      "NEAT1"
    )
}

safe_median <- function(x){

  if(!length(x)){
    return(NA_real_)
  }

  median(
    x,
    na.rm=TRUE
  )
}

# ============================================================
# Read discovery assignments
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
  )
)

assign$global_cell <- as.character(
  assign$global_cell
)

assign$project_id <- as.character(
  assign$project_id
)

assign$library_key <- as.character(
  assign$library_key
)

assign$native_cluster <- as.character(
  assign[[RESOLUTION]]
)

stopifnot(
  nrow(assign) == N_DISCOVERY_EXPECTED,
  !anyDuplicated(assign$global_cell),
  all(assign$project_id == DISCOVERY_PROJECT)
)

cluster_levels <- sort(
  unique(
    as.integer(
      assign$native_cluster
    )
  )
)

stopifnot(
  length(cluster_levels) ==
    N_CLUSTERS_EXPECTED,
  identical(
    cluster_levels,
    0:13
  )
)

cluster_levels <- as.character(
  cluster_levels
)

library_levels <- sort(
  unique(
    assign$library_key
  )
)

stopifnot(
  length(library_levels) ==
    N_LIBRARIES_EXPECTED
)

cat(
  "Discovery cells=",
  nrow(assign),
  "\n",
  sep=""
)

cat(
  "Libraries=",
  length(library_levels),
  "\n",
  sep=""
)

cat(
  "Clusters=",
  length(cluster_levels),
  "\n",
  sep=""
)

# ============================================================
# Cluster × library balance
# ============================================================

balance_tab <- table(
  cluster=assign$native_cluster,
  library_key=assign$library_key
)

cluster_balance <- do.call(
  rbind,
  lapply(
    cluster_levels,
    function(cl){

      x <- as.numeric(
        balance_tab[cl, ]
      )

      n <- sum(x)

      p <- x[x > 0] / n

      data.frame(
        cluster=cl,
        n_cells=n,
        n_libraries_present=sum(x > 0),
        max_library_fraction=max(p),
        effective_n_libraries=
          1 / sum(p^2),
        stringsAsFactors=FALSE
      )
    }
  )
)

cluster_balance$evidence_scope <- ifelse(
  cluster_balance$n_cells < 100L,
  "tiny_library_skewed_candidate",
  ifelse(
    cluster_balance$max_library_fraction >= 0.90,
    "library_confounded_candidate",
    "standard_native_discovery"
  )
)

stopifnot(
  setequal(
    cluster_balance$cluster[
      cluster_balance$evidence_scope ==
        "library_confounded_candidate"
    ],
    c("5","6","10","11")
  ),
  identical(
    cluster_balance$cluster[
      cluster_balance$evidence_scope ==
        "tiny_library_skewed_candidate"
    ],
    "13"
  )
)

# ============================================================
# Recover original_cell
# ============================================================

transfer_files <- sort(
  list.files(
    transfer_dir,
    pattern="\\.tsv\\.gz$",
    full.names=TRUE
  )
)

if(!length(transfer_files)){
  stop("No annotation-transfer files")
}

meta_list <- list()

for(f in transfer_files){

  d <- read_tsv(f)

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
    as.character(d$project_id) ==
      DISCOVERY_PROJECT &
      as.character(
        d$integration_compartment_primary_v1
      ) == "Progenitor",
    req,
    drop=FALSE
  ]

  if(nrow(d)){
    meta_list[[length(meta_list)+1L]] <- d
  }
}

meta <- do.call(
  rbind,
  meta_list
)

rm(meta_list)

for(nm in c(
  "global_cell",
  "original_cell",
  "project_id",
  "library_key"
)){
  meta[[nm]] <- as.character(
    meta[[nm]]
  )
}

stopifnot(
  nrow(meta) == N_DISCOVERY_EXPECTED,
  !anyDuplicated(meta$global_cell)
)

mi <- match(
  assign$global_cell,
  meta$global_cell
)

stopifnot(
  !anyNA(mi)
)

assign$original_cell <-
  meta$original_cell[mi]

stopifnot(
  assign$project_id ==
    meta$project_id[mi],
  assign$library_key ==
    meta$library_key[mi]
)

rm(meta, mi)
invisible(gc())

# ============================================================
# Resolved source libraries
# ============================================================

lib <- read_tsv(
  lib_file
)

stopifnot(
  all(
    c(
      "project_id",
      "library_key",
      "final_rds_resolved"
    ) %in% names(lib)
  )
)

lib$project_id <- as.character(
  lib$project_id
)

lib$library_key <- as.character(
  lib$library_key
)

lib$key <- paste(
  lib$project_id,
  lib$library_key,
  sep="|||"
)

# ============================================================
# Read native SoupX RNA once per discovery library
# ============================================================

counts_by_library <- list()

for(ii in seq_along(library_levels)){

  lk <- library_levels[[ii]]

  key <- paste(
    DISCOVERY_PROJECT,
    lk,
    sep="|||"
  )

  li <- match(
    key,
    lib$key
  )

  if(is.na(li)){
    stop(
      "Library lookup failed: ",
      key
    )
  }

  source_rds <- as.character(
    lib$final_rds_resolved[[li]]
  )

  if(!file.exists(source_rds)){
    stop(
      "Missing source RDS: ",
      source_rds
    )
  }

  d <- assign[
    assign$library_key == lk,
    ,
    drop=FALSE
  ]

  obj0 <- readRDS(
    source_rds
  )

  counts0 <- get_rna_counts(
    obj0
  )

  stopifnot(
    all(
      d$original_cell %in%
        colnames(counts0)
    )
  )

  m <- counts0[
    ,
    d$original_cell,
    drop=FALSE
  ]

  colnames(m) <- d$global_cell

  counts_by_library[[lk]] <- m

  cat(
    sprintf(
      "READ %d/%d %s cells=%d features=%d\n",
      ii,
      length(library_levels),
      lk,
      ncol(m),
      nrow(m)
    )
  )

  rm(
    obj0,
    counts0,
    m
  )

  invisible(gc())
}

# ============================================================
# Harmonize feature space
# ============================================================

feature_sets <- lapply(
  counts_by_library,
  rownames
)

common_features <- Reduce(
  intersect,
  feature_sets
)

first_features <- rownames(
  counts_by_library[[1]]
)

common_features <- first_features[
  first_features %in%
    common_features
]

if(length(common_features) < 10000L){
  stop(
    "Too few common features: ",
    length(common_features)
  )
}

counts_by_library <- lapply(
  counts_by_library,
  function(m){
    m[
      common_features,
      ,
      drop=FALSE
    ]
  }
)

cat(
  "Common native features=",
  length(common_features),
  "\n",
  sep=""
)

formal_gene_allowed <-
  !technical_gene_for_formal_marker(
    common_features
  )

# Globins MUST remain eligible.
globins_present <- intersect(
  c(
    "HBA1",
    "HBA2",
    "HBB",
    "HBD",
    "HBM"
  ),
  common_features
)

if(length(globins_present)){
  stopifnot(
    all(
      formal_gene_allowed[
        match(
          globins_present,
          common_features
        )
      ]
    )
  )
}

cat(
  "Globins retained for marker evidence=",
  paste(
    globins_present,
    collapse=","
  ),
  "\n",
  sep=""
)

# ============================================================
# Effect matrices
#
# One column per library for each cluster.
# ============================================================

ng <- length(common_features)
nl <- length(library_levels)

init_metric_list <- function(){

  setNames(
    lapply(
      cluster_levels,
      function(x){
        matrix(
          NA_real_,
          nrow=ng,
          ncol=nl,
          dimnames=list(
            common_features,
            library_levels
          )
        )
      }
    ),
    cluster_levels
  )
}

log2fc_mat <- init_metric_list()
delta_pct_mat <- init_metric_list()
target_pct_mat <- init_metric_list()
baseline_pct_mat <- init_metric_list()

library_support <- list()
library_top <- list()

# ============================================================
# Within-library cluster-vs-rest effects
# ============================================================

for(li in seq_along(library_levels)){

  lk <- library_levels[[li]]

  cat(
    "\n===== LIBRARY ",
    li,
    "/",
    nl,
    ": ",
    lk,
    " =====\n",
    sep=""
  )

  m <- counts_by_library[[lk]]

  ai <- match(
    colnames(m),
    assign$global_cell
  )

  stopifnot(
    !anyNA(ai)
  )

  cl <- assign$native_cluster[ai]

  cl_index <- match(
    cl,
    cluster_levels
  )

  stopifnot(
    !anyNA(cl_index)
  )

  design <- Matrix::sparseMatrix(
    i=seq_along(cl_index),
    j=cl_index,
    x=1,
    dims=c(
      length(cl_index),
      length(cluster_levels)
    ),
    dimnames=list(
      colnames(m),
      cluster_levels
    )
  )

  cluster_sum <- m %*% design
  cluster_detect <- (m > 0) %*% design

  total_gene_sum <- Matrix::rowSums(
    m
  )

  total_gene_detect <- Matrix::rowSums(
    m > 0
  )

  total_umi <- sum(
    total_gene_sum
  )

  cluster_n <- as.numeric(
    Matrix::colSums(
      design
    )
  )

  names(cluster_n) <- cluster_levels

  cluster_umi <- as.numeric(
    Matrix::colSums(
      cluster_sum
    )
  )

  names(cluster_umi) <- cluster_levels

  for(j in seq_along(cluster_levels)){

    clj <- cluster_levels[[j]]

    n_target <- cluster_n[[clj]]
    n_baseline <- ncol(m) - n_target

    eligible <-
      n_target >=
        MIN_TARGET_CELLS_PER_LIBRARY &&
      n_baseline >=
        MIN_BASELINE_CELLS_PER_LIBRARY

    library_support[[
      length(library_support)+1L
    ]] <- data.frame(
      cluster=clj,
      library_key=lk,
      n_target=n_target,
      n_baseline=n_baseline,
      eligible=eligible,
      stringsAsFactors=FALSE
    )

    if(!eligible){
      next
    }

    target_sum <- as.numeric(
      cluster_sum[, j]
    )

    target_detect <- as.numeric(
      cluster_detect[, j]
    )

    baseline_sum <-
      total_gene_sum -
      target_sum

    baseline_detect <-
      total_gene_detect -
      target_detect

    target_total_umi <-
      cluster_umi[[clj]]

    baseline_total_umi <-
      total_umi -
      target_total_umi

    if(
      target_total_umi <= 0 ||
      baseline_total_umi <= 0
    ){
      stop(
        "Invalid total UMI for ",
        lk,
        " cluster ",
        clj
      )
    }

    target_cpm <-
      target_sum /
      target_total_umi *
      1e6

    baseline_cpm <-
      baseline_sum /
      baseline_total_umi *
      1e6

    log2fc <- log2(
      (target_cpm + 1) /
        (baseline_cpm + 1)
    )

    target_pct <-
      target_detect /
      n_target

    baseline_pct <-
      baseline_detect /
      n_baseline

    delta_pct <-
      target_pct -
      baseline_pct

    log2fc_mat[[clj]][, li] <-
      log2fc

    delta_pct_mat[[clj]][, li] <-
      delta_pct

    target_pct_mat[[clj]][, li] <-
      target_pct

    baseline_pct_mat[[clj]][, li] <-
      baseline_pct

    # Descriptive per-library top markers.
    dd <- data.frame(
      cluster=clj,
      library_key=lk,
      gene=common_features,
      log2FC=log2fc,
      target_pct=target_pct,
      baseline_pct=baseline_pct,
      delta_pct=delta_pct,
      formal_gene_allowed=
        formal_gene_allowed,
      stringsAsFactors=FALSE
    )

    dd <- dd[
      dd$formal_gene_allowed &
        is.finite(dd$log2FC) &
        dd$log2FC > 0,
      ,
      drop=FALSE
    ]

    oo <- order(
      -dd$log2FC,
      -dd$delta_pct,
      dd$gene
    )

    dd <- dd[
      oo,
      ,
      drop=FALSE
    ]

    if(nrow(dd) > TOP_N){
      dd <- dd[
        seq_len(TOP_N),
        ,
        drop=FALSE
      ]
    }

    dd$rank_within_library <-
      seq_len(nrow(dd))

    library_top[[
      length(library_top)+1L
    ]] <- dd
  }

  rm(
    m,
    design,
    cluster_sum,
    cluster_detect,
    total_gene_sum,
    total_gene_detect
  )

  invisible(gc())
}

rm(counts_by_library)
invisible(gc())

library_support <- do.call(
  rbind,
  library_support
)

library_top <- do.call(
  rbind,
  library_top
)

rownames(library_support) <- NULL
rownames(library_top) <- NULL

# ============================================================
# Aggregate libraries with equal weight
# ============================================================

aggregate_list <- list()
top_list <- list()
cluster_summary <- list()

for(cl in cluster_levels){

  support_cl <- library_support[
    library_support$cluster == cl,
    ,
    drop=FALSE
  ]

  eligible_libraries <-
    support_cl$library_key[
      support_cl$eligible
    ]

  idx <- match(
    eligible_libraries,
    library_levels
  )

  idx <- idx[
    !is.na(idx)
  ]

  n_eligible <- length(idx)

  lfc <- log2fc_mat[[cl]][
    ,
    idx,
    drop=FALSE
  ]

  dp <- delta_pct_mat[[cl]][
    ,
    idx,
    drop=FALSE
  ]

  tp <- target_pct_mat[[cl]][
    ,
    idx,
    drop=FALSE
  ]

  bp <- baseline_pct_mat[[cl]][
    ,
    idx,
    drop=FALSE
  ]

  if(n_eligible > 0L){

    median_log2fc <- apply(
      lfc,
      1,
      median,
      na.rm=TRUE
    )

    mean_log2fc <- rowMeans(
      lfc,
      na.rm=TRUE
    )

    median_delta_pct <- apply(
      dp,
      1,
      median,
      na.rm=TRUE
    )

    median_target_pct <- apply(
      tp,
      1,
      median,
      na.rm=TRUE
    )

    median_baseline_pct <- apply(
      bp,
      1,
      median,
      na.rm=TRUE
    )

    positive_library_fraction <- rowMeans(
      (lfc > 0) &
        (dp > 0),
      na.rm=TRUE
    )

  } else {

    median_log2fc <- rep(
      NA_real_,
      ng
    )

    mean_log2fc <- rep(
      NA_real_,
      ng
    )

    median_delta_pct <- rep(
      NA_real_,
      ng
    )

    median_target_pct <- rep(
      NA_real_,
      ng
    )

    median_baseline_pct <- rep(
      NA_real_,
      ng
    )

    positive_library_fraction <- rep(
      NA_real_,
      ng
    )
  }

  scope <- cluster_balance$evidence_scope[
    match(
      cl,
      cluster_balance$cluster
    )
  ]

  marker_effect_pass <-
    n_eligible >=
      MIN_ELIGIBLE_LIBRARIES_FORMAL &
    formal_gene_allowed &
    median_log2fc >=
      MARKER_MIN_MEDIAN_LOG2FC &
    median_delta_pct >=
      MARKER_MIN_MEDIAN_DELTA_PCT &
    positive_library_fraction >=
      MARKER_MIN_POSITIVE_LIBRARY_FRACTION

  formal_marker <-
    marker_effect_pass &
    scope ==
      "standard_native_discovery"

  d <- data.frame(
    cluster=cl,
    gene=common_features,
    evidence_scope=scope,
    n_eligible_libraries=n_eligible,
    median_log2FC=median_log2fc,
    mean_log2FC=mean_log2fc,
    median_target_pct=median_target_pct,
    median_baseline_pct=median_baseline_pct,
    median_delta_pct=median_delta_pct,
    positive_library_fraction=
      positive_library_fraction,
    formal_gene_allowed=
      formal_gene_allowed,
    marker_effect_pass=
      marker_effect_pass,
    formal_marker=
      formal_marker,
    stringsAsFactors=FALSE
  )

  aggregate_list[[
    length(aggregate_list)+1L
  ]] <- d

  candidate <- d[
    d$formal_gene_allowed &
      is.finite(d$median_log2FC),
    ,
    drop=FALSE
  ]

  oo <- order(
    -as.integer(candidate$formal_marker),
    -as.integer(candidate$marker_effect_pass),
    -candidate$median_log2FC,
    -candidate$median_delta_pct,
    -candidate$positive_library_fraction,
    candidate$gene
  )

  candidate <- candidate[
    oo,
    ,
    drop=FALSE
  ]

  if(nrow(candidate) > TOP_N){
    candidate <- candidate[
      seq_len(TOP_N),
      ,
      drop=FALSE
    ]
  }

  candidate$rank <-
    seq_len(nrow(candidate))

  top_list[[
    length(top_list)+1L
  ]] <- candidate

  cluster_summary[[
    length(cluster_summary)+1L
  ]] <- data.frame(
    cluster=cl,
    n_cells=
      cluster_balance$n_cells[
        match(
          cl,
          cluster_balance$cluster
        )
      ],
    n_libraries_present=
      cluster_balance$n_libraries_present[
        match(
          cl,
          cluster_balance$cluster
        )
      ],
    max_library_fraction=
      cluster_balance$max_library_fraction[
        match(
          cl,
          cluster_balance$cluster
        )
      ],
    effective_n_libraries=
      cluster_balance$effective_n_libraries[
        match(
          cl,
          cluster_balance$cluster
        )
      ],
    n_eligible_marker_libraries=
      n_eligible,
    evidence_scope=
      scope,
    formal_marker_claim_allowed=
      scope ==
        "standard_native_discovery" &&
      n_eligible >=
        MIN_ELIGIBLE_LIBRARIES_FORMAL,
    n_marker_effect_pass=
      sum(
        marker_effect_pass,
        na.rm=TRUE
      ),
    n_formal_markers=
      sum(
        formal_marker,
        na.rm=TRUE
      ),
    stringsAsFactors=FALSE
  )
}

aggregate_df <- do.call(
  rbind,
  aggregate_list
)

top_df <- do.call(
  rbind,
  top_list
)

cluster_summary <- do.call(
  rbind,
  cluster_summary
)

rownames(aggregate_df) <- NULL
rownames(top_df) <- NULL
rownames(cluster_summary) <- NULL

formal_df <- aggregate_df[
  aggregate_df$formal_marker,
  ,
  drop=FALSE
]

# ============================================================
# Assertions on scope policy
# ============================================================

stopifnot(
  setequal(
    cluster_summary$cluster[
      cluster_summary$evidence_scope ==
        "library_confounded_candidate"
    ],
    c("5","6","10","11")
  ),
  identical(
    cluster_summary$cluster[
      cluster_summary$evidence_scope ==
        "tiny_library_skewed_candidate"
    ],
    "13"
  ),
  all(
    cluster_summary$n_formal_markers[
      cluster_summary$cluster %in%
        c("5","6","10","11","13")
    ] == 0L
  )
)

# ============================================================
# Write outputs
# ============================================================

aggregate_file <- file.path(
  out_dir,
  paste0(
    "progenitor_r0p4_marker_aggregate_v1__",
    tag,
    ".tsv.gz"
  )
)

write_gz_tsv(
  aggregate_df,
  aggregate_file
)

write.table(
  top_df,
  file.path(
    out_dir,
    "progenitor_r0p4_top_markers_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  formal_df,
  file.path(
    out_dir,
    "progenitor_r0p4_formal_markers_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  library_support,
  file.path(
    out_dir,
    "progenitor_r0p4_library_cluster_support_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  library_top,
  file.path(
    out_dir,
    "progenitor_r0p4_library_top_markers_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  cluster_summary,
  file.path(
    out_dir,
    "progenitor_r0p4_cluster_marker_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  cluster_balance,
  file.path(
    out_dir,
    "progenitor_r0p4_cluster_library_balance_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Small freeze object
# ============================================================

freeze <- list(
  marker_version=
    "Progenitor_native_marker_aggregation_v1",

  discovery_project=
    DISCOVERY_PROJECT,

  n_discovery_cells=
    N_DISCOVERY_EXPECTED,

  backbone_resolution=
    RESOLUTION,

  cluster_levels=
    cluster_levels,

  library_levels=
    library_levels,

  common_features=
    common_features,

  cluster_summary=
    cluster_summary,

  cluster_balance=
    cluster_balance,

  top_markers=
    top_df,

  formal_markers=
    formal_df,

  source_cluster_dir=
    cluster_dir,

  source_assignment_file=
    assignment_file
)

rds_file <- file.path(
  out_dir,
  paste0(
    "progenitor_native_marker_aggregation_v1__",
    tag,
    ".rds"
  )
)

tmp_rds <- paste0(
  rds_file,
  ".tmp"
)

saveRDS(
  freeze,
  tmp_rds,
  compress=TRUE
)

if(!file.rename(
  tmp_rds,
  rds_file
)){
  stop("Atomic RDS rename failed")
}

chk <- readRDS(
  rds_file
)

stopifnot(
  chk$n_discovery_cells ==
    N_DISCOVERY_EXPECTED,
  chk$backbone_resolution ==
    RESOLUTION,
  identical(
    chk$cluster_levels,
    cluster_levels
  )
)

# ============================================================
# gzip re-read
# ============================================================

aggregate_chk <- read_tsv(
  aggregate_file
)

stopifnot(
  nrow(aggregate_chk) ==
    nrow(aggregate_df),
  identical(
    as.character(
      aggregate_chk$cluster
    ),
    as.character(
      aggregate_df$cluster
    )
  ),
  identical(
    as.character(
      aggregate_chk$gene
    ),
    as.character(
      aggregate_df$gene
    )
  )
)

# ============================================================
# SHA256
# ============================================================

sha_files <- c(
  aggregate_file,
  file.path(
    out_dir,
    "progenitor_r0p4_top_markers_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_r0p4_formal_markers_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_r0p4_library_cluster_support_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_r0p4_library_top_markers_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_r0p4_cluster_marker_summary_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_r0p4_cluster_library_balance_v1.tsv"
  ),
  rds_file
)

sha_file <- file.path(
  out_dir,
  paste0(
    "SHA256SUMS_v1__",
    tag,
    ".txt"
  )
)

status <- system2(
  "sha256sum",
  sha_files,
  stdout=sha_file
)

if(status != 0L){
  stop("sha256sum creation failed")
}

status <- system2(
  "sha256sum",
  c("-c", sha_file),
  stdout=FALSE,
  stderr=FALSE
)

if(status != 0L){
  stop("sha256 verification failed")
}

# ============================================================
# Completion marker LAST
# ============================================================

writeLines(
  c(
    "PASS",
    "Progenitor native marker aggregation v1",
    "discovery_project=GSE216007",
    "n_discovery_cells=26962",
    "backbone_resolution=native_cluster_r0p4",
    "n_clusters=14",
    "n_discovery_libraries=5",
    "",
    "comparison=cluster-vs-rest within each library",
    "library aggregation=equal weight",
    "summary effect=median across eligible libraries",
    "minimum target cells/library=10",
    "minimum baseline cells/library=20",
    "minimum eligible libraries for formal marker=3",
    "formal median log2FC threshold=0.50",
    "formal median delta detection threshold=0.10",
    "formal positive-library fraction threshold=0.60",
    "",
    "c5,c6,c10,c11=library_confounded_candidate",
    "c13=tiny_library_skewed_candidate",
    "confounded/tiny clusters descriptive only",
    "",
    "globins retained for biological marker evidence",
    "MT/ribosomal/MALAT1/NEAT1 excluded only from formal markers",
    "condition not used",
    "non-discovery projects not used",
    "no cross-project integration",
    "native SoupX-corrected RNA",
    "cellranger_count not accessed",
    "",
    "gzip_integrity=PASS",
    "gzip_re_read=PASS",
    "freeze_rds_re_read=PASS",
    "sha256_manifest=PASS"
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

print(
  top_df[
    top_df$rank <= 12L,
    c(
      "cluster",
      "rank",
      "gene",
      "median_log2FC",
      "median_delta_pct",
      "positive_library_fraction",
      "formal_marker",
      "evidence_scope"
    ),
    drop=FALSE
  ],
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Progenitor native marker aggregation completed\n"
)
