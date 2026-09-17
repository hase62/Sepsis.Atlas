#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

in_rds <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_native_subclustering_v1_primary9997",
  "b_plasma_pilot_native_subclustering_v1.rds"
)

diag_file <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_native_subclustering_v1_primary9997",
  "b_plasma_cluster_diagnostics_all_resolutions_v1.tsv"
)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_project_balanced_marker_audit_r0p4_v1"
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

for(f in c(in_rds, diag_file)){
  if(!file.exists(f)){
    stop("Missing input: ", f)
  }
}

# ============================================================
# Load authoritative B/plasma discovery object
# ============================================================

bp <- readRDS(in_rds)

stopifnot(
  ncol(bp) == 9997,
  nrow(bp) == 38606
)

md <- bp[[]]

required_md <- c(
  "project_id",
  "library_key",
  "condition_binary",
  "cluster_res_0p4"
)

miss <- setdiff(
  required_md,
  names(md)
)

if(length(miss)){
  stop(
    "Missing metadata columns: ",
    paste(miss, collapse=", ")
  )
}

cluster <- as.character(
  md$cluster_res_0p4
)

project <- as.character(
  md$project_id
)

features <- rownames(bp)

cluster_levels <- sort(
  unique(cluster)
)

cluster_levels <- cluster_levels[
  order(
    suppressWarnings(
      as.integer(cluster_levels)
    )
  )
]

projects <- sort(
  unique(project)
)

cat(
  "cells     =", ncol(bp), "\n",
  "features  =", nrow(bp), "\n",
  "clusters  =", length(cluster_levels), "\n",
  "projects  =", length(projects), "\n"
)

# ============================================================
# Native SoupX-corrected counts
# ============================================================

counts <- tryCatch(
  SeuratObject::LayerData(
    bp[["RNA"]],
    layer="counts"
  ),
  error=function(e) NULL
)

if(is.null(counts)){
  counts <- SeuratObject::GetAssayData(
    bp,
    assay="RNA",
    layer="counts"
  )
}

stopifnot(
  nrow(counts) == 38606,
  ncol(counts) == 9997,
  identical(
    colnames(counts),
    rownames(md)
  )
)

# ============================================================
# Marker flags
# ============================================================

technical_gene <- grepl(
  paste0(
    "^MT-|",
    "^RPL[0-9]|",
    "^RPS[0-9]|",
    "^HBA[12]$|",
    "^HBB$|",
    "^HBD$|",
    "^HBG[12]$|",
    "^MALAT1$|",
    "^NEAT1$"
  ),
  features
)

# V(D)J/rearranged genes are retained in the statistics,
# but flagged because clonotype can dominate ranking.
rearranged_ig_gene <- grepl(
  paste0(
    "^IGHV|",
    "^IGKV|",
    "^IGLV|",
    "^IGHJ[0-9]|",
    "^IGKJ[0-9]|",
    "^IGLJ[0-9]"
  ),
  features
)

# Constant/isotype genes are NOT flagged:
# IGHM, IGHD, IGHG*, IGHA*, JCHAIN remain informative.

# ============================================================
# Project-balanced comparison
#
# For each r0.4 cluster:
#   within each project,
#   target cluster vs all other B/plasma cells.
#
# This prevents a project-dominated cluster from receiving
# a subtype label solely because of between-study differences.
# ============================================================

MIN_IN <- 30L
MIN_OUT <- 100L

all_gene_stats <- list()
top_marker_list <- list()
cluster_summary <- list()

for(k in cluster_levels){

  cat(
    "\n============================================================\n",
    "CLUSTER ", k, "\n",
    "============================================================\n",
    sep=""
  )

  per_project_fc <- list()
  per_project_dp <- list()
  eligible_projects <- character()
  project_sizes <- list()

  for(p in projects){

    in_idx <- which(
      project == p &
      cluster == k
    )

    out_idx <- which(
      project == p &
      cluster != k
    )

    n_in <- length(in_idx)
    n_out <- length(out_idx)

    if(
      n_in < MIN_IN ||
      n_out < MIN_OUT
    ){
      next
    }

    cin <- Matrix::rowSums(
      counts[,in_idx,drop=FALSE]
    )

    cout <- Matrix::rowSums(
      counts[,out_idx,drop=FALSE]
    )

    din <- Matrix::rowSums(
      counts[,in_idx,drop=FALSE] > 0
    )

    dout <- Matrix::rowSums(
      counts[,out_idx,drop=FALSE] > 0
    )

    lib_in <- sum(cin)
    lib_out <- sum(cout)

    if(
      lib_in <= 0 ||
      lib_out <= 0
    ){
      next
    }

    cpm_in <- cin / lib_in * 1e6
    cpm_out <- cout / lib_out * 1e6

    fc <- log2(cpm_in + 0.1) -
      log2(cpm_out + 0.1)

    pct_in <- din / n_in
    pct_out <- dout / n_out

    dp <- pct_in - pct_out

    eligible_projects <- c(
      eligible_projects,
      p
    )

    per_project_fc[[p]] <- fc
    per_project_dp[[p]] <- dp

    project_sizes[[p]] <- data.frame(
      project_id=p,
      n_in=n_in,
      n_out=n_out,
      stringsAsFactors=FALSE
    )

    cat(
      sprintf(
        "%-12s n_in=%4d n_out=%5d\n",
        p,
        n_in,
        n_out
      )
    )
  }

  nproj <- length(
    eligible_projects
  )

  total_cluster_cells <- sum(
    cluster == k
  )

  if(nproj == 0){

    cluster_summary[[
      length(cluster_summary)+1L
    ]] <- data.frame(
      cluster=k,
      n_cells=total_cluster_cells,
      n_eligible_projects=0,
      n_reproducible_markers=0,
      n_strong_reproducible_markers=0,
      marker_reproducibility_class=
        "project_specific_or_unresolved",
      stringsAsFactors=FALSE
    )

    next
  }

  fc_mat <- do.call(
    cbind,
    per_project_fc
  )

  dp_mat <- do.call(
    cbind,
    per_project_dp
  )

  rownames(fc_mat) <- features
  rownames(dp_mat) <- features

  median_fc <- apply(
    fc_mat,
    1,
    median,
    na.rm=TRUE
  )

  median_dp <- apply(
    dp_mat,
    1,
    median,
    na.rm=TRUE
  )

  positive_fraction <- rowMeans(
    fc_mat >= 0.5 &
    dp_mat >= 0.05,
    na.rm=TRUE
  )

  strong_fraction <- rowMeans(
    fc_mat >= 1 &
    dp_mat >= 0.10,
    na.rm=TRUE
  )

  # Two eligible projects:
  # both must support the marker.
  #
  # >=3 eligible projects:
  # at least 2/3 must support the marker.
  support_threshold <- if(
    nproj == 2
  ){
    1
  } else if(
    nproj >= 3
  ){
    2/3
  } else {
    1
  }

  reproducible <-
    nproj >= 2 &
    median_fc >= 0.5 &
    median_dp >= 0.05 &
    positive_fraction >=
      support_threshold

  strong_reproducible <-
    nproj >= 2 &
    median_fc >= 1 &
    median_dp >= 0.10 &
    strong_fraction >=
      support_threshold

  d <- data.frame(
    cluster=k,
    gene=features,
    n_eligible_projects=nproj,
    median_log2FC=median_fc,
    median_delta_pct=median_dp,
    positive_project_fraction=
      positive_fraction,
    strong_project_fraction=
      strong_fraction,
    technical_gene=
      technical_gene,
    rearranged_ig_gene=
      rearranged_ig_gene,
    reproducible_marker=
      reproducible,
    strong_reproducible_marker=
      strong_reproducible,
    stringsAsFactors=FALSE
  )

  all_gene_stats[[
    length(all_gene_stats)+1L
  ]] <- d

  good <- d[
    d$reproducible_marker &
    !d$technical_gene &
    !d$rearranged_ig_gene,
    ,
    drop=FALSE
  ]

  good <- good[
    order(
      -good$positive_project_fraction,
      -good$median_delta_pct,
      -good$median_log2FC
    ),
    ,
    drop=FALSE
  ]

  strong_good <- good[
    good$strong_reproducible_marker,
    ,
    drop=FALSE
  ]

  top_marker_list[[
    length(top_marker_list)+1L
  ]] <- head(
    good,
    100
  )

  cls <- if(
    nproj >= 3 &&
    nrow(good) >= 10
  ){
    "cross_project_subtype_candidate"
  } else if(
    nproj >= 2 &&
    nrow(good) >= 5
  ){
    "limited_cross_project_support"
  } else {
    "project_specific_or_unresolved"
  }

  cluster_summary[[
    length(cluster_summary)+1L
  ]] <- data.frame(
    cluster=k,
    n_cells=total_cluster_cells,
    n_eligible_projects=nproj,
    eligible_projects=paste(
      eligible_projects,
      collapse=","
    ),
    n_reproducible_markers=
      nrow(good),
    n_strong_reproducible_markers=
      nrow(strong_good),
    marker_reproducibility_class=
      cls,
    stringsAsFactors=FALSE
  )

  cat(
    "eligible projects       =", nproj, "\n",
    "reproducible markers    =", nrow(good), "\n",
    "strong reproducible     =", nrow(strong_good), "\n",
    "class                   =", cls, "\n"
  )
}

# ============================================================
# Combine
# ============================================================

gene_stats <- do.call(
  rbind,
  all_gene_stats
)

top_markers <- if(
  length(top_marker_list)
){
  do.call(
    rbind,
    top_marker_list
  )
} else {
  data.frame()
}

summary_df <- do.call(
  rbind,
  cluster_summary
)

# ============================================================
# Join existing clustering diagnostics
# ============================================================

diag <- read.delim(
  diag_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

diag <- diag[
  diag$resolution == 0.4,
  ,
  drop=FALSE
]

diag$cluster <- as.character(
  diag$cluster
)

summary_df$cluster <-
  as.character(
    summary_df$cluster
  )

keep_diag <- c(
  "cluster",
  "n_libraries",
  "n_projects",
  "top_project",
  "max_project_fraction",
  "project_entropy_normalized",
  "top_old_global_cluster",
  "old_global_cluster_purity",
  "healthy_fraction",
  "disease_fraction"
)

summary_df <- merge(
  summary_df,
  diag[,keep_diag,drop=FALSE],
  by="cluster",
  all.x=TRUE,
  sort=FALSE
)

summary_df$cluster_num <-
  suppressWarnings(
    as.integer(
      summary_df$cluster
    )
  )

summary_df <- summary_df[
  order(summary_df$cluster_num),
  ,
  drop=FALSE
]

summary_df$cluster_num <- NULL

# ============================================================
# Write outputs
# ============================================================

write.table(
  summary_df,
  file=file.path(
    out_dir,
    "b_plasma_r0p4_project_balanced_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

write.table(
  top_markers,
  file=file.path(
    out_dir,
    "b_plasma_r0p4_reproducible_top_markers_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

gene_file <- file.path(
  out_dir,
  "b_plasma_r0p4_project_balanced_gene_stats_v1.tsv.gz"
)

con <- gzfile(
  gene_file,
  "wt"
)

tryCatch(
  write.table(
    gene_stats,
    con,
    sep="\t",
    quote=1,
    qmethod="double",
    row.names=FALSE
  ),
  finally=close(con)
)

saveRDS(
  list(
    resolution=0.4,
    n_cells=9997,
    min_in=MIN_IN,
    min_out=MIN_OUT,
    comparison=
      "cluster vs all other B/plasma cells within the same project",
    summary=summary_df,
    top_reproducible_markers=
      top_markers
  ),
  file.path(
    out_dir,
    "b_plasma_project_balanced_marker_audit_r0p4_v1.rds"
  ),
  compress=FALSE
)

cat(
  "\n===== PROJECT-BALANCED SUMMARY =====\n"
)

print(
  summary_df[
    ,
    c(
      "cluster",
      "n_cells",
      "n_projects",
      "max_project_fraction",
      "n_eligible_projects",
      "n_reproducible_markers",
      "n_strong_reproducible_markers",
      "marker_reproducibility_class",
      "top_old_global_cluster",
      "healthy_fraction",
      "disease_fraction"
    )
  ],
  row.names=FALSE
)

writeLines(
  c(
    "PASS",
    "B/plasma project-balanced native marker audit",
    "resolution=0.4",
    "n_cells=9997",
    "comparison=cluster vs rest within project",
    "min_in=30",
    "min_out=100",
    "r0.4 used as discovery scaffold only"
  ),
  file.path(
    out_dir,
    "PROJECT_BALANCED_MARKER_AUDIT_COMPLETE.ok"
  )
)

cat(
  "\nPASS: B/plasma project-balanced marker audit completed\n"
)

