#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
})

base <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1"
)

cl_rds <- file.path(
  base,
  "integrated_compartment_clustering",
  "Monocyte_DC_combined",
  "integrated_clustering_v1.rds"
)

meta_file <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "scmerge2_primary_runs",
  "Monocyte_DC_combined__full_v1_cells.tsv.gz"
)

lib_file <- file.path(
  root,
  "pre_integration",
  "pilot_unintegrated_large_v1",
  "resolved_library_table.tsv"
)

prov_file <- file.path(
  base,
  "integrated_compartment_clustering",
  "Monocyte_DC_combined",
  "cluster_diagnostics_v1.tsv"
)

out_dir <- file.path(
  base,
  "monodc_project_balanced_marker_audit_v1"
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

for(f in c(cl_rds, meta_file, lib_file, prov_file)){
  if(!file.exists(f)){
    stop("Missing input: ", f)
  }
}

# ============================================================
# Clustering + cell metadata
# ============================================================

x <- readRDS(cl_rds)

cells <- as.character(x$cells)

stopifnot(
  "cluster_res_0p6" %in% names(x$clusters)
)

cluster <- as.character(
  x$clusters$cluster_res_0p6
)

stopifnot(
  length(cluster) == length(cells)
)

meta <- read.delim(
  gzfile(meta_file),
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

required <- c(
  "global_cell",
  "project_id",
  "library_key",
  "integration_celltype_full_v1",
  "condition_binary"
)

miss <- setdiff(
  required,
  names(meta)
)

if(length(miss)){
  stop(
    "Missing metadata columns: ",
    paste(miss, collapse=", ")
  )
}

m <- match(
  cells,
  as.character(meta$global_cell)
)

if(anyNA(m)){
  stop(
    "Cell metadata mismatch: ",
    sum(is.na(m))
  )
}

meta <- meta[m,,drop=FALSE]

stopifnot(
  identical(
    as.character(meta$global_cell),
    cells
  )
)

broad <- as.character(
  meta$integration_celltype_full_v1
)

allowed <- c(
  "Classical_monocyte_like",
  "Nonclassical_monocyte_like",
  "DC2_monocyte_boundary"
)

if(any(!broad %in% allowed)){
  bad <- sort(unique(
    broad[!broad %in% allowed]
  ))

  stop(
    "Unexpected Mono/DC backbone labels: ",
    paste(bad, collapse=", ")
  )
}

# ============================================================
# Determine dominant broad identity for each r0.6 cluster
# ============================================================

cluster_levels <- sort(
  unique(cluster)
)

cluster_identity <- do.call(
  rbind,
  lapply(
    cluster_levels,
    function(k){

      ii <- which(cluster == k)

      tt <- sort(
        table(broad[ii]),
        decreasing=TRUE
      )

      data.frame(
        cluster=k,
        n_cells=length(ii),
        dominant_backbone=names(tt)[1],
        dominant_backbone_fraction=
          as.numeric(tt[1]) / length(ii),
        stringsAsFactors=FALSE
      )
    }
  )
)

write.table(
  cluster_identity,
  file=file.path(
    out_dir,
    "cluster_backbone_identity_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# Library table
# ============================================================

lib <- read.delim(
  lib_file,
  sep="\t",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

req_lib <- c(
  "library_key",
  "project_id",
  "final_rds_resolved"
)

miss <- setdiff(
  req_lib,
  names(lib)
)

if(length(miss)){
  stop(
    "Missing library-table columns: ",
    paste(miss, collapse=", ")
  )
}

used_libs <- unique(
  as.character(meta$library_key)
)

li <- match(
  used_libs,
  as.character(lib$library_key)
)

if(anyNA(li)){
  stop(
    "Missing libraries: ",
    paste(
      used_libs[is.na(li)],
      collapse=", "
    )
  )
}

lib <- lib[li,,drop=FALSE]

# ============================================================
# RNA counts accessor
# ============================================================

get_rna_counts <- function(obj){

  if(!"RNA" %in% names(obj@assays)){
    stop("RNA assay missing")
  }

  z <- tryCatch(
    SeuratObject::LayerData(
      obj[["RNA"]],
      layer="counts"
    ),
    error=function(e) NULL
  )

  if(is.null(z)){
    z <- tryCatch(
      SeuratObject::GetAssayData(
        obj,
        assay="RNA",
        layer="counts"
      ),
      error=function(e) NULL
    )
  }

  if(is.null(z)){
    stop("Could not retrieve RNA count layer")
  }

  z
}

# ============================================================
# Feature reference
# ============================================================

first_obj <- readRDS(
  as.character(
    lib$final_rds_resolved[1]
  )
)

first_counts <- get_rna_counts(
  first_obj
)

features <- rownames(
  first_counts
)

nfeat <- length(features)

rm(
  first_obj,
  first_counts
)

gc(verbose=FALSE)

cat(
  "Cells=", length(cells),
  " Features=", nfeat,
  " Libraries=", nrow(lib),
  "\n",
  sep=""
)

# ============================================================
# Project × backbone × cluster aggregate
# ============================================================

projects <- sort(
  unique(
    as.character(meta$project_id)
  )
)

# We only need groups actually present.
group_info <- unique(
  data.frame(
    project=as.character(meta$project_id),
    backbone=broad,
    cluster=cluster,
    stringsAsFactors=FALSE
  )
)

group_info$key <- paste(
  group_info$project,
  group_info$backbone,
  group_info$cluster,
  sep="|||"
)

group_info <- group_info[
  order(group_info$key),
  ,
  drop=FALSE
]

group_keys <- group_info$key
ng <- nrow(group_info)

sum_counts <- matrix(
  0,
  nrow=nfeat,
  ncol=ng,
  dimnames=list(
    features,
    group_keys
  )
)

detect_counts <- matrix(
  0,
  nrow=nfeat,
  ncol=ng,
  dimnames=list(
    features,
    group_keys
  )
)

n_cells_group <- setNames(
  integer(ng),
  group_keys
)

cat("\n=== AGGREGATING NATIVE RNA ===\n")

for(i in seq_len(nrow(lib))){

  key <- as.character(
    lib$library_key[i]
  )

  f <- as.character(
    lib$final_rds_resolved[i]
  )

  ii <- which(
    as.character(meta$library_key) == key
  )

  if(!length(ii)){
    next
  }

  cat(
    sprintf(
      "[%3d/%3d] %s cells=%d\n",
      i,
      nrow(lib),
      key,
      length(ii)
    )
  )

  obj <- readRDS(f)

  counts <- get_rna_counts(
    obj
  )

  if(!identical(
    rownames(counts),
    features
  )){
    stop(
      "Feature order mismatch: ",
      key
    )
  }

  prefix <- paste0(
    key,
    "___"
  )

  global <- cells[ii]

  if(!all(startsWith(
    global,
    prefix
  ))){
    stop(
      "Cell-prefix mismatch: ",
      key
    )
  }

  original <- substring(
    global,
    nchar(prefix) + 1L
  )

  if(!all(
    original %in% colnames(counts)
  )){
    stop(
      "Cells absent from final RDS: ",
      key
    )
  }

  mat <- counts[
    ,
    original,
    drop=FALSE
  ]

  local_keys <- paste(
    as.character(meta$project_id[ii]),
    broad[ii],
    cluster[ii],
    sep="|||"
  )

  for(g in unique(local_keys)){

    jj <- which(
      local_keys == g
    )

    if(!g %in% group_keys){
      stop(
        "Unexpected aggregation group: ",
        g
      )
    }

    sum_counts[,g] <-
      sum_counts[,g] +
      as.numeric(
        Matrix::rowSums(
          mat[,jj,drop=FALSE]
        )
      )

    detect_counts[,g] <-
      detect_counts[,g] +
      as.numeric(
        Matrix::rowSums(
          mat[,jj,drop=FALSE] > 0
        )
      )

    n_cells_group[g] <-
      n_cells_group[g] +
      length(jj)
  }

  rm(
    obj,
    counts,
    mat
  )

  gc(verbose=FALSE)
}

stopifnot(
  sum(n_cells_group) ==
    length(cells)
)

# ============================================================
# Project-balanced cluster-vs-rest analysis
#
# Comparison:
#   same project
#   same transferred backbone
#   target r0.6 cluster
#       vs
#   all other cells in that project/backbone
# ============================================================

MIN_IN <- 50L
MIN_OUT <- 200L

per_project <- list()
cluster_summary <- list()
top_markers <- list()

technical_gene <- grepl(
  paste0(
    "^(",
    "MT-|",
    "RPL[0-9]|",
    "RPS[0-9]|",
    "HBA[12]$|",
    "HBB$|",
    "HBD$|",
    "HBG[12]$|",
    "MALAT1$|",
    "NEAT1$",
    ")"
  ),
  features
)

for(k in cluster_levels){

  ident <- cluster_identity[
    cluster_identity$cluster == k,
    "dominant_backbone"
  ][1]

  project_stats <- list()

  for(p in projects){

    target_key <- paste(
      p,
      ident,
      k,
      sep="|||"
    )

    if(!target_key %in% group_keys){
      next
    }

    same_pb <- which(
      group_info$project == p &
      group_info$backbone == ident
    )

    target_col <- match(
      target_key,
      group_keys
    )

    n_in <- n_cells_group[
      target_key
    ]

    n_total <- sum(
      n_cells_group[
        group_info$key[same_pb]
      ]
    )

    n_out <- n_total - n_in

    if(
      n_in < MIN_IN ||
      n_out < MIN_OUT
    ){
      next
    }

    cin <- sum_counts[
      ,
      target_col
    ]

    cout <- rowSums(
      sum_counts[
        ,
        same_pb,
        drop=FALSE
      ]
    ) - cin

    din <- detect_counts[
      ,
      target_col
    ]

    dout <- rowSums(
      detect_counts[
        ,
        same_pb,
        drop=FALSE
      ]
    ) - din

    lib_in <- sum(cin)
    lib_out <- sum(cout)

    cpm_in <-
      if(lib_in > 0)
        cin / lib_in * 1e6
      else
        rep(0, nfeat)

    cpm_out <-
      if(lib_out > 0)
        cout / lib_out * 1e6
      else
        rep(0, nfeat)

    log2fc <-
      log2(cpm_in + 0.1) -
      log2(cpm_out + 0.1)

    pct_in <- din / n_in
    pct_out <- dout / n_out

    delta_pct <-
      pct_in - pct_out

    project_stats[[
      length(project_stats) + 1L
    ]] <- data.frame(
      cluster=k,
      backbone=ident,
      project=p,
      n_in=n_in,
      n_out=n_out,
      gene=features,
      log2FC=log2fc,
      pct_in=pct_in,
      pct_out=pct_out,
      delta_pct=delta_pct,
      stringsAsFactors=FALSE
    )
  }

  nproj <- length(
    project_stats
  )

  if(!nproj){

    cluster_summary[[
      length(cluster_summary) + 1L
    ]] <- data.frame(
      cluster=k,
      backbone=ident,
      n_cells=
        cluster_identity$n_cells[
          cluster_identity$cluster == k
        ],
      dominant_backbone_fraction=
        cluster_identity$dominant_backbone_fraction[
          cluster_identity$cluster == k
        ],
      n_eligible_projects=0,
      n_reproducible_markers=0,
      stringsAsFactors=FALSE
    )

    next
  }

  z <- do.call(
    rbind,
    project_stats
  )

  per_project[[
    length(per_project) + 1L
  ]] <- z

  # gene × project matrices
  spl <- split(
    z,
    z$project
  )

  fc_mat <- do.call(
    cbind,
    lapply(
      spl,
      function(d){
        d$log2FC[
          match(features, d$gene)
        ]
      }
    )
  )

  dp_mat <- do.call(
    cbind,
    lapply(
      spl,
      function(d){
        d$delta_pct[
          match(features, d$gene)
        ]
      }
    )
  )

  if(is.null(dim(fc_mat))){
    fc_mat <- matrix(
      fc_mat,
      ncol=1
    )
  }

  if(is.null(dim(dp_mat))){
    dp_mat <- matrix(
      dp_mat,
      ncol=1
    )
  }

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
    fc_mat > 0.5 &
    dp_mat > 0.05,
    na.rm=TRUE
  )

  strong_fraction <- rowMeans(
    fc_mat > 1 &
    dp_mat > 0.10,
    na.rm=TRUE
  )

  # For 2 projects require both.
  # For >=3 projects require >=2/3 support.
  support_threshold <-
    if(nproj <= 2)
      1
    else
      2/3

  reproducible <-
    median_fc >= 0.5 &
    median_dp >= 0.05 &
    positive_fraction >=
      support_threshold

  d <- data.frame(
    cluster=k,
    backbone=ident,
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
    reproducible_marker=
      reproducible,
    stringsAsFactors=FALSE
  )

  d <- d[
    order(
      !d$reproducible_marker,
      -d$positive_project_fraction,
      -d$median_delta_pct,
      -d$median_log2FC
    ),
    ,
    drop=FALSE
  ]

  good <- d[
    d$reproducible_marker &
    !d$technical_gene,
    ,
    drop=FALSE
  ]

  top_markers[[
    length(top_markers) + 1L
  ]] <- head(
    good,
    100
  )

  cluster_summary[[
    length(cluster_summary) + 1L
  ]] <- data.frame(
    cluster=k,
    backbone=ident,
    n_cells=
      cluster_identity$n_cells[
        cluster_identity$cluster == k
      ],
    dominant_backbone_fraction=
      cluster_identity$dominant_backbone_fraction[
        cluster_identity$cluster == k
      ],
    n_eligible_projects=nproj,
    n_reproducible_markers=
      nrow(good),
    stringsAsFactors=FALSE
  )
}

cluster_summary <- do.call(
  rbind,
  cluster_summary
)

if(length(top_markers)){
  top_markers <- do.call(
    rbind,
    top_markers
  )
} else {
  top_markers <- data.frame()
}

if(length(per_project)){
  per_project <- do.call(
    rbind,
    per_project
  )
} else {
  per_project <- data.frame()
}

# ============================================================
# Add original project-dominance metrics
# ============================================================

diag <- read.delim(
  prov_file,
  sep="\t",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

diag <- diag[
  diag$resolution == 0.6,
  ,
  drop=FALSE
]

diag$cluster <- as.character(
  diag$cluster
)

cluster_summary$cluster <-
  as.character(
    cluster_summary$cluster
  )

keep <- c(
  "cluster",
  "n_projects",
  "max_project_fraction",
  "project_entropy_normalized",
  "transfer_celltype_purity",
  "healthy_fraction",
  "disease_fraction"
)

cluster_summary <- merge(
  cluster_summary,
  diag[,keep,drop=FALSE],
  by="cluster",
  all.x=TRUE,
  sort=FALSE
)

cluster_summary$cluster_num <-
  suppressWarnings(
    as.integer(
      cluster_summary$cluster
    )
  )

cluster_summary <- cluster_summary[
  order(
    cluster_summary$cluster_num
  ),
  ,
  drop=FALSE
]

cluster_summary$cluster_num <- NULL

# Classification is deliberately technical only.
cluster_summary$state_reproducibility_class <-
  ifelse(
    cluster_summary$n_eligible_projects >= 3 &
    cluster_summary$n_reproducible_markers >= 10,
    "cross_project_state_candidate",
    ifelse(
      cluster_summary$n_eligible_projects >= 2 &
      cluster_summary$n_reproducible_markers >= 5,
      "limited_cross_project_support",
      "project_specific_or_unresolved"
    )
  )

write.table(
  cluster_summary,
  file=file.path(
    out_dir,
    "monodc_r0p6_project_balanced_summary_v1.tsv"
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
    "monodc_r0p6_reproducible_top_markers_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# Per-project table can be large.
if(nrow(per_project)){

  con <- gzfile(
    file.path(
      out_dir,
      "monodc_r0p6_per_project_marker_stats_v1.tsv.gz"
    ),
    open="wt"
  )

  tryCatch(
    write.table(
      per_project,
      con,
      sep="\t",
      quote=TRUE,
      qmethod="double",
      row.names=FALSE
    ),
    finally=close(con)
  )
}

saveRDS(
  list(
    resolution=0.6,
    comparison=
      "cluster vs rest within same project and transferred backbone",
    minimum_cells_in=MIN_IN,
    minimum_cells_out=MIN_OUT,
    cluster_backbone=
      cluster_identity,
    cluster_summary=
      cluster_summary,
    top_reproducible_markers=
      top_markers
  ),
  file.path(
    out_dir,
    "monodc_project_balanced_marker_audit_v1.rds"
  ),
  compress=FALSE
)

writeLines(
  c(
    "PASS",
    "resolution=0.6",
    paste0(
      "n_cells=",
      length(cells)
    ),
    paste0(
      "n_clusters=",
      length(cluster_levels)
    ),
    paste0(
      "min_in=",
      MIN_IN
    ),
    paste0(
      "min_out=",
      MIN_OUT
    )
  ),
  file.path(
    out_dir,
    "PROJECT_BALANCED_MARKER_AUDIT_COMPLETE.ok"
  )
)

cat(
  "\n===== MONO/DC PROJECT-BALANCED SUMMARY =====\n"
)

print(
  cluster_summary[
    ,
    c(
      "cluster",
      "backbone",
      "n_cells",
      "max_project_fraction",
      "n_eligible_projects",
      "n_reproducible_markers",
      "state_reproducibility_class"
    )
  ],
  row.names=FALSE
)

cat(
  "\nPASS: Mono/DC project-balanced marker audit completed\n"
)
