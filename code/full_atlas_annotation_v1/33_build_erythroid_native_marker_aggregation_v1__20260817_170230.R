#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3L){
  stop("Usage: script <root> <tag> <discovery_dir>")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]
discovery_dir <- normalizePath(args[[3]], mustWork=TRUE)

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

MIN_TARGET_CELLS <- 5L
MIN_BASELINE_CELLS <- 10L
MIN_LIBRARIES_PER_PROJECT <- 2L
MIN_PROJECTS <- 3L
TOP_N <- 50L

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "erythroid_native_marker_aggregation_v1__",
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
  "ERYTHROID_NATIVE_MARKER_AGGREGATION_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop("Already completed: ", done_file)
}

assignment_file <- file.path(
  discovery_dir,
  "erythroid_balanced_discovery_cluster_assignments_v1.tsv.gz"
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

# ============================================================
# Helpers
# ============================================================

read_tsv <- function(path){

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

  con <- gzfile(path, "wt")

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

  if(
    system2(
      "gzip",
      c("-t", path),
      stdout=FALSE,
      stderr=FALSE
    ) != 0L
  ){
    stop("gzip integrity failure: ", path)
  }
}

logcpm <- function(count, total){

  log2(
    ((count + 0.5) /
       (total + 1)) *
      1e6
  )
}

row_median <- function(x){

  if(requireNamespace("matrixStats", quietly=TRUE)){
    matrixStats::rowMedians(x, na.rm=TRUE)
  } else {
    apply(x, 1L, median, na.rm=TRUE)
  }
}

technical_gene <- function(g){

  grepl(
    paste0(
      "^MT-|",
      "^RPL[0-9]|",
      "^RPS[0-9]|",
      "^MALAT1$|",
      "^NEAT1$|",
      "^IG[HKL]"
    ),
    g
  )
}

# ============================================================
# Discovery assignments
# ============================================================

a <- read_tsv(assignment_file)

stopifnot(
  nrow(a) == 8692L,
  all(c(
    "global_cell",
    "project_id",
    "library_key",
    "cluster_r0p4"
  ) %in% names(a)),
  !anyDuplicated(a$global_cell)
)

a$global_cell <- as.character(a$global_cell)
a$project_id <- as.character(a$project_id)
a$library_key <- as.character(a$library_key)
a$cluster <- as.character(a$cluster_r0p4)

clusters <- sort(
  unique(a$cluster)
)

projects <- sort(
  unique(a$project_id)
)

cat(
  "n_cells=", nrow(a),
  " n_projects=", length(projects),
  " n_clusters=", length(clusters),
  "\n",
  sep=""
)

# ============================================================
# Recover original cell IDs
# ============================================================

tf <- list.files(
  transfer_dir,
  pattern="\\.tsv\\.gz$",
  full.names=TRUE
)

meta_list <- list()

for(f in tf){

  z <- read_tsv(f)

  req <- c(
    "global_cell",
    "original_cell",
    "project_id",
    "library_key",
    "integration_compartment_primary_v1"
  )

  stopifnot(
    all(req %in% names(z))
  )

  z <- z[
    as.character(
      z$integration_compartment_primary_v1
    ) == "Erythroid",
    req,
    drop=FALSE
  ]

  if(nrow(z)){
    meta_list[[length(meta_list)+1L]] <- z
  }
}

meta <- do.call(
  rbind,
  meta_list
)

meta$global_cell <- as.character(meta$global_cell)
meta$original_cell <- as.character(meta$original_cell)

stopifnot(
  nrow(meta) == 12177L,
  !anyDuplicated(meta$global_cell)
)

mi <- match(
  a$global_cell,
  meta$global_cell
)

stopifnot(
  !anyNA(mi)
)

a$original_cell <- meta$original_cell[mi]

stopifnot(
  a$project_id ==
    as.character(meta$project_id[mi]),
  a$library_key ==
    as.character(meta$library_key[mi])
)

rm(meta, meta_list, mi)
gc()

# ============================================================
# Library source lookup
# ============================================================

lib <- read_tsv(lib_file)

stopifnot(
  all(c(
    "project_id",
    "library_key",
    "final_rds_resolved"
  ) %in% names(lib))
)

lib$project_id <- as.character(lib$project_id)
lib$library_key <- as.character(lib$library_key)

lib$key <- paste(
  lib$project_id,
  lib$library_key,
  sep="|||"
)

# ============================================================
# Native library × cluster pseudobulk
# ============================================================

feature_reference <- NULL
library_effects <- list()
composition <- list()

library_keys <- unique(
  paste(
    a$project_id,
    a$library_key,
    sep="|||"
  )
)

for(n in seq_along(library_keys)){

  key <- library_keys[[n]]

  sp <- strsplit(
    key,
    "\\|\\|\\|"
  )[[1]]

  p <- sp[[1]]
  lk <- sp[[2]]

  md <- a[
    a$project_id == p &
      a$library_key == lk,
    ,
    drop=FALSE
  ]

  li <- match(
    key,
    lib$key
  )

  if(is.na(li)){
    stop("Missing resolved library: ", key)
  }

  rds <- as.character(
    lib$final_rds_resolved[[li]]
  )

  if(!file.exists(rds)){
    stop("Missing source RDS: ", rds)
  }

  obj <- readRDS(rds)
  counts <- get_rna_counts(obj)

  if(is.null(feature_reference)){

    feature_reference <- rownames(counts)

  } else {

    if(!setequal(
      feature_reference,
      rownames(counts)
    )){
      stop("Feature mismatch: ", key)
    }

    counts <- counts[
      feature_reference,
      ,
      drop=FALSE
    ]
  }

  cells <- md$original_cell

  if(!all(
    cells %in%
      colnames(counts)
  )){
    stop("Missing source cells: ", key)
  }

  x <- as(
    counts[, cells, drop=FALSE],
    "dgCMatrix"
  )

  stopifnot(
    identical(
      colnames(x),
      cells
    )
  )

  cf <- factor(
    md$cluster,
    levels=clusters
  )

  design <- Matrix::sparse.model.matrix(
    ~ 0 + cf
  )

  colnames(design) <- clusters

  sums <- x %*% design

  detected <- x
  detected@x[] <- 1

  detect_sums <- detected %*% design

  n_by_cluster <- as.integer(
    table(cf)
  )

  names(n_by_cluster) <- clusters

  umi_by_cluster <- as.numeric(
    Matrix::colSums(sums)
  )

  names(umi_by_cluster) <- clusters

  all_count <- Matrix::rowSums(sums)
  all_detect <- Matrix::rowSums(detect_sums)

  all_n <- ncol(x)
  all_umi <- sum(umi_by_cluster)

  for(k in clusters){

    target_n <- n_by_cluster[[k]]

    if(is.na(target_n)){
      target_n <- 0L
    }

    baseline_n <- all_n - target_n

    composition[[length(composition)+1L]] <- data.frame(
      project_id=p,
      library_key=lk,
      cluster=k,
      target_cells=target_n,
      library_cells=all_n,
      target_fraction=
        if(all_n > 0)
          target_n / all_n
        else
          NA_real_,
      eligible=
        target_n >= MIN_TARGET_CELLS &&
        baseline_n >= MIN_BASELINE_CELLS,
      stringsAsFactors=FALSE
    )

    if(
      target_n < MIN_TARGET_CELLS ||
      baseline_n < MIN_BASELINE_CELLS
    ){
      next
    }

    target_count <- as.numeric(
      sums[, k]
    )

    target_detect <- as.numeric(
      detect_sums[, k]
    )

    target_umi <- umi_by_cluster[[k]]

    baseline_count <-
      all_count -
      target_count

    baseline_detect <-
      all_detect -
      target_detect

    baseline_umi <-
      all_umi -
      target_umi

    delta_logcpm <-
      logcpm(
        target_count,
        target_umi
      ) -
      logcpm(
        baseline_count,
        baseline_umi
      )

    delta_pct <-
      target_detect /
        target_n -
      baseline_detect /
        baseline_n

    library_effects[[length(library_effects)+1L]] <-
      data.frame(
        project_id=p,
        library_key=lk,
        cluster=k,
        gene=feature_reference,
        target_cells=target_n,
        baseline_cells=baseline_n,
        delta_logCPM=delta_logcpm,
        delta_pct=delta_pct,
        stringsAsFactors=FALSE
      )
  }

  rm(
    obj,
    counts,
    x,
    detected,
    design,
    sums,
    detect_sums
  )

  invisible(gc())

  cat(
    sprintf(
      "LIB %3d/%3d %s\n",
      n,
      length(library_keys),
      key
    )
  )
}

composition <- do.call(
  rbind,
  composition
)

effects <- do.call(
  rbind,
  library_effects
)

rm(library_effects)
gc()

# ============================================================
# First aggregate libraries within each project
# ============================================================

project_stats <- list()

group_key <- paste(
  effects$project_id,
  effects$cluster,
  effects$gene,
  sep="|||"
)

groups <- split(
  seq_len(nrow(effects)),
  group_key
)

for(ii in groups){

  d <- effects[ii, , drop=FALSE]

  project_stats[[length(project_stats)+1L]] <-
    data.frame(
      project_id=d$project_id[[1]],
      cluster=d$cluster[[1]],
      gene=d$gene[[1]],
      n_libraries=nrow(d),
      median_delta_logCPM=
        median(
          d$delta_logCPM,
          na.rm=TRUE
        ),
      median_delta_pct=
        median(
          d$delta_pct,
          na.rm=TRUE
        ),
      positive_library_fraction=
        mean(
          d$delta_logCPM > 0,
          na.rm=TRUE
        ),
      stringsAsFactors=FALSE
    )
}

project_stats <- do.call(
  rbind,
  project_stats
)

# ============================================================
# Then aggregate projects with equal project weight
# ============================================================

eligible_project_stats <- project_stats[
  project_stats$n_libraries >=
    MIN_LIBRARIES_PER_PROJECT,
  ,
  drop=FALSE
]

consensus <- list()

group_key <- paste(
  eligible_project_stats$cluster,
  eligible_project_stats$gene,
  sep="|||"
)

groups <- split(
  seq_len(nrow(eligible_project_stats)),
  group_key
)

for(ii in groups){

  d <- eligible_project_stats[
    ii,
    ,
    drop=FALSE
  ]

  consensus[[length(consensus)+1L]] <-
    data.frame(
      cluster=d$cluster[[1]],
      gene=d$gene[[1]],
      n_projects=nrow(d),
      projects=
        paste(
          sort(
            unique(
              d$project_id
            )
          ),
          collapse=","
        ),
      median_project_delta_logCPM=
        median(
          d$median_delta_logCPM,
          na.rm=TRUE
        ),
      median_project_delta_pct=
        median(
          d$median_delta_pct,
          na.rm=TRUE
        ),
      positive_project_fraction=
        mean(
          d$median_delta_logCPM > 0,
          na.rm=TRUE
        ),
      strong_project_fraction=
        mean(
          d$median_delta_logCPM >= 0.5 &
            d$median_delta_pct >= 0.05,
          na.rm=TRUE
        ),
      stringsAsFactors=FALSE
    )
}

consensus <- do.call(
  rbind,
  consensus
)

consensus$technical_gene <-
  technical_gene(
    consensus$gene
  )

consensus$marker_candidate <-
  !consensus$technical_gene &
  consensus$n_projects >=
    MIN_PROJECTS &
  consensus$median_project_delta_logCPM >=
    0.5 &
  consensus$median_project_delta_pct >=
    0.05 &
  consensus$positive_project_fraction >=
    0.75

consensus$evidence_score <-
  pmax(
    consensus$median_project_delta_logCPM,
    0
  ) *
  pmax(
    consensus$median_project_delta_pct,
    0
  ) *
  consensus$positive_project_fraction *
  log1p(
    consensus$n_projects
  )

# ============================================================
# Top markers
# ============================================================

top_list <- list()

for(k in clusters){

  d <- consensus[
    consensus$cluster == k &
      consensus$marker_candidate,
    ,
    drop=FALSE
  ]

  d <- d[
    order(
      d$evidence_score,
      decreasing=TRUE
    ),
    ,
    drop=FALSE
  ]

  d <- head(
    d,
    TOP_N
  )

  if(nrow(d)){

    d$rank_in_cluster <-
      seq_len(
        nrow(d)
      )

    top_list[[length(top_list)+1L]] <- d
  }
}

top <- if(length(top_list)){
  do.call(rbind, top_list)
} else {
  data.frame()
}

# ============================================================
# Cluster support summary
# ============================================================

support <- do.call(
  rbind,
  lapply(
    clusters,
    function(k){

      c0 <- composition[
        composition$cluster == k,
        ,
        drop=FALSE
      ]

      c1 <- c0[
        c0$eligible,
        ,
        drop=FALSE
      ]

      projects_supported <- unique(
        c1$project_id
      )

      data.frame(
        cluster=k,

        n_discovery_cells=
          sum(
            a$cluster == k
          ),

        n_projects_present=
          length(
            unique(
              a$project_id[
                a$cluster == k
              ]
            )
          ),

        n_libraries_present=
          length(
            unique(
              a$library_key[
                a$cluster == k
              ]
            )
          ),

        n_eligible_libraries=
          nrow(c1),

        n_projects_with_eligible_libraries=
          length(
            projects_supported
          ),

        projects_with_eligible_libraries=
          paste(
            sort(
              projects_supported
            ),
            collapse=","
          ),

        stringsAsFactors=FALSE
      )
    }
  )
)

# ============================================================
# Compact top marker table
# ============================================================

compact <- list()

if(nrow(top)){

  for(k in clusters){

    d <- top[
      top$cluster == k,
      ,
      drop=FALSE
    ]

    compact[[length(compact)+1L]] <-
      data.frame(
        cluster=k,
        n_top_markers=nrow(d),
        top_markers=
          paste(
            d$gene,
            collapse=","
          ),
        top_effects=
          paste(
            sprintf(
              "%.3f",
              d$median_project_delta_logCPM
            ),
            collapse=","
          ),
        top_positive_project_fraction=
          paste(
            sprintf(
              "%.2f",
              d$positive_project_fraction
            ),
            collapse=","
          ),
        stringsAsFactors=FALSE
      )
  }
}

compact <- do.call(
  rbind,
  compact
)

# ============================================================
# Outputs
# ============================================================

write.table(
  support,
  file=file.path(
    out_dir,
    "erythroid_r0p4_native_marker_support_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  project_stats,
  file=file.path(
    out_dir,
    "erythroid_r0p4_project_native_marker_stats_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write_gz_tsv(
  consensus,
  file.path(
    out_dir,
    "erythroid_r0p4_project_balanced_marker_stats_v1.tsv.gz"
  )
)

write.table(
  top,
  file=file.path(
    out_dir,
    "erythroid_r0p4_project_balanced_top_markers_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  compact,
  file=file.path(
    out_dir,
    "erythroid_r0p4_project_balanced_top_markers_compact_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Re-read / integrity
# ============================================================

chk <- read_tsv(
  file.path(
    out_dir,
    "erythroid_r0p4_project_balanced_marker_stats_v1.tsv.gz"
  )
)

stopifnot(
  nrow(chk) ==
    nrow(consensus)
)

rm(chk)

# ============================================================
# Completion
# ============================================================

writeLines(
  c(
    "PASS",
    "Erythroid native project-balanced marker aggregation v1",
    "discovery_cells=8692",
    "backbone_resolution=r0.4",
    "",
    "native SoupX-corrected RNA",
    "unit=library x integrated-discovery-cluster",
    "contrast=cluster versus rest within same library",
    "library effects aggregated within project first",
    "projects then receive equal weight",
    "",
    paste0(
      "minimum_target_cells_per_library=",
      MIN_TARGET_CELLS
    ),
    paste0(
      "minimum_baseline_cells_per_library=",
      MIN_BASELINE_CELLS
    ),
    paste0(
      "minimum_libraries_per_project=",
      MIN_LIBRARIES_PER_PROJECT
    ),
    paste0(
      "minimum_projects_for_consensus_marker=",
      MIN_PROJECTS
    ),
    "",
    "integrated assay NOT used for marker evidence",
    "cellranger_count not accessed",
    "gzip_integrity=PASS",
    "re_read_check=PASS"
  ),
  done_file
)

cat(
  "\n===== SUPPORT =====\n"
)

print(
  support,
  row.names=FALSE
)

cat(
  "\n===== COMPACT TOP MARKERS =====\n"
)

print(
  compact,
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Erythroid native marker aggregation completed\n"
)
