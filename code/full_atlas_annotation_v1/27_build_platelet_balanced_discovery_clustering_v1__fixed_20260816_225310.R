#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 2L){
  stop("Usage: script <root> <tag>")
}

root <- normalizePath(
  args[[1]],
  mustWork=TRUE
)

tag <- args[[2]]

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
})

source(
  file.path(
    root,
    "full_atlas_primary_integration_v1",
    "00_common.R"
  )
)

set.seed(20260816)

N_FULL_EXPECTED <- 49585L

MAX_CELLS_PER_PROJECT <- 4000L
MIN_DISCOVERY_PROJECT_CELLS <- 100L

N_VARIABLE_FEATURES <- 2500L
N_INTEGRATION_FEATURES <- 2000L
N_PCS <- 20L

RESOLUTIONS <- c(
  0.2,
  0.4,
  0.6
)

# ============================================================
# Paths
# ============================================================

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
    "platelet_balanced_discovery_clustering_v1__",
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
  "PLATELET_BALANCED_DISCOVERY_CLUSTERING_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop(
    "Already completed: ",
    done_file
  )
}

stopifnot(
  dir.exists(transfer_dir),
  file.exists(lib_file)
)

# ============================================================
# Helpers
# ============================================================

read_tsv <- function(path){

  con <- if(
    grepl("\\.gz$", path)
  ){
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }

  tryCatch(
    {
      read.delim(
        con,
        sep="\t",
        quote="\"",
        comment.char="",
        stringsAsFactors=FALSE,
        check.names=FALSE
      )
    },
    finally={
      close(con)
    }
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

  invisible(TRUE)
}

technical_gene <- function(g){

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

find_resolution_column <- function(
  md,
  resolution
){

  rs <- format(
    resolution,
    nsmall=1,
    trim=TRUE
  )

  pattern <- paste0(
    "res\\.",
    gsub(
      "\\.",
      "\\\\.",
      rs
    ),
    "$"
  )

  x <- grep(
    pattern,
    names(md),
    value=TRUE
  )

  if(length(x) != 1L){
    stop(
      "Could not identify clustering column for resolution ",
      resolution,
      ": ",
      paste(x, collapse=",")
    )
  }

  x
}

# Equal first-pass quota across libraries, then fill remainder
# from unused cells. This prevents one very large library from
# defining the discovery pilot.
sample_project_balanced <- function(
  d,
  cap,
  seed
){

  if(nrow(d) <= cap){
    return(
      seq_len(nrow(d))
    )
  }

  set.seed(seed)

  groups <- split(
    seq_len(nrow(d)),
    d$library_key
  )

  n_lib <- length(groups)

  base_quota <- max(
    1L,
    floor(
      cap / n_lib
    )
  )

  selected <- integer()

  for(g in groups){

    take <- min(
      length(g),
      base_quota
    )

    selected <- c(
      selected,
      sample(
        g,
        take,
        replace=FALSE
      )
    )
  }

  selected <- unique(
    selected
  )

  remaining_n <-
    cap -
    length(selected)

  if(remaining_n > 0L){

    pool <- setdiff(
      seq_len(nrow(d)),
      selected
    )

    take <- min(
      remaining_n,
      length(pool)
    )

    if(take > 0L){

      selected <- c(
        selected,
        sample(
          pool,
          take,
          replace=FALSE
        )
      )
    }
  }

  sort(
    unique(selected)
  )
}

entropy_normalized <- function(
  x,
  n_total_projects
){

  tt <- table(x)

  p <- as.numeric(tt) /
    sum(tt)

  h <- -sum(
    p *
      log(p)
  )

  if(n_total_projects <= 1L){
    return(0)
  }

  h /
    log(n_total_projects)
}

# ============================================================
# Full Platelet/megakaryocyte metadata
# ============================================================

transfer_files <- list.files(
  transfer_dir,
  pattern="\\.tsv\\.gz$",
  full.names=TRUE
)

stopifnot(
  length(transfer_files) == 158L
)

lst <- list()

for(f in transfer_files){

  d <- read_tsv(f)

  required <- c(
    "project_id",
    "library_key",
    "original_cell",
    "global_cell",
    "condition_binary",
    "integration_compartment_primary_v1",
    "integration_celltype_full_v1",
    "transfer_confidence"
  )

  stopifnot(
    all(
      required %in%
        names(d)
    )
  )

  d <- d[
    as.character(
      d$integration_compartment_primary_v1
    ) ==
      "Platelet_megakaryocyte",
    required,
    drop=FALSE
  ]

  if(nrow(d)){
    lst[[
      length(lst)+1L
    ]] <- d
  }
}

meta <- do.call(
  rbind,
  lst
)

rm(lst)

meta$project_id <-
  as.character(
    meta$project_id
  )

meta$library_key <-
  as.character(
    meta$library_key
  )

meta$original_cell <-
  as.character(
    meta$original_cell
  )

meta$global_cell <-
  as.character(
    meta$global_cell
  )

meta$condition_binary <-
  as.character(
    meta$condition_binary
  )

stopifnot(
  nrow(meta) ==
    N_FULL_EXPECTED,
  !anyDuplicated(
    meta$global_cell
  )
)

# ============================================================
# Project manifest
# ============================================================

project_counts <- sort(
  table(meta$project_id),
  decreasing=TRUE
)

discovery_projects <- names(
  project_counts[
    project_counts >=
      MIN_DISCOVERY_PROJECT_CELLS
  ]
)

deferred_projects <- names(
  project_counts[
    project_counts <
      MIN_DISCOVERY_PROJECT_CELLS
  ]
)

if(length(discovery_projects) < 2L){
  stop(
    "Fewer than two projects are eligible for discovery integration"
  )
}

# ============================================================
# Deterministic project/library-balanced sampling
# ============================================================

sampled_list <- list()
library_manifest_list <- list()

project_order <- sort(
  discovery_projects
)

for(pi in seq_along(
  project_order
)){

  p <- project_order[[pi]]

  d <- meta[
    meta$project_id == p,
    ,
    drop=FALSE
  ]

  idx <- sample_project_balanced(
    d,
    cap=MAX_CELLS_PER_PROJECT,
    seed=20260816L + pi
  )

  ds <- d[
    idx,
    ,
    drop=FALSE
  ]

  sampled_list[[
    length(sampled_list)+1L
  ]] <- ds

  all_libs <- sort(
    unique(
      d$library_key
    )
  )

  for(lk in all_libs){

    n_full <- sum(
      d$library_key == lk
    )

    n_sampled <- sum(
      ds$library_key == lk
    )

    library_manifest_list[[
      length(library_manifest_list)+1L
    ]] <- data.frame(
      project_id=p,
      library_key=lk,
      n_full=n_full,
      n_discovery=n_sampled,
      discovery_fraction=
        n_sampled /
        n_full,
      stringsAsFactors=FALSE
    )
  }
}

pilot_meta <- do.call(
  rbind,
  sampled_list
)

rm(sampled_list)

library_manifest <- do.call(
  rbind,
  library_manifest_list
)

rm(library_manifest_list)

stopifnot(
  !anyDuplicated(
    pilot_meta$global_cell
  )
)

pilot_project_counts <- table(
  pilot_meta$project_id
)

project_manifest <- data.frame(
  project_id=
    names(project_counts),

  n_full=
    as.integer(
      project_counts
    ),

  discovery_eligible=
    names(project_counts) %in%
      discovery_projects,

  n_discovery=
    as.integer(
      pilot_project_counts[
        names(project_counts)
      ]
    ),

  stringsAsFactors=FALSE
)

project_manifest$n_discovery[
  is.na(
    project_manifest$n_discovery
  )
] <- 0L

project_manifest$policy <- ifelse(
  project_manifest$discovery_eligible,
  "balanced_discovery_integration",
  "deferred_from_discovery"
)

cat(
  "\n===== PROJECT MANIFEST =====\n"
)

print(
  project_manifest,
  row.names=FALSE
)

cat(
  "\nDiscovery pilot cells = ",
  nrow(pilot_meta),
  "\n",
  sep=""
)

# ============================================================
# Resolved source-library lookup
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
    ) %in%
      names(lib)
  )
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

# ============================================================
# Build native RNA object for each project
# ============================================================

object_list <- list()

feature_reference <- NULL

for(p in project_order){

  cat(
    "\n============================================\n",
    "BUILD PROJECT: ",
    p,
    "\n",
    "============================================\n",
    sep=""
  )

  mp <- pilot_meta[
    pilot_meta$project_id == p,
    ,
    drop=FALSE
  ]

  project_libraries <- sort(
    unique(
      mp$library_key
    )
  )

  matrices <- list()

  for(lk in project_libraries){

    ml <- mp[
      mp$library_key == lk,
      ,
      drop=FALSE
    ]

    key <- paste(
      p,
      lk,
      sep="|||"
    )

    ii <- match(
      key,
      lib$key
    )

    if(is.na(ii)){
      stop(
        "Library not found: ",
        key
      )
    }

    source_rds <- as.character(
      lib$final_rds_resolved[[ii]]
    )

    if(!file.exists(source_rds)){
      stop(
        "Missing source RDS: ",
        source_rds
      )
    }

    obj0 <- readRDS(
      source_rds
    )

    counts0 <- get_rna_counts(
      obj0
    )

    cells <- as.character(
      ml$original_cell
    )

    if(
      !all(
        cells %in%
          colnames(counts0)
      )
    ){
      stop(
        "Missing cells in source RDS: ",
        key
      )
    }

    m <- counts0[
      ,
      cells,
      drop=FALSE
    ]

    if(is.null(
      feature_reference
    )){

      feature_reference <-
        rownames(m)

    } else {

      if(
        !setequal(
          feature_reference,
          rownames(m)
        )
      ){
        stop(
          "Feature-set mismatch: ",
          key
        )
      }

      m <- m[
        feature_reference,
        ,
        drop=FALSE
      ]
    }

    colnames(m) <-
      as.character(
        ml$global_cell
      )

    matrices[[
      length(matrices)+1L
    ]] <- m

    rm(
      obj0,
      counts0,
      m
    )

    invisible(
      gc()
    )
  }

  counts <- do.call(
    cbind,
    matrices
  )

  rm(matrices)

  stopifnot(
    ncol(counts) ==
      nrow(mp),
    setequal(
      colnames(counts),
      mp$global_cell
    )
  )

  counts <- counts[
    ,
    mp$global_cell,
    drop=FALSE
  ]

  md <- mp[
    match(
      colnames(counts),
      mp$global_cell
    ),
    ,
    drop=FALSE
  ]

  rownames(md) <-
    md$global_cell

  sobj <- CreateSeuratObject(
    counts=counts,
    project=p,
    meta.data=md,
    min.cells=0,
    min.features=0
  )

  sobj <- NormalizeData(
    sobj,
    normalization.method="LogNormalize",
    scale.factor=10000,
    verbose=FALSE
  )

  sobj <- FindVariableFeatures(
    sobj,
    selection.method="vst",
    nfeatures=N_VARIABLE_FEATURES,
    verbose=FALSE
  )

  object_list[[p]] <- sobj

  cat(
    "cells=",
    ncol(sobj),
    " HVG=",
    length(
      VariableFeatures(sobj)
    ),
    "\n",
    sep=""
  )

  rm(
    counts,
    sobj
  )

  invisible(
    gc()
  )
}

# ============================================================
# Project-balanced integration features
# ============================================================

integration_features <- SelectIntegrationFeatures(
  object.list=object_list,
  nfeatures=N_INTEGRATION_FEATURES
)

integration_features <-
  integration_features[
    !technical_gene(
      integration_features
    )
  ]

if(length(
  integration_features
) < 1200L){
  stop(
    "Too few nontechnical integration features: ",
    length(integration_features)
  )
}

cat(
  "\nIntegration features retained = ",
  length(integration_features),
  "\n",
  sep=""
)

# ============================================================
# RPCA preparation
#
# No SCTransform.
# No nCount regression.
# No percent.mt regression.
# No condition regression.
# ============================================================

for(p in names(
  object_list
)){

  object_list[[p]] <- ScaleData(
    object_list[[p]],
    features=integration_features,
    verbose=FALSE
  )

  object_list[[p]] <- RunPCA(
    object_list[[p]],
    features=integration_features,
    npcs=N_PCS,
    verbose=FALSE
  )
}

min_project_n <- min(
  vapply(
    object_list,
    ncol,
    numeric(1)
  )
)

if(min_project_n < 60L){
  stop(
    "Smallest discovery project unexpectedly too small: ",
    min_project_n
  )
}

anchors <- FindIntegrationAnchors(
  object.list=object_list,
  anchor.features=integration_features,
  reduction="rpca",
  dims=seq_len(N_PCS),
  k.anchor=5,
  k.filter=min(
    50L,
    min_project_n - 1L
  ),
  verbose=TRUE
)

integrated <- IntegrateData(
  anchorset=anchors,
  dims=seq_len(N_PCS),
  k.weight=min(
    50L,
    min_project_n - 1L
  ),
  new.assay.name="integrated",
  verbose=TRUE
)

DefaultAssay(
  integrated
) <- "integrated"

integrated <- ScaleData(
  integrated,
  verbose=FALSE
)

integrated <- RunPCA(
  integrated,
  npcs=N_PCS,
  verbose=FALSE
)

integrated <- FindNeighbors(
  integrated,
  reduction="pca",
  dims=seq_len(N_PCS),
  verbose=FALSE
)

integrated <- FindClusters(
  integrated,
  resolution=RESOLUTIONS,
  algorithm=1,
  random.seed=20260816,
  verbose=FALSE
)

integrated <- RunUMAP(
  integrated,
  reduction="pca",
  dims=seq_len(N_PCS),
  seed.use=20260816,
  verbose=FALSE
)

# ============================================================
# Recover clustering columns
# ============================================================

md <- integrated[[]]

r02 <- find_resolution_column(
  md,
  0.2
)

r04 <- find_resolution_column(
  md,
  0.4
)

r06 <- find_resolution_column(
  md,
  0.6
)

assignments <- data.frame(
  global_cell=
    rownames(md),

  project_id=
    as.character(
      md$project_id
    ),

  library_key=
    as.character(
      md$library_key
    ),

  condition_binary=
    as.character(
      md$condition_binary
    ),

  cluster_r0p2=
    as.character(
      md[[r02]]
    ),

  cluster_r0p4=
    as.character(
      md[[r04]]
    ),

  cluster_r0p6=
    as.character(
      md[[r06]]
    ),

  stringsAsFactors=FALSE
)

stopifnot(
  nrow(assignments) ==
    nrow(pilot_meta),
  !anyDuplicated(
    assignments$global_cell
  )
)

# ============================================================
# Cluster diagnostics
# ============================================================

resolution_map <- c(
  r0p2="cluster_r0p2",
  r0p4="cluster_r0p4",
  r0p6="cluster_r0p6"
)

diagnostic_rows <- list()

N_DISCOVERY_PROJECTS <-
  length(
    unique(
      assignments$project_id
    )
  )

for(res_name in names(
  resolution_map
)){

  col <- resolution_map[[res_name]]

  clusters <- sort(
    unique(
      assignments[[col]]
    )
  )

  for(k in clusters){

    d <- assignments[
      assignments[[col]] == k,
      ,
      drop=FALSE
    ]

    pt <- sort(
      table(d$project_id),
      decreasing=TRUE
    )

    lt <- sort(
      table(d$library_key),
      decreasing=TRUE
    )

    ct <- table(
      d$condition_binary
    )

    diagnostic_rows[[
      length(diagnostic_rows)+1L
    ]] <- data.frame(
      resolution=res_name,
      cluster=k,
      n_cells=nrow(d),

      n_projects=
        length(pt),

      max_project=
        names(pt)[1],

      max_project_fraction=
        as.numeric(
          pt[1]
        ) /
        nrow(d),

      project_entropy_normalized=
        entropy_normalized(
          d$project_id,
          N_DISCOVERY_PROJECTS
        ),

      n_libraries=
        length(lt),

      max_library=
        names(lt)[1],

      max_library_fraction=
        as.numeric(
          lt[1]
        ) /
        nrow(d),

      healthy_fraction=
        if(
          "healthy" %in%
            names(ct)
        ){
          as.numeric(
            ct[["healthy"]]
          ) /
            nrow(d)
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
            nrow(d)
        } else {
          NA_real_
        },

      stringsAsFactors=FALSE
    )
  }
}

cluster_diagnostics <- do.call(
  rbind,
  diagnostic_rows
)

# ============================================================
# Resolution-level diagnostics
# ============================================================

resolution_summary <- do.call(
  rbind,
  lapply(
    names(
      resolution_map
    ),
    function(res_name){

      col <- resolution_map[[res_name]]

      tt <- table(
        assignments[[col]]
      )

      d <- cluster_diagnostics[
        cluster_diagnostics$resolution ==
          res_name,
        ,
        drop=FALSE
      ]

      data.frame(
        resolution=res_name,
        n_cells=nrow(assignments),
        n_clusters=length(tt),
        smallest_cluster=
          min(
            as.integer(tt)
          ),
        median_cluster=
          median(
            as.integer(tt)
          ),
        largest_cluster=
          max(
            as.integer(tt)
          ),
        median_max_project_fraction=
          median(
            d$max_project_fraction
          ),
        max_cluster_project_fraction=
          max(
            d$max_project_fraction
          ),
        n_clusters_project_fraction_ge_0p8=
          sum(
            d$max_project_fraction >=
              0.8
          ),
        median_project_entropy=
          median(
            d$project_entropy_normalized
          ),
        stringsAsFactors=FALSE
      )
    }
  )
)

# ============================================================
# Save compact discovery object
#
# We intentionally do NOT save the full integrated assay.
# Native RNA remains in original per-library RDS files.
# ============================================================

compact <- list(

  analysis_version=
    "Platelet_balanced_discovery_clustering_v1",

  n_full_cells=
    N_FULL_EXPECTED,

  n_discovery_cells=
    nrow(assignments),

  discovery_projects=
    discovery_projects,

  deferred_projects=
    deferred_projects,

  integration_method=
    "Seurat RPCA LogNormalize balanced discovery",

  integration_features=
    integration_features,

  assignments=
    assignments,

  pca=
    Embeddings(
      integrated,
      reduction="pca"
    ),

  umap=
    Embeddings(
      integrated,
      reduction="umap"
    ),

  project_manifest=
    project_manifest,

  library_sampling_manifest=
    library_manifest,

  cluster_diagnostics=
    cluster_diagnostics,

  resolution_summary=
    resolution_summary,

  rules=list(
    max_cells_per_project=
      MAX_CELLS_PER_PROJECT,
    minimum_discovery_project_cells=
      MIN_DISCOVERY_PROJECT_CELLS,
    project_sampling=
      "library-balanced first-pass quota plus random remainder",
    normalization=
      "LogNormalize",
    integration=
      "RPCA",
    n_pcs=
      N_PCS,
    resolutions=
      RESOLUTIONS,
    integrated_coordinates_role=
      "cluster discovery only",
    biological_annotation_evidence=
      "native RNA library-balanced markers"
  )
)

rds_file <- file.path(
  out_dir,
  "platelet_balanced_discovery_clustering_v1.rds"
)

tmp_file <- paste0(
  rds_file,
  ".tmp.",
  Sys.getpid()
)

saveRDS(
  compact,
  tmp_file,
  compress=TRUE
)

if(
  !file.rename(
    tmp_file,
    rds_file
  )
){
  stop(
    "Failed atomic rename: ",
    rds_file
  )
}

# ============================================================
# Write outputs
# ============================================================

write_gz_tsv(
  assignments,
  file.path(
    out_dir,
    "platelet_balanced_discovery_cluster_assignments_v1.tsv.gz"
  )
)

write.table(
  project_manifest,
  file=file.path(
    out_dir,
    "platelet_balanced_discovery_project_manifest_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  library_manifest,
  file=file.path(
    out_dir,
    "platelet_balanced_discovery_library_sampling_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  cluster_diagnostics,
  file=file.path(
    out_dir,
    "platelet_balanced_discovery_cluster_diagnostics_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  resolution_summary,
  file=file.path(
    out_dir,
    "platelet_balanced_discovery_resolution_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Integrity re-read
# ============================================================

chk <- readRDS(
  rds_file
)

stopifnot(
  identical(
    chk$analysis_version,
    "Platelet_balanced_discovery_clustering_v1"
  ),
  chk$n_full_cells ==
    N_FULL_EXPECTED,
  chk$n_discovery_cells ==
    nrow(assignments),
  nrow(
    chk$assignments
  ) ==
    nrow(assignments)
)

rm(chk)

# ============================================================
# Completion
# ============================================================

writeLines(
  c(
    "PASS",
    "Platelet/megakaryocyte balanced discovery clustering v1",
    "n_full_cells=49585",
    paste0(
      "n_discovery_cells=",
      nrow(assignments)
    ),
    paste0(
      "n_discovery_projects=",
      length(discovery_projects)
    ),
    paste0(
      "n_deferred_projects=",
      length(deferred_projects)
    ),
    "",
    "max_cells_per_project=4000",
    "sampling=library-balanced within project",
    "normalization=LogNormalize",
    "integration=Seurat RPCA",
    "PCA=20",
    "resolutions=0.2,0.4,0.6",
    "",
    "NO SCTransform",
    "NO nCount regression",
    "NO mitochondrial regression",
    "NO condition regression",
    "",
    "integrated coordinates are discovery-only",
    "final biological evidence must come from native RNA",
    "cellranger_count not accessed"
  ),
  done_file
)

cat(
  "\n===== PROJECT MANIFEST =====\n"
)

print(
  project_manifest,
  row.names=FALSE
)

cat(
  "\n===== RESOLUTION SUMMARY =====\n"
)

print(
  resolution_summary,
  row.names=FALSE
)

cat(
  "\n===== CLUSTER DIAGNOSTICS =====\n"
)

print(
  cluster_diagnostics,
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Platelet balanced discovery clustering completed\n"
)
