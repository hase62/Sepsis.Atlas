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

set.seed(20260817)

DISCOVERY_PROJECT <- "GSE216007"

N_FULL_EXPECTED <- 27775L
N_DISCOVERY_EXPECTED <- 26962L
N_REPLICATION_SUPPORT_EXPECTED <- 813L

N_HVG <- 3000L
N_PCS <- 30L
RESOLUTIONS <- c(0.2, 0.4, 0.6)

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
    "progenitor_dominant_project_native_clustering_v1__",
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
  "PROGENITOR_DOMINANT_PROJECT_NATIVE_CLUSTERING_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop("Already completed: ", done_file)
}

stopifnot(
  dir.exists(transfer_dir),
  file.exists(lib_file)
)

# ============================================================
# Read full-primary Neutrophil metadata
# ============================================================

files <- list.files(
  transfer_dir,
  pattern="\\.tsv\\.gz$",
  full.names=TRUE
)

stopifnot(
  length(files) == 158L
)

lst <- list()

for(i in seq_along(files)){

  d <- read.delim(
    files[[i]],
    sep="\t",
    quote="\"",
    comment.char="",
    stringsAsFactors=FALSE,
    check.names=FALSE
  )

  req <- c(
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
    all(req %in% names(d))
  )

  d <- d[
    as.character(
      d$integration_compartment_primary_v1
    ) == "Progenitor",
    req,
    drop=FALSE
  ]

  if(nrow(d)){
    lst[[length(lst)+1L]] <- d
  }
}

neu <- do.call(
  rbind,
  lst
)

rm(lst)

neu$project_id <- as.character(neu$project_id)
neu$library_key <- as.character(neu$library_key)
neu$original_cell <- as.character(neu$original_cell)
neu$global_cell <- as.character(neu$global_cell)

stopifnot(
  nrow(neu) == N_FULL_EXPECTED,
  !anyDuplicated(neu$global_cell)
)

rownames(neu) <- neu$global_cell

project_counts <- sort(
  table(neu$project_id),
  decreasing=TRUE
)

stopifnot(
  DISCOVERY_PROJECT %in% names(project_counts),
  length(project_counts) == 9L,
  as.integer(
    project_counts[[DISCOVERY_PROJECT]]
  ) == N_DISCOVERY_EXPECTED
)

eligible_projects <- DISCOVERY_PROJECT

excluded_projects <- setdiff(
  names(project_counts),
  DISCOVERY_PROJECT
)

stopifnot(
  length(eligible_projects) == 1L,
  length(excluded_projects) == 8L,
  sum(
    project_counts[excluded_projects]
  ) == N_REPLICATION_SUPPORT_EXPECTED
)

# ============================================================
# Library lookup
# ============================================================

lib <- read.delim(
  lib_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

stopifnot(
  all(c(
    "project_id",
    "library_key",
    "final_rds_resolved"
  ) %in% names(lib))
)

lib$key <- paste(
  lib$project_id,
  lib$library_key,
  sep="|||"
)

neu$key <- paste(
  neu$project_id,
  neu$library_key,
  sep="|||"
)

# ============================================================
# Helpers
# ============================================================

safe_name <- function(x){
  gsub(
    "[^A-Za-z0-9_.-]",
    "_",
    x
  )
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
  res
){

  pat <- paste0(
    "res\\.",
    gsub(
      "\\.",
      "\\\\.",
      as.character(res)
    ),
    "$"
  )

  z <- grep(
    pat,
    names(md),
    value=TRUE
  )

  if(length(z) != 1L){
    stop(
      "Could not uniquely identify resolution ",
      res,
      "; candidates=",
      paste(z, collapse=",")
    )
  }

  z
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
    c("-t", path),
    stdout=FALSE,
    stderr=FALSE
  )

  if(status != 0L){
    stop(
      "gzip integrity failure: ",
      path
    )
  }
}

# ============================================================
# Dominant-project native discovery
#
# No integration, Harmony, CCA, RPCA or batch correction.
# GSE216007 alone is clustered in native RNA space.
# All other projects are withheld for independent
# replication/support evidence.
# ============================================================

assignment_list <- list()
summary_list <- list()
pca_files <- character()

for(p in eligible_projects){

  cat(
    "\n============================================\n",
    "PROJECT: ", p, "\n",
    "============================================\n",
    sep=""
  )

  mdp <- neu[
    neu$project_id == p,
    ,
    drop=FALSE
  ]

  cat(
    "n_cells=",
    nrow(mdp),
    "\n"
  )

  library_keys <- unique(
    mdp$library_key
  )

  matrices <- list()
  feature_ref <- NULL

  for(lk in library_keys){

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

    mdlib <- mdp[
      mdp$library_key == lk,
      ,
      drop=FALSE
    ]

    obj0 <- readRDS(
      source_rds
    )

    counts0 <- get_rna_counts(
      obj0
    )

    cells <- as.character(
      mdlib$original_cell
    )

    if(!all(
      cells %in%
        colnames(counts0)
    )){
      stop(
        "Missing source cells: ",
        key
      )
    }

    m <- counts0[
      ,
      cells,
      drop=FALSE
    ]

    if(is.null(feature_ref)){

      feature_ref <- rownames(m)

    } else {

      if(
        !setequal(
          feature_ref,
          rownames(m)
        )
      ){
        stop(
          "Feature set mismatch: ",
          key
        )
      }

      m <- m[
        feature_ref,
        ,
        drop=FALSE
      ]
    }

    colnames(m) <- as.character(
      mdlib$global_cell
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
      nrow(mdp),
    !anyDuplicated(
      colnames(counts)
    ),
    setequal(
      colnames(counts),
      mdp$global_cell
    )
  )

  # Restore deterministic metadata order.
  counts <- counts[
    ,
    mdp$global_cell,
    drop=FALSE
  ]

  # ----------------------------------------------------------
  # Native Seurat preprocessing
  # ----------------------------------------------------------

  sobj <- CreateSeuratObject(
    counts=counts,
    project=p,
    min.cells=0,
    min.features=0
  )

  sobj$project_id <- p

  sobj$library_key <- as.character(
    mdp[
      colnames(sobj),
      "library_key"
    ]
  )

  sobj$condition_binary <- as.character(
    mdp[
      colnames(sobj),
      "condition_binary"
    ]
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
    nfeatures=N_HVG,
    verbose=FALSE
  )

  hvg_raw <- VariableFeatures(
    sobj
  )

  hvg <- hvg_raw[
    !technical_gene(
      hvg_raw
    )
  ]

  if(length(hvg) < 1500L){
    stop(
      "Too few nontechnical HVGs for ",
      p,
      ": ",
      length(hvg)
    )
  }

  VariableFeatures(sobj) <- hvg

  cat(
    "HVG raw=",
    length(hvg_raw),
    " retained=",
    length(hvg),
    "\n"
  )

  sobj <- ScaleData(
    sobj,
    features=hvg,
    verbose=FALSE
  )

  sobj <- RunPCA(
    sobj,
    features=hvg,
    npcs=N_PCS,
    verbose=FALSE
  )

  sobj <- FindNeighbors(
    sobj,
    reduction="pca",
    dims=seq_len(N_PCS),
    verbose=FALSE
  )

  sobj <- FindClusters(
    sobj,
    resolution=RESOLUTIONS,
    algorithm=1,
    random.seed=20260817,
    verbose=FALSE
  )

  md_out <- sobj[[]]

  c02 <- find_resolution_column(
    md_out,
    0.2
  )

  c04 <- find_resolution_column(
    md_out,
    0.4
  )

  c06 <- find_resolution_column(
    md_out,
    0.6
  )

  assignment <- data.frame(
    global_cell=rownames(md_out),
    project_id=p,
    library_key=as.character(
      md_out$library_key
    ),
    condition_binary=as.character(
      md_out$condition_binary
    ),
    native_cluster_r0p2=
      as.character(
        md_out[[c02]]
      ),
    native_cluster_r0p4=
      as.character(
        md_out[[c04]]
      ),
    native_cluster_r0p6=
      as.character(
        md_out[[c06]]
      ),
    stringsAsFactors=FALSE
  )

  stopifnot(
    nrow(assignment) ==
      nrow(mdp),
    !anyDuplicated(
      assignment$global_cell
    )
  )

  assignment_list[[
    length(assignment_list)+1L
  ]] <- assignment

  # ----------------------------------------------------------
  # Resolution diagnostics
  # ----------------------------------------------------------

  for(res_name in c(
    "native_cluster_r0p2",
    "native_cluster_r0p4",
    "native_cluster_r0p6"
  )){

    tt <- table(
      assignment[[res_name]]
    )

    summary_list[[
      length(summary_list)+1L
    ]] <- data.frame(
      project_id=p,
      resolution=res_name,
      n_cells=nrow(assignment),
      n_libraries=
        length(
          unique(
            assignment$library_key
          )
        ),
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
      largest_cluster_fraction=
        max(
          as.integer(tt)
        ) /
        nrow(assignment),
      n_clusters_lt20=
        sum(
          as.integer(tt) < 20L
        ),
      n_clusters_lt50=
        sum(
          as.integer(tt) < 50L
        ),
      stringsAsFactors=FALSE
    )
  }

  # ----------------------------------------------------------
  # Save PCA only; do not duplicate full RNA matrix
  # ----------------------------------------------------------

  pca <- Embeddings(
    sobj,
    reduction="pca"
  )

  pca_file <- file.path(
    out_dir,
    paste0(
      safe_name(p),
      "_native_pca30_v1.rds"
    )
  )

  saveRDS(
    pca,
    pca_file,
    compress=FALSE
  )

  pca_files <- c(
    pca_files,
    pca_file
  )

  # Per-project assignments.
  write_gz_tsv(
    assignment,
    file.path(
      out_dir,
      paste0(
        safe_name(p),
        "_native_cluster_assignments_v1.tsv.gz"
      )
    )
  )

  rm(
    counts,
    sobj,
    pca,
    assignment,
    md_out
  )

  invisible(
    gc()
  )
}

# ============================================================
# Combined assignments
# ============================================================

assignments <- do.call(
  rbind,
  assignment_list
)

summary_df <- do.call(
  rbind,
  summary_list
)

n_discovery <- sum(
  project_counts[
    eligible_projects
  ]
)

stopifnot(
  nrow(assignments) ==
    n_discovery,
  n_discovery ==
    N_DISCOVERY_EXPECTED,
  !anyDuplicated(
    assignments$global_cell
  )
)

# ============================================================
# Non-discovery projects retained for independent replication/support
# ============================================================

tiny <- neu[
  neu$project_id %in%
    excluded_projects,
  c(
    "global_cell",
    "project_id",
    "library_key",
    "original_cell",
    "condition_binary",
    "integration_celltype_full_v1",
    "transfer_confidence"
  ),
  drop=FALSE
]

stopifnot(
  nrow(tiny) == N_REPLICATION_SUPPORT_EXPECTED
)

# ============================================================
# Cluster-size detail
# ============================================================

cluster_detail <- list()

for(p in eligible_projects){

  a <- assignments[
    assignments$project_id == p,
    ,
    drop=FALSE
  ]

  for(res in c(
    "native_cluster_r0p2",
    "native_cluster_r0p4",
    "native_cluster_r0p6"
  )){

    tt <- as.data.frame(
      table(
        cluster=a[[res]]
      ),
      stringsAsFactors=FALSE
    )

    tt$project_id <- p
    tt$resolution <- res
    tt$fraction <-
      tt$Freq /
      nrow(a)

    tt <- tt[
      ,
      c(
        "project_id",
        "resolution",
        "cluster",
        "Freq",
        "fraction"
      )
    ]

    cluster_detail[[
      length(cluster_detail)+1L
    ]] <- tt
  }
}

cluster_detail <- do.call(
  rbind,
  cluster_detail
)

# ============================================================
# Write
# ============================================================

write_gz_tsv(
  assignments,
  file.path(
    out_dir,
    "progenitor_dominant_project_native_cluster_assignments_v1.tsv.gz"
  )
)

write.table(
  summary_df,
  file=file.path(
    out_dir,
    "neutrophil_project_native_clustering_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  cluster_detail,
  file=file.path(
    out_dir,
    "progenitor_dominant_project_native_cluster_sizes_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  tiny,
  file=file.path(
    out_dir,
    "progenitor_non_discovery_projects_replication_support_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

project_manifest <- data.frame(
  project_id=names(project_counts),
  n_cells=as.integer(project_counts),
  discovery_eligible=
    names(project_counts) %in%
      eligible_projects,
  policy=ifelse(
    names(project_counts) %in%
      eligible_projects,
    "dominant_project_native_discovery",
    "reserved_for_independent_replication_support"
  ),
  stringsAsFactors=FALSE
)

write.table(
  project_manifest,
  file=file.path(
    out_dir,
    "neutrophil_project_native_clustering_manifest_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

writeLines(
  c(
    "PASS",
    "Progenitor dominant-project native clustering v1",
    paste0(
      "n_full_progenitor_cells=",
      N_FULL_EXPECTED
    ),
    paste0(
      "n_discovery_cells=",
      n_discovery
    ),
    paste0(
      "n_replication_support_cells=",
      sum(
        project_counts[excluded_projects]
      )
    ),
    paste0(
      "n_discovery_projects=",
      length(eligible_projects)
    ),
    paste0(
      "n_replication_support_projects=",
      length(excluded_projects)
    ),
    paste0(
      "discovery_project=",
      DISCOVERY_PROJECT
    ),
    paste0(
      "replication_support_projects=",
      paste(
        excluded_projects,
        collapse=","
      )
    ),
    "",
    "normalization=Seurat LogNormalize",
    "HVG=3000 before technical-gene exclusion",
    "PCA=30",
    "resolutions=0.2,0.4,0.6",
    "",
    "GSE216007 clustered independently in native RNA space",
    "all other projects withheld from discovery clustering",
    "non-discovery projects reserved for independent replication/support",
    "no cross-project integration",
    "no Harmony",
    "no CCA",
    "no RPCA",
    "no regression of biological state",
    "globins excluded from clustering HVGs but retained in source native RNA",
    "native SoupX-corrected RNA retained as source",
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
  "\n===== CLUSTERING SUMMARY =====\n"
)

print(
  summary_df,
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Neutrophil project-native clustering completed\n"
)

