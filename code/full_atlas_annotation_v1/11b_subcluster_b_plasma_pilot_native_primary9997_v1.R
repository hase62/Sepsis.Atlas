#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

set.seed(20260814)

# ============================================================
# Paths
# ============================================================

pilot_file <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_annotation_v1",
  "b_plasma_pilot_primary_intersection_v1.tsv.gz"
)

full_bp_file <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_annotation_v1",
  "b_plasma_primary_cells_v1.tsv.gz"
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
  "b_plasma_native_subclustering_v1_primary9997"
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

for(f in c(
  pilot_file,
  full_bp_file,
  lib_file
)){
  if(!file.exists(f)){
    stop("Missing input: ", f)
  }
}

raw_checkpoint <- file.path(
  out_dir,
  "b_plasma_pilot_native_raw_v1.rds"
)

final_rds <- file.path(
  out_dir,
  "b_plasma_pilot_native_subclustering_v1.rds"
)

if(file.exists(final_rds)){
  stop(
    "Final output already exists; refusing to overwrite: ",
    final_rds
  )
}

# ============================================================
# Pilot manifest
# ============================================================

pilot <- read.delim(
  gzfile(pilot_file),
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

stopifnot(
  nrow(pilot) == 9997,
  length(unique(pilot$global_cell)) == 9997,
  all(c(
    "global_cell",
    "project_id",
    "library_key",
    "original_cell",
    "seurat_cluster"
  ) %in% names(pilot))
)

# Bring condition from authoritative full B/plasma cell set.
full_bp <- read.delim(
  gzfile(full_bp_file),
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

mfull <- match(
  pilot$global_cell,
  full_bp$global_cell
)

if(anyNA(mfull)){
  stop(
    "Pilot B/plasma cells absent from full primary set: ",
    sum(is.na(mfull))
  )
}

pilot$condition_binary <-
  full_bp$condition_binary[mfull]

pilot$clinical_domain_std_full <-
  full_bp$clinical_domain_std[mfull]

pilot$transfer_confidence_full <-
  full_bp$transfer_confidence[mfull]

rm(full_bp)
gc(verbose=FALSE)

# ============================================================
# Library table
# ============================================================

lib <- read.delim(
  lib_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

required_lib <- c(
  "library_key",
  "project_id",
  "final_rds_resolved"
)

miss <- setdiff(
  required_lib,
  names(lib)
)

if(length(miss)){
  stop(
    "Missing resolved-library columns: ",
    paste(miss, collapse=", ")
  )
}

pilot_libs <- unique(
  as.character(pilot$library_key)
)

li <- match(
  pilot_libs,
  as.character(lib$library_key)
)

if(anyNA(li)){
  stop(
    "Pilot library missing from resolved library table: ",
    paste(
      pilot_libs[is.na(li)],
      collapse=", "
    )
  )
}

lib <- lib[li,,drop=FALSE]

cat(
  "Pilot cells     =", nrow(pilot), "\n",
  "Pilot libraries =", nrow(lib), "\n",
  "Pilot projects  =",
  length(unique(pilot$project_id)),
  "\n"
)

# ============================================================
# Helpers
# ============================================================

resolve_rds_path <- function(f){

  f <- as.character(f)

  if(file.exists(f)){
    return(
      normalizePath(
        f,
        mustWork=TRUE
      )
    )
  }

  g <- file.path(
    root,
    f
  )

  if(file.exists(g)){
    return(
      normalizePath(
        g,
        mustWork=TRUE
      )
    )
  }

  stop(
    "Final QC RDS not found: ",
    f
  )
}

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
    stop(
      "Could not retrieve native RNA counts"
    )
  }

  z
}

# ============================================================
# Reconstruct native SoupX RNA
#
# Reuse checkpoint if a previous attempt already completed
# the expensive per-library extraction.
# ============================================================

if(file.exists(raw_checkpoint)){

  cat(
    "\n===== REUSING NATIVE RAW CHECKPOINT =====\n"
  )

  bp <- readRDS(
    raw_checkpoint
  )

  stopifnot(
    ncol(bp) == 9997,
    nrow(bp) == 38606,
    identical(
      colnames(bp),
      pilot$global_cell
    )
  )

} else {

  cat(
    "\n===== RECONSTRUCTING NATIVE RNA =====\n"
  )

  mats <- vector(
    "list",
    nrow(lib)
  )

  features <- NULL

  for(i in seq_len(nrow(lib))){

    key <- as.character(
      lib$library_key[i]
    )

    d <- pilot[
      pilot$library_key == key,
      ,
      drop=FALSE
    ]

    if(!nrow(d)){
      next
    }

    f <- resolve_rds_path(
      lib$final_rds_resolved[i]
    )

    cat(
      sprintf(
        "[%3d/%3d] %-50s cells=%d\n",
        i,
        nrow(lib),
        key,
        nrow(d)
      )
    )

    obj <- readRDS(f)

    counts <- get_rna_counts(
      obj
    )

    if(is.null(features)){

      features <- rownames(
        counts
      )

      if(length(features) != 38606){
        stop(
          "Expected 38606 features; got ",
          length(features),
          " in ",
          key
        )
      }

    } else {

      if(!identical(
        rownames(counts),
        features
      )){
        stop(
          "Feature order mismatch: ",
          key
        )
      }
    }

    cells_local <- as.character(
      d$original_cell
    )

    if(!all(
      cells_local %in%
        colnames(counts)
    )){
      missing_cells <- cells_local[
        !cells_local %in%
          colnames(counts)
      ]

      stop(
        "Pilot cells absent from final QC RDS ",
        key,
        ": ",
        paste(
          head(missing_cells,20),
          collapse=", "
        )
      )
    }

    mat <- counts[
      ,
      cells_local,
      drop=FALSE
    ]

    colnames(mat) <-
      as.character(
        d$global_cell
      )

    mats[[i]] <- mat

    rm(
      obj,
      counts,
      mat
    )

    gc(verbose=FALSE)
  }

  mats <- mats[
    !vapply(
      mats,
      is.null,
      logical(1)
    )
  ]

  counts_all <- do.call(
    cbind,
    mats
  )

  rm(mats)
  gc(verbose=FALSE)

  stopifnot(
    nrow(counts_all) == 38606,
    ncol(counts_all) == 9997,
    length(unique(
      colnames(counts_all)
    )) == 9997
  )

  mi <- match(
    pilot$global_cell,
    colnames(counts_all)
  )

  if(anyNA(mi)){
    stop(
      "Failed to restore pilot cell order"
    )
  }

  counts_all <- counts_all[
    ,
    mi,
    drop=FALSE
  ]

  stopifnot(
    identical(
      colnames(counts_all),
      pilot$global_cell
    )
  )

  md <- pilot

  rownames(md) <-
    md$global_cell

  bp <- CreateSeuratObject(
    counts=counts_all,
    project="B_plasma_pilot_native_v1",
    assay="RNA",
    min.cells=0,
    min.features=0,
    meta.data=md
  )

  stopifnot(
    nrow(bp) == 38606,
    ncol(bp) == 9997,
    identical(
      colnames(bp),
      pilot$global_cell
    )
  )

  saveRDS(
    bp,
    raw_checkpoint,
    compress=FALSE
  )

  cat(
    "\nNative raw checkpoint written:\n",
    raw_checkpoint,
    "\n"
  )

  rm(counts_all)
  gc(verbose=FALSE)
}

# ============================================================
# B/plasma-specific normalization / HVG
# ============================================================

cat(
  "\n===== NORMALIZE / HVG =====\n"
)

DefaultAssay(bp) <- "RNA"

bp <- NormalizeData(
  bp,
  normalization.method="LogNormalize",
  scale.factor=10000,
  verbose=FALSE
)

# Start with a larger ranked candidate set, then remove
# technical / clonotype-heavy features and retain 3000.
bp <- FindVariableFeatures(
  bp,
  selection.method="vst",
  nfeatures=5000,
  verbose=FALSE
)

vf_initial <- VariableFeatures(
  bp
)

technical <- grepl(
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
  vf_initial
)

# Exclude rearranged immunoglobulin variable/joining genes
# from PCA only. Constant/isotype genes such as IGHM,
# IGHD, IGHA*, IGHG* and JCHAIN remain eligible.
ig_rearranged <- grepl(
  paste0(
    "^IGHV|",
    "^IGKV|",
    "^IGLV|",
    "^IGHJ[0-9]|",
    "^IGKJ[0-9]|",
    "^IGLJ[0-9]"
  ),
  vf_initial
)

eligible <- !technical &
  !ig_rearranged

vf <- head(
  vf_initial[eligible],
  3000
)

if(length(vf) < 3000){
  stop(
    "Fewer than 3000 usable HVGs after filtering: ",
    length(vf)
  )
}

VariableFeatures(bp) <- vf

hvg_table <- data.frame(
  gene=vf_initial,
  technical_gene=technical,
  rearranged_ig_gene=ig_rearranged,
  selected_for_pca=
    vf_initial %in% vf,
  stringsAsFactors=FALSE
)

write.table(
  hvg_table,
  file=file.path(
    out_dir,
    "b_plasma_hvg_selection_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

write.table(
  data.frame(
    gene=vf,
    stringsAsFactors=FALSE
  ),
  file=file.path(
    out_dir,
    "b_plasma_hvg3000_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# PCA
# ============================================================

cat(
  "\n===== SCALE / PCA =====\n"
)

bp <- ScaleData(
  bp,
  features=vf,
  verbose=FALSE
)

bp <- RunPCA(
  bp,
  features=vf,
  npcs=50,
  seed.use=20260814,
  verbose=FALSE
)

pca_stdev <- Stdev(
  bp[["pca"]]
)

write.table(
  data.frame(
    PC=seq_along(pca_stdev),
    stdev=pca_stdev
  ),
  file=file.path(
    out_dir,
    "b_plasma_pca_stdev_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# kNN / SNN
# ============================================================

cat(
  "\n===== NEIGHBORS =====\n"
)

bp <- FindNeighbors(
  bp,
  reduction="pca",
  dims=1:30,
  k.param=30,
  nn.method="annoy",
  n.trees=100,
  annoy.metric="euclidean",
  compute.SNN=TRUE,
  prune.SNN=1/15,
  verbose=FALSE
)

# ============================================================
# Leiden clustering over a fixed resolution grid
# ============================================================

cat(
  "\n===== LEIDEN =====\n"
)

resolutions <- c(
  0.2,
  0.4,
  0.6,
  0.8,
  1.0,
  1.2
)

cluster_columns <- character()

for(r in resolutions){

  nm <- paste0(
    "cluster_res_",
    gsub(
      "\\.",
      "p",
      sprintf("%.1f", r)
    )
  )

  cat(
    "resolution ",
    r,
    " -> ",
    nm,
    "\n",
    sep=""
  )

  bp <- FindClusters(
    bp,
    graph.name="RNA_snn",
    algorithm=4,
    leiden_method="igraph",
    leiden_objective_function="modularity",
    resolution=r,
    n.iter=20,
    random.seed=20260814,
    verbose=FALSE
  )

  bp[[nm]] <-
    as.character(
      Idents(bp)
    )

  cluster_columns <- c(
    cluster_columns,
    nm
  )
}

# ============================================================
# UMAP = visualization only
# ============================================================

cat(
  "\n===== UMAP =====\n"
)

bp <- RunUMAP(
  bp,
  reduction="pca",
  dims=1:30,
  n.neighbors=30,
  min.dist=0.3,
  metric="cosine",
  seed.use=20260814,
  reduction.name="umap",
  verbose=FALSE
)

# ============================================================
# Cluster diagnostics
# ============================================================

normalized_entropy <- function(x){

  tt <- table(x)

  if(length(tt) <= 1){
    return(0)
  }

  p <- as.numeric(tt) /
    sum(tt)

  -sum(
    p * log(p)
  ) / log(length(tt))
}

diagnostics <- list()

md <- bp[[]]

for(j in seq_along(
  cluster_columns
)){

  nm <- cluster_columns[j]

  resolution <- resolutions[j]

  cl <- as.character(
    md[[nm]]
  )

  for(k in sort(unique(cl))){

    ii <- which(
      cl == k
    )

    pp <- sort(
      table(
        md$project_id[ii]
      ),
      decreasing=TRUE
    )

    ll <- table(
      md$library_key[ii]
    )

    old <- sort(
      table(
        md$seurat_cluster[ii]
      ),
      decreasing=TRUE
    )

    condition <- table(
      md$condition_binary[ii]
    )

    healthy_fraction <-
      if("healthy" %in%
         names(condition)){
        as.numeric(
          condition["healthy"]
        ) / length(ii)
      } else {
        0
      }

    disease_fraction <-
      if("disease" %in%
         names(condition)){
        as.numeric(
          condition["disease"]
        ) / length(ii)
      } else {
        0
      }

    diagnostics[[
      length(diagnostics)+1L
    ]] <- data.frame(
      resolution=resolution,
      cluster=k,
      n_cells=length(ii),
      n_libraries=length(ll),
      n_projects=length(pp),
      top_project=names(pp)[1],
      max_project_fraction=
        as.numeric(pp[1]) /
        length(ii),
      project_entropy_normalized=
        normalized_entropy(
          md$project_id[ii]
        ),
      top_old_global_cluster=
        names(old)[1],
      old_global_cluster_purity=
        as.numeric(old[1]) /
        length(ii),
      healthy_fraction=
        healthy_fraction,
      disease_fraction=
        disease_fraction,
      stringsAsFactors=FALSE
    )
  }
}

diagnostics <- do.call(
  rbind,
  diagnostics
)

write.table(
  diagnostics,
  file=file.path(
    out_dir,
    "b_plasma_cluster_diagnostics_all_resolutions_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# Resolution-level overview
overview <- do.call(
  rbind,
  lapply(
    resolutions,
    function(r){

      d <- diagnostics[
        diagnostics$resolution == r,
        ,
        drop=FALSE
      ]

      data.frame(
        resolution=r,
        n_clusters=nrow(d),
        min_cluster_cells=
          min(d$n_cells),
        median_cluster_cells=
          median(d$n_cells),
        max_cluster_cells=
          max(d$n_cells),
        median_max_project_fraction=
          median(
            d$max_project_fraction
          ),
        median_project_entropy=
          median(
            d$project_entropy_normalized
          ),
        median_old_global_cluster_purity=
          median(
            d$old_global_cluster_purity
          ),
        stringsAsFactors=FALSE
      )
    }
  )
)

write.table(
  overview,
  file=file.path(
    out_dir,
    "b_plasma_resolution_overview_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

cat(
  "\n===== RESOLUTION OVERVIEW =====\n"
)

print(
  overview,
  row.names=FALSE
)

# ============================================================
# Cell-level export
# ============================================================

umap <- Embeddings(
  bp,
  "umap"
)

cell_export <- data.frame(
  global_cell=colnames(bp),
  project_id=md$project_id,
  library_key=md$library_key,
  original_cell=md$original_cell,
  condition_binary=
    md$condition_binary,
  old_global_cluster=
    md$seurat_cluster,
  md[,cluster_columns,drop=FALSE],
  UMAP_1=umap[,1],
  UMAP_2=umap[,2],
  stringsAsFactors=FALSE,
  check.names=FALSE
)

cell_file <- file.path(
  out_dir,
  "b_plasma_cluster_assignments_v1.tsv.gz"
)

con <- gzfile(
  cell_file,
  "wt"
)

tryCatch(
  write.table(
    cell_export,
    con,
    sep="\t",
    quote=1,
    qmethod="double",
    row.names=FALSE
  ),
  finally=close(con)
)

# Read-back exact validation.
check <- read.delim(
  gzfile(cell_file),
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

stopifnot(
  nrow(check) == 9997,
  identical(
    check$global_cell,
    colnames(bp)
  )
)

rm(check)

# ============================================================
# Save authoritative discovery object
# ============================================================

saveRDS(
  bp,
  final_rds,
  compress=FALSE
)

writeLines(
  c(
    "PASS",
    "B/plasma native pilot subclustering v1",
    "n_cells=9997",
    "n_features=38606",
    "native_RNA=SoupX-corrected RNA counts",
    "HVG_candidates=5000",
    "HVG_PCA=3000",
    "rearranged_IG_variable_genes_excluded_from_PCA_only",
    "PCA=50 calculated; PCs1-30 used",
    "kNN=30",
    "SNN_prune=1/15",
    "Leiden=igraph algorithm4 modularity",
    "resolutions=0.2,0.4,0.6,0.8,1.0,1.2",
    "UMAP=visualization_only",
    paste0(
      "Seurat=",
      as.character(
        packageVersion("Seurat")
      )
    ),
    paste0(
      "igraph=",
      as.character(
        packageVersion("igraph")
      )
    )
  ),
  file.path(
    out_dir,
    "SUBCLUSTERING_COMPLETE.ok"
  )
)

cat(
  "\nPASS: B/plasma native pilot subclustering completed\n"
)
