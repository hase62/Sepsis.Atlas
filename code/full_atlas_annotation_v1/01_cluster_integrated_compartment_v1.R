#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if(length(args) < 2L){
  stop("Usage: 01_cluster_integrated_compartment_v1.R ROOT COMPARTMENT")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
compartment <- as.character(args[[2]])

required <- c("Seurat","SeuratObject","igraph","uwot")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop("Missing packages: ", paste(missing, collapse=", "))

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

set.seed(20260813)

spec <- data.frame(
  compartment=c("T_NK_combined","Monocyte_DC_combined"),
  ruvK=c(5L,2L),
  kPB=c(5L,5L),
  stringsAsFactors=FALSE
)

hit <- which(spec$compartment == compartment)
if(length(hit) != 1L) stop("Unknown compartment: ", compartment)

in_dir <- file.path(
  root, "pre_integration", "full_atlas_primary_v1", "scmerge2_primary_runs"
)

in_rds <- file.path(
  in_dir,
  paste0(
    compartment,
    "__ruvK", spec$ruvK[[hit]],
    "__kPB", spec$kPB[[hit]],
    "__full_v1.rds"
  )
)

meta_file <- file.path(
  in_dir,
  paste0(compartment, "__full_v1_cells.tsv.gz")
)

out_dir <- file.path(
  root, "atlas", "full_atlas_annotation_v1",
  "integrated_compartment_clustering", compartment
)
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

if(!file.exists(in_rds)) stop("Missing input RDS: ", in_rds)
if(!file.exists(meta_file)) stop("Missing cell metadata: ", meta_file)

out_rds <- file.path(out_dir, "integrated_clustering_v1.rds")
graph_rds <- file.path(out_dir, "snn_graph_v1.rds")
complete_file <- file.path(out_dir, "CLUSTERING_COMPLETE.ok")

if(file.exists(complete_file) && file.exists(out_rds) && file.exists(graph_rds)){
  cat("SKIP: completed output already exists for ", compartment, "\n", sep="")
  quit(save="no", status=0)
}

z <- readRDS(in_rds)

pca <- z$pca_embeddings
if(is.null(dim(pca))) stop("Missing PCA embeddings")
if(ncol(pca) < 30L) stop("Expected at least 30 PCs")
pca <- pca[,1:30,drop=FALSE]

cells <- rownames(pca)
if(is.null(cells)) stop("PCA embeddings have no cell names")
if(length(cells) != length(unique(cells))) stop("Duplicate PCA cell names")
if(!all(is.finite(pca))) stop("Non-finite PCA values")

meta <- read.delim(
  gzfile(meta_file),
  sep="\t",
  header=TRUE,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

if(!"global_cell" %in% names(meta)) stop("global_cell missing from metadata")
if(anyDuplicated(meta$global_cell)) stop("Duplicate global_cell in metadata")
if(!setequal(cells, meta$global_cell)) stop("PCA/metadata cell mismatch")

meta <- meta[match(cells, meta$global_cell),,drop=FALSE]
stopifnot(identical(meta$global_cell, cells))

needed_meta <- c(
  "project_id",
  "condition_binary",
  "integration_celltype_full_v1",
  "transfer_confidence"
)
missing_meta <- setdiff(needed_meta, names(meta))
if(length(missing_meta)){
  stop("Missing metadata columns: ", paste(missing_meta, collapse=", "))
}

fc <- get("FindClusters.default", envir=asNamespace("Seurat"), inherits=FALSE)
fc_args <- names(formals(fc))
need_fc <- c("leiden_method","leiden_objective_function")
if(!all(need_fc %in% fc_args)){
  stop(
    "Installed Seurat FindClusters API lacks required igraph Leiden arguments. ",
    "Seurat version=", as.character(packageVersion("Seurat"))
  )
}

params <- list(
  dims=1:30,
  k_param=30L,
  nn_method="annoy",
  n_trees=100L,
  annoy_metric="euclidean",
  prune_SNN=1/15,
  resolutions=c(0.2,0.4,0.6,0.8,1.0,1.2),
  clustering_algorithm=4L,
  leiden_method="igraph",
  leiden_objective_function="modularity",
  n_iter=20L,
  umap_neighbors=30L,
  umap_metric="cosine",
  umap_min_dist=0.3,
  seed=20260813L
)

cat("=== INTEGRATED COMPARTMENT CLUSTERING ===\n")
cat("compartment=", compartment, "\n", sep="")
cat("cells=", nrow(pca), "; PCs=", ncol(pca), "\n", sep="")
cat("Seurat=", as.character(packageVersion("Seurat")), "\n", sep="")
cat("igraph=", as.character(packageVersion("igraph")), "\n", sep="")

cat("\n--- LOAD/BUILD KNN/SNN ---\n")

if(file.exists(graph_rds)){
  cat("Reusing existing SNN checkpoint: ", graph_rds, "\n", sep="")
  g <- readRDS(graph_rds)

  if(is.null(g$snn)) stop("Existing graph checkpoint has no SNN")
  if(!identical(g$cells, cells)) stop("Existing graph checkpoint cell mismatch")

  snn <- g$snn
  rm(g)
  gc(verbose=FALSE)

} else {

  ng <- Seurat::FindNeighbors(
    object=pca,
    k.param=params$k_param,
    return.neighbor=FALSE,
    compute.SNN=TRUE,
    prune.SNN=params$prune_SNN,
    nn.method=params$nn_method,
    n.trees=params$n_trees,
    annoy.metric=params$annoy_metric,
    verbose=TRUE
  )

  if(is.null(ng$snn)) stop("SNN graph was not returned")
  snn <- ng$snn
  rm(ng)
  gc(verbose=FALSE)

  saveRDS(
    list(
      compartment=compartment,
      cells=cells,
      snn=snn,
      params=params,
      source_integration_rds=in_rds,
      Seurat_version=as.character(packageVersion("Seurat")),
      igraph_version=as.character(packageVersion("igraph"))
    ),
    graph_rds,
    compress=FALSE
  )
}

if(!identical(rownames(snn), cells) || !identical(colnames(snn), cells)){
  stop("SNN cell order mismatch")
}

cat("\n--- LEIDEN RESOLUTION SWEEP ---\n")
cl <- Seurat::FindClusters(
  object=snn,
  resolution=params$resolutions,
  algorithm=params$clustering_algorithm,
  leiden_method=params$leiden_method,
  leiden_objective_function=params$leiden_objective_function,
  n.iter=params$n_iter,
  random.seed=params$seed,
  group.singletons=TRUE,
  verbose=TRUE
)

if(!identical(rownames(cl), cells)) stop("Cluster output cell order mismatch")

resolution_from_name <- function(x){
  as.numeric(sub("^res\\.", "", x))
}

cluster_col_name <- function(r){
  paste0("cluster_res_", gsub("\\.", "p", format(r, trim=TRUE, scientific=FALSE)))
}

old_names <- names(cl)
resolutions_observed <- vapply(old_names, resolution_from_name, numeric(1))

new_names <- vapply(resolutions_observed, cluster_col_name, character(1))
names(cl) <- new_names

cat("\n--- UMAP ---\n")
umap <- Seurat::RunUMAP(
  object=pca,
  n.neighbors=params$umap_neighbors,
  n.components=2L,
  metric=params$umap_metric,
  min.dist=params$umap_min_dist,
  spread=1,
  seed.use=params$seed,
  verbose=TRUE
)

if(inherits(umap, "DimReduc")){
  umap <- SeuratObject::Embeddings(umap)
} else {
  umap <- as.matrix(umap)
}

umap <- as.matrix(umap)

if(nrow(umap) != length(cells)) stop("UMAP row count mismatch")
if(is.null(rownames(umap))) rownames(umap) <- cells
if(!identical(rownames(umap), cells)) umap <- umap[cells,,drop=FALSE]
colnames(umap) <- c("UMAP_1","UMAP_2")
if(!all(is.finite(umap))) stop("Non-finite UMAP values")

entropy_norm <- function(x){
  x <- x[x > 0]
  if(length(x) <= 1L) return(0)
  p <- x/sum(x)
  -sum(p*log(p))/log(length(p))
}

diagnostics <- list()
res_summary <- list()

for(cc in names(cl)){
  labels <- as.character(cl[[cc]])
  r <- resolutions_observed[match(cc,new_names)]

  tab_size <- table(labels)
  clusters <- names(tab_size)

  rows <- lapply(clusters, function(k){
    idx <- which(labels == k)
    n <- length(idx)

    proj_tab <- table(meta$project_id[idx])
    ct_tab <- table(meta$integration_celltype_full_v1[idx])
    cond_tab <- table(meta$condition_binary[idx])
    conf_tab <- table(meta$transfer_confidence[idx])

    max_proj <- max(proj_tab)/n
    top_ct <- names(ct_tab)[which.max(ct_tab)]
    top_ct_frac <- max(ct_tab)/n

    healthy_frac <- if("healthy" %in% names(cond_tab)) unname(cond_tab[["healthy"]])/n else 0
    disease_frac <- if("disease" %in% names(cond_tab)) unname(cond_tab[["disease"]])/n else 0
    low_frac <- if("low" %in% names(conf_tab)) unname(conf_tab[["low"]])/n else 0

    data.frame(
      compartment=compartment,
      resolution=r,
      cluster=k,
      n_cells=n,
      n_projects=length(proj_tab),
      max_project_fraction=max_proj,
      project_entropy_normalized=entropy_norm(as.numeric(proj_tab)),
      top_transfer_celltype=top_ct,
      transfer_celltype_purity=top_ct_frac,
      healthy_fraction=healthy_frac,
      disease_fraction=disease_frac,
      low_transfer_fraction=low_frac,
      stringsAsFactors=FALSE
    )
  })

  d <- do.call(rbind, rows)
  diagnostics[[cc]] <- d

  res_summary[[cc]] <- data.frame(
    compartment=compartment,
    resolution=r,
    n_clusters=nrow(d),
    min_cluster_size=min(d$n_cells),
    median_cluster_size=median(d$n_cells),
    max_cluster_size=max(d$n_cells),
    n_clusters_lt_100=sum(d$n_cells < 100),
    n_clusters_lt_200=sum(d$n_cells < 200),
    median_max_project_fraction=median(d$max_project_fraction),
    p95_max_project_fraction=as.numeric(
      quantile(d$max_project_fraction, probs=0.95, names=FALSE)
    ),
    median_project_entropy=median(d$project_entropy_normalized),
    median_transfer_celltype_purity=median(d$transfer_celltype_purity),
    stringsAsFactors=FALSE
  )
}

diagnostics <- do.call(rbind, diagnostics)
res_summary <- do.call(rbind, res_summary)

assignment <- data.frame(
  global_cell=cells,
  cl,
  stringsAsFactors=FALSE,
  check.names=FALSE
)

umap_df <- data.frame(
  global_cell=cells,
  UMAP_1=umap[,1],
  UMAP_2=umap[,2],
  stringsAsFactors=FALSE
)

{
  con <- gzfile(
    file.path(out_dir,"cluster_assignments_v1.tsv.gz"),
    open="wt"
  )
  tryCatch(
    write.table(
      assignment,
      con,
      sep="\t",
      quote=1,
      qmethod="double",
      row.names=FALSE
    ),
    finally=close(con)
  )
}

{
  con <- gzfile(
    file.path(out_dir,"umap_v1.tsv.gz"),
    open="wt"
  )
  tryCatch(
    write.table(
      umap_df,
      con,
      sep="\t",
      quote=1,
      qmethod="double",
      row.names=FALSE
    ),
    finally=close(con)
  )
}

write.table(
  diagnostics,
  file=file.path(out_dir,"cluster_diagnostics_v1.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

write.table(
  res_summary,
  file=file.path(out_dir,"resolution_summary_v1.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

saveRDS(
  list(
    compartment=compartment,
    cells=cells,
    clusters=cl,
    umap=umap,
    params=params,
    resolution_summary=res_summary,
    source_integration_rds=in_rds,
    source_cell_metadata=meta_file,
    Seurat_version=as.character(packageVersion("Seurat")),
    igraph_version=as.character(packageVersion("igraph"))
  ),
  out_rds,
  compress=FALSE
)

writeLines(
  c(
    paste0("compartment=",compartment),
    paste0("n_cells=",length(cells)),
    paste0("n_resolutions=",length(params$resolutions)),
    paste0("resolutions=",paste(params$resolutions,collapse=",")),
    paste0("k=",params$k_param),
    paste0("algorithm=Leiden_igraph_modularity"),
    paste0("output=",out_rds)
  ),
  complete_file
)

cat("\n=== RESOLUTION SUMMARY ===\n")
print(res_summary,row.names=FALSE)
cat("\nPASS: integrated compartment clustering completed for ", compartment, "\n", sep="")
