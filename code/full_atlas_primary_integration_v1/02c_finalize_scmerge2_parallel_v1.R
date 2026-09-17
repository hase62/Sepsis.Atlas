#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2L) {
  stop("Usage: 02c_finalize_scmerge2_parallel_v1.R ROOT COMPARTMENT")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
compartment <- as.character(args[[2]])

source(file.path(root,"full_atlas_primary_integration_v1","00_common.R"))

required <- c("irlba","scMerge")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing)) stop("Missing packages: ",paste(missing,collapse=", "))

full_dir <- file.path(root,"pre_integration","full_atlas_primary_v1")
run_dir <- file.path(full_dir,"scmerge2_primary_runs")
work_dir <- file.path(
  full_dir,
  "scmerge2_parallel_v1",
  compartment
)
adjusted_dir <- file.path(work_dir,"adjusted_chunks")

dir.create(run_dir,recursive=TRUE,showWarnings=FALSE)

model <- readRDS(file.path(work_dir,"model_v1.rds"))
manifest <- read_tsv(file.path(work_dir,"chunk_manifest_v1.tsv"))

out_rds <- file.path(
  run_dir,
  paste0(
    compartment,
    "__ruvK",model$ruvK,
    "__kPB",model$k_pseudoBulk,
    "__full_v1.rds"
  )
)

if(file.exists(out_rds)){
  cat("SKIP existing final output: ",out_rds,"\n",sep="")
  quit(save="no",status=0)
}

run_pca <- function(mat,n_pcs=30L){
  n_pcs <- min(n_pcs,nrow(mat)-1L,ncol(mat)-1L)
  fit <- irlba::prcomp_irlba(
    t(mat),
    n=n_pcs,
    center=TRUE,
    scale.=FALSE
  )
  v <- fit$sdev^2
  v <- v/sum(v)
  list(embeddings=fit$x,variance=v)
}

eta_squared <- function(score,group){
  keep <- is.finite(score)&!is.na(group)
  score <- score[keep]
  group <- droplevels(factor(group[keep]))
  if(length(score)<3L||nlevels(group)<2L) return(NA_real_)
  total <- sum((score-mean(score))^2)
  if(total<=0) return(0)
  means <- tapply(score,group,mean)
  counts <- table(group)
  sum(counts*(means-mean(score))^2)/total
}

weighted_eta <- function(pca,variance,group,n_pc=20L){
  n <- min(n_pc,ncol(pca),length(variance))
  eta <- vapply(
    seq_len(n),
    function(i) eta_squared(pca[,i],group),
    numeric(1)
  )
  w <- variance[seq_len(n)]
  keep <- is.finite(eta)&is.finite(w)
  if(!any(keep)) return(NA_real_)
  sum(eta[keep]*w[keep])/sum(w[keep])
}

distance_preservation <- function(base_pca,adj_pca,max_cells=1500L){
  common <- intersect(rownames(base_pca),rownames(adj_pca))
  if(length(common)<10L) return(NA_real_)
  set.seed(20260811)
  use <- if(length(common)<=max_cells) common else sample(common,max_cells)
  suppressWarnings(
    cor(
      as.numeric(dist(base_pca[use,,drop=FALSE])),
      as.numeric(dist(adj_pca[use,,drop=FALSE])),
      method="spearman"
    )
  )
}

missing_chunks <- character()

for(i in seq_len(nrow(manifest))){
  f <- file.path(adjusted_dir,manifest$adjusted_file[[i]])
  if(!file.exists(f)) missing_chunks <- c(missing_chunks,f)
}

if(length(missing_chunks)){
  stop(
    "Missing ",length(missing_chunks),
    " adjusted chunks; first: ",missing_chunks[[1]]
  )
}

n_gene <- length(model$chosen_hvg)
n_cell <- length(model$cells)

cat(
  "Allocating final adjusted matrix: ",
  n_gene," x ",n_cell,
  "; RSS=",sprintf("%.2f",rss_gb())," GB\n",
  sep=""
)

adjusted <- matrix(
  NA_real_,
  nrow=n_gene,
  ncol=n_cell,
  dimnames=list(model$chosen_hvg,model$cells)
)

chunk_elapsed <- numeric(nrow(manifest))

for(i in seq_len(nrow(manifest))){
  f <- file.path(adjusted_dir,manifest$adjusted_file[[i]])
  z <- readRDS(f)

  a <- manifest$start_index[[i]]
  b <- manifest$end_index[[i]]
  idx <- a:b

  if(!identical(z$cells,model$cells[idx])) {
    stop("Adjusted chunk cell order mismatch: ",f)
  }

  block <- as.matrix(z$adjusted)

  expected_dim <- c(length(model$chosen_hvg),length(idx))
  if(!identical(dim(block),expected_dim)) {
    stop(
      "Adjusted chunk dimension mismatch: ",f,
      "; observed=",paste(dim(block),collapse="x"),
      "; expected=",paste(expected_dim,collapse="x")
    )
  }

  if(!identical(rownames(block),model$chosen_hvg)) {
    stop("Adjusted chunk feature mismatch: ",f)
  }

  if(!identical(colnames(block),model$cells[idx])) {
    stop("Adjusted chunk column/cell mismatch: ",f)
  }

  adjusted[,idx] <- block
  chunk_elapsed[[i]] <- z$elapsed_seconds

  rm(z,block)
  gc(verbose=FALSE)

  cat(
    "LOADED ",i,"/",nrow(manifest),
    "; RSS=",sprintf("%.2f",rss_gb())," GB\n",
    sep=""
  )
}

if(anyNA(adjusted)) stop("NA detected in assembled adjusted matrix")

cat("Running adjusted PCA...\n")
pca_adj <- run_pca(adjusted,30L)
rownames(pca_adj$embeddings) <- model$cells

base_pca <- model$baseline_pca
batch <- model$batch
condition <- model$condition
cell_types <- model$cell_types

base_proj <- weighted_eta(
  base_pca$embeddings,
  base_pca$variance,
  batch
)
adj_proj <- weighted_eta(
  pca_adj$embeddings,
  pca_adj$variance,
  batch
)

base_cond <- weighted_eta(
  base_pca$embeddings,
  base_pca$variance,
  condition
)
adj_cond <- weighted_eta(
  pca_adj$embeddings,
  pca_adj$variance,
  condition
)

base_ct <- weighted_eta(
  base_pca$embeddings,
  base_pca$variance,
  cell_types
)
adj_ct <- weighted_eta(
  pca_adj$embeddings,
  pca_adj$variance,
  cell_types
)

dist_cor <- distance_preservation(
  base_pca$embeddings,
  pca_adj$embeddings
)

total_scmerge_seconds <-
  model$model_elapsed_seconds +
  sum(chunk_elapsed,na.rm=TRUE)

met <- data.frame(
  compartment=compartment,
  n_cells=length(model$cells),
  n_projects=length(unique(batch)),
  n_celltypes=length(unique(cell_types)),
  ruvK=model$ruvK,
  k_pseudoBulk=model$k_pseudoBulk,
  baseline_project_eta2=base_proj,
  adjusted_project_eta2=adj_proj,
  project_eta2_reduction_fraction=1-adj_proj/base_proj,
  baseline_condition_eta2=base_cond,
  adjusted_condition_eta2=adj_cond,
  baseline_celltype_eta2=base_ct,
  adjusted_celltype_eta2=adj_ct,
  distance_spearman_native_vs_adjusted=dist_cor,
  elapsed_seconds=total_scmerge_seconds,
  rss_gb_after_run=rss_gb(),
  stringsAsFactors=FALSE
)

final_object <- list(
  compartment=compartment,
  ruvK=model$ruvK,
  k_pseudoBulk=model$k_pseudoBulk,
  cells=model$cells,
  pca_embeddings=pca_adj$embeddings,
  pca_variance=pca_adj$variance,
  fullalpha=model$fullalpha,
  M=model$M,
  controls=model$controls,
  chosen_hvg=model$chosen_hvg,
  metrics=met,
  scMerge_version=model$scMerge_version,
  note=paste(
    "Adjusted expression matrix intentionally not retained;",
    "native SoupX-corrected RNA remains basis for markers/DE.",
    "Adjustment generated by checkpointed parallel getAdjustedMat workflow."
  )
)

tmp <- paste0(out_rds,".tmp.",Sys.getpid())
saveRDS(final_object,tmp,compress=FALSE)

if(!file.rename(tmp,out_rds)){
  unlink(tmp)
  stop("Atomic rename failed: ",out_rds)
}

cells_path <- file.path(
  run_dir,
  paste0(compartment,"__full_v1_cells.tsv.gz")
)

if(!file.exists(cells_path)){
  write_tsv(model$meta,cells_path)
}

write_tsv(
  met,
  file.path(work_dir,"final_metrics_v1.tsv")
)

writeLines(
  c(
    paste0("compartment=",compartment),
    paste0("n_cells=",length(model$cells)),
    paste0("n_chunks=",nrow(manifest)),
    paste0("output=",out_rds)
  ),
  file.path(work_dir,"FINALIZED.ok")
)

rm(adjusted,pca_adj,base_pca)
gc(verbose=FALSE)

cat("\n=== FINAL METRICS ===\n")
print(met,row.names=FALSE)
cat("\nPASS: finalized ",compartment,"\n",sep="")
