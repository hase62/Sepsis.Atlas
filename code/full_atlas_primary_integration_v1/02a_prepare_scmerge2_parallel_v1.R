#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3L) {
  stop("Usage: 02a_prepare_scmerge2_parallel_v1.R ROOT COMPARTMENT CHUNK_SIZE")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
compartment <- as.character(args[[2]])
chunk_size <- as.integer(args[[3]])
if (!is.finite(chunk_size) || chunk_size < 1000L) stop("Invalid CHUNK_SIZE")

source(file.path(root,"full_atlas_primary_integration_v1","00_common.R"))

required <- c(
  "scMerge","BiocParallel","BiocSingular","BiocNeighbors",
  "irlba","batchelor","Matrix"
)
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) stop("Missing packages: ",paste(missing,collapse=", "))

suppressPackageStartupMessages({
  library(BiocParallel)
  library(BiocSingular)
  library(BiocNeighbors)
  library(irlba)
})

set.seed(20260811)

base <- file.path(root,"pre_integration")
pilot_dir <- file.path(base,"pilot_unintegrated_large_v1")
full_dir <- file.path(base,"full_atlas_primary_v1")
transfer_dir <- file.path(full_dir,"annotation_transfer_by_library")

parallel_root <- file.path(full_dir,"scmerge2_parallel_v1")
work_dir <- file.path(parallel_root,compartment)
prepared_dir <- file.path(work_dir,"prepared_chunks")
adjusted_dir <- file.path(work_dir,"adjusted_chunks")

dir.create(prepared_dir,recursive=TRUE,showWarnings=FALSE)
dir.create(adjusted_dir,recursive=TRUE,showWarnings=FALSE)

plan <- data.frame(
  compartment=c("T_NK_combined","Monocyte_DC_combined"),
  ruvK=c(5L,2L),
  k_pseudoBulk=c(5L,5L),
  stringsAsFactors=FALSE
)

jj <- match(compartment,plan$compartment)
if (is.na(jj)) stop("Unknown compartment: ",compartment)

ruvK <- plan$ruvK[[jj]]
kpb <- plan$k_pseudoBulk[[jj]]
seed <- 20260811L + jj

model_path <- file.path(work_dir,"model_v1.rds")
manifest_path <- file.path(work_dir,"chunk_manifest_v1.tsv")
ready_path <- file.path(work_dir,"PREPARED.ok")

atomic_save_rds <- function(x,path,compress=FALSE){
  if (file.exists(path)) return(invisible(FALSE))
  tmp <- paste0(path,".tmp.",Sys.getpid())
  saveRDS(x,tmp,compress=compress)
  if (!file.rename(tmp,path)) {
    unlink(tmp)
    stop("Atomic rename failed: ",path)
  }
  invisible(TRUE)
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

read_transfer <- function(project_id,library_key){
  tag <- safe_name(paste(project_id,library_key,sep="__"))
  f <- file.path(transfer_dir,paste0(tag,".tsv.gz"))
  if (!file.exists(f)) stop("Missing transfer: ",f)
  read_tsv(f)
}

cat("=== PREPARE PARALLEL SCMERGE2 ===\n")
cat("compartment=",compartment,
    "; ruvK=",ruvK,
    "; kPB=",kpb,
    "; chunk_size=",chunk_size,"\n",sep="")

lib <- read_tsv(file.path(pilot_dir,"resolved_library_table.tsv"))
hvg <- read_tsv(file.path(pilot_dir,"large_pilot_consensus_hvg_3000.tsv"))
consensus_hvg <- unique(as.character(hvg$feature))
if (length(consensus_hvg)!=3000L) stop("Expected 3000 consensus HVGs")

first_obj <- readRDS(as.character(lib$final_rds_resolved[[1]]))
feature_reference <- rownames(get_rna_counts(first_obj))
rm(first_obj)
gc(verbose=FALSE)

data("segList",package="scMerge",envir=environment())
human_seg <- intersect(
  unique(as.character(segList$human$human_scSEG)),
  feature_reference
)
if (length(human_seg)<100L) stop("Too few human SEG controls")

input_features <- unique(c(consensus_hvg,human_seg))

counts_list <- list()
meta_list <- list()

for(i in seq_len(nrow(lib))){
  project_id <- as.character(lib$project_id[[i]])
  library_key <- as.character(lib$library_key[[i]])

  tr <- read_transfer(project_id,library_key)
  tr <- tr[
    tr$integration_compartment_primary_v1==compartment,
    ,
    drop=FALSE
  ]
  if(!nrow(tr)) next

  obj <- readRDS(as.character(lib$final_rds_resolved[[i]]))
  counts <- get_rna_counts(obj)

  miss <- setdiff(input_features,rownames(counts))
  if(length(miss)) {
    stop(compartment,": missing features in ",library_key)
  }

  if(!all(tr$original_cell %in% colnames(counts))) {
    stop(compartment,": missing cells in ",library_key)
  }

  m <- counts[
    input_features,
    tr$original_cell,
    drop=FALSE
  ]
  colnames(m) <- tr$global_cell

  counts_list[[length(counts_list)+1L]] <- m
  meta_list[[length(meta_list)+1L]] <- tr[
    ,
    c(
      "global_cell","project_id","library_key",
      "clinical_domain_std","condition_binary",
      "integration_celltype_full_v1","transfer_confidence"
    ),
    drop=FALSE
  ]

  rm(obj,counts,m,tr)
  gc(verbose=FALSE)
}

counts_sub <- do.call(cbind,counts_list)
rm(counts_list)
gc(verbose=FALSE)

meta <- do.call(rbind,meta_list)
rm(meta_list)

rownames(meta) <- meta$global_cell
cells <- colnames(counts_sub)

if(!setequal(cells,rownames(meta))) {
  stop(compartment,": metadata mismatch")
}
meta <- meta[cells,,drop=FALSE]

expressed <- Matrix::rowSums(counts_sub)>0
counts_sub <- counts_sub[expressed,,drop=FALSE]

data_sub <- log_normalize_sparse(counts_sub)

chosen_hvg <- intersect(consensus_hvg,rownames(data_sub))
ctl <- intersect(human_seg,rownames(data_sub))

if(length(chosen_hvg)<2500L) stop("Too few HVGs")
if(length(ctl)<100L) stop("Too few controls")

batch <- as.character(meta$project_id)
cell_types <- as.character(meta$integration_celltype_full_v1)
condition <- as.character(meta$condition_binary)

if(length(unique(batch))<2L || length(unique(condition))<2L) {
  stop("Insufficient batch/condition levels")
}

cat(
  "Cells=",length(cells),
  "; projects=",length(unique(batch)),
  "; cellTypes=",length(unique(cell_types)),
  "; features=",nrow(data_sub),
  "; RSS=",sprintf("%.2f",rss_gb())," GB\n",
  sep=""
)

cat("Computing native baseline PCA...\n")
baseline_pca <- run_pca(
  data_sub[chosen_hvg,cells,drop=FALSE],
  30L
)
rownames(baseline_pca$embeddings) <- cells

cat("Estimating scMerge2 model with return_matrix=FALSE...\n")
start_model <- proc.time()[["elapsed"]]

result <- scMerge::scMerge2(
  exprsMat=data_sub,
  batch=batch,
  cellTypes=cell_types,
  condition=condition,
  ctl=ctl,
  chosen.hvg=chosen_hvg,
  ruvK=ruvK,
  use_bpparam=BiocParallel::SerialParam(),
  use_bsparam=BiocSingular::RandomParam(),
  use_bnparam=BiocNeighbors::AnnoyParam(),
  pseudoBulk_fn="create_pseudoBulk",
  k_pseudoBulk=kpb,
  k_celltype=10,
  exprsMat_counts=counts_sub,
  cosineNorm=TRUE,
  return_subset=TRUE,
  return_subset_genes=chosen_hvg,
  return_matrix=FALSE,
  byChunk=TRUE,
  chunkSize=5000,
  verbose=TRUE,
  seed=seed
)

model_elapsed <- proc.time()[["elapsed"]] - start_model

if(is.null(result$fullalpha)) stop("scMerge2 returned no fullalpha")

cat(
  "Model estimation complete; elapsed=",
  sprintf("%.1f",model_elapsed),
  " sec; RSS=",sprintf("%.2f",rss_gb())," GB\n",
  sep=""
)

cat("Cosine normalising matrix for independent adjustment chunks...\n")
cosine_mat <- batchelor::cosineNorm(data_sub)

if (!identical(colnames(cosine_mat),cells)) {
  cosine_mat <- cosine_mat[,cells,drop=FALSE]
}

adjusted_means <- if (inherits(cosine_mat,"Matrix")) {
  Matrix::rowMeans(cosine_mat)
} else {
  rowMeans(cosine_mat)
}
names(adjusted_means) <- rownames(cosine_mat)

model <- list(
  compartment=compartment,
  ruvK=ruvK,
  k_pseudoBulk=kpb,
  seed=seed,
  chunk_size=chunk_size,
  cells=cells,
  meta=meta,
  batch=batch,
  cell_types=cell_types,
  condition=condition,
  chosen_hvg=chosen_hvg,
  controls=ctl,
  input_features=rownames(cosine_mat),
  adjusted_means=adjusted_means,
  fullalpha=result$fullalpha,
  M=result$M,
  baseline_pca=baseline_pca,
  model_elapsed_seconds=model_elapsed,
  scMerge_version=as.character(packageVersion("scMerge"))
)

if (file.exists(model_path)) {
  old <- readRDS(model_path)
  if (!identical(old$cells,model$cells) ||
      !identical(old$chosen_hvg,model$chosen_hvg) ||
      old$ruvK != model$ruvK ||
      old$k_pseudoBulk != model$k_pseudoBulk) {
    stop("Existing model_v1.rds is incompatible; refusing overwrite")
  }
  cat("Existing compatible model retained: ",model_path,"\n",sep="")
} else {
  atomic_save_rds(model,model_path,compress=FALSE)
  cat("Saved model: ",model_path,"\n",sep="")
}

chunk_id <- ceiling(seq_along(cells)/chunk_size)
chunk_levels <- sort(unique(chunk_id))

manifest <- do.call(
  rbind,
  lapply(chunk_levels,function(k){
    idx <- which(chunk_id==k)
    data.frame(
      chunk_id=as.integer(k),
      start_index=min(idx),
      end_index=max(idx),
      n_cells=length(idx),
      prepared_file=sprintf("prepared_chunk_%04d.rds",k),
      adjusted_file=sprintf("adjusted_chunk_%04d.rds",k),
      stringsAsFactors=FALSE
    )
  })
)

if (file.exists(manifest_path)) {
  old_manifest <- read_tsv(manifest_path)
  if (!identical(
    old_manifest[,names(manifest),drop=FALSE],
    manifest
  )) {
    stop("Existing chunk manifest differs; refusing overwrite")
  }
} else {
  write_tsv(manifest,manifest_path)
}

cat("Preparing ",nrow(manifest)," independent chunks...\n",sep="")

for (ii in seq_len(nrow(manifest))) {
  k <- manifest$chunk_id[[ii]]
  a <- manifest$start_index[[ii]]
  b <- manifest$end_index[[ii]]
  idx <- a:b

  f <- file.path(prepared_dir,manifest$prepared_file[[ii]])

  if (file.exists(f)) {
    z <- readRDS(f)
    if (!identical(z$cells,cells[idx])) {
      stop("Existing prepared chunk has wrong cells: ",f)
    }
    cat("SKIP existing chunk ",k,"\n",sep="")
    next
  }

  chunk <- cosine_mat[,idx,drop=FALSE]

  z <- list(
    chunk_id=k,
    cells=cells[idx],
    mat=chunk
  )

  atomic_save_rds(z,f,compress=FALSE)
  rm(z,chunk)
  gc(verbose=FALSE)

  cat(
    "PREPARED chunk ",k,
    "/",nrow(manifest),
    " cells=",length(idx),
    " RSS=",sprintf("%.2f",rss_gb())," GB\n",
    sep=""
  )
}

writeLines(
  c(
    paste0("compartment=",compartment),
    paste0("n_cells=",length(cells)),
    paste0("n_chunks=",nrow(manifest)),
    paste0("chunk_size=",chunk_size),
    paste0("ruvK=",ruvK),
    paste0("k_pseudoBulk=",kpb)
  ),
  ready_path
)

cat("PASS: model + prepared chunks complete for ",compartment,"\n",sep="")
