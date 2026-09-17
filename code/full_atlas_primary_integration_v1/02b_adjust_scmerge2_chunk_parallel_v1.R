#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 3L) {
  stop("Usage: 02b_adjust_scmerge2_chunk_parallel_v1.R ROOT COMPARTMENT TASK_ID")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
compartment <- as.character(args[[2]])
task_id <- as.integer(args[[3]])

source(file.path(root,"full_atlas_primary_integration_v1","00_common.R"))

if (!requireNamespace("scMerge",quietly=TRUE)) {
  stop("Missing scMerge")
}

full_dir <- file.path(root,"pre_integration","full_atlas_primary_v1")
work_dir <- file.path(
  full_dir,
  "scmerge2_parallel_v1",
  compartment
)
prepared_dir <- file.path(work_dir,"prepared_chunks")
adjusted_dir <- file.path(work_dir,"adjusted_chunks")

model_path <- file.path(work_dir,"model_v1.rds")
manifest_path <- file.path(work_dir,"chunk_manifest_v1.tsv")

if (!file.exists(model_path)) stop("Missing model: ",model_path)
if (!file.exists(manifest_path)) stop("Missing manifest: ",manifest_path)

model <- readRDS(model_path)
manifest <- read_tsv(manifest_path)

hit <- which(manifest$chunk_id==task_id)
if(length(hit)!=1L) stop("Invalid TASK_ID: ",task_id)

row <- manifest[hit,,drop=FALSE]

prepared_file <- file.path(prepared_dir,row$prepared_file[[1]])
adjusted_file <- file.path(adjusted_dir,row$adjusted_file[[1]])

if (!file.exists(prepared_file)) {
  stop("Missing prepared chunk: ",prepared_file)
}

if (file.exists(adjusted_file)) {
  z <- readRDS(adjusted_file)
  expected <- model$cells[row$start_index[[1]]:row$end_index[[1]]]
  if (!identical(z$cells,expected)) {
    stop("Existing adjusted chunk incompatible: ",adjusted_file)
  }
  cat("SKIP existing adjusted chunk: ",adjusted_file,"\n",sep="")
  quit(save="no",status=0)
}

chunk <- readRDS(prepared_file)

expected_cells <- model$cells[
  row$start_index[[1]]:row$end_index[[1]]
]

if (!identical(chunk$cells,expected_cells)) {
  stop("Prepared chunk cell order mismatch")
}

if (!identical(rownames(chunk$mat),model$input_features)) {
  stop("Prepared chunk feature order mismatch")
}

cat(
  "=== ADJUST CHUNK ===\n",
  "compartment=",compartment,
  "; task=",task_id,
  "; cells=",length(chunk$cells),
  "; ruvK=",model$ruvK,
  "; RSS=",sprintf("%.2f",rss_gb())," GB\n",
  sep=""
)

start <- proc.time()[["elapsed"]]

adjusted <- scMerge::getAdjustedMat(
  exprsMat=chunk$mat,
  fullalpha=model$fullalpha,
  ctl=model$controls,
  adjusted_means=model$adjusted_means,
  ruvK=model$ruvK,
  return_subset_genes=model$chosen_hvg
)

elapsed <- proc.time()[["elapsed"]] - start

if (is.null(rownames(adjusted))) {
  rownames(adjusted) <- model$chosen_hvg
}
if (is.null(colnames(adjusted))) {
  colnames(adjusted) <- chunk$cells
}

if (!identical(rownames(adjusted),model$chosen_hvg) ||
    !identical(colnames(adjusted),chunk$cells)) {
  adjusted <- adjusted[
    model$chosen_hvg,
    chunk$cells,
    drop=FALSE
  ]
}

# getAdjustedMat() returns a DelayedMatrix; materialize before checkpointing.
adjusted <- as.matrix(adjusted)

expected_dim <- c(length(model$chosen_hvg),length(chunk$cells))
if(!identical(dim(adjusted),expected_dim)) {
  stop(
    "Adjusted matrix dimension mismatch; observed=",
    paste(dim(adjusted),collapse="x"),
    "; expected=",paste(expected_dim,collapse="x")
  )
}

out <- list(
  compartment=compartment,
  chunk_id=task_id,
  cells=chunk$cells,
  adjusted=adjusted,
  elapsed_seconds=elapsed,
  rss_gb_after=rss_gb()
)

tmp <- paste0(adjusted_file,".tmp.",Sys.getpid())
saveRDS(out,tmp,compress=FALSE)

if (!file.rename(tmp,adjusted_file)) {
  unlink(tmp)
  stop("Atomic rename failed: ",adjusted_file)
}

cat(
  "PASS chunk ",task_id,
  "; elapsed=",sprintf("%.1f",elapsed)," sec",
  "; RSS=",sprintf("%.2f",rss_gb())," GB\n",
  sep=""
)
