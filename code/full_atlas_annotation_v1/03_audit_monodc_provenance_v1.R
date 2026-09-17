#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".", mustWork=TRUE)

suppressPackageStartupMessages({
  library(SeuratObject)
})

base <- file.path(root, "atlas", "full_atlas_annotation_v1", "integrated_compartment_clustering")
comp <- "Monocyte_DC_combined"

cl_rds <- file.path(base, comp, "integrated_clustering_v1.rds")
meta_file <- file.path(root, "pre_integration", "full_atlas_primary_v1", "scmerge2_primary_runs", "Monocyte_DC_combined__full_v1_cells.tsv.gz")
lib_file <- file.path(root, "pre_integration", "pilot_unintegrated_large_v1", "resolved_library_table.tsv")

out_dir <- file.path(root, "atlas", "full_atlas_annotation_v1", "monodc_provenance_audit_v1")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

stopifnot(file.exists(cl_rds), file.exists(meta_file), file.exists(lib_file))

x <- readRDS(cl_rds)
cells <- as.character(x$cells)
cl <- x$clusters

stopifnot(length(cells) == nrow(cl), identical(rownames(cl), cells))

meta <- read.delim(gzfile(meta_file), sep="\t", quote="\"", comment.char="", stringsAsFactors=FALSE, check.names=FALSE)
required_meta <- c("global_cell","project_id","library_key","condition_binary")
miss <- setdiff(required_meta, names(meta))
if(length(miss)) stop("Missing full-cell metadata columns: ", paste(miss, collapse=", "))

m <- match(cells, as.character(meta$global_cell))
if(anyNA(m)) stop("Full-cell metadata mismatch: ", sum(is.na(m)), " cells")
meta <- meta[m,,drop=FALSE]
stopifnot(identical(as.character(meta$global_cell), cells))

lib <- read.delim(lib_file, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
required_lib <- c("project_id","library_key","final_rds_resolved")
miss <- setdiff(required_lib, names(lib))
if(length(miss)) stop("Missing resolved-library columns: ", paste(miss, collapse=", "))

used_libs <- unique(as.character(meta$library_key))
idx <- match(used_libs, as.character(lib$library_key))
if(anyNA(idx)) stop("Libraries missing from resolved_library_table.tsv: ", paste(used_libs[is.na(idx)], collapse=", "))
lib <- lib[idx,,drop=FALSE]

first_nonempty <- function(v){
  v <- as.character(v)
  v <- unique(v[!is.na(v) & nzchar(v)])
  if(!length(v)) return(NA_character_)
  paste(v, collapse="|")
}

extract_field <- function(obj, field){
  if(!field %in% colnames(obj@meta.data)) return(NA_character_)
  first_nonempty(obj@meta.data[[field]])
}

fields <- c("source_blood_fraction","tissue_label","chemistry_version","library_type","reference_genome","frozen_or_fresh")
lib_meta <- vector("list", nrow(lib))

cat("=== BUILD LIBRARY CONTEXT TABLE ===\n")
for(i in seq_len(nrow(lib))){
  key <- as.character(lib$library_key[[i]])
  f <- as.character(lib$final_rds_resolved[[i]])
  if(!file.exists(f)) stop("Missing final RDS: ", f)
  cat(sprintf("[%3d/%3d] %s\n", i, nrow(lib), key))
  obj <- readRDS(f)

  row <- data.frame(project_id=as.character(lib$project_id[[i]]), library_key=key, stringsAsFactors=FALSE)
  for(field in fields){
    val <- NA_character_
    if(field %in% names(lib)) val <- first_nonempty(lib[[field]][i])
    if(is.na(val) || !nzchar(val)) val <- extract_field(obj, field)
    row[[field]] <- val
  }
  lib_meta[[i]] <- row
  rm(obj)
  gc(verbose=FALSE)
}

lib_meta <- do.call(rbind, lib_meta)
write.table(lib_meta, file=file.path(out_dir,"library_context_v1.tsv"), sep="\t", quote=TRUE, qmethod="double", row.names=FALSE)

lm <- match(as.character(meta$library_key), lib_meta$library_key)
stopifnot(!anyNA(lm))
for(field in fields) meta[[field]] <- lib_meta[[field]][lm]

safe_level <- function(x){
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- "NA"
  x
}

category_columns <- c(
  project="project_id",
  condition="condition_binary",
  blood_fraction="source_blood_fraction",
  tissue="tissue_label",
  chemistry="chemistry_version",
  library_type="library_type",
  reference="reference_genome"
)

resolutions <- c(0.6,0.8,1.0)
res_col <- function(r) paste0("cluster_res_", gsub("\\.", "p", format(r, trim=TRUE, scientific=FALSE)))

long <- list()
summary_rows <- list()

for(r in resolutions){
  cc <- res_col(r)
  if(!cc %in% names(cl)) stop("Missing clustering column: ", cc)
  labels <- as.character(cl[[cc]])

  for(k in sort(unique(labels))){
    ii <- which(labels == k)
    n <- length(ii)
    srow <- data.frame(resolution=r, cluster=k, n_cells=n, stringsAsFactors=FALSE)

    for(type in names(category_columns)){
      col <- category_columns[[type]]
      vals <- safe_level(meta[[col]][ii])
      tt <- sort(table(vals), decreasing=TRUE)
      top <- names(tt)[1]
      top_fraction <- as.numeric(tt[1]) / n
      srow[[paste0(type,"_n_levels")]] <- length(tt)
      srow[[paste0("top_",type)]] <- top
      srow[[paste0("top_",type,"_fraction")]] <- top_fraction
      long[[length(long)+1L]] <- data.frame(
        resolution=r, cluster=k, category=type, level=names(tt),
        n_cells=as.integer(tt), fraction=as.numeric(tt)/n,
        stringsAsFactors=FALSE
      )
    }
    summary_rows[[length(summary_rows)+1L]] <- srow
  }
}

long <- do.call(rbind, long)
summary_df <- do.call(rbind, summary_rows)

write.table(summary_df, file=file.path(out_dir,"cluster_provenance_summary_v1.tsv"), sep="\t", quote=TRUE, qmethod="double", row.names=FALSE)
write.table(long, file=file.path(out_dir,"cluster_provenance_long_v1.tsv"), sep="\t", quote=TRUE, qmethod="double", row.names=FALSE)

proj_context <- unique(lib_meta[,c("project_id","source_blood_fraction","tissue_label","chemistry_version","library_type","reference_genome"),drop=FALSE])
proj_context <- proj_context[order(proj_context$project_id),,drop=FALSE]
write.table(proj_context, file=file.path(out_dir,"project_library_context_v1.tsv"), sep="\t", quote=TRUE, qmethod="double", row.names=FALSE)

flag <- summary_df[summary_df$top_project_fraction >= 0.95,,drop=FALSE]
flag <- flag[order(flag$resolution, -flag$top_project_fraction, -flag$n_cells),,drop=FALSE]
write.table(flag, file=file.path(out_dir,"project_dominated_clusters_ge95pct_v1.tsv"), sep="\t", quote=TRUE, qmethod="double", row.names=FALSE)

writeLines(c(
  "PASS: Mono/DC provenance audit completed",
  paste0("n_cells=",length(cells)),
  paste0("n_libraries=",length(unique(meta$library_key))),
  paste0("resolutions=",paste(resolutions,collapse=",")),
  paste0("n_project_dominated_ge95pct=",nrow(flag))
), file.path(out_dir,"AUDIT_COMPLETE.ok"))

cat("\n=== MONO/DC PROVENANCE AUDIT SUMMARY ===\n")
cat("cells: ",length(cells),"\n",sep="")
cat("libraries: ",length(unique(meta$library_key)),"\n",sep="")
cat("project-dominated clusters >=95%: ",nrow(flag),"\n",sep="")
cat("\nTop flagged clusters:\n")
print(head(flag,30), row.names=FALSE)
cat("\nPASS: Mono/DC provenance audit completed\n")
