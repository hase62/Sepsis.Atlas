#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".", mustWork=TRUE)

suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
})

comp <- "T_NK_combined"
cl_rds <- file.path(root, "atlas", "full_atlas_annotation_v1", "integrated_compartment_clustering", comp, "integrated_clustering_v1.rds")
meta_file <- file.path(root, "pre_integration", "full_atlas_primary_v1", "scmerge2_primary_runs", "T_NK_combined__full_v1_cells.tsv.gz")
lib_file <- file.path(root, "pre_integration", "pilot_unintegrated_large_v1", "resolved_library_table.tsv")
out_dir <- file.path(root, "atlas", "full_atlas_annotation_v1", "tnk_native_marker_aggregation_v1")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
stopifnot(file.exists(cl_rds), file.exists(meta_file), file.exists(lib_file))

x <- readRDS(cl_rds)
cells <- as.character(x$cells)
clusters <- x$clusters
stopifnot(nrow(clusters) == length(cells), identical(rownames(clusters), cells))

meta <- read.delim(gzfile(meta_file), sep="\t", quote="\"", comment.char="", stringsAsFactors=FALSE, check.names=FALSE)
required_meta <- c("global_cell","project_id","library_key")
miss <- setdiff(required_meta, names(meta))
if(length(miss)) stop("Missing cell metadata: ", paste(miss,collapse=", "))
mm <- match(cells, as.character(meta$global_cell))
if(anyNA(mm)) stop("Cell metadata mismatch: ",sum(is.na(mm)))
meta <- meta[mm,,drop=FALSE]
stopifnot(identical(as.character(meta$global_cell), cells))

lib <- read.delim(lib_file, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
required_lib <- c("project_id","library_key","final_rds_resolved")
miss <- setdiff(required_lib, names(lib))
if(length(miss)) stop("Missing library table columns: ",paste(miss,collapse=", "))

resolutions <- c(0.4,0.6)
res_col <- function(r) paste0("cluster_res_",gsub("\\.","p",format(r,trim=TRUE,scientific=FALSE)))
for(r in resolutions) if(!res_col(r) %in% names(clusters)) stop("Missing clustering column: ",res_col(r))

get_rna_counts <- function(obj){
  if(!"RNA" %in% names(obj@assays)) stop("RNA assay missing")
  ans <- tryCatch(SeuratObject::LayerData(obj[["RNA"]], layer="counts"), error=function(e) NULL)
  if(is.null(ans)) ans <- tryCatch(SeuratObject::GetAssayData(obj, assay="RNA", layer="counts"), error=function(e) NULL)
  if(is.null(ans)) ans <- tryCatch(SeuratObject::GetAssayData(obj, assay="RNA", slot="counts"), error=function(e) NULL)
  if(is.null(ans)) stop("Could not retrieve RNA counts")
  ans
}

used_libs <- unique(as.character(meta$library_key))
li <- match(used_libs, as.character(lib$library_key))
if(anyNA(li)) stop("Missing used libraries: ",paste(used_libs[is.na(li)],collapse=", "))
lib <- lib[li,,drop=FALSE]

first_obj <- readRDS(as.character(lib$final_rds_resolved[[1]]))
feature_ref <- rownames(get_rna_counts(first_obj))
rm(first_obj)
gc(verbose=FALSE)
if(is.null(feature_ref) || !length(feature_ref)) stop("No RNA feature names")
nfeat <- length(feature_ref)

cat("Features=",nfeat,"\n",sep="")
cat("T/NK cells=",length(cells),"\n",sep="")
cat("Libraries=",nrow(lib),"\n",sep="")

agg <- list()
for(r in resolutions){
  labs <- as.character(clusters[[res_col(r)]])
  lev <- sort(unique(labs))
  agg[[as.character(r)]] <- list(
    resolution=r,
    levels=lev,
    sum_counts=matrix(0,nrow=nfeat,ncol=length(lev),dimnames=list(feature_ref,lev)),
    detect_counts=matrix(0,nrow=nfeat,ncol=length(lev),dimnames=list(feature_ref,lev)),
    n_cells=setNames(integer(length(lev)),lev)
  )
}

manifest <- list()
cat("\n=== AGGREGATE NATIVE SOUPX-CORRECTED RNA ===\n")

for(i in seq_len(nrow(lib))){
  key <- as.character(lib$library_key[[i]])
  project <- as.character(lib$project_id[[i]])
  f <- as.character(lib$final_rds_resolved[[i]])
  ii <- which(as.character(meta$library_key) == key)
  if(!length(ii)) next
  cat(sprintf("[%3d/%3d] %s  cells=%d\n",i,nrow(lib),key,length(ii)))
  if(!file.exists(f)) stop("Missing final RDS: ",f)
  obj <- readRDS(f)
  counts <- get_rna_counts(obj)
  if(!identical(rownames(counts),feature_ref)) stop("RNA feature order mismatch in ",key)

  prefix <- paste0(key,"___")
  global <- cells[ii]
  if(!all(startsWith(global,prefix))) stop("Global-cell prefix mismatch in ",key)
  original <- substring(global,nchar(prefix)+1L)
  if(!all(original %in% colnames(counts))){
    bad <- original[!original %in% colnames(counts)]
    stop("Original cells missing in ",key,": ",paste(head(bad,10),collapse=", "))
  }

  m <- counts[,original,drop=FALSE]
  for(r in resolutions){
    rr <- as.character(r)
    labs <- as.character(clusters[[res_col(r)]][ii])
    lev <- unique(labs)
    for(k in lev){
      jj <- which(labs == k)
      nc <- length(jj)
      sc <- Matrix::rowSums(m[,jj,drop=FALSE])
      dc <- Matrix::rowSums(m[,jj,drop=FALSE] > 0)
      agg[[rr]]$sum_counts[,k] <- agg[[rr]]$sum_counts[,k] + as.numeric(sc)
      agg[[rr]]$detect_counts[,k] <- agg[[rr]]$detect_counts[,k] + as.numeric(dc)
      agg[[rr]]$n_cells[[k]] <- agg[[rr]]$n_cells[[k]] + nc
      manifest[[length(manifest)+1L]] <- data.frame(
        project_id=project, library_key=key, resolution=r, cluster=k,
        n_cells=nc, total_native_counts=sum(sc), stringsAsFactors=FALSE
      )
    }
  }
  rm(obj,counts,m)
  gc(verbose=FALSE)
}

manifest <- do.call(rbind,manifest)
write.table(manifest, file=file.path(out_dir,"library_cluster_manifest_v1.tsv"), sep="\t",quote=TRUE,qmethod="double",row.names=FALSE)

technical_flag <- function(g){
  grepl("^(MT-|RPL[0-9]|RPS[0-9]|HBA[12]$|HBB$|HBD$|HBG[12]$|MALAT1$|NEAT1$)",g)
}

rank_outputs <- list()
for(r in resolutions){
  rr <- as.character(r)
  a <- agg[[rr]]
  stopifnot(sum(a$n_cells) == length(cells))
  total_gene <- rowSums(a$sum_counts)
  total_detect <- rowSums(a$detect_counts)
  total_counts <- sum(total_gene)
  total_cells <- sum(a$n_cells)
  all_stats <- list()
  top_stats <- list()

  for(k in a$levels){
    nk <- a$n_cells[[k]]
    nout <- total_cells - nk
    cin <- a$sum_counts[,k]
    cout <- total_gene - cin
    uin <- sum(cin)
    uout <- total_counts - uin
    pct_in <- a$detect_counts[,k] / nk
    pct_out <- if(nout > 0) (total_detect-a$detect_counts[,k]) / nout else 0
    cpm_in <- if(uin > 0) cin/uin*1e6 else rep(0,nfeat)
    cpm_out <- if(uout > 0) cout/uout*1e6 else rep(0,nfeat)
    log2fc <- log2(cpm_in + 0.1) - log2(cpm_out + 0.1)
    delta_pct <- pct_in - pct_out

    d <- data.frame(
      resolution=r, cluster=k, gene=feature_ref, n_cells_cluster=nk,
      native_cpm_in=cpm_in, native_cpm_out=cpm_out,
      log2FC_cluster_vs_rest=log2fc, pct_in=pct_in, pct_out=pct_out,
      delta_pct=delta_pct, technical_gene=technical_flag(feature_ref),
      marker_candidate=(pct_in >= 0.10 & delta_pct >= 0.05 & log2fc >= 0.5),
      stringsAsFactors=FALSE
    )

    ord <- order(!d$marker_candidate,-d$delta_pct,-d$log2FC_cluster_vs_rest,-d$pct_in,d$gene)
    d <- d[ord,,drop=FALSE]
    d$rank_within_cluster <- seq_len(nrow(d))
    all_stats[[length(all_stats)+1L]] <- d
    cand <- d[d$marker_candidate & !d$technical_gene,,drop=FALSE]
    if(!nrow(cand)) cand <- d[!d$technical_gene,,drop=FALSE]
    top_stats[[length(top_stats)+1L]] <- head(cand,200)
  }

  all_stats <- do.call(rbind,all_stats)
  top_stats <- do.call(rbind,top_stats)
  tag <- gsub("\\.","p",format(r,trim=TRUE,scientific=FALSE))

  con <- gzfile(file.path(out_dir,paste0("native_marker_stats_res_",tag,"_v1.tsv.gz")),open="wt")
  tryCatch(write.table(all_stats,con,sep="\t",quote=TRUE,qmethod="double",row.names=FALSE),finally=close(con))
  write.table(top_stats,file=file.path(out_dir,paste0("top_native_markers_res_",tag,"_v1.tsv")),sep="\t",quote=TRUE,qmethod="double",row.names=FALSE)

  cs <- data.frame(
    resolution=r,
    cluster=a$levels,
    n_cells=as.integer(a$n_cells[a$levels]),
    n_libraries=vapply(a$levels,function(k) sum(manifest$resolution == r & manifest$cluster == k & manifest$n_cells > 0),integer(1)),
    n_projects=vapply(a$levels,function(k) length(unique(manifest$project_id[manifest$resolution == r & manifest$cluster == k & manifest$n_cells > 0])),integer(1)),
    n_marker_candidates=vapply(a$levels,function(k) sum(all_stats$cluster == k & all_stats$marker_candidate & !all_stats$technical_gene),integer(1)),
    stringsAsFactors=FALSE
  )
  write.table(cs,file=file.path(out_dir,paste0("cluster_marker_summary_res_",tag,"_v1.tsv")),sep="\t",quote=TRUE,qmethod="double",row.names=FALSE)
  rank_outputs[[rr]] <- list(resolution=r,cluster_summary=cs)
}

saveRDS(list(
  compartment=comp,
  resolutions=resolutions,
  feature_reference=feature_ref,
  aggregates=agg,
  library_cluster_manifest=manifest,
  ranking_definition=list(
    marker_candidate="pct_in>=0.10 & delta_pct>=0.05 & log2FC_cluster_vs_rest>=0.5",
    ranking="candidate first; delta_pct desc; log2FC desc; pct_in desc",
    technical_gene_flag="MT/RPL/RPS/HBA/HBB/HBD/HBG/MALAT1/NEAT1"
  ),
  source_clustering_rds=cl_rds,
  source_cell_metadata=meta_file,
  source_library_table=lib_file
), file.path(out_dir,"tnk_native_marker_aggregates_v1.rds"), compress=FALSE)

writeLines(c(
  "PASS: T/NK native marker aggregation completed",
  paste0("n_cells=",length(cells)),
  paste0("n_features=",nfeat),
  paste0("n_libraries=",length(unique(meta$library_key))),
  paste0("resolutions=",paste(resolutions,collapse=","))
), file.path(out_dir,"MARKER_AGGREGATION_COMPLETE.ok"))

cat("\n=== T/NK MARKER AGGREGATION COMPLETE ===\n")
for(r in resolutions) print(rank_outputs[[as.character(r)]]$cluster_summary,row.names=FALSE)
cat("\nPASS: T/NK native marker aggregation completed\n")
