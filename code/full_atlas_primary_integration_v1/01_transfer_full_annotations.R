#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork=TRUE)
source(file.path(root,"full_atlas_primary_integration_v1","00_common.R"))
suppressPackageStartupMessages(library(BiocNeighbors))
set.seed(20260811)

base <- file.path(root,"pre_integration")
pilot_dir <- file.path(base,"pilot_unintegrated_large_v1")
freeze_dir <- file.path(pilot_dir,"annotation_freeze_v2")
out_dir <- file.path(base,"full_atlas_primary_v1")
by_lib_dir <- file.path(out_dir,"annotation_transfer_by_library")
pca_dir <- file.path(out_dir,"native_pca30_by_library")
index_dir <- file.path(root,"tmp","full_atlas_annoy_index_v1")
dir.create(by_lib_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(pca_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(index_dir, recursive=TRUE, showWarnings=FALSE)

pilot_path <- file.path(freeze_dir,"pilot_unintegrated_large_v1__annotation_freeze_v2.rds")
lib_path <- file.path(pilot_dir,"resolved_library_table.tsv")
for (f in c(pilot_path,lib_path)) if (!file.exists(f)) stop("Missing input: ",f)

pilot <- readRDS(pilot_path)
lib <- read_tsv(lib_path)
if (nrow(lib)!=158L) stop("Expected 158 libraries; found ",nrow(lib))
if (!all(file.exists(lib$final_rds_resolved))) stop("At least one final RDS path does not exist")

ref_pca <- Embeddings(pilot,"pca")
ref_pca <- ref_pca[,seq_len(min(30L,ncol(ref_pca))),drop=FALSE]
ref_cells <- rownames(ref_pca)
required_ref_md <- c("annotation_large_v2_broad","annotation_large_v2_fine","integration_compartment_v2","integration_celltype_v2")
miss <- setdiff(required_ref_md,names(pilot[[]]))
if (length(miss)) stop("Pilot freeze metadata missing: ",paste(miss,collapse=", "))
ref_md <- pilot[[]][ref_cells,required_ref_md,drop=FALSE]

pilot_data <- get_rna_data(pilot)
projection_features <- rownames(pilot_data)
loadings <- Loadings(pilot,"pca")
loadings <- loadings[,colnames(ref_pca),drop=FALSE]
pca_genes <- rownames(loadings)
x <- pilot_data[pca_genes,,drop=FALSE]
mu <- Matrix::rowMeans(x)
mu2 <- Matrix::rowMeans(x*x)
sdv <- sqrt(pmax(mu2-mu^2,0))
if (any(!is.finite(sdv)) || any(sdv<=1e-8)) stop("Invalid SD among PCA genes")
coef <- sweep(loadings,1,sdv,FUN="/")
offset <- as.numeric(crossprod(mu/sdv,loadings))
names(offset) <- colnames(loadings)

set.seed(20260811)
check_cells <- if (ncol(x)<=5000L) colnames(x) else sample(colnames(x),5000L)
reproj <- as.matrix(Matrix::t(x[,check_cells,drop=FALSE]) %*% coef)
reproj <- sweep(reproj,2,offset,FUN="-")
stored <- ref_pca[check_cells,colnames(reproj),drop=FALSE]
projection_rmse <- sqrt(mean((reproj-stored)^2))
projection_cor <- cor(as.numeric(reproj),as.numeric(stored))
preflight <- data.frame(metric=c("n_reference_cells","n_reference_pca_dims","n_pca_genes","projection_check_cells","projection_rmse","projection_correlation"), value=c(nrow(ref_pca),ncol(ref_pca),length(pca_genes),length(check_cells),projection_rmse,projection_cor))
write_tsv(preflight,file.path(out_dir,"full_annotation_projection_preflight.tsv"))
if (!is.finite(projection_cor) || projection_cor<0.999) stop("Pilot PCA projection validation failed")

oldwd <- getwd(); setwd(index_dir); on.exit(setwd(oldwd),add=TRUE)
ann_index <- BiocNeighbors::buildIndex(ref_pca,BNPARAM=BiocNeighbors::AnnoyParam())

vote_one <- function(indices, query_cell, label_vec, within=NULL, k=30L) {
  ids <- indices[ref_cells[indices] != query_cell]
  if (!is.null(within)) ids <- ids[within[ids]]
  ids <- head(ids,k)
  labs <- as.character(label_vec[ids]); labs <- labs[!is.na(labs)&nzchar(labs)]
  if (!length(labs)) return(c(label=NA_character_,fraction=NA_character_))
  tab <- sort(table(labs),decreasing=TRUE)
  c(label=names(tab)[[1]],fraction=as.character(as.integer(tab[[1]])/sum(tab)))
}

nn_ref <- BiocNeighbors::queryKNN(X=ref_pca,query=ref_pca,k=31L,BNINDEX=ann_index)
loocv_comp <- character(nrow(ref_pca)); loocv_fine <- character(nrow(ref_pca))
loocv_comp_frac <- numeric(nrow(ref_pca)); loocv_fine_frac <- numeric(nrow(ref_pca))
for (i in seq_len(nrow(ref_pca))) {
  ids <- nn_ref$index[i,]
  vc <- vote_one(ids,ref_cells[[i]],ref_md$integration_compartment_v2,k=30L)
  comp <- vc[["label"]]; loocv_comp[[i]] <- comp; loocv_comp_frac[[i]] <- as.numeric(vc[["fraction"]])
  vf <- vote_one(ids,ref_cells[[i]],ref_md$integration_celltype_v2,within=(ref_md$integration_compartment_v2==comp),k=30L)
  loocv_fine[[i]] <- vf[["label"]]; loocv_fine_frac[[i]] <- as.numeric(vf[["fraction"]])
}
loocv <- data.frame(metric=c("compartment_accuracy","fine_label_accuracy","median_compartment_vote_fraction","median_fine_vote_fraction"), value=c(mean(loocv_comp==ref_md$integration_compartment_v2),mean(loocv_fine==ref_md$integration_celltype_v2),median(loocv_comp_frac,na.rm=TRUE),median(loocv_fine_frac,na.rm=TRUE)))
write_tsv(loocv,file.path(out_dir,"full_annotation_transfer_loocv.tsv"))
if (loocv$value[loocv$metric=="compartment_accuracy"]<0.95) stop("LOOCV compartment transfer accuracy <0.95")

project_new <- function(data_new) {
  z <- data_new[pca_genes,,drop=FALSE]
  score <- as.matrix(Matrix::t(z) %*% coef)
  score <- sweep(score,2,offset,FUN="-")
  colnames(score) <- colnames(ref_pca)
  score
}

vote_query <- function(nn_idx,query_cells) {
  n <- nrow(nn_idx)
  out <- data.frame(integration_compartment_full_v1=character(n),integration_compartment_vote_fraction=numeric(n),integration_celltype_full_v1=character(n),integration_celltype_vote_fraction=numeric(n),annotation_broad_full_v1=character(n),annotation_broad_vote_fraction=numeric(n),transfer_confidence=character(n),stringsAsFactors=FALSE)
  for (i in seq_len(n)) {
    ids <- nn_idx[i,]; ids <- ids[ref_cells[ids] != query_cells[[i]]]; ids <- head(ids,30L)
    tc <- sort(table(ref_md$integration_compartment_v2[ids]),decreasing=TRUE)
    comp <- names(tc)[[1]]; comp_frac <- as.integer(tc[[1]])/sum(tc)
    ids_comp <- ids[ref_md$integration_compartment_v2[ids]==comp]
    tf <- sort(table(ref_md$integration_celltype_v2[ids_comp]),decreasing=TRUE)
    fine <- names(tf)[[1]]; fine_frac <- as.integer(tf[[1]])/sum(tf)
    tb <- sort(table(ref_md$annotation_large_v2_broad[ids_comp]),decreasing=TRUE)
    broad <- names(tb)[[1]]; broad_frac <- as.integer(tb[[1]])/sum(tb)
    conf <- if (comp_frac>=0.80 && fine_frac>=0.60) "high" else if (comp_frac>=0.60 && fine_frac>=0.50) "medium" else "low"
    out[i,] <- list(comp,comp_frac,fine,fine_frac,broad,broad_frac,conf)
  }
  out
}

library_summaries <- vector("list",nrow(lib))
for (i in seq_len(nrow(lib))) {
  project_id <- as.character(lib$project_id[[i]]); library_key <- as.character(lib$library_key[[i]])
  tag <- safe_name(paste(project_id,library_key,sep="__"))
  meta_out <- file.path(by_lib_dir,paste0(tag,".tsv.gz")); pca_out <- file.path(pca_dir,paste0(tag,"__pca30.rds"))
  if (file.exists(meta_out) && file.exists(pca_out)) {
    z <- read_tsv(meta_out)
    library_summaries[[i]] <- data.frame(project_id=project_id,library_key=library_key,n_cells=nrow(z),n_high=sum(z$transfer_confidence=="high"),n_medium=sum(z$transfer_confidence=="medium"),n_low=sum(z$transfer_confidence=="low"),skipped_existing=TRUE)
    cat("SKIP existing ",i,"/",nrow(lib),": ",library_key,"\n",sep="")
    next
  }
  obj <- readRDS(as.character(lib$final_rds_resolved[[i]])); counts <- get_rna_counts(obj)
  if (length(setdiff(pca_genes,rownames(counts)))) stop("PCA genes missing in ",library_key)
  original_cells <- colnames(counts); global_cells <- paste0(safe_name(library_key),"___",original_cells)
  if (length(setdiff(projection_features,rownames(counts)))) stop("Projection features missing in ",library_key)
  data_reduced <- log_normalize_sparse(counts[projection_features,,drop=FALSE])
  pca_new <- project_new(data_reduced); rownames(pca_new) <- global_cells
  nn <- BiocNeighbors::queryKNN(X=ref_pca,query=pca_new,k=31L,BNINDEX=ann_index)
  votes <- vote_query(nn$index,global_cells)
  md_out <- data.frame(project_id=project_id,library_key=library_key,original_cell=original_cells,global_cell=global_cells,clinical_domain_std=as.character(lib$clinical_domain_std[[i]]),condition_binary=healthy_binary(lib$clinical_domain_std[[i]]),votes,stringsAsFactors=FALSE)
  md_out$integration_compartment_primary_v1 <- ifelse(md_out$transfer_confidence=="low","Deferred_transfer_low_confidence",md_out$integration_compartment_full_v1)
  write_tsv(md_out,meta_out); saveRDS(pca_new,pca_out,compress=FALSE)
  library_summaries[[i]] <- data.frame(project_id=project_id,library_key=library_key,n_cells=nrow(md_out),n_high=sum(md_out$transfer_confidence=="high"),n_medium=sum(md_out$transfer_confidence=="medium"),n_low=sum(md_out$transfer_confidence=="low"),skipped_existing=FALSE)
  rm(obj,counts,data_reduced,pca_new,nn,votes,md_out); gc(verbose=FALSE)
  cat("ANNOTATED ",i,"/",nrow(lib)," ",library_key,"; RSS=",sprintf("%.2f",rss_gb())," GB\n",sep="")
}

library_summary <- do.call(rbind,library_summaries)
write_tsv(library_summary,file.path(out_dir,"full_annotation_transfer_library_summary.tsv"))
files <- list.files(by_lib_dir,pattern="\\.tsv\\.gz$",full.names=TRUE)
if (length(files)!=158L) stop("Expected 158 transfer files; found ",length(files))
agg <- lapply(files,function(f){z<-read_tsv(f); z[,c("integration_compartment_primary_v1","integration_celltype_full_v1","transfer_confidence","project_id"),drop=FALSE]})
agg <- do.call(rbind,agg)
if (nrow(agg)!=665816L) stop("Expected 665816 full Atlas cells; found ",nrow(agg))
comp_summary <- as.data.frame(table(integration_compartment_primary_v1=agg$integration_compartment_primary_v1),stringsAsFactors=FALSE); names(comp_summary)[2] <- "n_cells"
conf_summary <- as.data.frame(table(transfer_confidence=agg$transfer_confidence),stringsAsFactors=FALSE); names(conf_summary)[2] <- "n_cells"
write_tsv(comp_summary,file.path(out_dir,"full_annotation_transfer_compartment_counts.tsv")); write_tsv(conf_summary,file.path(out_dir,"full_annotation_transfer_confidence_counts.tsv"))
summary <- data.frame(metric=c("n_libraries","n_cells","n_high_confidence","n_medium_confidence","n_low_confidence","low_confidence_fraction","pilot_loocv_compartment_accuracy","pilot_loocv_fine_accuracy"),value=c(158,nrow(agg),sum(agg$transfer_confidence=="high"),sum(agg$transfer_confidence=="medium"),sum(agg$transfer_confidence=="low"),mean(agg$transfer_confidence=="low"),loocv$value[loocv$metric=="compartment_accuracy"],loocv$value[loocv$metric=="fine_label_accuracy"]))
write_tsv(summary,file.path(out_dir,"full_annotation_transfer_summary.tsv"))
cat("\n=== FULL ANNOTATION TRANSFER SUMMARY ===\n"); print(summary,row.names=FALSE); cat("\n=== FULL PRIMARY COMPARTMENT COUNTS ===\n"); print(comp_summary,row.names=FALSE); cat("\nPASS: full Atlas annotation transfer completed\n")
