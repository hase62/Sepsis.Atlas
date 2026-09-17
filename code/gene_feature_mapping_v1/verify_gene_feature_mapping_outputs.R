#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if(length(args)>=1L) args[[1]] else ".", mustWork=TRUE)
sp <- file.path(root,"pre_integration","gene_feature_mapping_summary.tsv")
vp <- file.path(root,"pre_integration","gene_feature_verification_by_library.tsv")
mp <- file.path(root,"pre_integration","gene_feature_mapping.tsv.gz")
stopifnot(file.exists(sp),file.exists(vp),file.exists(mp))
s <- read.delim(sp,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
v <- read.delim(vp,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
getv <- function(m) s$value[s$metric==m][1]
stopifnot(as.integer(getv("n_final_rds_checked"))==nrow(v),as.integer(getv("n_final_rds_problem"))==0L,all(v$problem==""))
cat("PASS: gene feature mapping outputs are internally consistent\nRDS checked: ",nrow(v),"\nFeatures: ",getv("n_features"),"\nSource: ",getv("selected_feature_file"),"\n",sep="")
