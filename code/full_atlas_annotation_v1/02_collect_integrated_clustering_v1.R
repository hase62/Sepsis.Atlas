#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".", mustWork=TRUE)

base <- file.path(
  root, "atlas", "full_atlas_annotation_v1",
  "integrated_compartment_clustering"
)

comps <- c("T_NK_combined","Monocyte_DC_combined")
all <- list()

for(comp in comps){
  wd <- file.path(base,comp)
  ok <- file.path(wd,"CLUSTERING_COMPLETE.ok")
  f <- file.path(wd,"resolution_summary_v1.tsv")

  if(!file.exists(ok)) stop("Missing completion marker: ", ok)
  if(!file.exists(f)) stop("Missing resolution summary: ", f)

  x <- read.delim(f, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)
  all[[comp]] <- x
}

out <- do.call(rbind, all)

write.table(
  out,
  file=file.path(base,"integrated_clustering_overview.tsv"),
  sep="\t", quote=FALSE, row.names=FALSE
)

cat("=== INTEGRATED CLUSTERING OVERVIEW ===\n")
print(out,row.names=FALSE)
cat("\nPASS: integrated clustering overview collected\n")
