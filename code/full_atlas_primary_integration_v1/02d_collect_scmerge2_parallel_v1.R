#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".",mustWork=TRUE)

source(file.path(root,"full_atlas_primary_integration_v1","00_common.R"))

full_dir <- file.path(root,"pre_integration","full_atlas_primary_v1")
run_dir <- file.path(full_dir,"scmerge2_primary_runs")

plan <- data.frame(
  compartment=c("T_NK_combined","Monocyte_DC_combined"),
  ruvK=c(5L,2L),
  k_pseudoBulk=c(5L,5L),
  stringsAsFactors=FALSE
)

metrics <- vector("list",nrow(plan))

for(i in seq_len(nrow(plan))){
  f <- file.path(
    run_dir,
    paste0(
      plan$compartment[[i]],
      "__ruvK",plan$ruvK[[i]],
      "__kPB",plan$k_pseudoBulk[[i]],
      "__full_v1.rds"
    )
  )

  if(!file.exists(f)) stop("Missing final compartment RDS: ",f)

  z <- readRDS(f)

  if(!identical(as.character(z$compartment),plan$compartment[[i]])) {
    stop("Compartment mismatch: ",f)
  }

  metrics[[i]] <- z$metrics
}

metrics <- do.call(rbind,metrics)

out <- file.path(full_dir,"full_primary_scmerge2_metrics.tsv")

if(file.exists(out)){
  old <- read_tsv(out)
  same <- isTRUE(all.equal(
    old,
    metrics,
    check.attributes=FALSE
  ))
  if(!same){
    stop("Existing full_primary_scmerge2_metrics.tsv differs; refusing overwrite")
  }
  cat("Existing identical metrics retained\n")
} else {
  write_tsv(metrics,out)
}

cat("\n=== FULL PRIMARY SCMERGE2 METRICS ===\n")
print(metrics,row.names=FALSE)
cat("\nPASS: collected parallel scMerge2 metrics\n")
