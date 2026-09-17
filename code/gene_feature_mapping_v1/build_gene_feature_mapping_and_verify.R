#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1L) args[[1]] else "."
features_arg <- if (length(args) >= 2L) args[[2]] else ""
project_root <- normalizePath(project_root, mustWork = TRUE)
suppressPackageStartupMessages(library(SeuratObject))
truthy <- function(x) tolower(trimws(as.character(x))) %in% c("1","true","t","yes","y")
rbind_fill <- function(xs) {
  if (!length(xs)) return(data.frame())
  all_names <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    for (nm in setdiff(all_names, names(x))) x[[nm]] <- NA
    x[, all_names, drop = FALSE]
  })
  out <- do.call(rbind, xs); rownames(out) <- NULL; out
}
read_features <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  x <- read.delim(con, header = FALSE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(x) < 2L) stop("Feature file has fewer than two columns: ", path)
  names(x)[1:2] <- c("gene_id", "gene_name")
  if (ncol(x) >= 3L) names(x)[3] <- "feature_type" else x$feature_type <- "Gene Expression"
  x[, c("gene_id","gene_name","feature_type"), drop = FALSE]
}
candidate_seurat_names <- function(features) {
  gene_name <- as.character(features$gene_name); gene_id <- as.character(features$gene_id)
  filled <- gene_name
  bad <- is.na(filled) | !nzchar(filled)
  filled[bad] <- gene_id[bad]
  list(gene_name = gene_name, gene_name_filled = filled, make_unique_dot = make.unique(filled), make_unique_dash = make.unique(filled, sep = "-"), gene_id = gene_id)
}
find_feature_files <- function(roots) {
  roots <- unique(roots[nzchar(roots) & dir.exists(roots)])
  out <- character()
  for (root in roots) {
    cmd <- sprintf("find %s -type f \\( -name 'features.tsv.gz' -o -name 'features.tsv' -o -name 'genes.tsv.gz' -o -name 'genes.tsv' \\) -print", shQuote(root))
    found <- tryCatch(system(cmd, intern = TRUE, ignore.stderr = TRUE), error = function(e) character())
    out <- c(out, found)
  }
  unique(out[file.exists(out)])
}
manifest_path <- file.path(project_root, "pre_integration", "final_qc_library_manifest.tsv")
if (!file.exists(manifest_path)) stop("Missing final-QC manifest: ", manifest_path)
manifest <- read.delim(manifest_path, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
required_manifest <- c("project_id","library_key","final_rds")
miss <- setdiff(required_manifest, names(manifest)); if (length(miss)) stop("Manifest missing: ", paste(miss, collapse=", "))
manifest <- manifest[file.exists(manifest$final_rds), , drop = FALSE]
if (!nrow(manifest)) stop("No readable final RDS files")
first_rds <- manifest$final_rds[[1]]
message("Reading reference RDS: ", first_rds)
ref_obj <- readRDS(first_rds)
assays <- SeuratObject::Assays(ref_obj)
ref_assay <- if ("RNA" %in% assays) "RNA" else SeuratObject::DefaultAssay(ref_obj)
ref_features <- rownames(ref_obj[[ref_assay]])
message("Reference assay: ", ref_assay, "; features: ", length(ref_features), "; cells: ", ncol(ref_obj))
rm(ref_obj); gc(verbose = FALSE)
roots <- c(Sys.getenv("ATLAS_CELLRANGER_OUTPUT_ROOT", unset = ""), file.path(project_root, "cellranger_count", "output"))
feature_files <- if (nzchar(features_arg)) normalizePath(features_arg, mustWork = TRUE) else find_feature_files(roots)
if (!length(feature_files)) stop("No feature file found; pass one explicitly as second argument")
message("Candidate feature files found: ", length(feature_files))
selected <- NULL; selection_rule <- NULL; diagnostics <- list()
for (path in feature_files) {
  features <- tryCatch(read_features(path), error = function(e) NULL)
  if (is.null(features)) next
  if (nrow(features) != length(ref_features)) {
    diagnostics[[length(diagnostics)+1L]] <- data.frame(path=path,n_features=nrow(features),best_rule=NA,best_match_fraction=NA)
    next
  }
  candidates <- candidate_seurat_names(features)
  frac <- vapply(candidates, function(x) mean(x == ref_features), numeric(1))
  best <- names(which.max(frac))[[1]]
  diagnostics[[length(diagnostics)+1L]] <- data.frame(path=path,n_features=nrow(features),best_rule=best,best_match_fraction=max(frac))
  if (max(frac) == 1) { selected <- list(path=path,features=features,seurat_names=candidates[[best]]); selection_rule <- best; break }
}
diagnostic_table <- rbind_fill(diagnostics)
diagnostic_path <- file.path(project_root,"pre_integration","gene_feature_source_candidates.tsv")
write.table(diagnostic_table, diagnostic_path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
if (is.null(selected)) stop("No feature file reproduced Seurat row names exactly. See: ", diagnostic_path)
message("Selected feature file: ", selected$path)
features <- selected$features; seurat_names <- selected$seurat_names
canonical <- as.character(features$gene_id)
ver <- grepl("^ENSG[0-9]+\\.[0-9]+$", canonical)
canonical[ver] <- sub("\\.[0-9]+$", "", canonical[ver])
mapping <- data.frame(
  feature_index=seq_len(nrow(features)), gene_id=as.character(features$gene_id), canonical_gene_id=canonical,
  gene_name=as.character(features$gene_name), feature_type=as.character(features$feature_type), seurat_feature_name=seurat_names,
  gene_id_is_versioned_ensembl=ver,
  gene_name_is_ensembl=grepl("^ENSG[0-9]+(?:\\.[0-9]+)?$", as.character(features$gene_name)),
  gene_name_missing=is.na(features$gene_name)|!nzchar(as.character(features$gene_name)),
  gene_name_duplicate_n=ave(seq_len(nrow(features)), as.character(features$gene_name), FUN=length),
  canonical_gene_id_duplicate_n=ave(seq_len(nrow(features)), canonical, FUN=length), stringsAsFactors=FALSE)
map_path <- file.path(project_root,"pre_integration","gene_feature_mapping.tsv.gz")
con <- gzfile(map_path,"wt"); write.table(mapping,con,sep="\t",quote=FALSE,row.names=FALSE,na="NA"); close(con)
verification <- vector("list", nrow(manifest))
message("Verifying ", nrow(manifest), " final RDS files...")
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i,,drop=FALSE]; obj <- tryCatch(readRDS(row$final_rds), error=function(e)e)
  if (inherits(obj,"error")) {
    verification[[i]] <- data.frame(project_id=row$project_id,library_key=row$library_key,final_rds=row$final_rds,rds_readable=FALSE,rna_present=NA,raw_present=NA,rna_feature_identical=NA,raw_feature_identical=NA,rna_raw_feature_identical=NA,rna_raw_cell_identical=NA,n_features_rna=NA,n_features_raw=NA,n_cells=NA,problem=conditionMessage(obj))
    next
  }
  aa <- SeuratObject::Assays(obj); rp <- "RNA" %in% aa; wp <- "RAW" %in% aa
  rf <- if (rp) rownames(obj[["RNA"]]) else character(); wf <- if (wp) rownames(obj[["RAW"]]) else character()
  rc <- if (rp) colnames(obj[["RNA"]]) else character(); wc <- if (wp) colnames(obj[["RAW"]]) else character()
  probs <- character(); if (!rp) probs <- c(probs,"RNA_assay_missing"); if (!wp) probs <- c(probs,"RAW_assay_missing")
  rfi <- rp && identical(rf,seurat_names); wfi <- wp && identical(wf,seurat_names); rwfi <- rp&&wp&&identical(rf,wf); rwci <- rp&&wp&&identical(rc,wc)
  if (rp&&!rfi) probs <- c(probs,"RNA_feature_order_mismatch"); if (wp&&!wfi) probs <- c(probs,"RAW_feature_order_mismatch"); if (rp&&wp&&!rwfi) probs <- c(probs,"RNA_RAW_feature_mismatch"); if (rp&&wp&&!rwci) probs <- c(probs,"RNA_RAW_cell_mismatch")
  verification[[i]] <- data.frame(project_id=row$project_id,library_key=row$library_key,final_rds=row$final_rds,rds_readable=TRUE,rna_present=rp,raw_present=wp,rna_feature_identical=rfi,raw_feature_identical=wfi,rna_raw_feature_identical=rwfi,rna_raw_cell_identical=rwci,n_features_rna=if(rp)length(rf)else NA,n_features_raw=if(wp)length(wf)else NA,n_cells=ncol(obj),problem=paste(probs,collapse=";"))
  rm(obj); if (i %% 10L == 0L) { message("  verified ",i,"/",nrow(manifest)); gc(verbose=FALSE) }
}
verification <- do.call(rbind,verification); rownames(verification)<-NULL
ver_path <- file.path(project_root,"pre_integration","gene_feature_verification_by_library.tsv")
write.table(verification,ver_path,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
summary <- data.frame(metric=c("selected_feature_file","seurat_feature_name_rule","n_features","n_unique_gene_id","n_unique_canonical_gene_id","n_versioned_ensembl_gene_id","n_gene_name_is_ensembl","n_missing_gene_name","n_duplicated_gene_name_rows","n_duplicated_canonical_gene_id_rows","n_final_rds_checked","n_final_rds_pass","n_final_rds_problem"), value=c(selected$path,selection_rule,nrow(mapping),length(unique(mapping$gene_id)),length(unique(mapping$canonical_gene_id)),sum(mapping$gene_id_is_versioned_ensembl),sum(mapping$gene_name_is_ensembl),sum(mapping$gene_name_missing),sum(mapping$gene_name_duplicate_n>1L),sum(mapping$canonical_gene_id_duplicate_n>1L),nrow(verification),sum(verification$problem==""),sum(verification$problem!="")), stringsAsFactors=FALSE)
sum_path <- file.path(project_root,"pre_integration","gene_feature_mapping_summary.tsv")
write.table(summary,sum_path,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
print(summary,row.names=FALSE)
problem_rows <- verification[verification$problem!="",,drop=FALSE]
cat("\nChecked:",nrow(verification)," Pass:",sum(verification$problem=="")," Problems:",nrow(problem_rows),"\n")
if (nrow(problem_rows)) { print(problem_rows,row.names=FALSE); stop("Feature verification failed") }
cat("PASS: one global Cell Ranger feature mapping is valid for all final RDS files\n")
