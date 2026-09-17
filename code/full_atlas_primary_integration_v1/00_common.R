suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))

write_tsv <- function(x, path) {
  con <- if (grepl("\\.gz$", path, ignore.case=TRUE)) gzfile(path, "wt") else file(path, "wt")
  on.exit(close(con), add=TRUE)
  write.table(x, con, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
}

read_tsv <- function(path) read.delim(path, sep="\t", stringsAsFactors=FALSE, check.names=FALSE)

get_rna_counts <- function(obj) {
  if (!"RNA" %in% as.character(names(obj@assays))) stop("RNA assay missing")
  layers <- as.character(SeuratObject::Layers(obj[["RNA"]]))
  if (!"counts" %in% layers) stop("RNA counts layer missing")
  x <- SeuratObject::LayerData(obj[["RNA"]], layer="counts")
  if (!inherits(x, "sparseMatrix")) x <- as(x, "dgCMatrix")
  x
}

get_rna_data <- function(obj) {
  if (!"RNA" %in% as.character(names(obj@assays))) stop("RNA assay missing")
  layers <- as.character(SeuratObject::Layers(obj[["RNA"]]))
  if (!"data" %in% layers) stop("RNA data layer missing")
  x <- SeuratObject::LayerData(obj[["RNA"]], layer="data")
  if (!inherits(x, "sparseMatrix")) x <- as(x, "dgCMatrix")
  x
}

log_normalize_sparse <- function(counts, scale_factor=10000) {
  libsize <- Matrix::colSums(counts)
  if (any(libsize <= 0)) stop("Zero library-size cell detected")
  dn <- dimnames(counts)
  x <- counts %*% Diagonal(x=scale_factor/libsize)
  x@x <- log1p(x@x)
  x <- as(x, "dgCMatrix")
  dimnames(x) <- dn
  x
}

rss_gb <- function() {
  f <- "/proc/self/status"
  if (!file.exists(f)) return(NA_real_)
  x <- readLines(f, warn=FALSE)
  z <- grep("^VmRSS:", x, value=TRUE)
  if (!length(z)) return(NA_real_)
  kb <- suppressWarnings(as.numeric(gsub("[^0-9]", "", z[[1]])))
  kb/1024/1024
}

healthy_binary <- function(x) {
  y <- tolower(trimws(as.character(x)))
  ifelse(y %in% c("healthy","healthy_control","healthy control","control","hc"), "healthy", "disease")
}
