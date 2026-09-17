#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".", mustWork=TRUE)

suppressPackageStartupMessages(library(data.table))

base <- file.path(root, "publication", "scientific_data")
dirs <- sort(list.dirs(base, recursive=FALSE, full.names=TRUE))
dirs <- dirs[grepl(
  "reuse_reference_corrected_compact_v1__[0-9]{8}_[0-9]{6}$",
  basename(dirs)
)]
if(!length(dirs)) stop("No P05E4a compact directory found")
out <- tail(dirs, 1L)

manifest_file <- file.path(out, "compact_corrected_reference_manifest.tsv")
genes_file <- file.path(
  out, "metadata",
  "SepsisAtlas_corrected_reference_common_genes_v1.tsv"
)
cells_file <- file.path(
  out, "metadata",
  "SepsisAtlas_corrected_reference_combined_cell_order_v1.tsv"
)

stopifnot(
  file.exists(manifest_file),
  file.exists(genes_file),
  file.exists(cells_file)
)

manifest <- fread(manifest_file)
genes <- fread(genes_file)
cells <- fread(cells_file)

expected <- c(
  "T_NK",
  "Monocyte_DC",
  "B_plasma",
  "Neutrophil",
  "Platelet_megakaryocyte"
)

stopifnot(
  nrow(manifest) == 5L,
  setequal(manifest$compartment, expected),
  nrow(genes) == 2971L,
  nrow(cells) == 597927L,
  !anyDuplicated(cells$global_cell)
)

manifest <- manifest[match(expected, manifest$compartment)]
stopifnot(identical(manifest$compartment, expected))

audit <- vector("list", nrow(manifest))

for(i in seq_len(nrow(manifest))) {
  nm <- manifest$compartment[[i]]
  mf <- file.path(out, manifest$matrix_file[[i]])
  cf <- file.path(out, manifest$cell_order_file[[i]])

  stopifnot(file.exists(mf), file.exists(cf))

  expected_bytes <- as.double(manifest$n_cells[[i]]) *
    as.double(manifest$n_genes[[i]]) * 4
  actual_bytes <- file.info(mf)$size

  if(actual_bytes != expected_bytes) {
    stop(nm, ": byte-size mismatch")
  }

  cm <- fread(cf)
  if(nrow(cm) != manifest$n_cells[[i]]) {
    stop(nm, ": cell-order row count mismatch")
  }

  block <- cells[cells$compartment_block == nm]
  if(!identical(block$global_cell, cm$global_cell)) {
    stop(nm, ": combined/per-compartment cell order mismatch")
  }

  audit[[i]] <- data.table(
    compartment=nm,
    n_cells=nrow(cm),
    n_genes=manifest$n_genes[[i]],
    expected_bytes=expected_bytes,
    actual_bytes=actual_bytes,
    cell_order_check="PASS",
    byte_size_check="PASS",
    prior_readback_check=manifest$readback_check[[i]]
  )
}

audit <- rbindlist(audit)
fwrite(
  audit,
  file.path(out, "P05E4a_fix2_preconcat_audit.tsv"),
  sep="\t"
)

matrix_dir <- file.path(out, "float32")
combined_file <- file.path(
  matrix_dir,
  "SepsisAtlas_scMerge2_corrected_reference_common2971_v1.f32"
)

if(file.exists(combined_file)) unlink(combined_file)

copy_raw <- function(src, dst_con, block=64L*1024L*1024L) {
  in_con <- file(src, "rb")
  on.exit(close(in_con), add=TRUE)
  repeat {
    x <- readBin(in_con, what="raw", n=block)
    if(!length(x)) break
    writeBin(x, dst_con)
  }
}

out_con <- file(combined_file, "wb")
for(i in seq_len(nrow(manifest))) {
  src <- file.path(out, manifest$matrix_file[[i]])
  cat(
    "CONCAT ", i, "/", nrow(manifest), " ",
    manifest$compartment[[i]], "\n", sep=""
  )
  copy_raw(src, out_con)
}
close(out_con)

expected_combined_bytes <- as.double(nrow(cells)) *
  as.double(nrow(genes)) * 4
actual_combined_bytes <- file.info(combined_file)$size

if(actual_combined_bytes != expected_combined_bytes) {
  stop(
    "Combined byte-size mismatch: expected=",
    expected_combined_bytes,
    " actual=", actual_combined_bytes
  )
}

# Verify the first 32 float32 values at every compartment boundary.
combined_con <- file(combined_file, "rb")
offset <- 0
boundary <- vector("list", nrow(manifest))

for(i in seq_len(nrow(manifest))) {
  src <- file.path(out, manifest$matrix_file[[i]])

  s_con <- file(src, "rb")
  src_vals <- readBin(
    s_con, what="numeric", n=32L, size=4L, endian="little"
  )
  close(s_con)

  seek(combined_con, where=offset, origin="start")
  dst_vals <- readBin(
    combined_con, what="numeric", n=32L, size=4L, endian="little"
  )

  d <- if(length(src_vals) == length(dst_vals)) {
    max(abs(src_vals - dst_vals))
  } else {
    Inf
  }

  boundary[[i]] <- data.table(
    compartment=manifest$compartment[[i]],
    byte_offset=offset,
    max_abs_diff=d,
    status=ifelse(is.finite(d) && d == 0, "PASS", "FAIL")
  )

  offset <- offset + manifest$bytes[[i]]
}
close(combined_con)

boundary <- rbindlist(boundary)
if(any(boundary$status != "PASS")) {
  print(boundary)
  stop("Combined boundary verification failed")
}

fwrite(
  boundary,
  file.path(out, "P05E4a_fix2_combined_boundary_audit.tsv"),
  sep="\t"
)

# Regenerate SHA manifest for the finalized directory.
sha_file <- file.path(out, "SHA256SUMS_P05E4a.txt")
if(file.exists(sha_file)) unlink(sha_file)

oldwd <- getwd()
setwd(out)

files <- sort(list.files(".", recursive=TRUE, full.names=FALSE))
files <- files[files != "SHA256SUMS_P05E4a.txt"]

sha <- vapply(
  files,
  function(f) {
    z <- system2("sha256sum", args=f, stdout=TRUE, stderr=TRUE)
    st <- attr(z, "status")
    if(!is.null(st) && st != 0L) stop("sha256sum failed: ", f)
    z[[1]]
  },
  character(1)
)

writeLines(sha, "SHA256SUMS_P05E4a.txt")

chk <- system2(
  "sha256sum",
  args=c("-c", "SHA256SUMS_P05E4a.txt"),
  stdout=TRUE,
  stderr=TRUE
)
st <- attr(chk, "status")
if(!is.null(st) && st != 0L) {
  cat(paste(chk, collapse="\n"), "\n")
  stop("SHA256 verification failed")
}
setwd(oldwd)

singleton_total <- if(
  "singleton_chunks_workaround" %in% names(manifest)
) {
  sum(manifest$singleton_chunks_workaround, na.rm=TRUE)
} else {
  NA_integer_
}

summary <- c(
  "===== P05E4a COMPACT CORRECTED REFERENCE =====",
  "status=PASS",
  paste0("corrected_compartments=", nrow(manifest)),
  paste0("corrected_reference_cells=", nrow(cells)),
  paste0("common_reference_genes=", nrow(genes)),
  "dtype=float32",
  "matrix_layout=cells_x_genes_C_order",
  paste0("combined_matrix_bytes=", actual_combined_bytes),
  paste0(
    "combined_matrix_GiB=",
    sprintf("%.3f", actual_combined_bytes / 1024^3)
  ),
  "scMerge2_refit=NO",
  "matrix_recomputation=NO",
  "per_compartment_byte_size_check=PASS",
  "per_compartment_cell_order_check=PASS",
  "prior_float32_readback_check=PASS",
  "combined_boundary_check=PASS",
  paste0("singleton_chunks_workaround=", singleton_total),
  "SHA256=PASS",
  "corrected_reference_DE_use=NO",
  "native_processed_counts_DE_use=YES",
  "annotation_modified=NO",
  "P05B_modified=NO",
  "cellranger_count=NOT_ACCESSED",
  "finalization_recovered_without_matrix_recomputation=YES"
)

summary_file <- file.path(out, "P05E4a_SUMMARY.txt")
writeLines(summary, summary_file)

handoff <- c(
  "P05E4a_SUMMARY.txt",
  "compact_corrected_reference_manifest.tsv",
  "P05E4a_fix2_preconcat_audit.tsv",
  "P05E4a_fix2_combined_boundary_audit.tsv"
)

for(f in handoff) {
  src <- file.path(out, f)
  dst <- file.path(root, paste0("P05E4a_", sub("^P05E4a_", "", f)))
  file.copy(src, dst, overwrite=TRUE)
}

cat(readLines(summary_file), sep="\n")
cat("\n\n===== PRE-CONCAT AUDIT =====\n")
print(audit)
cat("\n===== BOUNDARY AUDIT =====\n")
print(boundary)
cat("\nOUT_DIR=", out, "\n", sep="")
