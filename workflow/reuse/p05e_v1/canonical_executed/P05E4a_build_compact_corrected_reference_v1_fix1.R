#!/usr/bin/env Rscript

# P05E4a fix1
# Handles the scMerge::getAdjustedMat() one-cell DelayedArray drop edge case.
# Compact the frozen P05E4 scMerge2 corrected reference without re-fitting models.
# Uses the scMerge-recommended large-data strategy:
#   saved fullalpha -> cosine-normalized native expression -> getAdjustedMat() by chunk
# Writes float32 binary matrices plus exact cell/gene order.
# Does NOT read the giant P05E4 adjusted-expression RDS files.
# Does NOT access cellranger_count/.

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".", mustWork=TRUE)

Sys.setenv(
  OMP_NUM_THREADS="1",
  OPENBLAS_NUM_THREADS="1",
  MKL_NUM_THREADS="1",
  NUMEXPR_NUM_THREADS="1"
)

source(file.path(
  root,
  "full_atlas_primary_integration_v1",
  "00_common.R"
))

required <- c(
  "scMerge",
  "batchelor",
  "DelayedMatrixStats",
  "BiocParallel",
  "data.table",
  "Matrix"
)

missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly=TRUE)
]

if(length(missing)) {
  stop("Missing packages: ", paste(missing, collapse=", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

out <- file.path(
  root,
  "publication",
  "scientific_data",
  paste0("reuse_reference_corrected_compact_v1__", tag)
)

matrix_dir <- file.path(out, "float32")
metadata_dir <- file.path(out, "metadata")

dir.create(matrix_dir, recursive=TRUE)
dir.create(metadata_dir, recursive=TRUE)

# ============================================================
# Locate canonical P05E4
# ============================================================

p4dirs <- list.dirs(
  file.path(root, "publication", "scientific_data"),
  recursive=FALSE,
  full.names=TRUE
)

p4dirs <- sort(
  p4dirs[
    grepl(
      "reuse_reference_corrected_expression_v1__[0-9]{8}_[0-9]{6}$",
      basename(p4dirs)
    )
  ]
)

stopifnot(length(p4dirs) >= 1L)
p4 <- tail(p4dirs, 1L)

p4_summary <- file.path(p4, "P05E4_SUMMARY.txt")
stopifnot(file.exists(p4_summary))

summary_text <- readLines(p4_summary, warn=FALSE)
stopifnot(any(summary_text == "status=PASS"))

manifest <- fread(
  file.path(
    p4,
    "corrected_reference_compartment_manifest.tsv"
  )
)

common_gene_tab <- fread(
  file.path(
    p4,
    "corrected_reference_common_genes.tsv"
  )
)

stopifnot(
  all(c("common_gene_order","feature") %in% names(common_gene_tab)),
  nrow(common_gene_tab) == 2971L,
  !anyDuplicated(common_gene_tab$feature)
)

setorder(common_gene_tab, common_gene_order)
common_genes <- common_gene_tab$feature

expected_compartments <- c(
  "T_NK",
  "Monocyte_DC",
  "B_plasma",
  "Neutrophil",
  "Platelet_megakaryocyte"
)

stopifnot(
  nrow(manifest) == 5L,
  setequal(manifest$compartment, expected_compartments)
)

manifest <- manifest[
  match(expected_compartments, compartment)
]

# ============================================================
# final-QC RDS map
# ============================================================

rds_files <- Sys.glob(
  file.path(
    root,
    "pre_integration",
    "GSE*",
    "rds_final_qc",
    "*.rds"
  )
)

stopifnot(length(rds_files) == 158L)

rds_map <- rbindlist(
  lapply(rds_files, function(f) {
    data.table(
      project_id=basename(dirname(dirname(f))),
      library_key=sub(
        "__final_qc\\.rds$",
        "",
        basename(f)
      ),
      rds=f
    )
  })
)

stopifnot(
  nrow(rds_map) == 158L,
  !anyDuplicated(
    paste(
      rds_map$project_id,
      rds_map$library_key,
      sep="|||"
    )
  )
)

# ============================================================
# Helpers
# ============================================================

chunk_indices <- function(n, chunk_size=5000L) {
  starts <- seq.int(1L, n, by=chunk_size)
  lapply(
    starts,
    function(a) {
      a:min(n, a + chunk_size - 1L)
    }
  )
}

cosine_normalize_chunk <- function(counts) {
  dat <- log_normalize_sparse(counts)

  batchelor::cosineNorm(
    dat,
    BPPARAM=BiocParallel::SerialParam()
  )
}

copy_file_raw <- function(src, dst_con, block_bytes=64L*1024L*1024L) {
  in_con <- file(src, open="rb")
  on.exit(close(in_con), add=TRUE)

  repeat {
    x <- readBin(
      in_con,
      what="raw",
      n=block_bytes
    )

    if(!length(x)) break

    writeBin(x, dst_con)
  }
}

manifest_out <- list()
cell_orders <- list()

# ============================================================
# Process each compartment independently
# ============================================================

for(ci in seq_len(nrow(manifest))) {

  compartment <- manifest$compartment[[ci]]

  cat(
    "\n============================================================\n",
    "P05E4a COMPACT CORRECTED REFERENCE: ", compartment, "\n",
    "============================================================\n",
    sep=""
  )

  model_file <- file.path(
    p4,
    manifest$model_file[[ci]]
  )

  meta_file <- file.path(
    p4,
    manifest$metadata_file[[ci]]
  )

  stopifnot(
    file.exists(model_file),
    file.exists(meta_file)
  )

  model <- readRDS(model_file)
  meta <- fread(meta_file)

  stopifnot(
    identical(model$compartment, compartment),
    nrow(meta) == length(model$cells),
    identical(meta$global_cell, model$cells),
    all(common_genes %in% model$chosen_hvg),
    !is.null(model$fullalpha)
  )

  model_features <- colnames(model$fullalpha)

  if(is.null(model_features)) {
    stop(compartment, ": fullalpha has no feature names")
  }

  if(anyDuplicated(model_features)) {
    stop(compartment, ": duplicated model feature names")
  }

  if(!all(common_genes %in% model_features)) {
    stop(
      compartment,
      ": common genes missing from fullalpha: ",
      sum(!common_genes %in% model_features)
    )
  }

  controls <- intersect(
    model$controls,
    model_features
  )

  if(length(controls) < 100L) {
    stop(compartment, ": fewer than 100 controls available")
  }

  # Preserve exact P05E4 cell order.
  meta[, compact_row_index :=
    seq_len(.N)
  ]

  lib_order <- unique(
    meta[, .(
      project_id,
      library_key
    )]
  )

  cat(
    "cells=", nrow(meta),
    "; model_features=", length(model_features),
    "; common_genes=", length(common_genes),
    "; controls=", length(controls),
    "; libraries=", nrow(lib_order),
    "\n",
    sep=""
  )

  # ----------------------------------------------------------
  # PASS 1: global row means of cosine-normalized expression
  # ----------------------------------------------------------

  gene_sum <- numeric(length(model_features))
  names(gene_sum) <- model_features

  n_seen <- 0L

  for(li in seq_len(nrow(lib_order))) {

    pid <- lib_order$project_id[[li]]
    key <- lib_order$library_key[[li]]

    mm <- meta[
      project_id == pid &
      library_key == key
    ]

    rr <- rds_map[
      project_id == pid &
      library_key == key
    ]

    if(nrow(rr) != 1L) {
      stop(
        compartment,
        ": RDS mapping failure for ",
        pid, "::", key
      )
    }

    obj <- readRDS(rr$rds[[1]])
    counts <- get_rna_counts(obj)

    if(!all(model_features %in% rownames(counts))) {
      stop(
        compartment,
        ": model features missing in ",
        key
      )
    }

    if(!all(mm$original_cell %in% colnames(counts))) {
      stop(
        compartment,
        ": model cells missing in ",
        key
      )
    }

    chunks <- chunk_indices(
      nrow(mm),
      chunk_size=5000L
    )

    for(ix in chunks) {

      cc <- mm$original_cell[ix]

      x <- counts[
        model_features,
        cc,
        drop=FALSE
      ]

      cn <- cosine_normalize_chunk(x)

      rs <- DelayedMatrixStats::rowSums2(cn)

      if(length(rs) != length(model_features)) {
        stop(compartment, ": row-sum length mismatch")
      }

      gene_sum <- gene_sum + as.numeric(rs)
      n_seen <- n_seen + length(ix)

      rm(x, cn, rs)
      invisible(gc())
    }

    rm(obj, counts)
    invisible(gc())

    cat(
      "MEAN PASS ",
      li, "/", nrow(lib_order),
      " ",
      pid, "::", key,
      " cells=", nrow(mm),
      "\n",
      sep=""
    )
  }

  stopifnot(n_seen == nrow(meta))

  adjusted_means <- gene_sum / n_seen
  names(adjusted_means) <- model_features

  if(any(!is.finite(adjusted_means))) {
    stop(compartment, ": non-finite adjusted means")
  }

  # ----------------------------------------------------------
  # PASS 2: getAdjustedMat by chunk -> float32 binary
  #
  # Binary layout:
  #   shape = [n_cells, n_common_genes]
  #   dtype = IEEE float32 little-endian
  #   order = C / row-major
  #
  # R's column-major gene x cell vectorization produces exactly
  # one contiguous gene vector per cell, i.e. row-major cells x genes.
  # ----------------------------------------------------------

  bin_file <- file.path(
    matrix_dir,
    paste0(
      compartment,
      "__scMerge2_adjusted_common2971_v1.f32"
    )
  )

  con <- file(bin_file, open="wb")

  first_expected <- NULL
  n_written_cells <- 0L
  n_singleton_chunks <- 0L

  for(li in seq_len(nrow(lib_order))) {

    pid <- lib_order$project_id[[li]]
    key <- lib_order$library_key[[li]]

    mm <- meta[
      project_id == pid &
      library_key == key
    ]

    rr <- rds_map[
      project_id == pid &
      library_key == key
    ]

    obj <- readRDS(rr$rds[[1]])
    counts <- get_rna_counts(obj)

    chunks <- chunk_indices(
      nrow(mm),
      chunk_size=5000L
    )

    for(ix in chunks) {

      cc <- mm$original_cell[ix]
      gcells <- mm$global_cell[ix]

      x <- counts[
        model_features,
        cc,
        drop=FALSE
      ]

      colnames(x) <- gcells

      cn <- cosine_normalize_chunk(x)

      # scMerge::getAdjustedMat() has a one-cell edge case:
      # after its internal transpose, subsetting one row can drop to a
      # vector and DelayedArray arithmetic fails. For a singleton chunk,
      # duplicate the same normalized cell only for the adjustment call,
      # then retain the first result. Because fullalpha and adjusted_means
      # are fixed and adjustment is cell-wise, the duplicate does not
      # alter the corrected value of the original cell.
      singleton_chunk <- length(gcells) == 1L

      if(singleton_chunk) {
        n_singleton_chunks <- n_singleton_chunks + 1L
        cn_for_adjust <- cbind(cn, cn)
        colnames(cn_for_adjust) <- c(
          gcells,
          paste0(gcells, "__P05E4a_singleton_duplicate")
        )
      } else {
        cn_for_adjust <- cn
      }

      adj <- scMerge::getAdjustedMat(
        exprsMat=cn_for_adjust,
        fullalpha=model$fullalpha,
        ctl=controls,
        adjusted_means=adjusted_means,
        ruvK=model$ruvK,
        return_subset_genes=common_genes
      )

      adj <- as.matrix(adj)

      if(singleton_chunk) {
        if(ncol(adj) != 2L) {
          stop(
            compartment,
            ": singleton workaround returned unexpected ncol=",
            ncol(adj)
          )
        }

        # Numerical equality of duplicated outputs is required.
        if(max(abs(adj[, 1L] - adj[, 2L])) > 1e-10) {
          stop(
            compartment,
            ": singleton duplicate outputs are not identical"
          )
        }

        adj <- adj[, 1L, drop=FALSE]
        colnames(adj) <- gcells
      }

      if(!setequal(rownames(adj), common_genes)) {
        close(con)
        stop(compartment, ": adjusted common-gene universe mismatch")
      }

      if(!identical(rownames(adj), common_genes)) {
        adj <- adj[
          match(common_genes, rownames(adj)),
          ,
          drop=FALSE
        ]
      }

      if(!identical(colnames(adj), gcells)) {
        adj <- adj[
          ,
          match(gcells, colnames(adj)),
          drop=FALSE
        ]
      }

      stopifnot(
        identical(rownames(adj), common_genes),
        identical(colnames(adj), gcells),
        nrow(adj) == length(common_genes),
        ncol(adj) == length(gcells),
        all(is.finite(adj))
      )

      if(is.null(first_expected)) {
        first_expected <- as.numeric(adj)[
          seq_len(min(100L, length(adj)))
        ]
      }

      # float32 little-endian
      writeBin(
        as.numeric(adj),
        con,
        size=4L,
        endian="little"
      )

      n_written_cells <- n_written_cells + length(ix)

      rm(x, cn, adj)
      invisible(gc())
    }

    rm(obj, counts)
    invisible(gc())

    cat(
      "WRITE PASS ",
      li, "/", nrow(lib_order),
      " ",
      pid, "::", key,
      " cells=", nrow(mm),
      "\n",
      sep=""
    )
  }

  close(con)

  stopifnot(n_written_cells == nrow(meta))

  expected_bytes <- as.double(
    nrow(meta)
  ) * as.double(
    length(common_genes)
  ) * 4

  actual_bytes <- file.info(bin_file)$size

  if(actual_bytes != expected_bytes) {
    stop(
      compartment,
      ": float32 size mismatch: expected=",
      expected_bytes,
      " actual=",
      actual_bytes
    )
  }

  # ----------------------------------------------------------
  # Read-back verification
  # ----------------------------------------------------------

  check_con <- file(bin_file, open="rb")

  got <- readBin(
    check_con,
    what="numeric",
    n=length(first_expected),
    size=4L,
    endian="little"
  )

  close(check_con)

  tol <- 5e-6

  if(
    length(got) != length(first_expected) ||
    max(abs(got - first_expected)) > tol
  ) {
    stop(
      compartment,
      ": float32 read-back verification failed"
    )
  }

  cell_order_file <- file.path(
    metadata_dir,
    paste0(
      compartment,
      "__corrected_reference_cell_order_v1.tsv"
    )
  )

  fwrite(
    meta[, .(
      compact_row_index,
      global_cell,
      original_cell,
      project_id,
      library_key,
      clinical_domain_std,
      condition_binary,
      final_compartment_v1,
      atlas_core_identity_v1,
      atlas_axis_v1,
      atlas_state_v1,
      integration_celltype
    )],
    cell_order_file,
    sep="\t",
    quote=TRUE,
    na="NA"
  )

  check_meta <- fread(cell_order_file)

  stopifnot(
    nrow(check_meta) == nrow(meta),
    identical(
      check_meta$global_cell,
      meta$global_cell
    )
  )

  cell_orders[[compartment]] <- copy(
    meta[, .(
      global_cell,
      project_id,
      library_key,
      final_compartment_v1
    )]
  )

  manifest_out[[length(manifest_out)+1L]] <- data.table(
    compartment=compartment,
    n_cells=nrow(meta),
    n_genes=length(common_genes),
    dtype="float32",
    endian="little",
    matrix_layout="cells_x_genes_C_order",
    matrix_file=file.path(
      "float32",
      basename(bin_file)
    ),
    cell_order_file=file.path(
      "metadata",
      basename(cell_order_file)
    ),
    bytes=actual_bytes,
    size_gib=actual_bytes / 1024^3,
    ruvK=model$ruvK,
    k_pseudoBulk=model$k_pseudoBulk,
    scMerge_version=model$scMerge_version,
    correction_scope="within_compartment",
    DE_use="NO",
    readback_check="PASS",
    singleton_chunks_workaround=n_singleton_chunks
  )

  rm(
    model,
    meta,
    adjusted_means,
    gene_sum,
    first_expected,
    got,
    check_meta
  )

  invisible(gc())
}

# ============================================================
# Combined cell order and compact manifest
# ============================================================

compact_manifest <- rbindlist(manifest_out)

fwrite(
  compact_manifest,
  file.path(
    out,
    "compact_corrected_reference_manifest.tsv"
  ),
  sep="\t"
)

combined_cells <- rbindlist(
  lapply(
    expected_compartments,
    function(compartment) {
      x <- cell_orders[[compartment]]
      x[, compartment_block :=
        compartment
      ]
      x
    }
  ),
  use.names=TRUE,
  fill=TRUE
)

combined_cells[, combined_row_index :=
  seq_len(.N)
]

fwrite(
  combined_cells,
  file.path(
    metadata_dir,
    "SepsisAtlas_corrected_reference_combined_cell_order_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  na="NA"
)

fwrite(
  common_gene_tab,
  file.path(
    metadata_dir,
    "SepsisAtlas_corrected_reference_common_genes_v1.tsv"
  ),
  sep="\t"
)

stopifnot(
  nrow(combined_cells) == sum(compact_manifest$n_cells),
  nrow(combined_cells) == 597927L,
  nrow(common_gene_tab) == 2971L,
  !anyDuplicated(combined_cells$global_cell)
)

# ============================================================
# Optional combined float32 file: concatenate compartment blocks
# ============================================================

combined_bin <- file.path(
  matrix_dir,
  "SepsisAtlas_scMerge2_corrected_reference_common2971_v1.f32"
)

combined_con <- file(
  combined_bin,
  open="wb"
)

for(compartment in expected_compartments) {

  row <- compact_manifest[
    compact_manifest$compartment == compartment
  ]

  stopifnot(nrow(row) == 1L)

  src <- file.path(
    out,
    row$matrix_file[[1]]
  )

  stopifnot(file.exists(src))

  copy_file_raw(src, combined_con)
}

close(combined_con)

expected_combined_bytes <- as.double(
  nrow(combined_cells)
) * as.double(
  length(common_genes)
) * 4

actual_combined_bytes <- file.info(combined_bin)$size

if(actual_combined_bytes != expected_combined_bytes) {
  stop(
    "Combined float32 size mismatch: expected=",
    expected_combined_bytes,
    " actual=",
    actual_combined_bytes
  )
}

# ============================================================
# SHA256
# ============================================================

oldwd <- getwd()
setwd(out)

files_for_sha <- sort(
  list.files(
    ".",
    recursive=TRUE,
    full.names=FALSE
  )
)

files_for_sha <- files_for_sha[
  files_for_sha != "SHA256SUMS_P05E4a.txt"
]

sha_lines <- vapply(
  files_for_sha,
  function(f) {
    z <- system2(
      "sha256sum",
      args=f,
      stdout=TRUE,
      stderr=TRUE
    )

    status <- attr(z, "status")

    if(!is.null(status) && status != 0L) {
      stop("sha256sum failed for ", f)
    }

    z[[1]]
  },
  character(1)
)

writeLines(
  sha_lines,
  "SHA256SUMS_P05E4a.txt"
)

check <- system2(
  "sha256sum",
  args=c("-c", "SHA256SUMS_P05E4a.txt"),
  stdout=TRUE,
  stderr=TRUE
)

check_status <- attr(check, "status")

if(!is.null(check_status) && check_status != 0L) {
  cat(paste(check, collapse="\n"), "\n")
  stop("P05E4a SHA256 verification failed")
}

setwd(oldwd)

# ============================================================
# Summary
# ============================================================

summary_lines <- c(
  "===== P05E4a COMPACT CORRECTED REFERENCE =====",
  "status=PASS",
  paste0(
    "corrected_compartments=",
    nrow(compact_manifest)
  ),
  paste0(
    "corrected_reference_cells=",
    nrow(combined_cells)
  ),
  paste0(
    "common_reference_genes=",
    length(common_genes)
  ),
  "dtype=float32",
  "matrix_layout=cells_x_genes_C_order",
  paste0(
    "combined_matrix_bytes=",
    actual_combined_bytes
  ),
  paste0(
    "combined_matrix_GiB=",
    sprintf("%.3f", actual_combined_bytes / 1024^3)
  ),
  "reconstructed_from_saved_P05E4_fullalpha=YES",
  "scMerge2_refit=NO",
  "giant_P05E4_adjusted_RDS_read=NO",
  "float32_readback_check=PASS",
  paste0(
    "singleton_chunks_workaround=",
    sum(compact_manifest$singleton_chunks_workaround)
  ),
  "SHA256=PASS",
  "corrected_reference_DE_use=NO",
  "native_processed_counts_DE_use=YES",
  "annotation_modified=NO",
  "P05B_modified=NO",
  "cellranger_count=NOT_ACCESSED"
)

writeLines(
  summary_lines,
  file.path(
    out,
    "P05E4a_SUMMARY.txt"
  )
)

# ============================================================
# Root handoff copies
# ============================================================

handoff <- c(
  "P05E4a_SUMMARY.txt",
  "compact_corrected_reference_manifest.tsv"
)

for(f in handoff) {

  src <- file.path(out, f)

  dst <- file.path(
    root,
    paste0(
      "P05E4a_",
      sub("^P05E4a_", "", f)
    )
  )

  file.copy(
    src,
    dst,
    overwrite=TRUE
  )
}

cat(
  readLines(
    file.path(out, "P05E4a_SUMMARY.txt")
  ),
  sep="\n"
)

cat("\n\n===== COMPACT MANIFEST =====\n")
print(compact_manifest)

cat(
  "\nOUT_DIR=",
  out,
  "\n",
  sep=""
)

cat(
  "Root handoff files copied: ",
  length(handoff),
  "\n",
  sep=""
)
