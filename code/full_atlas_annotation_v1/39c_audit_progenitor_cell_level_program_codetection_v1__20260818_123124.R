#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3L){
  stop("Usage: script <root> <tag> <cluster_dir>")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]
cluster_dir <- normalizePath(args[[3]], mustWork=TRUE)

suppressPackageStartupMessages({
  library(Matrix)
})

source(
  file.path(
    root,
    "full_atlas_primary_integration_v1",
    "00_common.R"
  )
)

DISCOVERY_PROJECT <- "GSE216007"
N_EXPECTED <- 26962L
RESOLUTION <- "native_cluster_r0p4"

# Descriptive co-detection thresholds only.
# They are NOT annotation thresholds.
programs <- list(

  strict_HSC=c(
    "HLF","AVP","MECOM","MLLT3",
    "HOPX","PROM1","CRHBP","PRDM16"
  ),

  lymphoid_priming=c(
    "FLT3","DNTT","MME","SELL",
    "IL7R","TCF7","CD7"
  ),

  B_commitment=c(
    "EBF1","PAX5","VPREB1","VPREB3",
    "IGLL1","CD79A","CD19","MS4A1"
  ),

  granulocytic=c(
    "MPO","ELANE","AZU1","PRTN3",
    "CTSG","CLEC12A","CEBPA","GFI1"
  ),

  megakaryocytic=c(
    "ITGA2B","GP1BB","GP9","PF4","PPBP",
    "VWF","NFE2","FLI1","TUBB1","TREML1"
  ),

  erythroid_nonglobin=c(
    "GATA1","KLF1","ALAS2","AHSP",
    "ANK1","SLC4A1","CA1","BPGM","TFRC"
  ),

  eos_baso_mast=c(
    "HDC","MS4A2","ENPP3","CLC","PRG2",
    "CPA3","TPSAB1","KIT","GATA2","CCR3"
  ),

  T_NK=c(
    "NKG7","GNLY","PRF1","GZMB",
    "CD3D","CD3E","TRAC","TRBC1"
  ),

  cycling=c(
    "MKI67","TOP2A","UBE2C","CENPF",
    "CDC20","CDC45","PCNA","TYMS"
  )
)

thresholds <- c(
  strict_HSC=3L,
  lymphoid_priming=3L,
  B_commitment=3L,
  granulocytic=3L,
  megakaryocytic=2L,
  erythroid_nonglobin=3L,
  eos_baso_mast=3L,
  T_NK=2L,
  cycling=3L
)

read_tsv <- function(path){

  if(!file.exists(path)){
    stop("Missing file: ", path)
  }

  con <- if(grepl("\\.gz$", path)){
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }

  tryCatch(
    read.delim(
      con,
      sep="\t",
      quote="\"",
      comment.char="",
      stringsAsFactors=FALSE,
      check.names=FALSE
    ),
    finally=close(con)
  )
}

assignment_file <- file.path(
  cluster_dir,
  "progenitor_dominant_project_native_cluster_assignments_v1.tsv.gz"
)

assign <- read_tsv(assignment_file)

stopifnot(
  nrow(assign) == N_EXPECTED,
  !anyDuplicated(assign$global_cell),
  RESOLUTION %in% names(assign)
)

assign$global_cell <- as.character(assign$global_cell)
assign$project_id <- as.character(assign$project_id)
assign$library_key <- as.character(assign$library_key)
assign$cluster <- as.character(assign[[RESOLUTION]])

stopifnot(
  all(assign$project_id == DISCOVERY_PROJECT),
  setequal(unique(assign$cluster), as.character(0:13))
)

# ============================================================
# Recover original cell barcodes
# ============================================================

transfer_dir <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "annotation_transfer_by_library"
)

transfer_files <- sort(
  list.files(
    transfer_dir,
    pattern="\\.tsv\\.gz$",
    full.names=TRUE
  )
)

meta_list <- list()

for(f in transfer_files){

  d <- read_tsv(f)

  req <- c(
    "global_cell",
    "original_cell",
    "project_id",
    "library_key",
    "integration_compartment_primary_v1"
  )

  stopifnot(all(req %in% names(d)))

  d <- d[
    as.character(d$project_id) == DISCOVERY_PROJECT &
      as.character(
        d$integration_compartment_primary_v1
      ) == "Progenitor",
    req,
    drop=FALSE
  ]

  if(nrow(d)){
    meta_list[[length(meta_list)+1L]] <- d
  }
}

meta <- do.call(rbind, meta_list)

for(nm in c(
  "global_cell",
  "original_cell",
  "project_id",
  "library_key"
)){
  meta[[nm]] <- as.character(meta[[nm]])
}

stopifnot(
  nrow(meta) == N_EXPECTED,
  !anyDuplicated(meta$global_cell)
)

mi <- match(
  assign$global_cell,
  meta$global_cell
)

stopifnot(!anyNA(mi))

assign$original_cell <- meta$original_cell[mi]

stopifnot(
  assign$library_key == meta$library_key[mi]
)

rm(meta, meta_list, mi)
invisible(gc())

# ============================================================
# Resolved source libraries
# ============================================================

lib <- read_tsv(
  file.path(
    root,
    "pre_integration",
    "pilot_unintegrated_large_v1",
    "resolved_library_table.tsv"
  )
)

stopifnot(
  all(
    c(
      "project_id",
      "library_key",
      "final_rds_resolved"
    ) %in% names(lib)
  )
)

lib$project_id <- as.character(lib$project_id)
lib$library_key <- as.character(lib$library_key)

lib$key <- paste(
  lib$project_id,
  lib$library_key,
  sep="|||"
)

libraries <- sort(
  unique(assign$library_key)
)

all_panel_genes <- unique(
  unlist(programs)
)

cell_rows <- list()

# ============================================================
# Native SoupX RNA: cell-level detection
# ============================================================

for(ii in seq_along(libraries)){

  lk <- libraries[[ii]]

  key <- paste(
    DISCOVERY_PROJECT,
    lk,
    sep="|||"
  )

  li <- match(
    key,
    lib$key
  )

  if(is.na(li)){
    stop("Library lookup failed: ", key)
  }

  source_rds <- as.character(
    lib$final_rds_resolved[[li]]
  )

  if(!file.exists(source_rds)){
    stop("Missing source RDS: ", source_rds)
  }

  a <- assign[
    assign$library_key == lk,
    ,
    drop=FALSE
  ]

  obj <- readRDS(source_rds)
  counts0 <- get_rna_counts(obj)

  stopifnot(
    all(a$original_cell %in% colnames(counts0))
  )

  panel_present <- intersect(
    all_panel_genes,
    rownames(counts0)
  )

  x <- counts0[
    panel_present,
    a$original_cell,
    drop=FALSE
  ]

  colnames(x) <- a$global_cell

  out <- data.frame(
    global_cell=a$global_cell,
    library_key=lk,
    cluster=a$cluster,
    stringsAsFactors=FALSE
  )

  for(program in names(programs)){

    genes <- intersect(
      programs[[program]],
      rownames(x)
    )

    if(!length(genes)){
      n_detected <- rep(
        0L,
        ncol(x)
      )
    } else {
      n_detected <- as.integer(
        Matrix::colSums(
          x[genes, , drop=FALSE] > 0
        )
      )
    }

    out[[paste0(program, "_n")]] <-
      n_detected

    out[[paste0(program, "_flag")]] <-
      n_detected >=
      thresholds[[program]]
  }

  cell_rows[[length(cell_rows)+1L]] <- out

  cat(
    sprintf(
      "READ %d/%d %s cells=%d\n",
      ii,
      length(libraries),
      lk,
      nrow(out)
    )
  )

  rm(obj, counts0, x, out)
  invisible(gc())
}

cell_df <- do.call(
  rbind,
  cell_rows
)

rownames(cell_df) <- NULL

stopifnot(
  nrow(cell_df) == N_EXPECTED,
  !anyDuplicated(cell_df$global_cell)
)

# ============================================================
# Explicit pairwise co-detection
# ============================================================

cell_df$HSC_and_lymphoid <-
  cell_df$strict_HSC_flag &
  cell_df$lymphoid_priming_flag

cell_df$HSC_and_granulocytic <-
  cell_df$strict_HSC_flag &
  cell_df$granulocytic_flag

cell_df$ery_and_meg <-
  cell_df$erythroid_nonglobin_flag &
  cell_df$megakaryocytic_flag

cell_df$ery_and_cycling <-
  cell_df$erythroid_nonglobin_flag &
  cell_df$cycling_flag

cell_df$meg_and_cycling <-
  cell_df$megakaryocytic_flag &
  cell_df$cycling_flag

cell_df$B_and_lymphoid <-
  cell_df$B_commitment_flag &
  cell_df$lymphoid_priming_flag

cell_df$B_and_TNK <-
  cell_df$B_commitment_flag &
  cell_df$T_NK_flag

cell_df$eos_and_granulocytic <-
  cell_df$eos_baso_mast_flag &
  cell_df$granulocytic_flag

# ============================================================
# Summaries
# ============================================================

flag_names <- grep(
  "_flag$",
  names(cell_df),
  value=TRUE
)

pair_names <- c(
  "HSC_and_lymphoid",
  "HSC_and_granulocytic",
  "ery_and_meg",
  "ery_and_cycling",
  "meg_and_cycling",
  "B_and_lymphoid",
  "B_and_TNK",
  "eos_and_granulocytic"
)

summarize_cells <- function(d){

  out <- data.frame(
    n_cells=nrow(d),
    stringsAsFactors=FALSE
  )

  for(nm in flag_names){
    out[[paste0("fraction_", sub("_flag$", "", nm))]] <-
      mean(d[[nm]])
  }

  for(nm in pair_names){
    out[[paste0("fraction_", nm)]] <-
      mean(d[[nm]])
  }

  out
}

cluster_summary <- do.call(
  rbind,
  lapply(
    sort(unique(as.integer(cell_df$cluster))),
    function(cl){

      d <- cell_df[
        as.integer(cell_df$cluster) == cl,
        ,
        drop=FALSE
      ]

      z <- summarize_cells(d)
      z$cluster <- as.character(cl)

      z[
        ,
        c(
          "cluster",
          setdiff(names(z), "cluster")
        ),
        drop=FALSE
      ]
    }
  )
)

library_cluster_summary <- do.call(
  rbind,
  lapply(
    split(
      cell_df,
      paste(
        cell_df$library_key,
        cell_df$cluster,
        sep="|||"
      )
    ),
    function(d){

      z <- summarize_cells(d)

      z$library_key <-
        d$library_key[[1]]

      z$cluster <-
        d$cluster[[1]]

      z[
        ,
        c(
          "library_key",
          "cluster",
          setdiff(
            names(z),
            c("library_key","cluster")
          )
        ),
        drop=FALSE
      ]
    }
  )
)

rownames(cluster_summary) <- NULL
rownames(library_cluster_summary) <- NULL

# ============================================================
# Equal-library medians
# ============================================================

metric_cols <- setdiff(
  names(library_cluster_summary),
  c(
    "library_key",
    "cluster",
    "n_cells"
  )
)

library_equal_summary <- do.call(
  rbind,
  lapply(
    as.character(0:13),
    function(cl){

      d <- library_cluster_summary[
        library_cluster_summary$cluster == cl &
          library_cluster_summary$n_cells >= 5L,
        ,
        drop=FALSE
      ]

      z <- data.frame(
        cluster=cl,
        n_eligible_libraries=nrow(d),
        stringsAsFactors=FALSE
      )

      for(nm in metric_cols){

        z[[paste0(
          "median_library_",
          nm
        )]] <- if(nrow(d)){
          median(
            d[[nm]],
            na.rm=TRUE
          )
        } else {
          NA_real_
        }
      }

      z
    }
  )
)

# ============================================================
# Output
# ============================================================

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "progenitor_cell_level_program_codetection_v1__",
    tag
  )
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

cell_file <- file.path(
  out_dir,
  paste0(
    "progenitor_cell_level_program_flags_v1__",
    tag,
    ".tsv.gz"
  )
)

con <- gzfile(
  cell_file,
  "wt"
)

tryCatch(
  write.table(
    cell_df,
    con,
    sep="\t",
    quote=TRUE,
    row.names=FALSE
  ),
  finally=close(con)
)

if(system2(
  "gzip",
  c("-t", cell_file),
  stdout=FALSE,
  stderr=FALSE
) != 0L){
  stop("gzip integrity failure")
}

write.table(
  cluster_summary,
  file.path(
    out_dir,
    "progenitor_r0p4_cell_program_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  library_cluster_summary,
  file.path(
    out_dir,
    "progenitor_r0p4_library_cell_program_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  library_equal_summary,
  file.path(
    out_dir,
    "progenitor_r0p4_library_equal_program_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# Re-read cell file.
chk <- read_tsv(cell_file)

stopifnot(
  nrow(chk) == N_EXPECTED,
  identical(
    as.character(chk$global_cell),
    as.character(cell_df$global_cell)
  )
)

sha_files <- c(
  cell_file,
  file.path(
    out_dir,
    "progenitor_r0p4_cell_program_summary_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_r0p4_library_cell_program_summary_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_r0p4_library_equal_program_summary_v1.tsv"
  )
)

sha_file <- file.path(
  out_dir,
  paste0(
    "SHA256SUMS_v1__",
    tag,
    ".txt"
  )
)

status <- system2(
  "sha256sum",
  sha_files,
  stdout=sha_file
)

if(status != 0L){
  stop("sha256 creation failed")
}

status <- system2(
  "sha256sum",
  c("-c", sha_file),
  stdout=FALSE,
  stderr=FALSE
)

if(status != 0L){
  stop("sha256 verification failed")
}

writeLines(
  c(
    "PASS",
    "Progenitor cell-level program co-detection audit v1",
    "discovery_project=GSE216007",
    "n_cells=26962",
    "backbone=native_cluster_r0p4",
    "",
    "strict HSC excludes generic CD34/MEIS1-only interpretation",
    "erythroid program excludes globins",
    "cell-level lineage co-detection retained",
    "library-level and library-equal summaries retained",
    "",
    "descriptive thresholds only",
    "no automatic annotation",
    "condition not used",
    "non-discovery projects not used",
    "native SoupX-corrected RNA",
    "cellranger_count not accessed",
    "",
    "gzip_integrity=PASS",
    "re_read_count_order=PASS",
    "sha256_manifest=PASS"
  ),
  file.path(
    out_dir,
    "PROGENITOR_CELL_LEVEL_PROGRAM_CODETECTION_COMPLETE.ok"
  )
)

cat(
  "\n===== CELL-WEIGHTED CLUSTER SUMMARY =====\n"
)

print(
  cluster_summary,
  row.names=FALSE
)

cat(
  "\n===== LIBRARY-EQUAL SUMMARY =====\n"
)

print(
  library_equal_summary,
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Progenitor cell-level program audit completed\n"
)
