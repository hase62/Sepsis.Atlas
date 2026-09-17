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

MIN_CELLS_PER_LIBRARY_CLUSTER <- 5L

programs <- list(

  primitive_HSC=c(
    "HLF","AVP","MECOM","MEIS1","MLLT3",
    "HOPX","PROM1","CRHBP","PRDM16","CD34"
  ),

  lymphoid_priming=c(
    "FLT3","IL7R","DNTT","SELL","SATB1",
    "MME","CD7","LTB","TCF7"
  ),

  B_commitment=c(
    "EBF1","PAX5","VPREB1","VPREB3",
    "IGLL1","CD79A","CD19","MS4A1"
  ),

  myeloid_granulocytic=c(
    "MPO","ELANE","AZU1","PRTN3","CTSG",
    "CLEC12A","CSF3R","CEBPA","GFI1","LYZ"
  ),

  monocyte=c(
    "LYZ","CTSS","FCER1G","LST1","TYMP",
    "CSF1R","FCGR3A","S100A8","S100A9","CTSD"
  ),

  megakaryocytic=c(
    "ITGA2B","GP1BB","GP9","PF4","PPBP",
    "VWF","NFE2","FLI1","TUBB1","TREML1"
  ),

  erythroid=c(
    "GATA1","KLF1","ALAS2","AHSP","ANK1",
    "SLC4A1","CA1","BPGM","TFRC",
    "HBB","HBA1","HBA2"
  ),

  eos_baso_mast=c(
    "HDC","MS4A2","ENPP3","CLC","PRG2",
    "IL3RA","GATA2","CPA3","TPSAB1",
    "KIT","CCR3"
  ),

  dendritic=c(
    "FCER1A","CD1C","CLEC10A",
    "CLEC9A","IRF8"
  ),

  T_NK=c(
    "NKG7","GNLY","PRF1","GZMB",
    "CD3D","CD3E","TRAC","TRBC1"
  ),

  cycling=c(
    "MKI67","TOP2A","UBE2C","CENPF",
    "CDC20","CDC45","PCNA","TYMS"
  ),

  IFN=c(
    "ISG15","IFI6","IFI27","IFIT1",
    "IFIT2","IFIT3","MX1","OAS1",
    "OASL","IFITM3"
  ),

  immediate_early=c(
    "FOS","FOSB","JUN","JUNB",
    "JUND","DUSP1","EGR1","IER2"
  )
)

read_tsv <- function(path){

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
  all(as.character(assign$project_id) == DISCOVERY_PROJECT),
  RESOLUTION %in% names(assign)
)

assign$global_cell <- as.character(assign$global_cell)
assign$library_key <- as.character(assign$library_key)
assign$cluster <- as.character(assign[[RESOLUTION]])

# ============================================================
# Recover original barcodes
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
    as.character(d$integration_compartment_primary_v1) == "Progenitor",
    req,
    drop=FALSE
  ]

  if(nrow(d)){
    meta_list[[length(meta_list)+1L]] <- d
  }
}

meta <- do.call(rbind, meta_list)

for(nm in c(
  "global_cell","original_cell",
  "project_id","library_key"
)){
  meta[[nm]] <- as.character(meta[[nm]])
}

stopifnot(
  nrow(meta) == N_EXPECTED,
  !anyDuplicated(meta$global_cell)
)

mi <- match(assign$global_cell, meta$global_cell)
stopifnot(!anyNA(mi))

assign$original_cell <- meta$original_cell[mi]

# ============================================================
# Source libraries
# ============================================================

lib <- read_tsv(
  file.path(
    root,
    "pre_integration",
    "pilot_unintegrated_large_v1",
    "resolved_library_table.tsv"
  )
)

lib$project_id <- as.character(lib$project_id)
lib$library_key <- as.character(lib$library_key)

lib$key <- paste(
  lib$project_id,
  lib$library_key,
  sep="|||"
)

libraries <- sort(unique(assign$library_key))
clusters <- as.character(0:13)

gene_rows <- list()
complexity_rows <- list()

all_panel_genes <- unique(unlist(programs))

# ============================================================
# Native absolute expression
# ============================================================

for(ii in seq_along(libraries)){

  lk <- libraries[[ii]]

  key <- paste(
    DISCOVERY_PROJECT,
    lk,
    sep="|||"
  )

  li <- match(key, lib$key)

  if(is.na(li)){
    stop("Library lookup failure: ", key)
  }

  source_rds <- as.character(
    lib$final_rds_resolved[[li]]
  )

  obj <- readRDS(source_rds)
  counts0 <- get_rna_counts(obj)

  a <- assign[
    assign$library_key == lk,
    ,
    drop=FALSE
  ]

  stopifnot(
    all(a$original_cell %in% colnames(counts0))
  )

  genes_present <- intersect(
    all_panel_genes,
    rownames(counts0)
  )

  m <- counts0[
    ,
    a$original_cell,
    drop=FALSE
  ]

  colnames(m) <- a$global_cell

  for(cl in clusters){

    cells <- a$global_cell[
      a$cluster == cl
    ]

    n <- length(cells)

    if(n < MIN_CELLS_PER_LIBRARY_CLUSTER){
      next
    }

    mm <- m[, cells, drop=FALSE]

    ncount <- Matrix::colSums(mm)
    nfeature <- Matrix::colSums(mm > 0)

    complexity_rows[[
      length(complexity_rows)+1L
    ]] <- data.frame(
      cluster=cl,
      library_key=lk,
      n_cells=n,
      median_nCount=median(ncount),
      median_nFeature=median(nfeature),
      stringsAsFactors=FALSE
    )

    total_umi <- sum(ncount)

    for(program in names(programs)){

      pg <- intersect(
        programs[[program]],
        genes_present
      )

      if(!length(pg)){
        next
      }

      x <- mm[
        pg,
        ,
        drop=FALSE
      ]

      detection <- Matrix::rowMeans(x > 0)

      cpm <- if(total_umi > 0){
        Matrix::rowSums(x) / total_umi * 1e6
      } else {
        rep(NA_real_, length(pg))
      }

      gene_rows[[
        length(gene_rows)+1L
      ]] <- data.frame(
        cluster=cl,
        library_key=lk,
        n_cells=n,
        program=program,
        gene=pg,
        detection_fraction=as.numeric(detection),
        CPM=as.numeric(cpm),
        stringsAsFactors=FALSE
      )
    }

    rm(mm)
  }

  rm(obj, counts0, m)
  invisible(gc())

  cat(
    sprintf(
      "READ %d/%d %s\n",
      ii,
      length(libraries),
      lk
    )
  )
}

gene_df <- do.call(rbind, gene_rows)
complexity_df <- do.call(rbind, complexity_rows)

# ============================================================
# Equal-library gene aggregation
# ============================================================

keys <- unique(
  gene_df[, c("cluster","program","gene")]
)

gene_agg <- do.call(
  rbind,
  lapply(
    seq_len(nrow(keys)),
    function(i){

      k <- keys[i, ]

      d <- gene_df[
        gene_df$cluster == k$cluster &
        gene_df$program == k$program &
        gene_df$gene == k$gene,
        ,
        drop=FALSE
      ]

      data.frame(
        cluster=k$cluster,
        program=k$program,
        gene=k$gene,
        n_eligible_libraries=nrow(d),
        median_detection_fraction=
          median(d$detection_fraction),
        median_CPM=
          median(d$CPM),
        stringsAsFactors=FALSE
      )
    }
  )
)

# ============================================================
# Panel summaries
# Descriptive thresholds only; no automatic annotation.
# ============================================================

panel_keys <- unique(
  gene_agg[, c("cluster","program")]
)

panel_summary <- do.call(
  rbind,
  lapply(
    seq_len(nrow(panel_keys)),
    function(i){

      k <- panel_keys[i, ]

      d <- gene_agg[
        gene_agg$cluster == k$cluster &
        gene_agg$program == k$program,
        ,
        drop=FALSE
      ]

      expected <- programs[[k$program]]

      data.frame(
        cluster=k$cluster,
        program=k$program,
        n_panel_genes=length(expected),
        n_genes_available=nrow(d),
        n_genes_detection_ge_0p05=
          sum(d$median_detection_fraction >= 0.05),
        n_genes_detection_ge_0p10=
          sum(d$median_detection_fraction >= 0.10),
        n_genes_detection_ge_0p25=
          sum(d$median_detection_fraction >= 0.25),
        median_gene_detection=
          median(d$median_detection_fraction),
        median_gene_CPM=
          median(d$median_CPM),
        max_gene_detection=
          max(d$median_detection_fraction),
        stringsAsFactors=FALSE
      )
    }
  )
)

# ============================================================
# Complexity equal-library aggregation
# ============================================================

complexity_summary <- do.call(
  rbind,
  lapply(
    clusters,
    function(cl){

      d <- complexity_df[
        complexity_df$cluster == cl,
        ,
        drop=FALSE
      ]

      data.frame(
        cluster=cl,
        n_eligible_libraries=nrow(d),
        median_library_median_nCount=
          median(d$median_nCount),
        median_library_median_nFeature=
          median(d$median_nFeature),
        stringsAsFactors=FALSE
      )
    }
  )
)

# ============================================================
# Output
# ============================================================

dir.create(
  file.path(
    root,
    "atlas",
    "full_atlas_annotation_v1"
  ),
  recursive=TRUE,
  showWarnings=FALSE
)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "progenitor_native_absolute_program_audit_v1__",
    tag
  )
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

write.table(
  gene_agg,
  file.path(
    out_dir,
    "progenitor_r0p4_absolute_gene_program_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  panel_summary,
  file.path(
    out_dir,
    "progenitor_r0p4_absolute_program_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  complexity_df,
  file.path(
    out_dir,
    "progenitor_r0p4_native_complexity_by_library_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  complexity_summary,
  file.path(
    out_dir,
    "progenitor_r0p4_native_complexity_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

writeLines(
  c(
    "PASS",
    "Progenitor native absolute program audit v1",
    "discovery_project=GSE216007",
    "n_discovery_cells=26962",
    "backbone=native_cluster_r0p4",
    "native SoupX RNA",
    "absolute expression audit",
    "library-equal median aggregation",
    "globins retained",
    "descriptive only",
    "no automatic annotation",
    "condition not used",
    "non-discovery projects not used",
    "cellranger_count not accessed"
  ),
  file.path(
    out_dir,
    "PROGENITOR_NATIVE_ABSOLUTE_PROGRAM_AUDIT_COMPLETE.ok"
  )
)

cat("\n===== COMPLEXITY =====\n")
print(complexity_summary, row.names=FALSE)

cat("\n===== PROGRAM SUMMARY =====\n")
print(panel_summary, row.names=FALSE)

cat("\nOUT_DIR=", out_dir, "\n", sep="")
cat("\nPASS: absolute program audit completed\n")
