#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3L){
  stop(
    "Usage: script <root> <tag> <discovery_dir>"
  )
}

root <- normalizePath(
  args[[1]],
  mustWork=TRUE
)

tag <- args[[2]]

discovery_dir <- normalizePath(
  args[[3]],
  mustWork=TRUE
)

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

N_DISCOVERY <- 8692L
N_FULL <- 12177L
MIN_CELLS_PER_LIBRARY_CLUSTER <- 5L

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "erythroid_native_absolute_program_audit_v1__",
    tag
  )
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

done_file <- file.path(
  out_dir,
  "ERYTHROID_NATIVE_ABSOLUTE_PROGRAM_AUDIT_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop(
    "Already completed: ",
    done_file
  )
}

# ============================================================
# Helpers
# ============================================================

read_tsv <- function(path){

  if(!file.exists(path)){
    stop(
      "Missing: ",
      path
    )
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

write_gz_tsv <- function(x, path){

  con <- gzfile(
    path,
    "wt"
  )

  tryCatch(
    write.table(
      x,
      con,
      sep="\t",
      quote=TRUE,
      qmethod="double",
      row.names=FALSE
    ),
    finally=close(con)
  )

  status <- system2(
    "gzip",
    c("-t", path),
    stdout=FALSE,
    stderr=FALSE
  )

  if(status != 0L){
    stop(
      "gzip integrity failure: ",
      path
    )
  }
}

# ============================================================
# Marker panels
#
# IMPORTANT:
# These panels are descriptive evidence only.
# They DO NOT automatically assign lineage/state/stage.
# ============================================================

panels <- list(

  erythroid_non_globin_identity=c(
    "GYPA",
    "GYPB",
    "ALAS2",
    "AHSP",
    "SLC4A1",
    "EPB42",
    "ANK1",
    "KLF1",
    "TFRC",
    "CA1",
    "BPGM",
    "BLVRB",
    "FECH",
    "HMBS",
    "TSPO2",
    "HEMGN"
  ),

  erythroid_globin=c(
    "HBA1",
    "HBA2",
    "HBB",
    "HBD",
    "HBM"
  ),

  erythroid_early_regulatory=c(
    "TFRC",
    "KLF1",
    "GATA1",
    "TAL1",
    "EPOR",
    "SLC2A1",
    "HMBS",
    "ANK1",
    "EPB42",
    "TSPO2"
  ),

  erythroid_terminal_hemoglobinization=c(
    "HBA1",
    "HBA2",
    "HBB",
    "HBD",
    "AHSP",
    "BPGM",
    "CA1",
    "SLC4A1",
    "PRDX2",
    "BLVRB",
    "UBE2O",
    "EIF2AK1"
  ),

  platelet_specific=c(
    "PPBP",
    "PF4",
    "GP1BA",
    "GP1BB",
    "GP9",
    "ITGA2B",
    "ITGB3",
    "TUBB1",
    "TREML1",
    "CLEC1B",
    "P2RY12",
    "RGS18"
  ),

  myeloid_specific=c(
    "LYZ",
    "S100A8",
    "S100A9",
    "S100A12",
    "FCN1",
    "VCAN",
    "LST1",
    "TYROBP",
    "CTSS",
    "CTSD",
    "FCER1G"
  ),

  T_NK_specific=c(
    "CD3D",
    "CD3E",
    "CD3G",
    "TRAC",
    "TRBC1",
    "TRBC2",
    "NKG7",
    "GNLY",
    "CCL5",
    "GIMAP4"
  ),

  B_specific=c(
    "CD79A",
    "MS4A1",
    "CD37",
    "CD22",
    "CD83"
  ),

  cycling=c(
    "MKI67",
    "TOP2A",
    "CENPF",
    "UBE2C",
    "STMN1",
    "HMGB2",
    "TUBA1B"
  ),

  immediate_early=c(
    "FOS",
    "FOSB",
    "JUN",
    "JUNB",
    "DUSP1",
    "IER2",
    "IER3",
    "ATF3"
  ),

  interferon=c(
    "IFI6",
    "IFI27",
    "IFI44",
    "IFI44L",
    "IFIT1",
    "IFIT2",
    "IFIT3",
    "ISG15",
    "MX1",
    "MX2",
    "OAS1",
    "OAS2",
    "OAS3",
    "OASL",
    "RSAD2",
    "HERC5",
    "XAF1",
    "LY6E"
  )
)

panel_map <- do.call(
  rbind,
  lapply(
    names(panels),
    function(p){
      data.frame(
        panel=p,
        gene=panels[[p]],
        stringsAsFactors=FALSE
      )
    }
  )
)

marker_genes <- unique(
  panel_map$gene
)

# ============================================================
# Inputs
# ============================================================

assignment_file <- file.path(
  discovery_dir,
  "erythroid_balanced_discovery_cluster_assignments_v1.tsv.gz"
)

transfer_dir <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "annotation_transfer_by_library"
)

lib_file <- file.path(
  root,
  "pre_integration",
  "pilot_unintegrated_large_v1",
  "resolved_library_table.tsv"
)

a <- read_tsv(
  assignment_file
)

stopifnot(
  nrow(a) == N_DISCOVERY,
  !anyDuplicated(a$global_cell)
)

required_a <- c(
  "global_cell",
  "project_id",
  "library_key",
  "cluster_r0p4"
)

stopifnot(
  all(
    required_a %in%
      names(a)
  )
)

for(nm in required_a){
  a[[nm]] <- as.character(
    a[[nm]]
  )
}

# ============================================================
# Full metadata -> original cell IDs
# ============================================================

transfer_files <- list.files(
  transfer_dir,
  pattern="\\.tsv\\.gz$",
  full.names=TRUE
)

meta_list <- list()

for(f in transfer_files){

  z <- read_tsv(f)

  req <- c(
    "global_cell",
    "original_cell",
    "project_id",
    "library_key",
    "integration_compartment_primary_v1"
  )

  stopifnot(
    all(req %in% names(z))
  )

  z <- z[
    as.character(
      z$integration_compartment_primary_v1
    ) ==
      "Erythroid",
    req,
    drop=FALSE
  ]

  if(nrow(z)){
    meta_list[[
      length(meta_list)+1L
    ]] <- z
  }
}

full_meta <- do.call(
  rbind,
  meta_list
)

rm(meta_list)

for(nm in c(
  "global_cell",
  "original_cell",
  "project_id",
  "library_key"
)){
  full_meta[[nm]] <- as.character(
    full_meta[[nm]]
  )
}

stopifnot(
  nrow(full_meta) == N_FULL,
  !anyDuplicated(
    full_meta$global_cell
  )
)

mi <- match(
  a$global_cell,
  full_meta$global_cell
)

stopifnot(
  !anyNA(mi)
)

a$original_cell <-
  full_meta$original_cell[mi]

stopifnot(
  a$project_id ==
    full_meta$project_id[mi],

  a$library_key ==
    full_meta$library_key[mi]
)

# ============================================================
# Library lookup
# ============================================================

lib <- read_tsv(
  lib_file
)

stopifnot(
  all(
    c(
      "project_id",
      "library_key",
      "final_rds_resolved"
    ) %in%
      names(lib)
  )
)

lib$project_id <- as.character(
  lib$project_id
)

lib$library_key <- as.character(
  lib$library_key
)

lib$key <- paste(
  lib$project_id,
  lib$library_key,
  sep="|||"
)

# ============================================================
# Native absolute measurements
#
# For every library x cluster:
#   gene count
#   total UMI denominator
#   number of cells detecting gene
#
# Then combine libraries WITHIN project.
# Projects receive equal weight in final summary.
# ============================================================

library_keys <- unique(
  paste(
    a$project_id,
    a$library_key,
    sep="|||"
  )
)

gene_lib_list <- list()
cell_qc_list <- list()

marker_genes_present <- NULL
marker_genes_missing <- NULL

for(ii in seq_along(
  library_keys
)){

  key <- library_keys[[ii]]

  sp <- strsplit(
    key,
    "\\|\\|\\|"
  )[[1]]

  p <- sp[[1]]
  lk <- sp[[2]]

  d <- a[
    a$project_id == p &
      a$library_key == lk,
    ,
    drop=FALSE
  ]

  li <- match(
    key,
    lib$key
  )

  if(is.na(li)){
    stop(
      "Library lookup failed: ",
      key
    )
  }

  source_rds <- as.character(
    lib$final_rds_resolved[[li]]
  )

  if(!file.exists(source_rds)){
    stop(
      "Missing source RDS: ",
      source_rds
    )
  }

  obj0 <- readRDS(
    source_rds
  )

  counts0 <- get_rna_counts(
    obj0
  )

  if(is.null(
    marker_genes_present
  )){

    marker_genes_present <-
      marker_genes[
        marker_genes %in%
          rownames(counts0)
      ]

    marker_genes_missing <-
      setdiff(
        marker_genes,
        marker_genes_present
      )

    if(length(marker_genes_present) < 50L){
      stop(
        "Too few audit marker genes present: ",
        length(marker_genes_present)
      )
    }

  } else {

    if(!all(
      marker_genes_present %in%
        rownames(counts0)
    )){
      stop(
        "Marker feature mismatch: ",
        key
      )
    }
  }

  if(!all(
    d$original_cell %in%
      colnames(counts0)
  )){
    stop(
      "Missing discovery cells: ",
      key
    )
  }

  # ------------------------------------------
  # Per-cell RNA complexity
  # ------------------------------------------

  m_all <- counts0[
    ,
    d$original_cell,
    drop=FALSE
  ]

  ncount <- Matrix::colSums(
    m_all
  )

  nfeature <- Matrix::colSums(
    m_all > 0
  )

  cell_qc_list[[
    length(cell_qc_list)+1L
  ]] <- data.frame(
    project_id=p,
    library_key=lk,
    cluster=d$cluster_r0p4,
    global_cell=d$global_cell,
    nCount_native=as.numeric(ncount),
    nFeature_native=as.numeric(nfeature),
    stringsAsFactors=FALSE
  )

  rm(
    m_all,
    ncount,
    nfeature
  )

  # ------------------------------------------
  # Library x cluster absolute expression
  # ------------------------------------------

  for(cl in sort(
    unique(
      d$cluster_r0p4
    )
  )){

    dc <- d[
      d$cluster_r0p4 == cl,
      ,
      drop=FALSE
    ]

    n_cells <- nrow(dc)

    if(n_cells <
       MIN_CELLS_PER_LIBRARY_CLUSTER){
      next
    }

    cells <- dc$original_cell

    total_umis <- sum(
      Matrix::colSums(
        counts0[
          ,
          cells,
          drop=FALSE
        ]
      )
    )

    m <- counts0[
      marker_genes_present,
      cells,
      drop=FALSE
    ]

    gene_count <- Matrix::rowSums(
      m
    )

    detected <- Matrix::rowSums(
      m > 0
    )

    gene_lib_list[[
      length(gene_lib_list)+1L
    ]] <- data.frame(
      project_id=p,
      library_key=lk,
      cluster=cl,
      gene=marker_genes_present,
      n_cells=n_cells,
      gene_count_sum=
        as.numeric(
          gene_count[
            marker_genes_present
          ]
        ),
      total_umis=
        as.numeric(
          total_umis
        ),
      detected_cells=
        as.numeric(
          detected[
            marker_genes_present
          ]
        ),
      stringsAsFactors=FALSE
    )

    rm(
      m,
      gene_count,
      detected
    )
  }

  rm(
    obj0,
    counts0
  )

  invisible(gc())

  cat(
    sprintf(
      "READ %3d/%3d %s\n",
      ii,
      length(library_keys),
      key
    )
  )
}

gene_lib <- do.call(
  rbind,
  gene_lib_list
)

cell_qc <- do.call(
  rbind,
  cell_qc_list
)

rm(
  gene_lib_list,
  cell_qc_list
)

stopifnot(
  nrow(cell_qc) == N_DISCOVERY,
  !anyDuplicated(
    cell_qc$global_cell
  )
)

# ============================================================
# Aggregate libraries within project
# ============================================================

project_gene <- aggregate(
  cbind(
    gene_count_sum,
    total_umis,
    detected_cells,
    n_cells
  ) ~
    project_id +
    cluster +
    gene,
  data=gene_lib,
  FUN=sum
)

project_gene$log2CPM <-
  log2(
    1 +
      1e6 *
      project_gene$gene_count_sum /
      pmax(
        project_gene$total_umis,
        1
      )
  )

project_gene$detection_fraction <-
  project_gene$detected_cells /
  pmax(
    project_gene$n_cells,
    1
  )

# ============================================================
# Equal-project cluster/gene summary
# ============================================================

split_gene <- split(
  project_gene,
  paste(
    project_gene$cluster,
    project_gene$gene,
    sep="|||"
  )
)

gene_summary <- do.call(
  rbind,
  lapply(
    split_gene,
    function(d){

      data.frame(
        cluster=
          as.character(
            d$cluster[[1]]
          ),

        gene=
          as.character(
            d$gene[[1]]
          ),

        n_projects_evaluable=
          length(
            unique(
              d$project_id
            )
          ),

        mean_project_log2CPM=
          mean(
            d$log2CPM
          ),

        median_project_log2CPM=
          median(
            d$log2CPM
          ),

        mean_project_detection=
          mean(
            d$detection_fraction
          ),

        median_project_detection=
          median(
            d$detection_fraction
          ),

        fraction_projects_detection_ge_0p05=
          mean(
            d$detection_fraction >=
              0.05
          ),

        fraction_projects_detection_ge_0p10=
          mean(
            d$detection_fraction >=
              0.10
          ),

        fraction_projects_detection_ge_0p25=
          mean(
            d$detection_fraction >=
              0.25
          ),

        stringsAsFactors=FALSE
      )
    }
  )
)

rownames(gene_summary) <- NULL

# ============================================================
# Add panel membership
# ============================================================

panel_summary_source <- merge(
  gene_summary,
  panel_map,
  by="gene",
  all.x=FALSE,
  all.y=FALSE
)

sp <- split(
  panel_summary_source,
  paste(
    panel_summary_source$cluster,
    panel_summary_source$panel,
    sep="|||"
  )
)

program_summary <- do.call(
  rbind,
  lapply(
    sp,
    function(d){

      reproducible <-
        d$mean_project_detection >=
          0.10 &
        d$fraction_projects_detection_ge_0p10 >=
          0.50

      data.frame(
        cluster=
          as.character(
            d$cluster[[1]]
          ),

        panel=
          as.character(
            d$panel[[1]]
          ),

        n_genes_evaluable=
          nrow(d),

        n_genes_reproducibly_detected=
          sum(
            reproducible
          ),

        reproducibly_detected_genes=
          paste(
            d$gene[
              reproducible
            ],
            collapse=","
          ),

        median_gene_mean_project_detection=
          median(
            d$mean_project_detection
          ),

        median_gene_mean_project_log2CPM=
          median(
            d$mean_project_log2CPM
          ),

        stringsAsFactors=FALSE
      )
    }
  )
)

rownames(program_summary) <- NULL

# ============================================================
# Project-balanced RNA complexity
# ============================================================

project_qc <- aggregate(
  cbind(
    nCount_native,
    nFeature_native
  ) ~
    project_id +
    cluster,
  data=cell_qc,
  FUN=median
)

qc_split <- split(
  project_qc,
  project_qc$cluster
)

qc_summary <- do.call(
  rbind,
  lapply(
    qc_split,
    function(d){

      data.frame(
        cluster=
          as.character(
            d$cluster[[1]]
          ),

        n_projects_evaluable=
          length(
            unique(
              d$project_id
            )
          ),

        median_project_median_nCount=
          median(
            d$nCount_native
          ),

        median_project_median_nFeature=
          median(
            d$nFeature_native
          ),

        min_project_median_nFeature=
          min(
            d$nFeature_native
          ),

        max_project_median_nFeature=
          max(
            d$nFeature_native
          ),

        stringsAsFactors=FALSE
      )
    }
  )
)

rownames(qc_summary) <- NULL

# Add discovery counts.
cluster_counts <- as.data.frame(
  table(
    a$cluster_r0p4
  ),
  stringsAsFactors=FALSE
)

names(cluster_counts) <- c(
  "cluster",
  "n_discovery_cells"
)

cluster_counts$cluster <-
  as.character(
    cluster_counts$cluster
  )

qc_summary <- merge(
  cluster_counts,
  qc_summary,
  by="cluster",
  all.x=TRUE
)

# ============================================================
# Key-gene compact output
# ============================================================

key_genes <- c(
  "GYPA",
  "ALAS2",
  "AHSP",
  "KLF1",
  "TFRC",
  "GATA1",
  "SLC2A1",
  "ANK1",
  "EPB42",
  "SLC4A1",
  "CA1",
  "BPGM",
  "BLVRB",
  "HBA1",
  "HBA2",
  "HBB",
  "HBD",
  "PPBP",
  "PF4",
  "ITGA2B",
  "LYZ",
  "S100A8",
  "S100A9",
  "FCER1G",
  "FOS",
  "MKI67",
  "TOP2A",
  "IFI6",
  "ISG15"
)

key_summary <- gene_summary[
  gene_summary$gene %in%
    key_genes,
  ,
  drop=FALSE
]

key_summary <- key_summary[
  order(
    as.integer(
      key_summary$cluster
    ),
    match(
      key_summary$gene,
      key_genes
    )
  ),
  ,
  drop=FALSE
]

# ============================================================
# Output
# ============================================================

write.table(
  panel_map,
  file.path(
    out_dir,
    "erythroid_absolute_marker_panel_definition_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  gene_summary,
  file.path(
    out_dir,
    "erythroid_r0p4_absolute_gene_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  program_summary,
  file.path(
    out_dir,
    "erythroid_r0p4_absolute_program_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  qc_summary,
  file.path(
    out_dir,
    "erythroid_r0p4_native_complexity_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  key_summary,
  file.path(
    out_dir,
    "erythroid_r0p4_key_gene_absolute_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write_gz_tsv(
  project_gene,
  file.path(
    out_dir,
    "erythroid_r0p4_project_absolute_gene_metrics_v1.tsv.gz"
  )
)

# ============================================================
# Re-read / integrity checks
# ============================================================

chk <- read_tsv(
  file.path(
    out_dir,
    "erythroid_r0p4_project_absolute_gene_metrics_v1.tsv.gz"
  )
)

stopifnot(
  nrow(chk) ==
    nrow(project_gene)
)

rm(chk)

# ============================================================
# Completion marker LAST
# ============================================================

writeLines(
  c(
    "PASS",
    "Erythroid native absolute program audit v1",
    "discovery_cells=8692",
    "backbone_resolution=r0.4",
    "",
    "native SoupX-corrected RNA",
    "absolute expression evidence",
    "library data combined within project",
    "projects receive equal weight in cluster summaries",
    "",
    "globin genes retained as biological features",
    "non-globin erythroid identity scored separately",
    "RNA complexity summarized project-balanced",
    "",
    "NO automatic lineage assignment",
    "NO automatic maturation assignment",
    "NO automatic state assignment",
    "NO relative module-delta-only classification",
    "integrated assay NOT used as biological evidence",
    "cellranger_count not accessed",
    "gzip_integrity=PASS",
    "re_read_check=PASS"
  ),
  done_file
)

cat(
  "\n===== RNA COMPLEXITY =====\n"
)

print(
  qc_summary[
    order(
      as.integer(
        qc_summary$cluster
      )
    ),
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat(
  "\n===== PROGRAM SUMMARY =====\n"
)

print(
  program_summary[
    order(
      as.integer(
        program_summary$cluster
      ),
      program_summary$panel
    ),
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Erythroid native absolute audit completed\n"
)
