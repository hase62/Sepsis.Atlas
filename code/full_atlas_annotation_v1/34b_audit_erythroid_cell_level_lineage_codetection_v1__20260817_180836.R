#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3L){
  stop("Usage: script <root> <tag> <discovery_dir>")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]
discovery_dir <- normalizePath(args[[3]], mustWork=TRUE)

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

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "erythroid_cell_level_lineage_codetection_v1__",
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
  "ERYTHROID_CELL_LEVEL_LINEAGE_CODETECTION_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop("Already completed: ", done_file)
}

read_tsv <- function(path){

  if(!file.exists(path)){
    stop("Missing: ", path)
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

  con <- gzfile(path, "wt")

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
    stop("gzip integrity failure: ", path)
  }
}

# ============================================================
# Specific lineage panels.
#
# Globins are kept SEPARATE from erythroid non-globin identity.
#
# Threshold flags below are descriptive co-detection flags only.
# They are NOT annotation rules.
# ============================================================

panels <- list(

  erythroid_non_globin=c(
    "ALAS2","AHSP","SLC4A1","CA1","BPGM",
    "BLVRB","EPB42","ANK1","KLF1","GYPA",
    "GYPB","FECH","HEMGN","TFRC","TSPO2"
  ),

  globin=c(
    "HBA1","HBA2","HBB","HBD","HBM"
  ),

  myeloid=c(
    "LYZ","S100A8","S100A9","S100A12","FCN1",
    "VCAN","LST1","TYROBP","CTSS","FCER1G"
  ),

  platelet=c(
    "PPBP","PF4","GP1BA","GP1BB","ITGA2B",
    "TUBB1","TREML1","P2RY12","RGS18"
  ),

  T_NK=c(
    "CD3D","CD3E","CD3G","TRAC","TRBC1",
    "TRBC2","NKG7","GNLY","GIMAP4"
  ),

  B=c(
    "CD79A","MS4A1","CD22","CD19","CD37"
  )
)

all_panel_genes <- unique(
  unlist(
    panels,
    use.names=FALSE
  )
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

for(nm in c(
  "global_cell",
  "project_id",
  "library_key",
  "cluster_r0p4"
)){
  a[[nm]] <- as.character(a[[nm]])
}

# ============================================================
# Recover original cell IDs
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
    ) == "Erythroid",
    req,
    drop=FALSE
  ]

  if(nrow(z)){
    meta_list[[length(meta_list)+1L]] <- z
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
  full_meta[[nm]] <- as.character(full_meta[[nm]])
}

stopifnot(
  nrow(full_meta) == N_FULL,
  !anyDuplicated(full_meta$global_cell)
)

mi <- match(
  a$global_cell,
  full_meta$global_cell
)

stopifnot(!anyNA(mi))

a$original_cell <- full_meta$original_cell[mi]

stopifnot(
  a$project_id == full_meta$project_id[mi],
  a$library_key == full_meta$library_key[mi]
)

# ============================================================
# Library lookup
# ============================================================

lib <- read_tsv(lib_file)

lib$project_id <- as.character(lib$project_id)
lib$library_key <- as.character(lib$library_key)

lib$key <- paste(
  lib$project_id,
  lib$library_key,
  sep="|||"
)

library_keys <- unique(
  paste(
    a$project_id,
    a$library_key,
    sep="|||"
  )
)

# ============================================================
# Cell-level panel detection
# ============================================================

cell_list <- list()

for(ii in seq_along(library_keys)){

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
    stop("Library lookup failed: ", key)
  }

  source_rds <- as.character(
    lib$final_rds_resolved[[li]]
  )

  if(!file.exists(source_rds)){
    stop("Missing source RDS: ", source_rds)
  }

  obj0 <- readRDS(source_rds)

  counts0 <- get_rna_counts(obj0)

  cells <- d$original_cell

  if(!all(cells %in% colnames(counts0))){
    stop("Missing discovery cells: ", key)
  }

  genes_present <- intersect(
    all_panel_genes,
    rownames(counts0)
  )

  m <- counts0[
    genes_present,
    cells,
    drop=FALSE
  ]

  total_count <- Matrix::colSums(
    counts0[
      ,
      cells,
      drop=FALSE
    ]
  )

  total_feature <- Matrix::colSums(
    counts0[
      ,
      cells,
      drop=FALSE
    ] > 0
  )

  out <- data.frame(
    global_cell=d$global_cell,
    project_id=p,
    library_key=lk,
    cluster=d$cluster_r0p4,
    nCount_native=as.numeric(total_count),
    nFeature_native=as.numeric(total_feature),
    stringsAsFactors=FALSE
  )

  for(panel in names(panels)){

    g <- intersect(
      panels[[panel]],
      rownames(m)
    )

    if(!length(g)){
      stop("No genes available for panel: ", panel)
    }

    mm <- m[
      g,
      ,
      drop=FALSE
    ]

    out[[paste0(panel, "_n_detected")]] <-
      as.numeric(
        Matrix::colSums(mm > 0)
      )

    out[[paste0(panel, "_umis")]] <-
      as.numeric(
        Matrix::colSums(mm)
      )

    out[[paste0(panel, "_umi_fraction")]] <-
      out[[paste0(panel, "_umis")]] /
      pmax(
        out$nCount_native,
        1
      )
  }

  # ----------------------------------------------------------
  # Descriptive flags only.
  # ----------------------------------------------------------

  out$ery3 <-
    out$erythroid_non_globin_n_detected >= 3L

  out$myeloid3 <-
    out$myeloid_n_detected >= 3L

  out$platelet2 <-
    out$platelet_n_detected >= 2L

  out$T_NK2 <-
    out$T_NK_n_detected >= 2L

  out$B2 <-
    out$B_n_detected >= 2L

  out$any_nonery_lineage <-
    out$myeloid3 |
    out$platelet2 |
    out$T_NK2 |
    out$B2

  out$ery_and_nonery <-
    out$ery3 &
    out$any_nonery_lineage

  cell_list[[
    length(cell_list)+1L
  ]] <- out

  rm(
    obj0,
    counts0,
    m,
    total_count,
    total_feature,
    out
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

cell <- do.call(
  rbind,
  cell_list
)

rm(cell_list)

stopifnot(
  nrow(cell) == N_DISCOVERY,
  !anyDuplicated(cell$global_cell)
)

# ============================================================
# Cluster summary
# ============================================================

summarize_cells <- function(d){

  data.frame(
    n_cells=nrow(d),

    median_nCount=
      median(d$nCount_native),

    median_nFeature=
      median(d$nFeature_native),

    fraction_ery3=
      mean(d$ery3),

    fraction_myeloid3=
      mean(d$myeloid3),

    fraction_platelet2=
      mean(d$platelet2),

    fraction_T_NK2=
      mean(d$T_NK2),

    fraction_B2=
      mean(d$B2),

    fraction_any_nonery_lineage=
      mean(d$any_nonery_lineage),

    fraction_ery_and_nonery=
      mean(d$ery_and_nonery),

    fraction_ery3_myeloid3=
      mean(
        d$ery3 &
        d$myeloid3
      ),

    fraction_ery3_platelet2=
      mean(
        d$ery3 &
        d$platelet2
      ),

    fraction_ery3_T_NK2=
      mean(
        d$ery3 &
        d$T_NK2
      ),

    fraction_ery3_B2=
      mean(
        d$ery3 &
        d$B2
      ),

    median_ery_n_detected=
      median(
        d$erythroid_non_globin_n_detected
      ),

    median_myeloid_n_detected=
      median(
        d$myeloid_n_detected
      ),

    median_platelet_n_detected=
      median(
        d$platelet_n_detected
      ),

    median_T_NK_n_detected=
      median(
        d$T_NK_n_detected
      ),

    median_B_n_detected=
      median(
        d$B_n_detected
      ),

    median_ery_umi_fraction=
      median(
        d$erythroid_non_globin_umi_fraction
      ),

    median_globin_umi_fraction=
      median(
        d$globin_umi_fraction
      ),

    median_myeloid_umi_fraction=
      median(
        d$myeloid_umi_fraction
      ),

    median_platelet_umi_fraction=
      median(
        d$platelet_umi_fraction
      ),

    median_T_NK_umi_fraction=
      median(
        d$T_NK_umi_fraction
      ),

    median_B_umi_fraction=
      median(
        d$B_umi_fraction
      ),

    stringsAsFactors=FALSE
  )
}

cluster_split <- split(
  cell,
  cell$cluster
)

cluster_summary <- do.call(
  rbind,
  lapply(
    names(cluster_split),
    function(cl){

      x <- summarize_cells(
        cluster_split[[cl]]
      )

      data.frame(
        cluster=cl,
        x,
        stringsAsFactors=FALSE
      )
    }
  )
)

rownames(cluster_summary) <- NULL

cluster_summary <- cluster_summary[
  order(
    as.integer(cluster_summary$cluster)
  ),
  ,
  drop=FALSE
]

# ============================================================
# Project x cluster summary
# ============================================================

pc_key <- paste(
  cell$project_id,
  cell$cluster,
  sep="|||"
)

pc_split <- split(
  cell,
  pc_key
)

project_cluster_summary <- do.call(
  rbind,
  lapply(
    pc_split,
    function(d){

      x <- summarize_cells(d)

      data.frame(
        project_id=
          as.character(
            d$project_id[[1]]
          ),

        cluster=
          as.character(
            d$cluster[[1]]
          ),

        x,
        stringsAsFactors=FALSE
      )
    }
  )
)

rownames(project_cluster_summary) <- NULL

project_cluster_summary <-
  project_cluster_summary[
    order(
      as.integer(
        project_cluster_summary$cluster
      ),
      project_cluster_summary$project_id
    ),
    ,
    drop=FALSE
  ]

# ============================================================
# Compact suspicious-cluster output
# ============================================================

target_clusters <- c(
  "2","5","7","8","9"
)

target <- cell[
  cell$cluster %in%
    target_clusters,
  ,
  drop=FALSE
]

target_summary <- cluster_summary[
  cluster_summary$cluster %in%
    target_clusters,
  ,
  drop=FALSE
]

# ============================================================
# Write
# ============================================================

write.table(
  cluster_summary,
  file.path(
    out_dir,
    "erythroid_r0p4_cell_lineage_codetection_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  project_cluster_summary,
  file.path(
    out_dir,
    "erythroid_r0p4_project_cell_lineage_codetection_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  target_summary,
  file.path(
    out_dir,
    "erythroid_r0p4_suspicious_cluster_codetection_compact_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write_gz_tsv(
  target,
  file.path(
    out_dir,
    "erythroid_r0p4_suspicious_cluster_cell_metrics_v1.tsv.gz"
  )
)

chk <- read_tsv(
  file.path(
    out_dir,
    "erythroid_r0p4_suspicious_cluster_cell_metrics_v1.tsv.gz"
  )
)

stopifnot(
  nrow(chk) == nrow(target),
  identical(
    as.character(chk$global_cell),
    as.character(target$global_cell)
  )
)

writeLines(
  c(
    "PASS",
    "Erythroid cell-level lineage co-detection audit v1",
    "discovery_cells=8692",
    "backbone_resolution=r0.4",
    "",
    "native SoupX-corrected RNA",
    "globin and non-globin erythroid signals separated",
    "",
    "descriptive flags only:",
    "ery3 >=3 non-globin erythroid genes detected",
    "myeloid3 >=3 myeloid genes detected",
    "platelet2 >=2 platelet genes detected",
    "T_NK2 >=2 T/NK genes detected",
    "B2 >=2 B genes detected",
    "",
    "NO automatic annotation",
    "NO cluster reassignment",
    "cellranger_count not accessed",
    "gzip_integrity=PASS",
    "re_read_check=PASS"
  ),
  done_file
)

cat(
  "\n===== ALL CLUSTERS =====\n"
)

print(
  cluster_summary,
  row.names=FALSE
)

cat(
  "\n===== SUSPICIOUS CLUSTERS =====\n"
)

print(
  target_summary,
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Erythroid cell-level co-detection audit completed\n"
)
