#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
})

in_rds <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_native_subclustering_v1_primary9997",
  "b_plasma_pilot_native_subclustering_v1.rds"
)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_annotation_diagnostic_v1"
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

stopifnot(file.exists(in_rds))

bp <- readRDS(in_rds)

stopifnot(
  nrow(bp) == 38606,
  ncol(bp) == 9997
)

md <- bp[[]]

stopifnot(
  "cluster_res_0p4" %in% names(md),
  "project_id" %in% names(md)
)

cl <- as.character(
  md$cluster_res_0p4
)

counts <- LayerData(
  bp[["RNA"]],
  layer="counts"
)

genes <- rownames(counts)

# ------------------------------------------------------------
# Flags
# ------------------------------------------------------------

technical <- grepl(
  paste0(
    "^MT-|",
    "^RPL[0-9]|",
    "^RPS[0-9]|",
    "^HBA[12]$|",
    "^HBB$|",
    "^HBD$|",
    "^HBG[12]$|",
    "^MALAT1$|",
    "^NEAT1$"
  ),
  genes
)

rearranged_ig <- grepl(
  "^IGHV|^IGKV|^IGLV|^IGHJ[0-9]|^IGKJ[0-9]|^IGLJ[0-9]",
  genes
)

# ------------------------------------------------------------
# Global marker diagnostic for previously unresolved clusters
# ------------------------------------------------------------

targets <- c(
  "6","8","9","10","11","12","14"
)

all_markers <- list()

for(k in targets){

  ii <- which(cl == k)
  jj <- which(cl != k)

  cin <- Matrix::rowSums(
    counts[,ii,drop=FALSE]
  )

  cout <- Matrix::rowSums(
    counts[,jj,drop=FALSE]
  )

  din <- Matrix::rowSums(
    counts[,ii,drop=FALSE] > 0
  )

  dout <- Matrix::rowSums(
    counts[,jj,drop=FALSE] > 0
  )

  cpm_in <- cin / sum(cin) * 1e6
  cpm_out <- cout / sum(cout) * 1e6

  d <- data.frame(
    cluster=k,
    gene=genes,
    log2FC=
      log2(cpm_in + 0.1) -
      log2(cpm_out + 0.1),
    pct_in=din / length(ii),
    pct_out=dout / length(jj),
    delta_pct=
      din / length(ii) -
      dout / length(jj),
    technical=technical,
    rearranged_ig=rearranged_ig,
    stringsAsFactors=FALSE
  )

  d <- d[
    !d$technical &
    !d$rearranged_ig &
    d$pct_in >= 0.10 &
    d$delta_pct >= 0.05 &
    d$log2FC >= 0.5,
    ,
    drop=FALSE
  ]

  d <- d[
    order(
      -d$delta_pct,
      -d$log2FC
    ),
    ,
    drop=FALSE
  ]

  all_markers[[
    length(all_markers)+1L
  ]] <- d

  cat(
    "\n============================================================\n",
    "CLUSTER ", k,
    "   n=", length(ii), "\n",
    "============================================================\n",
    sep=""
  )

  print(
    head(d,30),
    row.names=FALSE
  )

  cat("\nProjects:\n")

  print(
    sort(
      table(md$project_id[ii]),
      decreasing=TRUE
    )
  )
}

all_markers <- do.call(
  rbind,
  all_markers
)

write.table(
  all_markers,
  file=file.path(
    out_dir,
    "remaining_clusters_global_markers_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ------------------------------------------------------------
# Canonical targeted panel across ALL r0.4 clusters
# ------------------------------------------------------------

panel <- unique(c(

  # generic B lineage
  "CD19","MS4A1","CD79A","CD79B",
  "CD37","CD74","HLA-DRA","CD22",
  "CD83","CD86",

  # naive / transitional
  "TCL1A","IGHD","IGHM","FCER2",
  "IL4R","BACH2","MME","CD38",
  "VPREB3","IGLL5",

  # memory
  "CD27","TNFRSF13B","AIM2",
  "GPR183","BANK1",

  # atypical / ABC
  "ITGAX","FCRL5","FCRL3",
  "FCRL4","TBX21","ZEB2",
  "TNFRSF1B",

  # switched isotypes
  "IGHA1","IGHA2",
  "IGHG1","IGHG2","IGHG3","IGHG4",

  # antibody-secreting lineage
  "MZB1","JCHAIN","XBP1","PRDM1",
  "TNFRSF17","SDC1","FKBP11",
  "DERL3","SEC11C","IRF4",

  # interferon
  "ISG15","IFI6","IFI44L","MX1",
  "IFIT1","IFIT3","OAS1","OASL",

  # cycling
  "MKI67","TOP2A","STMN1",
  "TYMS","UBE2C",

  # platelet
  "PPBP","PF4","GP1BB","GP9",
  "ITGA2B","TUBB1","TREML1",

  # T/NK
  "CD3D","CD3E","CD3G","TRAC",
  "CD247","CD2","IL7R","TCF7",
  "NKG7","GNLY",

  # myeloid
  "LYZ","S100A8","S100A9",
  "FCN1","CTSS","LILRB1"
))

panel <- panel[
  panel %in% genes
]

cluster_levels <- sort(
  unique(cl),
  method="radix"
)

cluster_levels <- cluster_levels[
  order(
    suppressWarnings(
      as.integer(cluster_levels)
    )
  )
]

panel_out <- list()

for(k in cluster_levels){

  ii <- which(cl == k)

  mat <- counts[
    panel,
    ii,
    drop=FALSE
  ]

  total <- Matrix::rowSums(mat)
  detected <- Matrix::rowSums(mat > 0)

  libsize <- sum(
    Matrix::colSums(
      counts[,ii,drop=FALSE]
    )
  )

  panel_out[[
    length(panel_out)+1L
  ]] <- data.frame(
    cluster=k,
    n_cells=length(ii),
    gene=panel,
    pct_detected=
      detected / length(ii),
    pseudobulk_CPM=
      total / libsize * 1e6,
    stringsAsFactors=FALSE
  )
}

panel_out <- do.call(
  rbind,
  panel_out
)

write.table(
  panel_out,
  file=file.path(
    out_dir,
    "all_clusters_canonical_marker_panel_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# Print compact genes with useful detection
cat(
  "\n\n===== CANONICAL MARKER PANEL =====\n"
)

for(k in cluster_levels){

  d <- panel_out[
    panel_out$cluster == k &
    panel_out$pct_detected >= 0.10,
    ,
    drop=FALSE
  ]

  d <- d[
    order(
      -d$pct_detected,
      -d$pseudobulk_CPM
    ),
    ,
    drop=FALSE
  ]

  cat(
    "\n----------------------------------------\n",
    "CLUSTER ", k, "\n",
    "----------------------------------------\n",
    sep=""
  )

  print(
    d,
    row.names=FALSE
  )
}

writeLines(
  c(
    "PASS",
    "B/plasma remaining-cluster and canonical-marker diagnostic v1",
    "resolution=0.4",
    "n_cells=9997",
    "global markers are diagnostic only",
    "canonical panel evaluated across all 15 clusters"
  ),
  file.path(
    out_dir,
    "B_PLASMA_DIAGNOSTIC_COMPLETE.ok"
  )
)

cat(
  "\nPASS: B/plasma annotation diagnostic completed\n"
)

