#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

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

MIN_CELLS <- 10L

ann_file <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_full_primary_transfer_v1__repaired_20260814_1430",
  "b_plasma_full_primary_annotation_repaired_v1.rds"
)

lib_file <- file.path(
  root,
  "pre_integration",
  "pilot_unintegrated_large_v1",
  "resolved_library_table.tsv"
)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_full_native_module_audit_v1__20260814_1448"
)

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

done_file <- file.path(
  out_dir,
  "B_PLASMA_FULL_NATIVE_MODULE_AUDIT_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop("Audit already completed: ", done_file)
}

# ============================================================
# Canonical modules
# ============================================================

modules <- list(

  naive_identity=c(
    "TCL1A","IGHD","IGHM","FCER2",
    "IL4R","BACH2","FCRL1","VPREB3"
  ),

  memory_identity=c(
    "CD27","TNFRSF13B","AIM2",
    "GPR183","CD82","ITGB1","POU2AF1"
  ),

  class_switched_identity=c(
    "IGHA1","IGHA2",
    "IGHG1","IGHG2","IGHG3","IGHG4"
  ),

  interferon=c(
    "IFI6","IFI27","IFI35","IFI44","IFI44L",
    "IFIT1","IFIT2","IFIT3",
    "IFITM1","IFITM2","IFITM3",
    "ISG15","MX1","MX2",
    "OAS1","OAS2","OAS3","OASL",
    "XAF1","BST2","LY6E","EPSTI1"
  ),

  early_activation=c(
    "CD69","NFKBIA",
    "EGR1","EGR2",
    "FOS","FOSB","JUN","JUNB",
    "DUSP1","DUSP2",
    "NR4A1","NR4A2"
  ),

  activation_stress=c(
    "ATF3","PPP1R15A",
    "FOS","FOSB","JUN","JUNB",
    "DUSP1","DUSP2","NFKBIA"
  ),

  atypical_activation=c(
    "ITGAX","FCRL5","TBX21","CD86"
  ),

  secretory=c(
    "XBP1","PRDM1","JCHAIN","MZB1",
    "FKBP11","DERL3","SEC11C",
    "ELL2","IRF4","TNFRSF17","SDC1"
  ),

  platelet=c(
    "PPBP","PF4","PF4V1",
    "GP1BB","GP9","ITGA2B",
    "TUBB1","TREML1","MPIG6B",
    "CAVIN2","GNG11"
  ),

  T_cell=c(
    "CD3D","CD3E","CD3G",
    "TRAC","TRBC1","TRBC2",
    "CD247","CD2","LCK",
    "IL7R","TCF7","GIMAP7"
  )
)

genes <- unique(unlist(modules))

contrasts <- list(
  T01_naive=list(
    target="T01",
    baseline=c("T02","T03"),
    module="naive_identity"
  ),
  T02_memory=list(
    target="T02",
    baseline=c("T01","T03"),
    module="memory_identity"
  ),
  T03_switched=list(
    target="T03",
    baseline=c("T01","T02"),
    module="class_switched_identity"
  ),
  T04_atypical=list(
    target="T04",
    baseline=c("T01","T02","T03"),
    module="atypical_activation"
  ),
  T05_IFN=list(
    target="T05",
    baseline="T01",
    module="interferon"
  ),
  T06_early_activation=list(
    target="T06",
    baseline="T01",
    module="early_activation"
  ),
  T07_secretory=list(
    target="T07",
    baseline="T01",
    module="secretory"
  ),
  T08_activation_stress=list(
    target="T08",
    baseline="T01",
    module="activation_stress"
  ),
  T11_IFN_switched=list(
    target="T11",
    baseline="T03",
    module="interferon"
  ),
  T10_platelet=list(
    target="T10",
    baseline=c("T01","T02","T03"),
    module="platelet"
  ),
  T12_Tcell=list(
    target="T12",
    baseline=c("T01","T02","T03"),
    module="T_cell"
  )
)

# ============================================================
# Metadata
# ============================================================

ann <- readRDS(ann_file)
lib <- read_tsv(lib_file)

stopifnot(
  nrow(ann) == 45955L,
  all(c(
    "project_id",
    "library_key",
    "original_cell",
    "b_plasma_final_taxonomy_id_v1"
  ) %in% names(ann))
)

lib$key <- paste(
  lib$project_id,
  lib$library_key,
  sep="|||"
)

ann$key <- paste(
  ann$project_id,
  ann$library_key,
  sep="|||"
)

# ============================================================
# Aggregate targeted native counts by library x taxonomy
# ============================================================

agg <- list()

for(i in seq_len(nrow(lib))){

  key <- lib$key[[i]]

  md <- ann[
    ann$key == key,
    ,
    drop=FALSE
  ]

  if(!nrow(md))
    next

  obj <- readRDS(
    as.character(
      lib$final_rds_resolved[[i]]
    )
  )

  counts <- get_rna_counts(obj)

  missing <- setdiff(
    genes,
    rownames(counts)
  )

  if(length(missing)){
    stop(
      "Missing canonical genes in ",
      key,
      ": ",
      paste(missing, collapse=",")
    )
  }

  cells <- as.character(
    md$original_cell
  )

  stopifnot(
    all(cells %in% colnames(counts))
  )

  x <- counts[
    genes,
    cells,
    drop=FALSE
  ]

  total_umi <- Matrix::colSums(
    counts[
      ,
      cells,
      drop=FALSE
    ]
  )

  labels <- as.character(
    md$b_plasma_final_taxonomy_id_v1
  )

  for(k in unique(labels)){

    jj <- which(labels == k)

    gene_sum <- Matrix::rowSums(
      x[,jj,drop=FALSE]
    )

    gene_detect <- Matrix::rowSums(
      x[,jj,drop=FALSE] > 0
    )

    z <- data.frame(
      project_id=as.character(md$project_id[[1]]),
      library_key=as.character(md$library_key[[1]]),
      taxonomy=k,
      n_cells=length(jj),
      total_umi=sum(total_umi[jj]),
      stringsAsFactors=FALSE
    )

    for(g in genes){
      z[[paste0("count__",g)]] <-
        as.numeric(gene_sum[[g]])

      z[[paste0("detect__",g)]] <-
        as.numeric(gene_detect[[g]])
    }

    agg[[length(agg)+1L]] <- z
  }

  cat(
    sprintf(
      "AGG %3d/%3d  %s  cells=%d\n",
      i,
      nrow(lib),
      key,
      nrow(md)
    )
  )

  rm(obj, counts, x)
  gc(verbose=FALSE)
}

agg <- do.call(rbind, agg)

# ============================================================
# Contrast helpers
# ============================================================

logcpm <- function(count, total){
  log2(
    ((count + 0.5) / (total + 1)) *
      1e6
  )
}

lib_keys <- unique(
  paste(
    agg$project_id,
    agg$library_key,
    sep="|||"
  )
)

module_rows <- list()
gene_rows <- list()

for(cname in names(contrasts)){

  cc <- contrasts[[cname]]

  module_name <- cc$module
  module_genes <- modules[[module_name]]

  for(lkey in lib_keys){

    d <- agg[
      paste(
        agg$project_id,
        agg$library_key,
        sep="|||"
      ) == lkey,
      ,
      drop=FALSE
    ]

    dt <- d[
      d$taxonomy == cc$target,
      ,
      drop=FALSE
    ]

    db <- d[
      d$taxonomy %in% cc$baseline,
      ,
      drop=FALSE
    ]

    if(
      !nrow(dt) ||
      !nrow(db)
    ){
      next
    }

    nt <- sum(dt$n_cells)
    nb <- sum(db$n_cells)

    if(
      nt < MIN_CELLS ||
      nb < MIN_CELLS
    ){
      next
    }

    ut <- sum(dt$total_umi)
    ub <- sum(db$total_umi)

    deltas <- numeric(
      length(module_genes)
    )

    names(deltas) <- module_genes

    for(g in module_genes){

      ct <- sum(
        dt[[paste0("count__",g)]]
      )

      cb <- sum(
        db[[paste0("count__",g)]]
      )

      det_t <- sum(
        dt[[paste0("detect__",g)]]
      ) / nt

      det_b <- sum(
        db[[paste0("detect__",g)]]
      ) / nb

      lt <- logcpm(ct, ut)
      lb <- logcpm(cb, ub)

      delta <- lt - lb
      deltas[[g]] <- delta

      gene_rows[[
        length(gene_rows)+1L
      ]] <- data.frame(
        contrast=cname,
        module=module_name,
        gene=g,
        project_id=dt$project_id[[1]],
        library_key=dt$library_key[[1]],
        n_target=nt,
        n_baseline=nb,
        target_logCPM=lt,
        baseline_logCPM=lb,
        delta_logCPM=delta,
        target_detection=det_t,
        baseline_detection=det_b,
        delta_detection=det_t-det_b,
        stringsAsFactors=FALSE
      )
    }

    module_rows[[
      length(module_rows)+1L
    ]] <- data.frame(
      contrast=cname,
      target=cc$target,
      baseline=paste(
        cc$baseline,
        collapse="+"
      ),
      module=module_name,
      project_id=dt$project_id[[1]],
      library_key=dt$library_key[[1]],
      n_target=nt,
      n_baseline=nb,
      module_delta_logCPM=
        mean(deltas),
      median_gene_delta_logCPM=
        median(deltas),
      positive_gene_fraction=
        mean(deltas > 0),
      stringsAsFactors=FALSE
    )
  }
}

module_lib <- do.call(
  rbind,
  module_rows
)

gene_lib <- do.call(
  rbind,
  gene_rows
)

# ============================================================
# Project-balanced summaries
# ============================================================

project_rows <- list()

for(cname in unique(module_lib$contrast)){

  d <- module_lib[
    module_lib$contrast == cname,
    ,
    drop=FALSE
  ]

  for(p in unique(d$project_id)){

    z <- d[
      d$project_id == p,
      ,
      drop=FALSE
    ]

    project_rows[[
      length(project_rows)+1L
    ]] <- data.frame(
      contrast=cname,
      module=z$module[[1]],
      project_id=p,
      n_libraries=nrow(z),
      n_target=sum(z$n_target),
      median_module_delta_logCPM=
        median(z$module_delta_logCPM),
      median_gene_delta_logCPM=
        median(z$median_gene_delta_logCPM),
      positive_library_fraction=
        mean(z$module_delta_logCPM > 0),
      stringsAsFactors=FALSE
    )
  }
}

module_project <- do.call(
  rbind,
  project_rows
)

summary_rows <- list()

for(cname in unique(module_project$contrast)){

  d <- module_project[
    module_project$contrast == cname,
    ,
    drop=FALSE
  ]

  summary_rows[[
    length(summary_rows)+1L
  ]] <- data.frame(
    contrast=cname,
    module=d$module[[1]],
    n_projects=nrow(d),
    n_libraries=sum(d$n_libraries),
    project_balanced_median_delta_logCPM=
      median(d$median_module_delta_logCPM),
    positive_project_fraction=
      mean(d$median_module_delta_logCPM > 0),
    min_project_delta_logCPM=
      min(d$median_module_delta_logCPM),
    max_project_delta_logCPM=
      max(d$median_module_delta_logCPM),
    stringsAsFactors=FALSE
  )
}

summary_df <- do.call(
  rbind,
  summary_rows
)

# ============================================================
# Gene-level project-balanced summary
# ============================================================

gene_project <- aggregate(
  cbind(
    delta_logCPM,
    delta_detection
  ) ~ contrast + module + gene + project_id,
  data=gene_lib,
  FUN=median
)

gene_summary <- do.call(
  rbind,
  lapply(
    split(
      gene_project,
      interaction(
        gene_project$contrast,
        gene_project$gene,
        drop=TRUE
      )
    ),
    function(d){

      data.frame(
        contrast=d$contrast[[1]],
        module=d$module[[1]],
        gene=d$gene[[1]],
        n_projects=nrow(d),
        project_balanced_median_delta_logCPM=
          median(d$delta_logCPM),
        positive_project_fraction=
          mean(d$delta_logCPM > 0),
        project_balanced_median_delta_detection=
          median(d$delta_detection),
        stringsAsFactors=FALSE
      )
    }
  )
)

# ============================================================
# Full taxonomy project composition
# ============================================================

composition <- as.data.frame(
  table(
    project_id=ann$project_id,
    taxonomy=
      ann$b_plasma_final_taxonomy_id_v1
  ),
  stringsAsFactors=FALSE
)

composition <- composition[
  composition$Freq > 0,
  ,
  drop=FALSE
]

# ============================================================
# Write
# ============================================================

write.table(
  module_lib,
  file=file.path(
    out_dir,
    "b_plasma_native_module_contrasts_by_library_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  module_project,
  file=file.path(
    out_dir,
    "b_plasma_native_module_contrasts_by_project_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  summary_df,
  file=file.path(
    out_dir,
    "b_plasma_native_module_contrast_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  gene_summary,
  file=file.path(
    out_dir,
    "b_plasma_native_module_gene_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  composition,
  file=file.path(
    out_dir,
    "b_plasma_full_taxonomy_project_composition_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Print
# ============================================================

cat(
  "\n===== PROJECT-BALANCED NATIVE MODULE AUDIT =====\n"
)

print(
  summary_df,
  row.names=FALSE
)

cat(
  "\n===== PROJECT-LEVEL CONTRASTS =====\n"
)

print(
  module_project,
  row.names=FALSE
)

cat(
  "\n===== GENE-LEVEL SUMMARY =====\n"
)

for(cname in unique(gene_summary$contrast)){

  cat(
    "\n--- ",
    cname,
    " ---\n",
    sep=""
  )

  d <- gene_summary[
    gene_summary$contrast == cname,
    ,
    drop=FALSE
  ]

  d <- d[
    order(
      d$project_balanced_median_delta_logCPM,
      decreasing=TRUE
    ),
    ,
    drop=FALSE
  ]

  print(
    d,
    row.names=FALSE
  )
}

writeLines(
  c(
    "PASS",
    "B/plasma full native module audit v1",
    "native SoupX-corrected RNA counts",
    "aggregation=library x frozen taxonomy",
    "contrast=within-library matched taxonomy",
    "project summary=median across eligible libraries",
    "atlas summary=project-balanced median",
    paste0(
      "minimum_cells_per_target_or_baseline=",
      MIN_CELLS
    ),
    "",
    "T09 unresolved_state intentionally has no positive biological module claim"
  ),
  done_file
)

cat(
  "\nPASS: B/plasma full native module audit completed\n"
)

