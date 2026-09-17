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
  "b_plasma_state_project_baseline_audit_v1__20260814_1530"
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

done_file <- file.path(
  out_dir,
  "B_PLASMA_STATE_PROJECT_BASELINE_AUDIT_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop("Already completed: ", done_file)
}

# ============================================================
# Modules
# ============================================================

modules <- list(

  atypical_activation=c(
    "ITGAX","FCRL5","TBX21","CD86"
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

  secretory=c(
    "XBP1","PRDM1","JCHAIN","MZB1",
    "FKBP11","DERL3","SEC11C",
    "ELL2","IRF4","TNFRSF17","SDC1"
  ),

  activation_stress=c(
    "ATF3","PPP1R15A",
    "FOS","FOSB","JUN","JUNB",
    "DUSP1","DUSP2","NFKBIA"
  )
)

contrasts <- list(

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
  )
)

genes <- unique(
  unlist(
    modules,
    use.names=FALSE
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

ann$key <- paste(
  ann$project_id,
  ann$library_key,
  sep="|||"
)

lib$key <- paste(
  lib$project_id,
  lib$library_key,
  sep="|||"
)

# ============================================================
# library x taxonomy pseudobulk
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

  if(!all(genes %in% rownames(counts))){
    stop(
      "Missing audit genes: ",
      key
    )
  }

  cells <- as.character(
    md$original_cell
  )

  stopifnot(
    all(cells %in% colnames(counts))
  )

  taxonomy <- as.character(
    md$b_plasma_final_taxonomy_id_v1
  )

  # Exclude unresolved transfer from positive taxonomy audit.
  keep <- taxonomy != "TRANSFER_UNRESOLVED"

  cells <- cells[keep]
  taxonomy <- taxonomy[keep]

  if(!length(cells))
    next

  total_umi <- Matrix::colSums(
    counts[
      ,
      cells,
      drop=FALSE
    ]
  )

  x <- counts[
    genes,
    cells,
    drop=FALSE
  ]

  for(k in unique(taxonomy)){

    jj <- which(
      taxonomy == k
    )

    gs <- Matrix::rowSums(
      x[,jj,drop=FALSE]
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
        as.numeric(gs[[g]])
    }

    agg[[length(agg)+1L]] <- z
  }

  cat(
    sprintf(
      "AGG %3d/%3d  %s\n",
      i,
      nrow(lib),
      key
    )
  )

  rm(
    obj,
    counts,
    x
  )

  gc(verbose=FALSE)
}

agg <- do.call(
  rbind,
  agg
)

# ============================================================
# Helper
# ============================================================

logcpm <- function(count, total){
  log2(
    ((count + 0.5) /
       (total + 1)) *
      1e6
  )
}

combine_baseline <- function(d, taxa){

  x <- d[
    d$taxonomy %in% taxa,
    ,
    drop=FALSE
  ]

  if(!nrow(x))
    return(NULL)

  groups <- split(
    seq_len(nrow(x)),
    paste(
      x$project_id,
      x$library_key,
      sep="|||"
    )
  )

  out <- lapply(
    groups,
    function(ii){

      y <- x[ii,,drop=FALSE]

      z <- data.frame(
        project_id=y$project_id[[1]],
        library_key=y$library_key[[1]],
        taxonomy="BASELINE",
        n_cells=sum(y$n_cells),
        total_umi=sum(y$total_umi),
        stringsAsFactors=FALSE
      )

      for(g in genes){
        z[[paste0("count__",g)]] <-
          sum(
            y[[paste0("count__",g)]]
          )
      }

      z
    }
  )

  do.call(rbind, out)
}

module_score <- function(d, module_genes){

  m <- vapply(
    module_genes,
    function(g){
      logcpm(
        d[[paste0("count__",g)]],
        d$total_umi
      )
    },
    numeric(nrow(d))
  )

  if(is.null(dim(m))){
    m <- matrix(
      m,
      ncol=1
    )
  }

  rowMeans(m)
}

# ============================================================
# Project-matched comparisons
# ============================================================

summary_rows <- list()
gene_rows <- list()

for(cname in names(contrasts)){

  cc <- contrasts[[cname]]
  module_genes <- modules[[cc$module]]

  target <- agg[
    agg$taxonomy == cc$target &
      agg$n_cells >= MIN_CELLS,
    ,
    drop=FALSE
  ]

  baseline <- combine_baseline(
    agg,
    cc$baseline
  )

  baseline <- baseline[
    baseline$n_cells >= MIN_CELLS,
    ,
    drop=FALSE
  ]

  projects <- intersect(
    unique(target$project_id),
    unique(baseline$project_id)
  )

  for(p in sort(projects)){

    tt <- target[
      target$project_id == p,
      ,
      drop=FALSE
    ]

    bb <- baseline[
      baseline$project_id == p,
      ,
      drop=FALSE
    ]

    if(
      !nrow(tt) ||
      !nrow(bb)
    ){
      next
    }

    target_scores <- module_score(
      tt,
      module_genes
    )

    baseline_scores <- module_score(
      bb,
      module_genes
    )

    delta <-
      median(target_scores) -
      median(baseline_scores)

    # Pooled pseudobulk secondary estimate.
    pooled_gene_delta <- numeric(
      length(module_genes)
    )

    names(pooled_gene_delta) <-
      module_genes

    for(g in module_genes){

      target_library_values <- logcpm(
        tt[[paste0("count__",g)]],
        tt$total_umi
      )

      baseline_library_values <- logcpm(
        bb[[paste0("count__",g)]],
        bb$total_umi
      )

      library_delta <-
        median(target_library_values) -
        median(baseline_library_values)

      pooled_delta <-
        logcpm(
          sum(tt[[paste0("count__",g)]]),
          sum(tt$total_umi)
        ) -
        logcpm(
          sum(bb[[paste0("count__",g)]]),
          sum(bb$total_umi)
        )

      pooled_gene_delta[[g]] <-
        pooled_delta

      gene_rows[[
        length(gene_rows)+1L
      ]] <- data.frame(
        contrast=cname,
        module=cc$module,
        gene=g,
        project_id=p,
        n_target_libraries=nrow(tt),
        n_baseline_libraries=nrow(bb),
        n_target_cells=sum(tt$n_cells),
        n_baseline_cells=sum(bb$n_cells),
        median_library_delta_logCPM=
          library_delta,
        pooled_delta_logCPM=
          pooled_delta,
        stringsAsFactors=FALSE
      )
    }

    summary_rows[[
      length(summary_rows)+1L
    ]] <- data.frame(
      contrast=cname,
      target=cc$target,
      baseline=paste(
        cc$baseline,
        collapse="+"
      ),
      module=cc$module,
      project_id=p,

      n_target_libraries=nrow(tt),
      n_baseline_libraries=nrow(bb),

      n_target_cells=sum(tt$n_cells),
      n_baseline_cells=sum(bb$n_cells),

      target_median_library_module_score=
        median(target_scores),

      baseline_median_library_module_score=
        median(baseline_scores),

      delta_median_library_module_score=
        delta,

      positive_target_library_fraction=
        mean(
          target_scores >
            median(baseline_scores)
        ),

      pooled_module_delta_logCPM=
        mean(pooled_gene_delta),

      positive_pooled_gene_fraction=
        mean(
          pooled_gene_delta > 0
        ),

      stringsAsFactors=FALSE
    )
  }
}

summary_df <- do.call(
  rbind,
  summary_rows
)

gene_df <- do.call(
  rbind,
  gene_rows
)

# ============================================================
# Hard sanity
# ============================================================

stopifnot(
  nrow(summary_df) > 0,
  nrow(gene_df) > 0
)

# ============================================================
# Write
# ============================================================

write.table(
  summary_df,
  file=file.path(
    out_dir,
    "b_plasma_state_project_baseline_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  gene_df,
  file=file.path(
    out_dir,
    "b_plasma_state_project_baseline_gene_detail_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Print
# ============================================================

cat(
  "\n===== PROJECT-MATCHED STATE AUDIT =====\n"
)

print(
  summary_df,
  row.names=FALSE
)

cat(
  "\n===== GENE DETAILS =====\n"
)

for(cname in unique(gene_df$contrast)){

  cat(
    "\n--- ",
    cname,
    " ---\n",
    sep=""
  )

  d <- gene_df[
    gene_df$contrast == cname,
    ,
    drop=FALSE
  ]

  d <- d[
    order(
      d$project_id,
      -d$pooled_delta_logCPM
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
    "B/plasma state project-baseline audit v1",
    "native SoupX-corrected RNA",
    "unit=library x taxonomy pseudobulk",
    "comparison=target libraries versus baseline libraries within the same project",
    paste0(
      "minimum_cells_per_library_group=",
      MIN_CELLS
    ),
    "",
    "designed as secondary audit for project-concentrated state labels",
    "does not replace the stricter within-library audit"
  ),
  done_file
)

cat(
  "\nPASS: project-baseline state audit completed\n"
)

