#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

suppressPackageStartupMessages({
  library(SeuratObject)
  library(Matrix)
})

# ============================================================
# Paths
# ============================================================

in_rds <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_native_subclustering_v1_primary9997",
  "b_plasma_pilot_native_subclustering_v1.rds"
)

marker_file <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_project_balanced_marker_audit_r0p4_v1",
  "b_plasma_r0p4_reproducible_top_markers_v1.tsv"
)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_underlying_identity_assignment_v1"
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

done_file <- file.path(
  out_dir,
  "B_PLASMA_IDENTITY_ASSIGNMENT_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop(
    "Output already completed; refusing to overwrite: ",
    done_file
  )
}

for(f in c(in_rds, marker_file)){
  if(!file.exists(f)){
    stop("Missing input: ", f)
  }
}

# ============================================================
# Load
# ============================================================

bp <- readRDS(in_rds)

stopifnot(
  nrow(bp) == 38606,
  ncol(bp) == 9997
)

md <- bp[[]]

stopifnot(
  all(c(
    "project_id",
    "cluster_res_0p4"
  ) %in% names(md))
)

cl <- as.character(
  md$cluster_res_0p4
)

project <- as.character(
  md$project_id
)

counts <- tryCatch(
  LayerData(
    bp[["RNA"]],
    layer="counts"
  ),
  error=function(e) NULL
)

if(is.null(counts)){
  counts <- GetAssayData(
    bp,
    assay="RNA",
    layer="counts"
  )
}

stopifnot(
  nrow(counts) == 38606,
  ncol(counts) == 9997,
  identical(
    colnames(counts),
    rownames(md)
  )
)

all_genes <- rownames(counts)

cat(
  "cells    =", ncol(bp), "\n",
  "features =", nrow(bp), "\n"
)

# ============================================================
# Frozen high-confidence anchor identities
#
# 1 = Naive B-like
# 2 = Memory B-like
# 3 = Class-switched memory B-like
#
# Clusters 4 and 7 are deliberately NOT used as references.
# They are evaluated as independent positive-control targets.
# ============================================================

anchor_map <- c(
  "1"="Naive_B_like",
  "2"="Memory_B_like",
  "3"="Class_switched_memory_B_like"
)

target_clusters <- c(
  "4","5","6","7","8",
  "9","10","11","12","14"
)

non_b_clusters <- c(
  "13","15"
)

stopifnot(
  all(
    c(
      names(anchor_map),
      target_clusters,
      non_b_clusters
    ) %in% unique(cl)
  )
)

# ============================================================
# Identity feature set
#
# Start from cross-project reproducible markers of anchor
# clusters 1/2/3.
# ============================================================

mk <- read.delim(
  marker_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

mk$cluster <- as.character(
  mk$cluster
)

mk <- mk[
  mk$cluster %in%
    names(anchor_map),
  ,
  drop=FALSE
]

# Top 80 reproducible markers per anchor.
marker_genes <- unlist(
  lapply(
    names(anchor_map),
    function(k){

      d <- mk[
        mk$cluster == k,
        ,
        drop=FALSE
      ]

      if("technical_gene" %in% names(d)){
        d <- d[
          !d$technical_gene,
          ,
          drop=FALSE
        ]
      }

      if("rearranged_ig_gene" %in% names(d)){
        d <- d[
          !d$rearranged_ig_gene,
          ,
          drop=FALSE
        ]
      }

      head(
        as.character(d$gene),
        80
      )
    }
  ),
  use.names=FALSE
)

# Ensure canonical identity axes are represented.
canonical_identity <- c(

  # naive
  "TCL1A",
  "IGHD",
  "IGHM",
  "FCER2",
  "IL4R",
  "BACH2",
  "FCRL1",
  "VPREB3",

  # memory
  "CD27",
  "TNFRSF13B",
  "AIM2",
  "GPR183",
  "CD82",
  "ITGB1",
  "POU2AF1",

  # class-switched memory
  "IGHA1",
  "IGHA2",
  "IGHG1",
  "IGHG2",
  "IGHG3",
  "IGHG4"
)

candidate_genes <- unique(
  c(
    marker_genes,
    canonical_identity
  )
)

candidate_genes <- intersect(
  candidate_genes,
  all_genes
)

# ============================================================
# Remove orthogonal state / technical axes
#
# Constant immunoglobulin isotypes are RETAINED.
# Rearranged V/J immunoglobulin genes are removed.
# ============================================================

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
  candidate_genes
)

rearranged_ig <- grepl(
  paste0(
    "^IGHV|",
    "^IGKV|",
    "^IGLV|",
    "^IGHJ[0-9]|",
    "^IGKJ[0-9]|",
    "^IGLJ[0-9]"
  ),
  candidate_genes
)

interferon_genes <- c(
  "IFI6","IFI27","IFI35","IFI44","IFI44L",
  "IFIT1","IFIT2","IFIT3","IFIT5",
  "IFITM1","IFITM2","IFITM3",
  "ISG15","MX1","MX2",
  "OAS1","OAS2","OAS3","OASL",
  "XAF1","IRF7","STAT1","STAT2",
  "GBP1","GBP2","GBP4","GBP5",
  "BST2","LY6E","PLSCR1","UBE2L6",
  "TRIM22","EIF2AK2","EPSTI1"
)

activation_stress_genes <- c(
  "FOS","FOSB",
  "JUN","JUNB","JUND",
  "EGR1","EGR2","EGR3",
  "ATF3",
  "NFKBIA",
  "DUSP1","DUSP2","DUSP5",
  "NR4A1","NR4A2","NR4A3",
  "CD69",
  "PPP1R15A",
  "HSPA1A","HSPA1B","HSPA8",
  "DNAJB1",
  "KLF4"
)

cycling_genes <- c(
  "MKI67","TOP2A","STMN1",
  "TYMS","UBE2C","CDC20",
  "CDK1","CCNB1","CCNB2",
  "TUBA1B"
)

secretory_genes <- c(
  "MZB1","JCHAIN",
  "XBP1","PRDM1",
  "TNFRSF17","SDC1",
  "FKBP11","DERL3",
  "SEC11C","ELL2",
  "IRF4"
)

contamination_genes <- c(

  # platelet
  "PPBP","PF4","PF4V1",
  "GP1BB","GP9","ITGA2B",
  "TUBB1","TREML1","MPIG6B",

  # T/NK
  "CD3D","CD3E","CD3G",
  "TRAC","TRBC1","TRBC2",
  "CD247","CD2",
  "NKG7","GNLY","CCL5",

  # myeloid
  "LYZ","S100A8","S100A9",
  "FCN1","VCAN",
  "LILRB1","CTSS"
)

state_or_contam <- candidate_genes %in%
  unique(c(
    interferon_genes,
    activation_stress_genes,
    cycling_genes,
    secretory_genes,
    contamination_genes
  ))

identity_genes <- candidate_genes[
  !technical &
  !rearranged_ig &
  !state_or_contam
]

identity_genes <- unique(
  identity_genes
)

if(length(identity_genes) < 50){
  stop(
    "Too few identity genes after filtering: ",
    length(identity_genes)
  )
}

cat(
  "identity genes =", length(identity_genes), "\n"
)

# ============================================================
# Feature provenance
# ============================================================

feature_table <- data.frame(
  gene=candidate_genes,
  source_cluster1=
    candidate_genes %in%
      mk$gene[mk$cluster == "1"],
  source_cluster2=
    candidate_genes %in%
      mk$gene[mk$cluster == "2"],
  source_cluster3=
    candidate_genes %in%
      mk$gene[mk$cluster == "3"],
  canonical_identity=
    candidate_genes %in%
      canonical_identity,
  technical=technical,
  rearranged_ig=rearranged_ig,
  interferon=
    candidate_genes %in%
      interferon_genes,
  activation_stress=
    candidate_genes %in%
      activation_stress_genes,
  cycling=
    candidate_genes %in%
      cycling_genes,
  secretory=
    candidate_genes %in%
      secretory_genes,
  contamination=
    candidate_genes %in%
      contamination_genes,
  retained_identity_feature=
    candidate_genes %in%
      identity_genes,
  stringsAsFactors=FALSE
)

write.table(
  feature_table,
  file=file.path(
    out_dir,
    "b_plasma_identity_feature_set_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# Helpers
# ============================================================

libsize <- Matrix::colSums(
  counts
)

names(libsize) <- colnames(counts)

logcpm_profile <- function(idx){

  if(!length(idx)){
    stop("Empty cell index")
  }

  gene_counts <- Matrix::rowSums(
    counts[
      identity_genes,
      idx,
      drop=FALSE
    ]
  )

  total <- sum(
    libsize[idx]
  )

  if(total <= 0){
    stop("Zero pseudobulk library size")
  }

  log2(
    gene_counts /
      total *
      1e6 +
      1
  )
}

safe_cor <- function(x,y,method){

  ok <- is.finite(x) &
    is.finite(y)

  if(sum(ok) < 20){
    return(NA_real_)
  }

  suppressWarnings(
    cor(
      x[ok],
      y[ok],
      method=method
    )
  )
}

# ============================================================
# Project baselines
#
# Exclude known non-B clusters 13 and 15.
# ============================================================

projects <- sort(
  unique(project)
)

project_baseline <- list()
project_baseline_n <- integer()

for(p in projects){

  ii <- which(
    project == p &
    !cl %in% non_b_clusters
  )

  if(length(ii) < 50){
    next
  }

  project_baseline[[p]] <-
    logcpm_profile(ii)

  project_baseline_n[p] <-
    length(ii)
}

cat(
  "projects with baseline =",
  length(project_baseline),
  "\n"
)

# ============================================================
# Anchor project-specific profiles
#
# profile = anchor pseudobulk - project B-cell baseline
# ============================================================

MIN_ANCHOR_CELLS <- 30L

anchor_profiles <- list()
anchor_meta <- list()

for(k in names(anchor_map)){

  identity <- unname(
    anchor_map[k]
  )

  for(p in projects){

    if(is.null(
      project_baseline[[p]]
    )){
      next
    }

    ii <- which(
      project == p &
      cl == k
    )

    if(length(ii) <
       MIN_ANCHOR_CELLS){
      next
    }

    prof <- logcpm_profile(ii) -
      project_baseline[[p]]

    key <- paste(
      identity,
      p,
      sep="|||"
    )

    anchor_profiles[[key]] <- prof

    anchor_meta[[
      length(anchor_meta)+1L
    ]] <- data.frame(
      cluster=k,
      identity=identity,
      project_id=p,
      n_cells=length(ii),
      stringsAsFactors=FALSE
    )
  }
}

anchor_meta <- do.call(
  rbind,
  anchor_meta
)

cat(
  "\n===== ANCHOR PROJECT COUNTS =====\n"
)

print(
  table(anchor_meta$identity)
)

# ============================================================
# Build project-balanced reference signatures
#
# Equal project weight via gene-wise median.
# exclude_project is used to avoid residual study matching.
# ============================================================

build_reference <- function(
  identity,
  exclude_project=NULL
){

  keys <- names(
    anchor_profiles
  )

  prefix <- paste0(
    identity,
    "|||"
  )

  keys <- keys[
    startsWith(
      keys,
      prefix
    )
  ]

  if(!is.null(
    exclude_project
  )){

    suffix <- paste0(
      "|||",
      exclude_project
    )

    keys <- keys[
      !endsWith(
        keys,
        suffix
      )
    ]
  }

  # Need at least two independent projects.
  if(length(keys) < 2){
    return(NULL)
  }

  m <- do.call(
    cbind,
    anchor_profiles[keys]
  )

  apply(
    m,
    1,
    median,
    na.rm=TRUE
  )
}

reference_identities <- unname(
  anchor_map
)

# ============================================================
# Leave-one-project-out validation on anchors
# ============================================================

anchor_validation <- list()

for(i in seq_len(
  nrow(anchor_meta)
)){

  true_identity <-
    anchor_meta$identity[i]

  p <-
    anchor_meta$project_id[i]

  key <- paste(
    true_identity,
    p,
    sep="|||"
  )

  target <- anchor_profiles[[key]]

  sims <- lapply(
    reference_identities,
    function(ref_identity){

      ref <- build_reference(
        ref_identity,
        exclude_project=p
      )

      if(is.null(ref)){
        return(
          data.frame(
            identity=ref_identity,
            pearson=NA_real_,
            spearman=NA_real_,
            combined=NA_real_
          )
        )
      }

      pe <- safe_cor(
        target,
        ref,
        "pearson"
      )

      sp <- safe_cor(
        target,
        ref,
        "spearman"
      )

      data.frame(
        identity=ref_identity,
        pearson=pe,
        spearman=sp,
        combined=mean(
          c(pe,sp),
          na.rm=TRUE
        )
      )
    }
  )

  sims <- do.call(
    rbind,
    sims
  )

  valid <- is.finite(
    sims$combined
  )

  if(!any(valid)){
    next
  }

  winner <- sims$identity[
    which.max(
      sims$combined
    )
  ]

  for(j in seq_len(
    nrow(sims)
  )){

    anchor_validation[[
      length(anchor_validation)+1L
    ]] <- data.frame(
      true_identity=true_identity,
      project_id=p,
      n_cells=
        anchor_meta$n_cells[i],
      compared_identity=
        sims$identity[j],
      pearson=
        sims$pearson[j],
      spearman=
        sims$spearman[j],
      combined=
        sims$combined[j],
      winner=winner,
      correct=
        winner ==
        true_identity,
      stringsAsFactors=FALSE
    )
  }
}

anchor_validation <- do.call(
  rbind,
  anchor_validation
)

write.table(
  anchor_validation,
  file=file.path(
    out_dir,
    "b_plasma_anchor_loocv_validation_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

anchor_calls <- unique(
  anchor_validation[
    ,
    c(
      "true_identity",
      "project_id",
      "n_cells",
      "winner",
      "correct"
    )
  ]
)

validation_accuracy <- mean(
  anchor_calls$correct
)

cat(
  "\n===== ANCHOR LOPO VALIDATION =====\n"
)

print(
  anchor_calls,
  row.names=FALSE
)

cat(
  "\naccuracy = ",
  validation_accuracy,
  "\n",
  sep=""
)

# ============================================================
# Target project profiles and similarity
# ============================================================

MIN_TARGET_CELLS <- 10L

target_similarity <- list()

for(k in target_clusters){

  for(p in projects){

    if(is.null(
      project_baseline[[p]]
    )){
      next
    }

    ii <- which(
      project == p &
      cl == k
    )

    if(length(ii) <
       MIN_TARGET_CELLS){
      next
    }

    target <- logcpm_profile(ii) -
      project_baseline[[p]]

    sim_rows <- list()

    for(ref_identity in
        reference_identities){

      # Critical:
      # exclude the target's own project from the
      # reference signature.
      ref <- build_reference(
        ref_identity,
        exclude_project=p
      )

      if(is.null(ref)){
        next
      }

      pe <- safe_cor(
        target,
        ref,
        "pearson"
      )

      sp <- safe_cor(
        target,
        ref,
        "spearman"
      )

      comb <- mean(
        c(pe,sp),
        na.rm=TRUE
      )

      sim_rows[[
        length(sim_rows)+1L
      ]] <- data.frame(
        cluster=k,
        project_id=p,
        n_cells=length(ii),
        compared_identity=
          ref_identity,
        pearson=pe,
        spearman=sp,
        combined=comb,
        stringsAsFactors=FALSE
      )
    }

    if(!length(sim_rows)){
      next
    }

    sim_rows <- do.call(
      rbind,
      sim_rows
    )

    valid <- is.finite(
      sim_rows$combined
    )

    if(!any(valid)){
      next
    }

    ord <- order(
      sim_rows$combined,
      decreasing=TRUE,
      na.last=NA
    )

    winner <-
      sim_rows$compared_identity[
        ord[1]
      ]

    best <-
      sim_rows$combined[
        ord[1]
      ]

    second <- if(
      length(ord) >= 2
    ){
      sim_rows$combined[
        ord[2]
      ]
    } else {
      NA_real_
    }

    margin <- best - second

    sim_rows$winner <-
      winner

    sim_rows$winner_combined <-
      best

    sim_rows$winner_margin <-
      margin

    target_similarity[[
      length(target_similarity)+1L
    ]] <- sim_rows
  }
}

target_similarity <- do.call(
  rbind,
  target_similarity
)

write.table(
  target_similarity,
  file=file.path(
    out_dir,
    "b_plasma_target_project_similarity_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# Aggregate project-level evidence equally across projects
# ============================================================

summary_list <- list()

for(k in target_clusters){

  d <- target_similarity[
    target_similarity$cluster == k,
    ,
    drop=FALSE
  ]

  if(!nrow(d)){

    summary_list[[
      length(summary_list)+1L
    ]] <- data.frame(
      cluster=k,
      n_cells=sum(cl == k),
      n_project_profiles=0,
      candidate_identity=NA_character_,
      median_best_similarity=NA_real_,
      median_margin=NA_real_,
      project_winner_fraction=NA_real_,
      assignment_strength="unresolved",
      stringsAsFactors=FALSE
    )

    next
  }

  # One winner per project.
  proj_calls <- unique(
    d[
      ,
      c(
        "project_id",
        "winner",
        "winner_combined",
        "winner_margin"
      )
    ]
  )

  # Median similarity for each candidate identity.
  agg <- aggregate(
    combined ~ compared_identity,
    data=d,
    FUN=median,
    na.rm=TRUE
  )

  agg <- agg[
    order(
      agg$combined,
      decreasing=TRUE
    ),
    ,
    drop=FALSE
  ]

  candidate <-
    agg$compared_identity[1]

  best_similarity <-
    agg$combined[1]

  second_similarity <-
    if(nrow(agg) >= 2){
      agg$combined[2]
    } else {
      NA_real_
    }

  median_margin <-
    best_similarity -
    second_similarity

  winner_fraction <-
    mean(
      proj_calls$winner ==
        candidate
    )

  n_profiles <-
    nrow(proj_calls)

  # Diagnostic evidence strength only.
  # This is NOT a final biological freeze decision.
  strength <- if(
    n_profiles >= 2 &&
    is.finite(best_similarity) &&
    best_similarity >= 0.15 &&
    is.finite(median_margin) &&
    median_margin >= 0.08 &&
    winner_fraction >= 2/3
  ){
    "strong"
  } else if(
    is.finite(best_similarity) &&
    best_similarity >= 0.10 &&
    is.finite(median_margin) &&
    median_margin >= 0.04 &&
    winner_fraction >= 0.5
  ){
    "moderate"
  } else {
    "weak"
  }

  summary_list[[
    length(summary_list)+1L
  ]] <- data.frame(
    cluster=k,
    n_cells=sum(cl == k),
    n_project_profiles=
      n_profiles,
    candidate_identity=
      candidate,
    median_best_similarity=
      best_similarity,
    median_margin=
      median_margin,
    project_winner_fraction=
      winner_fraction,
    assignment_strength=
      strength,
    stringsAsFactors=FALSE
  )
}

summary_df <- do.call(
  rbind,
  summary_list
)

summary_df$cluster_num <-
  as.integer(
    summary_df$cluster
  )

summary_df <- summary_df[
  order(
    summary_df$cluster_num
  ),
  ,
  drop=FALSE
]

summary_df$cluster_num <- NULL

# ============================================================
# Add fixed anchors and known non-B rows for complete r0.4 map
# ============================================================

anchor_rows <- data.frame(
  cluster=names(anchor_map),
  n_cells=vapply(
    names(anchor_map),
    function(k)
      sum(cl == k),
    numeric(1)
  ),
  n_project_profiles=NA_integer_,
  candidate_identity=
    unname(anchor_map),
  median_best_similarity=NA_real_,
  median_margin=NA_real_,
  project_winner_fraction=NA_real_,
  assignment_strength=
    "anchor_high_confidence",
  stringsAsFactors=FALSE
)

nonb_rows <- data.frame(
  cluster=c("13","15"),
  n_cells=c(
    sum(cl == "13"),
    sum(cl == "15")
  ),
  n_project_profiles=NA_integer_,
  candidate_identity=c(
    "Deferred_non_B_platelet_like",
    "Deferred_non_B_T_cell_like"
  ),
  median_best_similarity=NA_real_,
  median_margin=NA_real_,
  project_winner_fraction=NA_real_,
  assignment_strength=
    "non_B_by_marker",
  stringsAsFactors=FALSE
)

complete <- rbind(
  anchor_rows,
  summary_df,
  nonb_rows
)

complete$cluster_num <-
  as.integer(
    complete$cluster
  )

complete <- complete[
  order(
    complete$cluster_num
  ),
  ,
  drop=FALSE
]

complete$cluster_num <- NULL

stopifnot(
  nrow(complete) == 15,
  sum(complete$n_cells) == 9997
)

write.table(
  summary_df,
  file=file.path(
    out_dir,
    "b_plasma_target_identity_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

write.table(
  complete,
  file=file.path(
    out_dir,
    "b_plasma_r0p4_complete_identity_diagnostic_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# Print
# ============================================================

cat(
  "\n===== TARGET IDENTITY SUMMARY =====\n"
)

print(
  summary_df,
  row.names=FALSE
)

cat(
  "\n===== COMPLETE r0.4 DIAGNOSTIC MAP =====\n"
)

print(
  complete,
  row.names=FALSE
)

writeLines(
  c(
    "PASS",
    "B/plasma state-excluded underlying identity assignment v1",
    "n_cells=9997",
    paste0(
      "n_identity_genes=",
      length(identity_genes)
    ),
    paste0(
      "anchor_LOPO_accuracy=",
      sprintf(
        "%.6f",
        validation_accuracy
      )
    ),
    "anchors=cluster1 Naive_B_like; cluster2 Memory_B_like; cluster3 Class_switched_memory_B_like",
    "reference signatures=project-centered pseudobulk median across projects",
    "target project excluded from its reference signatures",
    "interferon/activation/cycling/secretory/contamination genes excluded",
    "assignments are diagnostic and not yet frozen"
  ),
  done_file
)

cat(
  "\nPASS: B/plasma underlying identity diagnostic completed\n"
)

