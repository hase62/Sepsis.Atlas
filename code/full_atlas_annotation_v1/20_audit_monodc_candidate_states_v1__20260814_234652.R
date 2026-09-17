#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 2L){
  stop("Usage: script <root> <tag>")
}

root <- normalizePath(
  args[[1]],
  mustWork=TRUE
)

tag <- args[[2]]

in_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "monodc_project_balanced_marker_audit_v1"
)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "monodc_candidate_state_audit_v1__",
    tag
  )
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

summary_file <- file.path(
  in_dir,
  "monodc_r0p6_project_balanced_summary_v1.tsv"
)

if(!file.exists(summary_file)){
  stop(
    "Missing summary: ",
    summary_file
  )
}

# ============================================================
# Load r0.6 cluster summary
# ============================================================

cs <- read.delim(
  summary_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

cs$cluster <- as.character(
  cs$cluster
)

targets <- c(
  "3",
  "14",
  "21",
  "30"
)

stopifnot(
  all(targets %in% cs$cluster)
)

# ============================================================
# Auto-detect project-balanced gene-level marker TSV
#
# Deliberately search ONLY the existing Monocyte/DC audit dir.
# ============================================================

required_cols <- c(
  "cluster",
  "gene",
  "n_eligible_projects",
  "median_log2FC",
  "median_delta_pct",
  "positive_project_fraction",
  "strong_project_fraction"
)

files <- list.files(
  in_dir,
  pattern="\\.tsv$",
  full.names=TRUE
)

candidate_files <- character()

for(f in files){

  x <- tryCatch(
    read.delim(
      f,
      sep="\t",
      quote="\"",
      comment.char="",
      stringsAsFactors=FALSE,
      check.names=FALSE,
      nrows=3
    ),
    error=function(e) NULL
  )

  if(
    !is.null(x) &&
    all(required_cols %in% names(x))
  ){
    candidate_files <- c(
      candidate_files,
      f
    )
  }
}

if(!length(candidate_files)){
  stop(
    "Could not auto-detect gene-level ",
    "project-balanced marker TSV in ",
    in_dir
  )
}

# Prefer an r0p6-labelled file if multiple tables qualify.
ord <- order(
  !grepl(
    "r0p6",
    basename(candidate_files),
    ignore.case=TRUE
  ),
  nchar(
    basename(candidate_files)
  )
)

marker_file <- candidate_files[
  ord[1]
]

cat(
  "marker_file =",
  marker_file,
  "\n"
)

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
  mk$cluster %in% targets,
  ,
  drop=FALSE
]

stopifnot(
  all(targets %in% mk$cluster)
)

for(v in c(
  "n_eligible_projects",
  "median_log2FC",
  "median_delta_pct",
  "positive_project_fraction",
  "strong_project_fraction"
)){
  mk[[v]] <- suppressWarnings(
    as.numeric(
      mk[[v]]
    )
  )
}

# ============================================================
# Evidence weight
#
# This is an audit score, NOT an automated annotation rule.
# Effect size × project reproducibility.
# ============================================================

mk$gene_weight <-
  pmax(
    mk$median_delta_pct,
    0
  ) *
  pmax(
    mk$positive_project_fraction,
    0
  ) *
  (
    0.5 +
    0.5 *
      pmax(
        mk$strong_project_fraction,
        0
      )
  )

# ============================================================
# Prespecified biological modules
# ============================================================

modules <- list(

  DC2_identity=c(
    "FCER1A","CD1C","CD1E",
    "CLEC10A","NAPSB",
    "CD74",
    "HLA-DRA","HLA-DRB1",
    "HLA-DPA1","HLA-DPB1",
    "HLA-DQA1","HLA-DQB1"
  ),

  interferon=c(
    "IFI6","IFI27","IFI35",
    "IFI44","IFI44L",
    "IFIT1","IFIT2","IFIT3","IFIT5",
    "IFITM1","IFITM2","IFITM3",
    "ISG15",
    "MX1","MX2",
    "OAS1","OAS2","OAS3","OASL",
    "XAF1","IRF7","STAT1","STAT2",
    "BST2","EPSTI1","SIGLEC1",
    "PLSCR1","UBE2L6","EIF2AK2"
  ),

  MS1_emergency_myeloid=c(
    "RETN","ALOX5AP","IL1R2",
    "S100A8","S100A9","S100A12",
    "VCAN","PLAC8",
    "FCN1","MNDA",
    "LILRB1","LILRB3",
    "CTSD","CTSB",
    "SELL","CTSS"
  ),

  inflammatory_NFkB=c(
    "NFKBIA","NFKBIZ",
    "TNF","IL1B","CXCL8",
    "CCL3","CCL4",
    "PTGS2",
    "DUSP1","DUSP2",
    "EGR1",
    "FOS","FOSB",
    "JUN","JUNB",
    "PELI1","PELI2",
    "PDE4B",
    "TLR2","TLR4",
    "AREG"
  ),

  nonclassical_monocyte=c(
    "FCGR3A","MS4A7",
    "CDKN1C","LST1",
    "IFITM3","LILRB1",
    "SERPINA1","SAT1",
    "FCER1G","TYROBP"
  ),

  platelet_contamination=c(
    "PPBP","PF4","PF4V1",
    "GP1BB","GP9",
    "ITGA2B","TUBB1",
    "TREML1","MPIG6B"
  ),

  cycling=c(
    "MKI67","TOP2A","STMN1",
    "PCLAF","PCNA",
    "MCM2","MCM3","MCM4",
    "MCM5","MCM6","MCM7",
    "TYMS","RRM2"
  )
)

# ============================================================
# Module audit
# ============================================================

res <- list()

for(k in targets){

  d <- mk[
    mk$cluster == k,
    ,
    drop=FALSE
  ]

  # Keep one record per gene if necessary.
  d <- d[
    order(
      d$gene,
      -d$gene_weight
    ),
    ,
    drop=FALSE
  ]

  d <- d[
    !duplicated(
      d$gene
    ),
    ,
    drop=FALSE
  ]

  for(m in names(modules)){

    mod <- unique(
      modules[[m]]
    )

    h <- d[
      d$gene %in% mod,
      ,
      drop=FALSE
    ]

    if(nrow(h)){
      h <- h[
        order(
          h$gene_weight,
          decreasing=TRUE
        ),
        ,
        drop=FALSE
      ]
    }

    res[[
      length(res)+1L
    ]] <- data.frame(
      cluster=k,
      module=m,
      module_size=length(mod),
      n_hits=nrow(h),
      hit_fraction=
        nrow(h) / length(mod),
      score_sum=
        if(nrow(h)){
          sum(
            h$gene_weight,
            na.rm=TRUE
          )
        } else {
          0
        },
      score_mean_hit=
        if(nrow(h)){
          mean(
            h$gene_weight,
            na.rm=TRUE
          )
        } else {
          0
        },
      n_strong_hits=
        if(nrow(h)){
          sum(
            h$strong_project_fraction >= 0.5,
            na.rm=TRUE
          )
        } else {
          0
        },
      hit_genes=
        if(nrow(h)){
          paste(
            h$gene,
            collapse=","
          )
        } else {
          ""
        },
      stringsAsFactors=FALSE
    )
  }
}

module_scores <- do.call(
  rbind,
  res
)

# ============================================================
# Top weighted native markers for audit trail
# ============================================================

top_list <- list()

for(k in targets){

  d <- mk[
    mk$cluster == k,
    ,
    drop=FALSE
  ]

  d <- d[
    order(
      d$gene_weight,
      decreasing=TRUE
    ),
    ,
    drop=FALSE
  ]

  d <- head(
    d,
    40
  )

  top_list[[
    length(top_list)+1L
  ]] <- d
}

top_markers <- do.call(
  rbind,
  top_list
)

# ============================================================
# Provisional biological interpretation
#
# These remain conservative:
# core identity and orthogonal state are kept separate.
# ============================================================

decision <- data.frame(

  cluster=c(
    "3",
    "14",
    "21",
    "30"
  ),

  proposed_core=c(
    "Classical_monocyte_like",
    "DC2_like",
    "Classical_monocyte_like",
    "Classical_monocyte_like"
  ),

  proposed_state=c(
    "RETN_PADI4_myeloid_state_candidate",
    "none",
    "interferon_stimulated_candidate",
    "LPAR1_PDE4B_state_candidate"
  ),

  proposed_confidence=c(
    "low_medium",
    "high",
    "low_medium",
    "medium"
  ),

  rationale=c(
    paste0(
      "RETN/PADI4 and activation markers; ",
      "only limited cross-project support; ",
      "do not force MS1 identity"
    ),
    paste0(
      "FCER1A/CD1C/HLA-DQA1/NAPSB program; ",
      "interpret as DC2 core-identity refinement, ",
      "not an orthogonal state"
    ),
    paste0(
      "strong interferon-response program; ",
      "limited number of eligible projects, ",
      "therefore candidate state"
    ),
    paste0(
      "reproducible LPAR1/PDE4B/PELI2/NFKBIZ-axis; ",
      "retain data-driven candidate name rather than ",
      "forcing a canonical monocyte state"
    )
  ),

  stringsAsFactors=FALSE
)

# ============================================================
# Outputs
# ============================================================

write.table(
  cs[
    cs$cluster %in% targets,
    ,
    drop=FALSE
  ],
  file=file.path(
    out_dir,
    "monodc_target_cluster_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

write.table(
  module_scores,
  file=file.path(
    out_dir,
    "monodc_candidate_state_module_scores_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

write.table(
  top_markers,
  file=file.path(
    out_dir,
    "monodc_candidate_state_top_markers_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

write.table(
  decision,
  file=file.path(
    out_dir,
    "monodc_candidate_state_provisional_decision_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

writeLines(
  c(
    "PASS",
    "Monocyte/DC candidate state audit v1",
    paste0(
      "marker_file=",
      marker_file
    ),
    "targets=3,14,21,30",
    "module scores are descriptive audit evidence only",
    "core identity and orthogonal state remain separate",
    "no raw Cell Ranger output accessed"
  ),
  file.path(
    out_dir,
    "MONODC_CANDIDATE_STATE_AUDIT_COMPLETE.ok"
  )
)

cat(
  "\n===== TARGET CLUSTER SUMMARY =====\n"
)

print(
  cs[
    cs$cluster %in% targets,
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat(
  "\n===== MODULE SCORES =====\n"
)

print(
  module_scores[
    order(
      as.integer(
        module_scores$cluster
      ),
      -module_scores$score_sum
    ),
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat(
  "\n===== PROVISIONAL DECISION =====\n"
)

print(
  decision,
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Monocyte/DC candidate-state audit completed\n"
)
