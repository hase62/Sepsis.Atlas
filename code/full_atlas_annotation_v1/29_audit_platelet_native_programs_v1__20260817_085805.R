#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3L){
  stop("Usage: script <root> <tag> <marker_dir>")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]
marker_dir <- normalizePath(args[[3]], mustWork=TRUE)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "platelet_native_program_audit_v1__",
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
  "PLATELET_NATIVE_PROGRAM_AUDIT_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop("Already completed: ", done_file)
}

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

project_file <- file.path(
  marker_dir,
  "platelet_r0p4_project_native_marker_stats_v1.tsv"
)

support_file <- file.path(
  marker_dir,
  "platelet_r0p4_native_marker_support_v1.tsv"
)

consensus_file <- file.path(
  marker_dir,
  "platelet_r0p4_project_balanced_marker_stats_v1.tsv.gz"
)

for(f in c(
  project_file,
  support_file,
  consensus_file
)){
  if(!file.exists(f)){
    stop("Missing input: ", f)
  }
}

# ============================================================
# Diagnostic programs
#
# These are evidence axes, not automatic biological labels.
# ============================================================

programs <- list(

  platelet_identity=c(
    "PPBP","PF4","PF4V1",
    "GP1BA","GP1BB","GP9",
    "ITGA2B","ITGB3",
    "TUBB1","GP6","CLEC1B",
    "P2RY12","TREML1","RGS18",
    "FERMT3","TLN1"
  ),

  megakaryocytic_transcription=c(
    "NFE2","GATA1","GATA2",
    "FLI1","LYL1","HEMGN",
    "ESAM","RASGRP2",
    "DAB2","MPP1",
    "RAB27B","PTGIR",
    "P2RY1"
  ),

  activation_degranulation=c(
    "SELP","THBS1",
    "CD40LG","CD36",
    "GP6","TREML1",
    "F13A1","MMRN1",
    "PF4V1","ITGB3",
    "FERMT3","TLN1"
  ),

  cytoskeletal_contractile=c(
    "MYL9","MYL6",
    "TPM1","TPM4",
    "FLNA","ACTN1",
    "PFN1","GSN",
    "ZYX","PARVB",
    "CFL1","TUBA1B",
    "TUBA1C"
  ),

  interferon=c(
    "IFI6","IFI27",
    "IFI44","IFI44L",
    "IFIT1","IFIT2","IFIT3",
    "IFITM1","IFITM2","IFITM3",
    "ISG15","MX1","MX2",
    "OAS1","OAS2","OAS3","OASL",
    "RSAD2","HERC5",
    "XAF1","LY6E",
    "STAT1","IRF7"
  ),

  myeloid=c(
    "S100A8","S100A9",
    "S100A12",
    "LYZ","FCN1",
    "VCAN","LST1",
    "CTSS","CTSD",
    "TYROBP",
    "HLA-DRA","CD74",
    "PTPRC","ITGB2"
  ),

  T_NK=c(
    "CD3D","CD3E","CD3G",
    "TRAC","TRBC1","TRBC2",
    "LCK","CD247",
    "NKG7","GNLY",
    "CCL5","GIMAP4",
    "PTPRC","ITGB2"
  ),

  erythroid=c(
    "HBA1","HBA2",
    "HBB","HBD",
    "AHSP","GYPA",
    "GYPB","ALAS2",
    "CA1","SLC4A1"
  ),

  basophil=c(
    "HDC","CLC",
    "IL3RA","FCER1A",
    "MS4A2","ENPP3",
    "CPA3","GATA2",
    "HRH4"
  )
)

# ============================================================
# Load project-level native differential statistics
# ============================================================

cat("Reading project native marker statistics...\n")

ps <- read_tsv(project_file)

required <- c(
  "project_id",
  "cluster",
  "gene",
  "n_libraries",
  "median_delta_logCPM",
  "median_delta_pct",
  "positive_library_fraction"
)

stopifnot(
  all(required %in% names(ps))
)

ps$project_id <- as.character(ps$project_id)
ps$cluster <- as.character(ps$cluster)
ps$gene <- as.character(ps$gene)

support <- read_tsv(support_file)

support$cluster <- as.character(
  support$cluster
)

clusters <- sort(
  unique(
    support$cluster
  )
)

stopifnot(
  length(clusters) == 10L
)

all_program_genes <- unique(
  unlist(
    programs,
    use.names=FALSE
  )
)

# Only retain genes required for this audit.
p <- ps[
  ps$gene %in%
    all_program_genes,
  ,
  drop=FALSE
]

rm(ps)
invisible(gc())

# ============================================================
# Project x cluster x program scores
# ============================================================

project_program <- list()

for(k in clusters){

  pk <- p[
    p$cluster == k,
    ,
    drop=FALSE
  ]

  projects <- sort(
    unique(
      pk$project_id
    )
  )

  for(prj in projects){

    pp <- pk[
      pk$project_id == prj,
      ,
      drop=FALSE
    ]

    for(program_name in names(programs)){

      genes <- programs[[program_name]]

      d <- pp[
        pp$gene %in% genes,
        ,
        drop=FALSE
      ]

      if(!nrow(d)){
        next
      }

      project_program[[
        length(project_program)+1L
      ]] <- data.frame(
        cluster=k,
        project_id=prj,
        program=program_name,

        n_program_genes_total=
          length(
            unique(genes)
          ),

        n_program_genes_observed=
          length(
            unique(
              d$gene
            )
          ),

        median_gene_delta_logCPM=
          median(
            d$median_delta_logCPM,
            na.rm=TRUE
          ),

        mean_gene_delta_logCPM=
          mean(
            d$median_delta_logCPM,
            na.rm=TRUE
          ),

        positive_gene_fraction=
          mean(
            d$median_delta_logCPM > 0,
            na.rm=TRUE
          ),

        strong_gene_fraction=
          mean(
            d$median_delta_logCPM >= 0.5,
            na.rm=TRUE
          ),

        median_gene_delta_pct=
          median(
            d$median_delta_pct,
            na.rm=TRUE
          ),

        stringsAsFactors=FALSE
      )
    }
  }
}

project_program <- do.call(
  rbind,
  project_program
)

# ============================================================
# Equal-project-weight cluster summary
# ============================================================

cluster_program <- list()

for(k in clusters){

  for(program_name in names(programs)){

    d <- project_program[
      project_program$cluster == k &
        project_program$program == program_name,
      ,
      drop=FALSE
    ]

    if(!nrow(d)){
      next
    }

    cluster_program[[
      length(cluster_program)+1L
    ]] <- data.frame(
      cluster=k,
      program=program_name,

      n_projects=nrow(d),

      projects=
        paste(
          sort(
            unique(
              d$project_id
            )
          ),
          collapse=","
        ),

      median_project_score=
        median(
          d$median_gene_delta_logCPM,
          na.rm=TRUE
        ),

      median_project_positive_gene_fraction=
        median(
          d$positive_gene_fraction,
          na.rm=TRUE
        ),

      positive_project_fraction=
        mean(
          d$median_gene_delta_logCPM > 0,
          na.rm=TRUE
        ),

      project_fraction_score_ge_0p25=
        mean(
          d$median_gene_delta_logCPM >= 0.25,
          na.rm=TRUE
        ),

      project_fraction_score_ge_0p5=
        mean(
          d$median_gene_delta_logCPM >= 0.5,
          na.rm=TRUE
        ),

      median_project_delta_pct=
        median(
          d$median_gene_delta_pct,
          na.rm=TRUE
        ),

      stringsAsFactors=FALSE
    )
  }
}

cluster_program <- do.call(
  rbind,
  cluster_program
)

# ============================================================
# Program genes among formal consensus marker candidates
# ============================================================

cons <- read_tsv(
  consensus_file
)

cons$cluster <- as.character(
  cons$cluster
)

cons$gene <- as.character(
  cons$gene
)

cons$marker_candidate <-
  as.logical(
    cons$marker_candidate
  )

marker_program <- list()

for(k in clusters){

  d <- cons[
    cons$cluster == k &
      cons$marker_candidate,
    ,
    drop=FALSE
  ]

  for(program_name in names(programs)){

    hits <- d[
      d$gene %in%
        programs[[program_name]],
      ,
      drop=FALSE
    ]

    hits <- hits[
      order(
        hits$evidence_score,
        decreasing=TRUE
      ),
      ,
      drop=FALSE
    ]

    marker_program[[
      length(marker_program)+1L
    ]] <- data.frame(
      cluster=k,
      program=program_name,
      n_consensus_marker_hits=nrow(hits),
      consensus_marker_hits=
        paste(
          hits$gene,
          collapse=","
        ),
      stringsAsFactors=FALSE
    )
  }
}

marker_program <- do.call(
  rbind,
  marker_program
)

# ============================================================
# Merge and produce compact wide diagnostic table
# ============================================================

summary_long <- merge(
  cluster_program,
  marker_program,
  by=c(
    "cluster",
    "program"
  ),
  all.x=TRUE,
  sort=FALSE
)

summary_long <- summary_long[
  order(
    as.numeric(
      summary_long$cluster
    ),
    summary_long$program
  ),
  ,
  drop=FALSE
]

compact <- support[
  ,
  c(
    "cluster",
    "n_discovery_cells",
    "n_projects_present",
    "n_libraries_present",
    "n_eligible_libraries",
    "n_projects_with_eligible_libraries"
  ),
  drop=FALSE
]

for(program_name in names(programs)){

  d <- summary_long[
    summary_long$program ==
      program_name,
    ,
    drop=FALSE
  ]

  mi <- match(
    compact$cluster,
    d$cluster
  )

  compact[[
    paste0(
      program_name,
      "_score"
    )
  ]] <- d$median_project_score[
    mi
  ]

  compact[[
    paste0(
      program_name,
      "_positive_projects"
    )
  ]] <- d$positive_project_fraction[
    mi
  ]

  compact[[
    paste0(
      program_name,
      "_marker_hits"
    )
  ]] <- d$consensus_marker_hits[
    mi
  ]
}

# ============================================================
# Conservative diagnostic screen
#
# Only obvious non-platelet and IFN cases are automatically
# screened. Megakaryocytic/activation axes are intentionally
# NOT converted to final labels here.
# ============================================================

compact$diagnostic_scope_screen <-
  "Platelet_megakaryocyte"

compact$diagnostic_state_screen <-
  "none"

compact$diagnostic_reason <- ""

myeloid_hit <-
  compact$myeloid_score >= 0.5 &
  compact$myeloid_positive_projects >= 0.75

tnk_hit <-
  compact$T_NK_score >= 0.5 &
  compact$T_NK_positive_projects >= 0.75

ifn_hit <-
  compact$interferon_score >= 0.25 &
  compact$interferon_positive_projects >= 0.75

compact$diagnostic_scope_screen[
  myeloid_hit
] <- "deferred_non_platelet_myeloid_like"

compact$diagnostic_reason[
  myeloid_hit
] <- "project-balanced native myeloid program"

compact$diagnostic_scope_screen[
  tnk_hit
] <- "deferred_non_platelet_T_NK_like"

compact$diagnostic_reason[
  tnk_hit
] <- "project-balanced native T/NK program"

platelet_rows <-
  compact$diagnostic_scope_screen ==
    "Platelet_megakaryocyte"

compact$diagnostic_state_screen[
  platelet_rows &
    ifn_hit
] <- "interferon_stimulated_candidate"

compact$diagnostic_reason[
  platelet_rows &
    ifn_hit
] <- "project-balanced native interferon program"

# ============================================================
# Outputs
# ============================================================

write.table(
  project_program,
  file=file.path(
    out_dir,
    "platelet_r0p4_project_program_scores_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  summary_long,
  file=file.path(
    out_dir,
    "platelet_r0p4_project_balanced_program_summary_long_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  compact,
  file=file.path(
    out_dir,
    "platelet_r0p4_program_audit_compact_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

writeLines(
  c(
    "PASS",
    "Platelet native program audit v1",
    "source=step28 project-native marker statistics",
    "backbone_resolution=r0.4",
    "",
    "projects receive equal weight",
    "no raw RNA reread",
    "no integrated-assay marker evidence",
    "no cellranger_count access",
    "",
    "megakaryocytic transcription is diagnostic only",
    "megakaryocytic transcript enrichment does not by itself imply intact megakaryocyte identity",
    "activation/degranulation is diagnostic only",
    "candidate screen is NOT final annotation"
  ),
  done_file
)

cat(
  "\n===== COMPACT PROGRAM AUDIT =====\n"
)

print(
  compact,
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Platelet native program audit completed\n"
)
