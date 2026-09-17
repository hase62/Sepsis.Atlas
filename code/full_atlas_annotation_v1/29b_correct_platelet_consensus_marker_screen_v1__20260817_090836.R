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
    "platelet_consensus_marker_screen_v1__",
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
  "PLATELET_CONSENSUS_MARKER_SCREEN_COMPLETE.ok"
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

to_logical <- function(x){

  if(is.logical(x)){
    return(x)
  }

  if(is.numeric(x)){
    return(x != 0)
  }

  tolower(as.character(x)) %in%
    c("true", "t", "1", "yes")
}

cons_file <- file.path(
  marker_dir,
  "platelet_r0p4_project_balanced_marker_stats_v1.tsv.gz"
)

support_file <- file.path(
  marker_dir,
  "platelet_r0p4_native_marker_support_v1.tsv"
)

stopifnot(
  file.exists(cons_file),
  file.exists(support_file)
)

cons <- read_tsv(cons_file)
support <- read_tsv(support_file)

cons$cluster <- as.character(cons$cluster)
cons$gene <- as.character(cons$gene)
cons$marker_candidate <- to_logical(cons$marker_candidate)

support$cluster <- as.character(support$cluster)

markers <- cons[
  cons$marker_candidate,
  ,
  drop=FALSE
]

# ============================================================
# Specific marker sets.
#
# Important:
# - PTPRC and ITGB2 are deliberately NOT T/NK-defining.
# - GATA2 is deliberately NOT basophil-defining.
# - IFITM genes alone are deliberately insufficient for IFN.
# ============================================================

sets <- list(

  platelet_identity=c(
    "PPBP","PF4","PF4V1",
    "GP1BA","GP1BB","GP9",
    "ITGA2B","ITGB3",
    "TUBB1","GP6","CLEC1B",
    "P2RY12","TREML1",
    "RGS18","FERMT3"
  ),

  megakaryocytic_transcription=c(
    "NFE2","GATA1","GATA2",
    "FLI1","LYL1","HEMGN",
    "ESAM","RASGRP2",
    "DAB2","MPP1",
    "RAB27B","PTGIR",
    "P2RY1"
  ),

  activation_degranulation_specific=c(
    "SELP",
    "THBS1",
    "CD36",
    "F13A1",
    "MMRN1",
    "CD40LG"
  ),

  interferon_specific=c(
    "IFI6","IFI27",
    "IFI44","IFI44L",
    "IFIT1","IFIT2","IFIT3",
    "ISG15",
    "MX1","MX2",
    "OAS1","OAS2","OAS3","OASL",
    "RSAD2","HERC5",
    "XAF1","LY6E",
    "STAT1","IRF7"
  ),

  myeloid_specific=c(
    "S100A8","S100A9","S100A12",
    "LYZ","FCN1","VCAN",
    "LST1","CTSS","CTSD",
    "TYROBP",
    "HLA-DRA","CD74",
    "FCER1G"
  ),

  T_NK_specific=c(
    "CD3D","CD3E","CD3G",
    "TRAC","TRBC1","TRBC2",
    "LCK","CD247",
    "NKG7","GNLY",
    "GIMAP4","CCL5"
  ),

  erythroid_specific=c(
    "HBA1","HBA2","HBB","HBD",
    "AHSP","GYPA","GYPB",
    "ALAS2","SLC4A1","CA1"
  ),

  basophil_specific=c(
    "HDC","CLC",
    "IL3RA","FCER1A",
    "MS4A2","ENPP3",
    "CPA3","HRH4"
  )
)

clusters <- sort(
  unique(support$cluster)
)

out <- support

get_hits <- function(k, genes){

  d <- markers[
    markers$cluster == k &
      markers$gene %in% genes,
    ,
    drop=FALSE
  ]

  if(!nrow(d)){
    return(character())
  }

  d <- d[
    order(
      d$evidence_score,
      decreasing=TRUE
    ),
    ,
    drop=FALSE
  ]

  unique(d$gene)
}

for(nm in names(sets)){

  out[[paste0("n_", nm, "_hits")]] <- 0L
  out[[paste0(nm, "_hits")]] <- ""
}

for(i in seq_len(nrow(out))){

  k <- out$cluster[[i]]

  for(nm in names(sets)){

    hits <- get_hits(
      k,
      sets[[nm]]
    )

    out[[paste0(
      "n_",
      nm,
      "_hits"
    )]][i] <- length(hits)

    out[[paste0(
      nm,
      "_hits"
    )]][i] <- paste(
      hits,
      collapse=","
    )
  }
}

# ============================================================
# Corrected conservative scope screen
# ============================================================

out$corrected_scope_screen <-
  "Platelet_like"

out$corrected_lineage_evidence <-
  "retained_platelet_compartment_without_specific_contaminant_program"

myeloid <- out$n_myeloid_specific_hits >= 4L
tnk <- out$n_T_NK_specific_hits >= 4L
ery <- out$n_erythroid_specific_hits >= 4L
baso <- out$n_basophil_specific_hits >= 4L

# Ambiguous multiple contamination programs are not silently
# overwritten. They become unresolved.
n_contam <- (
  as.integer(myeloid) +
  as.integer(tnk) +
  as.integer(ery) +
  as.integer(baso)
)

out$corrected_scope_screen[
  n_contam > 1L
] <- "Deferred_non_platelet_lineage_unresolved"

out$corrected_lineage_evidence[
  n_contam > 1L
] <- "multiple lineage-specific consensus programs"

out$corrected_scope_screen[
  myeloid &
    n_contam == 1L
] <- "Deferred_non_platelet_myeloid_like"

out$corrected_lineage_evidence[
  myeloid &
    n_contam == 1L
] <- ">=4 specific myeloid consensus markers"

out$corrected_scope_screen[
  tnk &
    n_contam == 1L
] <- "Deferred_non_platelet_T_NK_like"

out$corrected_lineage_evidence[
  tnk &
    n_contam == 1L
] <- ">=4 specific T/NK consensus markers"

out$corrected_scope_screen[
  ery &
    n_contam == 1L
] <- "Deferred_non_platelet_erythroid_like"

out$corrected_lineage_evidence[
  ery &
    n_contam == 1L
] <- ">=4 specific erythroid consensus markers"

out$corrected_scope_screen[
  baso &
    n_contam == 1L
] <- "Deferred_non_platelet_basophil_like"

out$corrected_lineage_evidence[
  baso &
    n_contam == 1L
] <- ">=4 specific basophil consensus markers"

# ============================================================
# Separate platelet transcriptional axis
# ============================================================

retained <- out$corrected_scope_screen ==
  "Platelet_like"

out$megakaryocytic_transcription_screen <-
  "none"

out$megakaryocytic_transcription_confidence <-
  "not_applicable"

strong_mk <-
  retained &
  out$n_megakaryocytic_transcription_hits >= 6L

candidate_mk <-
  retained &
  out$n_megakaryocytic_transcription_hits >= 3L &
  !strong_mk

out$megakaryocytic_transcription_screen[
  strong_mk
] <- "megakaryocytic_transcription_high"

out$megakaryocytic_transcription_confidence[
  strong_mk
] <- "high"

out$megakaryocytic_transcription_screen[
  candidate_mk
] <- "megakaryocytic_transcription_candidate"

out$megakaryocytic_transcription_confidence[
  candidate_mk
] <- "medium"

# ============================================================
# Orthogonal state
#
# IFN takes priority over activation only when there are >=4
# specific consensus IFN genes.
#
# IFITM-only enrichment is NOT sufficient.
# ============================================================

out$corrected_state_screen <- "none"
out$corrected_state_confidence <- "not_applicable"

ifn <-
  retained &
  out$n_interferon_specific_hits >= 4L

activation <-
  retained &
  out$n_activation_degranulation_specific_hits >= 4L

out$corrected_state_screen[
  activation
] <- "activation_degranulation_candidate"

out$corrected_state_confidence[
  activation
] <- "medium"

out$corrected_state_screen[
  ifn
] <- "interferon_stimulated_candidate"

out$corrected_state_confidence[
  ifn
] <- "medium"

# ============================================================
# Compact decision-oriented output
# ============================================================

keep <- c(
  "cluster",
  "n_discovery_cells",
  "n_projects_present",
  "n_libraries_present",
  "n_eligible_libraries",
  "n_projects_with_eligible_libraries",

  "n_platelet_identity_hits",
  "platelet_identity_hits",

  "n_megakaryocytic_transcription_hits",
  "megakaryocytic_transcription_hits",

  "n_activation_degranulation_specific_hits",
  "activation_degranulation_specific_hits",

  "n_interferon_specific_hits",
  "interferon_specific_hits",

  "n_myeloid_specific_hits",
  "myeloid_specific_hits",

  "n_T_NK_specific_hits",
  "T_NK_specific_hits",

  "n_erythroid_specific_hits",
  "erythroid_specific_hits",

  "n_basophil_specific_hits",
  "basophil_specific_hits",

  "corrected_scope_screen",
  "corrected_lineage_evidence",

  "megakaryocytic_transcription_screen",
  "megakaryocytic_transcription_confidence",

  "corrected_state_screen",
  "corrected_state_confidence"
)

compact <- out[
  ,
  keep,
  drop=FALSE
]

write.table(
  compact,
  file=file.path(
    out_dir,
    "platelet_r0p4_corrected_consensus_marker_screen_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

writeLines(
  c(
    "PASS",
    "Platelet corrected consensus-marker screen v1",
    "source=step28 formal project-balanced consensus marker candidates",
    "backbone_resolution=r0.4",
    "",
    "relative module-score-only classification is NOT used",
    "absolute expression is not inferred from differential scores",
    "PTPRC/ITGB2 excluded from T/NK-specific definition",
    "GATA2 excluded from basophil-specific definition",
    "IFITM-only enrichment insufficient for IFN state",
    "",
    "scope, megakaryocytic transcription, and orthogonal state kept separate",
    "no RNA reread",
    "no integrated assay marker evidence",
    "no cellranger_count access",
    "",
    "screen is decision evidence, not yet final 49,585-cell freeze"
  ),
  done_file
)

cat(
  "\n===== CORRECTED CONSENSUS MARKER SCREEN =====\n"
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

cat("\nPASS\n")
