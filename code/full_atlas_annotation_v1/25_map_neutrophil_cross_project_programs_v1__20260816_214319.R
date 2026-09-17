#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3L){
  stop("Usage: script <root> <tag> <marker_dir>")
}

root <- normalizePath(
  args[[1]],
  mustWork=TRUE
)

tag <- args[[2]]

marker_dir <- normalizePath(
  args[[3]],
  mustWork=TRUE
)

suppressPackageStartupMessages({
  library(Matrix)
})

MIN_CELLS_PER_LIBRARY_CLUSTER <- 10L
TOP_MARKERS_FOR_OVERLAP <- 100L

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "neutrophil_cross_project_program_mapping_v1__",
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
  "NEUTROPHIL_CROSS_PROJECT_PROGRAM_MAPPING_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop("Already completed: ", done_file)
}

# ============================================================
# Helpers
# ============================================================

read_tsv <- function(path){

  con <- if(
    grepl("\\.gz$", path)
  ){
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }

  on.exit(
    close(con)
  )

  read.delim(
    con,
    sep="\t",
    quote="\"",
    comment.char="",
    stringsAsFactors=FALSE,
    check.names=FALSE
  )
}

logcpm <- function(count, total){

  log2(
    ((count + 0.5) /
       (total + 1)) *
      1e6
  )
}

# ============================================================
# Prespecified programs
#
# Core maturation and orthogonal states are intentionally
# separated.
# ============================================================

modules <- list(

  mature_circulating=c(
    "FCGR3B","CXCR2","CSF3R",
    "FPR1","FPR2",
    "VNN2","SELL",
    "SERPINA1","CMTM2",
    "IL1R2"
  ),

  early_granulopoiesis=c(
    "MPO","ELANE","PRTN3",
    "AZU1","CTSG",
    "DEFA1","DEFA1B",
    "DEFA3","DEFA4",
    "MS4A3","CEBPE"
  ),

  late_immature_granulopoiesis=c(
    "LTF","LCN2","MMP8",
    "CAMP","CEACAM8",
    "BPI","PGLYRP1",
    "RETN","ARG1",
    "OLR1","CD24","TCN1"
  ),

  interferon=c(
    "IFI6","IFI27","IFI35",
    "IFI44","IFI44L",
    "IFIT1","IFIT2","IFIT3","IFIT5",
    "IFITM1","IFITM2","IFITM3",
    "ISG15",
    "MX1","MX2",
    "OAS1","OAS2","OAS3","OASL",
    "RSAD2","HERC5",
    "XAF1","EPSTI1",
    "STAT1","STAT2",
    "IRF7","IRF9",
    "UBE2L6","EIF2AK2"
  ),

  inflammatory_NFkB=c(
    "TNFAIP3","TNFAIP6",
    "NFKBIA","NFKBIZ",
    "IRAK2","IL1RN",
    "PTGS2","ICAM1",
    "CD83","IER3",
    "ZFP36","BCL2A1",
    "PDE4B"
  ),

  immediate_early_stress=c(
    "FOS","FOSB",
    "JUN","JUNB",
    "DUSP1","DUSP2",
    "EGR1","EGR2",
    "ATF3",
    "PPP1R15A",
    "GADD45B"
  ),

  cycling=c(
    "MKI67","TOP2A",
    "STMN1","PCLAF",
    "PCNA","TYMS","RRM2",
    "MCM2","MCM3","MCM4",
    "MCM5","MCM6","MCM7",
    "CENPF","NUSAP1",
    "KIF11","DTL"
  ),

  platelet=c(
    "PPBP","PF4","PF4V1",
    "GP1BB","GP9",
    "ITGA2B","TUBB1",
    "TREML1","MPIG6B",
    "CAVIN2","GNG11",
    "NRGN"
  ),

  basophil=c(
    "HDC","GATA2",
    "IL3RA","FCER1A",
    "MS4A2","ENPP3",
    "CLC","GCSAML",
    "HRH4","PTGER3",
    "CD200R1"
  ),

  T_NK=c(
    "CD3D","CD3E","CD3G",
    "TRAC","TRBC1","TRBC2",
    "LCK","BCL11B",
    "IL7R","TCF7",
    "NKG7","GNLY",
    "NCAM1","PRF1",
    "CCL5"
  ),

  eosinophil=c(
    "SIGLEC8","IL5RA",
    "PRG2","PRG3",
    "RNASE2","RNASE3",
    "EPX","CLC",
    "CCR3","ALOX15"
  )
)

# Marker-support sets used only for conservative screening.
marker_programs <- list(

  platelet=c(
    "PPBP","PF4","GP1BB",
    "GP9","ITGA2B","TUBB1",
    "TREML1","MPIG6B"
  ),

  basophil=c(
    "HDC","GATA2","IL3RA",
    "FCER1A","MS4A2",
    "ENPP3","CLC"
  ),

  T_NK=c(
    "CD3D","CD3E","CD3G",
    "TRAC","TRBC1","TRBC2",
    "LCK","IL7R","TCF7",
    "NKG7","GNLY","NCAM1"
  ),

  early_immature=c(
    "MPO","ELANE","PRTN3",
    "AZU1",
    "DEFA1","DEFA1B",
    "DEFA3","DEFA4",
    "MS4A3"
  ),

  late_immature=c(
    "LTF","LCN2","MMP8",
    "CAMP","CEACAM8",
    "BPI","PGLYRP1",
    "RETN","ARG1",
    "OLR1","CD24"
  ),

  mature=c(
    "FCGR3B","CXCR2","CSF3R",
    "FPR1","FPR2",
    "VNN2","SELL",
    "SERPINA1","IL1R2"
  ),

  IFN=c(
    "IFI6","IFI27","IFI44","IFI44L",
    "IFIT1","IFIT2","IFIT3",
    "ISG15","MX1","MX2",
    "RSAD2","HERC5",
    "XAF1","EPSTI1",
    "OAS1","OAS2","OAS3","OASL",
    "STAT1"
  ),

  inflammatory=c(
    "TNFAIP3","TNFAIP6",
    "NFKBIA","NFKBIZ",
    "IRAK2","IL1RN",
    "PTGS2","ICAM1",
    "CD83","IER3"
  ),

  cycling=c(
    "MKI67","TOP2A","STMN1",
    "PCLAF","PCNA",
    "TYMS","RRM2",
    "MCM2","MCM3","MCM4",
    "MCM5","MCM6","MCM7",
    "CENPF","KIF11","DTL"
  )
)

# ============================================================
# Cluster summary from step 24
# ============================================================

cluster_summary_file <- file.path(
  marker_dir,
  "neutrophil_native_cluster_summary_r0p4_v1.tsv"
)

cluster_summary <- read_tsv(
  cluster_summary_file
)

cluster_summary$project_id <-
  as.character(cluster_summary$project_id)

cluster_summary$cluster <-
  as.character(cluster_summary$cluster)

cluster_summary$key <- paste(
  cluster_summary$project_id,
  cluster_summary$cluster,
  sep="|||"
)

stopifnot(
  !anyDuplicated(cluster_summary$key)
)

projects <- sort(
  unique(cluster_summary$project_id)
)

stopifnot(
  length(projects) == 6L
)

# ============================================================
# Absolute native pseudobulk module scores
# ============================================================

library_module_rows <- list()

for(p in projects){

  agg_file <- file.path(
    marker_dir,
    paste0(
      p,
      "_native_library_cluster_pseudobulk_r0p4_v1.rds"
    )
  )

  if(!file.exists(agg_file)){
    stop(
      "Missing pseudobulk checkpoint: ",
      agg_file
    )
  }

  z <- readRDS(
    agg_file
  )

  features <- as.character(
    z$features
  )

  aggs <- z$aggregates

  module_idx <- lapply(
    modules,
    function(g){
      match(
        intersect(
          g,
          features
        ),
        features
      )
    }
  )

  for(lk in names(aggs)){

    a <- aggs[[lk]]

    for(k in a$clusters){

      n <- as.integer(
        a$n_by_cluster[[k]]
      )

      if(
        is.na(n) ||
        n <
          MIN_CELLS_PER_LIBRARY_CLUSTER
      ){
        next
      }

      total <- as.numeric(
        a$umi_by_cluster[[k]]
      )

      if(
        is.na(total) ||
        total <= 0
      ){
        next
      }

      for(m in names(modules)){

        idx <- module_idx[[m]]

        if(!length(idx)){
          next
        }

        counts <- as.numeric(
          a$count_by_cluster[
            idx,
            k
          ]
        )

        detected <- as.numeric(
          a$detect_by_cluster[
            idx,
            k
          ]
        )

        library_module_rows[[
          length(library_module_rows)+1L
        ]] <- data.frame(
          project_id=p,
          library_key=lk,
          cluster=as.character(k),
          module=m,
          n_cells=n,
          total_umi=total,
          n_module_genes=length(idx),
          module_logCPM=
            mean(
              logcpm(
                counts,
                total
              )
            ),
          module_detect_fraction=
            mean(
              detected /
                n
            ),
          stringsAsFactors=FALSE
        )
      }
    }
  }

  rm(
    z,
    aggs
  )

  invisible(
    gc()
  )
}

library_module <- do.call(
  rbind,
  library_module_rows
)

# ============================================================
# project x cluster module summary
# ============================================================

groups <- split(
  seq_len(nrow(library_module)),
  paste(
    library_module$project_id,
    library_module$cluster,
    library_module$module,
    sep="|||"
  )
)

module_summary_list <- lapply(
  groups,
  function(ii){

    d <- library_module[
      ii,
      ,
      drop=FALSE
    ]

    data.frame(
      project_id=d$project_id[[1]],
      cluster=d$cluster[[1]],
      module=d$module[[1]],
      n_libraries=nrow(d),
      n_cells=sum(d$n_cells),
      median_module_logCPM=
        median(
          d$module_logCPM
        ),
      median_module_detect_fraction=
        median(
          d$module_detect_fraction
        ),
      stringsAsFactors=FALSE
    )
  }
)

module_summary <- do.call(
  rbind,
  module_summary_list
)

module_summary$key <- paste(
  module_summary$project_id,
  module_summary$cluster,
  sep="|||"
)

# ============================================================
# Within-project relative module enrichment
# ============================================================

module_summary$delta_vs_project_cluster_median <-
  NA_real_

module_summary$project_module_rank_fraction <-
  NA_real_

module_summary$robust_z_within_project <-
  NA_real_

for(p in projects){

  for(m in names(modules)){

    ii <- which(
      module_summary$project_id == p &
        module_summary$module == m
    )

    if(!length(ii)){
      next
    }

    v <- module_summary$median_module_logCPM[
      ii
    ]

    module_summary$delta_vs_project_cluster_median[
      ii
    ] <-
      v -
      median(v)

    module_summary$project_module_rank_fraction[
      ii
    ] <-
      rank(
        v,
        ties.method="average"
      ) /
      length(v)

    if(
      length(v) >= 4L
    ){

      s <- mad(
        v,
        center=median(v),
        constant=1
      )

      if(
        is.finite(s) &&
        s > 0
      ){

        module_summary$robust_z_within_project[
          ii
        ] <-
          (
            v -
              median(v)
          ) /
          s
      }
    }
  }
}

# ============================================================
# Marker candidate sets
# ============================================================

stat_files <- list.files(
  marker_dir,
  pattern="_native_marker_stats_r0p4_v1\\.tsv\\.gz$",
  full.names=TRUE
)

marker_stats <- list()

for(f in stat_files){

  d <- read_tsv(f)

  if(nrow(d)){
    marker_stats[[
      length(marker_stats)+1L
    ]] <- d
  }
}

if(length(marker_stats)){

  marker_stats <- do.call(
    rbind,
    marker_stats
  )

} else {

  marker_stats <- data.frame()
}

marker_sets <- list()

if(nrow(marker_stats)){

  marker_stats$project_id <-
    as.character(
      marker_stats$project_id
    )

  marker_stats$cluster <-
    as.character(
      marker_stats$cluster
    )

  marker_stats$gene <-
    as.character(
      marker_stats$gene
    )

  marker_stats$key <- paste(
    marker_stats$project_id,
    marker_stats$cluster,
    sep="|||"
  )

  # Defensive normalization of logical column.
  marker_stats$marker_candidate <-
    as.logical(
      marker_stats$marker_candidate
    )

  keys <- unique(
    marker_stats$key
  )

  for(key in keys){

    d <- marker_stats[
      marker_stats$key == key &
        marker_stats$marker_candidate &
        !marker_stats$technical_gene,
      ,
      drop=FALSE
    ]

    d <- d[
      order(
        d$evidence_score,
        decreasing=TRUE
      ),
      ,
      drop=FALSE
    ]

    marker_sets[[key]] <-
      head(
        unique(d$gene),
        TOP_MARKERS_FOR_OVERLAP
      )
  }
}

# ============================================================
# Conservative candidate screen
#
# This is NOT the final annotation.
# ============================================================

support_count <- function(
  genes,
  program
){

  sum(
    unique(genes) %in%
      marker_programs[[program]]
  )
}

screen <- cluster_summary[
  ,
  c(
    "project_id",
    "cluster",
    "n_cells",
    "n_libraries",
    "max_library_fraction",
    "healthy_fraction",
    "disease_fraction",
    "n_eligible_marker_libraries",
    "marker_support_status"
  ),
  drop=FALSE
]

screen$key <- paste(
  screen$project_id,
  screen$cluster,
  sep="|||"
)

program_names <- names(
  marker_programs
)

for(pr in program_names){

  screen[[
    paste0(
      "n_marker_",
      pr
    )
  ]] <- 0L
}

screen$top_marker_set <- ""

for(i in seq_len(
  nrow(screen)
)){

  key <- screen$key[[i]]

  genes <- marker_sets[[key]]

  if(is.null(genes)){
    genes <- character()
  }

  screen$top_marker_set[[i]] <-
    paste(
      genes,
      collapse=","
    )

  for(pr in program_names){

    screen[[paste0("n_marker_", pr)]][i] <-
      support_count(
        genes,
        pr
      )
  }
}

# Discovery role:
# library-dominated clusters are retained, but they do not
# define the cross-project taxonomy.
screen$taxonomy_discovery_role <-
  ifelse(
    screen$marker_support_status == "eligible" &
      screen$max_library_fraction < 0.90,
    "taxonomy_discovery_eligible",
    "assignment_only"
  )

screen$provisional_core_screen <-
  "Neutrophil_maturation_unresolved"

screen$provisional_state_screen <-
  "none"

screen$provisional_scope_screen <-
  "Neutrophil"

screen$screen_reason <- ""

for(i in seq_len(
  nrow(screen)
)){

  # ----------------------------------------------------------
  # Non-neutrophil contamination first
  # ----------------------------------------------------------

  if(
    screen$n_marker_platelet[i] >= 4L
  ){

    screen$provisional_core_screen[i] <-
      "Deferred_non_neutrophil_platelet_like"

    screen$provisional_state_screen[i] <-
      "not_applicable"

    screen$provisional_scope_screen[i] <-
      "deferred_non_neutrophil"

    screen$screen_reason[i] <-
      ">=4 platelet marker candidates"

    next
  }

  if(
    screen$n_marker_basophil[i] >= 4L
  ){

    screen$provisional_core_screen[i] <-
      "Deferred_non_neutrophil_basophil_like"

    screen$provisional_state_screen[i] <-
      "not_applicable"

    screen$provisional_scope_screen[i] <-
      "deferred_non_neutrophil"

    screen$screen_reason[i] <-
      ">=4 basophil marker candidates"

    next
  }

  if(
    screen$n_marker_T_NK[i] >= 4L
  ){

    screen$provisional_core_screen[i] <-
      "Deferred_non_neutrophil_T_NK_like"

    screen$provisional_state_screen[i] <-
      "not_applicable"

    screen$provisional_scope_screen[i] <-
      "deferred_non_neutrophil"

    screen$screen_reason[i] <-
      ">=4 T/NK marker candidates"

    next
  }

  # ----------------------------------------------------------
  # Neutrophil maturation
  # ----------------------------------------------------------

  early <- screen$n_marker_early_immature[i]
  late  <- screen$n_marker_late_immature[i]
  mature <- screen$n_marker_mature[i]

  if(
    early >= 4L &&
    late >= 4L
  ){

    screen$provisional_core_screen[i] <-
      "Immature_neutrophil_mixed_granulopoiesis_candidate"

    screen$screen_reason[i] <-
      "early and late granulopoiesis marker support"

  } else if(
    early >= 4L
  ){

    screen$provisional_core_screen[i] <-
      "Early_immature_neutrophil_like_candidate"

    screen$screen_reason[i] <-
      "early granulopoiesis marker support"

  } else if(
    late >= 4L
  ){

    screen$provisional_core_screen[i] <-
      "Late_immature_neutrophil_like_candidate"

    screen$screen_reason[i] <-
      "LTF/LCN2/MMP8-class late immature marker support"

  } else if(
    mature >= 4L
  ){

    screen$provisional_core_screen[i] <-
      "Mature_neutrophil_like_candidate"

    screen$screen_reason[i] <-
      "mature circulating neutrophil marker support"
  }

  # ----------------------------------------------------------
  # Orthogonal state
  # ----------------------------------------------------------

  if(
    screen$n_marker_IFN[i] >= 5L
  ){

    screen$provisional_state_screen[i] <-
      "interferon_stimulated_candidate"

  } else if(
    screen$n_marker_inflammatory[i] >= 4L
  ){

    screen$provisional_state_screen[i] <-
      "inflammatory_activation_candidate"

  } else if(
    screen$n_marker_cycling[i] >= 4L
  ){

    screen$provisional_state_screen[i] <-
      "cycling_candidate"
  }
}

# ============================================================
# Recurrence of candidate programs across projects
#
# Only discovery-eligible clusters count as evidence defining
# a taxonomy family.
# ============================================================

disc <- screen[
  screen$taxonomy_discovery_role ==
    "taxonomy_discovery_eligible",
  ,
  drop=FALSE
]

family_core <- split(
  disc,
  disc$provisional_core_screen
)

core_recurrence <- do.call(
  rbind,
  lapply(
    names(family_core),
    function(x){

      d <- family_core[[x]]

      data.frame(
        dimension="core",
        candidate=x,
        n_projects=
          length(
            unique(
              d$project_id
            )
          ),
        n_clusters=nrow(d),
        n_cells=sum(d$n_cells),
        projects=
          paste(
            sort(
              unique(
                d$project_id
              )
            ),
            collapse=","
          ),
        stringsAsFactors=FALSE
      )
    }
  )
)

family_state <- split(
  disc,
  disc$provisional_state_screen
)

state_recurrence <- do.call(
  rbind,
  lapply(
    names(family_state),
    function(x){

      d <- family_state[[x]]

      data.frame(
        dimension="state",
        candidate=x,
        n_projects=
          length(
            unique(
              d$project_id
            )
          ),
        n_clusters=nrow(d),
        n_cells=sum(d$n_cells),
        projects=
          paste(
            sort(
              unique(
                d$project_id
              )
            ),
            collapse=","
          ),
        stringsAsFactors=FALSE
      )
    }
  )
)

recurrence <- rbind(
  core_recurrence,
  state_recurrence
)

# ============================================================
# Cross-project similarity
#
# Similarity combines:
#   1. module-profile concordance
#   2. top native-marker Jaccard overlap
#
# It is diagnostic only.
# ============================================================

get_module_vector <- function(key){

  d <- module_summary[
    module_summary$key == key,
    ,
    drop=FALSE
  ]

  v <- d$delta_vs_project_cluster_median

  names(v) <- d$module

  v[
    names(modules)
  ]
}

jaccard <- function(a, b){

  if(
    length(a) < 5L ||
    length(b) < 5L
  ){
    return(NA_real_)
  }

  u <- union(a, b)

  if(!length(u)){
    return(NA_real_)
  }

  length(
    intersect(a, b)
  ) /
    length(u)
}

keys <- screen$key

pair_rows <- list()

for(i in seq_len(
  length(keys) - 1L
)){

  for(j in seq.int(
    i + 1L,
    length(keys)
  )){

    p1 <- screen$project_id[
      match(
        keys[[i]],
        screen$key
      )
    ]

    p2 <- screen$project_id[
      match(
        keys[[j]],
        screen$key
      )
    ]

    if(p1 == p2){
      next
    }

    v1 <- get_module_vector(
      keys[[i]]
    )

    v2 <- get_module_vector(
      keys[[j]]
    )

    ok <- is.finite(v1) &
      is.finite(v2)

    rho <- if(
      sum(ok) >= 5L
    ){
      suppressWarnings(
        cor(
          v1[ok],
          v2[ok],
          method="spearman"
        )
      )
    } else {
      NA_real_
    }

    module_similarity <- if(
      is.finite(rho)
    ){
      (rho + 1) / 2
    } else {
      NA_real_
    }

    m1 <- marker_sets[[keys[[i]]]]

    m2 <- marker_sets[[keys[[j]]]]

    if(is.null(m1)){
      m1 <- character()
    }

    if(is.null(m2)){
      m2 <- character()
    }

    jac <- jaccard(
      m1,
      m2
    )

    combined <- if(
      is.finite(module_similarity) &&
      is.finite(jac)
    ){

      0.70 *
        module_similarity +
        0.30 *
        jac

    } else if(
      is.finite(module_similarity)
    ){

      module_similarity

    } else if(
      is.finite(jac)
    ){

      jac

    } else {

      NA_real_
    }

    pair_rows[[
      length(pair_rows)+1L
    ]] <- data.frame(
      key1=keys[[i]],
      project1=p1,
      cluster1=
        screen$cluster[
          match(
            keys[[i]],
            screen$key
          )
        ],
      key2=keys[[j]],
      project2=p2,
      cluster2=
        screen$cluster[
          match(
            keys[[j]],
            screen$key
          )
        ],
      module_spearman=rho,
      module_similarity_0to1=
        module_similarity,
      marker_jaccard=jac,
      combined_similarity=
        combined,
      stringsAsFactors=FALSE
    )
  }
}

pairs <- do.call(
  rbind,
  pair_rows
)

# ============================================================
# Nearest cross-project neighbours for every cluster
# ============================================================

neighbor_rows <- list()

for(key in keys){

  d1 <- pairs[
    pairs$key1 == key,
    ,
    drop=FALSE
  ]

  if(nrow(d1)){

    a <- data.frame(
      source_key=d1$key1,
      source_project=d1$project1,
      source_cluster=d1$cluster1,
      neighbor_key=d1$key2,
      neighbor_project=d1$project2,
      neighbor_cluster=d1$cluster2,
      module_spearman=d1$module_spearman,
      marker_jaccard=d1$marker_jaccard,
      combined_similarity=
        d1$combined_similarity,
      stringsAsFactors=FALSE
    )

    neighbor_rows[[
      length(neighbor_rows)+1L
    ]] <- a
  }

  d2 <- pairs[
    pairs$key2 == key,
    ,
    drop=FALSE
  ]

  if(nrow(d2)){

    b <- data.frame(
      source_key=d2$key2,
      source_project=d2$project2,
      source_cluster=d2$cluster2,
      neighbor_key=d2$key1,
      neighbor_project=d2$project1,
      neighbor_cluster=d2$cluster1,
      module_spearman=d2$module_spearman,
      marker_jaccard=d2$marker_jaccard,
      combined_similarity=
        d2$combined_similarity,
      stringsAsFactors=FALSE
    )

    neighbor_rows[[
      length(neighbor_rows)+1L
    ]] <- b
  }
}

neighbors <- do.call(
  rbind,
  neighbor_rows
)

# Keep best neighbour from each other project.
best_rows <- list()

groups <- split(
  seq_len(nrow(neighbors)),
  paste(
    neighbors$source_key,
    neighbors$neighbor_project,
    sep="|||"
  )
)

for(ii in groups){

  d <- neighbors[
    ii,
    ,
    drop=FALSE
  ]

  d <- d[
    order(
      d$combined_similarity,
      decreasing=TRUE,
      na.last=TRUE
    ),
    ,
    drop=FALSE
  ]

  best_rows[[
    length(best_rows)+1L
  ]] <- d[1,,drop=FALSE]
}

nearest <- do.call(
  rbind,
  best_rows
)

nearest <- nearest[
  order(
    nearest$source_project,
    as.numeric(
      nearest$source_cluster
    ),
    -nearest$combined_similarity
  ),
  ,
  drop=FALSE
]

# Add candidate labels.
source_idx <- match(
  nearest$source_key,
  screen$key
)

neighbor_idx <- match(
  nearest$neighbor_key,
  screen$key
)

nearest$source_core_screen <-
  screen$provisional_core_screen[
    source_idx
  ]

nearest$source_state_screen <-
  screen$provisional_state_screen[
    source_idx
  ]

nearest$neighbor_core_screen <-
  screen$provisional_core_screen[
    neighbor_idx
  ]

nearest$neighbor_state_screen <-
  screen$provisional_state_screen[
    neighbor_idx
  ]

# ============================================================
# Write outputs
# ============================================================

write.table(
  library_module,
  file=file.path(
    out_dir,
    "neutrophil_library_cluster_module_scores_r0p4_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  module_summary[
    ,
    setdiff(
      names(module_summary),
      "key"
    ),
    drop=FALSE
  ],
  file=file.path(
    out_dir,
    "neutrophil_cluster_program_scores_r0p4_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  screen[
    ,
    setdiff(
      names(screen),
      "key"
    ),
    drop=FALSE
  ],
  file=file.path(
    out_dir,
    "neutrophil_cluster_candidate_screen_r0p4_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  recurrence,
  file=file.path(
    out_dir,
    "neutrophil_candidate_family_recurrence_r0p4_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  pairs,
  file=file.path(
    out_dir,
    "neutrophil_cross_project_cluster_similarity_r0p4_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  nearest,
  file=file.path(
    out_dir,
    "neutrophil_cross_project_nearest_neighbors_r0p4_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Completion
# ============================================================

writeLines(
  c(
    "PASS",
    "Neutrophil cross-project program mapping v1",
    "source=step24 native library x cluster pseudobulk",
    "primary_resolution=r0.4",
    "",
    "core maturation and orthogonal state evaluated separately",
    "non-neutrophil contamination screened conservatively",
    "cross-project similarity uses module profile + native marker overlap",
    "library-dominated clusters retained but excluded from taxonomy discovery role",
    "",
    "NO raw RNA reread",
    "NO cross-project integration",
    "NO Cell Ranger access",
    "",
    "candidate screen is diagnostic and is NOT final annotation"
  ),
  done_file
)

cat(
  "\n===== CANDIDATE SCREEN =====\n"
)

print(
  screen[
    ,
    c(
      "project_id",
      "cluster",
      "n_cells",
      "max_library_fraction",
      "taxonomy_discovery_role",
      "provisional_core_screen",
      "provisional_state_screen",
      "provisional_scope_screen",
      "screen_reason"
    ),
    drop=FALSE
  ],
  row.names=FALSE
)

cat(
  "\n===== CROSS-PROJECT RECURRENCE =====\n"
)

print(
  recurrence[
    order(
      recurrence$dimension,
      -recurrence$n_projects,
      -recurrence$n_cells
    ),
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat(
  "\n===== NEAREST CROSS-PROJECT NEIGHBOURS =====\n"
)

print(
  nearest,
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Neutrophil cross-project program mapping completed\n"
)
