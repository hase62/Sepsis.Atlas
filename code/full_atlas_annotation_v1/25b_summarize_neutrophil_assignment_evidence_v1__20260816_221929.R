#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]
map_dir <- normalizePath(args[[3]], mustWork=TRUE)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "neutrophil_assignment_evidence_v1__",
    tag
  )
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

read_tsv <- function(f){
  read.delim(
    f,
    sep="\t",
    quote="\"",
    comment.char="",
    stringsAsFactors=FALSE,
    check.names=FALSE
  )
}

screen <- read_tsv(
  file.path(
    map_dir,
    "neutrophil_cluster_candidate_screen_r0p4_v1.tsv"
  )
)

prog <- read_tsv(
  file.path(
    map_dir,
    "neutrophil_cluster_program_scores_r0p4_v1.tsv"
  )
)

nn <- read_tsv(
  file.path(
    map_dir,
    "neutrophil_cross_project_nearest_neighbors_r0p4_v1.tsv"
  )
)

screen$project_id <- as.character(screen$project_id)
screen$cluster <- as.character(screen$cluster)

prog$project_id <- as.character(prog$project_id)
prog$cluster <- as.character(prog$cluster)

screen$key <- paste(
  screen$project_id,
  screen$cluster,
  sep="|||"
)

prog$key <- paste(
  prog$project_id,
  prog$cluster,
  sep="|||"
)

stopifnot(
  !anyDuplicated(screen$key)
)

# ============================================================
# Add selected module evidence in wide form
# ============================================================

wanted <- c(
  mature="mature_circulating",
  early="early_granulopoiesis",
  late="late_immature_granulopoiesis",
  IFN="interferon",
  inflammatory="inflammatory_NFkB",
  stress="immediate_early_stress",
  cycling="cycling",
  platelet="platelet",
  basophil="basophil",
  T_NK="T_NK",
  eosinophil="eosinophil"
)

res <- screen

for(prefix in names(wanted)){

  module_name <- wanted[[prefix]]

  d <- prog[
    prog$module == module_name,
    ,
    drop=FALSE
  ]

  mi <- match(
    res$key,
    d$key
  )

  res[[paste0(prefix, "_median_logCPM")]] <-
    d$median_module_logCPM[mi]

  res[[paste0(prefix, "_delta")]] <-
    d$delta_vs_project_cluster_median[mi]

  res[[paste0(prefix, "_rank")]] <-
    d$project_module_rank_fraction[mi]

  res[[paste0(prefix, "_robust_z")]] <-
    d$robust_z_within_project[mi]
}

# ============================================================
# Dominant maturation evidence
# ============================================================

maturation_cols <- c(
  mature="mature_delta",
  early="early_delta",
  late="late_delta"
)

res$dominant_maturation_module <- NA_character_
res$dominant_maturation_delta <- NA_real_
res$dominant_maturation_rank <- NA_real_

for(i in seq_len(nrow(res))){

  vv <- c(
    mature=res$mature_delta[i],
    early=res$early_delta[i],
    late=res$late_delta[i]
  )

  vv[!is.finite(vv)] <- -Inf

  if(all(vv == -Inf)){
    next
  }

  winner <- names(vv)[which.max(vv)]

  res$dominant_maturation_module[i] <- winner
  res$dominant_maturation_delta[i] <- vv[[winner]]

  rank_col <- paste0(
    winner,
    "_rank"
  )

  res$dominant_maturation_rank[i] <-
    res[[rank_col]][i]
}

# ============================================================
# Top 3 cross-project neighbours
# ============================================================

nn$source_key <- as.character(nn$source_key)
nn$neighbor_key <- as.character(nn$neighbor_key)

neighbor_screen_idx <- match(
  nn$neighbor_key,
  screen$key
)

nn$neighbor_role <-
  screen$taxonomy_discovery_role[
    neighbor_screen_idx
  ]

nn$neighbor_core <-
  screen$provisional_core_screen[
    neighbor_screen_idx
  ]

nn$neighbor_state <-
  screen$provisional_state_screen[
    neighbor_screen_idx
  ]

res$top_cross_project_neighbors <- ""

for(i in seq_len(nrow(res))){

  d <- nn[
    nn$source_key == res$key[i],
    ,
    drop=FALSE
  ]

  if(!nrow(d)){
    next
  }

  # Prefer taxonomy-discovery reference clusters.
  d$role_priority <- ifelse(
    d$neighbor_role ==
      "taxonomy_discovery_eligible",
    1,
    0
  )

  d <- d[
    order(
      -d$role_priority,
      -d$combined_similarity
    ),
    ,
    drop=FALSE
  ]

  d <- head(
    d,
    3
  )

  txt <- paste0(
    d$neighbor_key,
    ":",
    d$neighbor_core,
    ":",
    d$neighbor_state,
    ":sim=",
    sprintf("%.3f", d$combined_similarity)
  )

  res$top_cross_project_neighbors[i] <-
    paste(
      txt,
      collapse=" | "
    )
}

# ============================================================
# Compact decision evidence
# ============================================================

keep <- c(
  "project_id",
  "cluster",
  "n_cells",
  "n_libraries",
  "max_library_fraction",
  "marker_support_status",
  "taxonomy_discovery_role",
  "provisional_core_screen",
  "provisional_state_screen",
  "provisional_scope_screen",

  "mature_delta",
  "mature_rank",
  "early_delta",
  "early_rank",
  "late_delta",
  "late_rank",

  "IFN_delta",
  "IFN_rank",
  "inflammatory_delta",
  "inflammatory_rank",
  "cycling_delta",
  "cycling_rank",

  "platelet_delta",
  "platelet_rank",
  "basophil_delta",
  "basophil_rank",
  "T_NK_delta",
  "T_NK_rank",

  "dominant_maturation_module",
  "dominant_maturation_delta",
  "dominant_maturation_rank",

  "top_cross_project_neighbors"
)

compact <- res[
  ,
  keep,
  drop=FALSE
]

# Most important unresolved / assignment-only set.
needs_assignment <- compact[
  compact$taxonomy_discovery_role == "assignment_only" |
    compact$provisional_core_screen ==
      "Neutrophil_maturation_unresolved",
  ,
  drop=FALSE
]

needs_assignment <- needs_assignment[
  order(
    needs_assignment$project_id,
    -needs_assignment$n_cells
  ),
  ,
  drop=FALSE
]

# ============================================================
# State evidence across projects
# ============================================================

state_modules <- c(
  "interferon",
  "inflammatory_NFkB",
  "cycling"
)

state_top <- list()

for(p in sort(unique(prog$project_id))){

  for(m in state_modules){

    d <- prog[
      prog$project_id == p &
        prog$module == m,
      ,
      drop=FALSE
    ]

    if(!nrow(d)){
      next
    }

    d <- d[
      order(
        -d$delta_vs_project_cluster_median,
        -d$project_module_rank_fraction
      ),
      ,
      drop=FALSE
    ]

    d <- head(
      d,
      3
    )

    state_top[[
      length(state_top)+1L
    ]] <- d[
      ,
      c(
        "project_id",
        "cluster",
        "module",
        "n_libraries",
        "n_cells",
        "median_module_logCPM",
        "delta_vs_project_cluster_median",
        "project_module_rank_fraction",
        "robust_z_within_project"
      ),
      drop=FALSE
    ]
  }
}

state_top <- do.call(
  rbind,
  state_top
)

# ============================================================
# Counts
# ============================================================

role_counts <- aggregate(
  list(
    n_cells=compact$n_cells
  ),
  by=list(
    taxonomy_discovery_role=
      compact$taxonomy_discovery_role
  ),
  FUN=sum
)

# ============================================================
# Write
# ============================================================

write.table(
  compact,
  file=file.path(
    out_dir,
    "neutrophil_all_cluster_assignment_evidence_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  needs_assignment,
  file=file.path(
    out_dir,
    "neutrophil_unresolved_assignment_evidence_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  state_top,
  file=file.path(
    out_dir,
    "neutrophil_state_top_clusters_by_project_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  role_counts,
  file=file.path(
    out_dir,
    "neutrophil_discovery_role_cell_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

writeLines(
  c(
    "PASS",
    "Neutrophil assignment evidence summary v1",
    "source=step25 only",
    "no RNA reread",
    "no integration",
    "maturation and orthogonal state remain separate",
    "this is diagnostic evidence before final freeze"
  ),
  file.path(
    out_dir,
    "NEUTROPHIL_ASSIGNMENT_EVIDENCE_COMPLETE.ok"
  )
)

cat(
  "\n===== DISCOVERY ROLE COUNTS =====\n"
)

print(
  role_counts,
  row.names=FALSE
)

cat(
  "\n===== UNRESOLVED / ASSIGNMENT-ONLY =====\n"
)

print(
  needs_assignment,
  row.names=FALSE
)

cat(
  "\n===== TOP STATE MODULES BY PROJECT =====\n"
)

print(
  state_top,
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS\n"
)
