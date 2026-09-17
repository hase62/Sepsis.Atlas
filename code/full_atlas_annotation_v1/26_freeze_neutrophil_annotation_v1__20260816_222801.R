#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 5L){
  stop(
    paste(
      "Usage: script",
      "<root> <tag> <cluster_dir> <map_dir> <evidence_dir>"
    )
  )
}

root <- normalizePath(
  args[[1]],
  mustWork=TRUE
)

tag <- args[[2]]

cluster_dir <- normalizePath(
  args[[3]],
  mustWork=TRUE
)

map_dir <- normalizePath(
  args[[4]],
  mustWork=TRUE
)

evidence_dir <- normalizePath(
  args[[5]],
  mustWork=TRUE
)

N_FULL <- 53501L
N_DISCOVERY <- 53460L
N_TINY <- 41L

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "neutrophil_annotation_freeze_v1__",
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
  "NEUTROPHIL_ANNOTATION_FREEZE_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop(
    "Already completed: ",
    done_file
  )
}

# ============================================================
# Helpers
# ============================================================

read_tsv <- function(path){

  if(!file.exists(path)){
    stop(
      "Missing file: ",
      path
    )
  }

  if(grepl("\\.gz$", path)){

    con <- gzfile(
      path,
      open="rt"
    )

  } else {

    con <- file(
      path,
      open="rt"
    )
  }

  tryCatch(
    {
      read.delim(
        con,
        sep="\t",
        quote="\"",
        comment.char="",
        stringsAsFactors=FALSE,
        check.names=FALSE
      )
    },
    finally={
      close(con)
    }
  )
}

write_gz_tsv <- function(
  x,
  path
){

  con <- gzfile(
    path,
    open="wt"
  )

  tryCatch(
    {
      write.table(
        x,
        file=con,
        sep="\t",
        quote=TRUE,
        qmethod="double",
        row.names=FALSE
      )
    },
    finally={
      close(con)
    }
  )

  status <- system2(
    "gzip",
    c(
      "-t",
      path
    ),
    stdout=FALSE,
    stderr=FALSE
  )

  if(status != 0L){
    stop(
      "gzip integrity check failed: ",
      path
    )
  }

  invisible(TRUE)
}

# ============================================================
# Inputs
# ============================================================

assignment_file <- file.path(
  cluster_dir,
  "neutrophil_project_native_cluster_assignments_v1.tsv.gz"
)

tiny_file <- file.path(
  cluster_dir,
  "neutrophil_tiny_projects_deferred_v1.tsv"
)

cluster_done <- file.path(
  cluster_dir,
  "NEUTROPHIL_PROJECT_NATIVE_CLUSTERING_COMPLETE.ok"
)

map_done <- file.path(
  map_dir,
  "NEUTROPHIL_CROSS_PROJECT_PROGRAM_MAPPING_COMPLETE.ok"
)

evidence_done <- file.path(
  evidence_dir,
  "NEUTROPHIL_ASSIGNMENT_EVIDENCE_COMPLETE.ok"
)

evidence_file <- file.path(
  evidence_dir,
  "neutrophil_all_cluster_assignment_evidence_v1.tsv"
)

for(f in c(
  assignment_file,
  tiny_file,
  cluster_done,
  map_done,
  evidence_done,
  evidence_file
)){
  if(!file.exists(f)){
    stop(
      "Missing required input: ",
      f
    )
  }
}

assign <- read_tsv(
  assignment_file
)

tiny <- read_tsv(
  tiny_file
)

evidence <- read_tsv(
  evidence_file
)

stopifnot(
  nrow(assign) == N_DISCOVERY,
  nrow(tiny) == N_TINY,
  nrow(assign) + nrow(tiny) == N_FULL,
  !anyDuplicated(
    assign$global_cell
  ),
  !anyDuplicated(
    tiny$global_cell
  ),
  !any(
    assign$global_cell %in%
      tiny$global_cell
  )
)

assign$project_id <-
  as.character(
    assign$project_id
  )

assign$native_cluster_r0p4 <-
  as.character(
    assign$native_cluster_r0p4
  )

assign$key <- paste(
  assign$project_id,
  assign$native_cluster_r0p4,
  sep="|||"
)

evidence$project_id <-
  as.character(
    evidence$project_id
  )

evidence$cluster <-
  as.character(
    evidence$cluster
  )

evidence$key <- paste(
  evidence$project_id,
  evidence$cluster,
  sep="|||"
)

stopifnot(
  !anyDuplicated(
    evidence$key
  )
)

observed_keys <- sort(
  unique(
    assign$key
  )
)

stopifnot(
  length(
    observed_keys
  ) == 41L
)

# ============================================================
# Final frozen cluster decision
#
# Core identity, maturation and orthogonal state are separate.
#
# "none" means:
#   no frozen orthogonal state claim.
# It does NOT mean biological absence of state.
# ============================================================

decision <- data.frame(
  project_id=character(),
  cluster=character(),
  neutrophil_core_v1=character(),
  neutrophil_maturation_v1=character(),
  neutrophil_state_v1=character(),
  maturation_evidence_confidence_v1=character(),
  state_evidence_confidence_v1=character(),
  decision_rationale_v1=character(),
  stringsAsFactors=FALSE
)

add <- function(
  project,
  cluster,
  core="Neutrophil_like",
  maturation="unresolved",
  state="none",
  maturation_confidence="unresolved",
  state_confidence="not_applicable",
  rationale=""
){

  data.frame(
    project_id=project,
    cluster=as.character(cluster),
    neutrophil_core_v1=core,
    neutrophil_maturation_v1=maturation,
    neutrophil_state_v1=state,
    maturation_evidence_confidence_v1=
      maturation_confidence,
    state_evidence_confidence_v1=
      state_confidence,
    decision_rationale_v1=
      rationale,
    stringsAsFactors=FALSE
  )
}

d <- list()

# ============================================================
# GSE163668
# ============================================================

d[[length(d)+1L]] <- add(
  "GSE163668","0",
  maturation="unresolved",
  rationale=
    "library-dominated; no reproducible maturation program"
)

d[[length(d)+1L]] <- add(
  "GSE163668","1",
  maturation="mature",
  maturation_confidence="high",
  rationale=
    "native mature marker program"
)

d[[length(d)+1L]] <- add(
  "GSE163668","2",
  maturation="mature",
  state="interferon_stimulated",
  maturation_confidence="medium",
  state_confidence="low_medium",
  rationale=
    paste0(
      "mature module support; ",
      "strong IFN module and cross-project similarity; ",
      "library-dominated cluster"
    )
)

d[[length(d)+1L]] <- add(
  "GSE163668","3",
  maturation="mature",
  state="interferon_stimulated",
  maturation_confidence="high",
  state_confidence="medium",
  rationale=
    paste0(
      "native mature markers; ",
      "project-leading IFN module"
    )
)

d[[length(d)+1L]] <- add(
  "GSE163668","4",
  maturation="mature",
  maturation_confidence="high",
  rationale=
    "native mature marker program"
)

d[[length(d)+1L]] <- add(
  "GSE163668","5",
  core="Deferred_non_neutrophil_basophil_like",
  maturation="not_applicable",
  state="not_applicable",
  maturation_confidence="not_applicable",
  state_confidence="not_applicable",
  rationale=
    "CLC/IL3RA/FCER1A/GATA2/HDC/MS4A2/ENPP3 program"
)

d[[length(d)+1L]] <- add(
  "GSE163668","6",
  maturation="immature_transition_candidate",
  maturation_confidence="low_medium",
  rationale=
    paste0(
      "combined early and late granulopoiesis markers; ",
      "cross-project immature-neutrophil neighbors"
    )
)

d[[length(d)+1L]] <- add(
  "GSE163668","7",
  core="Deferred_non_neutrophil_T_NK_like",
  maturation="not_applicable",
  state="not_applicable",
  maturation_confidence="not_applicable",
  state_confidence="not_applicable",
  rationale=
    "CD3/TCR/T-cell marker program"
)

# ============================================================
# GSE167363
# ============================================================

d[[length(d)+1L]] <- add(
  "GSE167363","0",
  maturation="unresolved",
  rationale=
    "no decisive neutrophil maturation program"
)

d[[length(d)+1L]] <- add(
  "GSE167363","1",
  maturation="mature",
  maturation_confidence="medium",
  rationale=
    "mature module is project-leading; marker evidence weaker"
)

d[[length(d)+1L]] <- add(
  "GSE167363","2",
  maturation="early_immature",
  maturation_confidence="medium",
  rationale=
    paste0(
      "strong early-granulopoiesis module across libraries; ",
      "cross-project early/mixed immature neighbors"
    )
)

d[[length(d)+1L]] <- add(
  "GSE167363","3",
  maturation="late_immature",
  maturation_confidence="high",
  rationale=
    "LTF/LCN2/MMP8 late-immature program"
)

# ============================================================
# GSE216020
# ============================================================

d[[length(d)+1L]] <- add(
  "GSE216020","0",
  maturation="mature",
  maturation_confidence="high",
  rationale=
    "CXCR2/FCGR3B/CSF3R mature marker program"
)

d[[length(d)+1L]] <- add(
  "GSE216020","1",
  maturation="mature",
  state="interferon_stimulated",
  maturation_confidence="medium",
  state_confidence="high",
  rationale=
    paste0(
      "IFIT/ISG15/RSAD2/HERC5 IFN program; ",
      "mature neutrophil background"
    )
)

d[[length(d)+1L]] <- add(
  "GSE216020","2",
  maturation="mature",
  state="inflammatory_activation_candidate",
  maturation_confidence="low_medium",
  state_confidence="medium",
  rationale=
    paste0(
      "activation-dominant cluster without immature markers; ",
      "TNFAIP6/IRAK2/TNFAIP3/ICAM1/NFKBIA program"
    )
)

d[[length(d)+1L]] <- add(
  "GSE216020","3",
  maturation="unresolved",
  state="interferon_stimulated",
  maturation_confidence="unresolved",
  state_confidence="high",
  rationale=
    paste0(
      "very strong IFN program; ",
      "mature and late-immature maturation modules nearly tied"
    )
)

d[[length(d)+1L]] <- add(
  "GSE216020","4",
  maturation="immature_transition_candidate",
  maturation_confidence="low",
  rationale=
    paste0(
      "early and late maturation modules elevated; ",
      "strong single-library dominance prevents promotion"
    )
)

d[[length(d)+1L]] <- add(
  "GSE216020","5",
  maturation="mature",
  maturation_confidence="high",
  rationale=
    "MME/CXCR2/FCGR3B/VNN2 mature marker program"
)

d[[length(d)+1L]] <- add(
  "GSE216020","6",
  maturation="unresolved",
  state="inflammatory_activation_candidate",
  maturation_confidence="unresolved",
  state_confidence="low_medium",
  rationale=
    paste0(
      "strong inflammatory module; ",
      "cluster is strongly library-dominated"
    )
)

d[[length(d)+1L]] <- add(
  "GSE216020","7",
  maturation="late_immature",
  maturation_confidence="high",
  rationale=
    "LTF/MMP8/RETN/LCN2/CEACAM8/CAMP late-immature program"
)

d[[length(d)+1L]] <- add(
  "GSE216020","8",
  maturation="mature",
  maturation_confidence="high",
  rationale=
    "mature neutrophil marker program"
)

d[[length(d)+1L]] <- add(
  "GSE216020","9",
  maturation="late_immature",
  state="inflammatory_activation_candidate",
  maturation_confidence="low_medium",
  state_confidence="medium",
  rationale=
    paste0(
      "late-immature module and cross-project neighbor support; ",
      "strong inflammatory activation; ",
      "library dominance lowers maturation confidence"
    )
)

d[[length(d)+1L]] <- add(
  "GSE216020","10",
  core="Deferred_non_neutrophil_platelet_like",
  maturation="not_applicable",
  state="not_applicable",
  maturation_confidence="not_applicable",
  state_confidence="not_applicable",
  rationale=
    "GP1BB/PPBP/PF4/TUBB1/GP9/ITGA2B platelet program"
)

d[[length(d)+1L]] <- add(
  "GSE216020","11",
  maturation="unresolved",
  rationale=
    "strong single-library dominance and no decisive maturation markers"
)

d[[length(d)+1L]] <- add(
  "GSE216020","12",
  maturation="early_immature",
  state="cycling",
  maturation_confidence="high",
  state_confidence="high",
  rationale=
    paste0(
      "MPO/ELANE/AZU1/DEFA/MS4A3 early granulopoiesis; ",
      "RRM2/KIF11/DTL/MKI67 cycling program"
    )
)

d[[length(d)+1L]] <- add(
  "GSE216020","13",
  maturation="unresolved",
  rationale=
    paste0(
      "library-dominated; IFN module elevated but ",
      "insufficient evidence for frozen state assignment"
    )
)

d[[length(d)+1L]] <- add(
  "GSE216020","14",
  core="Deferred_non_neutrophil_basophil_like",
  maturation="not_applicable",
  state="not_applicable",
  maturation_confidence="not_applicable",
  state_confidence="not_applicable",
  rationale=
    "HDC/GATA2/IL3RA/MS4A2/ENPP3/CLC basophil program"
)

d[[length(d)+1L]] <- add(
  "GSE216020","15",
  core="Deferred_non_neutrophil_T_NK_like",
  maturation="not_applicable",
  state="not_applicable",
  maturation_confidence="not_applicable",
  state_confidence="not_applicable",
  rationale=
    "NKG7/GNLY/CD3/TRBC/T-cell-NK program"
)

# ============================================================
# GSE217906
# ============================================================

d[[length(d)+1L]] <- add(
  "GSE217906","0",
  maturation="mature",
  state="interferon_stimulated",
  maturation_confidence="high",
  state_confidence="medium",
  rationale=
    paste0(
      "mature marker program plus project-leading ",
      "IFN module"
    )
)

d[[length(d)+1L]] <- add(
  "GSE217906","1",
  maturation="mature",
  maturation_confidence="high",
  rationale=
    "native mature marker program"
)

d[[length(d)+1L]] <- add(
  "GSE217906","2",
  maturation="late_immature",
  maturation_confidence="high",
  rationale=
    "LCN2/LTF/RETN/CAMP/MMP8 late-immature program"
)

d[[length(d)+1L]] <- add(
  "GSE217906","3",
  maturation="mature",
  maturation_confidence="low_medium",
  rationale=
    paste0(
      "mature module and cross-project neighbor support; ",
      "single-library dominance limits confidence"
    )
)

d[[length(d)+1L]] <- add(
  "GSE217906","4",
  maturation="early_immature",
  maturation_confidence="high",
  rationale=
    "AZU1/ELANE/DEFA4/MPO/PRTN3 early-granulopoiesis program"
)

d[[length(d)+1L]] <- add(
  "GSE217906","5",
  maturation="mature",
  state="interferon_stimulated",
  maturation_confidence="medium",
  state_confidence="low_medium",
  rationale=
    paste0(
      "mature module support; ",
      "IFN module and cross-project IFN neighbor support"
    )
)

d[[length(d)+1L]] <- add(
  "GSE217906","6",
  maturation="unresolved",
  state="inflammatory_activation_candidate",
  maturation_confidence="unresolved",
  state_confidence="low_medium",
  rationale=
    "project-high inflammatory module; maturation remains unresolved"
)

d[[length(d)+1L]] <- add(
  "GSE217906","7",
  maturation="late_immature",
  maturation_confidence="low",
  rationale=
    paste0(
      "very strong late-immature module and immature ",
      "cross-project neighbors; single-library dominance"
    )
)

# ============================================================
# GSE242127
# ============================================================

d[[length(d)+1L]] <- add(
  "GSE242127","0",
  maturation="unresolved",
  rationale=
    "insufficient independent library support"
)

d[[length(d)+1L]] <- add(
  "GSE242127","1",
  maturation="mature",
  maturation_confidence="low",
  rationale=
    paste0(
      "mature module ranks highest; ",
      "insufficient independent library support"
    )
)

d[[length(d)+1L]] <- add(
  "GSE242127","2",
  maturation="immature_transition_candidate",
  maturation_confidence="low",
  rationale=
    paste0(
      "very strong early and late granulopoiesis modules; ",
      "only one eligible library and mixed module profile"
    )
)

# ============================================================
# GSE252331
# ============================================================

d[[length(d)+1L]] <- add(
  "GSE252331","0",
  maturation="early_immature",
  maturation_confidence="medium",
  rationale=
    "DEFA/AZU1/MPO/MS4A3 early-granulopoiesis marker program"
)

d[[length(d)+1L]] <- add(
  "GSE252331","1",
  maturation="late_immature",
  maturation_confidence="medium",
  rationale=
    "LTF/CAMP/LCN2/RETN late-immature marker program"
)

decision <- do.call(
  rbind,
  d
)

decision$key <- paste(
  decision$project_id,
  decision$cluster,
  sep="|||"
)

stopifnot(
  nrow(decision) == 41L,
  !anyDuplicated(
    decision$key
  ),
  setequal(
    decision$key,
    observed_keys
  )
)

# ============================================================
# Attach evidence / derive core confidence
# ============================================================

ei <- match(
  decision$key,
  evidence$key
)

stopifnot(
  !anyNA(ei)
)

decision$taxonomy_discovery_role_v1 <-
  as.character(
    evidence$taxonomy_discovery_role[
      ei
    ]
  )

decision$marker_support_status_v1 <-
  as.character(
    evidence$marker_support_status[
      ei
    ]
  )

decision$max_library_fraction_v1 <-
  as.numeric(
    evidence$max_library_fraction[
      ei
    ]
  )

decision$mature_delta_v1 <-
  as.numeric(
    evidence$mature_delta[
      ei
    ]
  )

decision$early_delta_v1 <-
  as.numeric(
    evidence$early_delta[
      ei
    ]
  )

decision$late_delta_v1 <-
  as.numeric(
    evidence$late_delta[
      ei
    ]
  )

decision$IFN_delta_v1 <-
  as.numeric(
    evidence$IFN_delta[
      ei
    ]
  )

decision$inflammatory_delta_v1 <-
  as.numeric(
    evidence$inflammatory_delta[
      ei
    ]
  )

decision$cycling_delta_v1 <-
  as.numeric(
    evidence$cycling_delta[
      ei
    ]
  )

decision$neutrophil_core_evidence_confidence_v1 <-
  ifelse(
    grepl(
      "^Deferred_non_neutrophil_",
      decision$neutrophil_core_v1
    ),
    "high",
    ifelse(
      decision$taxonomy_discovery_role_v1 ==
        "taxonomy_discovery_eligible",
      "high",
      "medium"
    )
  )

decision$neutrophil_annotation_scope_v1 <-
  ifelse(
    grepl(
      "^Deferred_non_neutrophil_",
      decision$neutrophil_core_v1
    ),
    "deferred_non_neutrophil",
    "Neutrophil"
  )

decision$annotation_basis_v1 <-
  ifelse(
    decision$neutrophil_annotation_scope_v1 ==
      "deferred_non_neutrophil",
    "native_marker_non_neutrophil_defer",
    ifelse(
      decision$taxonomy_discovery_role_v1 ==
        "taxonomy_discovery_eligible",
      "native_marker_and_project_relative_module",
      "cross_project_module_assignment_low_support"
    )
  )

# ============================================================
# Cell-level map: 53,460 discovery cells
# ============================================================

di <- match(
  assign$key,
  decision$key
)

stopifnot(
  !anyNA(di)
)

ann_discovery <- data.frame(
  global_cell=
    as.character(
      assign$global_cell
    ),
  project_id=
    as.character(
      assign$project_id
    ),
  library_key=
    as.character(
      assign$library_key
    ),
  condition_binary=
    as.character(
      assign$condition_binary
    ),
  native_cluster_r0p4=
    as.character(
      assign$native_cluster_r0p4
    ),

  neutrophil_core_v1=
    decision$neutrophil_core_v1[
      di
    ],

  neutrophil_maturation_v1=
    decision$neutrophil_maturation_v1[
      di
    ],

  neutrophil_state_v1=
    decision$neutrophil_state_v1[
      di
    ],

  neutrophil_annotation_scope_v1=
    decision$neutrophil_annotation_scope_v1[
      di
    ],

  neutrophil_core_evidence_confidence_v1=
    decision$neutrophil_core_evidence_confidence_v1[
      di
    ],

  neutrophil_maturation_evidence_confidence_v1=
    decision$maturation_evidence_confidence_v1[
      di
    ],

  neutrophil_state_evidence_confidence_v1=
    decision$state_evidence_confidence_v1[
      di
    ],

  neutrophil_annotation_basis_v1=
    decision$annotation_basis_v1[
      di
    ],

  stringsAsFactors=FALSE,
  check.names=FALSE
)

# ============================================================
# 41 tiny-project cells
#
# They do not define taxonomy.
# Preserve Neutrophil core but do not force maturation/state.
# ============================================================

required_tiny <- c(
  "global_cell",
  "project_id",
  "library_key",
  "condition_binary"
)

stopifnot(
  all(
    required_tiny %in%
      names(tiny)
  )
)

ann_tiny <- data.frame(
  global_cell=
    as.character(
      tiny$global_cell
    ),
  project_id=
    as.character(
      tiny$project_id
    ),
  library_key=
    as.character(
      tiny$library_key
    ),
  condition_binary=
    as.character(
      tiny$condition_binary
    ),

  native_cluster_r0p4=
    "deferred_tiny_project",

  neutrophil_core_v1=
    "Neutrophil_like",

  neutrophil_maturation_v1=
    "unresolved",

  neutrophil_state_v1=
    "none",

  neutrophil_annotation_scope_v1=
    "Neutrophil",

  neutrophil_core_evidence_confidence_v1=
    "low_medium",

  neutrophil_maturation_evidence_confidence_v1=
    "unresolved",

  neutrophil_state_evidence_confidence_v1=
    "not_applicable",

  neutrophil_annotation_basis_v1=
    "tiny_project_deferred_unresolved",

  stringsAsFactors=FALSE,
  check.names=FALSE
)

ann <- rbind(
  ann_discovery,
  ann_tiny
)

stopifnot(
  nrow(ann) == N_FULL,
  !anyDuplicated(
    ann$global_cell
  ),
  !anyNA(
    ann$neutrophil_core_v1
  ),
  !anyNA(
    ann$neutrophil_maturation_v1
  ),
  !anyNA(
    ann$neutrophil_state_v1
  )
)

# ============================================================
# Strong invariants
# ============================================================

stopifnot(
  all(
    ann$neutrophil_maturation_v1[
      ann$neutrophil_annotation_scope_v1 ==
        "deferred_non_neutrophil"
    ] ==
      "not_applicable"
  ),

  all(
    ann$neutrophil_state_v1[
      ann$neutrophil_annotation_scope_v1 ==
        "deferred_non_neutrophil"
    ] ==
      "not_applicable"
  )
)

# ============================================================
# Summary tables
# ============================================================

make_count <- function(
  x,
  name
){

  z <- as.data.frame(
    table(x),
    stringsAsFactors=FALSE
  )

  names(z) <- c(
    name,
    "n_cells"
  )

  z
}

core_counts <- make_count(
  ann$neutrophil_core_v1,
  "neutrophil_core_v1"
)

maturation_counts <- make_count(
  ann$neutrophil_maturation_v1,
  "neutrophil_maturation_v1"
)

state_counts <- make_count(
  ann$neutrophil_state_v1,
  "neutrophil_state_v1"
)

joint_counts <- aggregate(
  list(
    n_cells=
      rep(
        1L,
        nrow(ann)
      )
  ),
  by=list(
    neutrophil_core_v1=
      ann$neutrophil_core_v1,
    neutrophil_maturation_v1=
      ann$neutrophil_maturation_v1,
    neutrophil_state_v1=
      ann$neutrophil_state_v1,
    maturation_confidence=
      ann$neutrophil_maturation_evidence_confidence_v1,
    state_confidence=
      ann$neutrophil_state_evidence_confidence_v1
  ),
  FUN=sum
)

stopifnot(
  sum(core_counts$n_cells) ==
    N_FULL,
  sum(maturation_counts$n_cells) ==
    N_FULL,
  sum(state_counts$n_cells) ==
    N_FULL,
  sum(joint_counts$n_cells) ==
    N_FULL
)

# ============================================================
# Write cluster decision
# ============================================================

decision_out <- decision[
  ,
  setdiff(
    names(decision),
    "key"
  ),
  drop=FALSE
]

decision_file <- file.path(
  out_dir,
  paste0(
    "neutrophil_r0p4_annotation_decision_v1__",
    tag,
    ".tsv"
  )
)

write.table(
  decision_out,
  file=decision_file,
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# Cell annotation TSV
# ============================================================

annotation_tsv <- file.path(
  out_dir,
  paste0(
    "neutrophil_cell_annotation_freeze_v1__",
    tag,
    ".tsv.gz"
  )
)

write_gz_tsv(
  ann,
  annotation_tsv
)

# ============================================================
# Summary outputs
# ============================================================

write.table(
  core_counts,
  file=file.path(
    out_dir,
    paste0(
      "neutrophil_core_counts_v1__",
      tag,
      ".tsv"
    )
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  maturation_counts,
  file=file.path(
    out_dir,
    paste0(
      "neutrophil_maturation_counts_v1__",
      tag,
      ".tsv"
    )
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  state_counts,
  file=file.path(
    out_dir,
    paste0(
      "neutrophil_state_counts_v1__",
      tag,
      ".tsv"
    )
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  joint_counts,
  file=file.path(
    out_dir,
    paste0(
      "neutrophil_core_maturation_state_counts_v1__",
      tag,
      ".tsv"
    )
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Freeze RDS
# ============================================================

freeze <- list(

  annotation_version=
    "Neutrophil_annotation_v1",

  n_cells=
    N_FULL,

  n_project_native_clustered=
    N_DISCOVERY,

  n_tiny_project_deferred=
    N_TINY,

  annotation=
    ann,

  cluster_decision=
    decision_out,

  source_cluster_directory=
    cluster_dir,

  source_cross_project_mapping_directory=
    map_dir,

  source_assignment_evidence_directory=
    evidence_dir,

  rules=list(

    primary_resolution=
      "project-native r0.4",

    cross_project_integration_primary=
      FALSE,

    core_state_separation=
      TRUE,

    maturation_separate_from_core=
      TRUE,

    canonical_core=
      "Neutrophil_like",

    maturation_levels=
      c(
        "mature",
        "late_immature",
        "early_immature",
        "immature_transition_candidate",
        "unresolved"
      ),

    orthogonal_states=
      c(
        "none",
        "interferon_stimulated",
        "inflammatory_activation_candidate",
        "cycling"
      ),

    deferred_non_neutrophil=
      c(
        "Deferred_non_neutrophil_platelet_like",
        "Deferred_non_neutrophil_basophil_like",
        "Deferred_non_neutrophil_T_NK_like"
      ),

    tiny_project_policy=
      paste0(
        "retain Neutrophil_like core; ",
        "maturation unresolved; no state claim"
      ),

    none_state_definition=
      paste0(
        "no frozen reproducible orthogonal state claim; ",
        "not biological proof of state absence"
      )
  )
)

freeze_rds <- file.path(
  out_dir,
  paste0(
    "neutrophil_annotation_freeze_v1__",
    tag,
    ".rds"
  )
)

tmp_rds <- paste0(
  freeze_rds,
  ".tmp.",
  Sys.getpid()
)

saveRDS(
  freeze,
  tmp_rds,
  compress=TRUE
)

if(
  !file.rename(
    tmp_rds,
    freeze_rds
  )
){
  stop(
    "Atomic rename failed for freeze RDS"
  )
}

# ============================================================
# Re-read integrity checks
# ============================================================

check_rds <- readRDS(
  freeze_rds
)

stopifnot(
  identical(
    check_rds$annotation_version,
    "Neutrophil_annotation_v1"
  ),
  check_rds$n_cells ==
    N_FULL,
  nrow(
    check_rds$annotation
  ) ==
    N_FULL
)

check_tsv <- read_tsv(
  annotation_tsv
)

stopifnot(
  nrow(check_tsv) ==
    N_FULL,
  !anyDuplicated(
    check_tsv$global_cell
  ),
  identical(
    as.character(
      check_tsv$global_cell
    ),
    as.character(
      ann$global_cell
    )
  )
)

rm(
  check_rds,
  check_tsv
)

invisible(
  gc()
)

# ============================================================
# Completion marker
# ============================================================

writeLines(
  c(
    "PASS",
    "Neutrophil annotation freeze v1",
    "annotation_version=Neutrophil_annotation_v1",
    "n_cells=53501",
    "n_project_native_clustered=53460",
    "n_tiny_project_deferred=41",
    "primary_resolution=project-native-r0.4",
    "",
    "core=Neutrophil_like",
    "maturation is stored separately from core identity",
    "state is stored separately from maturation",
    "",
    "cross-project integration NOT used as primary evidence",
    "native project-level clustering used",
    "library-balanced native marker evidence used",
    "cross-project module/neighbor evidence used only for assignment",
    "",
    "IFN state supported across multiple projects",
    "inflammatory activation retained as candidate state",
    "cycling frozen only where direct marker evidence is strong",
    "",
    "tiny-project cells remain maturation-unresolved",
    "gzip_integrity=PASS",
    "cell_count_recheck=PASS",
    paste0(
      "freeze_rds=",
      freeze_rds
    ),
    paste0(
      "freeze_tsv=",
      annotation_tsv
    )
  ),
  done_file
)

cat(
  "\n===== CORE COUNTS =====\n"
)

print(
  core_counts,
  row.names=FALSE
)

cat(
  "\n===== MATURATION COUNTS =====\n"
)

print(
  maturation_counts,
  row.names=FALSE
)

cat(
  "\n===== STATE COUNTS =====\n"
)

print(
  state_counts,
  row.names=FALSE
)

cat(
  "\n===== CORE x MATURATION x STATE =====\n"
)

print(
  joint_counts[
    order(
      joint_counts$neutrophil_core_v1,
      joint_counts$neutrophil_maturation_v1,
      joint_counts$neutrophil_state_v1
    ),
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Neutrophil annotation freeze v1 completed\n"
)
