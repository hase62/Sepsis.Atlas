#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3L){
  stop(
    "Usage: script <root> <tag> <candidate_audit_dir>"
  )
}

root <- normalizePath(
  args[[1]],
  mustWork=TRUE
)

tag <- args[[2]]

audit_dir <- normalizePath(
  args[[3]],
  mustWork=TRUE
)

suppressPackageStartupMessages({
  library(SeuratObject)
})

# ============================================================
# Paths
# ============================================================

summary_file <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "monodc_project_balanced_marker_audit_v1",
  "monodc_r0p6_project_balanced_summary_v1.tsv"
)

audit_decision_file <- file.path(
  audit_dir,
  "monodc_candidate_state_provisional_decision_v1.tsv"
)

audit_done_file <- file.path(
  audit_dir,
  "MONODC_CANDIDATE_STATE_AUDIT_COMPLETE.ok"
)

clust_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "integrated_compartment_clustering",
  "Monocyte_DC_combined"
)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "monodc_annotation_freeze_v1__",
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
  paste0(
    "MONODC_ANNOTATION_FREEZE_COMPLETE__",
    tag,
    ".ok"
  )
)

if(file.exists(done_file)){
  stop(
    "Freeze already completed: ",
    done_file
  )
}

for(f in c(
  summary_file,
  audit_decision_file,
  audit_done_file
)){
  if(!file.exists(f)){
    stop(
      "Missing required input: ",
      f
    )
  }
}

if(!dir.exists(clust_dir)){
  stop(
    "Missing clustering directory: ",
    clust_dir
  )
}

# ============================================================
# Load cluster-level evidence
# ============================================================

sm <- read.delim(
  summary_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

sm$cluster <- as.character(
  sm$cluster
)

required_summary <- c(
  "cluster",
  "backbone",
  "n_cells",
  "dominant_backbone_fraction",
  "n_eligible_projects",
  "n_reproducible_markers",
  "n_projects",
  "max_project_fraction",
  "project_entropy_normalized",
  "transfer_celltype_purity",
  "healthy_fraction",
  "disease_fraction",
  "state_reproducibility_class"
)

stopifnot(
  all(
    required_summary %in%
      names(sm)
  )
)

stopifnot(
  nrow(sm) == 37L,
  sum(
    as.integer(
      sm$n_cells
    )
  ) == 192040L,
  !anyDuplicated(
    sm$cluster
  )
)

audit <- read.delim(
  audit_decision_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

audit$cluster <- as.character(
  audit$cluster
)

required_audit <- c(
  "cluster",
  "proposed_core",
  "proposed_state",
  "proposed_confidence",
  "rationale"
)

stopifnot(
  all(
    required_audit %in%
      names(audit)
  ),
  setequal(
    audit$cluster,
    c(
      "3",
      "14",
      "21",
      "30"
    )
  ),
  !anyDuplicated(
    audit$cluster
  )
)

# ============================================================
# Build frozen cluster decision table
# ============================================================

decision <- sm

decision$monodc_core_v1 <-
  as.character(
    decision$backbone
  )

decision$monodc_state_v1 <-
  "none"

decision$monodc_core_support_fraction_v1 <-
  as.numeric(
    decision$dominant_backbone_fraction
  )

decision$monodc_core_support_level_v1 <-
  ifelse(
    decision$monodc_core_support_fraction_v1 >= 0.95,
    "high",
    ifelse(
      decision$monodc_core_support_fraction_v1 >= 0.85,
      "medium",
      ifelse(
        decision$monodc_core_support_fraction_v1 >= 0.70,
        "low_medium",
        "low"
      )
    )
  )

decision$monodc_state_evidence_confidence_v1 <-
  "not_applicable"

decision$monodc_annotation_basis_v1 <-
  ifelse(
    decision$monodc_core_support_level_v1 == "low",
    "r0.6_backbone_low_support_no_reproducible_state_claim",
    "r0.6_backbone_no_reproducible_state_claim"
  )

decision$audit_rationale_v1 <-
  ""

# ============================================================
# Apply the four audited decisions
# ============================================================

for(i in seq_len(
  nrow(audit)
)){

  k <- audit$cluster[i]

  j <- match(
    k,
    decision$cluster
  )

  stopifnot(
    !is.na(j)
  )

  decision$monodc_core_v1[j] <-
    audit$proposed_core[i]

  decision$monodc_state_v1[j] <-
    audit$proposed_state[i]

  decision$audit_rationale_v1[j] <-
    audit$rationale[i]

  if(
    audit$proposed_state[i] != "none"
  ){

    decision$monodc_state_evidence_confidence_v1[j] <-
      audit$proposed_confidence[i]

    decision$monodc_annotation_basis_v1[j] <-
      if(
        decision$state_reproducibility_class[j] ==
          "cross_project_state_candidate"
      ){
        "project_balanced_cross_project_candidate_state"
      } else {
        "project_balanced_limited_support_candidate_state"
      }

  } else {

    # cluster 14:
    # reproducible DC2 program is an identity refinement,
    # not an orthogonal state.
    decision$monodc_core_support_level_v1[j] <-
      audit$proposed_confidence[i]

    decision$monodc_state_evidence_confidence_v1[j] <-
      "not_applicable"

    decision$monodc_annotation_basis_v1[j] <-
      "project_balanced_native_DC2_identity_refinement"
  }
}

# Explicit invariants.
stopifnot(
  decision$monodc_core_v1[
    decision$cluster == "14"
  ] == "DC2_like",
  decision$monodc_state_v1[
    decision$cluster == "14"
  ] == "none",
  decision$monodc_state_v1[
    decision$cluster == "3"
  ] ==
    "RETN_PADI4_myeloid_state_candidate",
  decision$monodc_state_v1[
    decision$cluster == "21"
  ] ==
    "interferon_stimulated_candidate",
  decision$monodc_state_v1[
    decision$cluster == "30"
  ] ==
    "LPAR1_PDE4B_state_candidate"
)

# cluster 26 is intentionally NOT promoted to a new subtype.
# Its r0.6 backbone is retained, with low support encoded separately.
if("26" %in% decision$cluster){

  j26 <- match(
    "26",
    decision$cluster
  )

  stopifnot(
    decision$monodc_core_support_fraction_v1[j26] <
      0.70
  )

  decision$monodc_annotation_basis_v1[j26] <-
    "r0.6_backbone_low_support_no_state_claim"
}

# ============================================================
# Write decision sheet before cell-level freeze
# ============================================================

decision_file <- file.path(
  out_dir,
  paste0(
    "monodc_r0p6_annotation_decision_v1__",
    tag,
    ".tsv"
  )
)

write.table(
  decision,
  file=decision_file,
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# Find the cell-level r0.6 clustering metadata
#
# Search ONLY:
# atlas/full_atlas_annotation_v1/
#   integrated_compartment_clustering/Monocyte_DC_combined
#
# cellranger_count is never searched.
#
# We identify r0.6 by requiring its cluster counts to match the
# frozen 37-cluster summary exactly.
# ============================================================

target_counts <- setNames(
  as.integer(
    sm$n_cells
  ),
  sm$cluster
)

N_EXPECTED <- 192040L

is_matching_cluster_vector <- function(v){

  if(
    !(is.atomic(v) ||
      is.factor(v))
  ){
    return(FALSE)
  }

  if(
    length(v) !=
      N_EXPECTED
  ){
    return(FALSE)
  }

  if(anyNA(v)){
    return(FALSE)
  }

  x <- as.character(
    v
  )

  tt <- table(x)

  if(
    length(tt) !=
      length(target_counts)
  ){
    return(FALSE)
  }

  if(
    !setequal(
      names(tt),
      names(target_counts)
    )
  ){
    return(FALSE)
  }

  all(
    as.integer(
      tt[
        names(target_counts)
      ]
    ) ==
      as.integer(
        target_counts
      )
  )
}

extract_data_frames <- function(obj){

  out <- list()

  if(
    inherits(
      obj,
      "Seurat"
    )
  ){
    out[[
      "Seurat_meta.data"
    ]] <- obj[[]]
  }

  if(
    is.data.frame(
      obj
    )
  ){
    out[[
      "object_data.frame"
    ]] <- obj
  }

  if(
    is.list(
      obj
    ) &&
    length(
      names(obj)
    )
  ){

    for(nm in names(obj)){

      z <- obj[[nm]]

      if(
        is.data.frame(
          z
        )
      ){
        out[[
          paste0(
            "list$",
            nm
          )
        ]] <- z
      }
    }
  }

  out
}

rds_files <- list.files(
  clust_dir,
  pattern="\\.rds$",
  full.names=TRUE,
  recursive=TRUE
)

if(!length(rds_files)){
  stop(
    "No RDS files found under ",
    clust_dir
  )
}

# Prefer filenames that look like clustering objects.
priority <-
  10L *
    grepl(
      "cluster",
      basename(
        rds_files
      ),
      ignore.case=TRUE
    ) +
  5L *
    grepl(
      "object|seurat|integrated",
      basename(
        rds_files
      ),
      ignore.case=TRUE
    )

sizes <- file.info(
  rds_files
)$size

rds_files <- rds_files[
  order(
    -priority,
    sizes,
    na.last=TRUE
  )
]

source_rds <- NULL
source_frame_name <- NULL
source_cluster_column <- NULL
md <- NULL

cat(
  "\n===== SEARCHING CLUSTERING RDS =====\n"
)

for(f in rds_files){

  cat(
    "checking: ",
    f,
    "\n",
    sep=""
  )

  obj <- tryCatch(
    readRDS(f),
    error=function(e){
      cat(
        "  read failed: ",
        conditionMessage(e),
        "\n",
        sep=""
      )
      NULL
    }
  )

  if(is.null(obj)){
    next
  }

  frames <- extract_data_frames(
    obj
  )

  if(length(frames)){

    for(fn in names(frames)){

      d <- frames[[fn]]

      if(
        nrow(d) !=
          N_EXPECTED
      ){
        next
      }

      good_cols <- names(d)[
        vapply(
          d,
          is_matching_cluster_vector,
          logical(1)
        )
      ]

      if(length(good_cols)){

        # If multiple columns match, they must represent the
        # same cell-wise clustering or we stop.
        if(
          length(good_cols) >
            1L
        ){

          vv <- lapply(
            good_cols,
            function(nm)
              as.character(
                d[[nm]]
              )
          )

          same <- all(
            vapply(
              vv[-1],
              identical,
              logical(1),
              vv[[1]]
            )
          )

          if(!same){
            stop(
              "Multiple non-identical columns match the ",
              "r0.6 cluster-count distribution in ",
              f,
              ": ",
              paste(
                good_cols,
                collapse=", "
              )
            )
          }

          preferred <- good_cols[
            grepl(
              "0p6|0\\.6|r0p6",
              good_cols,
              ignore.case=TRUE
            )
          ]

          if(length(preferred)){
            good_cols <- c(
              preferred,
              setdiff(
                good_cols,
                preferred
              )
            )
          }
        }

        source_rds <- f
        source_frame_name <- fn
        source_cluster_column <-
          good_cols[1]

        md <- d

        break
      }
    }
  }

  rm(
    obj,
    frames
  )

  invisible(
    gc()
  )

  if(!is.null(source_rds)){
    break
  }
}

if(is.null(source_rds)){
  stop(
    paste0(
      "Could not locate a 192,040-cell metadata frame ",
      "containing a column whose cluster counts exactly ",
      "match the frozen r0.6 summary.\n",
      "RDS files checked:\n",
      paste(
        rds_files,
        collapse="\n"
      )
    )
  )
}

cat(
  "\nsource_rds = ",
  source_rds,
  "\n",
  "metadata_frame = ",
  source_frame_name,
  "\n",
  "cluster_column = ",
  source_cluster_column,
  "\n",
  sep=""
)

stopifnot(
  nrow(md) ==
    N_EXPECTED
)

cluster_r0p6 <- as.character(
  md[[
    source_cluster_column
  ]]
)

stopifnot(
  is_matching_cluster_vector(
    cluster_r0p6
  )
)

# ============================================================
# Cell identity
# ============================================================

if(
  "global_cell" %in%
    names(md)
){

  global_cell <- as.character(
    md$global_cell
  )

} else {

  global_cell <- rownames(
    md
  )
}

if(
  is.null(global_cell) ||
  length(global_cell) !=
    N_EXPECTED ||
  anyNA(global_cell) ||
  any(
    !nzchar(
      global_cell
    )
  ) ||
  anyDuplicated(
    global_cell
  )
){
  stop(
    "Could not establish a unique global_cell identifier"
  )
}

get_char_col <- function(
  x,
  nm,
  n
){

  if(
    nm %in%
      names(x)
  ){
    return(
      as.character(
        x[[nm]]
      )
    )
  }

  rep(
    NA_character_,
    n
  )
}

# ============================================================
# Map r0.6 cluster decisions onto all 192,040 cells
# ============================================================

mi <- match(
  cluster_r0p6,
  decision$cluster
)

if(anyNA(mi)){
  stop(
    "Some cell clusters are absent from decision sheet"
  )
}

ann <- data.frame(
  global_cell=
    global_cell,

  project_id=
    get_char_col(
      md,
      "project_id",
      N_EXPECTED
    ),

  library_key=
    get_char_col(
      md,
      "library_key",
      N_EXPECTED
    ),

  condition_binary=
    get_char_col(
      md,
      "condition_binary",
      N_EXPECTED
    ),

  integration_celltype_full_v1=
    get_char_col(
      md,
      "integration_celltype_full_v1",
      N_EXPECTED
    ),

  upstream_transfer_confidence=
    get_char_col(
      md,
      "transfer_confidence",
      N_EXPECTED
    ),

  cluster_r0p6=
    cluster_r0p6,

  monodc_backbone_v1=
    decision$backbone[
      mi
    ],

  monodc_core_v1=
    decision$monodc_core_v1[
      mi
    ],

  monodc_state_v1=
    decision$monodc_state_v1[
      mi
    ],

  monodc_core_support_fraction_v1=
    decision$monodc_core_support_fraction_v1[
      mi
    ],

  monodc_core_support_level_v1=
    decision$monodc_core_support_level_v1[
      mi
    ],

  monodc_state_evidence_confidence_v1=
    decision$monodc_state_evidence_confidence_v1[
      mi
    ],

  monodc_state_reproducibility_class_v1=
    decision$state_reproducibility_class[
      mi
    ],

  monodc_annotation_basis_v1=
    decision$monodc_annotation_basis_v1[
      mi
    ],

  stringsAsFactors=FALSE,
  check.names=FALSE
)

stopifnot(
  nrow(ann) ==
    N_EXPECTED,
  !anyDuplicated(
    ann$global_cell
  ),
  !anyNA(
    ann$cluster_r0p6
  ),
  !anyNA(
    ann$monodc_core_v1
  ),
  !anyNA(
    ann$monodc_state_v1
  )
)

# Verify cell counts against the frozen r0.6 summary.
tab_check <- table(
  ann$cluster_r0p6
)

stopifnot(
  all(
    as.integer(
      tab_check[
        names(target_counts)
      ]
    ) ==
      as.integer(
        target_counts
      )
  )
)

# ============================================================
# Safe gzip writer
#
# Explicitly close gzip connection before gzip -t.
# ============================================================

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

  if(!identical(
    status,
    0L
  )){
    stop(
      "gzip integrity check failed: ",
      path
    )
  }

  invisible(
    TRUE
  )
}

annotation_tsv <- file.path(
  out_dir,
  paste0(
    "monodc_annotation_freeze_v1__",
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

core_counts <- as.data.frame(
  table(
    ann$monodc_core_v1
  ),
  stringsAsFactors=FALSE
)

names(core_counts) <- c(
  "monodc_core_v1",
  "n_cells"
)

state_counts <- as.data.frame(
  table(
    ann$monodc_state_v1
  ),
  stringsAsFactors=FALSE
)

names(state_counts) <- c(
  "monodc_state_v1",
  "n_cells"
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
    monodc_core_v1=
      ann$monodc_core_v1,
    monodc_state_v1=
      ann$monodc_state_v1,
    monodc_core_support_level_v1=
      ann$monodc_core_support_level_v1,
    monodc_state_evidence_confidence_v1=
      ann$monodc_state_evidence_confidence_v1
  ),
  FUN=sum
)

write.table(
  core_counts,
  file=file.path(
    out_dir,
    paste0(
      "monodc_core_counts_v1__",
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
      "monodc_state_counts_v1__",
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
      "monodc_core_state_counts_v1__",
      tag,
      ".tsv"
    )
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Primary frozen RDS
# ============================================================

freeze <- list(

  annotation_version=
    "Monocyte_DC_annotation_v1",

  n_cells=
    N_EXPECTED,

  annotation=
    ann,

  cluster_decision=
    decision,

  source_clustering_rds=
    normalizePath(
      source_rds,
      mustWork=TRUE
    ),

  source_clustering_metadata_frame=
    source_frame_name,

  source_r0p6_cluster_column=
    source_cluster_column,

  source_project_balanced_summary=
    normalizePath(
      summary_file,
      mustWork=TRUE
    ),

  source_candidate_audit=
    normalizePath(
      audit_decision_file,
      mustWork=TRUE
    ),

  rules=list(
    backbone_resolution=0.6,
    core_state_separation=TRUE,
    cluster14_core_refinement=
      "DC2_like",
    candidate_state_clusters=
      c(
        "3",
        "21",
        "30"
      ),
    no_MS1_promotion=
      "cluster 3 remains RETN_PADI4_myeloid_state_candidate",
    low_core_support_policy=
      paste0(
        "retain r0.6 backbone identity and encode ",
        "support separately; do not invent a new subtype"
      ),
    no_reproducible_state_policy=
      paste0(
        "state=none means no frozen reproducible ",
        "orthogonal state claim"
      )
  )
)

freeze_rds <- file.path(
  out_dir,
  paste0(
    "monodc_annotation_freeze_v1__",
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
    "Failed atomic rename of freeze RDS"
  )
}

# ============================================================
# Re-read final outputs
# ============================================================

check_rds <- readRDS(
  freeze_rds
)

stopifnot(
  identical(
    check_rds$annotation_version,
    "Monocyte_DC_annotation_v1"
  ),
  check_rds$n_cells ==
    N_EXPECTED,
  nrow(
    check_rds$annotation
  ) ==
    N_EXPECTED,
  !anyDuplicated(
    check_rds$annotation$global_cell
  )
)

check_tsv <- read.delim(
  annotation_tsv,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

stopifnot(
  nrow(check_tsv) ==
    N_EXPECTED,
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
  check_tsv,
  check_rds
)

invisible(
  gc()
)

# ============================================================
# Completion
# ============================================================

writeLines(
  c(
    "PASS",
    "Monocyte/DC annotation freeze v1",
    "annotation_version=Monocyte_DC_annotation_v1",
    paste0(
      "n_cells=",
      N_EXPECTED
    ),
    "n_clusters=37",
    "backbone_resolution=0.6",
    "cluster14=DC2_like + none",
    "cluster3=Classical_monocyte_like + RETN_PADI4_myeloid_state_candidate",
    "cluster21=Classical_monocyte_like + interferon_stimulated_candidate",
    "cluster30=Classical_monocyte_like + LPAR1_PDE4B_state_candidate",
    paste0(
      "source_clustering_rds=",
      source_rds
    ),
    paste0(
      "source_cluster_column=",
      source_cluster_column
    ),
    paste0(
      "freeze_rds=",
      freeze_rds
    ),
    paste0(
      "freeze_tsv=",
      annotation_tsv
    ),
    "gzip_integrity=PASS",
    "cell_count_recheck=PASS"
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
  "\n===== STATE COUNTS =====\n"
)

print(
  state_counts,
  row.names=FALSE
)

cat(
  "\n===== CORE x STATE =====\n"
)

print(
  joint_counts,
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Monocyte/DC annotation freeze v1 completed\n"
)

