#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 6L){
  stop(
    "Usage: script <root> <tag> <assembly_v2_dir> <assembly_v1_dir> <recon_dir> <out_dir>"
  )
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]
assembly_v2_dir <- normalizePath(args[[3]], mustWork=TRUE)
assembly_v1_dir <- normalizePath(args[[4]], mustWork=TRUE)
recon_dir <- normalizePath(args[[5]], mustWork=TRUE)
out_dir <- args[[6]]

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

N_FINAL <- 665816L
N_V1 <- 662541L
N_UPSTREAM_DEFERRED <- 3275L
N_STEP44_DEFERRED <- 13279L
N_FINAL_DEFERRED <- 9547L

N_LIBRARIES_ALL <- 159L
N_LIBRARIES_INCLUDED <- 158L
N_LIBRARIES_EXCLUDED <- 1L

MIN_FINAL_CELLS <- 200L

done_file <- file.path(
  out_dir,
  "ATLAS_GLOBAL_QC_INTEGRITY_AUDIT_V2_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop("Already completed: ", done_file)
}

# ============================================================
# Helpers
# ============================================================

read_tsv <- function(path){

  if(!file.exists(path)){
    stop("Missing: ", path)
  }

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

single_file <- function(dir, pattern){

  f <- list.files(
    dir,
    pattern=pattern,
    full.names=TRUE
  )

  if(length(f) != 1L){
    stop(
      "Expected one file in ",
      dir,
      " pattern=",
      pattern,
      " found=",
      length(f)
    )
  }

  f[[1]]
}

logical_safe <- function(x){

  if(is.logical(x)){
    return(x)
  }

  y <- toupper(trimws(as.character(x)))

  z <- rep(NA, length(y))

  z[y %in% c("TRUE","T","1","YES")] <- TRUE
  z[y %in% c("FALSE","F","0","NO")] <- FALSE

  z
}

same_na_safe <- function(a, b){

  a <- as.character(a)
  b <- as.character(b)

  (is.na(a) & is.na(b)) |
    (
      !is.na(a) &
      !is.na(b) &
      a == b
    )
}

checks <- list()

check <- function(
  name,
  pass,
  observed,
  expected,
  detail=""
){

  checks[[length(checks)+1L]] <<-
    data.frame(
      check=name,
      status=if(isTRUE(pass)) "PASS" else "FAIL",
      observed=as.character(observed),
      expected=as.character(expected),
      detail=as.character(detail),
      stringsAsFactors=FALSE
    )
}

# ============================================================
# Canonical files
# ============================================================

atlas_v2_file <- single_file(
  assembly_v2_dir,
  "^atlas_final_cell_annotation_v2__.*\\.tsv\\.gz$"
)

atlas_v1_file <- single_file(
  assembly_v1_dir,
  "^atlas_final_cell_annotation_v1__.*\\.tsv\\.gz$"
)

ledger_file <- single_file(
  recon_dir,
  "^atlas_deferred_reconciliation_cell_ledger_v1__.*\\.tsv\\.gz$"
)

fqc_file <- file.path(
  root,
  "pre_integration",
  "final_qc_library_manifest.tsv"
)

decision_file <- file.path(
  root,
  "pre_integration",
  "final_qc_decision_summary.tsv"
)

resolved_file <- file.path(
  root,
  "pre_integration",
  "pilot_unintegrated_large_v1",
  "resolved_library_table.tsv"
)

transfer_summary_file <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "full_annotation_transfer_library_summary.tsv"
)

transfer_dir <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "annotation_transfer_by_library"
)

fov_file <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "freeze_v1__20260813_1312",
  "final_object_validation_v1.tsv"
)

# ============================================================
# Load final Atlas v2
# ============================================================

atlas <- read_tsv(atlas_v2_file)

required <- c(
  "global_cell",
  "project_id",
  "library_key",
  "source_compartment",
  "frozen_core",
  "frozen_axis",
  "frozen_state",
  "final_compartment_v1",
  "reconciliation_action_v1",
  "reconciliation_status_v1",
  "atlas_core_identity_v1",
  "atlas_axis_v1",
  "atlas_state_v1",
  "annotation_resolution_policy_v1",
  "upstream_primary_compartment_v1",
  "upstream_transfer_confidence_v1"
)

stopifnot(
  all(required %in% names(atlas))
)

for(nm in required){
  atlas[[nm]] <- as.character(atlas[[nm]])
}

check(
  "atlas_n_cells",
  nrow(atlas) == N_FINAL,
  nrow(atlas),
  N_FINAL
)

check(
  "atlas_unique_global_cell",
  length(unique(atlas$global_cell)) == N_FINAL,
  length(unique(atlas$global_cell)),
  N_FINAL
)

check(
  "atlas_duplicate_global_cell",
  !anyDuplicated(atlas$global_cell),
  sum(duplicated(atlas$global_cell)),
  0
)

check(
  "atlas_missing_global_cell",
  !anyNA(atlas$global_cell) &&
    !any(atlas$global_cell == ""),
  sum(
    is.na(atlas$global_cell) |
      atlas$global_cell == ""
  ),
  0
)

check(
  "atlas_missing_project_id",
  !anyNA(atlas$project_id) &&
    !any(atlas$project_id == ""),
  sum(
    is.na(atlas$project_id) |
      atlas$project_id == ""
  ),
  0
)

check(
  "atlas_missing_library_key",
  !anyNA(atlas$library_key) &&
    !any(atlas$library_key == ""),
  sum(
    is.na(atlas$library_key) |
      atlas$library_key == ""
  ),
  0
)

# ============================================================
# Exact final-compartment counts
# ============================================================

expected_final_counts <- c(
  B_plasma=45154L,
  Deferred_unresolved=9547L,
  Erythroid=11526L,
  Monocyte_DC=192040L,
  Neutrophil=51469L,
  Platelet_megakaryocyte=47470L,
  Progenitor=27071L,
  T_NK=281539L
)

final_tab <- table(
  atlas$final_compartment_v1
)

observed_final_counts <- setNames(
  as.integer(final_tab),
  names(final_tab)
)

check(
  "final_compartment_vocabulary",
  setequal(
    names(observed_final_counts),
    names(expected_final_counts)
  ),
  paste(
    sort(names(observed_final_counts)),
    collapse=","
  ),
  paste(
    sort(names(expected_final_counts)),
    collapse=","
  )
)

check(
  "final_compartment_exact_counts",
  setequal(
    names(observed_final_counts),
    names(expected_final_counts)
  ) &&
    all(
      observed_final_counts[
        names(expected_final_counts)
      ] ==
        expected_final_counts
    ),
  paste(
    observed_final_counts[
      names(expected_final_counts)
    ],
    collapse=","
  ),
  paste(
    expected_final_counts,
    collapse=","
  )
)

# ============================================================
# Step45-v1 exact preservation
# ============================================================

v1 <- read_tsv(atlas_v1_file)

stopifnot(
  nrow(v1) == N_V1,
  !anyDuplicated(v1$global_cell)
)

vi <- match(
  v1$global_cell,
  atlas$global_cell
)

check(
  "v1_all_cells_present_in_v2",
  !anyNA(vi),
  sum(is.na(vi)),
  0
)

if(!anyNA(vi)){

  check(
    "v1_core_preserved_in_v2",
    all(
      same_na_safe(
        v1$atlas_core_identity_v1,
        atlas$atlas_core_identity_v1[vi]
      )
    ),
    sum(
      !same_na_safe(
        v1$atlas_core_identity_v1,
        atlas$atlas_core_identity_v1[vi]
      )
    ),
    0
  )

  check(
    "v1_axis_preserved_in_v2",
    all(
      same_na_safe(
        v1$atlas_axis_v1,
        atlas$atlas_axis_v1[vi]
      )
    ),
    sum(
      !same_na_safe(
        v1$atlas_axis_v1,
        atlas$atlas_axis_v1[vi]
      )
    ),
    0
  )

  check(
    "v1_state_preserved_in_v2",
    all(
      same_na_safe(
        v1$atlas_state_v1,
        atlas$atlas_state_v1[vi]
      )
    ),
    sum(
      !same_na_safe(
        v1$atlas_state_v1,
        atlas$atlas_state_v1[vi]
      )
    ),
    0
  )

  check(
    "v1_final_compartment_preserved_in_v2",
    all(
      v1$final_compartment_v1 ==
        atlas$final_compartment_v1[vi]
    ),
    sum(
      v1$final_compartment_v1 !=
        atlas$final_compartment_v1[vi]
    ),
    0
  )

  check(
    "v1_reconciliation_status_preserved_in_v2",
    all(
      v1$reconciliation_status_v1 ==
        atlas$reconciliation_status_v1[vi]
    ),
    sum(
      v1$reconciliation_status_v1 !=
        atlas$reconciliation_status_v1[vi]
    ),
    0
  )
}

# ============================================================
# Upstream-deferred cells
# ============================================================

ud <- atlas[
  atlas$reconciliation_status_v1 ==
    "upstream_deferred",
  ,
  drop=FALSE
]

check(
  "upstream_deferred_n",
  nrow(ud) == N_UPSTREAM_DEFERRED,
  nrow(ud),
  N_UPSTREAM_DEFERRED
)

check(
  "upstream_deferred_source",
  nrow(ud) == N_UPSTREAM_DEFERRED &&
    all(
      ud$source_compartment ==
        "Upstream_deferred"
    ),
  sum(
    ud$source_compartment !=
      "Upstream_deferred"
  ),
  0
)

check(
  "upstream_deferred_final_compartment",
  nrow(ud) == N_UPSTREAM_DEFERRED &&
    all(
      ud$final_compartment_v1 ==
        "Deferred_unresolved"
    ),
  sum(
    ud$final_compartment_v1 !=
      "Deferred_unresolved"
  ),
  0
)

check(
  "upstream_deferred_no_core_taxonomy",
  nrow(ud) == N_UPSTREAM_DEFERRED &&
    all(is.na(ud$atlas_core_identity_v1)),
  sum(!is.na(ud$atlas_core_identity_v1)),
  0
)

check(
  "upstream_deferred_no_axis_taxonomy",
  nrow(ud) == N_UPSTREAM_DEFERRED &&
    all(is.na(ud$atlas_axis_v1)),
  sum(!is.na(ud$atlas_axis_v1)),
  0
)

check(
  "upstream_deferred_no_state_taxonomy",
  nrow(ud) == N_UPSTREAM_DEFERRED &&
    all(is.na(ud$atlas_state_v1)),
  sum(!is.na(ud$atlas_state_v1)),
  0
)

ud_primary <- table(
  ud$upstream_primary_compartment_v1
)

check(
  "upstream_deferred_primary_counts",
  identical(
    as.integer(
      ud_primary["Deferred_ambiguous"]
    ),
    2266L
  ) &&
    identical(
      as.integer(
        ud_primary[
          "Deferred_transfer_low_confidence"
        ]
      ),
      1009L
    ),
  paste(
    names(ud_primary),
    as.integer(ud_primary),
    sep="=",
    collapse=","
  ),
  "Deferred_ambiguous=2266,Deferred_transfer_low_confidence=1009"
)

ud_conf <- table(
  ud$upstream_primary_compartment_v1,
  ud$upstream_transfer_confidence_v1
)

check(
  "upstream_deferred_confidence_structure",
  ud_conf[
    "Deferred_ambiguous",
    "high"
  ] == 2236L &&
    ud_conf[
      "Deferred_ambiguous",
      "medium"
    ] == 30L &&
    ud_conf[
      "Deferred_transfer_low_confidence",
      "low"
    ] == 1009L,
  "high_ambiguous=2236;medium_ambiguous=30;low_transfer=1009",
  "high_ambiguous=2236;medium_ambiguous=30;low_transfer=1009"
)

# ============================================================
# Annotation-resolution policies
# ============================================================

policy_tab <- table(
  atlas$annotation_resolution_policy_v1
)

check(
  "frozen_taxonomy_policy_n",
  as.integer(
    policy_tab[
      "frozen_compartment_taxonomy_preserved"
    ]
  ) == N_V1,
  as.integer(
    policy_tab[
      "frozen_compartment_taxonomy_preserved"
    ]
  ),
  N_V1
)

check(
  "upstream_deferred_policy_n",
  as.integer(
    policy_tab[
      "upstream_primary_deferred_no_compartment_taxonomy"
    ]
  ) == N_UPSTREAM_DEFERRED,
  as.integer(
    policy_tab[
      "upstream_primary_deferred_no_compartment_taxonomy"
    ]
  ),
  N_UPSTREAM_DEFERRED
)

# ============================================================
# Reconciliation status exact counts
# ============================================================

expected_status <- c(
  compartment_reassigned=2712L,
  global_deferred=6272L,
  source_retained=649262L,
  source_retained_ambiguous=4295L,
  upstream_deferred=3275L
)

status_tab <- table(
  atlas$reconciliation_status_v1
)

observed_status <- setNames(
  as.integer(status_tab),
  names(status_tab)
)

check(
  "reconciliation_status_exact_counts",
  setequal(
    names(observed_status),
    names(expected_status)
  ) &&
    all(
      observed_status[
        names(expected_status)
      ] ==
        expected_status
    ),
  paste(
    observed_status[
      names(expected_status)
    ],
    collapse=","
  ),
  paste(
    expected_status,
    collapse=","
  )
)

is_deferred_final <-
  atlas$final_compartment_v1 ==
    "Deferred_unresolved"

is_deferred_status <-
  atlas$reconciliation_status_v1 %in%
    c(
      "global_deferred",
      "upstream_deferred"
    )

check(
  "final_deferred_status_consistency",
  all(
    is_deferred_final ==
      is_deferred_status
  ),
  sum(
    is_deferred_final !=
      is_deferred_status
  ),
  0
)

check(
  "final_deferred_total",
  sum(is_deferred_final) ==
    N_FINAL_DEFERRED,
  sum(is_deferred_final),
  N_FINAL_DEFERRED
)

# ============================================================
# Step44 ledger still exact
# ============================================================

ledger <- read_tsv(
  ledger_file
)

stopifnot(
  nrow(ledger) ==
    N_STEP44_DEFERRED,
  !anyDuplicated(ledger$global_cell)
)

li <- match(
  ledger$global_cell,
  atlas$global_cell
)

check(
  "step44_ledger_all_cells_in_v2",
  !anyNA(li),
  sum(is.na(li)),
  0
)

if(!anyNA(li)){

  check(
    "step44_ledger_core_preserved",
    all(
      ledger$frozen_core ==
        atlas$frozen_core[li]
    ),
    sum(
      ledger$frozen_core !=
        atlas$frozen_core[li]
    ),
    0
  )

  expected_target <- ifelse(
    ledger$reconciliation_action ==
      "retain_source_compartment",
    ledger$source_compartment,
    ledger$reconciled_compartment
  )

  check(
    "step44_ledger_final_compartment_exact",
    all(
      expected_target ==
        atlas$final_compartment_v1[li]
    ),
    sum(
      expected_target !=
        atlas$final_compartment_v1[li]
    ),
    0
  )
}

# ============================================================
# Final-QC library manifest
# ============================================================

fqc <- read_tsv(fqc_file)

fqc$project_id <- as.character(
  fqc$project_id
)

fqc$library_key <- as.character(
  fqc$library_key
)

fqc$n_final <- as.integer(
  fqc$n_final
)

fqc$include_bool <- logical_safe(
  fqc$include_in_atlas
)

fqc$key <- paste(
  fqc$project_id,
  fqc$library_key,
  sep="|||"
)

check(
  "final_qc_manifest_n_libraries",
  nrow(fqc) == N_LIBRARIES_ALL,
  nrow(fqc),
  N_LIBRARIES_ALL
)

check(
  "final_qc_manifest_unique_libraries",
  !anyDuplicated(fqc$key),
  sum(duplicated(fqc$key)),
  0
)

check(
  "final_qc_include_parse",
  !anyNA(fqc$include_bool),
  sum(is.na(fqc$include_bool)),
  0
)

included <- fqc[
  fqc$include_bool,
  ,
  drop=FALSE
]

excluded <- fqc[
  !fqc$include_bool,
  ,
  drop=FALSE
]

check(
  "included_library_count",
  nrow(included) ==
    N_LIBRARIES_INCLUDED,
  nrow(included),
  N_LIBRARIES_INCLUDED
)

check(
  "excluded_library_count",
  nrow(excluded) ==
    N_LIBRARIES_EXCLUDED,
  nrow(excluded),
  N_LIBRARIES_EXCLUDED
)

check(
  "included_final_cells_ge_200",
  all(
    included$n_final >=
      MIN_FINAL_CELLS
  ),
  min(included$n_final),
  ">=200"
)

check(
  "excluded_final_cells_lt_200",
  nrow(excluded) == 1L &&
    all(
      excluded$n_final <
        MIN_FINAL_CELLS
    ),
  paste(
    excluded$n_final,
    collapse=","
  ),
  "<200"
)

check(
  "excluded_decision_low_final_cells",
  nrow(excluded) == 1L &&
    all(
      excluded$atlas_decision ==
        "exclude_low_final_cells"
    ),
  paste(
    excluded$atlas_decision,
    collapse=","
  ),
  "exclude_low_final_cells"
)

check(
  "included_final_cell_sum",
  sum(included$n_final) ==
    N_FINAL,
  sum(included$n_final),
  N_FINAL
)

# ============================================================
# Decision summary
# ============================================================

decision <- read_tsv(
  decision_file
)

decision$n_libraries <- as.integer(
  decision$n_libraries
)

check(
  "decision_summary_n_libraries",
  sum(decision$n_libraries) ==
    N_LIBRARIES_ALL,
  sum(decision$n_libraries),
  N_LIBRARIES_ALL
)

low_n <- sum(
  decision$n_libraries[
    decision$atlas_decision ==
      "exclude_low_final_cells"
  ]
)

check(
  "decision_summary_low_final_exclusion",
  low_n == 1L,
  low_n,
  1
)

# ============================================================
# Resolved library table
# ============================================================

resolved <- read_tsv(
  resolved_file
)

resolved$project_id <- as.character(
  resolved$project_id
)

resolved$library_key <- as.character(
  resolved$library_key
)

resolved$n_cells <- as.integer(
  resolved$n_cells
)

resolved$key <- paste(
  resolved$project_id,
  resolved$library_key,
  sep="|||"
)

check(
  "resolved_library_count",
  nrow(resolved) ==
    N_LIBRARIES_INCLUDED,
  nrow(resolved),
  N_LIBRARIES_INCLUDED
)

check(
  "resolved_library_set",
  setequal(
    resolved$key,
    included$key
  ),
  length(
    intersect(
      resolved$key,
      included$key
    )
  ),
  N_LIBRARIES_INCLUDED
)

ri <- match(
  included$key,
  resolved$key
)

check(
  "resolved_n_cells_equals_final_qc",
  !anyNA(ri) &&
    all(
      resolved$n_cells[ri] ==
        included$n_final
    ),
  if(anyNA(ri)){
    paste0(
      "missing=",
      sum(is.na(ri))
    )
  } else {
    sum(
      resolved$n_cells[ri] !=
        included$n_final
    )
  },
  0
)

check(
  "resolved_final_rds_paths_exist",
  all(
    file.exists(
      resolved$final_rds_resolved
    )
  ),
  sum(
    file.exists(
      resolved$final_rds_resolved
    )
  ),
  N_LIBRARIES_INCLUDED
)

# ============================================================
# Annotation-transfer summary
# ============================================================

ts <- read_tsv(
  transfer_summary_file
)

ts$project_id <- as.character(
  ts$project_id
)

ts$library_key <- as.character(
  ts$library_key
)

ts$n_cells <- as.integer(
  ts$n_cells
)

ts$key <- paste(
  ts$project_id,
  ts$library_key,
  sep="|||"
)

check(
  "transfer_summary_library_count",
  nrow(ts) ==
    N_LIBRARIES_INCLUDED,
  nrow(ts),
  N_LIBRARIES_INCLUDED
)

check(
  "transfer_summary_library_set",
  setequal(
    ts$key,
    included$key
  ),
  length(
    intersect(
      ts$key,
      included$key
    )
  ),
  N_LIBRARIES_INCLUDED
)

check(
  "transfer_summary_cell_sum",
  sum(ts$n_cells) ==
    N_FINAL,
  sum(ts$n_cells),
  N_FINAL
)

ti <- match(
  included$key,
  ts$key
)

check(
  "transfer_summary_per_library_counts",
  !anyNA(ti) &&
    all(
      ts$n_cells[ti] ==
        included$n_final
    ),
  if(anyNA(ti)){
    paste0(
      "missing=",
      sum(is.na(ti))
    )
  } else {
    sum(
      ts$n_cells[ti] !=
        included$n_final
    )
  },
  0
)

# ============================================================
# Atlas per-library counts
# ============================================================

alc <- aggregate(
  global_cell ~ project_id + library_key,
  data=atlas,
  FUN=length
)

names(alc)[
  names(alc) == "global_cell"
] <- "atlas_n_cells"

alc$key <- paste(
  alc$project_id,
  alc$library_key,
  sep="|||"
)

check(
  "atlas_library_count",
  nrow(alc) ==
    N_LIBRARIES_INCLUDED,
  nrow(alc),
  N_LIBRARIES_INCLUDED
)

check(
  "atlas_library_set",
  setequal(
    alc$key,
    included$key
  ),
  length(
    intersect(
      alc$key,
      included$key
    )
  ),
  N_LIBRARIES_INCLUDED
)

ai <- match(
  alc$key,
  included$key
)

alc$final_qc_n_final <-
  included$n_final[ai]

alc$count_delta <-
  alc$atlas_n_cells -
    alc$final_qc_n_final

check(
  "atlas_per_library_cell_counts",
  !anyNA(ai) &&
    all(
      alc$count_delta == 0L
    ),
  if(anyNA(ai)){
    paste0(
      "missing=",
      sum(is.na(ai))
    )
  } else {
    sum(
      alc$count_delta != 0L
    )
  },
  0
)

check(
  "atlas_min_library_cells_ge_200",
  min(alc$atlas_n_cells) >=
    MIN_FINAL_CELLS,
  min(alc$atlas_n_cells),
  ">=200"
)

# ============================================================
# Exact cell universe from annotation-transfer tables
# ============================================================

transfer_files <- sort(
  list.files(
    transfer_dir,
    pattern="\\.tsv\\.gz$",
    full.names=TRUE
  )
)

check(
  "transfer_file_count",
  length(transfer_files) ==
    N_LIBRARIES_INCLUDED,
  length(transfer_files),
  N_LIBRARIES_INCLUDED
)

lst <- vector(
  "list",
  length(transfer_files)
)

for(i in seq_along(transfer_files)){

  z <- read_tsv(
    transfer_files[[i]]
  )

  stopifnot(
    all(
      c(
        "global_cell",
        "project_id",
        "library_key"
      ) %in% names(z)
    )
  )

  lst[[i]] <- data.frame(
    global_cell=
      as.character(z$global_cell),
    project_id=
      as.character(z$project_id),
    library_key=
      as.character(z$library_key),
    stringsAsFactors=FALSE
  )

  if(
    i %% 20L == 0L ||
    i == length(transfer_files)
  ){
    cat(
      sprintf(
        "Read transfer table %d/%d\n",
        i,
        length(transfer_files)
      )
    )
  }
}

tr <- do.call(
  rbind,
  lst
)

rm(lst)
invisible(gc())

check(
  "transfer_cell_rows",
  nrow(tr) == N_FINAL,
  nrow(tr),
  N_FINAL
)

check(
  "transfer_unique_global_cell",
  !anyDuplicated(tr$global_cell),
  sum(
    duplicated(tr$global_cell)
  ),
  0
)

atlas_only <- setdiff(
  atlas$global_cell,
  tr$global_cell
)

transfer_only <- setdiff(
  tr$global_cell,
  atlas$global_cell
)

check(
  "exact_cell_universe",
  length(atlas_only) == 0L &&
    length(transfer_only) == 0L,
  paste0(
    "atlas_only=",
    length(atlas_only),
    ";transfer_only=",
    length(transfer_only)
  ),
  "atlas_only=0;transfer_only=0"
)

ui <- match(
  atlas$global_cell,
  tr$global_cell
)

check(
  "cell_level_project_id_concordance",
  !anyNA(ui) &&
    all(
      atlas$project_id ==
        tr$project_id[ui]
    ),
  if(anyNA(ui)){
    paste0(
      "missing=",
      sum(is.na(ui))
    )
  } else {
    sum(
      atlas$project_id !=
        tr$project_id[ui]
    )
  },
  0
)

check(
  "cell_level_library_key_concordance",
  !anyNA(ui) &&
    all(
      atlas$library_key ==
        tr$library_key[ui]
    ),
  if(anyNA(ui)){
    paste0(
      "missing=",
      sum(is.na(ui))
    )
  } else {
    sum(
      atlas$library_key !=
        tr$library_key[ui]
    )
  },
  0
)

# ============================================================
# Large integrated-object sanity check
# ============================================================

if(file.exists(fov_file)){

  fov <- read_tsv(fov_file)

  fov$n_cells <- as.integer(
    fov$n_cells
  )

  fov$n_unique_cells <- as.integer(
    fov$n_unique_cells
  )

  exp_large <- c(
    T_NK_combined=281508L,
    Monocyte_DC_combined=192040L
  )

  obs_large <- setNames(
    fov$n_cells,
    fov$compartment
  )

  check(
    "large_object_counts",
    all(
      obs_large[
        names(exp_large)
      ] ==
        exp_large
    ),
    paste(
      obs_large[
        names(exp_large)
      ],
      collapse=","
    ),
    paste(
      exp_large,
      collapse=","
    )
  )

  check(
    "large_object_unique_cells",
    all(
      fov$n_cells ==
        fov$n_unique_cells
    ),
    sum(
      fov$n_cells !=
        fov$n_unique_cells
    ),
    0
  )
}

# ============================================================
# Audit output
# ============================================================

assertions <- do.call(
  rbind,
  checks
)

n_fail <- sum(
  assertions$status == "FAIL"
)

universe <- data.frame(
  metric=c(
    "atlas_cells",
    "annotation_transfer_cells",
    "intersection_cells",
    "atlas_only_cells",
    "transfer_only_cells"
  ),
  value=c(
    nrow(atlas),
    nrow(tr),
    length(
      intersect(
        atlas$global_cell,
        tr$global_cell
      )
    ),
    length(atlas_only),
    length(transfer_only)
  ),
  stringsAsFactors=FALSE
)

write.table(
  assertions,
  file.path(
    out_dir,
    "atlas_global_qc_assertions_v2.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

write.table(
  universe,
  file.path(
    out_dir,
    "atlas_cell_universe_concordance_v2.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

write.table(
  alc,
  file.path(
    out_dir,
    "atlas_library_cell_count_concordance_v2.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

write.table(
  excluded[
    ,
    c(
      "project_id",
      "library_key",
      "atlas_decision",
      "n_final"
    ),
    drop=FALSE
  ],
  file.path(
    out_dir,
    "atlas_excluded_low_final_library_v2.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

final_counts <- as.data.frame(
  final_tab,
  stringsAsFactors=FALSE
)

names(final_counts) <- c(
  "final_compartment_v1",
  "n_cells"
)

write.table(
  final_counts,
  file.path(
    out_dir,
    "atlas_final_compartment_counts_qc_v2.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

status_counts <- as.data.frame(
  status_tab,
  stringsAsFactors=FALSE
)

names(status_counts) <- c(
  "reconciliation_status_v1",
  "n_cells"
)

write.table(
  status_counts,
  file.path(
    out_dir,
    "atlas_reconciliation_status_counts_qc_v2.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

manifest <- data.frame(
  source=c(
    "Step45_v2_final_atlas",
    "Step45_v1_preservation_reference",
    "Step44_reconciliation",
    "final_qc_library_manifest",
    "final_qc_decision_summary",
    "resolved_library_table",
    "full_annotation_transfer_library_summary",
    "annotation_transfer_by_library",
    "final_object_validation"
  ),
  path=c(
    atlas_v2_file,
    atlas_v1_file,
    ledger_file,
    fqc_file,
    decision_file,
    resolved_file,
    transfer_summary_file,
    transfer_dir,
    fov_file
  ),
  stringsAsFactors=FALSE
)

write.table(
  manifest,
  file.path(
    out_dir,
    "atlas_global_qc_source_manifest_v2.tsv"
  ),
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)

audit <- list(
  audit_version=
    "Atlas_global_QC_integrity_audit_v2",
  n_cells=N_FINAL,
  n_libraries=N_LIBRARIES_INCLUDED,
  n_upstream_deferred=N_UPSTREAM_DEFERRED,
  n_final_deferred=N_FINAL_DEFERRED,
  assertions=assertions,
  universe=universe,
  library_counts=alc,
  final_counts=final_counts,
  status_counts=status_counts,
  manifest=manifest
)

rds_file <- file.path(
  out_dir,
  paste0(
    "atlas_global_qc_integrity_audit_v2__",
    tag,
    ".rds"
  )
)

tmp <- paste0(
  rds_file,
  ".tmp"
)

saveRDS(
  audit,
  tmp,
  compress=TRUE
)

if(!file.rename(tmp, rds_file)){
  stop("RDS rename failed")
}

rr <- readRDS(rds_file)

stopifnot(
  rr$n_cells == N_FINAL,
  rr$n_libraries == N_LIBRARIES_INCLUDED
)

sha_files <- c(
  file.path(
    out_dir,
    "atlas_global_qc_assertions_v2.tsv"
  ),
  file.path(
    out_dir,
    "atlas_cell_universe_concordance_v2.tsv"
  ),
  file.path(
    out_dir,
    "atlas_library_cell_count_concordance_v2.tsv"
  ),
  file.path(
    out_dir,
    "atlas_excluded_low_final_library_v2.tsv"
  ),
  file.path(
    out_dir,
    "atlas_final_compartment_counts_qc_v2.tsv"
  ),
  file.path(
    out_dir,
    "atlas_reconciliation_status_counts_qc_v2.tsv"
  ),
  file.path(
    out_dir,
    "atlas_global_qc_source_manifest_v2.tsv"
  ),
  rds_file
)

sha_file <- file.path(
  out_dir,
  paste0(
    "SHA256SUMS_v2__",
    tag,
    ".txt"
  )
)

status <- system2(
  "sha256sum",
  sha_files,
  stdout=sha_file
)

if(status != 0L){
  stop("SHA256 creation failed")
}

status <- system2(
  "sha256sum",
  c("-c", sha_file),
  stdout=FALSE,
  stderr=FALSE
)

if(status != 0L){
  stop("SHA256 verification failed")
}

if(n_fail > 0L){

  cat("\n===== FAILED ASSERTIONS =====\n")

  print(
    assertions[
      assertions$status == "FAIL",
      ,
      drop=FALSE
    ],
    row.names=FALSE
  )

  stop(
    "Step46 v2 failed: ",
    n_fail,
    " assertion(s)"
  )
}

writeLines(
  c(
    "PASS",
    "Atlas global QC / integrity audit v2",
    "",
    "n_final_qc_cells=665816",
    "n_unique_global_cells=665816",
    "n_included_libraries=158",
    "n_excluded_low_final_cells_libraries=1",
    "n_compartment_annotated_cells=662541",
    "n_upstream_primary_deferred_cells=3275",
    "n_final_deferred_unresolved=9547",
    "",
    "exact final-QC cell universe=PASS",
    "per-library final-QC counts=PASS",
    "200-cell exclusion rule=PASS",
    "Step45-v1 annotation preservation=PASS",
    "Step44 reconciliation preservation=PASS",
    "upstream deferred handling=PASS",
    "project/library cell metadata concordance=PASS",
    "global_cell uniqueness=PASS",
    "annotation granularity unchanged=PASS",
    "",
    "no clustering",
    "no marker discovery",
    "no expression reanalysis",
    "no annotation refinement",
    "cellranger_count not accessed",
    "",
    "audit RDS reread=PASS",
    "SHA256 manifest=PASS",
    "all assertions=PASS"
  ),
  done_file
)

cat("\n===== ASSERTIONS =====\n")
print(assertions, row.names=FALSE)

cat("\n===== CELL UNIVERSE =====\n")
print(universe, row.names=FALSE)

cat("\n===== EXCLUDED LOW-FINAL LIBRARY =====\n")
print(
  excluded[
    ,
    c(
      "project_id",
      "library_key",
      "atlas_decision",
      "n_final"
    ),
    drop=FALSE
  ],
  row.names=FALSE
)

cat("\nOUT_DIR=", out_dir, "\n", sep="")
cat("\nPASS: Atlas global QC / integrity audit v2 completed\n")
