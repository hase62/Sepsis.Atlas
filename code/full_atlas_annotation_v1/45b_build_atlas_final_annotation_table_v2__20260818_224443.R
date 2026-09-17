#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3L){
  stop("Usage: script <root> <tag> <assembly_v1_dir>")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]
assembly_v1_dir <- normalizePath(args[[3]], mustWork=TRUE)

N_V1 <- 662541L
N_UPSTREAM_DEFERRED <- 3275L
N_FINAL <- 665816L

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "atlas_final_annotation_assembly_v2__",
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
  "ATLAS_FINAL_ANNOTATION_ASSEMBLY_V2_COMPLETE.ok"
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

write_gz <- function(x, path){

  con <- gzfile(path, "wt")

  tryCatch(
    write.table(
      x,
      con,
      sep="\t",
      quote=TRUE,
      qmethod="double",
      row.names=FALSE
    ),
    finally=close(con)
  )

  status <- system2(
    "gzip",
    c("-t", path),
    stdout=FALSE,
    stderr=FALSE
  )

  if(status != 0L){
    stop("gzip integrity failure: ", path)
  }
}

v1_file <- list.files(
  assembly_v1_dir,
  pattern=
    "^atlas_final_cell_annotation_v1__.*\\.tsv\\.gz$",
  full.names=TRUE
)

stopifnot(length(v1_file) == 1L)

atlas <- read_tsv(v1_file)

stopifnot(
  nrow(atlas) == N_V1,
  !anyDuplicated(atlas$global_cell)
)

# ============================================================
# Read exact upstream annotation-transfer universe
# ============================================================

transfer_dir <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "annotation_transfer_by_library"
)

transfer_files <- sort(
  list.files(
    transfer_dir,
    pattern="\\.tsv\\.gz$",
    full.names=TRUE
  )
)

stopifnot(length(transfer_files) == 158L)

lst <- vector(
  "list",
  length(transfer_files)
)

required_transfer <- c(
  "global_cell",
  "project_id",
  "library_key",
  "integration_compartment_primary_v1",
  "integration_compartment_full_v1",
  "integration_celltype_full_v1",
  "annotation_broad_full_v1",
  "transfer_confidence"
)

optional_transfer <- c(
  "integration_compartment_vote_fraction",
  "integration_celltype_vote_fraction",
  "annotation_broad_vote_fraction"
)

for(i in seq_along(transfer_files)){

  z <- read_tsv(
    transfer_files[[i]]
  )

  stopifnot(
    all(
      required_transfer %in%
        names(z)
    )
  )

  cols <- c(
    required_transfer,
    intersect(
      optional_transfer,
      names(z)
    )
  )

  lst[[i]] <- z[
    ,
    cols,
    drop=FALSE
  ]

  if(
    i %% 20L == 0L ||
    i == length(transfer_files)
  ){
    cat(
      sprintf(
        "Read transfer %d/%d\n",
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

for(nm in c(
  "global_cell",
  "project_id",
  "library_key",
  "integration_compartment_primary_v1",
  "integration_compartment_full_v1",
  "integration_celltype_full_v1",
  "annotation_broad_full_v1",
  "transfer_confidence"
)){
  tr[[nm]] <- as.character(
    tr[[nm]]
  )
}

stopifnot(
  nrow(tr) == N_FINAL,
  !anyDuplicated(tr$global_cell)
)

# ============================================================
# Identify cells intentionally deferred BEFORE the seven
# compartment-specific annotation workflows
# ============================================================

extra <- tr[
  !(tr$global_cell %in% atlas$global_cell),
  ,
  drop=FALSE
]

stopifnot(
  nrow(extra) ==
    N_UPSTREAM_DEFERRED,
  !anyDuplicated(extra$global_cell)
)

primary_counts <- table(
  extra$integration_compartment_primary_v1
)

stopifnot(
  identical(
    as.integer(
      primary_counts[
        "Deferred_ambiguous"
      ]
    ),
    2266L
  ),
  identical(
    as.integer(
      primary_counts[
        "Deferred_transfer_low_confidence"
      ]
    ),
    1009L
  ),
  sum(primary_counts) ==
    N_UPSTREAM_DEFERRED
)

# ============================================================
# Add explicit upstream evidence columns.
#
# Existing 662,541 cells are NOT reannotated.
# ============================================================

new_cols <- c(
  "upstream_primary_compartment_v1",
  "upstream_full_compartment_v1",
  "upstream_celltype_full_v1",
  "upstream_annotation_broad_full_v1",
  "upstream_transfer_confidence_v1",
  "upstream_compartment_vote_fraction_v1",
  "upstream_celltype_vote_fraction_v1",
  "upstream_annotation_broad_vote_fraction_v1"
)

for(nm in new_cols){
  atlas[[nm]] <- NA
}

# Create typed empty rows from existing table schema.
extra_out <- atlas[
  rep(
    NA_integer_,
    nrow(extra)
  ),
  ,
  drop=FALSE
]

extra_out$global_cell <-
  extra$global_cell

extra_out$project_id <-
  extra$project_id

extra_out$library_key <-
  extra$library_key

# This is provenance, not a biological compartment taxonomy.
extra_out$source_compartment <-
  "Upstream_deferred"

# No compartment-specific frozen annotation exists for these
# cells, so do NOT manufacture a frozen core/axis/state.
extra_out$frozen_core <-
  NA_character_

extra_out$frozen_axis <-
  NA_character_

extra_out$frozen_state <-
  NA_character_

extra_out$frozen_confidence <-
  NA_character_

extra_out$frozen_scope <-
  "upstream_primary_deferred"

if(
  "inventory_deferred_candidate" %in%
    names(extra_out)
){
  extra_out$inventory_deferred_candidate <-
    TRUE
}

# Final top-level placement only.
extra_out$final_compartment_v1 <-
  "Deferred_unresolved"

extra_out$reconciliation_action_v1 <-
  "retain_global_deferred"

extra_out$reconciliation_status_v1 <-
  "upstream_deferred"

extra_out$reconciliation_rationale_v1 <-
  ifelse(
    extra$integration_compartment_primary_v1 ==
      "Deferred_transfer_low_confidence",
    "upstream primary compartment transfer low confidence; no forced compartment assignment",
    "upstream primary compartment ambiguous; no forced compartment assignment"
  )

# No compartment-specific biological taxonomy is invented.
extra_out$atlas_core_identity_v1 <-
  NA_character_

extra_out$atlas_axis_v1 <-
  NA_character_

extra_out$atlas_state_v1 <-
  NA_character_

extra_out$atlas_annotation_confidence_v1 <-
  NA_character_

extra_out$annotation_resolution_policy_v1 <-
  "upstream_primary_deferred_no_compartment_taxonomy"

# Preserve upstream evidence separately.
extra_out$upstream_primary_compartment_v1 <-
  extra$integration_compartment_primary_v1

extra_out$upstream_full_compartment_v1 <-
  extra$integration_compartment_full_v1

extra_out$upstream_celltype_full_v1 <-
  extra$integration_celltype_full_v1

extra_out$upstream_annotation_broad_full_v1 <-
  extra$annotation_broad_full_v1

extra_out$upstream_transfer_confidence_v1 <-
  extra$transfer_confidence

if(
  "integration_compartment_vote_fraction" %in%
    names(extra)
){
  extra_out$upstream_compartment_vote_fraction_v1 <-
    extra$integration_compartment_vote_fraction
}

if(
  "integration_celltype_vote_fraction" %in%
    names(extra)
){
  extra_out$upstream_celltype_vote_fraction_v1 <-
    extra$integration_celltype_vote_fraction
}

if(
  "annotation_broad_vote_fraction" %in%
    names(extra)
){
  extra_out$upstream_annotation_broad_vote_fraction_v1 <-
    extra$annotation_broad_vote_fraction
}

# ============================================================
# Full corrected Atlas
# ============================================================

full <- rbind(
  atlas,
  extra_out
)

stopifnot(
  nrow(full) == N_FINAL,
  !anyDuplicated(full$global_cell),
  setequal(
    full$global_cell,
    tr$global_cell
  )
)

# Existing compartment-specific annotation is bit-for-bit
# preserved for the original 662,541 cells.
stopifnot(
  identical(
    full$global_cell[
      seq_len(N_V1)
    ],
    atlas$global_cell
  ),
  identical(
    full$atlas_core_identity_v1[
      seq_len(N_V1)
    ],
    atlas$atlas_core_identity_v1
  ),
  identical(
    full$atlas_axis_v1[
      seq_len(N_V1)
    ],
    atlas$atlas_axis_v1
  ),
  identical(
    full$atlas_state_v1[
      seq_len(N_V1)
    ],
    atlas$atlas_state_v1
  )
)

# ============================================================
# Expected final compartment counts
# ============================================================

final_counts <- as.data.frame(
  table(
    final_compartment_v1=
      full$final_compartment_v1
  ),
  stringsAsFactors=FALSE
)

expected <- c(
  B_plasma=45154L,
  Deferred_unresolved=9547L,
  Erythroid=11526L,
  Monocyte_DC=192040L,
  Neutrophil=51469L,
  Platelet_megakaryocyte=47470L,
  Progenitor=27071L,
  T_NK=281539L
)

observed <- setNames(
  as.integer(final_counts$Freq),
  as.character(
    final_counts$final_compartment_v1
  )
)

stopifnot(
  setequal(
    names(observed),
    names(expected)
  ),
  all(
    observed[
      names(expected)
    ] ==
      expected
  ),
  sum(observed) == N_FINAL
)

status_counts <- as.data.frame(
  table(
    reconciliation_status_v1=
      full$reconciliation_status_v1
  ),
  stringsAsFactors=FALSE
)

source_counts <- as.data.frame(
  table(
    source_compartment=
      full$source_compartment
  ),
  stringsAsFactors=FALSE
)

upstream_deferred_counts <- as.data.frame(
  table(
    upstream_primary_compartment_v1=
      extra_out$upstream_primary_compartment_v1,
    upstream_transfer_confidence_v1=
      extra_out$upstream_transfer_confidence_v1
  ),
  stringsAsFactors=FALSE
)

upstream_deferred_counts <-
  upstream_deferred_counts[
    upstream_deferred_counts$Freq > 0,
    ,
    drop=FALSE
  ]

# ============================================================
# Write
# ============================================================

full_file <- file.path(
  out_dir,
  paste0(
    "atlas_final_cell_annotation_v2__",
    tag,
    ".tsv.gz"
  )
)

write_gz(
  full,
  full_file
)

extra_file <- file.path(
  out_dir,
  paste0(
    "atlas_upstream_deferred_cells_v1__",
    tag,
    ".tsv.gz"
  )
)

write_gz(
  extra_out,
  extra_file
)

write.table(
  final_counts,
  file.path(
    out_dir,
    "atlas_final_compartment_counts_v2.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  status_counts,
  file.path(
    out_dir,
    "atlas_reconciliation_status_counts_v2.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  source_counts,
  file.path(
    out_dir,
    "atlas_source_compartment_counts_v2.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  upstream_deferred_counts,
  file.path(
    out_dir,
    "atlas_upstream_deferred_reason_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

source_manifest <- data.frame(
  source=c(
    "Step45_v1_compartment_resolved_assembly",
    "full_annotation_transfer_universe"
  ),
  path=c(
    assembly_v1_dir,
    transfer_dir
  ),
  stringsAsFactors=FALSE
)

write.table(
  source_manifest,
  file.path(
    out_dir,
    "atlas_final_annotation_v2_source_manifest.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# RDS
# ============================================================

rds_file <- file.path(
  out_dir,
  paste0(
    "atlas_final_cell_annotation_v2__",
    tag,
    ".rds"
  )
)

tmp <- paste0(
  rds_file,
  ".tmp"
)

saveRDS(
  full,
  tmp,
  compress=TRUE
)

if(!file.rename(
  tmp,
  rds_file
)){
  stop("RDS rename failed")
}

rr <- readRDS(rds_file)

stopifnot(
  nrow(rr) == N_FINAL,
  identical(
    rr$global_cell,
    full$global_cell
  )
)

chk <- read_tsv(
  full_file
)

stopifnot(
  nrow(chk) == N_FINAL,
  identical(
    as.character(chk$global_cell),
    as.character(full$global_cell)
  )
)

# ============================================================
# SHA256
# ============================================================

sha_files <- c(
  full_file,
  extra_file,
  rds_file,
  file.path(
    out_dir,
    "atlas_final_compartment_counts_v2.tsv"
  ),
  file.path(
    out_dir,
    "atlas_reconciliation_status_counts_v2.tsv"
  ),
  file.path(
    out_dir,
    "atlas_source_compartment_counts_v2.tsv"
  ),
  file.path(
    out_dir,
    "atlas_upstream_deferred_reason_counts_v1.tsv"
  ),
  file.path(
    out_dir,
    "atlas_final_annotation_v2_source_manifest.tsv"
  )
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
  stop("SHA creation failed")
}

status <- system2(
  "sha256sum",
  c("-c", sha_file),
  stdout=FALSE,
  stderr=FALSE
)

if(status != 0L){
  stop("SHA verification failed")
}

writeLines(
  c(
    "PASS",
    "Atlas final annotation assembly v2",
    "",
    "n_final_qc_cells=665816",
    "n_compartment_annotated_cells=662541",
    "n_upstream_primary_deferred_cells=3275",
    "n_final_deferred_unresolved=9547",
    "",
    "upstream Deferred_ambiguous=2266",
    "upstream Deferred_transfer_low_confidence=1009",
    "",
    "existing 662541 compartment-specific annotations preserved exactly",
    "upstream deferred cells not forced into biological compartments",
    "no new core subtype assigned",
    "no new axis assigned",
    "no new state assigned",
    "annotation granularity unchanged",
    "",
    "exact final-QC cell universe=665816",
    "duplicate global_cell=0",
    "cellranger_count not accessed",
    "",
    "gzip_integrity=PASS",
    "RDS_reread=PASS",
    "SHA256_manifest=PASS"
  ),
  done_file
)

cat("\n===== FINAL COMPARTMENTS V2 =====\n")
print(final_counts, row.names=FALSE)

cat("\n===== RECONCILIATION STATUS V2 =====\n")
print(status_counts, row.names=FALSE)

cat("\n===== UPSTREAM DEFERRED =====\n")
print(upstream_deferred_counts, row.names=FALSE)

cat("\nOUT_DIR=", out_dir, "\n", sep="")
cat("\nPASS: Atlas final annotation assembly v2 completed\n")
