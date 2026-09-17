#!/usr/bin/env Rscript

# P05E5 fix1: explicit data.table logical-column filtering.
# P05E5: build rich public metadata for the full native Atlas and the
# corrected-reference subset.
#
# Inputs:
#   - frozen P05B cell/library metadata
#   - 158 final-QC Seurat RDS objects
#   - P05E1 metadata schema classification
#   - P05E4/P05E4a corrected-reference provenance
#
# Outputs:
#   - full 665,816-cell rich metadata
#   - 597,927-cell corrected-reference rich metadata in exact matrix row order
#   - rich library metadata
#   - metadata dictionary and provenance/audit tables
#
# No expression values are modified.
# cellranger_count/ is not accessed.

args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(if(length(args)) args[[1]] else ".", mustWork=TRUE)

suppressPackageStartupMessages({
  library(data.table)
})

tag <- format(Sys.time(), "%Y%m%d_%H%M%S")

out <- file.path(
  root,
  "publication",
  "scientific_data",
  paste0("reuse_reference_rich_metadata_v1__", tag)
)

dir.create(out, recursive=TRUE)

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

latest_dir <- function(base, pattern) {
  x <- list.dirs(base, recursive=FALSE, full.names=TRUE)
  x <- sort(x[grepl(pattern, basename(x))])
  if(!length(x)) {
    stop("No directory matching: ", pattern)
  }
  tail(x, 1L)
}

clean_atomic <- function(x) {
  if(is.factor(x)) return(as.character(x))
  if(inherits(x, "Date")) return(as.character(x))
  if(inherits(x, "POSIXt")) return(as.character(x))
  if(is.atomic(x)) return(x)

  # Non-atomic metadata are not suitable for a flat public table.
  rep(NA_character_, length(x))
}

nonmissing <- function(x) {
  if(is.character(x)) {
    !is.na(x) & nzchar(trimws(x))
  } else {
    !is.na(x)
  }
}

one_constant_value <- function(x) {
  x <- clean_atomic(x)
  keep <- nonmissing(x)
  u <- unique(x[keep])

  if(length(u) == 0L) {
    return(NA)
  }

  if(length(u) == 1L) {
    return(u[[1]])
  }

  structure(
    NA,
    conflict=TRUE,
    n_unique=length(u)
  )
}

is_internal_column_name <- function(x) {
  grepl(
    paste(
      c(
        "(^|_)path($|_)",
        "(^|_)dir($|_)",
        "(^|_)root($|_)",
        "final_rds",
        "input_root",
        "output_root",
        "fallback_filtered_dir",
        "command_hash",
        "tmp",
        "tempfile",
        "log_file",
        "working_dir"
      ),
      collapse="|"
    ),
    x,
    ignore.case=TRUE
  )
}

path_value_pattern <- paste(
  c(
    "/home/",
    "/hshare",
    "/Users/",
    "[A-Za-z]:\\\\"
  ),
  collapse="|"
)

write_gzip_tsv <- function(dt, gzfile_path) {
  plain <- sub("\\.gz$", "", gzfile_path)

  if(identical(plain, gzfile_path)) {
    stop("Output must end in .gz: ", gzfile_path)
  }

  fwrite(
    dt,
    plain,
    sep="\t",
    quote=TRUE,
    na="NA"
  )

  status <- system2(
    "gzip",
    args=c("-f", plain),
    stdout=FALSE,
    stderr=FALSE
  )

  if(status != 0L) {
    stop("gzip failed for ", plain)
  }

  if(!file.exists(gzfile_path)) {
    stop("gzip output missing: ", gzfile_path)
  }

  status <- system2(
    "gzip",
    args=c("-t", gzfile_path),
    stdout=FALSE,
    stderr=FALSE
  )

  if(status != 0L) {
    stop("gzip -t failed for ", gzfile_path)
  }
}

# -------------------------------------------------------------------------
# Canonical paths
# -------------------------------------------------------------------------

sd_base <- file.path(
  root,
  "publication",
  "scientific_data"
)

p05b <- file.path(
  sd_base,
  "processed_data_release_candidate_v1__20260819_194718"
)

p05b_cell_file <- file.path(
  p05b,
  "metadata",
  "SepsisAtlas_cell_metadata_v1__processed_data_release_candidate_v1__20260819_194718.tsv.gz"
)

p05b_lib_file <- file.path(
  p05b,
  "metadata",
  "SepsisAtlas_library_metadata_with_source_accessions_v1__processed_data_release_candidate_v1__20260819_194718.tsv"
)

stopifnot(
  file.exists(p05b_cell_file),
  file.exists(p05b_lib_file)
)

p05e1 <- latest_dir(
  sd_base,
  "reuse_reference_design_v1__[0-9]{8}_[0-9]{6}$"
)

p05e4 <- latest_dir(
  sd_base,
  "reuse_reference_corrected_expression_v1__[0-9]{8}_[0-9]{6}$"
)

p05e4a <- latest_dir(
  sd_base,
  "reuse_reference_corrected_compact_v1__[0-9]{8}_[0-9]{6}$"
)

schema_file <- file.path(
  p05e1,
  "rich_metadata_schema_candidates.tsv"
)

corrected_cells_file <- file.path(
  p05e4a,
  "metadata",
  "SepsisAtlas_corrected_reference_combined_cell_order_v1.tsv"
)

p05e4_excluded_file <- file.path(
  p05e4,
  "corrected_reference_excluded_cells.tsv"
)

stopifnot(
  file.exists(schema_file),
  file.exists(corrected_cells_file)
)

# -------------------------------------------------------------------------
# Frozen P05B metadata
# -------------------------------------------------------------------------

p05b_cells <- fread(p05b_cell_file)
p05b_libs <- fread(p05b_lib_file)

stopifnot(
  nrow(p05b_cells) == 665816L,
  nrow(p05b_libs) == 158L,
  !anyDuplicated(p05b_cells$global_cell),
  all(c("project_id","library_key") %in% names(p05b_cells)),
  all(c("project_id","library_key") %in% names(p05b_libs))
)

p05b_cell_order <- p05b_cells$global_cell

# -------------------------------------------------------------------------
# P05E1 candidate classification
# -------------------------------------------------------------------------

schema <- fread(schema_file)

required_schema <- c(
  "column",
  "source",
  "metadata_location"
)

stopifnot(all(required_schema %in% names(schema)))

fq_schema <- schema[
  source == "final_qc_rds"
]

cell_candidates <- unique(
  fq_schema[
    metadata_location == "CELL_RICH",
    column
  ]
)

library_candidates <- unique(
  fq_schema[
    metadata_location == "LIBRARY_RICH",
    column
  ]
)

cell_candidates <- cell_candidates[
  !is_internal_column_name(cell_candidates)
]

library_candidates <- library_candidates[
  !is_internal_column_name(library_candidates)
]

# P05B frozen fields take precedence if names overlap.
cell_candidates <- setdiff(
  cell_candidates,
  names(p05b_cells)
)

library_candidates <- setdiff(
  library_candidates,
  names(p05b_libs)
)

# -------------------------------------------------------------------------
# final-QC RDS map
# -------------------------------------------------------------------------

rds_files <- Sys.glob(
  file.path(
    root,
    "pre_integration",
    "GSE*",
    "rds_final_qc",
    "*.rds"
  )
)

stopifnot(length(rds_files) == 158L)

rds_map <- rbindlist(
  lapply(rds_files, function(f) {
    data.table(
      project_id=basename(dirname(dirname(f))),
      library_key=sub(
        "__final_qc\\.rds$",
        "",
        basename(f)
      ),
      rds=f
    )
  })
)

stopifnot(
  nrow(rds_map) == 158L,
  !anyDuplicated(
    paste(
      rds_map$project_id,
      rds_map$library_key,
      sep="|||"
    )
  )
)

# -------------------------------------------------------------------------
# Extract cell-varying metadata and library-constant metadata
# -------------------------------------------------------------------------

cell_rows <- vector("list", nrow(rds_map))
lib_rows <- vector("list", nrow(rds_map))
conflicts <- list()
unsupported_columns <- list()

for(i in seq_len(nrow(rds_map))) {

  pid <- rds_map$project_id[[i]]
  key <- rds_map$library_key[[i]]
  f <- rds_map$rds[[i]]

  obj <- readRDS(f)
  md <- obj@meta.data

  if(is.null(rownames(md)) || any(!nzchar(rownames(md)))) {
    stop("Missing cell barcodes in metadata: ", key)
  }

  global_cell <- paste0(
    key,
    "___",
    rownames(md)
  )

  # ----- Cell-varying metadata -----
  present_cell <- intersect(
    cell_candidates,
    names(md)
  )

  ct <- data.table(
    global_cell=global_cell
  )

  for(nm in present_cell) {
    x <- md[[nm]]

    if(!is.atomic(x) && !is.factor(x) &&
       !inherits(x, "Date") && !inherits(x, "POSIXt")) {
      unsupported_columns[[length(unsupported_columns)+1L]] <-
        data.table(
          project_id=pid,
          library_key=key,
          column=nm,
          scope="CELL_RICH",
          class=paste(class(x), collapse=";")
        )
      next
    }

    set(
      ct,
      j=nm,
      value=clean_atomic(x)
    )
  }

  cell_rows[[i]] <- ct

  # ----- Library-constant metadata -----
  lr <- data.table(
    project_id=pid,
    library_key=key
  )

  present_lib <- intersect(
    library_candidates,
    names(md)
  )

  for(nm in present_lib) {
    x <- md[[nm]]

    if(!is.atomic(x) && !is.factor(x) &&
       !inherits(x, "Date") && !inherits(x, "POSIXt")) {
      unsupported_columns[[length(unsupported_columns)+1L]] <-
        data.table(
          project_id=pid,
          library_key=key,
          column=nm,
          scope="LIBRARY_RICH",
          class=paste(class(x), collapse=";")
        )
      next
    }

    xx <- clean_atomic(x)
    keep <- nonmissing(xx)
    u <- unique(xx[keep])

    if(length(u) > 1L) {
      conflicts[[length(conflicts)+1L]] <- data.table(
        project_id=pid,
        library_key=key,
        column=nm,
        n_unique_nonmissing=length(u),
        example_values=paste(
          head(as.character(u), 5L),
          collapse=" || "
        )
      )
      next
    }

    val <- if(length(u) == 1L) {
      u[[1]]
    } else {
      NA
    }

    set(
      lr,
      j=nm,
      value=val
    )
  }

  lib_rows[[i]] <- lr

  rm(obj, md, ct, lr)
  invisible(gc())

  cat(
    "READ ",
    i, "/", nrow(rds_map),
    " ",
    pid, "::", key,
    "\n",
    sep=""
  )
}

cell_extra <- rbindlist(
  cell_rows,
  use.names=TRUE,
  fill=TRUE
)

lib_extra <- rbindlist(
  lib_rows,
  use.names=TRUE,
  fill=TRUE
)

stopifnot(
  nrow(cell_extra) == 665816L,
  !anyDuplicated(cell_extra$global_cell),
  nrow(lib_extra) == 158L
)

conflict_dt <- if(length(conflicts)) {
  rbindlist(conflicts, fill=TRUE)
} else {
  data.table(
    project_id=character(),
    library_key=character(),
    column=character(),
    n_unique_nonmissing=integer(),
    example_values=character()
  )
}

unsupported_dt <- if(length(unsupported_columns)) {
  unique(rbindlist(unsupported_columns, fill=TRUE))
} else {
  data.table(
    project_id=character(),
    library_key=character(),
    column=character(),
    scope=character(),
    class=character()
  )
}

# -------------------------------------------------------------------------
# Assemble rich library metadata
#
# Precedence:
#   P05B canonical public library metadata > final-QC library-constant fields
# -------------------------------------------------------------------------

rich_lib <- merge(
  p05b_libs,
  lib_extra,
  by=c("project_id","library_key"),
  all.x=TRUE,
  sort=FALSE
)

stopifnot(
  nrow(rich_lib) == 158L,
  !anyDuplicated(
    paste(
      rich_lib$project_id,
      rich_lib$library_key,
      sep="|||"
    )
  )
)

# -------------------------------------------------------------------------
# Assemble full rich cell metadata
#
# P05B frozen annotation fields remain authoritative.
# -------------------------------------------------------------------------

rich_cell <- merge(
  p05b_cells,
  cell_extra,
  by="global_cell",
  all.x=TRUE,
  sort=FALSE
)

idx <- match(
  p05b_cell_order,
  rich_cell$global_cell
)

stopifnot(!anyNA(idx))

rich_cell <- rich_cell[idx]

stopifnot(
  identical(
    rich_cell$global_cell,
    p05b_cell_order
  )
)

# -------------------------------------------------------------------------
# Propagate useful canonical library-level fields into cell metadata
# -------------------------------------------------------------------------

preferred_library_to_cell <- c(
  "atlas_library_id",
  "database_accession",
  "sample_id",
  "analysis_sample_id",
  "analysis_subject_id",
  "analysis_subject_id_source",
  "clinical_domain_std",
  "study_group",
  "study_group_std",
  "is_sepsis",
  "is_sepsis_std",
  "cov19",
  "cov19_std",
  "ards",
  "ards_std",
  "septic_shock",
  "septic_shock_std",
  "severity",
  "severity_label_std",
  "days",
  "sampling_day",
  "timepoint_label",
  "sex",
  "sex_std",
  "age",
  "age_years",
  "sofa",
  "sofa_score",
  "peep",
  "peep_cm_h2o",
  "tissue_label",
  "tissue_std",
  "source_blood_fraction",
  "blood_fraction_std",
  "frozen_or_fresh",
  "preservation_std",
  "chemistry_version",
  "chemistry_std",
  "cellranger_or_equivalent_version",
  "reference_genome",
  "reference_genome_std",
  "analysis_batch_project",
  "analysis_batch_chemistry",
  "analysis_batch_source",
  "source_GSE",
  "source_GSM",
  "source_PRJNA",
  "source_SAMN",
  "source_SRX",
  "source_SRR"
)

propagate <- intersect(
  preferred_library_to_cell,
  names(rich_lib)
)

# Do not overwrite existing cell-level fields.
propagate <- setdiff(
  propagate,
  names(rich_cell)
)

if(length(propagate)) {
  lib_propagate <- rich_lib[
    ,
    c(
      "project_id",
      "library_key",
      propagate
    ),
    with=FALSE
  ]

  rich_cell <- merge(
    rich_cell,
    lib_propagate,
    by=c("project_id","library_key"),
    all.x=TRUE,
    sort=FALSE
  )

  idx <- match(
    p05b_cell_order,
    rich_cell$global_cell
  )

  stopifnot(!anyNA(idx))
  rich_cell <- rich_cell[idx]
}

stopifnot(
  nrow(rich_cell) == 665816L,
  identical(
    rich_cell$global_cell,
    p05b_cell_order
  )
)

# -------------------------------------------------------------------------
# Corrected-reference membership/provenance
# -------------------------------------------------------------------------

corrected_cells <- fread(corrected_cells_file)

stopifnot(
  nrow(corrected_cells) == 597927L,
  !anyDuplicated(corrected_cells$global_cell)
)

corrected_set <- corrected_cells$global_cell

rich_cell[, corrected_reference_included :=
  global_cell %in% corrected_set
]

rich_cell[, corrected_reference_method :=
  fifelse(
    corrected_reference_included,
    "scMerge2",
    NA_character_
  )
]

rich_cell[, corrected_reference_scope :=
  fifelse(
    corrected_reference_included,
    "within_compartment",
    NA_character_
  )
]

ruvK_map <- c(
  T_NK=5L,
  Monocyte_DC=2L,
  B_plasma=5L,
  Neutrophil=3L,
  Platelet_megakaryocyte=3L
)

rich_cell[, corrected_reference_ruvK :=
  as.integer(
    ruvK_map[
      final_compartment_v1
    ]
  )
]

rich_cell[
  corrected_reference_included == FALSE,
  corrected_reference_ruvK := NA_integer_
]

rich_cell[, corrected_reference_status :=
  fifelse(
    corrected_reference_included,
    "included_scMerge2",
    fifelse(
      final_compartment_v1 == "Erythroid",
      "native_only_condition_not_identifiable",
      fifelse(
        final_compartment_v1 == "Progenitor",
        "native_only_structurally_confound",
        fifelse(
          final_compartment_v1 == "Deferred_unresolved",
          "excluded_deferred_unresolved",
          "excluded_from_corrected_reference"
        )
      )
    )
  )
]

# Add exact P05E4 exclusion reason where available.
if(file.exists(p05e4_excluded_file)) {

  exc <- fread(p05e4_excluded_file)

  if(all(c(
    "global_cell",
    "exclusion_reason"
  ) %in% names(exc))) {

    exc <- unique(
      exc[, .(
        global_cell,
        corrected_reference_exclusion_reason=
          exclusion_reason
      )],
      by="global_cell"
    )

    rich_cell <- merge(
      rich_cell,
      exc,
      by="global_cell",
      all.x=TRUE,
      sort=FALSE
    )

    idx <- match(
      p05b_cell_order,
      rich_cell$global_cell
    )

    stopifnot(!anyNA(idx))
    rich_cell <- rich_cell[idx]
  }
}

if(!"corrected_reference_exclusion_reason" %in% names(rich_cell)) {
  rich_cell[, corrected_reference_exclusion_reason := NA_character_]
}

rich_cell[
  corrected_reference_included == TRUE,
  corrected_reference_exclusion_reason := NA_character_
]

# -------------------------------------------------------------------------
# Path-value audit BEFORE export
# -------------------------------------------------------------------------

audit_path_values <- function(dt, scope) {
  char_cols <- names(dt)[
    vapply(dt, is.character, logical(1))
  ]

  hits <- list()

  for(nm in char_cols) {
    x <- dt[[nm]]
    idx <- which(
      !is.na(x) &
      grepl(
        path_value_pattern,
        x,
        perl=TRUE
      )
    )

    if(length(idx)) {
      hits[[length(hits)+1L]] <- data.table(
        scope=scope,
        column=nm,
        n_hits=length(idx),
        example=head(x[idx], 1L)
      )
    }
  }

  if(length(hits)) {
    rbindlist(hits)
  } else {
    data.table(
      scope=character(),
      column=character(),
      n_hits=integer(),
      example=character()
    )
  }
}

path_audit_lib <- audit_path_values(
  rich_lib,
  "library"
)

path_audit_cell <- audit_path_values(
  rich_cell,
  "cell"
)

path_audit <- rbindlist(
  list(
    path_audit_lib,
    path_audit_cell
  ),
  fill=TRUE
)

fwrite(
  path_audit,
  file.path(
    out,
    "rich_metadata_private_path_audit.tsv"
  ),
  sep="\t",
  quote=TRUE
)

if(nrow(path_audit)) {
  print(path_audit)
  stop(
    "Private/local path values remain in rich metadata. ",
    "Review rich_metadata_private_path_audit.tsv"
  )
}

# -------------------------------------------------------------------------
# Corrected-reference metadata in exact float32/h5ad row order
# -------------------------------------------------------------------------

rich_corrected <- rich_cell[
  match(
    corrected_cells$global_cell,
    global_cell
  )
]

stopifnot(
  nrow(rich_corrected) == 597927L,
  identical(
    rich_corrected$global_cell,
    corrected_cells$global_cell
  ),
  all(rich_corrected$corrected_reference_included)
)

# -------------------------------------------------------------------------
# Metadata dictionary
# -------------------------------------------------------------------------

make_dictionary <- function(dt, scope, source_note) {

  rbindlist(
    lapply(names(dt), function(nm) {
      x <- dt[[nm]]
      keep <- nonmissing(x)

      data.table(
        scope=scope,
        column=nm,
        storage_class=paste(class(x), collapse=";"),
        n_rows=length(x),
        n_nonmissing=sum(keep),
        missing_fraction=
          1 - sum(keep)/length(x),
        n_unique_nonmissing=
          length(unique(x[keep])),
        provenance=source_note
      )
    })
  )
}

dict <- rbindlist(
  list(
    make_dictionary(
      rich_cell,
      "cell",
      paste(
        "P05B frozen cell metadata + selected final-QC",
        "cell-varying metadata + selected canonical library fields"
      )
    ),
    make_dictionary(
      rich_lib,
      "library",
      paste(
        "P05B canonical library metadata + selected final-QC",
        "library-constant metadata"
      )
    )
  ),
  fill=TRUE
)

fwrite(
  dict,
  file.path(
    out,
    "SepsisAtlas_rich_metadata_dictionary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  na="NA"
)

# -------------------------------------------------------------------------
# Export
# -------------------------------------------------------------------------

cell_out <- file.path(
  out,
  "SepsisAtlas_rich_cell_metadata_v1.tsv.gz"
)

corrected_out <- file.path(
  out,
  "SepsisAtlas_corrected_reference_rich_cell_metadata_v1.tsv.gz"
)

library_out <- file.path(
  out,
  "SepsisAtlas_rich_library_metadata_v1.tsv.gz"
)

write_gzip_tsv(
  rich_cell,
  cell_out
)

write_gzip_tsv(
  rich_corrected,
  corrected_out
)

write_gzip_tsv(
  rich_lib,
  library_out
)

# -------------------------------------------------------------------------
# Re-read and exact count/order verification
# -------------------------------------------------------------------------

cell_check <- fread(cell_out)
corrected_check <- fread(corrected_out)
lib_check <- fread(library_out)

stopifnot(
  nrow(cell_check) == 665816L,
  nrow(corrected_check) == 597927L,
  nrow(lib_check) == 158L,
  identical(
    cell_check$global_cell,
    rich_cell$global_cell
  ),
  identical(
    corrected_check$global_cell,
    corrected_cells$global_cell
  ),
  identical(
    paste(
      lib_check$project_id,
      lib_check$library_key,
      sep="|||"
    ),
    paste(
      rich_lib$project_id,
      rich_lib$library_key,
      sep="|||"
    )
  )
)

# -------------------------------------------------------------------------
# Coverage / provenance outputs
# -------------------------------------------------------------------------

coverage <- dict[
  ,
  .(
    scope,
    column,
    storage_class,
    n_rows,
    n_nonmissing,
    missing_fraction,
    n_unique_nonmissing
  )
]

fwrite(
  coverage,
  file.path(
    out,
    "rich_metadata_column_coverage.tsv"
  ),
  sep="\t"
)

fwrite(
  conflict_dt,
  file.path(
    out,
    "rich_library_metadata_conflicts.tsv"
  ),
  sep="\t",
  quote=TRUE
)

fwrite(
  unsupported_dt,
  file.path(
    out,
    "rich_metadata_unsupported_columns.tsv"
  ),
  sep="\t",
  quote=TRUE
)

status_counts <- rich_cell[
  ,
  .N,
  by=corrected_reference_status
][order(-N)]

fwrite(
  status_counts,
  file.path(
    out,
    "corrected_reference_status_counts.tsv"
  ),
  sep="\t"
)

# -------------------------------------------------------------------------
# SHA256
# -------------------------------------------------------------------------

oldwd <- getwd()
setwd(out)

files <- sort(
  list.files(
    ".",
    recursive=TRUE,
    full.names=FALSE
  )
)

files <- files[
  files != "SHA256SUMS_P05E5.txt"
]

sha <- vapply(
  files,
  function(f) {
    z <- system2(
      "sha256sum",
      args=f,
      stdout=TRUE,
      stderr=TRUE
    )

    st <- attr(z, "status")

    if(!is.null(st) && st != 0L) {
      stop("sha256sum failed for ", f)
    }

    z[[1]]
  },
  character(1)
)

writeLines(
  sha,
  "SHA256SUMS_P05E5.txt"
)

chk <- system2(
  "sha256sum",
  args=c("-c", "SHA256SUMS_P05E5.txt"),
  stdout=TRUE,
  stderr=TRUE
)

st <- attr(chk, "status")

if(!is.null(st) && st != 0L) {
  cat(paste(chk, collapse="\n"), "\n")
  stop("P05E5 SHA256 verification failed")
}

setwd(oldwd)

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------

summary <- c(
  "===== P05E5 RICH METADATA =====",
  "status=PASS",
  "full_atlas_cells=665816",
  "corrected_reference_cells=597927",
  "libraries=158",
  paste0(
    "rich_cell_metadata_columns=",
    ncol(rich_cell)
  ),
  paste0(
    "rich_library_metadata_columns=",
    ncol(rich_lib)
  ),
  paste0(
    "final_qc_cell_rich_candidate_columns=",
    length(cell_candidates)
  ),
  paste0(
    "final_qc_library_rich_candidate_columns=",
    length(library_candidates)
  ),
  paste0(
    "library_constant_conflicts=",
    nrow(conflict_dt)
  ),
  paste0(
    "unsupported_nonatomic_metadata_records=",
    nrow(unsupported_dt)
  ),
  "P05B_cell_metadata_precedence=YES",
  "P05B_library_metadata_precedence=YES",
  "corrected_reference_row_order_match=PASS",
  "gzip_integrity=PASS",
  "full_cell_metadata_reread_count_order=PASS",
  "corrected_cell_metadata_reread_count_order=PASS",
  "library_metadata_reread_count_order=PASS",
  "private_path_audit=0",
  "SHA256=PASS",
  "expression_modified=NO",
  "annotation_modified=NO",
  "P05B_modified=NO",
  "cellranger_count=NOT_ACCESSED"
)

summary_file <- file.path(
  out,
  "P05E5_SUMMARY.txt"
)

writeLines(
  summary,
  summary_file
)

# -------------------------------------------------------------------------
# Root handoff copies
# -------------------------------------------------------------------------

handoff <- c(
  "P05E5_SUMMARY.txt",
  "rich_metadata_column_coverage.tsv",
  "corrected_reference_status_counts.tsv",
  "rich_library_metadata_conflicts.tsv"
)

for(f in handoff) {
  src <- file.path(out, f)
  dst <- file.path(
    root,
    paste0(
      "P05E5_",
      sub("^P05E5_", "", f)
    )
  )
  file.copy(src, dst, overwrite=TRUE)
}

cat(readLines(summary_file), sep="\n")

cat("\n\n===== CORRECTED REFERENCE STATUS =====\n")
print(status_counts)

cat("\nOUT_DIR=", out, "\n", sep="")
cat("Root handoff files copied: ", length(handoff), "\n", sep="")
