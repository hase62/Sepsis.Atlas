#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

base <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_annotation_v1"
)

pilot_file <- file.path(
  base,
  "b_plasma_pilot_cells_v1.tsv.gz"
)

primary_file <- file.path(
  base,
  "b_plasma_primary_cells_v1.tsv.gz"
)

transfer_dir <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "annotation_transfer_by_library"
)

read_gz <- function(f) {
  read.delim(
    gzfile(f),
    sep="\t",
    quote="\"",
    comment.char="",
    stringsAsFactors=FALSE,
    check.names=FALSE
  )
}

pilot <- read_gz(pilot_file)
primary <- read_gz(primary_file)

stopifnot(
  nrow(pilot) == 10018,
  nrow(primary) == 45955
)

# ============================================================
# Frozen-pilot B/plasma cells absent from current primary B/plasma
# ============================================================

in_primary <- pilot$global_cell %in%
  primary$global_cell

keep <- pilot[
  in_primary,
  ,
  drop=FALSE
]

missing <- pilot[
  !in_primary,
  ,
  drop=FALSE
]

stopifnot(
  nrow(keep) == 9997,
  nrow(missing) == 21
)

cat("\n===== PILOT / PRIMARY =====\n")
cat("pilot frozen B/plasma =", nrow(pilot), "\n")
cat("current primary       =", nrow(keep), "\n")
cat("not current primary   =", nrow(missing), "\n")

# ============================================================
# Resolve all 21 against ALL 158 transfer tables
# ============================================================

files <- sort(list.files(
  transfer_dir,
  pattern="\\.tsv\\.gz$",
  full.names=TRUE
))

stopifnot(length(files) == 158)

target <- as.character(
  missing$global_cell
)

resolved_list <- list()

for(i in seq_along(files)) {

  z <- read_gz(files[i])

  ii <- which(
    as.character(z$global_cell) %in%
      target
  )

  if(!length(ii)) next

  d <- z[ii,,drop=FALSE]
  d$source_transfer_file <-
    basename(files[i])

  resolved_list[[
    length(resolved_list) + 1L
  ]] <- d
}

if(!length(resolved_list)) {
  stop("None of the 21 cells could be resolved")
}

resolved <- do.call(
  rbind,
  resolved_list
)

# One authoritative transfer row per global cell expected.
if(anyDuplicated(resolved$global_cell)) {
  dup <- unique(
    resolved$global_cell[
      duplicated(resolved$global_cell)
    ]
  )

  stop(
    "Duplicated global_cell in transfer tables: ",
    paste(dup, collapse=", ")
  )
}

mr <- match(
  target,
  resolved$global_cell
)

cat("\n===== RESOLUTION =====\n")
cat(
  "resolved in all transfer tables = ",
  sum(!is.na(mr)),
  " / 21\n",
  sep=""
)

if(anyNA(mr)) {
  cat("\nUNRESOLVED CELLS:\n")
  print(target[is.na(mr)])
  stop(
    "Some pilot cells are absent from current transfer tables"
  )
}

resolved <- resolved[
  mr,
  ,
  drop=FALSE
]

stopifnot(
  identical(
    as.character(resolved$global_cell),
    target
  )
)

# ============================================================
# Diagnostics
# ============================================================

cat("\n===== integration_compartment_full_v1 =====\n")
print(
  sort(
    table(
      resolved$integration_compartment_full_v1,
      useNA="ifany"
    ),
    decreasing=TRUE
  )
)

cat("\n===== integration_compartment_primary_v1 =====\n")
print(
  sort(
    table(
      resolved$integration_compartment_primary_v1,
      useNA="ifany"
    ),
    decreasing=TRUE
  )
)

cat("\n===== integration_celltype_full_v1 =====\n")
print(
  sort(
    table(
      resolved$integration_celltype_full_v1,
      useNA="ifany"
    ),
    decreasing=TRUE
  )
)

cat("\n===== transfer_confidence =====\n")
print(
  sort(
    table(
      resolved$transfer_confidence,
      useNA="ifany"
    ),
    decreasing=TRUE
  )
)

cat("\n===== project =====\n")
print(
  sort(
    table(
      resolved$project_id,
      useNA="ifany"
    ),
    decreasing=TRUE
  )
)

cat("\n===== CELL-BY-CELL =====\n")

show_cols <- c(
  "global_cell",
  "project_id",
  "library_key",
  "integration_compartment_full_v1",
  "integration_compartment_vote_fraction",
  "integration_celltype_full_v1",
  "integration_celltype_vote_fraction",
  "transfer_confidence",
  "integration_compartment_primary_v1"
)

show_cols <- intersect(
  show_cols,
  names(resolved)
)

print(
  resolved[,show_cols,drop=FALSE],
  row.names=FALSE
)

# ============================================================
# Save 21-cell provenance
# ============================================================

write.table(
  resolved,
  file=file.path(
    base,
    "b_plasma_pilot_21_nonprimary_resolution_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

# ============================================================
# Create authoritative 9,997-cell subtype-discovery pilot
# ============================================================

mk <- match(
  keep$global_cell,
  primary$global_cell
)

stopifnot(!anyNA(mk))

# Append authoritative current full-atlas metadata.
for(nm in c(
  "clinical_domain_std",
  "condition_binary",
  "transfer_confidence",
  "integration_compartment_full_v1",
  "integration_compartment_vote_fraction",
  "integration_celltype_full_v1",
  "integration_celltype_vote_fraction",
  "integration_compartment_primary_v1"
)) {

  if(nm %in% names(primary)) {

    keep[[
      paste0(nm, "_full")
    ]] <- primary[[nm]][mk]
  }
}

stopifnot(
  nrow(keep) == 9997,
  length(unique(keep$global_cell)) ==
    9997,
  all(
    keep$integration_compartment_primary_v1_full ==
      "B_plasma"
  )
)

out_file <- file.path(
  base,
  "b_plasma_pilot_primary_intersection_v1.tsv.gz"
)

con <- gzfile(
  out_file,
  "wt"
)

tryCatch(
  write.table(
    keep,
    con,
    sep="\t",
    quote=1,
    qmethod="double",
    row.names=FALSE
  ),
  finally=close(con)
)

# Read-back validation.
chk <- read_gz(out_file)

stopifnot(
  nrow(chk) == 9997,
  identical(
    as.character(chk$global_cell),
    as.character(keep$global_cell)
  )
)

write.table(
  data.frame(
    category=c(
      "frozen_pilot_B_plasma",
      "current_primary_B_plasma_intersection",
      "excluded_from_current_primary_B_plasma"
    ),
    n_cells=c(
      10018,
      9997,
      21
    )
  ),
  file=file.path(
    base,
    "b_plasma_pilot_primary_intersection_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

writeLines(
  c(
    "PASS",
    "Frozen B/plasma pilot intersected with authoritative current primary B/plasma",
    "frozen_pilot=10018",
    "current_primary_intersection=9997",
    "excluded_nonprimary=21",
    "all 21 excluded cells resolved against complete annotation_transfer_by_library set"
  ),
  file.path(
    base,
    "B_PLASMA_PILOT_PRIMARY_INTERSECTION_COMPLETE_v1.ok"
  )
)

cat(
  "\nPASS: authoritative B/plasma subtype-discovery pilot = ",
  nrow(keep),
  " cells\n",
  sep=""
)

