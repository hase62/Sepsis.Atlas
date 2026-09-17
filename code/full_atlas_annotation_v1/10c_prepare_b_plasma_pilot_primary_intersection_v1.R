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

fullvote_file <- file.path(
  base,
  "b_plasma_transfer_cells_v1.tsv.gz"
)

for(f in c(
  pilot_file,
  primary_file,
  fullvote_file
)){
  if(!file.exists(f)){
    stop("Missing input: ", f)
  }
}

read_gz <- function(f){
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
fullvote <- read_gz(fullvote_file)

stopifnot(
  nrow(pilot) == 10018,
  nrow(primary) == 45955,
  nrow(fullvote) == 46029
)

# ============================================================
# Pilot cells absent from authoritative primary B/plasma set
# ============================================================

in_primary <- pilot$global_cell %in%
  primary$global_cell

missing <- pilot[
  !in_primary,
  ,
  drop=FALSE
]

keep <- pilot[
  in_primary,
  ,
  drop=FALSE
]

cat("\n===== PILOT / PRIMARY INTERSECTION =====\n")
cat("pilot total      =", nrow(pilot), "\n")
cat("in primary       =", nrow(keep), "\n")
cat("absent primary   =", nrow(missing), "\n")

# ============================================================
# Resolve the missing cells against full-transfer B/plasma vote
# ============================================================

mm <- match(
  missing$global_cell,
  fullvote$global_cell
)

cat("\n===== MISSING CELL RESOLUTION =====\n")
cat(
  "found in full-vote B/plasma table =",
  sum(!is.na(mm)),
  "/",
  nrow(missing),
  "\n"
)

resolved <- fullvote[
  mm[!is.na(mm)],
  ,
  drop=FALSE
]

if(nrow(resolved)){

  cat("\nPrimary compartment:\n")
  print(
    table(
      resolved$integration_compartment_primary_v1,
      useNA="ifany"
    )
  )

  cat("\nTransfer confidence:\n")
  print(
    table(
      resolved$transfer_confidence,
      useNA="ifany"
    )
  )

  cat("\nProject:\n")
  print(
    sort(
      table(resolved$project_id),
      decreasing=TRUE
    )
  )
}

# ============================================================
# Safety checks
# ============================================================

stopifnot(
  nrow(missing) == 21,
  nrow(keep) == 9997,
  !anyNA(mm),
  all(
    resolved$integration_compartment_primary_v1 ==
      "Deferred_transfer_low_confidence"
  ),
  length(unique(keep$global_cell)) == 9997
)

# ============================================================
# Add authoritative full-primary metadata to retained cells
# ============================================================

mk <- match(
  keep$global_cell,
  primary$global_cell
)

stopifnot(!anyNA(mk))

for(nm in c(
  "clinical_domain_std",
  "condition_binary",
  "transfer_confidence",
  "integration_compartment_full_v1",
  "integration_compartment_primary_v1",
  "integration_celltype_full_v1"
)){
  if(nm %in% names(primary)){
    keep[[paste0(nm, "_full")]] <-
      primary[[nm]][mk]
  }
}

# ============================================================
# Save authoritative B/plasma discovery pilot
# ============================================================

out_file <- file.path(
  base,
  "b_plasma_pilot_primary_intersection_v1.tsv.gz"
)

con <- gzfile(out_file, "wt")

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

excluded <- cbind(
  missing,
  resolved[
    ,
    c(
      "integration_compartment_full_v1",
      "integration_compartment_primary_v1",
      "transfer_confidence"
    ),
    drop=FALSE
  ]
)

write.table(
  excluded,
  file=file.path(
    base,
    "b_plasma_pilot_excluded_primary_deferred_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

summary <- data.frame(
  category=c(
    "pilot_frozen_B_plasma",
    "pilot_primary_intersection",
    "pilot_primary_deferred_excluded"
  ),
  n_cells=c(
    10018,
    9997,
    21
  )
)

write.table(
  summary,
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
    "B/plasma discovery pilot intersected with authoritative full primary compartment",
    "pilot_frozen_B_plasma=10018",
    "pilot_primary_intersection=9997",
    "excluded_primary_deferred=21"
  ),
  file.path(
    base,
    "B_PLASMA_PILOT_PRIMARY_INTERSECTION_COMPLETE.ok"
  )
)

cat(
  "\nPASS: authoritative B/plasma discovery pilot = ",
  nrow(keep),
  " cells\n",
  sep=""
)
