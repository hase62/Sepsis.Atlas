#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

infile <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_annotation_v1",
  "b_plasma_transfer_cells_v1.tsv.gz"
)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_annotation_v1"
)

stopifnot(file.exists(infile))

x <- read.delim(
  gzfile(infile),
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

cat("\n===== INPUT =====\n")
cat("cells =", nrow(x), "\n")

stopifnot(
  nrow(x) == 46029,
  "integration_compartment_primary_v1" %in% names(x),
  "global_cell" %in% names(x),
  "project_id" %in% names(x),
  "library_key" %in% names(x)
)

cat("\n===== PRIMARY COMPARTMENT =====\n")
print(
  table(
    x$integration_compartment_primary_v1,
    useNA="ifany"
  )
)

bp <- x[
  x$integration_compartment_primary_v1 ==
    "B_plasma",
  ,
  drop=FALSE
]

deferred <- x[
  x$integration_compartment_primary_v1 ==
    "Deferred_transfer_low_confidence",
  ,
  drop=FALSE
]

stopifnot(
  nrow(bp) == 45955,
  nrow(deferred) == 74,
  nrow(bp) + nrow(deferred) == nrow(x),
  length(unique(bp$global_cell)) == 45955,
  length(unique(deferred$global_cell)) == 74,
  length(intersect(
    bp$global_cell,
    deferred$global_cell
  )) == 0
)

cat("\n===== PRIMARY B/PLASMA =====\n")
cat("cells     =", nrow(bp), "\n")
cat("libraries =", length(unique(bp$library_key)), "\n")
cat("projects  =", length(unique(bp$project_id)), "\n")

cat("\nCells by project:\n")
print(
  sort(
    table(bp$project_id),
    decreasing=TRUE
  )
)

cat("\nTransfer confidence:\n")
print(
  table(
    bp$transfer_confidence,
    useNA="ifany"
  )
)

cat("\nCondition:\n")
print(
  table(
    bp$condition_binary,
    useNA="ifany"
  )
)

# ------------------------------------------------------------
# Authoritative downstream B/plasma cell set
# ------------------------------------------------------------

bp_file <- file.path(
  out_dir,
  "b_plasma_primary_cells_v1.tsv.gz"
)

con <- gzfile(bp_file, "wt")

tryCatch(
  write.table(
    bp,
    con,
    sep="\t",
    quote=1,
    qmethod="double",
    row.names=FALSE
  ),
  finally=close(con)
)

# ------------------------------------------------------------
# Preserve the 74 excluded cells explicitly for later
# Deferred annotation.
# ------------------------------------------------------------

deferred_file <- file.path(
  out_dir,
  "b_plasma_full_vote_but_primary_deferred_v1.tsv.gz"
)

con <- gzfile(deferred_file, "wt")

tryCatch(
  write.table(
    deferred,
    con,
    sep="\t",
    quote=1,
    qmethod="double",
    row.names=FALSE
  ),
  finally=close(con)
)

# ------------------------------------------------------------
# Project summary
# ------------------------------------------------------------

projects <- sort(unique(bp$project_id))

project_summary <- do.call(
  rbind,
  lapply(
    projects,
    function(p){

      d <- bp[
        bp$project_id == p,
        ,
        drop=FALSE
      ]

      data.frame(
        project_id=p,
        n_cells=nrow(d),
        n_libraries=
          length(unique(d$library_key)),
        healthy_fraction=
          mean(
            d$condition_binary ==
              "healthy",
            na.rm=TRUE
          ),
        disease_fraction=
          mean(
            d$condition_binary ==
              "disease",
            na.rm=TRUE
          ),
        high_transfer_fraction=
          mean(
            d$transfer_confidence ==
              "high",
            na.rm=TRUE
          ),
        stringsAsFactors=FALSE
      )
    }
  )
)

write.table(
  project_summary,
  file=file.path(
    out_dir,
    "b_plasma_primary_project_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

writeLines(
  c(
    "PASS",
    "B/plasma primary downstream cell set",
    "n_cells=45955",
    "selection=integration_compartment_primary_v1 == B_plasma",
    "74 full-transfer B_plasma cells remain Deferred_transfer_low_confidence"
  ),
  file.path(
    out_dir,
    "B_PLASMA_PRIMARY_CELLSET_COMPLETE.ok"
  )
)

cat(
  "\nPASS: authoritative B/plasma primary cell set = ",
  nrow(bp),
  " cells\n",
  sep=""
)
