#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

pilot_rds <- file.path(
  root,
  "pre_integration",
  "pilot_unintegrated_large_v1",
  "pilot_unintegrated_large_v1__reduced3000.rds"
)

freeze_file <- file.path(
  root,
  "pre_integration",
  "pilot_unintegrated_large_v1",
  "annotation_freeze_v2",
  "large_pilot_annotation_freeze_v2.tsv"
)

sampled_file <- file.path(
  root,
  "pre_integration",
  "pilot_unintegrated_large_v1",
  "large_pilot_sampled_cells.tsv.gz"
)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_annotation_v1"
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

for(f in c(
  pilot_rds,
  freeze_file,
  sampled_file
)){
  if(!file.exists(f)){
    stop("Missing input: ", f)
  }
}

# ============================================================
# Pilot Seurat object
# ============================================================

x <- readRDS(pilot_rds)
md <- x[[]]

stopifnot(
  ncol(x) == 154496,
  nrow(md) == 154496,
  "seurat_clusters" %in% names(md),
  "project_id" %in% names(md),
  "library_key" %in% names(md),
  "original_cell" %in% names(md)
)

# ============================================================
# Frozen global-cluster annotation
# ============================================================

fr <- read.delim(
  freeze_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

required_fr <- c(
  "cluster",
  "final_broad",
  "final_fine",
  "final_confidence",
  "integration_compartment_v2",
  "integration_celltype_v2"
)

miss <- setdiff(
  required_fr,
  names(fr)
)

if(length(miss)){
  stop(
    "Missing freeze columns: ",
    paste(miss, collapse=", ")
  )
}

fr$cluster <- as.character(fr$cluster)

bp_fr <- fr[
  fr$integration_compartment_v2 ==
    "B_plasma",
  ,
  drop=FALSE
]

cat(
  "\n===== B/PLASMA FROZEN GLOBAL CLUSTERS =====\n"
)

print(
  bp_fr,
  row.names=FALSE
)

stopifnot(
  nrow(bp_fr) == 2
)

bp_clusters <- as.character(
  bp_fr$cluster
)

cat(
  "\nFrozen B/plasma clusters: ",
  paste(bp_clusters, collapse=", "),
  "\n",
  sep=""
)

# ============================================================
# Determine B/plasma cells directly from Seurat clustering
# ============================================================

cl <- as.character(
  md$seurat_clusters
)

bp_idx <- which(
  cl %in% bp_clusters
)

cat(
  "\nB/plasma pilot cells from cluster membership = ",
  length(bp_idx),
  "\n",
  sep=""
)

# Authoritative expected pilot count from frozen compartment count
stopifnot(
  length(bp_idx) == 10018
)

# ============================================================
# Resolve global cell IDs
# ============================================================

sampled <- read.delim(
  gzfile(sampled_file),
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

stopifnot(
  nrow(sampled) == 154496,
  all(c(
    "project_id",
    "library_key",
    "original_cell",
    "global_cell"
  ) %in% names(sampled))
)

# Build exact key in both tables so row order does not matter.
object_key <- paste(
  as.character(md$project_id),
  as.character(md$library_key),
  as.character(md$original_cell),
  sep="|||"
)

sampled_key <- paste(
  as.character(sampled$project_id),
  as.character(sampled$library_key),
  as.character(sampled$original_cell),
  sep="|||"
)

m <- match(
  object_key,
  sampled_key
)

if(anyNA(m)){
  stop(
    "Could not map ",
    sum(is.na(m)),
    " pilot cells to sampled-cell manifest"
  )
}

global_cell <- as.character(
  sampled$global_cell[m]
)

stopifnot(
  length(global_cell) == 154496,
  length(unique(global_cell)) == 154496
)

# ============================================================
# Join frozen cluster annotation
# ============================================================

mi <- match(
  cl,
  fr$cluster
)

if(anyNA(mi)){
  stop(
    "Pilot cluster(s) absent from annotation freeze: ",
    paste(
      sort(unique(
        cl[is.na(mi)]
      )),
      collapse=", "
    )
  )
}

cell_annot <- data.frame(
  global_cell=global_cell,
  project_id=
    as.character(md$project_id),
  library_key=
    as.character(md$library_key),
  original_cell=
    as.character(md$original_cell),
  seurat_cluster=cl,
  final_broad=
    fr$final_broad[mi],
  final_fine=
    fr$final_fine[mi],
  final_confidence=
    fr$final_confidence[mi],
  integration_compartment_v2=
    fr$integration_compartment_v2[mi],
  integration_celltype_v2=
    fr$integration_celltype_v2[mi],
  stringsAsFactors=FALSE
)

for(nm in c(
  "clinical_domain_std",
  "study_group_std",
  "severity_label_std",
  "blood_fraction_std",
  "chemistry_std"
)){
  if(nm %in% names(md)){
    cell_annot[[nm]] <-
      md[[nm]]
  }
}

bp <- cell_annot[
  cell_annot$integration_compartment_v2 ==
    "B_plasma",
  ,
  drop=FALSE
]

stopifnot(
  nrow(bp) == 10018,
  length(unique(bp$global_cell)) == 10018,
  setequal(
    unique(bp$seurat_cluster),
    bp_clusters
  ),
  all(bp$final_fine == "B_cell")
)

# ============================================================
# Summary
# ============================================================

summary_list <- lapply(
  bp_clusters,
  function(k){

    d <- bp[
      bp$seurat_cluster == k,
      ,
      drop=FALSE
    ]

    pp <- sort(
      table(d$project_id),
      decreasing=TRUE
    )

    ff <- bp_fr[
      bp_fr$cluster == k,
      ,
      drop=FALSE
    ]

    data.frame(
      cluster=k,
      n_cells=nrow(d),
      final_broad=
        ff$final_broad[1],
      final_fine=
        ff$final_fine[1],
      final_confidence=
        ff$final_confidence[1],
      n_projects=
        length(unique(
          d$project_id
        )),
      n_libraries=
        length(unique(
          d$library_key
        )),
      top_project=
        names(pp)[1],
      top_project_fraction=
        as.numeric(pp[1]) /
        nrow(d),
      stringsAsFactors=FALSE
    )
  }
)

summary_df <- do.call(
  rbind,
  summary_list
)

stopifnot(
  sum(summary_df$n_cells) ==
    10018
)

cat(
  "\n===== B/PLASMA PILOT SUMMARY =====\n"
)

print(
  summary_df,
  row.names=FALSE
)

cat(
  "\n===== PROJECT DISTRIBUTION BY GLOBAL CLUSTER =====\n"
)

for(k in bp_clusters){

  d <- bp[
    bp$seurat_cluster == k,
    ,
    drop=FALSE
  ]

  cat(
    "\n----------------------------------------\n",
    "cluster ",
    k,
    "  n=",
    nrow(d),
    "\n",
    sep=""
  )

  print(
    sort(
      table(d$project_id),
      decreasing=TRUE
    )
  )
}

# ============================================================
# Outputs
# ============================================================

cell_file <- file.path(
  out_dir,
  "b_plasma_pilot_cells_v1.tsv.gz"
)

con <- gzfile(
  cell_file,
  "wt"
)

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

summary_file <- file.path(
  out_dir,
  "b_plasma_pilot_global_cluster_summary_v1.tsv"
)

write.table(
  summary_df,
  summary_file,
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

writeLines(
  c(
    "PASS",
    "B/plasma pilot cell manifest",
    "n_cells=10018",
    paste0(
      "global_clusters=",
      paste(
        bp_clusters,
        collapse=","
      )
    ),
    "cell count derived directly from pilot Seurat cluster membership",
    "freeze n_cells column intentionally not used because it is NA"
  ),
  file.path(
    out_dir,
    "B_PLASMA_PILOT_CELLSET_COMPLETE_v1.ok"
  )
)

cat(
  "\nPASS: B/plasma pilot manifest = ",
  nrow(bp),
  " cells\n",
  sep=""
)

