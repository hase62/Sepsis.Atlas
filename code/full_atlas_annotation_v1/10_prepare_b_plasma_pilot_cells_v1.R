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

for(f in c(pilot_rds, freeze_file, sampled_file)){
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
# Frozen cluster annotation
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
  "integration_celltype_v2",
  "n_cells"
)

miss <- setdiff(required_fr, names(fr))

if(length(miss)){
  stop(
    "Missing freeze columns: ",
    paste(miss, collapse=", ")
  )
}

fr$cluster <- as.character(fr$cluster)

bp_fr <- fr[
  fr$integration_compartment_v2 == "B_plasma",
  ,
  drop=FALSE
]

cat("\n===== B/PLASMA FROZEN GLOBAL CLUSTERS =====\n")
print(bp_fr, row.names=FALSE)

stopifnot(
  nrow(bp_fr) == 2,
  sum(bp_fr$n_cells) == 10018
)

bp_clusters <- as.character(bp_fr$cluster)

# ============================================================
# Resolve global cell IDs robustly
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
  "global_cell" %in% names(sampled)
)

object_cells <- colnames(x)

if(setequal(
  object_cells,
  sampled$global_cell
)){
  global_cell <- object_cells

} else {

  reconstructed <- paste(
    md$library_key,
    md$original_cell,
    sep="___"
  )

  if(!setequal(
    reconstructed,
    sampled$global_cell
  )){
    stop(
      "Could not reconcile pilot object cells with sampled-cell manifest"
    )
  }

  global_cell <- reconstructed
}

# ============================================================
# Join frozen cluster annotation to every pilot cell
# ============================================================

cl <- as.character(md$seurat_clusters)

mi <- match(cl, fr$cluster)

if(anyNA(mi)){
  stop(
    "Pilot cluster(s) absent from annotation freeze: ",
    paste(
      sort(unique(cl[is.na(mi)])),
      collapse=", "
    )
  )
}

cell_annot <- data.frame(
  global_cell=global_cell,
  project_id=as.character(md$project_id),
  library_key=as.character(md$library_key),
  original_cell=as.character(md$original_cell),
  seurat_cluster=cl,
  final_broad=fr$final_broad[mi],
  final_fine=fr$final_fine[mi],
  final_confidence=fr$final_confidence[mi],
  integration_compartment_v2=
    fr$integration_compartment_v2[mi],
  integration_celltype_v2=
    fr$integration_celltype_v2[mi],
  stringsAsFactors=FALSE
)

# bring clinical metadata when available
for(nm in c(
  "clinical_domain_std",
  "study_group_std",
  "severity_label_std"
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
  )
)

# ============================================================
# Summary by frozen global cluster
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
      freeze_n_cells=ff$n_cells[1],
      final_broad=ff$final_broad[1],
      final_fine=ff$final_fine[1],
      final_confidence=
        ff$final_confidence[1],
      n_projects=
        length(unique(d$project_id)),
      n_libraries=
        length(unique(d$library_key)),
      top_project=
        names(pp)[1],
      top_project_fraction=
        as.numeric(pp[1]) / nrow(d),
      stringsAsFactors=FALSE
    )
  }
)

summary_df <- do.call(
  rbind,
  summary_list
)

stopifnot(
  all(
    summary_df$n_cells ==
      summary_df$freeze_n_cells
  )
)

cat("\n===== B/PLASMA PILOT SUMMARY =====\n")
print(summary_df, row.names=FALSE)

cat("\n===== PROJECT DISTRIBUTION BY GLOBAL CLUSTER =====\n")

for(k in bp_clusters){

  cat(
    "\n----------------------------------------\n"
  )
  cat("cluster ", k, "\n", sep="")

  d <- bp[
    bp$seurat_cluster == k,
    ,
    drop=FALSE
  ]

  print(
    sort(
      table(d$project_id),
      decreasing=TRUE
    )
  )
}

# ============================================================
# Write authoritative pilot B/plasma manifest
# ============================================================

cell_file <- file.path(
  out_dir,
  "b_plasma_pilot_cells_v1.tsv.gz"
)

con <- gzfile(cell_file, "wt")

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

write.table(
  summary_df,
  file=file.path(
    out_dir,
    "b_plasma_pilot_global_cluster_summary_v1.tsv"
  ),
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
      paste(bp_clusters, collapse=",")
    ),
    "source=large pilot frozen cluster annotation v2"
  ),
  file.path(
    out_dir,
    "B_PLASMA_PILOT_CELLSET_COMPLETE.ok"
  )
)

cat(
  "\nPASS: B/plasma pilot manifest = ",
  nrow(bp),
  " cells\n",
  sep=""
)
