#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

base <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1"
)

cl_file <- file.path(
  base,
  "integrated_compartment_clustering",
  "T_NK_combined",
  "integrated_clustering_v1.rds"
)

diag_file <- file.path(
  base,
  "integrated_compartment_clustering",
  "T_NK_combined",
  "cluster_diagnostics_v1.tsv"
)

marker_file <- file.path(
  base,
  "tnk_native_marker_aggregation_v1",
  "top_native_markers_res_0p4_v1.tsv"
)

summary_file <- file.path(
  base,
  "tnk_native_marker_aggregation_v1",
  "cluster_marker_summary_res_0p4_v1.tsv"
)

out_dir <- file.path(
  base,
  "tnk_annotation_decision_v1"
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

for (f in c(cl_file, diag_file, marker_file, summary_file)) {
  if (!file.exists(f)) {
    stop("Missing input: ", f)
  }
}

cat("Reading clustering RDS...\n")
x <- readRDS(cl_file)

stopifnot(
  "cluster_res_0p4" %in% names(x$clusters),
  "cluster_res_0p6" %in% names(x$clusters)
)

c04 <- as.character(x$clusters$cluster_res_0p4)
c06 <- as.character(x$clusters$cluster_res_0p6)

stopifnot(
  length(c04) == length(c06),
  length(c04) == length(x$cells)
)

cat("Building r0.4 -> r0.6 mapping...\n")

tab <- table(
  parent_r0p4=c04,
  child_r0p6=c06
)

mapping <- do.call(
  rbind,
  lapply(
    rownames(tab),
    function(parent) {

      z <- tab[parent,]
      z <- sort(
        z[z > 0],
        decreasing=TRUE
      )

      data.frame(
        cluster=parent,
        child_clusters=paste(
          names(z),
          collapse=","
        ),
        child_cells=paste(
          as.integer(z),
          collapse=","
        ),
        child_fractions=paste(
          sprintf(
            "%.3f",
            as.integer(z) / sum(z)
          ),
          collapse=","
        ),
        n_children=length(z),
        largest_child_fraction=
          max(z) / sum(z),
        stringsAsFactors=FALSE
      )
    }
  )
)

cat("Reading diagnostics...\n")

diag <- read.delim(
  diag_file,
  sep="\t",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

diag <- diag[
  diag$resolution == 0.4,
  ,
  drop=FALSE
]

diag$cluster <- as.character(
  diag$cluster
)

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

cat("Reading native markers...\n")

mk <- read.delim(
  marker_file,
  sep="\t",
  quote="\"",
  comment.char="",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

mk$cluster <- as.character(
  mk$cluster
)

marker_rows <- lapply(
  sort(unique(sm$cluster)),
  function(k) {

    d <- mk[
      mk$cluster == k &
      !mk$technical_gene &
      mk$marker_candidate,
      ,
      drop=FALSE
    ]

    d <- head(d, 15)

    data.frame(
      cluster=k,
      top_markers=
        if (nrow(d))
          paste(d$gene, collapse=",")
        else
          NA_character_,
      top_marker_log2FC=
        if (nrow(d))
          paste(
            sprintf(
              "%.2f",
              d$log2FC_cluster_vs_rest
            ),
            collapse=","
          )
        else
          NA_character_,
      top_marker_delta_pct=
        if (nrow(d))
          paste(
            sprintf(
              "%.2f",
              d$delta_pct
            ),
            collapse=","
          )
        else
          NA_character_,
      stringsAsFactors=FALSE
    )
  }
)

marker_rows <- do.call(
  rbind,
  marker_rows
)

keep_diag <- c(
  "cluster",
  "n_cells",
  "n_projects",
  "max_project_fraction",
  "project_entropy_normalized",
  "top_transfer_celltype",
  "transfer_celltype_purity",
  "healthy_fraction",
  "disease_fraction"
)

out <- merge(
  sm,
  diag[, keep_diag, drop=FALSE],
  by="cluster",
  all.x=TRUE,
  sort=FALSE
)

out <- merge(
  out,
  mapping,
  by="cluster",
  all.x=TRUE,
  sort=FALSE
)

out <- merge(
  out,
  marker_rows,
  by="cluster",
  all.x=TRUE,
  sort=FALSE
)

out$cluster_num <- suppressWarnings(
  as.integer(out$cluster)
)

out <- out[
  order(out$cluster_num),
  ,
  drop=FALSE
]

out$cluster_num <- NULL

out_file <- file.path(
  out_dir,
  "tnk_r0p4_annotation_decision_sheet_v1.tsv"
)

write.table(
  out,
  file=out_file,
  sep="\t",
  quote=TRUE,
  qmethod="double",
  row.names=FALSE
)

cat("\n===== VALIDATION =====\n")
cat("rows              =", nrow(out), "\n")
cat("unique clusters   =", length(unique(out$cluster)), "\n")
cat("expected clusters =", length(unique(c04)), "\n")
cat("output             =", out_file, "\n")

stopifnot(
  nrow(out) == length(unique(c04)),
  length(unique(out$cluster)) ==
    length(unique(c04)),
  file.exists(out_file)
)

cat("\nPASS: T/NK annotation decision sheet created\n")
