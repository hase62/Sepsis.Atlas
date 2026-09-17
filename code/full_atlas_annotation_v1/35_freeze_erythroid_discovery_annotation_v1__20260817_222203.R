#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 4L){
  stop("Usage: script <root> <tag> <discovery_dir> <decision_file>")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]
discovery_dir <- normalizePath(args[[3]], mustWork=TRUE)
decision_file <- normalizePath(args[[4]], mustWork=TRUE)

N_DISCOVERY <- 8692L

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0("erythroid_discovery_annotation_freeze_v1__", tag)
)

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

done_file <- file.path(
  out_dir,
  "ERYTHROID_DISCOVERY_ANNOTATION_FREEZE_COMPLETE.ok"
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

write_gz_tsv <- function(x, path){

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

  if(system2("gzip", c("-t", path)) != 0L){
    stop("gzip test failed: ", path)
  }
}

assignment_file <- file.path(
  discovery_dir,
  "erythroid_balanced_discovery_cluster_assignments_v1.tsv.gz"
)

a <- read_tsv(assignment_file)
decision <- read_tsv(decision_file)

stopifnot(
  nrow(a) == N_DISCOVERY,
  !anyDuplicated(a$global_cell),
  nrow(decision) == 10L,
  setequal(as.character(decision$cluster), as.character(0:9))
)

a$cluster_r0p4 <- as.character(a$cluster_r0p4)
decision$cluster <- as.character(decision$cluster)

mi <- match(
  a$cluster_r0p4,
  decision$cluster
)

stopifnot(!anyNA(mi))

cols <- c(
  "erythroid_core_v1",
  "erythroid_maturation_axis_v1",
  "erythroid_state_v1",
  "erythroid_evidence_confidence_v1",
  "erythroid_evidence_flag_v1",
  "rationale"
)

for(nm in cols){
  a[[nm]] <- decision[[nm]][mi]
}

# Deferred populations must not carry maturation/state claims.
deferred <- grepl("^Deferred_", a$erythroid_core_v1)

stopifnot(
  all(a$erythroid_maturation_axis_v1[deferred] == "not_applicable"),
  all(a$erythroid_state_v1[deferred] == "not_applicable")
)

stopifnot(
  all(a$erythroid_maturation_axis_v1[!deferred] != "not_applicable"),
  all(a$erythroid_state_v1[!deferred] == "none")
)

a$erythroid_transfer_class_v1 <- paste(
  a$erythroid_core_v1,
  a$erythroid_maturation_axis_v1,
  a$erythroid_state_v1,
  sep="|||"
)

count1 <- function(x, name){

  z <- as.data.frame(
    table(x),
    stringsAsFactors=FALSE
  )

  names(z) <- c(name, "n_cells")
  z
}

core_counts <- count1(
  a$erythroid_core_v1,
  "erythroid_core_v1"
)

maturation_counts <- count1(
  a$erythroid_maturation_axis_v1,
  "erythroid_maturation_axis_v1"
)

state_counts <- count1(
  a$erythroid_state_v1,
  "erythroid_state_v1"
)

transfer_counts <- count1(
  a$erythroid_transfer_class_v1,
  "erythroid_transfer_class_v1"
)

stopifnot(
  sum(core_counts$n_cells) == N_DISCOVERY,
  sum(maturation_counts$n_cells) == N_DISCOVERY,
  sum(state_counts$n_cells) == N_DISCOVERY,
  sum(transfer_counts$n_cells) == N_DISCOVERY
)

# Expected discovery semantics.
stopifnot(
  sum(a$erythroid_core_v1 == "Erythroid_like") == 8094L,
  sum(a$erythroid_core_v1 == "Deferred_non_erythroid_myeloid_like") == 352L,
  sum(a$erythroid_core_v1 == "Deferred_mixed_lineage_like") == 246L,

  sum(a$erythroid_maturation_axis_v1 ==
        "late_maturation_globin_dominant_candidate") == 3744L,

  sum(a$erythroid_maturation_axis_v1 ==
        "terminal_hemoglobinization_high") == 2491L,

  sum(a$erythroid_maturation_axis_v1 ==
        "regulatory_rich_erythroid") == 1081L,

  sum(a$erythroid_maturation_axis_v1 ==
        "unresolved") == 778L,

  sum(a$erythroid_maturation_axis_v1 ==
        "not_applicable") == 598L
)

annotation_file <- file.path(
  out_dir,
  paste0(
    "erythroid_discovery_cell_annotation_v1__",
    tag,
    ".tsv.gz"
  )
)

write_gz_tsv(a, annotation_file)

decision_out <- file.path(
  out_dir,
  "erythroid_r0p4_annotation_decision_v1.tsv"
)

write.table(
  decision,
  decision_out,
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  core_counts,
  file.path(out_dir, "erythroid_discovery_core_counts_v1.tsv"),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  maturation_counts,
  file.path(out_dir, "erythroid_discovery_maturation_counts_v1.tsv"),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  state_counts,
  file.path(out_dir, "erythroid_discovery_state_counts_v1.tsv"),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  transfer_counts,
  file.path(out_dir, "erythroid_discovery_transfer_class_counts_v1.tsv"),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

freeze <- list(
  annotation_version="Erythroid_discovery_annotation_v1",
  backbone_resolution="r0.4",
  n_cells=N_DISCOVERY,
  annotation=a,
  decision=decision,
  source_discovery_dir=discovery_dir,
  source_assignment_file=assignment_file,
  source_decision_file=decision_file
)

rds_file <- file.path(
  out_dir,
  paste0(
    "erythroid_discovery_annotation_freeze_v1__",
    tag,
    ".rds"
  )
)

tmp_rds <- paste0(rds_file, ".tmp")

saveRDS(
  freeze,
  tmp_rds,
  compress=TRUE
)

if(!file.rename(tmp_rds, rds_file)){
  stop("Atomic RDS rename failed")
}

chk_rds <- readRDS(rds_file)

stopifnot(
  chk_rds$n_cells == N_DISCOVERY,
  nrow(chk_rds$annotation) == N_DISCOVERY,
  identical(
    as.character(chk_rds$annotation$global_cell),
    as.character(a$global_cell)
  )
)

chk_tsv <- read_tsv(annotation_file)

stopifnot(
  nrow(chk_tsv) == N_DISCOVERY,
  identical(
    as.character(chk_tsv$global_cell),
    as.character(a$global_cell)
  ),
  identical(
    as.character(chk_tsv$erythroid_core_v1),
    as.character(a$erythroid_core_v1)
  ),
  identical(
    as.character(chk_tsv$erythroid_maturation_axis_v1),
    as.character(a$erythroid_maturation_axis_v1)
  )
)

sha_files <- c(
  annotation_file,
  decision_out,
  file.path(out_dir, "erythroid_discovery_core_counts_v1.tsv"),
  file.path(out_dir, "erythroid_discovery_maturation_counts_v1.tsv"),
  file.path(out_dir, "erythroid_discovery_state_counts_v1.tsv"),
  file.path(out_dir, "erythroid_discovery_transfer_class_counts_v1.tsv"),
  rds_file
)

sha_out <- file.path(
  out_dir,
  paste0("SHA256SUMS_v1__", tag, ".txt")
)

status <- system2(
  "sha256sum",
  sha_files,
  stdout=sha_out
)

if(status != 0L){
  stop("sha256sum failed")
}

writeLines(
  c(
    "PASS",
    "Erythroid discovery annotation freeze v1",
    "annotation_version=Erythroid_discovery_annotation_v1",
    "n_cells=8692",
    "backbone_resolution=r0.4",
    "",
    "Erythroid_like=8094",
    "Deferred_non_erythroid_myeloid_like=352",
    "Deferred_mixed_lineage_like=246",
    "",
    "maturation axis is program-based; canonical BasoE/PolyE/OrthoE stages not asserted",
    "c5 retained as Erythroid_like with unresolved maturation",
    "c5 evidence flag=low_level_nonerythroid_signal",
    "core / maturation / state / evidence flag remain separate",
    "",
    "gzip_integrity=PASS",
    "re_read_check=PASS",
    "freeze_rds_re_read=PASS",
    "sha256_manifest=PASS",
    "cellranger_count not accessed"
  ),
  done_file
)

cat("\n===== CORE =====\n")
print(core_counts, row.names=FALSE)

cat("\n===== MATURATION =====\n")
print(maturation_counts, row.names=FALSE)

cat("\n===== STATE =====\n")
print(state_counts, row.names=FALSE)

cat("\n===== TRANSFER CLASSES =====\n")
print(transfer_counts, row.names=FALSE)

cat("\nOUT_DIR=", out_dir, "\n", sep="")
cat("\nPASS: Erythroid discovery annotation freeze completed\n")
