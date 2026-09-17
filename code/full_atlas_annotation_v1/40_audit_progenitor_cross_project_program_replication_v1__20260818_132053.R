#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3L){
  stop("Usage: script <root> <tag> <cluster_dir>")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]
cluster_dir <- normalizePath(args[[3]], mustWork=TRUE)

suppressPackageStartupMessages({
  library(Matrix)
})

source(
  file.path(
    root,
    "full_atlas_primary_integration_v1",
    "00_common.R"
  )
)

DISCOVERY_PROJECT <- "GSE216007"
N_SUPPORT_EXPECTED <- 813L
N_SUPPORT_PROJECTS_EXPECTED <- 8L

# Identical definitions to Step39c.
programs <- list(

  strict_HSC=c(
    "HLF","AVP","MECOM","MLLT3",
    "HOPX","PROM1","CRHBP","PRDM16"
  ),

  lymphoid_priming=c(
    "FLT3","DNTT","MME","SELL",
    "IL7R","TCF7","CD7"
  ),

  B_commitment=c(
    "EBF1","PAX5","VPREB1","VPREB3",
    "IGLL1","CD79A","CD19","MS4A1"
  ),

  granulocytic=c(
    "MPO","ELANE","AZU1","PRTN3",
    "CTSG","CLEC12A","CEBPA","GFI1"
  ),

  megakaryocytic=c(
    "ITGA2B","GP1BB","GP9","PF4","PPBP",
    "VWF","NFE2","FLI1","TUBB1","TREML1"
  ),

  erythroid_nonglobin=c(
    "GATA1","KLF1","ALAS2","AHSP",
    "ANK1","SLC4A1","CA1","BPGM","TFRC"
  ),

  eos_baso_mast=c(
    "HDC","MS4A2","ENPP3","CLC","PRG2",
    "CPA3","TPSAB1","KIT","GATA2","CCR3"
  ),

  T_NK=c(
    "NKG7","GNLY","PRF1","GZMB",
    "CD3D","CD3E","TRAC","TRBC1"
  ),

  cycling=c(
    "MKI67","TOP2A","UBE2C","CENPF",
    "CDC20","CDC45","PCNA","TYMS"
  )
)

thresholds <- c(
  strict_HSC=3L,
  lymphoid_priming=3L,
  B_commitment=3L,
  granulocytic=3L,
  megakaryocytic=2L,
  erythroid_nonglobin=3L,
  eos_baso_mast=3L,
  T_NK=2L,
  cycling=3L
)

# Descriptive replication rule.
# A project must contain >=10 Progenitor cells to be assessable.
# Program support = >=3 matching cells AND >=5% of its Progenitor cells.
MIN_PROJECT_CELLS <- 10L
MIN_SUPPORT_CELLS <- 3L
MIN_SUPPORT_FRACTION <- 0.05

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

  if(system2(
    "gzip",
    c("-t", path),
    stdout=FALSE,
    stderr=FALSE
  ) != 0L){
    stop("gzip integrity failure: ", path)
  }
}

# ============================================================
# The 813 cells withheld from GSE216007 discovery
# ============================================================

support_file <- file.path(
  cluster_dir,
  "progenitor_non_discovery_projects_replication_support_v1.tsv"
)

support <- read_tsv(
  support_file
)

required <- c(
  "global_cell",
  "project_id",
  "library_key",
  "original_cell"
)

stopifnot(
  all(required %in% names(support))
)

for(nm in required){
  support[[nm]] <- as.character(
    support[[nm]]
  )
}

stopifnot(
  nrow(support) == N_SUPPORT_EXPECTED,
  !anyDuplicated(support$global_cell),
  !any(support$project_id == DISCOVERY_PROJECT),
  length(unique(support$project_id)) ==
    N_SUPPORT_PROJECTS_EXPECTED
)

cat(
  "Replication-support cells=",
  nrow(support),
  "\nProjects=",
  length(unique(support$project_id)),
  "\n",
  sep=""
)

# ============================================================
# Resolved native-RNA sources
# ============================================================

lib <- read_tsv(
  file.path(
    root,
    "pre_integration",
    "pilot_unintegrated_large_v1",
    "resolved_library_table.tsv"
  )
)

stopifnot(
  all(
    c(
      "project_id",
      "library_key",
      "final_rds_resolved"
    ) %in% names(lib)
  )
)

lib$project_id <- as.character(lib$project_id)
lib$library_key <- as.character(lib$library_key)

lib$key <- paste(
  lib$project_id,
  lib$library_key,
  sep="|||"
)

support$key <- paste(
  support$project_id,
  support$library_key,
  sep="|||"
)

all_panel_genes <- unique(
  unlist(programs)
)

cell_rows <- list()

keys <- sort(
  unique(support$key)
)

# ============================================================
# Cell-level native program detection
# ============================================================

for(ii in seq_along(keys)){

  key <- keys[[ii]]

  d <- support[
    support$key == key,
    ,
    drop=FALSE
  ]

  li <- match(
    key,
    lib$key
  )

  if(is.na(li)){
    stop("Library lookup failed: ", key)
  }

  source_rds <- as.character(
    lib$final_rds_resolved[[li]]
  )

  if(!file.exists(source_rds)){
    stop("Missing source RDS: ", source_rds)
  }

  obj <- readRDS(source_rds)
  counts0 <- get_rna_counts(obj)

  stopifnot(
    all(
      d$original_cell %in%
        colnames(counts0)
    )
  )

  full_x <- counts0[
    ,
    d$original_cell,
    drop=FALSE
  ]

  panel_present <- intersect(
    all_panel_genes,
    rownames(counts0)
  )

  x <- counts0[
    panel_present,
    d$original_cell,
    drop=FALSE
  ]

  out <- data.frame(
    global_cell=d$global_cell,
    project_id=d$project_id,
    library_key=d$library_key,
    nCount=as.numeric(
      Matrix::colSums(full_x)
    ),
    nFeature=as.numeric(
      Matrix::colSums(full_x > 0)
    ),
    stringsAsFactors=FALSE
  )

  for(program in names(programs)){

    genes <- intersect(
      programs[[program]],
      rownames(x)
    )

    if(length(genes)){

      n_detected <- as.integer(
        Matrix::colSums(
          x[
            genes,
            ,
            drop=FALSE
          ] > 0
        )
      )

    } else {

      n_detected <- rep(
        0L,
        ncol(x)
      )
    }

    out[[paste0(program, "_n")]] <-
      n_detected

    out[[paste0(program, "_flag")]] <-
      n_detected >=
        thresholds[[program]]
  }

  # ----------------------------------------------------------
  # Candidate biological combinations.
  # Descriptive only: these are NOT canonical annotations.
  # ----------------------------------------------------------

  out$primitive_HSC_low_commitment_flag <-
    out$strict_HSC_flag &
    !out$B_commitment_flag &
    !out$granulocytic_flag &
    !out$eos_baso_mast_flag &
    !out$T_NK_flag

  out$lymphoid_primed_nonB_flag <-
    out$lymphoid_priming_flag &
    !out$B_commitment_flag &
    !out$T_NK_flag

  out$MegE_flag <-
    out$megakaryocytic_flag &
    out$erythroid_nonglobin_flag

  out$meg_only_flag <-
    out$megakaryocytic_flag &
    !out$erythroid_nonglobin_flag

  out$HSC_granulocytic_flag <-
    out$strict_HSC_flag &
    out$granulocytic_flag

  out$B_and_lymphoid_flag <-
    out$B_commitment_flag &
    out$lymphoid_priming_flag

  out$B_and_TNK_flag <-
    out$B_commitment_flag &
    out$T_NK_flag

  out$eos_and_meg_flag <-
    out$eos_baso_mast_flag &
    out$megakaryocytic_flag

  cell_rows[[
    length(cell_rows)+1L
  ]] <- out

  cat(
    sprintf(
      "READ %3d/%3d %s cells=%d\n",
      ii,
      length(keys),
      key,
      nrow(out)
    )
  )

  rm(
    obj,
    counts0,
    full_x,
    x,
    out
  )

  invisible(gc())
}

cell_df <- do.call(
  rbind,
  cell_rows
)

rownames(cell_df) <- NULL

stopifnot(
  nrow(cell_df) ==
    N_SUPPORT_EXPECTED,
  !anyDuplicated(
    cell_df$global_cell
  )
)

# ============================================================
# Long-form summaries
# ============================================================

flag_cols <- grep(
  "_flag$",
  names(cell_df),
  value=TRUE
)

signal_names <- sub(
  "_flag$",
  "",
  flag_cols
)

summarize_group <- function(d, group_name){

  do.call(
    rbind,
    lapply(
      seq_along(flag_cols),
      function(i){

        flag <- flag_cols[[i]]

        data.frame(
          group=group_name,
          signal=signal_names[[i]],
          n_cells=nrow(d),
          n_flagged=sum(
            d[[flag]],
            na.rm=TRUE
          ),
          fraction_flagged=mean(
            d[[flag]],
            na.rm=TRUE
          ),
          stringsAsFactors=FALSE
        )
      }
    )
  )
}

project_split <- split(
  cell_df,
  cell_df$project_id
)

project_summary <- do.call(
  rbind,
  lapply(
    names(project_split),
    function(p){

      z <- summarize_group(
        project_split[[p]],
        p
      )

      names(z)[1] <- "project_id"
      z
    }
  )
)

rownames(project_summary) <- NULL

library_split <- split(
  cell_df,
  paste(
    cell_df$project_id,
    cell_df$library_key,
    sep="|||"
  )
)

library_summary <- do.call(
  rbind,
  lapply(
    names(library_split),
    function(k){

      d <- library_split[[k]]

      z <- summarize_group(
        d,
        k
      )

      names(z)[1] <- "project_library"

      z$project_id <-
        d$project_id[[1]]

      z$library_key <-
        d$library_key[[1]]

      z[
        ,
        c(
          "project_id",
          "library_key",
          "signal",
          "n_cells",
          "n_flagged",
          "fraction_flagged"
        ),
        drop=FALSE
      ]
    }
  )
)

rownames(library_summary) <- NULL

# ============================================================
# Project replication assessment
# ============================================================

replication_summary <- do.call(
  rbind,
  lapply(
    signal_names,
    function(signal){

      d <- project_summary[
        project_summary$signal == signal,
        ,
        drop=FALSE
      ]

      eligible <-
        d$n_cells >=
          MIN_PROJECT_CELLS

      supported <-
        eligible &
        d$n_flagged >=
          MIN_SUPPORT_CELLS &
        d$fraction_flagged >=
          MIN_SUPPORT_FRACTION

      flag_col <- paste0(
        signal,
        "_flag"
      )

      n_total_flagged <- sum(
        cell_df[[flag_col]]
      )

      data.frame(
        signal=signal,

        n_projects_total=nrow(d),

        n_projects_eligible=
          sum(eligible),

        n_projects_with_any=
          sum(d$n_flagged > 0),

        n_projects_supported=
          sum(supported),

        supported_projects=
          paste(
            d$project_id[supported],
            collapse=","
          ),

        total_flagged_cells=
          n_total_flagged,

        total_fraction=
          n_total_flagged /
          nrow(cell_df),

        replicated_ge2_projects=
          sum(supported) >= 2L,

        stringsAsFactors=FALSE
      )
    }
  )
)

# ============================================================
# Project cell counts
# ============================================================

project_counts <- as.data.frame(
  table(
    project_id=cell_df$project_id
  ),
  stringsAsFactors=FALSE
)

names(project_counts)[2] <- "n_cells"

project_counts$n_cells <- as.integer(
  project_counts$n_cells
)

# ============================================================
# Outputs
# ============================================================

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "progenitor_cross_project_program_replication_v1__",
    tag
  )
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

cell_file <- file.path(
  out_dir,
  paste0(
    "progenitor_replication_support_cell_programs_v1__",
    tag,
    ".tsv.gz"
  )
)

write_gz_tsv(
  cell_df,
  cell_file
)

write.table(
  project_counts,
  file.path(
    out_dir,
    "progenitor_replication_project_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  project_summary,
  file.path(
    out_dir,
    "progenitor_replication_project_program_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  library_summary,
  file.path(
    out_dir,
    "progenitor_replication_library_program_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  replication_summary,
  file.path(
    out_dir,
    "progenitor_program_cross_project_replication_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Re-read / SHA
# ============================================================

chk <- read_tsv(
  cell_file
)

stopifnot(
  nrow(chk) ==
    N_SUPPORT_EXPECTED,
  identical(
    as.character(chk$global_cell),
    as.character(cell_df$global_cell)
  )
)

sha_files <- c(
  cell_file,
  file.path(
    out_dir,
    "progenitor_replication_project_counts_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_replication_project_program_summary_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_replication_library_program_summary_v1.tsv"
  ),
  file.path(
    out_dir,
    "progenitor_program_cross_project_replication_summary_v1.tsv"
  )
)

sha_file <- file.path(
  out_dir,
  paste0(
    "SHA256SUMS_v1__",
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
  stop("sha256 creation failed")
}

status <- system2(
  "sha256sum",
  c("-c", sha_file),
  stdout=FALSE,
  stderr=FALSE
)

if(status != 0L){
  stop("sha256 verification failed")
}

# Completion LAST.
writeLines(
  c(
    "PASS",
    "Progenitor cross-project program replication audit v1",
    "discovery_project excluded=GSE216007",
    "n_replication_support_cells=813",
    "n_replication_support_projects=8",
    "",
    "program definitions identical to Step39c",
    "native SoupX-corrected RNA",
    "no cross-project integration",
    "no cluster-label transfer",
    "",
    "project assessable threshold=10 cells",
    "project program support threshold=3 cells and 5 percent",
    "replicated criterion=support in at least 2 projects",
    "",
    "descriptive replication criterion only",
    "no automatic annotation",
    "condition not used",
    "cellranger_count not accessed",
    "",
    "gzip_integrity=PASS",
    "re_read_count_order=PASS",
    "sha256_manifest=PASS"
  ),
  file.path(
    out_dir,
    "PROGENITOR_CROSS_PROJECT_PROGRAM_REPLICATION_COMPLETE.ok"
  )
)

cat(
  "\n===== PROJECT COUNTS =====\n"
)

print(
  project_counts,
  row.names=FALSE
)

cat(
  "\n===== REPLICATION SUMMARY =====\n"
)

print(
  replication_summary,
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS: Progenitor cross-project program replication completed\n"
)
