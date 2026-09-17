#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if(length(args) < 3L){
  stop("Usage: script <root> <tag> <calibration_dir>")
}

root <- normalizePath(args[[1]], mustWork=TRUE)
tag <- args[[2]]
cal_dir <- normalizePath(args[[3]], mustWork=TRUE)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "platelet_transfer_precision_audit_v1__",
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
  "PLATELET_TRANSFER_PRECISION_AUDIT_COMPLETE.ok"
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

pred_file <- file.path(
  cal_dir,
  "platelet_transfer_calibration_predictions_v1.tsv.gz"
)

if(!file.exists(pred_file)){
  stop("Missing: ", pred_file)
}

p <- read_tsv(pred_file)

required <- c(
  "predicted.id",
  "prediction.score.max",
  "actual_class",
  "project_id",
  "actual_core",
  "predicted_core",
  "actual_axis",
  "predicted_axis",
  "actual_state",
  "predicted_state"
)

stopifnot(
  all(required %in% names(p))
)

p$predicted.id <- as.character(p$predicted.id)
p$actual_class <- as.character(p$actual_class)

p$actual_core <- as.character(p$actual_core)
p$predicted_core <- as.character(p$predicted_core)

p$actual_axis <- as.character(p$actual_axis)
p$predicted_axis <- as.character(p$predicted_axis)

p$actual_state <- as.character(p$actual_state)
p$predicted_state <- as.character(p$predicted_state)

p$project_id <- as.character(p$project_id)

p$prediction.score.max <- as.numeric(
  p$prediction.score.max
)

stopifnot(
  nrow(p) == 5312L,
  all(is.finite(p$prediction.score.max))
)

p$exact_correct <-
  p$predicted.id == p$actual_class

p$core_correct <-
  p$predicted_core == p$actual_core

p$axis_correct <-
  p$predicted_axis == p$actual_axis

p$state_correct <-
  p$predicted_state == p$actual_state

# ============================================================
# Threshold grid
# ============================================================

thresholds <- c(
  0.00,
  0.60,
  0.70,
  0.75,
  0.80,
  0.85,
  0.90,
  0.95
)

# ============================================================
# Overall conditional accuracy / coverage
# ============================================================

overall <- do.call(
  rbind,
  lapply(
    thresholds,
    function(thr){

      d <- p[
        p$prediction.score.max >= thr,
        ,
        drop=FALSE
      ]

      data.frame(
        threshold=thr,
        n_accepted=nrow(d),
        coverage=nrow(d) / nrow(p),

        exact_accuracy=
          if(nrow(d))
            mean(d$exact_correct)
          else
            NA_real_,

        core_accuracy=
          if(nrow(d))
            mean(d$core_correct)
          else
            NA_real_,

        axis_accuracy=
          if(nrow(d))
            mean(d$axis_correct)
          else
            NA_real_,

        state_accuracy=
          if(nrow(d))
            mean(d$state_correct)
          else
            NA_real_,

        stringsAsFactors=FALSE
      )
    }
  )
)

# ============================================================
# Predicted-class precision by threshold
#
# This is what we need before full transfer.
# ============================================================

classes <- sort(
  unique(
    p$predicted.id
  )
)

class_precision <- list()

for(cl in classes){

  actual_total <- sum(
    p$actual_class == cl
  )

  for(thr in thresholds){

    d <- p[
      p$predicted.id == cl &
        p$prediction.score.max >= thr,
      ,
      drop=FALSE
    ]

    n_correct <- if(nrow(d)){
      sum(d$actual_class == cl)
    } else {
      0L
    }

    class_precision[[
      length(class_precision)+1L
    ]] <- data.frame(
      predicted_class=cl,
      threshold=thr,

      n_predicted_accepted=
        nrow(d),

      exact_precision=
        if(nrow(d))
          mean(d$actual_class == cl)
        else
          NA_real_,

      correct_accepted=
        n_correct,

      actual_class_total=
        actual_total,

      accepted_recall=
        if(actual_total > 0L)
          n_correct / actual_total
        else
          NA_real_,

      median_score=
        if(nrow(d))
          median(d$prediction.score.max)
        else
          NA_real_,

      stringsAsFactors=FALSE
    )
  }
}

class_precision <- do.call(
  rbind,
  class_precision
)

# ============================================================
# Core precision by predicted core
# ============================================================

core_values <- sort(
  unique(
    p$predicted_core
  )
)

core_precision <- list()

for(value in core_values){

  for(thr in thresholds){

    d <- p[
      p$predicted_core == value &
        p$prediction.score.max >= thr,
      ,
      drop=FALSE
    ]

    core_precision[[
      length(core_precision)+1L
    ]] <- data.frame(
      predicted_core=value,
      threshold=thr,
      n_predicted_accepted=nrow(d),

      precision=
        if(nrow(d))
          mean(d$actual_core == value)
        else
          NA_real_,

      stringsAsFactors=FALSE
    )
  }
}

core_precision <- do.call(
  rbind,
  core_precision
)

# ============================================================
# Transcriptional-axis precision
# ============================================================

axis_values <- sort(
  unique(
    p$predicted_axis
  )
)

axis_precision <- list()

for(value in axis_values){

  for(thr in thresholds){

    d <- p[
      p$predicted_axis == value &
        p$prediction.score.max >= thr,
      ,
      drop=FALSE
    ]

    axis_precision[[
      length(axis_precision)+1L
    ]] <- data.frame(
      predicted_axis=value,
      threshold=thr,
      n_predicted_accepted=nrow(d),

      precision=
        if(nrow(d))
          mean(d$actual_axis == value)
        else
          NA_real_,

      stringsAsFactors=FALSE
    )
  }
}

axis_precision <- do.call(
  rbind,
  axis_precision
)

# ============================================================
# State precision
# ============================================================

state_values <- sort(
  unique(
    p$predicted_state
  )
)

state_precision <- list()

for(value in state_values){

  for(thr in thresholds){

    d <- p[
      p$predicted_state == value &
        p$prediction.score.max >= thr,
      ,
      drop=FALSE
    ]

    state_precision[[
      length(state_precision)+1L
    ]] <- data.frame(
      predicted_state=value,
      threshold=thr,
      n_predicted_accepted=nrow(d),

      precision=
        if(nrow(d))
          mean(d$actual_state == value)
        else
          NA_real_,

      stringsAsFactors=FALSE
    )
  }
}

state_precision <- do.call(
  rbind,
  state_precision
)

# ============================================================
# High-value compact table at 0.70 / 0.80 / 0.90
# ============================================================

compact_thresholds <- c(
  0.70,
  0.80,
  0.90
)

compact <- class_precision[
  class_precision$threshold %in%
    compact_thresholds,
  ,
  drop=FALSE
]

compact <- compact[
  order(
    compact$predicted_class,
    compact$threshold
  ),
  ,
  drop=FALSE
]

# ============================================================
# Score distribution for correct vs incorrect predictions
# ============================================================

score_distribution <- do.call(
  rbind,
  lapply(
    sort(
      unique(
        p$actual_class
      )
    ),
    function(cl){

      d <- p[
        p$actual_class == cl,
        ,
        drop=FALSE
      ]

      data.frame(
        actual_class=cl,
        n_cells=nrow(d),

        correct_fraction=
          mean(d$exact_correct),

        median_score_correct=
          if(any(d$exact_correct))
            median(
              d$prediction.score.max[
                d$exact_correct
              ]
            )
          else
            NA_real_,

        median_score_incorrect=
          if(any(!d$exact_correct))
            median(
              d$prediction.score.max[
                !d$exact_correct
              ]
            )
          else
            NA_real_,

        q25_score_correct=
          if(any(d$exact_correct))
            unname(
              quantile(
                d$prediction.score.max[
                  d$exact_correct
                ],
                0.25
              )
            )
          else
            NA_real_,

        q75_score_incorrect=
          if(any(!d$exact_correct))
            unname(
              quantile(
                d$prediction.score.max[
                  !d$exact_correct
                ],
                0.75
              )
            )
          else
            NA_real_,

        stringsAsFactors=FALSE
      )
    }
  )
)

# ============================================================
# Write
# ============================================================

write.table(
  overall,
  file=file.path(
    out_dir,
    "platelet_transfer_precision_overall_by_threshold_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  class_precision,
  file=file.path(
    out_dir,
    "platelet_transfer_class_precision_by_threshold_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  core_precision,
  file=file.path(
    out_dir,
    "platelet_transfer_core_precision_by_threshold_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  axis_precision,
  file=file.path(
    out_dir,
    "platelet_transfer_axis_precision_by_threshold_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  state_precision,
  file=file.path(
    out_dir,
    "platelet_transfer_state_precision_by_threshold_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  compact,
  file=file.path(
    out_dir,
    "platelet_transfer_class_precision_compact_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  score_distribution,
  file=file.path(
    out_dir,
    "platelet_transfer_score_correct_incorrect_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

writeLines(
  c(
    "PASS",
    "Platelet transfer precision audit v1",
    "source=step30 held-out validation predictions only",
    "n_validation=5312",
    "",
    "no RNA reread",
    "no remapping",
    "no Cell Ranger access",
    "",
    "purpose=freeze full-transfer acceptance thresholds",
    "class precision is evaluated as a function of prediction score"
  ),
  done_file
)

cat(
  "\n===== OVERALL BY THRESHOLD =====\n"
)

print(
  overall,
  row.names=FALSE
)

cat(
  "\n===== CLASS PRECISION: 0.70 / 0.80 / 0.90 =====\n"
)

print(
  compact,
  row.names=FALSE
)

cat(
  "\n===== CORE PRECISION =====\n"
)

print(
  core_precision[
    core_precision$threshold %in%
      compact_thresholds,
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat(
  "\n===== AXIS PRECISION =====\n"
)

print(
  axis_precision[
    axis_precision$threshold %in%
      compact_thresholds,
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat(
  "\n===== STATE PRECISION =====\n"
)

print(
  state_precision[
    state_precision$threshold %in%
      compact_thresholds,
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat(
  "\n===== CORRECT vs INCORRECT SCORE =====\n"
)

print(
  score_distribution,
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat("\nPASS\n")
