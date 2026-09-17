#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

root <- normalizePath(
  args[[1]],
  mustWork=TRUE
)

tag <- args[[2]]

transfer_dir <- file.path(
  root,
  "pre_integration",
  "full_atlas_primary_v1",
  "annotation_transfer_by_library"
)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  paste0(
    "remaining_compartment_profile_v1__",
    tag
  )
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

targets <- c(
  "Neutrophil_granulocyte",
  "Platelet_megakaryocyte",
  "Erythroid",
  "Progenitor"
)

expected <- c(
  Neutrophil_granulocyte=53501L,
  Platelet_megakaryocyte=49585L,
  Erythroid=12177L,
  Progenitor=27775L
)

files <- list.files(
  transfer_dir,
  pattern="\\.tsv\\.gz$",
  full.names=TRUE
)

stopifnot(
  length(files) == 158L
)

z <- list()

for(i in seq_along(files)){

  d <- read.delim(
    files[[i]],
    sep="\t",
    quote="\"",
    comment.char="",
    stringsAsFactors=FALSE,
    check.names=FALSE
  )

  req <- c(
    "project_id",
    "library_key",
    "condition_binary",
    "integration_compartment_primary_v1"
  )

  stopifnot(
    all(req %in% names(d))
  )

  d <- d[
    d$integration_compartment_primary_v1 %in%
      targets,
    req,
    drop=FALSE
  ]

  if(nrow(d)){
    z[[length(z)+1L]] <- d
  }
}

x <- do.call(
  rbind,
  z
)

rm(z)

x$integration_compartment_primary_v1 <-
  as.character(
    x$integration_compartment_primary_v1
  )

x$project_id <-
  as.character(
    x$project_id
  )

x$library_key <-
  as.character(
    x$library_key
  )

x$condition_binary <-
  as.character(
    x$condition_binary
  )

# ============================================================
# Hard count check
# ============================================================

observed <- table(
  x$integration_compartment_primary_v1
)

for(k in names(expected)){

  if(
    !k %in% names(observed)
  ){
    stop(
      "Missing compartment: ",
      k
    )
  }

  if(
    as.integer(observed[[k]]) !=
      expected[[k]]
  ){
    stop(
      paste0(
        "Cell count mismatch for ",
        k,
        ": observed=",
        observed[[k]],
        " expected=",
        expected[[k]]
      )
    )
  }
}

# ============================================================
# Summary
# ============================================================

summary_rows <- list()
project_rows <- list()

for(k in targets){

  d <- x[
    x$integration_compartment_primary_v1 == k,
    ,
    drop=FALSE
  ]

  pt <- sort(
    table(d$project_id),
    decreasing=TRUE
  )

  p <- as.numeric(pt) /
    sum(pt)

  effective_n_projects <-
    1 / sum(p^2)

  max_project <-
    names(pt)[1]

  max_project_fraction <-
    as.numeric(pt[1]) /
    nrow(d)

  cond <- table(
    d$condition_binary
  )

  healthy_fraction <-
    if("healthy" %in% names(cond)){
      as.numeric(cond[["healthy"]]) /
        nrow(d)
    } else {
      NA_real_
    }

  disease_fraction <-
    if("disease" %in% names(cond)){
      as.numeric(cond[["disease"]]) /
        nrow(d)
    } else {
      NA_real_
    }

  summary_rows[[
    length(summary_rows)+1L
  ]] <- data.frame(
    compartment=k,
    n_cells=nrow(d),
    n_projects=length(
      unique(d$project_id)
    ),
    n_libraries=length(
      unique(d$library_key)
    ),
    effective_n_projects=
      effective_n_projects,
    max_project=max_project,
    max_project_fraction=
      max_project_fraction,
    healthy_fraction=
      healthy_fraction,
    disease_fraction=
      disease_fraction,
    stringsAsFactors=FALSE
  )

  for(pr in names(pt)){

    q <- d[
      d$project_id == pr,
      ,
      drop=FALSE
    ]

    project_rows[[
      length(project_rows)+1L
    ]] <- data.frame(
      compartment=k,
      project_id=pr,
      n_cells=nrow(q),
      cell_fraction=
        nrow(q) / nrow(d),
      n_libraries=
        length(
          unique(q$library_key)
        ),
      healthy_fraction=
        mean(
          q$condition_binary ==
            "healthy",
          na.rm=TRUE
        ),
      disease_fraction=
        mean(
          q$condition_binary ==
            "disease",
          na.rm=TRUE
        ),
      stringsAsFactors=FALSE
    )
  }
}

summary_df <- do.call(
  rbind,
  summary_rows
)

project_df <- do.call(
  rbind,
  project_rows
)

# ============================================================
# Neutrophil-specific decision flags
#
# These are descriptive, not an automatic integration decision.
# ============================================================

neu <- summary_df[
  summary_df$compartment ==
    "Neutrophil_granulocyte",
  ,
  drop=FALSE
]

neu$pilot_policy <-
  "cautious_diagnostic_only"

neu$full_reassessment_needed <-
  TRUE

neu$interpretation <- if(
  neu$max_project_fraction >= 0.70
){

  paste0(
    "strong_project_skew; ",
    "do_not_use_global_integration_as_primary ",
    "annotation evidence"
  )

} else if(
  neu$effective_n_projects < 3
){

  paste0(
    "low_effective_project_diversity; ",
    "prefer native/project-balanced annotation"
  )

} else {

  paste0(
    "full cohort diversity improved; ",
    "cautious integration benchmark may be evaluated, ",
    "but native RNA remains primary evidence"
  )
}

# ============================================================
# Write
# ============================================================

write.table(
  summary_df,
  file=file.path(
    out_dir,
    "remaining_compartment_summary_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  project_df,
  file=file.path(
    out_dir,
    "remaining_compartment_by_project_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  neu,
  file=file.path(
    out_dir,
    "neutrophil_full_design_reassessment_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

writeLines(
  c(
    "PASS",
    "remaining compartment profile v1",
    "n_Neutrophil_granulocyte=53501",
    "n_Platelet_megakaryocyte=49585",
    "n_Erythroid=12177",
    "n_Progenitor=27775",
    "source=full primary annotation transfer metadata",
    "no Cell Ranger output accessed"
  ),
  file.path(
    out_dir,
    "REMAINING_COMPARTMENT_PROFILE_COMPLETE.ok"
  )
)

cat(
  "\n===== REMAINING COMPARTMENTS =====\n"
)

print(
  summary_df,
  row.names=FALSE
)

cat(
  "\n===== NEUTROPHIL BY PROJECT =====\n"
)

print(
  project_df[
    project_df$compartment ==
      "Neutrophil_granulocyte",
    ,
    drop=FALSE
  ],
  row.names=FALSE
)

cat(
  "\n===== NEUTROPHIL DESIGN =====\n"
)

print(
  neu,
  row.names=FALSE
)

cat(
  "\nOUT_DIR=",
  out_dir,
  "\n",
  sep=""
)

cat(
  "\nPASS\n"
)
