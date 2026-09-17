#!/usr/bin/env Rscript

root <- normalizePath(".", mustWork=TRUE)

in_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_full_primary_transfer_v1__repaired_20260814_1430"
)

ann_file <- file.path(
  in_dir,
  "b_plasma_full_primary_annotation_repaired_v1.rds"
)

out_dir <- file.path(
  root,
  "atlas",
  "full_atlas_annotation_v1",
  "b_plasma_full_annotation_freeze_v1__20260814"
)

dir.create(
  out_dir,
  recursive=TRUE,
  showWarnings=FALSE
)

done_file <- file.path(
  out_dir,
  "B_PLASMA_FULL_ANNOTATION_FREEZE_COMPLETE.ok"
)

if(file.exists(done_file)){
  stop(
    "Freeze already completed; refusing to overwrite: ",
    done_file
  )
}

stopifnot(
  file.exists(ann_file)
)

ann <- readRDS(ann_file)

stopifnot(
  nrow(ann) == 45955L,
  !anyDuplicated(ann$cell_id)
)

required <- c(
  "b_plasma_final_taxonomy_id_v1",
  "b_plasma_transfer_confidence_v1",
  "b_plasma_transfer_accepted_v1",
  "b_plasma_prediction_score_v1",
  "b_plasma_prediction_margin_v1",
  "b_plasma_annotation_source_v1"
)

stopifnot(
  all(required %in% names(ann))
)

# ============================================================
# Final biological interpretation after native-RNA audits
# ============================================================

taxonomy <- data.frame(

  taxonomy_id=c(
    "T01","T02","T03","T04",
    "T05","T06","T07","T08",
    "T09","T10","T11","T12"
  ),

  final_core_identity=c(
    "Naive_B_like",
    "Memory_B_like",
    "Class_switched_memory_B_like",
    "B_cell_unresolved",
    "Naive_B_like",
    "Naive_B_like",
    "Naive_B_like",
    "Naive_B_like",
    "Naive_B_like",
    "Deferred_non_B_platelet_like",
    "Class_switched_memory_B_like",
    "Deferred_non_B_T_cell_like"
  ),

  final_state=c(
    "none",
    "none",
    "none",
    "activated_atypical_B_candidate",
    "interferon_stimulated",
    "early_activation_candidate",
    "secretory_differentiation_candidate",
    "activation_stress",
    "unresolved_state",
    "not_applicable",
    "interferon_stimulated_candidate",
    "not_applicable"
  ),

  final_annotation_scope=c(
    "B",
    "B",
    "B",
    "B_unresolved",
    "B",
    "B",
    "B",
    "B",
    "B",
    "deferred_non_B",
    "B",
    "deferred_non_B"
  ),

  biological_evidence_confidence=c(
    "high",
    "high",
    "high",
    "medium",
    "high",
    "low_medium",
    "low",
    "medium",
    "low",
    "high",
    "low_medium",
    "high"
  ),

  evidence_basis=c(
    "8-project project-balanced native naive identity program",
    "4-project project-balanced native memory identity program",
    "4-of-5 projects strongly support native class-switched program; single 10-cell GSE167363 exception",
    "native atypical activation program supported by strict within-library audit; secondary project audit mixed in GSE167363; retained as candidate",
    "strong native interferon program in GSE163668 and GSE167363",
    "strict within-library audit slightly negative; secondary project-matched audit weakly positive; downgraded to candidate",
    "pilot secretory differentiation evidence retained but full project-matched module effect weak; remains candidate with low confidence",
    "native activation/stress program positive in strict and secondary GSE220189 audits; additional positive GSE151263 signal",
    "no positive biological state claim",
    "strong native platelet program across four projects",
    "pilot IFN evidence present but full independent project-matched audit unavailable; downgraded to candidate",
    "strong native T-cell program across four projects"
  ),

  stringsAsFactors=FALSE
)

stopifnot(
  !anyDuplicated(taxonomy$taxonomy_id)
)

# ============================================================
# Map final taxonomy
# ============================================================

id <- as.character(
  ann$b_plasma_final_taxonomy_id_v1
)

resolved <- id != "TRANSFER_UNRESOLVED"

ii <- match(
  id,
  taxonomy$taxonomy_id
)

stopifnot(
  all(!is.na(ii[resolved])),
  all(is.na(ii[!resolved]))
)

ann$b_plasma_atlas_core_identity_v1 <-
  "B_plasma_transfer_unresolved"

ann$b_plasma_atlas_state_v1 <-
  "unresolved"

ann$b_plasma_atlas_annotation_scope_v1 <-
  "unresolved"

ann$b_plasma_atlas_biological_evidence_confidence_v1 <-
  "unresolved"

ann$b_plasma_atlas_annotation_decision_basis_v1 <-
  "transfer_score_below_frozen_accept_threshold"

ann$b_plasma_atlas_core_identity_v1[resolved] <-
  taxonomy$final_core_identity[
    ii[resolved]
  ]

ann$b_plasma_atlas_state_v1[resolved] <-
  taxonomy$final_state[
    ii[resolved]
  ]

ann$b_plasma_atlas_annotation_scope_v1[resolved] <-
  taxonomy$final_annotation_scope[
    ii[resolved]
  ]

ann$b_plasma_atlas_biological_evidence_confidence_v1[resolved] <-
  taxonomy$biological_evidence_confidence[
    ii[resolved]
  ]

ann$b_plasma_atlas_annotation_decision_basis_v1[resolved] <-
  taxonomy$evidence_basis[
    ii[resolved]
  ]

# ============================================================
# Preserve transfer confidence independently
# ============================================================

ann$b_plasma_atlas_transfer_confidence_v1 <-
  as.character(
    ann$b_plasma_transfer_confidence_v1
  )

ann$b_plasma_atlas_prediction_score_v1 <-
  ann$b_plasma_prediction_score_v1

ann$b_plasma_atlas_prediction_margin_v1 <-
  ann$b_plasma_prediction_margin_v1

ann$b_plasma_atlas_annotation_source_v1 <-
  as.character(
    ann$b_plasma_annotation_source_v1
  )

# Frozen pilot cells have no prediction score by design.
stopifnot(
  all(
    ann$b_plasma_atlas_annotation_source_v1 %in%
      c(
        "pilot_freeze_v1",
        "full_transfer_v1"
      )
  )
)

# ============================================================
# Integrity
# ============================================================

stopifnot(
  nrow(ann) == 45955L,
  !anyNA(
    ann$b_plasma_atlas_core_identity_v1
  ),
  !anyNA(
    ann$b_plasma_atlas_state_v1
  ),
  !anyNA(
    ann$b_plasma_atlas_annotation_scope_v1
  ),
  !anyNA(
    ann$b_plasma_atlas_biological_evidence_confidence_v1
  )
)

n_unresolved <- sum(
  ann$b_plasma_atlas_core_identity_v1 ==
    "B_plasma_transfer_unresolved"
)

stopifnot(
  n_unresolved == 2370L
)

# ============================================================
# Summary tables
# ============================================================

summary_taxonomy <- as.data.frame(
  table(
    core_identity=
      ann$b_plasma_atlas_core_identity_v1,
    state=
      ann$b_plasma_atlas_state_v1
  ),
  stringsAsFactors=FALSE
)

summary_taxonomy <- summary_taxonomy[
  summary_taxonomy$Freq > 0,
  ,
  drop=FALSE
]

summary_evidence <- as.data.frame(
  table(
    biological_evidence_confidence=
      ann$b_plasma_atlas_biological_evidence_confidence_v1
  ),
  stringsAsFactors=FALSE
)

summary_transfer <- as.data.frame(
  table(
    transfer_confidence=
      ann$b_plasma_atlas_transfer_confidence_v1
  ),
  stringsAsFactors=FALSE
)

# ============================================================
# Write lossless primary object
# ============================================================

saveRDS(
  ann,
  file=file.path(
    out_dir,
    "b_plasma_full_annotation_freeze_v1.rds"
  ),
  compress=FALSE
)

write.table(
  taxonomy,
  file=file.path(
    out_dir,
    "b_plasma_final_taxonomy_definition_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  summary_taxonomy,
  file=file.path(
    out_dir,
    "b_plasma_final_annotation_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  summary_evidence,
  file=file.path(
    out_dir,
    "b_plasma_final_biological_evidence_confidence_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

write.table(
  summary_transfer,
  file=file.path(
    out_dir,
    "b_plasma_final_transfer_confidence_counts_v1.tsv"
  ),
  sep="\t",
  quote=TRUE,
  row.names=FALSE
)

# ============================================================
# Re-read frozen RDS
# ============================================================

chk <- readRDS(
  file.path(
    out_dir,
    "b_plasma_full_annotation_freeze_v1.rds"
  )
)

stopifnot(
  nrow(chk) == 45955L,
  identical(
    chk$cell_id,
    ann$cell_id
  ),
  identical(
    chk$b_plasma_atlas_core_identity_v1,
    ann$b_plasma_atlas_core_identity_v1
  ),
  identical(
    chk$b_plasma_atlas_state_v1,
    ann$b_plasma_atlas_state_v1
  )
)

rm(chk)

# ============================================================
# Print
# ============================================================

cat(
  "\n===== FINAL B/PLASMA ATLAS TAXONOMY =====\n"
)

print(
  taxonomy,
  row.names=FALSE
)

cat(
  "\n===== FINAL CELL COUNTS =====\n"
)

print(
  summary_taxonomy,
  row.names=FALSE
)

cat(
  "\n===== BIOLOGICAL EVIDENCE CONFIDENCE =====\n"
)

print(
  summary_evidence,
  row.names=FALSE
)

cat(
  "\n===== TRANSFER CONFIDENCE =====\n"
)

print(
  summary_transfer,
  row.names=FALSE
)

# ============================================================
# Completion marker
# ============================================================

writeLines(
  c(
    "PASS",
    "B/plasma full Atlas annotation freeze v1",
    "n_full_cells=45955",
    "n_transfer_unresolved=2370",
    "",
    "core identity and orthogonal state stored separately",
    "transfer confidence and biological evidence confidence stored separately",
    "",
    "post-audit revisions:",
    "T06 early_activation -> early_activation_candidate",
    "T07 secretory_differentiation_candidate retained but biological confidence downgraded to low",
    "T11 interferon_stimulated -> interferon_stimulated_candidate",
    "",
    "T09 remains unresolved_state",
    "T10/T12 remain deferred non-B",
    "",
    "native evidence sources:",
    "strict within-library library-pseudobulk audit",
    "secondary within-project library-pseudobulk audit",
    "",
    "this freeze does not alter native RNA expression"
  ),
  done_file
)

cat(
  "\nPASS: B/plasma full annotation freeze v1 completed\n"
)

