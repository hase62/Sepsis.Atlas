#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$1" && pwd)"
TAG="$2"
OUT_DIR="$3"

ATLAS_ROOT="$ROOT/atlas/full_atlas_annotation_v1"

mkdir -p \
  "$OUT_DIR/release" \
  "$OUT_DIR/audit" \
  "$OUT_DIR/provenance"

DONE="$OUT_DIR/ATLAS_ANNOTATION_RELEASE_CANDIDATE_V1_COMPLETE.ok"

if [[ -e "$DONE" ]]; then
  echo "ERROR: release candidate already completed:"
  echo "$DONE"
  exit 1
fi

latest_dir_by_marker() {
  local marker="$1"
  local result

  result="$(
    find "$ATLAS_ROOT" \
      -maxdepth 2 \
      -type f \
      -name "$marker" \
      -printf '%T@\t%h\n' \
    | sort -nr \
    | head -n 1 \
    | cut -f2-
  )"

  if [[ -z "$result" ]]; then
    echo "ERROR: marker not found: $marker" >&2
    exit 1
  fi

  printf '%s\n' "$result"
}

single_file() {
  local dir="$1"
  local pattern="$2"

  local -a files=()

  while IFS= read -r -d '' f; do
    files+=("$f")
  done < <(
    find "$dir" \
      -maxdepth 1 \
      -type f \
      -name "$pattern" \
      -print0
  )

  if [[ "${#files[@]}" -ne 1 ]]; then
    echo "ERROR: expected exactly one file" >&2
    echo "DIR=$dir" >&2
    echo "PATTERN=$pattern" >&2
    echo "FOUND=${#files[@]}" >&2
    exit 1
  fi

  printf '%s\n' "${files[0]}"
}

# ============================================================
# Resolve canonical inputs
# ============================================================

B_DIR="$(
  latest_dir_by_marker \
    'B_PLASMA_FULL_ANNOTATION_FREEZE_COMPLETE.ok'
)"

TNK_DIR="$ATLAS_ROOT/tnk_annotation_freeze_v1"
MONODC_DIR="$ATLAS_ROOT/monodc_annotation_freeze_v1"

NEUT_DIR="$(
  latest_dir_by_marker \
    'NEUTROPHIL_ANNOTATION_FREEZE_COMPLETE.ok'
)"

PLATELET_DIR="$(
  latest_dir_by_marker \
    'PLATELET_FULL_ANNOTATION_FREEZE_COMPLETE.ok'
)"

ERY_DIR="$(
  latest_dir_by_marker \
    'ERYTHROID_FULL_ANNOTATION_FREEZE_COMPLETE.ok'
)"

PROG_DIR="$(
  latest_dir_by_marker \
    'PROGENITOR_FULL_ANNOTATION_FREEZE_COMPLETE.ok'
)"

RECON_DIR="$(
  latest_dir_by_marker \
    'GLOBAL_DEFERRED_RECONCILIATION_COMPLETE.ok'
)"

ASSEMBLY_V1_DIR="$(
  latest_dir_by_marker \
    'ATLAS_FINAL_ANNOTATION_ASSEMBLY_COMPLETE.ok'
)"

ASSEMBLY_V2_DIR="$(
  latest_dir_by_marker \
    'ATLAS_FINAL_ANNOTATION_ASSEMBLY_V2_COMPLETE.ok'
)"

QC_V2_DIR="$(
  latest_dir_by_marker \
    'ATLAS_GLOBAL_QC_INTEGRITY_AUDIT_V2_COMPLETE.ok'
)"

for d in \
  "$B_DIR" \
  "$TNK_DIR" \
  "$MONODC_DIR" \
  "$NEUT_DIR" \
  "$PLATELET_DIR" \
  "$ERY_DIR" \
  "$PROG_DIR" \
  "$RECON_DIR" \
  "$ASSEMBLY_V1_DIR" \
  "$ASSEMBLY_V2_DIR" \
  "$QC_V2_DIR"
do
  if [[ ! -d "$d" ]]; then
    echo "ERROR: missing canonical directory:"
    echo "$d"
    exit 1
  fi
done

# ============================================================
# Completion markers
# ============================================================

test -f \
  "$B_DIR/B_PLASMA_FULL_ANNOTATION_FREEZE_COMPLETE.ok"

test -f \
  "$TNK_DIR/ANNOTATION_FREEZE_COMPLETE.ok"

test -f \
  "$MONODC_DIR/ANNOTATION_FREEZE_COMPLETE.ok"

test -f \
  "$NEUT_DIR/NEUTROPHIL_ANNOTATION_FREEZE_COMPLETE.ok"

test -f \
  "$PLATELET_DIR/PLATELET_FULL_ANNOTATION_FREEZE_COMPLETE.ok"

test -f \
  "$ERY_DIR/ERYTHROID_FULL_ANNOTATION_FREEZE_COMPLETE.ok"

test -f \
  "$PROG_DIR/PROGENITOR_FULL_ANNOTATION_FREEZE_COMPLETE.ok"

test -f \
  "$RECON_DIR/GLOBAL_DEFERRED_RECONCILIATION_COMPLETE.ok"

test -f \
  "$ASSEMBLY_V1_DIR/ATLAS_FINAL_ANNOTATION_ASSEMBLY_COMPLETE.ok"

test -f \
  "$ASSEMBLY_V2_DIR/ATLAS_FINAL_ANNOTATION_ASSEMBLY_V2_COMPLETE.ok"

test -f \
  "$QC_V2_DIR/ATLAS_GLOBAL_QC_INTEGRITY_AUDIT_V2_COMPLETE.ok"

# ============================================================
# Locate final release files
# ============================================================

FINAL_TSV="$(
  single_file \
    "$ASSEMBLY_V2_DIR" \
    'atlas_final_cell_annotation_v2__*.tsv.gz'
)"

FINAL_RDS="$(
  single_file \
    "$ASSEMBLY_V2_DIR" \
    'atlas_final_cell_annotation_v2__*.rds'
)"

AUDIT_RDS="$(
  single_file \
    "$QC_V2_DIR" \
    'atlas_global_qc_integrity_audit_v2__*.rds'
)"

LEDGER="$(
  single_file \
    "$RECON_DIR" \
    'atlas_deferred_reconciliation_cell_ledger_v1__*.tsv.gz'
)"

# ============================================================
# Copy helper with exact-byte verification
# ============================================================

COPY_MANIFEST="$OUT_DIR/release_copy_manifest_v1.tsv"

printf \
  'role\tsource_path\trelease_path\tbytes\tsha256\n' \
  > "$COPY_MANIFEST"

copy_verified() {
  local role="$1"
  local src="$2"
  local subdir="$3"
  local dstname="$4"

  if [[ ! -f "$src" ]]; then
    echo "ERROR: missing source: $src"
    exit 1
  fi

  local dst="$OUT_DIR/$subdir/$dstname"

  cp -p -- "$src" "$dst"

  if ! cmp -s -- "$src" "$dst"; then
    echo "ERROR: copied file differs from source:"
    echo "$src"
    exit 1
  fi

  local sha
  local bytes

  sha="$(
    sha256sum "$dst" |
      awk '{print $1}'
  )"

  bytes="$(
    stat -c '%s' "$dst"
  )"

  printf \
    '%s\t%s\t%s\t%s\t%s\n' \
    "$role" \
    "$src" \
    "${dst#"$OUT_DIR/"}" \
    "$bytes" \
    "$sha" \
    >> "$COPY_MANIFEST"
}

# ============================================================
# Release payload
# ============================================================

copy_verified \
  "final_cell_annotation_tsv" \
  "$FINAL_TSV" \
  "release" \
  "atlas_final_cell_annotation_v2__RC_${TAG}.tsv.gz"

copy_verified \
  "final_cell_annotation_rds" \
  "$FINAL_RDS" \
  "release" \
  "atlas_final_cell_annotation_v2__RC_${TAG}.rds"

copy_verified \
  "final_compartment_counts" \
  "$ASSEMBLY_V2_DIR/atlas_final_compartment_counts_v2.tsv" \
  "release" \
  "atlas_final_compartment_counts_v2__RC_${TAG}.tsv"

copy_verified \
  "reconciliation_status_counts" \
  "$ASSEMBLY_V2_DIR/atlas_reconciliation_status_counts_v2.tsv" \
  "release" \
  "atlas_reconciliation_status_counts_v2__RC_${TAG}.tsv"

copy_verified \
  "upstream_deferred_counts" \
  "$ASSEMBLY_V2_DIR/atlas_upstream_deferred_reason_counts_v1.tsv" \
  "release" \
  "atlas_upstream_deferred_reason_counts_v1__RC_${TAG}.tsv"

copy_verified \
  "assembly_source_manifest" \
  "$ASSEMBLY_V2_DIR/atlas_final_annotation_v2_source_manifest.tsv" \
  "release" \
  "atlas_final_annotation_v2_source_manifest__RC_${TAG}.tsv"

# ============================================================
# Audit payload
# ============================================================

copy_verified \
  "global_qc_assertions" \
  "$QC_V2_DIR/atlas_global_qc_assertions_v2.tsv" \
  "audit" \
  "atlas_global_qc_assertions_v2__RC_${TAG}.tsv"

copy_verified \
  "cell_universe_concordance" \
  "$QC_V2_DIR/atlas_cell_universe_concordance_v2.tsv" \
  "audit" \
  "atlas_cell_universe_concordance_v2__RC_${TAG}.tsv"

copy_verified \
  "library_cell_count_concordance" \
  "$QC_V2_DIR/atlas_library_cell_count_concordance_v2.tsv" \
  "audit" \
  "atlas_library_cell_count_concordance_v2__RC_${TAG}.tsv"

copy_verified \
  "excluded_low_final_library" \
  "$QC_V2_DIR/atlas_excluded_low_final_library_v2.tsv" \
  "audit" \
  "atlas_excluded_low_final_library_v2__RC_${TAG}.tsv"

copy_verified \
  "qc_final_compartment_counts" \
  "$QC_V2_DIR/atlas_final_compartment_counts_qc_v2.tsv" \
  "audit" \
  "atlas_final_compartment_counts_qc_v2__RC_${TAG}.tsv"

copy_verified \
  "qc_reconciliation_status_counts" \
  "$QC_V2_DIR/atlas_reconciliation_status_counts_qc_v2.tsv" \
  "audit" \
  "atlas_reconciliation_status_counts_qc_v2__RC_${TAG}.tsv"

copy_verified \
  "global_qc_source_manifest" \
  "$QC_V2_DIR/atlas_global_qc_source_manifest_v2.tsv" \
  "audit" \
  "atlas_global_qc_source_manifest_v2__RC_${TAG}.tsv"

copy_verified \
  "global_qc_audit_rds" \
  "$AUDIT_RDS" \
  "audit" \
  "atlas_global_qc_integrity_audit_v2__RC_${TAG}.rds"

copy_verified \
  "deferred_reconciliation_ledger" \
  "$LEDGER" \
  "audit" \
  "atlas_deferred_reconciliation_cell_ledger_v1__RC_${TAG}.tsv.gz"

copy_verified \
  "deferred_reconciliation_rules" \
  "$RECON_DIR/atlas_deferred_reconciliation_rules_v1.tsv" \
  "audit" \
  "atlas_deferred_reconciliation_rules_v1__RC_${TAG}.tsv"

copy_verified \
  "deferred_reconciliation_label_counts" \
  "$RECON_DIR/atlas_deferred_reconciliation_label_counts_v1.tsv" \
  "audit" \
  "atlas_deferred_reconciliation_label_counts_v1__RC_${TAG}.tsv"

copy_verified \
  "deferred_reconciliation_action_counts" \
  "$RECON_DIR/atlas_deferred_reconciliation_action_counts_v1.tsv" \
  "audit" \
  "atlas_deferred_reconciliation_action_counts_v1__RC_${TAG}.tsv"

# ============================================================
# Upstream QC provenance
# ============================================================

copy_verified \
  "final_qc_library_manifest" \
  "$ROOT/pre_integration/final_qc_library_manifest.tsv" \
  "provenance" \
  "final_qc_library_manifest__RC_${TAG}.tsv"

copy_verified \
  "final_qc_decision_summary" \
  "$ROOT/pre_integration/final_qc_decision_summary.tsv" \
  "provenance" \
  "final_qc_decision_summary__RC_${TAG}.tsv"

copy_verified \
  "resolved_library_table" \
  "$ROOT/pre_integration/pilot_unintegrated_large_v1/resolved_library_table.tsv" \
  "provenance" \
  "resolved_library_table__RC_${TAG}.tsv"

copy_verified \
  "annotation_transfer_library_summary" \
  "$ROOT/pre_integration/full_atlas_primary_v1/full_annotation_transfer_library_summary.tsv" \
  "provenance" \
  "full_annotation_transfer_library_summary__RC_${TAG}.tsv"

copy_verified \
  "large_object_validation" \
  "$ROOT/pre_integration/full_atlas_primary_v1/freeze_v1__20260813_1312/final_object_validation_v1.tsv" \
  "provenance" \
  "final_object_validation_v1__RC_${TAG}.tsv"

# ============================================================
# Canonical source hash manifest
# ============================================================

INPUT_HASHES="$OUT_DIR/canonical_input_sha256_v1.tsv"

printf \
  'source_role\tfile_path\tbytes\tsha256\n' \
  > "$INPUT_HASHES"

hash_dir_top() {
  local role="$1"
  local dir="$2"

  while IFS= read -r -d '' f; do

    local sha
    local bytes

    sha="$(
      sha256sum "$f" |
        awk '{print $1}'
    )"

    bytes="$(
      stat -c '%s' "$f"
    )"

    printf \
      '%s\t%s\t%s\t%s\n' \
      "$role" \
      "$f" \
      "$bytes" \
      "$sha" \
      >> "$INPUT_HASHES"

  done < <(
    find "$dir" \
      -maxdepth 1 \
      -type f \
      -print0 |
    sort -z
  )
}

hash_dir_top "B_plasma_freeze" "$B_DIR"
hash_dir_top "T_NK_freeze" "$TNK_DIR"
hash_dir_top "Monocyte_DC_freeze" "$MONODC_DIR"
hash_dir_top "Neutrophil_freeze" "$NEUT_DIR"
hash_dir_top "Platelet_megakaryocyte_freeze" "$PLATELET_DIR"
hash_dir_top "Erythroid_freeze" "$ERY_DIR"
hash_dir_top "Progenitor_freeze" "$PROG_DIR"

hash_dir_top "Step44_reconciliation" "$RECON_DIR"
hash_dir_top "Step45_v1_assembly" "$ASSEMBLY_V1_DIR"
hash_dir_top "Step45_v2_assembly" "$ASSEMBLY_V2_DIR"
hash_dir_top "Step46_v2_global_QC" "$QC_V2_DIR"

TRANSFER_DIR="$ROOT/pre_integration/full_atlas_primary_v1/annotation_transfer_by_library"

TRANSFER_N="$(
  find "$TRANSFER_DIR" \
    -maxdepth 1 \
    -type f \
    -name '*.tsv.gz' |
  wc -l
)"

if [[ "$TRANSFER_N" -ne 158 ]]; then
  echo "ERROR: expected 158 annotation-transfer files, found $TRANSFER_N"
  exit 1
fi

hash_dir_top \
  "annotation_transfer_by_library" \
  "$TRANSFER_DIR"

# ============================================================
# Canonical directory manifest
# ============================================================

DIR_MANIFEST="$OUT_DIR/canonical_directory_manifest_v1.tsv"

cat > "$DIR_MANIFEST" <<EOF
rolepath
B_plasma_freeze$B_DIR
T_NK_freeze$TNK_DIR
Monocyte_DC_freeze$MONODC_DIR
Neutrophil_freeze$NEUT_DIR
Platelet_megakaryocyte_freeze$PLATELET_DIR
Erythroid_freeze$ERY_DIR
Progenitor_freeze$PROG_DIR
Step44_reconciliation$RECON_DIR
Step45_v1_assembly$ASSEMBLY_V1_DIR
Step45_v2_assembly$ASSEMBLY_V2_DIR
Step46_v2_global_QC$QC_V2_DIR
annotation_transfer_by_library$TRANSFER_DIR
EOF

# ============================================================
# Release metadata
# ============================================================

cat > "$OUT_DIR/RELEASE_CANDIDATE_METADATA_v1.tsv" <<EOF
fieldvalue
release_candidateatlas_annotation_release_candidate_v1
release_tag$TAG
n_final_qc_cells665816
n_unique_global_cells665816
n_included_libraries158
n_excluded_low_final_cells_libraries1
excluded_low_final_libraryGSE163668::GSM4995445__SAMN17184537
excluded_low_final_library_n_final132
n_compartment_annotated_cells662541
n_upstream_primary_deferred_cells3275
n_step44_global_deferred_cells6272
n_final_deferred_unresolved9547
n_compartment_reassigned2712
n_source_retained_ambiguous4295
annotation_granularityfrozen
new_clustering_after_compartment_freezeno
new_marker_discovery_after_compartment_freezeno
core_taxonomy_changedno
axis_taxonomy_changedno
state_taxonomy_changedno
exact_final_qc_cell_universePASS
global_cell_duplicates0
minimum_final_cells_rule200
step46_all_assertionsPASS
EOF

cat > "$OUT_DIR/README_ANNOTATION_RELEASE_CANDIDATE_v1.txt" <<EOF
Sepsis.Atlas annotation release candidate v1
Release tag: $TAG

This release candidate freezes the completed cell-annotation and
global-QC state of the Atlas.

Final Atlas universe:
  665816 final-QC cells
  158 included libraries
  0 duplicated global_cell identifiers

Seven compartment-specific annotation workflows contain 662541 cells.
An additional 3275 cells were intentionally deferred during upstream
primary-compartment assignment and are retained as Deferred_unresolved
without inventing compartment-specific core, axis, or state labels.

Final Deferred_unresolved:
  6272 cells from Step44 global reconciliation
  3275 upstream primary-deferred cells
  9547 cells total

The existing compartment-specific annotation granularity is frozen.
No new clustering, marker discovery, subtype refinement, taxonomy
splitting, or taxonomy collapsing is performed by this release step.

Step46 v2 verified exact concordance between the 665816-cell Atlas and
the 665816-cell final annotation-transfer universe.

The final-QC rule excluding libraries with fewer than 200 retained
cells was verified. One library was excluded:
  GSE163668::GSM4995445__SAMN17184537
  n_final = 132

This directory is a release-candidate freeze, not a new analysis.
EOF

# ============================================================
# Deep validation of copied release object
# ============================================================

RC_TSV="$OUT_DIR/release/atlas_final_cell_annotation_v2__RC_${TAG}.tsv.gz"
RC_RDS="$OUT_DIR/release/atlas_final_cell_annotation_v2__RC_${TAG}.rds"
RC_LEDGER="$OUT_DIR/audit/atlas_deferred_reconciliation_cell_ledger_v1__RC_${TAG}.tsv.gz"

gzip -t "$RC_TSV"
gzip -t "$RC_LEDGER"

conda run \
  -n sepsis-scrna-r43-seurat551 \
  --no-capture-output \
  Rscript - \
  "$RC_RDS" \
  "$RC_TSV" <<'RS'

args <- commandArgs(trailingOnly=TRUE)

rds_file <- args[[1]]
tsv_file <- args[[2]]

r <- readRDS(rds_file)

stopifnot(
  is.data.frame(r),
  nrow(r) == 665816L,
  !anyDuplicated(r$global_cell)
)

expected <- c(
  B_plasma=45154L,
  Deferred_unresolved=9547L,
  Erythroid=11526L,
  Monocyte_DC=192040L,
  Neutrophil=51469L,
  Platelet_megakaryocyte=47470L,
  Progenitor=27071L,
  T_NK=281539L
)

obs <- table(
  r$final_compartment_v1
)

stopifnot(
  setequal(
    names(obs),
    names(expected)
  ),
  all(
    as.integer(
      obs[names(expected)]
    ) ==
      expected
  )
)

con <- gzfile(
  tsv_file,
  "rt"
)

t <- tryCatch(
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

stopifnot(
  nrow(t) == 665816L,
  !anyDuplicated(t$global_cell),
  identical(
    as.character(t$global_cell),
    as.character(r$global_cell)
  ),
  identical(
    as.character(t$final_compartment_v1),
    as.character(r$final_compartment_v1)
  )
)

cat("PASS: release RDS/TSV deep validation\n")
RS

# ============================================================
# Release file manifest
# ============================================================

RELEASE_MANIFEST="$OUT_DIR/RELEASE_FILE_MANIFEST_v1.tsv"

printf \
  'relative_path\tbytes\tsha256\n' \
  > "$RELEASE_MANIFEST"

while IFS= read -r -d '' f; do

  rel="${f#"$OUT_DIR/"}"

  case "$rel" in
    SHA256SUMS_v1__*.txt)
      continue
      ;;
    ATLAS_ANNOTATION_RELEASE_CANDIDATE_V1_COMPLETE.ok)
      continue
      ;;
  esac

  bytes="$(
    stat -c '%s' "$f"
  )"

  sha="$(
    sha256sum "$f" |
      awk '{print $1}'
  )"

  printf \
    '%s\t%s\t%s\n' \
    "$rel" \
    "$bytes" \
    "$sha" \
    >> "$RELEASE_MANIFEST"

done < <(
  find "$OUT_DIR" \
    -type f \
    -print0 |
  sort -z
)

# ============================================================
# Final SHA256 manifest
# ============================================================

SHA_FILE="$OUT_DIR/SHA256SUMS_v1__${TAG}.txt"

(
  cd "$OUT_DIR"

  find . \
    -type f \
    ! -name "SHA256SUMS_v1__${TAG}.txt" \
    ! -name "ATLAS_ANNOTATION_RELEASE_CANDIDATE_V1_COMPLETE.ok" \
    -print0 |
  sort -z |
  xargs -0 sha256sum \
    > "SHA256SUMS_v1__${TAG}.txt"

  sha256sum \
    -c \
    "SHA256SUMS_v1__${TAG}.txt" \
    >/dev/null
)

# ============================================================
# Completion marker LAST
# ============================================================

cat > "$DONE" <<EOF
PASS
Atlas annotation release candidate v1

release_tag=$TAG

n_final_qc_cells=665816
n_unique_global_cells=665816
n_included_libraries=158
n_excluded_low_final_cells_libraries=1

n_compartment_annotated_cells=662541
n_upstream_primary_deferred_cells=3275
n_final_deferred_unresolved=9547

exact final-QC cell universe=PASS
Step44 reconciliation preservation=PASS
Step45-v1 annotation preservation=PASS
Step45-v2 full-universe assembly=PASS
Step46-v2 all assertions=PASS

annotation granularity frozen
no new clustering
no new marker discovery
no subtype refinement
no taxonomy split
no taxonomy collapse

release TSV/RDS deep validation=PASS
gzip integrity=PASS
copy byte identity=PASS
canonical input SHA256 provenance=PASS
release SHA256 manifest=PASS

CELL ANNOTATION / GLOBAL QC PIPELINE COMPLETE
EOF

# ============================================================
# Freeze release candidate read-only
# ============================================================

find "$OUT_DIR" \
  -type f \
  -exec chmod 0444 {} +

find "$OUT_DIR" \
  -type d \
  -exec chmod 0555 {} +

echo
echo "PASS: Atlas annotation release candidate frozen"
echo "OUT_DIR=$OUT_DIR"
