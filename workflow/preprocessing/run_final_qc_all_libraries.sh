#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd -P)"
ROOT="${1:-$REPO_ROOT}"
ROOT="$(cd "$ROOT" && pwd -P)"
cd "$ROOT" || exit 1

PROJECTS=(
  GSE151263 GSE163668 GSE167363 GSE216007 GSE216020
  GSE217906 GSE220189 GSE242127 GSE252331
)
if [[ -n "${ATLAS_FINAL_QC_PROJECTS:-}" ]]; then
  PROJECT_TEXT="${ATLAS_FINAL_QC_PROJECTS//,/ }"
  read -r -a PROJECTS <<< "$PROJECT_TEXT"
fi

export TMPDIR="$ROOT/tmp/r_tmp_final_qc"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"
mkdir -p "$TMPDIR" "$ROOT/tmp/logs"

export ATLAS_MIQC_POSTERIOR_CUTOFF="${ATLAS_MIQC_POSTERIOR_CUTOFF:-0.90}"
export ATLAS_MIQC_NSTARTS="${ATLAS_MIQC_NSTARTS:-3}"
export ATLAS_MIQC_MIN_CELLS="${ATLAS_MIQC_MIN_CELLS:-200}"
export ATLAS_SCDOUBLET_DBR_PER1K="${ATLAS_SCDOUBLET_DBR_PER1K:-0.008}"
export ATLAS_FINAL_MIN_CELLS_REVIEW="${ATLAS_FINAL_MIN_CELLS_REVIEW:-200}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"

STAMP="$(date +%Y%m%d_%H%M%S)"
LOG_ROOT="$ROOT/tmp/logs/final_qc_${STAMP}"
mkdir -p "$LOG_ROOT"
SUMMARY="$LOG_ROOT/batch_summary.tsv"
printf 'project_id\tlibrary_key\texit_code\tstatus\treason\tinput_rds\toutput_rds\tlog\n' > "$SUMMARY"

package_check() {
  Rscript - <<'RS'
pkgs <- c("Seurat", "SeuratObject", "Matrix", "SingleCellExperiment", "scDblFinder", "BiocParallel", "scater", "miQC", "flexmix")
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing final-QC packages: ", paste(missing, collapse = ", "))
cat("PASS: final-QC package preflight\n")
RS
}

if ! package_check; then
  echo "Package preflight failed. Run workflow/environment/install_final_qc_packages_r43.R first." >&2
  exit 1
fi

for P in "${PROJECTS[@]}"; do
  INPUT_DIR="$ROOT/pre_integration/$P/rds_preqc_raw"
  PLOG="$LOG_ROOT/$P"
  mkdir -p "$PLOG"
  if [[ ! -d "$INPUT_DIR" ]]; then
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$P" "NA" "1" "failed" "input_directory_missing" "$INPUT_DIR" "NA" "$PLOG" >> "$SUMMARY"
    continue
  fi

  while IFS= read -r -d '' INPUT; do
    BN="$(basename "$INPUT")"
    KEY="${BN%__preqc_raw.rds}"
    if [[ "$KEY" == "$BN" ]]; then KEY="${BN%.rds}"; fi
    LOG="$PLOG/${KEY}.log"
    SUM="$ROOT/pre_integration/$P/qc/final_qc/${KEY}__final_qc_summary.tsv"
    OUT="$ROOT/pre_integration/$P/rds_final_qc/${KEY}__final_qc.rds"

    if [[ "${ATLAS_FINAL_QC_OVERWRITE:-false}" != "true" && -s "$SUM" ]]; then
      STATUS="$(awk -F '\t' 'NR==2 {for(i=1;i<=NF;i++) if(h[i]=="status") print $i} NR==1 {for(i=1;i<=NF;i++) h[i]=$i}' "$SUM")"
      REASON="$(awk -F '\t' 'NR==2 {for(i=1;i<=NF;i++) if(h[i]=="reason") print $i} NR==1 {for(i=1;i<=NF;i++) h[i]=$i}' "$SUM")"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$P" "$KEY" "0" "skipped_${STATUS:-existing}" "${REASON:-existing_summary}" "$INPUT" "$OUT" "$LOG" >> "$SUMMARY"
      continue
    fi

    if [[ "${ATLAS_FINAL_QC_OVERWRITE:-false}" == "true" ]]; then
      rm -f "$SUM" "$OUT" "$ROOT/pre_integration/$P/qc/final_qc/${KEY}__final_qc_cells.csv.gz"
    fi

    echo "START: $P / $KEY"
    Rscript "$ROOT/workflow/preprocessing/run_final_qc_one_library.R" --project-root "$ROOT" --input "$INPUT" > "$LOG" 2>&1
    EC=$?
    if [[ -s "$SUM" ]]; then
      STATUS="$(awk -F '\t' 'NR==2 {for(i=1;i<=NF;i++) if(h[i]=="status") print $i} NR==1 {for(i=1;i<=NF;i++) h[i]=$i}' "$SUM")"
      REASON="$(awk -F '\t' 'NR==2 {for(i=1;i<=NF;i++) if(h[i]=="reason") print $i} NR==1 {for(i=1;i<=NF;i++) h[i]=$i}' "$SUM")"
    else
      STATUS="failed"
      REASON="no_summary_written"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$P" "$KEY" "$EC" "$STATUS" "$REASON" "$INPUT" "$OUT" "$LOG" >> "$SUMMARY"
    echo "END: $P / $KEY status=$STATUS exit=$EC"
  done < <(find "$INPUT_DIR" -maxdepth 1 -type f -name '*.rds' -print0 | sort -z)
done

printf '\n================ FINAL QC BATCH SUMMARY ================\n'
column -t -s $'\t' "$SUMMARY" || cat "$SUMMARY"
printf 'Summary: %s\n' "$SUMMARY"

NFAIL=$(awk -F '\t' 'NR>1 && $4=="failed" {n++} END{print n+0}' "$SUMMARY")
NREVIEW=$(awk -F '\t' 'NR>1 && ($4=="review" || $4 ~ /^skipped_review/) {n++} END{print n+0}' "$SUMMARY")
NPASS=$(awk -F '\t' 'NR>1 && ($4=="pass" || $4 ~ /^skipped_pass/) {n++} END{print n+0}' "$SUMMARY")
echo "pass=$NPASS review=$NREVIEW failed=$NFAIL"
[[ "$NFAIL" -eq 0 ]]
