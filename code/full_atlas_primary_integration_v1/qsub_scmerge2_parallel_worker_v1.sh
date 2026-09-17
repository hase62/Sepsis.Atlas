#!/usr/bin/env bash
set -euo pipefail

MODE="${1:?mode required}"
ROOT="${2:?root required}"
COMPARTMENT="${3:-}"
CHUNK_SIZE="${4:-10000}"

ROOT="$(cd "${ROOT}" && pwd -P)"
SCRIPT_DIR="${ROOT}/full_atlas_primary_integration_v1"

TASK="${SGE_TASK_ID:-0}"
JOB="${JOB_ID:-manual}"

TMPDIR="${ROOT}/tmp/scmerge2_parallel_v1/${COMPARTMENT:-verify}/${JOB}.${TASK}"
mkdir -p "${TMPDIR}"

export TMPDIR
export TMP="${TMPDIR}"
export TEMP="${TMPDIR}"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

CONDA_BIN="${CONDA_EXE:-}"
if [[ -z "${CONDA_BIN}" || ! -x "${CONDA_BIN}" ]]; then
  CONDA_BIN="$(command -v conda || true)"
fi

if [[ -z "${CONDA_BIN}" ]]; then
  echo "ERROR: conda executable not found" >&2
  exit 2
fi

run_r () {
  "${CONDA_BIN}" run \
    -n sepsis-scrna-r43 \
    --no-capture-output \
    Rscript "$@"
}

echo "=== SCMERGE2 PARALLEL WORKER ==="
date
hostname
echo "MODE=${MODE}"
echo "COMPARTMENT=${COMPARTMENT:-NA}"
echo "JOB_ID=${JOB}"
echo "SGE_TASK_ID=${TASK}"
echo "TMPDIR=${TMPDIR}"

case "${MODE}" in
  prepare)
    run_r \
      "${SCRIPT_DIR}/02a_prepare_scmerge2_parallel_v1.R" \
      "${ROOT}" \
      "${COMPARTMENT}" \
      "${CHUNK_SIZE}"
    ;;

  adjust)
    if [[ "${TASK}" -lt 1 ]]; then
      echo "ERROR: adjust mode requires SGE_TASK_ID" >&2
      exit 3
    fi

    run_r \
      "${SCRIPT_DIR}/02b_adjust_scmerge2_chunk_parallel_v1.R" \
      "${ROOT}" \
      "${COMPARTMENT}" \
      "${TASK}"
    ;;

  finalize)
    run_r \
      "${SCRIPT_DIR}/02c_finalize_scmerge2_parallel_v1.R" \
      "${ROOT}" \
      "${COMPARTMENT}"
    ;;

  verify)
    run_r \
      "${SCRIPT_DIR}/02d_collect_scmerge2_parallel_v1.R" \
      "${ROOT}"

    run_r \
      "${SCRIPT_DIR}/03_verify_full_primary.R" \
      "${ROOT}"
    ;;

  *)
    echo "ERROR: unknown mode: ${MODE}" >&2
    exit 4
    ;;
esac

rm -rf "${TMPDIR}"

date
echo "PASS MODE=${MODE} COMPARTMENT=${COMPARTMENT:-NA} TASK=${TASK}"
