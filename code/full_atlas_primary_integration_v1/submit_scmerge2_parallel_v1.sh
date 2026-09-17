#!/usr/bin/env bash
set -euo pipefail

ROOT="${1:-$PWD}"
ROOT="$(cd "${ROOT}" && pwd -P)"

SCRIPT_DIR="${ROOT}/full_atlas_primary_integration_v1"
FULL_DIR="${ROOT}/pre_integration/full_atlas_primary_v1"
WORKER="${SCRIPT_DIR}/qsub_scmerge2_parallel_worker_v1.sh"

CHUNK_SIZE="${CHUNK_SIZE:-10000}"
MAX_CONCURRENT_PER_COMPARTMENT="${MAX_CONCURRENT_PER_COMPARTMENT:-8}"

LOG_DIR="${FULL_DIR}/qsub_logs/parallel_v1"
mkdir -p "${LOG_DIR}"

COUNTS="${FULL_DIR}/full_annotation_transfer_compartment_counts.tsv"

if [[ ! -f "${COUNTS}" ]]; then
  echo "ERROR: missing ${COUNTS}" >&2
  exit 2
fi

if [[ ! -x "${WORKER}" ]]; then
  echo "ERROR: worker is not executable: ${WORKER}" >&2
  exit 2
fi

get_n_cells () {
  local compartment="$1"
  awk -F '\t' -v c="${compartment}" '
    NR > 1 && $1 == c {print $2}
  ' "${COUNTS}"
}

n_chunks () {
  local n="$1"
  echo $(( (n + CHUNK_SIZE - 1) / CHUNK_SIZE ))
}

submit_and_get_id () {
  local out id
  out="$(qsub "$@")"
  echo "${out}" >&2

  id="$(
    printf '%s\n' "${out}" |
      sed -nE 's/.*job(-array)?[[:space:]]+([0-9]+).*/\2/p' |
      head -1
  )"

  if [[ -z "${id}" ]]; then
    echo "ERROR: could not parse qsub job ID" >&2
    return 1
  fi

  printf '%s\n' "${id}"
}

T_NK_CELLS="$(get_n_cells T_NK_combined)"
MONO_CELLS="$(get_n_cells Monocyte_DC_combined)"

if [[ -z "${T_NK_CELLS}" || -z "${MONO_CELLS}" ]]; then
  echo "ERROR: could not resolve compartment counts" >&2
  exit 3
fi

T_NK_CHUNKS="$(n_chunks "${T_NK_CELLS}")"
MONO_CHUNKS="$(n_chunks "${MONO_CELLS}")"

echo "=== PARALLEL SCMERGE2 SUBMISSION ==="
echo "ROOT=${ROOT}"
echo "CHUNK_SIZE=${CHUNK_SIZE}"
echo "MAX_CONCURRENT_PER_COMPARTMENT=${MAX_CONCURRENT_PER_COMPARTMENT}"
echo "T_NK cells=${T_NK_CELLS}; chunks=${T_NK_CHUNKS}"
echo "Mono/DC cells=${MONO_CELLS}; chunks=${MONO_CHUNKS}"
echo

TNK_PREP="$(
  submit_and_get_id \
    -V -cwd -S /bin/bash \
    -N s2p_tnk \
    -l s_vmem=40G \
    -l mem_req=40G \
    -o "${LOG_DIR}" \
    -e "${LOG_DIR}" \
    "${WORKER}" \
    prepare "${ROOT}" T_NK_combined "${CHUNK_SIZE}"
)"

MONO_PREP="$(
  submit_and_get_id \
    -V -cwd -S /bin/bash \
    -N s2p_mdc \
    -l s_vmem=40G \
    -l mem_req=40G \
    -o "${LOG_DIR}" \
    -e "${LOG_DIR}" \
    "${WORKER}" \
    prepare "${ROOT}" Monocyte_DC_combined "${CHUNK_SIZE}"
)"

TNK_ARRAY="$(
  submit_and_get_id \
    -V -cwd -S /bin/bash \
    -N s2a_tnk \
    -hold_jid "${TNK_PREP}" \
    -t "1-${T_NK_CHUNKS}" \
    -tc "${MAX_CONCURRENT_PER_COMPARTMENT}" \
    -l s_vmem=12G \
    -l mem_req=12G \
    -o "${LOG_DIR}" \
    -e "${LOG_DIR}" \
    "${WORKER}" \
    adjust "${ROOT}" T_NK_combined "${CHUNK_SIZE}"
)"

MONO_ARRAY="$(
  submit_and_get_id \
    -V -cwd -S /bin/bash \
    -N s2a_mdc \
    -hold_jid "${MONO_PREP}" \
    -t "1-${MONO_CHUNKS}" \
    -tc "${MAX_CONCURRENT_PER_COMPARTMENT}" \
    -l s_vmem=12G \
    -l mem_req=12G \
    -o "${LOG_DIR}" \
    -e "${LOG_DIR}" \
    "${WORKER}" \
    adjust "${ROOT}" Monocyte_DC_combined "${CHUNK_SIZE}"
)"

TNK_FINAL="$(
  submit_and_get_id \
    -V -cwd -S /bin/bash \
    -N s2f_tnk \
    -hold_jid "${TNK_ARRAY}" \
    -l s_vmem=64G \
    -l mem_req=64G \
    -o "${LOG_DIR}" \
    -e "${LOG_DIR}" \
    "${WORKER}" \
    finalize "${ROOT}" T_NK_combined "${CHUNK_SIZE}"
)"

MONO_FINAL="$(
  submit_and_get_id \
    -V -cwd -S /bin/bash \
    -N s2f_mdc \
    -hold_jid "${MONO_ARRAY}" \
    -l s_vmem=64G \
    -l mem_req=64G \
    -o "${LOG_DIR}" \
    -e "${LOG_DIR}" \
    "${WORKER}" \
    finalize "${ROOT}" Monocyte_DC_combined "${CHUNK_SIZE}"
)"

VERIFY="$(
  submit_and_get_id \
    -V -cwd -S /bin/bash \
    -N s2_verify \
    -hold_jid "${TNK_FINAL},${MONO_FINAL}" \
    -l s_vmem=16G \
    -l mem_req=16G \
    -o "${LOG_DIR}" \
    -e "${LOG_DIR}" \
    "${WORKER}" \
    verify "${ROOT}"
)"

cat > "${FULL_DIR}/parallel_v1_submitted_jobs.tsv" <<EOF
stagecompartmentjob_idn_tasks
prepareT_NK_combined${TNK_PREP}1
prepareMonocyte_DC_combined${MONO_PREP}1
adjustT_NK_combined${TNK_ARRAY}${T_NK_CHUNKS}
adjustMonocyte_DC_combined${MONO_ARRAY}${MONO_CHUNKS}
finalizeT_NK_combined${TNK_FINAL}1
finalizeMonocyte_DC_combined${MONO_FINAL}1
verifyall${VERIFY}1
EOF

echo
echo "PASS: dependency chain submitted"
cat "${FULL_DIR}/parallel_v1_submitted_jobs.tsv"
echo
echo "Monitor with:"
echo "  qstat -u ${USER}"
echo "Logs:"
echo "  ${LOG_DIR}"
