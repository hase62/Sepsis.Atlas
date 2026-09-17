#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
ENV_YML="${1:-$SCRIPT_DIR/seurat_r43_conda_environment.yml}"
ENV_NAME="${2:-sepsis-scrna-r43}"

if command -v mamba >/dev/null 2>&1; then
  mamba env create -f "$ENV_YML" -n "$ENV_NAME"
elif command -v conda >/dev/null 2>&1; then
  conda env create -f "$ENV_YML" -n "$ENV_NAME"
else
  echo "Neither mamba nor conda was found in PATH." >&2
  exit 1
fi

echo "Activate with: conda activate $ENV_NAME"
