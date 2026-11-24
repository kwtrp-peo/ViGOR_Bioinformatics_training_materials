#!/usr/bin/env bash
#
# env_prep.sh — create the conda environment from the shared YAML
#
# Usage:
#   bash scripts/env_prep.sh
#   mamba activate metagenomics
# -----------------------------------------------------------------------------
set -euo pipefail
export PATH="$HOME/miniforge3/bin:$PATH"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if command -v mamba >/dev/null 2>&1; then CONDA=mamba; else CONDA=conda; fi
echo ">> Creating environment with: $CONDA"
"$CONDA" env create -f "$ROOT/envs/metagenomics.yml"

echo
echo ">> Done. Activate it with:"
echo "     mamba activate metagenomics"
