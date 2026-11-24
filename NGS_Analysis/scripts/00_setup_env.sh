#!/usr/bin/env bash
#
# 00_setup_env.sh — create the conda environment from the shared YAML
# -----------------------------------------------------------------------------
# Run once. Creates the "ebola-training" environment with every tool the
# workflow needs (mafft, iqtree, seqkit, plus the read-mapping tools).
#
# Usage:
#   bash scripts/00_setup_env.sh
#   conda activate ebola-training
# -----------------------------------------------------------------------------
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if command -v mamba >/dev/null 2>&1; then CONDA=mamba; else CONDA=conda; fi
echo ">> Creating environment with: $CONDA"
"$CONDA" env create -f "$ROOT/env/ebola-training.yml"

echo
echo ">> Done. Activate it with:"
echo "     conda activate ebola-training"
