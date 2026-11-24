#!/usr/bin/env bash
#
# 02_tree.sh — maximum-likelihood trees with IQ-TREE
# -----------------------------------------------------------------------------
# Builds an ML tree for each aligned dataset:
#   -m MFP     ModelFinder Plus picks the best substitution model
#   -B 1000    1000 ultrafast bootstrap replicates (branch support)
#   -T AUTO    auto-detect thread count
# Outputs land in results/trees/<dataset>/ with prefix <dataset>.
#
# Prereq:  conda activate ebola-training  (after scripts/01_align.sh)
# Usage:   bash scripts/02_tree.sh
# -----------------------------------------------------------------------------
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

# IQ-TREE ships under different binary names across versions; pick what's present.
if   command -v iqtree3 >/dev/null 2>&1; then IQ=iqtree3
elif command -v iqtree2 >/dev/null 2>&1; then IQ=iqtree2
elif command -v iqtree  >/dev/null 2>&1; then IQ=iqtree
else echo "ERROR: no iqtree binary on PATH (activate the env first)"; exit 1; fi
echo ">> Using $IQ"

for DS in genotype bdbv; do
    echo ">> Building tree for ${DS}"
    mkdir -p "results/trees/${DS}"
    "$IQ" -s "phylo/aligned/${DS}.aln.fasta" \
          -m MFP -B 1000 -T AUTO \
          --prefix "results/trees/${DS}/${DS}" -redo
done

echo
echo ">> Trees written. Key files per dataset:"
echo "     results/trees/<ds>/<ds>.treefile     ML tree (Newick, with UFBoot support)"
echo "     results/trees/<ds>/<ds>.iqtree       full report (model, likelihoods)"
echo ">> Next: root & visualise — see README (iTOL section)."
