#!/usr/bin/env bash
#
# 01_align.sh — multiple sequence alignment with MAFFT
# -----------------------------------------------------------------------------
# Aligns both datasets (genotype = cross-species panel, bdbv = Bundibugyo
# outbreak panel), reading from phylo/unaligned/ and writing to phylo/aligned/.
#
# Prereq:  conda activate ebola-training
# Usage:   bash scripts/01_align.sh          (or: THREADS=8 bash scripts/01_align.sh)
# -----------------------------------------------------------------------------
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
THREADS="${THREADS:-4}"

mkdir -p phylo/aligned
for DS in genotype bdbv; do
    echo ">> Aligning ${DS} with MAFFT (--auto)"
    mafft --auto --thread "$THREADS" \
        "phylo/unaligned/${DS}.fasta" > "phylo/aligned/${DS}.aln.fasta"
    seqkit stats "phylo/aligned/${DS}.aln.fasta"
done
echo ">> Alignments written to phylo/aligned/"
