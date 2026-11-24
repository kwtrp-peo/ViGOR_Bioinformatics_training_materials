#!/usr/bin/env bash
#
# qc.sh — paired-end read QC, trimming and human/PhiX dehosting
#
# Usage:
#   mamba activate metagenomics
#   bash scripts/qc.sh
#
set -Eeuo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

THREADS="${THREADS:-8}"
INPUT_DIR="${INPUT_DIR:-$ROOT/data/raw}"
REF_FASTA="${REF_FASTA:-$ROOT/references/human_phix.fna}"
OUT="$ROOT/results/qc"

# Output directories
RAW_QC="$OUT/fastqc_raw"
TRIM="$OUT/trimmed"
TRIM_QC="$OUT/fastqc_trimmed"
DEHOST="$OUT/dehosted"
REPORTS="$OUT/reports"
LOGS="$OUT/logs"
REF_INDEX="$ROOT/databases/bowtie2_human_phix/human_phix"

mkdir -p "$RAW_QC" "$TRIM" "$TRIM_QC" \
         "$DEHOST" "$REPORTS" "$LOGS" \
         "$(dirname "$REF_INDEX")"

# Check required programs
for tool in fastqc fastp multiqc bowtie2 bowtie2-build samtools; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: Required tool not found: $tool"
        exit 1
    fi
done

# Check input and reference files
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

if [[ ! -s "$REF_FASTA" ]]; then
    echo "ERROR: Combined human + PhiX FASTA not found: $REF_FASTA"
    echo "Provide a FASTA containing both reference sequences."
    exit 1
fi

# Build Bowtie2 index if not already available
if [[ ! -s "${REF_INDEX}.1.bt2" &&
      ! -s "${REF_INDEX}.1.bt2l" ]]; then
    echo ">> Building human + PhiX Bowtie2 index"
    bowtie2-build --threads "$THREADS" \
        "$REF_FASTA" "$REF_INDEX" \
        > "$LOGS/bowtie2_build.log" 2>&1
fi

# Discover paired-end R1 FASTQs
shopt -s nullglob
R1_FILES=("$INPUT_DIR"/*_R1*.fastq.gz)

if [[ ${#R1_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No files matching *_R1*.fastq.gz in $INPUT_DIR"
    exit 1
fi

echo ">> Found ${#R1_FILES[@]} R1 FASTQ files"
echo ">> Input directory: $INPUT_DIR"
echo ">> Threads: $THREADS"

for R1 in "${R1_FILES[@]}"; do

    BASE="$(basename "$R1")"
    SAMPLE="${BASE%%_R1*}"
    R2="${R1/_R1/_R2}"

    if [[ ! -s "$R2" ]]; then
        echo "ERROR: Matching R2 file not found for $R1"
        exit 1
    fi

    echo
    echo "=================================================="
    echo ">> Processing sample: $SAMPLE"
    echo "=================================================="

    # 1. FastQC on raw reads
    echo ">> [1/5] FastQC on raw reads"

    fastqc --threads "$THREADS" \
        --outdir "$RAW_QC" \
        "$R1" "$R2" \
        > "$LOGS/${SAMPLE}.fastqc_raw.log" 2>&1

    # 2. Adapter trimming and quality filtering
    echo ">> [2/5] Trimming with fastp"

    TRIM_R1="$TRIM/${SAMPLE}_R1.trimmed.fastq.gz"
    TRIM_R2="$TRIM/${SAMPLE}_R2.trimmed.fastq.gz"

    fastp \
        --in1 "$R1" \
        --in2 "$R2" \
        --out1 "$TRIM_R1" \
        --out2 "$TRIM_R2" \
        --thread "$THREADS" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --length_required 50 \
        --html "$REPORTS/${SAMPLE}.fastp.html" \
        --json "$REPORTS/${SAMPLE}.fastp.json" \
        2> "$LOGS/${SAMPLE}.fastp.log"

    # 3. FastQC after trimming
    echo ">> [3/5] FastQC on trimmed reads"

    fastqc --threads "$THREADS" \
        --outdir "$TRIM_QC" \
        "$TRIM_R1" "$TRIM_R2" \
        > "$LOGS/${SAMPLE}.fastqc_trimmed.log" 2>&1

    # 4. Dehost: retain only pairs where BOTH reads are
    # unmapped to the combined human + PhiX reference.
    echo ">> [4/5] Removing human/PhiX-mapped read pairs"

    OUT_R1="$DEHOST/${SAMPLE}_R1.clean.fastq"
    OUT_R2="$DEHOST/${SAMPLE}_R2.clean.fastq"

    bowtie2 \
        --very-sensitive \
        --threads "$THREADS" \
        -x "$REF_INDEX" \
        -1 "$TRIM_R1" \
        -2 "$TRIM_R2" \
        2> "$LOGS/${SAMPLE}.dehost_bowtie2.log" \
    | samtools view \
        -@ "$THREADS" \
        -b -f 12 -F 2304 \
    | samtools fastq \
        -@ "$THREADS" \
        -n -f 12 \
        -1 "$OUT_R1" \
        -2 "$OUT_R2" \
        -0 /dev/null \
        -s /dev/null \
        - \
        2> "$LOGS/${SAMPLE}.dehost_fastq.log"

    # Compress the dehosted FASTQ files
      gzip "$OUT_R1" "$OUT_R2"
      OUT_R1="${OUT_R1}.gz"
      OUT_R2="${OUT_R2}.gz"

    # 5. Report read counts and final output files
    echo ">> [5/5] Counting retained reads"

    {
        echo "Sample: $SAMPLE"
        echo "Trimmed R1:"
        gzip -cd "$TRIM_R1" | awk 'END {print "Reads:", NR/4}'
        echo "Trimmed R2:"
        gzip -cd "$TRIM_R2" | awk 'END {print "Reads:", NR/4}'
        echo "Dehosted R1:"
        gzip -cd "$OUT_R1" | awk 'END {print "Reads:", NR/4}'
        echo "Dehosted R2:"
        gzip -cd "$OUT_R2" | awk 'END {print "Reads:", NR/4}'
    } > "$REPORTS/${SAMPLE}.read_counts.txt"

    echo ">> Completed: $SAMPLE"
done

# Summarise raw FastQC, trimmed FastQC and fastp reports
echo
echo ">> Generating MultiQC reports"

multiqc "$RAW_QC" \
    --outdir "$REPORTS" \
    --filename multiqc_raw.html \
    --force

multiqc "$TRIM_QC" "$REPORTS" \
    --outdir "$REPORTS" \
    --filename multiqc_processed.html \
    --ignore "*multiqc*" \
    --force

echo
echo ">> All samples processed."
echo ">> Trimmed reads:  $TRIM"
echo ">> Dehosted reads: $DEHOST"
echo ">> Reports:        $REPORTS"
