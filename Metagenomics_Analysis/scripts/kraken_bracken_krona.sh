
#!/usr/bin/env bash
#
# kraken_bracken_krona.sh
#
# Paired-end FASTQ pipeline:
#   Kraken 2 -> Bracken -> Krona
#
# Examples:
#   bash scripts/kraken_bracken_krona.sh
#   THREADS=4 bash scripts/kraken_bracken_krona.sh
#   SAMPLE=S29 bash scripts/kraken_bracken_krona.sh
#   SAMPLE=S29 FORCE=1 bash scripts/kraken_bracken_krona.sh
#   KRAKEN_DB=/path/to/another/db READ_LENGTH=150 bash scripts/kraken_bracken_krona.sh
#
# Set FORCE=1 to rerun all three tools for the selected samples.

set -Eeuo pipefail
shopt -s nullglob

# --------------------------- Configuration ---------------------------

PROJECT_DIR="${PROJECT_DIR:-$HOME/Documents/metagenomics}"
READ_DIR="${READ_DIR:-$PROJECT_DIR/results/qc/dehosted}"
KRAKEN_DB="${KRAKEN_DB:-$PROJECT_DIR/databases/kraken2_standard_08_20260626}"

THREADS="${THREADS:-8}"
READ_LENGTH="${READ_LENGTH:-300}"
BRACKEN_LEVEL="${BRACKEN_LEVEL:-S}"
BRACKEN_THRESHOLD="${BRACKEN_THRESHOLD:-0}"

SAMPLE="${SAMPLE:-}"
FORCE="${FORCE:-0}"

KRAKEN_CLASS_DIR="$PROJECT_DIR/results/kraken/classifications"
KRAKEN_REPORT_DIR="$PROJECT_DIR/results/kraken/reports"
KRAKEN_LOG_DIR="$PROJECT_DIR/results/kraken/logs"
BRACKEN_DIR="$PROJECT_DIR/results/bracken"
BRACKEN_LOG_DIR="$PROJECT_DIR/results/bracken/logs"
KRONA_DIR="$PROJECT_DIR/results/krona"
KRONA_LOG_DIR="$PROJECT_DIR/results/krona/logs"
KRONA_COMPARE_OUT="$KRONA_DIR/all_samples_comparison.krona.html"
KRONA_COMPARE_LOG="$KRONA_LOG_DIR/all_samples_comparison.log"
TAXPASTA_DIR="$PROJECT_DIR/results/taxpasta"
TAXPASTA_OUT="$TAXPASTA_DIR/bracken_merged.tsv"
TAXPASTA_LOG="$TAXPASTA_DIR/taxpasta_merge.log"

mkdir -p \
    "$KRAKEN_CLASS_DIR" "$KRAKEN_REPORT_DIR" "$KRAKEN_LOG_DIR" \
    "$BRACKEN_DIR" "$BRACKEN_LOG_DIR" \
    "$KRONA_DIR" "$KRONA_LOG_DIR" \
    "$TAXPASTA_DIR"

# ------------------------------ Helpers ------------------------------

log() {
    printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 ||
        die "Required program not found in PATH: $1"
}

# --------------------------- Pre-flight checks -----------------------

[[ -d "$PROJECT_DIR" ]] || die "Project directory not found: $PROJECT_DIR"
[[ -d "$READ_DIR" ]] || die "Input directory not found: $READ_DIR"

[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] ||
    die "THREADS must be a positive integer."
[[ "$READ_LENGTH" =~ ^[1-9][0-9]*$ ]] ||
    die "READ_LENGTH must be a positive integer."
[[ "$BRACKEN_THRESHOLD" =~ ^[0-9]+$ ]] ||
    die "BRACKEN_THRESHOLD must be a non-negative integer."
[[ "$FORCE" == 0 || "$FORCE" == 1 ]] ||
    die "FORCE must be 0 or 1."

for cmd in kraken2 bracken ktImportTaxonomy taxpasta awk; do
    require_command "$cmd"
done

for file in hash.k2d opts.k2d taxo.k2d; do
    [[ -s "$KRAKEN_DB/$file" ]] ||
        die "Missing Kraken database file: $KRAKEN_DB/$file"
done


[[ -s "$KRAKEN_DB/database${READ_LENGTH}mers.kmer_distrib" ]] ||
    die "Missing Bracken k-mer distribution. Run bracken-build."

log "Project: $PROJECT_DIR"
log "Database: $KRAKEN_DB"
log "Threads: $THREADS | Read length: $READ_LENGTH | Force: $FORCE"

# ----------------------------- Find samples --------------------------

R1_FILES=("$READ_DIR"/*_R1.clean.fastq.gz)
SAMPLES=()

for R1 in "${R1_FILES[@]}"; do
    BASE="$(basename "$R1")"
    ID="${BASE%_R1.clean.fastq.gz}"

    [[ -z "$SAMPLE" || "$ID" == "$SAMPLE" ]] || continue

    R2="$READ_DIR/${ID}_R2.clean.fastq.gz"
    [[ -s "$R2" ]] || die "Missing or empty R2 file for $ID: $R2"

    SAMPLES+=("$ID")
done

[[ ${#SAMPLES[@]} -gt 0 ]] ||
    die "No matching samples found in $READ_DIR. Check filenames or SAMPLE=$SAMPLE."

# ----------------------------- Main pipeline -------------------------

FAILED=0

for ID in "${SAMPLES[@]}"; do
    R1="$READ_DIR/${ID}_R1.clean.fastq.gz"
    R2="$READ_DIR/${ID}_R2.clean.fastq.gz"

    CLASS_OUT="$KRAKEN_CLASS_DIR/${ID}.kraken"
    REPORT="$KRAKEN_REPORT_DIR/${ID}.report"
    KRAKEN_LOG="$KRAKEN_LOG_DIR/${ID}.kraken.log"

    BRACKEN_OUT="$BRACKEN_DIR/${ID}.bracken"
    BRACKEN_LOG="$BRACKEN_LOG_DIR/${ID}.bracken.log"

    KRONA_INPUT="$KRONA_DIR/${ID}.krona_input.tsv"
    KRONA_HTML="$KRONA_DIR/${ID}.krona.html"
    KRONA_LOG="$KRONA_LOG_DIR/${ID}.krona.log"

    log "============================================================"
    log "Processing sample: $ID"

    KRAKEN_RERAN=0
    KRONA_INPUT_REFRESHED=0

    # 1. Kraken 2
    # Do not use --use-names: column 3 must contain numeric taxids.

    KRAKEN_VALID=0

    if [[ "$FORCE" == 0 && -s "$CLASS_OUT" && -s "$REPORT" ]]; then
        if awk -F '\t' '
            NF < 3 { bad=1; exit }
            $3 !~ /^[0-9]+$/ { bad=1; exit }
            END { exit bad }
        ' "$CLASS_OUT"; then
            KRAKEN_VALID=1
        fi
    fi

    if [[ "$KRAKEN_VALID" == 1 ]]; then
        log "Kraken outputs exist and have numeric taxids; skipping Kraken."
    else
        log "Running Kraken 2..."
        KRAKEN_RERAN=1

        rm -f "$CLASS_OUT" "$REPORT"

        if kraken2 \
            --db "$KRAKEN_DB" \
            --threads "$THREADS" \
            --paired \
            --gzip-compressed \
            --report "$REPORT" \
            --output "$CLASS_OUT" \
            "$R1" "$R2" > "$KRAKEN_LOG" 2>&1; then
            :
        else
            log "Kraken 2 failed. See $KRAKEN_LOG"
            rm -f "$CLASS_OUT" "$REPORT"
            FAILED=1
            continue
        fi

        if [[ ! -s "$CLASS_OUT" || ! -s "$REPORT" ]]; then
            log "Kraken did not produce both expected outputs for $ID."
            FAILED=1
            continue
        fi

        if ! awk -F '\t' '
            NF < 3 { bad=1; exit }
            $3 !~ /^[0-9]+$/ { bad=1; exit }
            END { exit bad }
        ' "$CLASS_OUT"; then
            log "Kraken output does not contain numeric taxids: $ID"
            FAILED=1
            continue
        fi

        log "Kraken 2 completed."
    fi

    # 2. Bracken
    # Rerun if forced, Kraken reran, or the output is missing.

    if [[ "$FORCE" == 1 || "$KRAKEN_RERAN" == 1 || ! -s "$BRACKEN_OUT" ]]; then
        log "Running Bracken..."

        rm -f "$BRACKEN_OUT"
        rm -f "$BRACKEN_LOG"

        if bracken \
            -d "$KRAKEN_DB" \
            -i "$REPORT" \
            -o "$BRACKEN_OUT" \
            -r "$READ_LENGTH" \
            -l "$BRACKEN_LEVEL" \
            -t "$BRACKEN_THRESHOLD" > "$BRACKEN_LOG" 2>&1; then
            :
        else
            log "Bracken failed. See $BRACKEN_LOG"
            rm -f "$BRACKEN_OUT"
            FAILED=1
            continue
        fi

        if [[ ! -s "$BRACKEN_OUT" ]]; then
            log "Bracken output is empty for $ID. See $BRACKEN_LOG"
            FAILED=1
            continue
        fi

        log "Bracken completed."
    else
        log "Bracken output exists; skipping Bracken."
    fi

    # 3. Krona input
    # Two columns: sequence/read ID and numeric taxonomy ID.
    # Refresh after Kraken reruns, or when forced/missing.

    if [[ "$FORCE" == 1 || "$KRAKEN_RERAN" == 1 || ! -s "$KRONA_INPUT" ]]; then
        log "Preparing Krona input..."

        rm -f "$KRONA_INPUT"

        if ! awk -F '\t' 'BEGIN { OFS="\t" }
            $1 == "C" && $3 ~ /^[0-9]+$/ && $3 > 0 {
                print $2, $3
            }
        ' "$CLASS_OUT" > "$KRONA_INPUT"; then
            log "Failed to create Krona input for $ID."
            FAILED=1
            continue
        fi

        if [[ ! -s "$KRONA_INPUT" ]]; then
            log "Krona input contains no classified reads for $ID."
            FAILED=1
            continue
        fi

        KRONA_INPUT_REFRESHED=1
    else
        log "Krona input exists; skipping input preparation."
    fi

    # 4. Krona HTML
    # Regenerate if forced, Kraken/input changed, or HTML is missing.

    if [[ "$FORCE" == 1 || "$KRAKEN_RERAN" == 1 ||
          "$KRONA_INPUT_REFRESHED" == 1 || ! -s "$KRONA_HTML" ]]; then
        log "Generating Krona HTML..."

        rm -f "$KRONA_HTML"

        if ktImportTaxonomy \
            -q 1 \
            -t 2 \
            "$KRONA_INPUT" \
            -o "$KRONA_HTML" > "$KRONA_LOG" 2>&1; then
            :
        else
            log "Krona failed. See $KRONA_LOG"
            rm -f "$KRONA_HTML"
            FAILED=1
            continue
        fi

        if [[ ! -s "$KRONA_HTML" ]]; then
            log "Krona HTML is empty for $ID."
            FAILED=1
            continue
        fi

        log "Krona completed."
    else
        log "Krona HTML exists; skipping Krona."
    fi

    log "Finished sample: $ID"
    log "Kraken report: $REPORT"
    log "Bracken table: $BRACKEN_OUT"
    log "Krona report:  $KRONA_HTML"
done


if [[ "$FAILED" -ne 0 ]]; then
    die "One or more samples failed. Check the per-sample logs."
fi


# ---------------------- Comparative Krona aggregation -----------------

log "Building comparative Krona chart..."

KRONA_INPUT_FILES=("$KRONA_DIR"/*.krona_input.tsv)

if [[ ${#KRONA_INPUT_FILES[@]} -eq 0 ]]; then
    die "No Krona input files found in $KRONA_DIR"
fi

# Each file becomes an independently selectable dataset.
# The file,name syntax labels each dataset with its sample ID.
KRONA_COMPARE_ARGS=()

for file in "${KRONA_INPUT_FILES[@]}"; do
    sample=$(basename "$file" .krona_input.tsv)
    KRONA_COMPARE_ARGS+=("${file},${sample}")
done

if ktImportTaxonomy \
    -q 1 \
    -t 2 \
    "${KRONA_COMPARE_ARGS[@]}" \
    -o "$KRONA_COMPARE_OUT" \
    > "$KRONA_COMPARE_LOG" 2>&1; then

    if [[ -s "$KRONA_COMPARE_OUT" ]]; then
        log "Comparative Krona chart created: $KRONA_COMPARE_OUT"
    else
        die "Comparative Krona output is missing or empty. Check $KRONA_COMPARE_LOG"
    fi
else
    log "Comparative Krona failed. See $KRONA_COMPARE_LOG"
    die "Comparative Krona aggregation failed."
fi

# ------------------------- Taxpasta aggregation -----------------------

mkdir -p "$TAXPASTA_DIR"

BRACKEN_FILES=("$BRACKEN_DIR"/*.bracken)

if [[ ${#BRACKEN_FILES[@]} -eq 0 ]]; then
    die "No Bracken profiles found in $BRACKEN_DIR"
fi

if [[ ${#BRACKEN_FILES[@]} -lt 2 ]]; then
    log "Skipping Taxpasta merge: found ${#BRACKEN_FILES[@]} Bracken profile; at least two are required."
else
    log "Merging ${#BRACKEN_FILES[@]} Bracken profiles with Taxpasta..."

if taxpasta merge \
    --profiler bracken \
    "${BRACKEN_FILES[@]}" \
    --output "$TAXPASTA_OUT" \
    > "$TAXPASTA_LOG" 2>&1; then

    [[ -s "$TAXPASTA_OUT" ]] ||
        die "Taxpasta completed but the merged table is missing or empty."

    log "Taxpasta merge completed."
    log "Merged table: $TAXPASTA_OUT"
    log "Taxpasta log: $TAXPASTA_LOG"
else
    log "Taxpasta merge failed. See $TAXPASTA_LOG"
    die "Taxpasta aggregation failed."
fi
fi

log "All requested pipeline stages completed successfully."