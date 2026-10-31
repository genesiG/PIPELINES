#!/bin/bash
# markdup.sh -- Picard MarkDuplicates (duplicates marked, not removed).
# Called by MARKDUP (modules/local/markdup.nf):
#   markdup.sh --sample S --bam S.filtered.bam --mem-gb G --threads N
# Arguments: --sample, --bam (coordinate-sorted), --mem-gb (task memory in GB; JVM -Xmx = max(1, G-2) g),
#   --threads (unused by Picard, accepted for a uniform interface); other arguments go to picard.
# Legacy parameters: OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 REMOVE_DUPLICATES=false.
# CREATE_INDEX=true writes S.markdup.bai, renamed to S.markdup.bam.bai (no samtools in the image).
# Outputs: S.markdup.bam, S.markdup.bam.bai, S.markdup_metrics.txt
set -euo pipefail

sample="" bam="" mem_gb=4 threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)  sample="$2"; shift 2 ;;
        --bam)     bam="$2"; shift 2 ;;
        --mem-gb)  mem_gb="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        *)         extra+=("$1"); shift ;;
    esac
done
[[ -n "${sample}" && -s "${bam}" ]] || { echo "markdup.sh: missing --sample or empty --bam" >&2; exit 3; }

xmx=$(( mem_gb - 2 )); (( xmx >= 1 )) || xmx=1
mkdir -p picard_tmp

picard -Xmx"${xmx}g" MarkDuplicates \
    INPUT="${bam}" \
    OUTPUT="${sample}.markdup.bam" \
    METRICS_FILE="${sample}.markdup_metrics.txt" \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
    REMOVE_DUPLICATES=false \
    ASSUME_SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=SILENT \
    CREATE_INDEX=true \
    TMP_DIR=picard_tmp \
    ${extra[@]+"${extra[@]}"}

mv "${sample}.markdup.bai" "${sample}.markdup.bam.bai"
rm -rf picard_tmp
