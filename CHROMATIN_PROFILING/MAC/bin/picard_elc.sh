#!/bin/bash
# picard_elc.sh -- Picard EstimateLibraryComplexity on the pre-duplicate-removal filtered BAM.
# Called by PICARD_ELC (modules/local/picard_elc.nf):
#   picard_elc.sh --sample S --bam S.filtered.bam --mem-gb G --threads N [picard args...]
# Arguments: --sample, --bam (coordinate-sorted, duplicates still present), --mem-gb (task memory in GB; JVM -Xmx =
#   max(1, G-2) g), --threads (accepted for a uniform interface, unused by Picard). Other arguments go to Picard.
# Output: S.elc.txt (Picard metrics: LIBRARY, READ_PAIRS_EXAMINED, UNIQUE_READ_PAIRS, READ_PAIR_DUPLICATES,
#   PERCENT_DUPLICATION, ESTIMATED_LIBRARY_SIZE, ...)
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
[[ -n "${sample}" && -s "${bam}" ]] || { echo "picard_elc.sh: missing --sample or empty --bam" >&2; exit 3; }

xmx=$(( mem_gb - 2 )); (( xmx >= 1 )) || xmx=1
mkdir -p picard_tmp
picard -Xmx"${xmx}g" EstimateLibraryComplexity \
    INPUT="${bam}" \
    OUTPUT="${sample}.elc.txt" \
    TMP_DIR=picard_tmp \
    VALIDATION_STRINGENCY=SILENT \
    ${extra[@]+"${extra[@]}"}
rm -rf picard_tmp
