#!/bin/bash
# qc_fingerprint.sh -- deepTools plotFingerprint for the BAMs of one group.
# Called by FINGERPRINT (modules/local/fingerprint.nf):
#   qc_fingerprint.sh --group G --bams a.final.bam,b.final.bam --labels a,b --threads N [plotFingerprint args...]
# Arguments:
#   --group    group name (used in the plot title and the output file names)
#   --bams     comma-separated list of indexed final BAMs (the .bai files must sit next to them)
#   --labels   comma-separated labels, same number and order as --bams (the sample ids)
#   --threads  processors
#   Other arguments go to plotFingerprint.
# Settings (legacy): --skipZeros --numberOfSamples 500000 --binSize 500; new: --extendReads (paired-end fragments
#   instead of single reads). Output format SVG.
# Outputs: fingerprint_<group>.svg (plot), fingerprint_<group>.tsv (raw per-bin counts, --outRawCounts)
# Exit 3 on a label/BAM count mismatch or fewer than 1 BAM.
set -euo pipefail

group="" bams="" labels="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --group)   group="$2"; shift 2 ;;
        --bams)    bams="$2"; shift 2 ;;
        --labels)  labels="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        *)         extra+=("$1"); shift ;;
    esac
done
[[ -n "${group}" && -n "${bams}" && -n "${labels}" ]] || { echo "qc_fingerprint.sh: need --group, --bams and --labels" >&2; exit 3; }
IFS=',' read -r -a bam_arr <<< "${bams}"
IFS=',' read -r -a label_arr <<< "${labels}"
if (( ${#bam_arr[@]} != ${#label_arr[@]} )); then
    echo "qc_fingerprint.sh: ${#bam_arr[@]} BAMs but ${#label_arr[@]} labels" >&2; exit 3
fi
for i in "${!bam_arr[@]}"; do
    [[ "${bam_arr[$i]}" == "${label_arr[$i]}."* ]] || echo "qc_fingerprint.sh: WARNING label '${label_arr[$i]}' does not prefix BAM '${bam_arr[$i]}' (check the order)" >&2
done

plotFingerprint \
    --bamfiles "${bam_arr[@]}" \
    --labels "${label_arr[@]}" \
    --extendReads \
    --skipZeros \
    --numberOfSamples 500000 \
    --binSize 500 \
    --numberOfProcessors "${threads}" \
    --plotFileFormat svg \
    --plotFile "fingerprint_${group}.svg" \
    --outRawCounts "fingerprint_${group}.tsv" \
    --plotTitle "${group} fingerprint" \
    ${extra[@]+"${extra[@]}"}
