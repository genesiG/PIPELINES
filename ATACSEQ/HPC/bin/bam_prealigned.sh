#!/bin/bash
# bam_prealigned.sh -- bring a user-supplied pre-aligned BAM into the final-BAM layout.
# Called by BAM_PREALIGNED (modules/local/bam_prealigned.nf):
#   bam_prealigned.sh --sample S --bam input.bam --threads N
# Arguments: --sample, --bam, --threads; other arguments go to `samtools sort` (only if sorting is needed).
# The input is assumed final (filtered, duplicate-removed, blacklist-filtered: it gets no BAM_FILTER, DEDUP or BLACKLIST_FILTER step).
# The input is copied (cp -L) to S.final.bam; if its header is not SO:coordinate it is sorted.
# No raw stats (BAM_FILTER) and no orphan check (BAM_FINALIZE) are produced for pre-aligned samples.
# Outputs: S.final.bam, S.final.bam.bai, S.final.flagstat.txt, S.final.idxstats.txt
set -euo pipefail

sample="" bam="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)  sample="$2"; shift 2 ;;
        --bam)     bam="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        *)         extra+=("$1"); shift ;;
    esac
done
[[ -n "${sample}" && -s "${bam}" ]] || { echo "bam_prealigned.sh: missing --sample or empty --bam" >&2; exit 3; }

hd="$(samtools view -H "${bam}" | grep '^@HD' || true)"
if [[ "${hd}" == *SO:coordinate* ]]; then
    cp -L "${bam}" "${sample}.final.bam"
else
    echo "Input BAM is not coordinate-sorted; sorting" >&2
    samtools sort -m 1G -@ "${threads}" ${extra[@]+"${extra[@]}"} -o "${sample}.final.bam" "${bam}"
fi
samtools index -@ "${threads}" "${sample}.final.bam"
samtools flagstat -@ "${threads}" "${sample}.final.bam" > "${sample}.final.flagstat.txt"
samtools idxstats "${sample}.final.bam" > "${sample}.final.idxstats.txt"
