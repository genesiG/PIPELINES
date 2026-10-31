#!/bin/bash
# bam_filter.sh -- sort, stats of the unfiltered alignment, MAPQ + flag filter.
# Called by BAM_FILTER (modules/local/bam_filter.nf):
#   bam_filter.sh --sample S --bam S.aligned.bam --mapq 30 --threads N
# Arguments: --sample, --bam (unsorted aligned BAM), --mapq (samtools -q), --threads;
#   other arguments are appended to the filtering `samtools view`.
# Steps: coordinate sort + index; flagstat/idxstats of that unfiltered sorted BAM;
#   samtools view -b -q MAPQ -f 3 -F 2828
# Flags copied from legacy step_3_process_bam.py (paired-end):
#   "-f 1 -f 2 -F 4 -F 8 -F 256 -F 512 -F 1024 -F 2048". samtools ORs repeated options, so
#   -f 1 -f 2 = -f 3 (paired + proper pair); -F 4,8,256,512,2048 = -F 2828.
#   The legacy -F 1024 is omitted: no duplicates are marked at this stage (MARKDUP follows),
#   so it was a no-op there. chrM and spike contigs are kept; duplicates are NOT removed.
# Outputs: S.filtered.bam, S.filtered.bam.bai, S.raw.flagstat.txt, S.raw.idxstats.txt
set -euo pipefail

sample="" bam="" mapq="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)  sample="$2"; shift 2 ;;
        --bam)     bam="$2"; shift 2 ;;
        --mapq)    mapq="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        *)         extra+=("$1"); shift ;;
    esac
done
for v in sample bam mapq; do
    [[ -n "${!v}" ]] || { echo "bam_filter.sh: missing --${v}" >&2; exit 3; }
done
[[ -s "${bam}" ]] || { echo "bam_filter.sh: empty or missing input BAM ${bam}" >&2; exit 3; }

sorted="${sample}.sorted_unfiltered.tmp.bam"
samtools sort -m 1G -@ "${threads}" -o "${sorted}" "${bam}"
samtools index -@ "${threads}" "${sorted}"
samtools flagstat -@ "${threads}" "${sorted}" > "${sample}.raw.flagstat.txt"
samtools idxstats "${sorted}" > "${sample}.raw.idxstats.txt"

samtools view -b -@ "${threads}" -q "${mapq}" -f 3 -F 2828 ${extra[@]+"${extra[@]}"} \
    -o "${sample}.filtered.bam" "${sorted}"
samtools index -@ "${threads}" "${sample}.filtered.bam"
rm -f "${sorted}" "${sorted}.bai"
