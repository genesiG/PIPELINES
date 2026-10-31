#!/bin/bash
# seacr_bedgraph.sh - fragment-level bedGraph of a final paired-end BAM for SEACR.
# Called by SEACR_BEDGRAPH (modules/local/seacr_bedgraph.nf):
#   seacr_bedgraph.sh --sample S --bam S.final.bam
# Arguments:
#   --sample  sample id
#   --bam     final BAM (coordinate-sorted, mate-consistent: BAM_FINALIZE orphan check)
# Command (as legacy step_6b): bedtools genomecov -ibam BAM -bg -pc | LC_ALL=C sort -k1,1 -k2,2n
# Note: -pc builds each fragment from the mate position / template length stored in the
# record, so no name sorting is needed and only bedtools is used (the image has no samtools).
# Legacy passed -g chrom.sizes next to -ibam; bedtools ignores -g for BAM input (header is used).
# Output: <sample>.fragments.bedgraph (sorted)
set -euo pipefail

sample="" bam="" extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample) sample="$2"; shift 2 ;;
        --bam)    bam="$2"; shift 2 ;;
        *)        extra+=("$1"); shift ;;
    esac
done
[[ -n "$sample" && -n "$bam" ]] || { echo "seacr_bedgraph.sh: --sample and --bam are required" >&2; exit 2; }

bedtools genomecov -ibam "$bam" -bg -pc ${extra[@]+"${extra[@]}"} | LC_ALL=C sort -k1,1 -k2,2n > "${sample}.fragments.bedgraph"
