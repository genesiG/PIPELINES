#!/bin/bash
# blacklist_filter.sh -- pair-aware ENCODE blacklist removal.
# Called by BLACKLIST_FILTER (modules/local/blacklist_filter.nf):
#   blacklist_filter.sh --sample S --bam S.target.namesorted.bam --blacklist blacklist.bed
# Arguments: --sample, --bam (NAME-sorted, so mates are adjacent), --blacklist (BED);
#   other arguments are appended to bedtools pairtobed.
# `bedtools pairtobed -abam in -b bl.bed -type neither` keeps a pair only when NEITHER mate overlaps
# the blacklist: BOTH mates are dropped when either one overlaps (the legacy `intersect -v` dropped
# single reads and left orphan mates). With -abam, pairtobed reads the BAM directly (no BEDPE
# conversion needed) and writes a BAM, which stays name-sorted.
# Output: S.noblacklist.namesorted.bam
# Exit 3 when the blacklist is missing or empty.
set -euo pipefail

sample="" bam="" blacklist="" extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)    sample="$2"; shift 2 ;;
        --bam)       bam="$2"; shift 2 ;;
        --blacklist) blacklist="$2"; shift 2 ;;
        *)           extra+=("$1"); shift ;;
    esac
done
[[ -n "${sample}" && -s "${bam}" ]] || { echo "blacklist_filter.sh: missing --sample or empty --bam" >&2; exit 3; }
[[ -n "${blacklist}" && -s "${blacklist}" ]] || { echo "blacklist_filter.sh: blacklist missing or empty: '${blacklist}'" >&2; exit 3; }

bedtools pairtobed -abam "${bam}" -b "${blacklist}" -type neither ${extra[@]+"${extra[@]}"} \
    > "${sample}.noblacklist.namesorted.bam"
