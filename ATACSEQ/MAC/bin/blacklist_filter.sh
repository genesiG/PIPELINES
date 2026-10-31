#!/bin/bash
# blacklist_filter.sh -- pair-aware blacklist removal (both mates dropped when either one overlaps the blacklist).
# Called by BLACKLIST_FILTER (modules/local/blacklist_filter.nf):
#   blacklist_filter.sh --sample S --bam S.dedup.bam --blacklist blacklist.bed
# Arguments:
#   --sample     sample id
#   --bam        dedup BAM (DEDUP output; any sort order on input)
#   --blacklist  BED file of blacklisted regions (config blacklist; required here, the process is skipped when empty)
#   other arguments are appended to bedtools pairtobed.
# Steps: samtools sort -n (name order, so the two mates of a pair are adjacent) into S.namesort.tmp.bam;
#   bedtools pairtobed -abam S.namesort.tmp.bam -b blacklist -type neither > S.noblacklist.namesorted.bam;
#   the temporary BAM is removed.
# `-type neither` keeps a pair only when NEITHER mate overlaps the blacklist, so BOTH mates are dropped when
# either overlaps (a plain `intersect -v` drops single reads and leaves orphan mates). MAPQ and sequences are
# kept: no BEDPE round trip (legacy bedpetobam set MAPQ 255 and dropped sequences).
# Output: S.noblacklist.namesorted.bam (name-sorted; BAM_FINALIZE coordinate-sorts it)
# Exit 3 when the input BAM or the blacklist is missing or empty.
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

namesorted="${sample}.namesort.tmp.bam"
samtools sort -n -m "${ATAC_SORT_MEM:-1G}" -o "${namesorted}" "${bam}"
bedtools pairtobed -abam "${namesorted}" -b "${blacklist}" -type neither ${extra[@]+"${extra[@]}"} \
    > "${sample}.noblacklist.namesorted.bam"
rm -f "${namesorted}"
