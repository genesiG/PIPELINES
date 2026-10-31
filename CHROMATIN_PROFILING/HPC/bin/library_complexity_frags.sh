#!/bin/bash
# library_complexity_frags.sh -- NRF / PBC1 / PBC2 from fragment coordinates of one pre-duplicate-removal BAM.
# Called by LIBCOMPLEXITY_FRAGS (modules/local/libcomplexity_frags.nf):
#   library_complexity_frags.sh --sample S --bam S.filtered.bam --threads N [--sort-mem 2G]
# Arguments:
#   --sample    sample id (written to the output)
#   --bam       coordinate-sorted, MAPQ/flag-filtered BAM in which duplicates are still present (S.filtered.bam)
#   --threads   threads for samtools and sort (default 1)
#   --sort-mem  memory per `sort` (default 2G)
# Method: read1 of every proper pair (samtools view -f 66 -F 2308 = read1 + proper pair, excluding unmapped,
#   secondary and supplementary). Fragment = chrom, min(pos, pnext), min(pos, pnext) + |tlen|. Each distinct
#   fragment coordinate is counted with `sort | uniq -c | awk`.
#   total_fragments    = number of fragments (read pairs)
#   distinct_fragments = number of distinct (chrom, start, end)
#   once / twice       = distinct fragments observed exactly once / exactly twice
#   NRF  = distinct / total ; PBC1 = once / distinct ; PBC2 = once / twice (NA when twice = 0)
# The BAM still holds chrM and spike-in contigs, so these fragments are counted too.
# Output: S.complexity.tsv  (header + one row: sample_id total_fragments distinct_fragments once twice NRF PBC1 PBC2)
# Exit 3 when the BAM contains no proper-pair read1 alignments.
set -euo pipefail

sample="" bam="" threads=1 sort_mem="2G"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)   sample="$2"; shift 2 ;;
        --bam)      bam="$2"; shift 2 ;;
        --threads)  threads="$2"; shift 2 ;;
        --sort-mem) sort_mem="$2"; shift 2 ;;
        *) echo "library_complexity_frags.sh: unknown argument $1" >&2; exit 2 ;;
    esac
done
[[ -n "${sample}" && -s "${bam}" ]] || { echo "library_complexity_frags.sh: missing --sample or empty --bam" >&2; exit 3; }

mkdir -p sort_tmp
# GNU sort options only where GNU sort exists (minimal images may ship busybox sort)
sort_opts=()
if sort --version 2>/dev/null | grep -q GNU; then
    sort_opts=(-T sort_tmp -S "${sort_mem}" --parallel="${threads}")
fi
samtools view -@ "${threads}" -f 66 -F 2308 "${bam}" \
  | awk 'BEGIN{OFS="\t"} {s = ($4 < $8) ? $4 : $8; l = ($9 < 0) ? -$9 : $9; print $3, s, s + l}' \
  | LC_ALL=C sort ${sort_opts[@]+"${sort_opts[@]}"} \
  | uniq -c \
  | awk -v s="${sample}" 'BEGIN{OFS="\t"}
      {tot += $1; dist++; if ($1 == 1) once++; if ($1 == 2) twice++}
      END{
        if (tot == 0) { print "no proper-pair read1 alignments in the BAM" > "/dev/stderr"; exit 3 }
        print "sample_id", "total_fragments", "distinct_fragments", "once", "twice", "NRF", "PBC1", "PBC2"
        print s, tot, dist, once + 0, twice + 0, sprintf("%.4f", dist / tot), sprintf("%.4f", (once + 0) / dist), (twice > 0 ? sprintf("%.4f", once / twice) : "NA")
      }' > "${sample}.complexity.tsv"
rm -rf sort_tmp
