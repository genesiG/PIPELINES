#!/bin/bash
# count_pairs.sh -- number of read pairs of one FASTQ pair, for the chunk list of ALIGN (align_chunk_pairs > 0).
# Called by COUNT_PAIRS (modules/local/count_pairs.nf):
#   count_pairs.sh --sample S --r1 R1 --r2 R2
# Arguments: --sample, --r1, --r2 (FASTQ, gzipped or plain; the reads ALIGN reads, i.e. trimmed when trim is on).
#   Other arguments are collected and ignored.
# Output: S.pair_count.txt, one line with one integer: the number of pairs (records of R1). Nothing else.
# Exit 3: missing arguments or files; a line count not divisible by 4 (the file is named); R1 and R2 with different
#   numbers of records; zero pairs; a FASTQ that cannot be read.
set -euo pipefail

sample="" r1="" r2="" extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample) sample="$2"; shift 2 ;;
        --r1)     r1="$2"; shift 2 ;;
        --r2)     r2="$2"; shift 2 ;;
        *)        extra+=("$1"); shift ;;
    esac
done

[[ -n "${sample}" && -n "${r1}" && -n "${r2}" ]] || { echo "count_pairs.sh: missing --sample, --r1 or --r2" >&2; exit 3; }
for f in "${r1}" "${r2}"; do
    [[ -f "${f}" ]] || { echo "count_pairs.sh: FASTQ not found: ${f}" >&2; exit 3; }
done

# Lines (NR) of a FASTQ, gzipped or plain; zcat -f passes plain files through.
n1=$(zcat -f "${r1}" | awk 'END{print NR}') || { echo "count_pairs.sh: cannot read ${r1}" >&2; exit 3; }
n2=$(zcat -f "${r2}" | awk 'END{print NR}') || { echo "count_pairs.sh: cannot read ${r2}" >&2; exit 3; }
# a FASTQ record is 4 lines
if (( n1 % 4 != 0 )); then
    echo "count_pairs.sh: ${r1} has ${n1} lines, not a multiple of 4 (not a FASTQ file?)" >&2
    exit 3
fi
if (( n2 % 4 != 0 )); then
    echo "count_pairs.sh: ${r2} has ${n2} lines, not a multiple of 4 (not a FASTQ file?)" >&2
    exit 3
fi
if (( n1 != n2 )); then
    echo "count_pairs.sh: R1 has $((n1 / 4)) records and R2 has $((n2 / 4)) (${r1} vs ${r2})" >&2
    exit 3
fi
if (( n1 == 0 )); then
    echo "count_pairs.sh: no read pairs in ${r1}" >&2
    exit 3
fi

echo $(( n1 / 4 )) > "${sample}.pair_count.txt"
