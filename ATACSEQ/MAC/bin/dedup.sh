#!/bin/bash
# dedup.sh -- sambamba markdup on the coordinate-sorted filtered BAM: duplicates removed, or only marked.
# Called by DEDUP (modules/local/dedup.nf):
#   dedup.sh --sample S --bam S.filtered.bam --remove true --threads N
# Arguments:
#   --sample   sample id
#   --bam      coordinate-sorted BAM (BAM_FILTER output)
#   --remove   true: sambamba markdup -r (config remove_duplicates, default true);
#              false: duplicates are marked but kept
#   --threads  sambamba -t (task.cpus)
#   other arguments are appended to sambamba markdup.
# Command: sambamba markdup [-r] -t N --tmpdir . in.bam S.dedup.bam 2> S.sambamba_markdup.log
# Legacy (step_3_process_bam.py): `sambamba markdup -r` when removing, `sambamba markdup` when marking.
# Exit 3 when the input BAM is missing or empty, or the output contains no reads.
# Outputs: S.dedup.bam (not published; BLACKLIST_FILTER reads it), S.sambamba_markdup.log (sambamba stderr)
set -euo pipefail

sample="" bam="" remove="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)  sample="$2"; shift 2 ;;
        --bam)     bam="$2"; shift 2 ;;
        --remove)  remove="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        *)         extra+=("$1"); shift ;;
    esac
done
for v in sample bam remove; do
    [[ -n "${!v}" ]] || { echo "dedup.sh: missing --${v}" >&2; exit 3; }
done
[[ -s "${bam}" ]] || { echo "dedup.sh: empty or missing input BAM ${bam}" >&2; exit 3; }

case "${remove}" in
    true)  rm_flag=(-r) ;;
    false) rm_flag=() ;;
    *) echo "dedup.sh: --remove must be true or false (got '${remove}')" >&2; exit 3 ;;
esac

sambamba markdup ${rm_flag[@]+"${rm_flag[@]}"} -t "${threads}" --tmpdir . ${extra[@]+"${extra[@]}"} \
    "${bam}" "${sample}.dedup.bam" 2> "${sample}.sambamba_markdup.log"
cat "${sample}.sambamba_markdup.log" >&2

reads=$(sambamba view -c "${sample}.dedup.bam")
[[ "${reads}" -gt 0 ]] || { echo "dedup.sh: no reads in ${sample}.dedup.bam after markdup" >&2; exit 3; }
