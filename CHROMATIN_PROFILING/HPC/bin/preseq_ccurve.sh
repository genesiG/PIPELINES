#!/bin/bash
# preseq_ccurve.sh -- preseq c_curve (distinct fragments vs fragments sequenced) for one filtered BAM.
# Called by PRESEQ (modules/local/preseq.nf):
#   preseq_ccurve.sh --sample S --bam S.filtered.bam [--step N] [preseq args...]
# Arguments:
#   --sample  sample id
#   --bam     coordinate-sorted filtered BAM, duplicates still present (S.filtered.bam)
#   --step    c_curve step size in fragments. Default: 1,000,000; when that yields fewer than 10 curve points
#             (small libraries) it is re-run with step = max(1000, last_total / 20) so every curve has ~20 points.
#             (The legacy script took the step from samtools flagstat; the preseq image has no samtools.)
# Invocation retained from the legacy pipeline: `preseq c_curve -B -P -v -s STEP -o OUT BAM` (BAM input, paired-end).
# Output: S.preseq_ccurve.txt (preseq output: total_reads <tab> distinct_reads, one row per step)
# Exit 3 if the BAM is empty or preseq produces no curve.
set -euo pipefail

sample="" bam="" step=1000000 step_given=0 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample) sample="$2"; shift 2 ;;
        --bam)    bam="$2"; shift 2 ;;
        --step)   step="$2"; step_given=1; shift 2 ;;
        *)        extra+=("$1"); shift ;;
    esac
done
[[ -n "${sample}" && -s "${bam}" ]] || { echo "preseq_ccurve.sh: missing --sample or empty --bam" >&2; exit 3; }
out="${sample}.preseq_ccurve.txt"

# -v prints "TOTAL READS = N" on stderr; keep it so the retry step can be sized from the real depth
run() { preseq c_curve -B -P -v -s "$1" -o "${out}" ${extra[@]+"${extra[@]}"} "${bam}" 2> preseq.verbose.log; local rc=$?; cat preseq.verbose.log >&2; return $rc; }
npoints() { grep -c -E '^[0-9]' "${out}" || true; }

run "${step}"
if (( step_given == 0 )) && (( $(npoints) < 10 )); then
    # Size the step from the library depth preseq reported (a step larger than the
    # library gives a single "0 0" row, so the curve's own maximum cannot be used).
    last=$(awk -F'=' '/TOTAL READS/ {gsub(/ /, "", $2); print $2; exit}' preseq.verbose.log)
    [[ "${last}" =~ ^[0-9]+$ ]] && (( last > 0 )) || { echo "preseq_ccurve.sh: preseq reported no reads for ${bam}" >&2; exit 3; }
    newstep=$(( last / 20 )); (( newstep >= 1000 )) || newstep=1000
    echo "preseq_ccurve.sh: only $(npoints) points at step ${step}; re-running with step ${newstep}" >&2
    run "${newstep}"
fi
(( $(npoints) >= 2 )) || { echo "preseq_ccurve.sh: no curve produced for ${bam}" >&2; exit 3; }
