#!/bin/bash
# infer_strandedness.sh -- RSeQC infer_experiment.py: fraction of reads per strand-class, one BAM per task.
# Called by INFER_STRANDEDNESS (modules/local/infer_strandedness.nf):
#   infer_strandedness.sh --sample S --bam S.final.bam --bed genes.bed12 --reads 200000 [extra infer_experiment args...]
# Arguments:
#   --sample  sample id (output prefix)
#   --bam     final BAM (its .bai must sit next to it)
#   --bed     gene model, BED12 (GENEPRED_TO_BED12 output)
#   --reads   infer_experiment -s, number of reads sampled (config infer_experiment_reads, default 200000)
#   other arguments are appended to infer_experiment.py
# Command (the strandedness check; not in the legacy pipeline, which never checked it):
#   infer_experiment.py -r genes.bed12 -i S.final.bam -s N > S.infer_experiment.txt
# infer_experiment.py stderr goes to S.infer_experiment.log and is repeated on the task stderr.
# Result lines in S.infer_experiment.txt (paired-end data):
#   This is PairEnd Data
#   Fraction of reads failed to determine: F_undetermined
#   Fraction of reads explained by "1++,1--,2+-,2-+": F_forward   -> forward: read 1 is on the same strand as the gene
#   Fraction of reads explained by "1+-,1-+,2++,2--": F_reverse   -> reverse: read 1 is on the opposite strand
#   (the reverse class is the usual one for dUTP-based stranded libraries; an unstranded library gives both near 0.5)
# A single-end layout ("++,--") is not expected for the lab's paired-end data; qc_strandedness_table.py reads only the
#   paired-end classes above.
# Exit 3 when an input is missing or empty, or when the output lacks the line "Fraction of reads failed to determine".
# Outputs: S.infer_experiment.txt, S.infer_experiment.log
set -euo pipefail

sample="" bam="" bed="" reads="" extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample) sample="$2"; shift 2 ;;
        --bam)    bam="$2"; shift 2 ;;
        --bed)    bed="$2"; shift 2 ;;
        --reads)  reads="$2"; shift 2 ;;
        *)        extra+=("$1"); shift ;;
    esac
done
for v in sample bam bed reads; do
    [[ -n "${!v}" ]] || { echo "infer_strandedness.sh: missing --${v}" >&2; exit 3; }
done
[[ -s "${bam}" ]] || { echo "infer_strandedness.sh: empty or missing BAM ${bam}" >&2; exit 3; }
[[ -s "${bed}" ]] || { echo "infer_strandedness.sh: empty or missing BED12 ${bed}" >&2; exit 3; }
[[ "${reads}" =~ ^[1-9][0-9]*$ ]] || { echo "infer_strandedness.sh: --reads must be a positive integer (got '${reads}')" >&2; exit 3; }

out="${sample}.infer_experiment.txt"
log="${sample}.infer_experiment.log"
if ! infer_experiment.py -r "${bed}" -i "${bam}" -s "${reads}" ${extra[@]+"${extra[@]}"} > "${out}" 2> "${log}"; then
    cat "${log}" >&2
    echo "infer_strandedness.sh: infer_experiment.py failed for ${sample}" >&2
    exit 1
fi
cat "${log}" >&2

grep -qF "Fraction of reads failed to determine" "${out}" \
    || { echo "infer_strandedness.sh: ${out} has no 'Fraction of reads failed to determine' line" >&2; exit 3; }
