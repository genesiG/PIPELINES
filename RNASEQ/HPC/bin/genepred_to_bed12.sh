#!/bin/bash
# genepred_to_bed12.sh -- genePred to BED12 (UCSC genePredToBed), the gene model read by RSeQC infer_experiment.py.
# Called by GENEPRED_TO_BED12 (modules/local/gtf_to_bed12.nf), run once per pipeline run:
#   genepred_to_bed12.sh --genepred genes.genepred --out genes.bed12 [extra genePredToBed options...]
# Arguments:
#   --genepred  input genePred file (GTF_TO_GENEPRED output)
#   --out       output BED12 file (genes.bed12)
#   other arguments are placed before the file arguments of genePredToBed
# Command (from the spec):
#   genePredToBed genes.genepred genes.bed12
# Exit 3 when the genePred input is missing or empty, or when the BED12 output is empty.
# Outputs: genes.bed12 (or the --out name)
set -euo pipefail

genepred="" out="" extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --genepred) genepred="$2"; shift 2 ;;
        --out)      out="$2"; shift 2 ;;
        *)          extra+=("$1"); shift ;;
    esac
done
for v in genepred out; do
    [[ -n "${!v}" ]] || { echo "genepred_to_bed12.sh: missing --${v}" >&2; exit 3; }
done
[[ -s "${genepred}" ]] || { echo "genepred_to_bed12.sh: empty or missing genePred ${genepred}" >&2; exit 3; }

genePredToBed ${extra[@]+"${extra[@]}"} "${genepred}" "${out}"

[[ -s "${out}" ]] || { echo "genepred_to_bed12.sh: empty output ${out}" >&2; exit 3; }
