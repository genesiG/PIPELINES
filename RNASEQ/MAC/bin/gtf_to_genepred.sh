#!/bin/bash
# gtf_to_genepred.sh -- GTF to genePred (UCSC gtfToGenePred), the gene model for RSeQC.
# Called by GTF_TO_GENEPRED (modules/local/gtf_to_bed12.nf), run once per pipeline run:
#   gtf_to_genepred.sh --gtf GTF --out genes.genepred [extra gtfToGenePred options...]
# Arguments:
#   --gtf     annotation GTF (config gtf)
#   --out     output genePred file (genes.genepred)
#   other arguments are placed before the file arguments of gtfToGenePred
# Command (legacy-equivalent, from the spec):
#   gtfToGenePred -ignoreGroupsWithoutExons GTF genes.genepred
# Exit 3 when the GTF is missing or empty, or when the genePred output is empty.
# Outputs: genes.genepred (or the --out name)
set -euo pipefail

gtf="" out="" extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --gtf)  gtf="$2"; shift 2 ;;
        --out)  out="$2"; shift 2 ;;
        *)      extra+=("$1"); shift ;;
    esac
done
for v in gtf out; do
    [[ -n "${!v}" ]] || { echo "gtf_to_genepred.sh: missing --${v}" >&2; exit 3; }
done
[[ -s "${gtf}" ]] || { echo "gtf_to_genepred.sh: empty or missing GTF ${gtf}" >&2; exit 3; }

gtfToGenePred -ignoreGroupsWithoutExons ${extra[@]+"${extra[@]}"} "${gtf}" "${out}"

[[ -s "${out}" ]] || { echo "gtf_to_genepred.sh: empty output ${out}" >&2; exit 3; }
