#!/bin/bash
# featurecounts.sh -- featureCounts (subread) gene counts for one final BAM, paired-end read-pair counting.
# Called by FEATURECOUNTS (modules/local/featurecounts.nf):
#   featurecounts.sh --sample S --bam S.final.bam --gtf GTF --strand 0 --feature-type exon --attribute gene_id \
#       --threads N [extra featureCounts args...]
# Arguments:
#   --sample        sample id (output prefix)
#   --bam           final BAM (MARK_DUPLICATES output; all reads counted, duplicates are marked, not removed by default)
#   --gtf           annotation GTF (config gtf)
#   --strand        featureCounts -s: 0 unstranded, 1 forward, 2 reverse (main.nf maps config strandedness to this code)
#   --feature-type  featureCounts -t (config featurecounts_feature_type, default exon)
#   --attribute     featureCounts -g (config featurecounts_attribute, default gene_id)
#   --threads       featureCounts -T (task.cpus)
#   other arguments are placed before -a (config featurecounts_args, default empty)
# Command (spec; legacy step_4_count_reads.py had no -T and no -s, i.e. -s 0):
#   featureCounts -T N -p --countReadPairs -B -C --primary -F GTF -t TYPE -g ATTR -s CODE [args] -a GTF -o S_featurecounts.txt S.final.bam
#   tail -n+3 S_featurecounts.txt | cut -f1,7 > S_counts.txt
# The summary S_featurecounts.txt.summary is written by featureCounts.
# Exit 3 when an input is missing or empty, --strand is not 0, 1 or 2, or S_counts.txt is empty.
# Outputs: S_featurecounts.txt, S_featurecounts.txt.summary, S_counts.txt
set -euo pipefail

sample="" bam="" gtf="" strand="" feature_type="" attribute="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)        sample="$2"; shift 2 ;;
        --bam)           bam="$2"; shift 2 ;;
        --gtf)           gtf="$2"; shift 2 ;;
        --strand)        strand="$2"; shift 2 ;;
        --feature-type)  feature_type="$2"; shift 2 ;;
        --attribute)     attribute="$2"; shift 2 ;;
        --threads)       threads="$2"; shift 2 ;;
        *)               extra+=("$1"); shift ;;
    esac
done
for v in sample bam gtf strand feature_type attribute; do
    [[ -n "${!v}" ]] || { echo "featurecounts.sh: missing --${v//_/-}" >&2; exit 3; }
done
case "${strand}" in
    0|1|2) ;;
    *) echo "featurecounts.sh: --strand must be 0, 1 or 2 (got '${strand}')" >&2; exit 3 ;;
esac
[[ -s "${bam}" ]] || { echo "featurecounts.sh: empty or missing BAM ${bam}" >&2; exit 3; }
[[ -s "${gtf}" ]] || { echo "featurecounts.sh: empty or missing GTF ${gtf}" >&2; exit 3; }

featureCounts -T "${threads}" -p --countReadPairs -B -C --primary \
    -F GTF -t "${feature_type}" -g "${attribute}" -s "${strand}" \
    ${extra[@]+"${extra[@]}"} \
    -a "${gtf}" -o "${sample}_featurecounts.txt" "${bam}"

tail -n+3 "${sample}_featurecounts.txt" | cut -f1,7 > "${sample}_counts.txt"

[[ -s "${sample}_counts.txt" ]] || { echo "featurecounts.sh: empty count table ${sample}_counts.txt" >&2; exit 3; }
