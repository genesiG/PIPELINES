#!/bin/bash
# qc_correlation.sh -- genome-wide read-count correlation between all samples (deepTools).
# Called by CORRELATION (modules/local/correlation.nf):
#   qc_correlation.sh --bams a.final.bam,b.final.bam,... --labels a,b,... --threads N [multiBamSummary args...]
# Arguments:
#   --bams      comma-separated indexed final BAMs (.bai next to them)
#   --labels    comma-separated labels, same number and order as --bams (explicit; never parsed from file names)
#   --threads   processors
#   --bin-size  bin size in bp (default 5000, legacy)
#   Other arguments go to multiBamSummary.
# Steps (legacy settings): multiBamSummary bins --binSize 5000 --distanceBetweenBins 0 --extendReads (paired-end
#   fragments, new), then plotCorrelation --corMethod spearman --whatToPlot heatmap --skipZeros --removeOutliers
#   --colorMap RdYlBu --plotNumbers.
# Outputs: correlation_matrix.npz, correlation_heatmap.svg, correlation_matrix.tsv (Spearman matrix)
# Exit 3 on a label/BAM count mismatch or fewer than 2 BAMs.
set -euo pipefail

bams="" labels="" threads=1 bin_size=5000 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bams)     bams="$2"; shift 2 ;;
        --labels)   labels="$2"; shift 2 ;;
        --threads)  threads="$2"; shift 2 ;;
        --bin-size) bin_size="$2"; shift 2 ;;
        *)          extra+=("$1"); shift ;;
    esac
done
[[ -n "${bams}" && -n "${labels}" ]] || { echo "qc_correlation.sh: need --bams and --labels" >&2; exit 3; }
IFS=',' read -r -a bam_arr <<< "${bams}"
IFS=',' read -r -a label_arr <<< "${labels}"
if (( ${#bam_arr[@]} != ${#label_arr[@]} )); then
    echo "qc_correlation.sh: ${#bam_arr[@]} BAMs but ${#label_arr[@]} labels" >&2; exit 3
fi
if (( ${#bam_arr[@]} < 2 )); then
    echo "qc_correlation.sh: a correlation needs at least 2 BAMs (got ${#bam_arr[@]})" >&2; exit 3
fi
for i in "${!bam_arr[@]}"; do
    [[ "${bam_arr[$i]}" == "${label_arr[$i]}."* ]] || echo "qc_correlation.sh: WARNING label '${label_arr[$i]}' does not prefix BAM '${bam_arr[$i]}' (check the order)" >&2
done

multiBamSummary bins \
    --bamfiles "${bam_arr[@]}" \
    --labels "${label_arr[@]}" \
    --binSize "${bin_size}" \
    --distanceBetweenBins 0 \
    --extendReads \
    --numberOfProcessors "${threads}" \
    --outFileName correlation_matrix.npz \
    ${extra[@]+"${extra[@]}"}

plotCorrelation \
    --corData correlation_matrix.npz \
    --corMethod spearman \
    --whatToPlot heatmap \
    --skipZeros \
    --removeOutliers \
    --colorMap RdYlBu \
    --plotNumbers \
    --plotTitle "Spearman correlation of read counts (${bin_size} bp bins)" \
    --plotFileFormat svg \
    --outFileCorMatrix correlation_matrix.tsv \
    --plotFile correlation_heatmap.svg
