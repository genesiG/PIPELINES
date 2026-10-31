#!/bin/bash
# qc_frip.sh -- FRiP of one sample against ITS OWN peak set.
# Called by FRIP (modules/local/frip.nf):
#   qc_frip.sh --sample S --bam S.final.bam --peaks peaks.bed --peak-source LABEL --flagstat S.final.flagstat.txt
# Arguments:
#   --sample       sample id
#   --bam          final BAM (alignments to test)
#   --peaks        the peak set of this sample (BED/narrowPeak/broadPeak, >= 3 columns); passed explicitly, never globbed
#   --peak-source  label of the peak set (e.g. the peakset name or caller); goes into the file name and the output row
#   --flagstat     S.final.flagstat.txt of the same sample (must start with "<sample>."; exit 3 otherwise)
# FRiP = alignments overlapping at least one peak / mapped alignments in the final BAM (flagstat "mapped" line).
#   Alignments are counted once each (bedtools intersect -u on the BAM, converted to BED records, one line per
#   alignment); both mates are counted in numerator and denominator.
# Output: S.<peak_source>.frip.tsv  (header + one row: sample_id peak_source peaks_file n_peaks
#   total_mapped_alignments reads_in_peaks frip)
# An empty peak file gives FRiP 0 with a warning (not an error).
# Exit 3: empty BAM, missing peak file, flagstat of another sample, or zero mapped alignments.
set -euo pipefail

sample="" bam="" peaks="" source_label="" flagstat=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)      sample="$2"; shift 2 ;;
        --bam)         bam="$2"; shift 2 ;;
        --peaks)       peaks="$2"; shift 2 ;;
        --peak-source) source_label="$2"; shift 2 ;;
        --flagstat)    flagstat="$2"; shift 2 ;;
        *) echo "qc_frip.sh: unknown argument $1" >&2; exit 2 ;;
    esac
done
[[ -n "${sample}" && -n "${source_label}" && -s "${bam}" && -s "${flagstat}" ]] || { echo "qc_frip.sh: missing/empty --sample, --peak-source, --bam or --flagstat" >&2; exit 3; }
[[ -e "${peaks}" ]] || { echo "qc_frip.sh: peak file '${peaks}' of sample ${sample} (${source_label}) is missing" >&2; exit 3; }
[[ "$(basename "${flagstat}")" == "${sample}."* ]] || { echo "qc_frip.sh: flagstat '${flagstat}' does not belong to sample ${sample}" >&2; exit 3; }

total=$(awk '/^[0-9]+ \+ [0-9]+ mapped \(/ {print $1; exit}' "${flagstat}")
[[ -n "${total}" && "${total}" -gt 0 ]] || { echo "qc_frip.sh: no mapped alignments in ${flagstat}" >&2; exit 3; }

n_peaks=$(grep -c -v -E '^(#|track|browser)' "${peaks}" || true)
# An empty peak set is a valid result (e.g. a shallow library or a group with no peaks):
# report n_peaks 0 and FRiP 0 with a warning instead of failing the run.
if (( n_peaks == 0 )); then
    echo "qc_frip.sh: WARNING peak set '${peaks}' of sample ${sample} (${source_label}) has no peaks; FRiP reported as 0" >&2
    in_peaks=0
else
    in_peaks=$(bedtools intersect -u -abam "${bam}" -b "${peaks}" -bed | wc -l)
fi

{
    printf 'sample_id\tpeak_source\tpeaks_file\tn_peaks\ttotal_mapped_alignments\treads_in_peaks\tfrip\n'
    awk -v s="${sample}" -v src="${source_label}" -v pk="$(basename "${peaks}")" -v np="${n_peaks}" -v t="${total}" -v r="${in_peaks}" \
        'BEGIN {printf "%s\t%s\t%s\t%d\t%d\t%d\t%.6f\n", s, src, pk, np, t, r, r / t}'
} > "${sample}.${source_label}.frip.tsv"
