#!/bin/bash
# guide_anchor_trim.sh — positions every sample's guide exactly, regardless
# of the (possibly variable) length of whatever comes before or after the
# anchor, by searching for the anchor sequence per read (cutadapt, real
# alignment) rather than assuming a fixed offset from the read start.
# Called by GUIDE_ANCHOR_TRIM (modules/local/guide_anchor_trim.nf):
#   guide_anchor_trim.sh --anchor-file lib_anchor.txt --fastq "s1.fastq.gz s2.fastq.gz" \
#                        --labels "s1,s2" --outdir . [cutadapt args...]
# --fastq and --labels are in the same order; each output is named
# <label>_trimmed.fastq.gz.
# Any other arguments (config.yaml `guide_anchor_trim_args`, e.g.
# "--error-rate=0.05") are passed to cutadapt.
# Must write, per sample: <label>_trimmed.fastq.gz, and one combined
# <label>_anchor_trim_stats.csv per invocation (sample,reads,anchor_found_reads).
#
# No --discard-untrimmed: a read where the anchor genuinely can't be found
# (even with the configured error tolerance) is passed through UNCHANGED,
# not dropped. Its guide-length slice at position 0 is then essentially
# arbitrary sequence, which will almost certainly not match any real guide
# (or its 1-mismatch neighborhood) — so it still gets correctly classified
# as "far-miss" downstream, keeping the read accounting complete: the total
# read count into mageck count matches the total that came out of
# demultiplexing, with nothing silently dropped from the denominator.
set -euo pipefail

anchor_file="" fastqs="" labels="" outdir="." extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --anchor-file) anchor_file="$2"; shift 2 ;;
        --fastq)       fastqs="$2"; shift 2 ;;
        --labels)      labels="$2"; shift 2 ;;
        --outdir)      outdir="$2"; shift 2 ;;
        *)             extra+=("$1"); shift ;;
    esac
done
mkdir -p "${outdir}"

anchor=$(tr -d '[:space:]' < "${anchor_file}")
if [[ -z "${anchor}" ]]; then
    echo "ERROR: ${anchor_file} is empty — no anchor to search for." >&2
    exit 1
fi
# Require a near-full-length anchor match (a couple of bases' slack) so a
# short coincidental match doesn't trigger a spurious trim.
overlap=$(( ${#anchor} - 2 ))
[[ "${overlap}" -lt 4 ]] && overlap=4

read -ra fastq_files <<< "${fastqs}"
IFS=',' read -ra sample_labels <<< "${labels}"
if [[ "${#fastq_files[@]}" -ne "${#sample_labels[@]}" ]]; then
    echo "ERROR: ${#fastq_files[@]} --fastq file(s) but ${#sample_labels[@]} --labels." >&2
    exit 1
fi

stats="${outdir}/anchor_trim_stats.csv"
echo "sample,reads,anchor_found_reads,pct_anchor_found" > "${stats}"

for i in "${!fastq_files[@]}"; do
    fq="${fastq_files[$i]}"
    label="${sample_labels[$i]}"
    out="${outdir}/${label}_trimmed.fastq.gz"
    log="${outdir}/${label}_anchor_trim.log"

    cutadapt \
        -g "${anchor}" \
        -e 0.1 \
        --overlap "${overlap}" \
        --cores 0 \
        ${extra[@]+"${extra[@]}"} \
        -o "${out}" \
        "${fq}" \
        > "${log}" 2>&1

    # cutadapt reports "Total reads processed:  N" and "Reads with adapters:
    # N (P%)". awk, not grep -P: the cutadapt image's grep is BusyBox and
    # has no -P (Perl regex) support.
    total=$(awk -F: '/^Total reads processed:/ { v=$2; gsub(/[, ]/, "", v); print v; exit }' "${log}")
    found=$(awk -F: '/^Reads with adapters:/ { v=$2; sub(/\(.*/, "", v); gsub(/[, ]/, "", v); print v; exit }' "${log}")
    total="${total:-0}"
    found="${found:-0}"
    pct="0.00"
    if [[ "${total}" -gt 0 ]]; then
        pct=$(awk -v f="${found}" -v t="${total}" 'BEGIN { printf "%.2f", f / t * 100 }')
    fi
    echo "${label},${total},${found},${pct}" >> "${stats}"
done
