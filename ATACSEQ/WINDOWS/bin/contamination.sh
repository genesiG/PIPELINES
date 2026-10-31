#!/bin/bash
# contamination.sh -- screen one sample's unmapped read pairs against a contaminant bowtie2 index.
# Port of legacy step_3b_contamination_check.py without bsub/conda. Called by CONTAMINATION (modules/local/contamination.nf):
#   contamination.sh --sample S --bam S.aligned.bam --index-base PREFIX --threads N
# Arguments:
#   --sample      sample id
#   --bam         unsorted aligned BAM (ALIGN output, S.aligned.bam)
#   --index-base  bare file-name prefix of the contaminant bowtie2 index (index files are staged in the work dir)
#   --threads     bowtie2 -p (task.cpus)
#   other arguments are appended to bowtie2.
# Steps:
#   1. samtools view -u -f 12 -F 256: pairs with both mates unmapped (no secondary records) -> bedtools bamtofastq
#      -> temporary R1/R2 FASTQ (paired FASTQ must go to files; they are deleted at the end).
#      Unmapped individual reads = R1 records + R2 records (= 2 x unmapped pairs).
#   2. bowtie2 -p N --no-unal -x PREFIX -1 R1 -2 R2 | samtools view -F 4 | cut -f3: only the SAM text is piped, never a
#      SAM file. Aligned records are counted per reference name. organism = the reference name up to the first '|'
#      (the whole name when there is no '|'), as in the legacy script. Records count mates separately.
#   3. pct_of_unmapped = 100 * aligned_reads / total unmapped individual reads.
# Outputs:
#   S.contamination.tsv          header: organism, aligned_reads, pct_of_unmapped; one row per organism, most reads first
#   S.contamination_summary.txt  unmapped read pairs and reads, contaminant-aligned reads and %, the organism table,
#                                bowtie2 alignment statistics
# No unmapped pairs: both files are still written (TSV with the header only; summary says "No unmapped reads") and exit 0.
# Exit 3 when --bam is missing or empty, or --index-base has no bowtie2 index files.
set -euo pipefail

sample="" bam="" index_base="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)     sample="$2"; shift 2 ;;
        --bam)        bam="$2"; shift 2 ;;
        --index-base) index_base="$2"; shift 2 ;;
        --threads)    threads="$2"; shift 2 ;;
        *)            extra+=("$1"); shift ;;
    esac
done
for v in sample bam index_base; do
    [[ -n "${!v}" ]] || { echo "contamination.sh: missing --${v//_/-}" >&2; exit 3; }
done
[[ -s "${bam}" ]] || { echo "contamination.sh: empty or missing input BAM ${bam}" >&2; exit 3; }
ls "${index_base}"*.bt2* >/dev/null 2>&1 || { echo "contamination.sh: no bowtie2 index files for prefix '${index_base}'" >&2; exit 3; }

unmapped_bam="${sample}.unmapped_pairs.tmp.bam"
r1="${sample}.unmapped_R1.tmp.fastq"
r2="${sample}.unmapped_R2.tmp.fastq"
org_counts="${sample}.organism_counts.tmp.tsv"
bt2_log="${sample}.contamination_bowtie2.log"
tsv="${sample}.contamination.tsv"
summary="${sample}.contamination_summary.txt"

# ---- 1. unmapped pairs -> FASTQ ----
samtools view -u -f 12 -F 256 "${bam}" > "${unmapped_bam}"
bedtools bamtofastq -i "${unmapped_bam}" -fq "${r1}" -fq2 "${r2}"
n_r1=$(( $(wc -l < "${r1}") / 4 ))
n_r2=$(( $(wc -l < "${r2}") / 4 ))
total_unmapped=$(( n_r1 + n_r2 ))
printf 'organism\taligned_reads\tpct_of_unmapped\n' > "${tsv}"

if (( total_unmapped == 0 )); then
    {
        echo "Contamination report for ${sample}"
        echo "No unmapped reads"
        echo "No unmapped read pairs were found, so the contaminant screen was skipped."
    } > "${summary}"
    rm -f "${unmapped_bam}" "${r1}" "${r2}"
    exit 0
fi

# ---- 2. screen against the contaminant index; count aligned records per reference name ----
bowtie2 -p "${threads}" --no-unal -x "${index_base}" -1 "${r1}" -2 "${r2}" ${extra[@]+"${extra[@]}"} \
    2> "${bt2_log}" \
  | samtools view -F 4 - \
  | cut -f3 \
  | awk -F'|' '{ n[$1]++ } END { for (o in n) printf "%s\t%d\n", o, n[o] }' \
  | sort -t$'\t' -k2,2nr -k1,1 > "${org_counts}"

# ---- 3. percentages ----
aligned_reads=$(awk -F'\t' '{ s += $2 } END { print s + 0 }' "${org_counts}")
pct_total=$(awk -v a="${aligned_reads}" -v t="${total_unmapped}" 'BEGIN { printf "%.4f", a * 100.0 / t }')
awk -F'\t' -v t="${total_unmapped}" 'BEGIN { OFS = "\t" } { printf "%s\t%d\t%.4f\n", $1, $2, $2 * 100.0 / t }' \
    "${org_counts}" >> "${tsv}"

{
    echo "================================================================================"
    echo "Contamination Report for ${sample}"
    echo "================================================================================"
    echo
    echo "OVERALL STATISTICS:"
    echo "------------------"
    echo "Total unmapped read pairs (R1): ${n_r1}"
    echo "Total unmapped individual reads (R1 + R2): ${total_unmapped}"
    echo "Total contaminant-aligned reads: ${aligned_reads} (${pct_total}%)"
    echo
    echo "CONTAMINATION BY ORGANISM:"
    echo "-------------------------"
    printf "%-40s %15s %15s\n" "Organism" "Aligned Reads" "% of Unmapped"
    printf "%-40s %15s %15s\n" "--------" "-------------" "-------------"
    awk -F'\t' -v t="${total_unmapped}" '{ printf "%-40s %15d %14.4f%%\n", $1, $2, $2 * 100.0 / t }' "${org_counts}"
    echo
    echo "BOWTIE2 ALIGNMENT STATISTICS:"
    echo "----------------------------"
    cat "${bt2_log}"
} > "${summary}"

rm -f "${unmapped_bam}" "${r1}" "${r2}" "${org_counts}"
