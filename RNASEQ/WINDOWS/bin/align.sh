#!/bin/bash
# align.sh -- HISAT2 paired-end alignment piped into a coordinate-sorted, indexed BAM.
# Called by ALIGN (modules/local/align.nf):
#   align.sh --sample S --r1 R1.fastq.gz --r2 R2.fastq.gz --index-base BASE --threads N [-- hisat2 args...]
# Arguments:
#   --sample      sample id (output prefix)
#   --r1, --r2    paired FASTQ files (trimmed, or raw when trim is off)
#   --index-base  bare file-name prefix of the HISAT2 index (<prefix>.N.ht2 or .ht2l; the index files are staged in
#                 the work dir; config hisat2_index)
#   --threads     hisat2 -p and samtools -@ (task.cpus)
#   --sort-mem    samtools sort -m, memory per sort buffer (default 1G; ALIGN passes task.ext.sort_mem, config ext.sort_mem)
#   --            everything after it goes to hisat2 unchanged (config hisat2_args, default empty; legacy none)
# Command (legacy step_2_align_hisat2.py; the SAM is piped instead of written to disk):
#   hisat2 -p N -q --dta -x BASE -1 R1 -2 R2 --no-discordant --no-mixed --summary-file S.hisat2_summary.txt [args] \
#     | samtools view -b -@ N - \
#     | samtools sort -m SORTMEM -@ N -o S.aligned.bam -
#   samtools index S.aligned.bam
# Legacy differences: SAM written to disk (-S) then converted; samtools sort -m 2G; HISAT2 taken from the ~/.bashrc PATH.
#   No --no-unal, as legacy: unaligned records stay in the BAM for the statistics. --rna-strandness is not set here
#   (legacy none; use hisat2_args).
# The hisat2 summary goes to S.hisat2_summary.txt; hisat2's stderr stays on the task stderr.
# Exit 3 on missing arguments, an empty or missing FASTQ, no HISAT2 index files for --index-base, or a BAM with no reads.
# Outputs: S.aligned.bam, S.aligned.bam.bai, S.hisat2_summary.txt
set -euo pipefail

sample="" r1="" r2="" index_base="" threads=1 sort_mem="1G" extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)     sample="$2"; shift 2 ;;
        --r1)         r1="$2"; shift 2 ;;
        --r2)         r2="$2"; shift 2 ;;
        --index-base) index_base="$2"; shift 2 ;;
        --threads)    threads="$2"; shift 2 ;;
        --sort-mem)   sort_mem="$2"; shift 2 ;;
        --)           shift; extra+=("$@"); break ;;
        *)            extra+=("$1"); shift ;;
    esac
done

for v in sample r1 r2 index_base; do
    [[ -n "${!v}" ]] || { echo "align.sh: missing --${v//_/-}" >&2; exit 3; }
done
for v in r1 r2; do
    [[ -s "${!v}" ]] || { echo "align.sh: empty or missing FASTQ ${!v}" >&2; exit 3; }
done
ls "${index_base}".*.ht2* >/dev/null 2>&1 || { echo "align.sh: no HISAT2 index files for prefix '${index_base}'" >&2; exit 3; }

hisat2 -p "${threads}" -q --dta -x "${index_base}" -1 "${r1}" -2 "${r2}" \
    --no-discordant --no-mixed \
    --summary-file "${sample}.hisat2_summary.txt" \
    ${extra[@]+"${extra[@]}"} \
  | samtools view -b -@ "${threads}" - \
  | samtools sort -m "${sort_mem}" -@ "${threads}" -o "${sample}.aligned.bam" -

samtools index "${sample}.aligned.bam"

reads=$(samtools view -c -@ "${threads}" "${sample}.aligned.bam")
[[ "${reads}" -gt 0 ]] || { echo "align.sh: no reads in ${sample}.aligned.bam" >&2; exit 3; }
