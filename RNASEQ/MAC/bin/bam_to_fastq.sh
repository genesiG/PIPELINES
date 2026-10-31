#!/bin/bash
# bam_to_fastq.sh -- paired-end FASTQ files from one BAM (legacy step_0d_bam_to_fastq.py, paired-end branch).
# Called by BAM_TO_FASTQ (modules/local/bam_to_fastq.nf):
#   bam_to_fastq.sh --sample S --bam IN.bam --threads N [extra samtools fastq args...]
# Arguments:
#   --sample   sample id (output prefix)
#   --bam      input BAM (name-collated here, so it need not be sorted by name)
#   --threads  samtools -@ (task.cpus)
#   other arguments are appended to samtools fastq (none by default)
# Command (legacy command, plus -@ and -n as the lab spec asks):
#   samtools collate -u -O -@ N IN.bam S.collate.tmp \
#     | samtools fastq -@ N -1 S.R1.fastq -2 S.R2.fastq -0 /dev/null -s /dev/null -n -
#   gzip S.R1.fastq S.R2.fastq
# Legacy differences: step_0d and the batch files that ran (Original_Data/fastq_files/log/*.batch) have no -@ and
#   no -n (read names keep the /1 /2 suffix the legacy FASTQ had). The legacy also wrote uncompressed FASTQ and
#   gzipped afterwards; the same is done here, so the output does not depend on samtools' compression-by-extension.
# Exit 3 when the BAM is missing or empty, or when a FASTQ comes out empty (no read pairs, e.g. a single-end BAM).
# Outputs: S.R1.fastq.gz, S.R2.fastq.gz
set -euo pipefail

sample="" bam="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)  sample="$2"; shift 2 ;;
        --bam)     bam="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        *)         extra+=("$1"); shift ;;
    esac
done

for v in sample bam; do
    [[ -n "${!v}" ]] || { echo "bam_to_fastq.sh: missing --${v}" >&2; exit 3; }
done
[[ -s "${bam}" ]] || { echo "bam_to_fastq.sh: empty or missing input BAM ${bam}" >&2; exit 3; }

samtools collate -u -O -@ "${threads}" "${bam}" "${sample}.collate.tmp" \
    | samtools fastq -@ "${threads}" -1 "${sample}.R1.fastq" -2 "${sample}.R2.fastq" \
        -0 /dev/null -s /dev/null -n ${extra[@]+"${extra[@]}"} -

for r in R1 R2; do
    [[ -s "${sample}.${r}.fastq" ]] || { echo "bam_to_fastq.sh: no read pairs in ${bam} (${sample}.${r}.fastq is empty; single-end BAM?)" >&2; exit 3; }
done

gzip "${sample}.R1.fastq" "${sample}.R2.fastq"
