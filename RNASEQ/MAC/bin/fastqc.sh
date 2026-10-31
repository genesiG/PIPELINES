#!/bin/bash
# fastqc.sh -- FastQC on one FASTQ file (one task per FASTQ; both reads of every pair, raw and trimmed).
# Called by FASTQC_RAW and FASTQC_TRIM (modules/local/fastqc.nf):
#   fastqc.sh --fastq reads.fastq.gz --threads N [extra fastqc args...]
# Arguments:
#   --fastq    input FASTQ (.fastq/.fastq.gz)
#   --threads  fastqc threads (task.cpus)
#   any other arguments (config.yaml `fastqc_args`) are passed to fastqc.
# Legacy: fastqc --extract on the trimmed R1 only (spec: both reads, raw and trimmed).
# Outputs: <name>_fastqc.html, <name>_fastqc.zip (and the extracted <name>_fastqc/ directory)
set -euo pipefail

fastq="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq)   fastq="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        *)         extra+=("$1"); shift ;;
    esac
done

[[ -n "${fastq}" && -s "${fastq}" ]] || { echo "fastqc.sh: missing or empty --fastq '${fastq}'" >&2; exit 3; }

fastqc -o . --extract --threads "${threads}" ${extra[@]+"${extra[@]}"} "${fastq}"
