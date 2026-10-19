#!/bin/bash
# fastqc.sh — FastQC on one per-sample FASTQ.
# Called by FASTQC (modules/local/fastqc.nf):
#   fastqc.sh --fastq sample.fastq.gz --threads N [fastqc args...]
# Any other arguments (config.yaml `fastqc_args`) are passed to fastqc.
# Must write: *_fastqc.html, *_fastqc.zip
set -euo pipefail

fastq="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq)   fastq="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        *)         extra+=("$1"); shift ;;
    esac
done

fastqc -o . --extract --threads "${threads}" ${extra[@]+"${extra[@]}"} "${fastq}"
