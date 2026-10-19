#!/bin/bash
# multiqc.sh — one MultiQC report over every QC file staged in the task
# directory (FastQC zips, the cutadapt log).
# Called by MULTIQC (modules/local/multiqc.nf):
#   multiqc.sh [multiqc args...]
# Any arguments (config.yaml `multiqc_args`) are passed to multiqc.
# Must write: multiqc_report.html, multiqc_report_data/
set -euo pipefail

multiqc . --filename multiqc_report.html --force "$@"
