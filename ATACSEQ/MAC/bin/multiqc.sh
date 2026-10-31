#!/bin/bash
# multiqc.sh -- one MultiQC report over every QC file staged in the task directory.
# Called by MULTIQC (modules/local/multiqc.nf):
#   multiqc.sh [multiqc args...]
# Inputs (staged flat by the module; MultiQC detects each by content):
#   FastQC zips (raw and trimmed), Trimmomatic logs, bowtie2 logs, sambamba markdup logs,
#   raw and final samtools flagstat/idxstats, MACS3 *_peaks.xls
# Any arguments (config.yaml `multiqc_args`) are passed to multiqc.
# Outputs: multiqc_report.html, multiqc_report_data/
set -euo pipefail

multiqc . --filename multiqc_report.html --force "$@"
