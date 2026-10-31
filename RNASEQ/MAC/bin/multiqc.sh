#!/bin/bash
# multiqc.sh -- one MultiQC report over every QC file staged in the task directory.
# Called by MULTIQC (modules/local/multiqc.nf):
#   multiqc.sh [multiqc args...]
# Inputs (staged flat by the module; MultiQC detects each by content):
#   FastQC zips (raw and trimmed, both reads), Trimmomatic logs, HISAT2 summaries, samtools flagstat and idxstats
#   (raw and filtered), Picard MarkDuplicates metrics, featureCounts summaries, RSeQC infer_experiment outputs.
# Writes multiqc_config.yaml in the task directory, so sample names are cleaned to sample_id: the
#   extra_fn_clean_exts list removes the stage suffixes (added to MultiQC's default cleaning list).
#   FastQC entries keep their read label (R1 or R2, and .trim for trimmed reads), because the read label is part of
#   the FASTQ file name and is not a stage suffix; those rows are per read, not per sample.
# Any arguments (config.yaml `multiqc_args`) are passed to multiqc.
# Outputs: multiqc_report.html, multiqc_report_data/, multiqc_config.yaml
set -euo pipefail

cat > multiqc_config.yaml <<'EOF'
# Written by bin/multiqc.sh. Stage suffixes removed from file names so that sample names equal sample_id.
extra_fn_clean_exts:
  - "_featurecounts.txt.summary"
  - "_featurecounts.txt"
  - ".trimmomatic.log"
  - ".hisat2_summary.txt"
  - ".raw.flagstat.txt"
  - ".raw.idxstats.txt"
  - ".filtered.flagstat.txt"
  - ".filtered.idxstats.txt"
  - ".picard_metrics.txt"
  - ".infer_experiment.txt"
EOF

multiqc . --filename multiqc_report.html --force -c multiqc_config.yaml "$@"
