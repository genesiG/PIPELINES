/*
 * Parses the cutadapt demux log into per-sample read-count / uniformity
 * metrics (Gini, CV, min/max fold) — ported from step_1qc.py's
 * run_cutadapt_stats(). Runs in the cutadapt container since no extra
 * dependencies are needed beyond the Python stdlib.
 */
include { runner; extraArgs } from './functions'

process CUTADAPT_STATS {
    tag "${params.project_name}"
    container 'quay.io/biocontainers/cutadapt:5.2--py39hbcbf7aa_0'
    // command.sh excluded here, published separately below -- see fastqc.nf
    // for why (no pattern here would otherwise also match it).
    publishDir "${params.outdir}/demultiplexing/qc/cutadapt_stats", mode: params.publish_mode, pattern: "${params.project_name}_step1_*"
    publishDir "${params.outdir}/demultiplexing/qc/cutadapt_stats/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${params.project_name}_step1_command.sh" }

    input:
    path cutadapt_log
    path samplesheet
    path helper_script   // bin/cutadapt_stats.py, or the project's `cutadapt_stats:` replacement

    output:
    path "${params.project_name}_step1_demux_metrics.json", emit: json
    path "${params.project_name}_step1_sample_counts.csv",  emit: csv
    path "command.sh",                                        emit: command_sh

    script:
    """
    ${runner(helper_script)} "${helper_script}" \\
        --log "${cutadapt_log}" \\
        --samplesheet "${samplesheet}" \\
        --project-name "${params.project_name}" \\
        --outdir . \\
        ${extraArgs('cutadapt_stats')}
    """

    stub:
    """
    touch "${params.project_name}_step1_demux_metrics.json" "${params.project_name}_step1_sample_counts.csv"
    """
}
