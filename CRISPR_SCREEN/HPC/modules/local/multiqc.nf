include { runner } from './functions'

process MULTIQC {
    tag "${params.project_name}"
    container 'quay.io/biocontainers/multiqc:1.34--pyhdfd78af_0'
    // command.sh excluded here, published separately below -- see fastqc.nf
    // for why (no pattern here would otherwise also match it).
    publishDir "${params.outdir}/demultiplexing/qc/multiqc", mode: params.publish_mode, pattern: "multiqc_report*"
    publishDir "${params.outdir}/demultiplexing/qc/multiqc/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${params.project_name}_command.sh" }

    input:
    path '*'           // any mix of fastqc zips and the cutadapt log; multiqc auto-detects
    path step_script   // bin/multiqc.sh, or the project's `multiqc:` replacement

    output:
    path "multiqc_report.html",      emit: report
    path "multiqc_report_data",       emit: data   // MultiQC derives the data dir name from --filename
    path "command.sh",                emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" ${params.multiqc_args ?: ''}
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_report_data
    """
}
