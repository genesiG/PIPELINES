include { runner } from './functions'

process FASTQC {
    tag "${sample_id} ${fastq.simpleName}"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    // Specific pattern: matches the report files only, never command.sh.
    publishDir "${params.outdir}/qc/fastqc", mode: params.publish_mode, pattern: "*_fastqc.{html,zip}"
    publishDir "${params.outdir}/qc/fastqc/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_${fastq.simpleName}_fastqc_command.sh" }

    input:
    tuple val(sample_id), path(fastq)
    path step_script   // bin/fastqc.sh, or the project's `fastqc:` replacement

    output:
    tuple val(sample_id), path("*_fastqc.html"), emit: html
    tuple val(sample_id), path("*_fastqc.zip"),  emit: zip
    path "command.sh",                           emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --fastq "${fastq}" \\
        --threads ${task.cpus} \\
        ${params.fastqc_args ?: ''}
    """

    stub:
    """
    touch "${fastq.simpleName}_fastqc.html" "${fastq.simpleName}_fastqc.zip"
    """
}
