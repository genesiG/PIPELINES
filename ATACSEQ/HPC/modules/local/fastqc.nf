include { runner; extraArgs } from './functions'

process FASTQC {
    tag "${sample_id} ${fastq.simpleName}"
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir { "${params.outdir}/qc/fastqc/${stage}" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/qc/fastqc/${stage}/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_${fastq.name.replace('.fastq.gz', '')}_fastqc_command.sh" }

    input:
    tuple val(sample_id), path(fastq)
    val stage          // raw | trimmed (publish directory)
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
        ${extraArgs('fastqc')}
    """

    stub:
    def base = fastq.name.replace('.fastq.gz', '')
    """
    touch "${base}_fastqc.html" "${base}_fastqc.zip"
    """
}
