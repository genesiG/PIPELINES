include { runner; extraArgs } from './functions'

process FASTQC {
    tag "${fastq.simpleName}"
    // command.sh (from conf/base.config's afterScript) is excluded here and
    // published separately below — this publishDir has no pattern, so
    // without the exclusion it would also publish command.sh unrenamed,
    // colliding across every sample's task.
    publishDir "${params.outdir}/demultiplexing/qc/fastqc", mode: params.publish_mode, pattern: "*_fastqc.*"
    publishDir "${params.outdir}/demultiplexing/qc/fastqc/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${fastq.simpleName}_command.sh" }

    input:
    path fastq
    path step_script   // bin/fastqc.sh, or the project's `fastqc:` replacement

    output:
    path "*_fastqc.html", emit: html
    path "*_fastqc.zip",  emit: zip
    path "command.sh",    emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --fastq "${fastq}" \\
        --threads ${task.cpus} \\
        ${extraArgs('fastqc')}
    """

    stub:
    """
    touch "${fastq.simpleName}_fastqc.html" "${fastq.simpleName}_fastqc.zip"
    """
}
