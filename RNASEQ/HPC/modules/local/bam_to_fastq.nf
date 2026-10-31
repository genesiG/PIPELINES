include { runner; extraArgs } from './functions'

process BAM_TO_FASTQ {
    tag "${sample_id}"
    container 'quay.io/biocontainers/samtools:1.23--h96c455f_0'
    publishDir "${params.intermediate_dir}/fastq_from_bam", mode: params.publish_mode, pattern: "*.fastq.gz"
    publishDir "${params.intermediate_dir}/fastq_from_bam/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_bam_to_fastq_command.sh" }

    input:
    tuple val(sample_id), path(bam)
    path step_script   // bin/bam_to_fastq.sh, or the project's `bam_to_fastq:` replacement

    output:
    tuple val(sample_id), path("${sample_id}.R1.fastq.gz"), path("${sample_id}.R2.fastq.gz"), emit: reads
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --threads ${task.cpus} \\
        ${extraArgs('bam_to_fastq')}
    """

    stub:
    """
    printf '@stub/1\\nA\\n+\\nI\\n' | gzip -c > "${sample_id}.R1.fastq.gz"
    printf '@stub/2\\nA\\n+\\nI\\n' | gzip -c > "${sample_id}.R2.fastq.gz"
    """
}
