include { runner } from './functions'

process MARKDUP {
    tag "${sample_id}"
    container 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
    publishDir "${params.intermediate_dir}/markdup", mode: params.publish_mode, pattern: "*.markdup.bam*"
    publishDir "${params.outdir}/qc/markdup", mode: params.publish_mode, pattern: "*.markdup_metrics.txt"
    publishDir "${params.outdir}/qc/markdup/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_markdup_command.sh" }

    input:
    tuple val(sample_id), path(bam), path(bai)
    path step_script   // bin/markdup.sh, or the project's `markdup:` replacement

    output:
    tuple val(sample_id), path("${sample_id}.markdup.bam"), path("${sample_id}.markdup.bam.bai"), emit: bam
    path "${sample_id}.markdup_metrics.txt", emit: metrics
    path "command.sh",                       emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --mem-gb ${task.memory.toGiga()} \\
        --threads ${task.cpus} \\
        ${params.markdup_args ?: ''}
    """

    stub:
    """
    touch "${sample_id}.markdup.bam" "${sample_id}.markdup.bam.bai"
    printf '## METRICS CLASS\\tpicard.sam.DuplicationMetrics\\nLIBRARY\\tUNPAIRED_READS_EXAMINED\\tREAD_PAIRS_EXAMINED\\n' > "${sample_id}.markdup_metrics.txt"
    """
}
