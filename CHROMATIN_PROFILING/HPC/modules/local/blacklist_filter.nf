include { runner } from './functions'

process BLACKLIST_FILTER {
    tag "${sample_id}"
    container 'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3'
    // The name-sorted BAM is an intermediate and is not published; only command.sh is.
    publishDir "${params.outdir}/bam/final/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_blacklist_filter_command.sh" }

    input:
    tuple val(sample_id), path(target_bam)
    path blacklist
    path step_script   // bin/blacklist_filter.sh, or the project's `blacklist_filter:` replacement

    output:
    tuple val(sample_id), path("${sample_id}.noblacklist.namesorted.bam"), emit: bam
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${target_bam}" \\
        --blacklist "${blacklist}" \\
        ${params.blacklist_filter_args ?: ''}
    """

    stub:
    """
    touch "${sample_id}.noblacklist.namesorted.bam"
    """
}
