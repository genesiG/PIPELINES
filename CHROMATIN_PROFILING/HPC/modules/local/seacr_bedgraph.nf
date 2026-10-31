include { runner } from './functions'

process SEACR_BEDGRAPH {
    tag "${sample_id}"
    container 'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3'
    publishDir { "${params.intermediate_dir}/seacr" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.intermediate_dir}/seacr/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_seacr_bedgraph_command.sh" }

    input:
    tuple val(sample_id), path(bam), path(bai)
    path step_script        // bin/seacr_bedgraph.sh, or the project's replacement

    output:
    tuple val(sample_id), path("${sample_id}.fragments.bedgraph"), emit: bedgraph
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        ${params.seacr_bedgraph_args ?: ''}
    """

    stub:
    """
    touch "${sample_id}.fragments.bedgraph"
    """
}
