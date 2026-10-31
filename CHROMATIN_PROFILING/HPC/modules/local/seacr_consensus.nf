include { runner } from './functions'

process SEACR_CONSENSUS {
    tag "${group}"
    container 'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3'
    publishDir { "${params.outdir}/peaks/seacr" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/peaks/seacr/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${group}_seacr_consensus_command.sh" }

    input:
    tuple val(group), path(beds)
    val min_replicates
    val reciprocal_overlap
    path step_script        // bin/seacr_consensus.sh, or the project's replacement

    output:
    tuple val(group), path("${group}.seacr_consensus.bed"), emit: bed
    path "command.sh", emit: command_sh

    script:
    def names = [beds].flatten().collect { it.name }.join(',')
    """
    ${runner(step_script)} "${step_script}" \\
        --group "${group}" \\
        --beds "${names}" \\
        --min-replicates ${min_replicates} \\
        --reciprocal-overlap ${reciprocal_overlap} \\
        ${params.seacr_consensus_args ?: ''}
    """

    stub:
    """
    touch "${group}.seacr_consensus.bed"
    """
}
