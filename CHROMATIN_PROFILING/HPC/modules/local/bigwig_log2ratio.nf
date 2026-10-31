include { runner } from './functions'

process BIGWIG_LOG2RATIO {
    tag "${sample_id}"
    container 'quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0'
    publishDir { "${params.outdir}/bigwig" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/bigwig/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}.vs_${control_id}_log2ratio_command.sh" }

    input:
    tuple val(sample_id), path(target_bw), val(control_id), path(control_bw)
    val bin_size
    path step_script        // bin/bigwig_log2ratio.sh, or the project's replacement

    output:
    tuple val(sample_id), path("${sample_id}.vs_${control_id}.log2ratio.bw"), emit: bw
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --target "${target_bw}" \\
        --control "${control_bw}" \\
        --control-id "${control_id}" \\
        --bin-size ${bin_size} \\
        --threads ${task.cpus} \\
        ${params.bigwig_log2ratio_args ?: ''}
    """

    stub:
    """
    touch "${sample_id}.vs_${control_id}.log2ratio.bw"
    """
}
