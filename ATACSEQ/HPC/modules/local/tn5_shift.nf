include { runner; extraArgs } from './functions'

process TN5_SHIFT {
    tag "${sample_id}"
    container 'quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0'
    publishDir "${params.outdir}/bam/tn5", mode: params.publish_mode, pattern: "*.tn5.bam*"
    publishDir "${params.outdir}/bam/tn5/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_tn5_shift_command.sh" }

    input:
    tuple val(sample_id), path(bam), path(bai)
    path step_script   // bin/tn5_shift.sh, or the project's `tn5_shift:` replacement

    output:
    tuple val(sample_id), path("${sample_id}.tn5.bam"), path("${sample_id}.tn5.bam.bai"), emit: bam
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --threads ${task.cpus} \\
        ${extraArgs('tn5_shift')}
    """

    stub:
    """
    touch "${sample_id}.tn5.bam" "${sample_id}.tn5.bam.bai"
    """
}
