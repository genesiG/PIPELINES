include { runner } from './functions'

process BAM_PREALIGNED {
    tag "${sample_id}"
    container 'quay.io/biocontainers/samtools:1.23--h96c455f_0'
    publishDir "${params.outdir}/bam/final", mode: params.publish_mode, pattern: "*.final.bam*"
    publishDir "${params.outdir}/bam/stats", mode: params.publish_mode, pattern: "*.flagstat.txt"
    publishDir "${params.outdir}/bam/stats", mode: params.publish_mode, pattern: "*.idxstats.txt"
    publishDir "${params.outdir}/bam/final/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_bam_prealigned_command.sh" }

    input:
    tuple val(sample_id), path(bam, stageAs: "input/*")   // own folder: an input named S.final.bam must not collide with the output
    path step_script   // bin/bam_prealigned.sh, or the project's `bam_prealigned:` replacement

    output:
    tuple val(sample_id), path("${sample_id}.final.bam"), path("${sample_id}.final.bam.bai"), emit: bam
    path "${sample_id}.final.flagstat.txt", emit: flagstat
    path "${sample_id}.final.idxstats.txt", emit: idxstats
    path "command.sh",                      emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --threads ${task.cpus} \\
        ${params.bam_prealigned_args ?: ''}
    """

    stub:
    """
    touch "${sample_id}.final.bam" "${sample_id}.final.bam.bai"
    echo "0 + 0 in total (QC-passed reads + QC-failed reads)" > "${sample_id}.final.flagstat.txt"
    printf 'stub_chr\\t1000\\t0\\t0\\n' > "${sample_id}.final.idxstats.txt"
    """
}
