include { runner; extraArgs } from './functions'

process BAM_FINALIZE {
    tag "${sample_id}"
    container 'quay.io/biocontainers/samtools:1.23--h96c455f_0'
    publishDir "${params.outdir}/bam/final", mode: params.publish_mode, pattern: "*.final.bam*"
    publishDir "${params.outdir}/bam/stats", mode: params.publish_mode, pattern: "*.flagstat.txt"
    publishDir "${params.outdir}/bam/stats", mode: params.publish_mode, pattern: "*.idxstats.txt"
    publishDir "${params.outdir}/bam/stats", mode: params.publish_mode, pattern: "*.orphan_check.txt"
    publishDir "${params.outdir}/bam/final/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_bam_finalize_command.sh" }

    input:
    tuple val(sample_id), path(bam)
    path step_script   // bin/bam_finalize.sh, or the project's `bam_finalize:` replacement

    output:
    tuple val(sample_id), path("${sample_id}.final.bam"), path("${sample_id}.final.bam.bai"), emit: bam
    path "${sample_id}.final.flagstat.txt", emit: flagstat
    path "${sample_id}.final.idxstats.txt", emit: idxstats
    path "${sample_id}.orphan_check.txt",   emit: orphan_check
    path "command.sh",                      emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --threads ${task.cpus} \\
        ${extraArgs('bam_finalize')}
    """

    stub:
    """
    touch "${sample_id}.final.bam" "${sample_id}.final.bam.bai"
    echo "0 + 0 in total (QC-passed reads + QC-failed reads)" > "${sample_id}.final.flagstat.txt"
    printf 'stub_chr\\t1000\\t0\\t0\\n' > "${sample_id}.final.idxstats.txt"
    echo OK > "${sample_id}.orphan_check.txt"
    """
}
