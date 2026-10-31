include { runner } from './functions'

process BAM_FILTER {
    tag "${sample_id}"
    container 'quay.io/biocontainers/samtools:1.23--h96c455f_0'
    publishDir "${params.intermediate_dir}/filtered", mode: params.publish_mode, pattern: "*.filtered.bam*", enabled: params.publish_intermediates.toString() == 'true'
    publishDir "${params.outdir}/bam/stats", mode: params.publish_mode, pattern: "*.flagstat.txt"
    publishDir "${params.outdir}/bam/stats", mode: params.publish_mode, pattern: "*.idxstats.txt"
    publishDir "${params.outdir}/bam/stats/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_bam_filter_command.sh" }

    input:
    tuple val(sample_id), path(aligned_bam)
    val mapq
    path step_script   // bin/bam_filter.sh, or the project's `bam_filter:` replacement

    output:
    tuple val(sample_id), path("${sample_id}.filtered.bam"), path("${sample_id}.filtered.bam.bai"), emit: bam
    path "${sample_id}.raw.flagstat.txt", emit: raw_flagstat
    path "${sample_id}.raw.idxstats.txt", emit: raw_idxstats
    path "command.sh",                    emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${aligned_bam}" \\
        --mapq ${mapq} \\
        --threads ${task.cpus} \\
        ${params.bam_filter_args ?: ''}
    """

    stub:
    """
    touch "${sample_id}.filtered.bam" "${sample_id}.filtered.bam.bai"
    echo "0 + 0 in total (QC-passed reads + QC-failed reads)" > "${sample_id}.raw.flagstat.txt"
    printf 'stub_chr\\t1000\\t0\\t0\\n' > "${sample_id}.raw.idxstats.txt"
    """
}
