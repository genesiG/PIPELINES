include { runner } from './functions'

process DEDUP_SPLIT {
    tag "${sample_id}"
    container 'quay.io/biocontainers/samtools:1.23--h96c455f_0'
    publishDir "${params.outdir}/bam/spikein", mode: params.publish_mode, pattern: "*.spikein.bam*"
    publishDir "${params.outdir}/bam/stats", mode: params.publish_mode, pattern: "*.dedup_counts.tsv"
    publishDir "${params.outdir}/bam/stats", mode: params.publish_mode, pattern: "*.spikein.idxstats.txt"
    publishDir "${params.outdir}/bam/stats/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_dedup_split_command.sh" }

    input:
    tuple val(sample_id), path(bam), path(bai), val(is_control)
    val remove_duplicates
    val spikein_dedup
    val spike_prefix
    val mito_chrom
    path step_script   // bin/dedup_split.sh, or the project's `dedup_split:` replacement

    output:
    tuple val(sample_id), path("${sample_id}.target.namesorted.bam"), emit: target
    tuple val(sample_id), path("${sample_id}.spikein.bam"), path("${sample_id}.spikein.bam.bai"), emit: spikein
    path "${sample_id}.spikein.idxstats.txt", emit: spike_idxstats
    path "${sample_id}.dedup_counts.tsv", emit: counts
    path "command.sh",                    emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --is-control ${is_control} \\
        --remove-duplicates ${remove_duplicates} \\
        --spikein-dedup ${spikein_dedup} \\
        --spike-prefix "${spike_prefix}" \\
        --mito-chrom "${mito_chrom}" \\
        --threads ${task.cpus} \\
        ${params.dedup_split_args ?: ''}
    """

    stub:
    """
    touch "${sample_id}.target.namesorted.bam" "${sample_id}.spikein.bam" "${sample_id}.spikein.bam.bai"
    printf 'ecoli_chr1\\t1000\\t0\\t0\\n*\\t0\\t0\\t0\\n' > "${sample_id}.spikein.idxstats.txt"
    printf 'sample_id\\tstage\\treads\\n${sample_id}\\ttarget_before_dedup\\t0\\n${sample_id}\\ttarget_after_dedup\\t0\\n${sample_id}\\tspikein_before_dedup\\t0\\n${sample_id}\\tspikein_after_dedup\\t0\\n' > "${sample_id}.dedup_counts.tsv"
    """
}
