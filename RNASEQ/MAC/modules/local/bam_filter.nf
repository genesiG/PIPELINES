include { runner; extraArgs } from './functions'

process BAM_FILTER {
    tag "${sample_id}"
    container 'quay.io/biocontainers/samtools:1.23--h96c455f_0'
    publishDir "${params.intermediate_dir}/filtered", mode: params.publish_mode, pattern: "*.filtered.bam*", enabled: params.publish_intermediates.toString() == 'true'
    publishDir "${params.outdir}/qc/samtools", mode: params.publish_mode, pattern: "*.flagstat.txt"
    publishDir "${params.outdir}/qc/samtools", mode: params.publish_mode, pattern: "*.idxstats.txt"
    publishDir "${params.outdir}/qc/samtools/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_bam_filter_command.sh" }

    input:
    tuple val(sample_id), path(aligned_bam), path(aligned_bai)
    val mapq
    val remove_chroms   // comma list of contigs to drop (params.remove_chroms); '' = keep all contigs
    path step_script    // bin/bam_filter.sh, or the project's `bam_filter` script replacement

    output:
    tuple val(sample_id), path("${sample_id}.filtered.bam"), path("${sample_id}.filtered.bam.bai"), emit: bam
    path "${sample_id}.raw.flagstat.txt",      emit: raw_flagstat
    path "${sample_id}.raw.idxstats.txt",      emit: raw_idxstats
    path "${sample_id}.filtered.flagstat.txt", emit: filtered_flagstat
    path "${sample_id}.filtered.idxstats.txt", emit: filtered_idxstats
    path "command.sh",                         emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${aligned_bam}" \\
        --mapq ${mapq} \\
        --remove-chroms "${remove_chroms ?: ''}" \\
        --threads ${task.cpus} \\
        ${extraArgs('bam_filter')}
    """

    stub:
    """
    touch "${sample_id}.filtered.bam" "${sample_id}.filtered.bam.bai"
    echo "0 + 0 in total (QC-passed reads + QC-failed reads)" > "${sample_id}.raw.flagstat.txt"
    echo "0 + 0 in total (QC-passed reads + QC-failed reads)" > "${sample_id}.filtered.flagstat.txt"
    printf 'stub_chr\\t1000\\t0\\t0\\n' > "${sample_id}.raw.idxstats.txt"
    printf 'stub_chr\\t1000\\t0\\t0\\n' > "${sample_id}.filtered.idxstats.txt"
    """
}
