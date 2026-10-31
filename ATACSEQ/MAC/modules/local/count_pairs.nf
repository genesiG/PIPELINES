/*
 * Counts the read pairs of the reads ALIGN reads (trimmed when trim is on). Used only to plan the ALIGN chunks
 * (align_chunk_pairs > 0); bin/count_pairs.sh exits 3 on a data problem (R1 and R2 differ, or a line count not divisible by 4).
 */
include { runner; extraArgs } from './functions'

process COUNT_PAIRS {
    tag "${sample_id}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    publishDir "${params.outdir}/qc/bowtie2/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_count_pairs_command.sh" }

    input:
    tuple val(sample_id), path(r1), path(r2)
    path step_script   // bin/count_pairs.sh, or the project's count_pairs: replacement

    output:
    tuple val(sample_id), path("${sample_id}.pair_count.txt"), emit: count
    path "command.sh",                                         emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --r1 "${r1}" \\
        --r2 "${r2}" \\
        ${extraArgs('count_pairs')}
    """

    stub:
    """
    echo 25000 > "${sample_id}.pair_count.txt"
    """
}
