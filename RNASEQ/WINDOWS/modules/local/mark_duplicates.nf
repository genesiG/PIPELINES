include { runner; extraArgs } from './functions'

process MARK_DUPLICATES {
    tag "${sample_id}"
    container 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
    publishDir "${params.outdir}/bam/final", mode: params.publish_mode, pattern: "*.final.bam*"
    publishDir "${params.outdir}/qc/picard", mode: params.publish_mode, pattern: "*.picard_metrics.txt"
    publishDir "${params.outdir}/qc/picard/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_mark_duplicates_command.sh" }

    input:
    tuple val(sample_id), path(bam), path(bai)
    val mark            // 'true' = Picard MarkDuplicates; 'false' = S.final.bam is a copy of the filtered BAM
    val remove          // 'true' = REMOVE_DUPLICATES (params.remove_duplicates); 'false' = mark only
    val pixel_distance  // OPTICAL_DUPLICATE_PIXEL_DISTANCE (params.optical_duplicate_pixel_distance)
    path step_script    // bin/mark_duplicates.sh, or the project's `mark_duplicates` script replacement

    output:
    tuple val(sample_id), path("${sample_id}.final.bam"), path("${sample_id}.final.bam.bai"), emit: bam
    path "${sample_id}.picard_metrics.txt", optional: true, emit: metrics   // only when mark is 'true'
    path "command.sh", emit: command_sh

    script:
    def java_gb = Math.max(1, (task.memory.toGiga() * 0.75) as int)   // Java heap from the task's memory
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --mark ${mark} \\
        --remove ${remove} \\
        --pixel-distance ${pixel_distance} \\
        --java-gb ${java_gb} \\
        --threads ${task.cpus} \\
        ${extraArgs('mark_duplicates')}
    """

    stub:
    def metrics_cmd = mark.toString() == 'true' ? "echo 'stub picard metrics' > '${sample_id}.picard_metrics.txt'" : ':'
    """
    touch "${sample_id}.final.bam" "${sample_id}.final.bam.bai"
    ${metrics_cmd}
    """
}
