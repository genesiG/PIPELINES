/*
 * Joins the chunk BAMs of one library (samtools cat, in chunk order) and sums the chunk bowtie2 logs into one
 * S.bowtie2.log (bin/align_merge.sh calls bin/merge_bowtie2_logs.py). Runs as soon as the library's chunks are done.
 * chunk_bams and chunk_logs: staged flat, in chunk order.
 */
include { runner; extraArgs } from './functions'

process ALIGN_MERGE {
    tag "${sample_id}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    publishDir "${params.intermediate_dir}/aligned", mode: params.publish_mode, pattern: "*.aligned.bam", enabled: params.publish_intermediates.toString() == 'true'
    publishDir "${params.outdir}/qc/bowtie2", mode: params.publish_mode, pattern: "*.bowtie2.log"
    publishDir "${params.outdir}/qc/bowtie2/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_align_merge_command.sh" }

    input:
    tuple val(sample_id), path(chunk_bams), path(chunk_logs)   // in chunk order
    path log_merger    // bin/merge_bowtie2_logs.py, or the project's merge_bowtie2_logs: replacement
    path step_script   // bin/align_merge.sh, or the project's align_merge: replacement

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam"), emit: bam
    path "${sample_id}.bowtie2.log",                        emit: log
    path "command.sh",                                      emit: command_sh

    script:
    def bnames = [chunk_bams].flatten().collect { it.name }.join(',')
    def lnames = [chunk_logs].flatten().collect { it.name }.join(',')
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bams "${bnames}" \\
        --logs "${lnames}" \\
        --log-merger "${log_merger}" \\
        ${extraArgs('align_merge')}
    """

    stub:
    """
    touch "${sample_id}.aligned.bam"
    cat > "${sample_id}.bowtie2.log" <<'EOF'
10000 reads; of these:
  10000 (100.00%) were paired; of these:
    1000 (10.00%) aligned concordantly 0 times
    8000 (80.00%) aligned concordantly exactly 1 time
    1000 (10.00%) aligned concordantly >1 times
90.00% overall alignment rate
EOF
    """
}
