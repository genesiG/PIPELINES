include { runner; extraArgs } from './functions'

process ALIGN {
    tag "${sample_id}"
    container '/project/rephilllab/nf_containers/rnaseq-align_1.0.sif'   // lab image: HISAT2 2.2.1 + samtools 1.23
    publishDir "${params.intermediate_dir}/aligned", mode: params.publish_mode, pattern: "*.aligned.bam*", enabled: params.publish_intermediates.toString() == 'true'
    publishDir "${params.outdir}/qc/hisat2", mode: params.publish_mode, pattern: "*.hisat2_summary.txt"
    publishDir "${params.outdir}/qc/hisat2/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_align_command.sh" }

    input:
    tuple val(sample_id), path(r1), path(r2)
    path index_files   // all HISAT2 index files of the prefix (<prefix>.N.ht2 / .ht2l), staged flat
    val index_base     // bare file-name prefix of the index
    val hisat2_args    // params.hisat2_args, passed to hisat2 after `--`
    path step_script   // bin/align.sh, or the project's `align` script replacement

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam"), path("${sample_id}.aligned.bam.bai"), emit: bam
    path "${sample_id}.hisat2_summary.txt", emit: summary
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --r1 "${r1}" \\
        --r2 "${r2}" \\
        --index-base "${index_base}" \\
        --threads ${task.cpus} \\
        --sort-mem ${task.ext.sort_mem ?: '1G'} \\
        -- ${hisat2_args ?: ''} ${extraArgs('align')}
    """

    stub:
    """
    touch "${sample_id}.aligned.bam" "${sample_id}.aligned.bam.bai"
    cat > "${sample_id}.hisat2_summary.txt" <<'EOF'
10000 reads; of these:
  10000 (100.00%) were paired; of these:
    1000 (10.00%) aligned concordantly 0 times
    8000 (80.00%) aligned concordantly exactly 1 time
    1000 (10.00%) aligned concordantly >1 times
90.00% overall alignment rate
EOF
    """
}
