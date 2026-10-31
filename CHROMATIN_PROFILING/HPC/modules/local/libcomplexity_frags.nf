include { runner } from './functions'

process LIBCOMPLEXITY_FRAGS {
    tag "${sample_id}"
    container 'quay.io/biocontainers/samtools:1.23--h96c455f_0'
    // BAM input: no cache 'deep'
    publishDir "${params.outdir}/qc/library_complexity", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/library_complexity/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_complexity_command.sh" }

    input:
    tuple val(sample_id), path(bam), path(bai)   // S.filtered.bam (pre-duplicate-removal)
    path step_script                              // bin/library_complexity_frags.sh, or the project's `library_complexity_frags:` replacement

    output:
    path "${sample_id}.complexity.tsv", emit: tsv
    path "command.sh",                  emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --threads ${task.cpus} \\
        --sort-mem ${Math.max(1, (task.memory.toGiga() / 2) as int)}G \\
        ${params.library_complexity_frags_args ?: ''}
    """

    stub:
    """
    printf 'sample_id\\ttotal_fragments\\tdistinct_fragments\\tonce\\ttwice\\tNRF\\tPBC1\\tPBC2\\n${sample_id}\\t1000\\t900\\t850\\t40\\t0.9000\\t0.9444\\t21.2500\\n' > "${sample_id}.complexity.tsv"
    """
}
