include { runner } from './functions'

process PICARD_ELC {
    tag "${sample_id}"
    container 'quay.io/biocontainers/picard:2.27.5--hdfd78af_0'
    // BAM input: no cache 'deep'
    publishDir "${params.outdir}/qc/library_complexity", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/library_complexity/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_elc_command.sh" }

    input:
    tuple val(sample_id), path(bam), path(bai)   // S.filtered.bam (pre-duplicate-removal)
    path step_script                              // bin/picard_elc.sh, or the project's `picard_elc:` replacement

    output:
    path "${sample_id}.elc.txt", emit: metrics
    path "command.sh",           emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --mem-gb ${task.memory.toGiga()} \\
        --threads ${task.cpus} \\
        ${params.picard_elc_args ?: ''}
    """

    stub:
    """
    printf '## METRICS CLASS\\tpicard.sam.DuplicationMetrics\\nLIBRARY\\tREAD_PAIRS_EXAMINED\\tUNIQUE_READ_PAIRS\\tPERCENT_DUPLICATION\\tESTIMATED_LIBRARY_SIZE\\nstub\\t1000\\t900\\t0.1\\t10000\\n' > "${sample_id}.elc.txt"
    """
}
