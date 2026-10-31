include { runner } from './functions'

process PRESEQ {
    tag "${sample_id}"
    container 'quay.io/biocontainers/preseq:2.0.2--gsl1.16_0'
    // BAM input: no cache 'deep'
    publishDir "${params.outdir}/qc/preseq", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/preseq/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_preseq_command.sh" }

    input:
    tuple val(sample_id), path(bam), path(bai)   // S.filtered.bam (pre-duplicate-removal)
    path step_script                              // bin/preseq_ccurve.sh, or the project's `preseq_ccurve:` replacement

    output:
    path "${sample_id}.preseq_ccurve.txt", emit: curve
    path "command.sh",                     emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        ${params.preseq_ccurve_args ?: ''}
    """

    stub:
    """
    printf 'TOTAL_READS\\tDISTINCT_READS\\n0\\t0\\n1000\\t900\\n2000\\t1700\\n' > "${sample_id}.preseq_ccurve.txt"
    """
}
