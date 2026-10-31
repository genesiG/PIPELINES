include { runner } from './functions'

process CORRELATION {
    tag "correlation"
    container 'quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0'
    // BAM inputs: no cache 'deep'
    publishDir "${params.outdir}/qc/correlation", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/correlation/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "correlation_command.sh" }

    input:
    path bams           // all non-skipped S.final.bam
    path bais           // their S.final.bam.bai (not in the contract's input list; deepTools needs the indexes next to the BAMs)
    val labels          // comma-separated labels (sample ids), same order as bams
    path step_script    // bin/qc_correlation.sh, or the project's `qc_correlation:` replacement

    output:
    path "correlation_matrix.npz",    emit: npz
    path "correlation_heatmap.svg",   emit: svg
    path "correlation_matrix.tsv",    emit: tsv
    path "command.sh",                emit: command_sh

    script:
    def bam_list = (bams instanceof Collection ? bams : [bams]).collect { it.name }.join(',')
    """
    ${runner(step_script)} "${step_script}" \\
        --bams "${bam_list}" \\
        --labels "${labels}" \\
        --threads ${task.cpus} \\
        ${params.qc_correlation_args ?: ''}
    """

    stub:
    """
    touch correlation_matrix.npz
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > correlation_heatmap.svg
    printf "#plotCorrelation --outFileCorMatrix\\n\\tstub1\\tstub2\\nstub1\\t1.0\\t0.9\\nstub2\\t0.9\\t1.0\\n" > correlation_matrix.tsv
    """
}
