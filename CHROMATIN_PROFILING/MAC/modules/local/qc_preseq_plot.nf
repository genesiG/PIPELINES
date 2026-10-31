include { runner } from './functions'

process QC_PRESEQ_PLOT {
    tag "preseq_ccurve"
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    cache 'deep'
    publishDir "${params.outdir}/qc/preseq", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/preseq/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "preseq_plot_command.sh" }

    input:
    path curves         // all S.preseq_ccurve.txt, staged flat
    path theme          // bin/my_theme.R
    path step_script    // bin/qc_preseq_plot.R, or the project's `qc_preseq_plot:` replacement

    output:
    path "preseq_ccurve.svg", emit: svg
    path "command.sh",        emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --curves-dir=. \\
        --theme="${theme}" \\
        ${params.qc_preseq_plot_args ?: ''}
    """

    stub:
    """
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > preseq_ccurve.svg
    """
}
