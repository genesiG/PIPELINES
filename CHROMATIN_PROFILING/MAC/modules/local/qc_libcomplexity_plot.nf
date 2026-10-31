include { runner } from './functions'

process QC_LIBCOMPLEXITY_PLOT {
    tag "library_complexity_plot"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    cache 'deep'
    publishDir "${params.outdir}/qc/library_complexity", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/library_complexity/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "library_complexity_plot_command.sh" }

    input:
    path table          // library_complexity.tsv
    path theme          // bin/my_theme_py.py
    path step_script    // bin/qc_library_complexity_plot.py, or the project's `qc_library_complexity_plot:` replacement

    output:
    path "library_complexity.svg", emit: svg
    path "command.sh",            emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --table "${table}" \\
        --theme "${theme}" \\
        --out library_complexity.svg \\
        ${params.qc_library_complexity_plot_args ?: ''}
    """

    stub:
    """
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > library_complexity.svg
    """
}
