include { runner } from './functions'

process QC_FRAGMENT_SIZES_PLOT {
    tag "fragment_sizes_plots"
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    cache 'deep'
    publishDir "${params.outdir}/qc/fragment_sizes", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/fragment_sizes/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "fragment_sizes_plot_command.sh" }

    input:
    path summary        // fragment_sizes.tsv
    path distribution   // fragment_size_distribution.tsv
    path theme          // bin/my_theme.R
    path step_script    // bin/qc_fragment_sizes_plot.R, or the project's `qc_fragment_sizes_plot:` replacement

    output:
    path "fragment_sizes_*.svg", emit: svg
    path "command.sh",           emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --summary="${summary}" \\
        --distribution="${distribution}" \\
        --theme="${theme}" \\
        ${params.qc_fragment_sizes_plot_args ?: ''}
    """

    stub:
    """
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > fragment_sizes_density.svg
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > fragment_sizes_boxplot.svg
    """
}
