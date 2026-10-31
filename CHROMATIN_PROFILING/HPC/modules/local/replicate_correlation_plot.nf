/*
 * Plots of the replicate concordance step (bin/qc_replicate_correlation_plot.R): replicate_correlation.svg and
 * replicate_scatter.svg (+ .p2.svg, ...). Reads the two tables of REPLICATE_CORRELATION; the BIN_COUNTS counts only for the scatter.
 */
include { runner } from './functions'

process REPLICATE_CORRELATION_PLOT {
    tag "replicate_correlation_plot"
    cache 'deep'
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    publishDir "${params.outdir}/qc/replicate_correlation", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/replicate_correlation/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "replicate_correlation_plot_command.sh" }

    input:
    path table              // replicate_correlation.tsv
    path pairs_all          // replicate_correlation_pairs_all.tsv
    tuple val(bin_size), path(counts), path(libsizes)
    path theme              // bin/my_theme.R
    path step_script        // bin/qc_replicate_correlation_plot.R, or the project's replacement

    output:
    path "replicate_*.svg", emit: svg
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --table=${table} \\
        --pairs-all=${pairs_all} \\
        --counts=${counts} \\
        --libsizes=${libsizes} \\
        --theme=${theme} \\
        --outdir=. \\
        ${params.qc_replicate_correlation_plot_args ?: ''}
    """

    stub:
    """
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > replicate_correlation.svg
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > replicate_scatter.svg
    """
}
