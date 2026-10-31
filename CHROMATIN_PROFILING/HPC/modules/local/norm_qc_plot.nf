/*
 * Normalisation QC MA plots for one norm_run (bin/norm_qc_plot.R). Reads only tables: the pairs table of NORM_QC (the line
 * of every method is read from it), the wide scale-factor table (reason text of unavailable methods) and the BIN_COUNTS
 * counts + libsizes (the cloud). Pages: <norm_run>.ma.svg, <norm_run>.ma.p2.svg, ... (norm_qc_panels_per_page panels each).
 */
include { runner } from './functions'

process NORM_QC_PLOT {
    tag "${norm_run}"
    cache 'deep'
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    publishDir { "${params.outdir}/norm_qc" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/norm_qc/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${norm_run}_norm_qc_plot_command.sh" }

    input:
    tuple val(norm_run), path(pairs), path(scale_factors_all), path(counts), path(libsizes)
    path step_script   // bin/norm_qc_plot.R

    output:
    path "${norm_run}.ma*.svg", emit: svg
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --pairs=${pairs} \\
        --scale-factors-all=${scale_factors_all} \\
        --counts=${counts} \\
        --libsizes=${libsizes} \\
        --panels-per-page=${params.norm_qc_panels_per_page} \\
        --out-prefix=${norm_run} \\
        ${params.norm_qc_plot_args ?: ''}
    """

    stub:
    """
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "${norm_run}.ma.svg"
    """
}
