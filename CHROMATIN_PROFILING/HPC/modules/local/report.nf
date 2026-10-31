/*
 * Self-contained HTML QC report assembled from the QC tables and SVG plots
 * only (no recomputation). Downstream of every table and plot process.
 */
include { runner } from './functions'

process REPORT {
    tag "${project_name}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    cache 'deep'
    publishDir "${params.outdir}/qc/report", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/report/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${project_name}_report_command.sh" }

    input:
    path tables         // all QC tables (*.tsv), staged flat; absent optional tables simply drop their section
    path plots          // all QC SVGs, staged flat
    val project_name
    path step_script    // bin/build_report.py, or the project's `build_report:` replacement

    output:
    path "${project_name}_qc_report.html", emit: report
    path "command.sh",                     emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --project "${project_name}" \\
        --tables-dir . \\
        --plots-dir . \\
        --outdir . \\
        ${params.build_report_args ?: ''}
    """

    stub:
    """
    touch "${project_name}_qc_report.html"
    """
}
