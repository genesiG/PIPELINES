/*
 * MA plot of one classification (bin/classify_regions_plot.R): logCPM vs logFC coloured by class, with threshold lines.
 * Reads only the per_region table and the summary written by CLASSIFY_REGIONS; one <id>.ma.svg.
 */
include { runner } from './functions'

process CLASSIFY_REGIONS_PLOT {
    tag "${classification}"
    cache 'deep'
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    publishDir { "${params.outdir}/classify/${classification}" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/classify/${classification}/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${classification}_classify_regions_plot_command.sh" }

    input:
    tuple val(classification), path(per_region), path(summary)
    path step_script   // bin/classify_regions_plot.R

    output:
    path "${classification}.ma.svg", emit: svg
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --per-region=${per_region} \\
        --summary=${summary} \\
        --id=${classification} \\
        --out=${classification}.ma.svg \\
        ${params.classify_regions_plot_args ?: ''}
    """

    stub:
    """
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "${classification}.ma.svg"
    """
}
