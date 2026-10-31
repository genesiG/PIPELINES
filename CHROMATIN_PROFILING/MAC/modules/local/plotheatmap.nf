/*
 * deepTools plotHeatmap for one heatmap (bin/heatmap_plot.sh). Lists are comma-separated; optional arguments are
 * passed only when non-empty. labels may contain a literal backslash-n, so they go through a single-quoted shell
 * variable (the script turns \n into a newline). Ref/x labels, size, plot type, summary average and legend location
 * come from params.heatmap_*; y_label stays a task input.
 */
include { runner } from './functions'

process PLOTHEATMAP {
    tag "${heatmap}"
    cache 'deep'
    container 'quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0'
    publishDir { "${params.outdir}/heatmaps/plots" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/heatmaps/plots/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${heatmap}_plotheatmap_command.sh" }

    input:
    tuple val(heatmap), path(matrix), val(labels), val(colors), val(z_min), val(z_max), val(y_max), val(region_labels), val(sort_idx)
    val y_label
    path step_script   // bin/heatmap_plot.sh

    output:
    path "${heatmap}.heatmap.svg", emit: svg
    path "${heatmap}.sorted_regions.bed", emit: sorted
    path "command.sh", emit: command_sh

    script:
    """
    LABELS='${labels.toString().replace("'", "'\\''")}'
    ${runner(step_script)} "${step_script}" \\
        --matrix "${matrix}" \\
        --out-prefix "${heatmap}" \\
        --labels "\$LABELS" \\
        --colors "${colors}" \\
        ${z_min ? "--z-min '${z_min}'" : ''} \\
        ${z_max ? "--z-max '${z_max}'" : ''} \\
        ${y_max ? "--y-max '${y_max}'" : ''} \\
        ${region_labels ? "--region-labels '${region_labels}'" : ''} \\
        ${sort_idx ? "--sort-using '${sort_idx}'" : ''} \\
        --y-label "${y_label}" \\
        --ref-label "${params.heatmap_ref_label}" \\
        --x-label "${params.heatmap_x_label}" \\
        --height ${params.heatmap_height} \\
        --width ${params.heatmap_width} \\
        --plot-type ${params.heatmap_plot_type} \\
        --average-type-summary ${params.heatmap_average_type_summary} \\
        --legend-location ${params.heatmap_legend_location} \\
        --threads ${task.cpus} \\
        ${params.heatmap_plot_args ?: ''}
    """

    stub:
    """
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "${heatmap}.heatmap.svg"
    touch "${heatmap}.sorted_regions.bed"
    """
}
