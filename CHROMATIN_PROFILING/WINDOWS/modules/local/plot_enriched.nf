/*
 * EnrichedHeatmap-style plot for one heatmap, in R (bin/heatmap_plot_enriched.R). Non-clustered when assignments is
 * empty; otherwise rows are split by cluster (assignments = NAME.cluster_assignments.csv from HEATMAP_CLUSTER).
 * labels, legend_titles and cluster_labels go through single-quoted shell variables (the script turns \n into a newline).
 */
include { runner } from './functions'

process PLOT_ENRICHED {
    tag "${heatmap}"
    cache 'deep'
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    publishDir { "${params.outdir}/heatmaps/plots" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/heatmaps/plots/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${heatmap}_plot_enriched_command.sh" }

    input:
    tuple val(heatmap), path(norm_mats), path(targets), val(sample_ids), val(labels), val(colors), val(legend_titles), val(z_min), val(z_max), val(y_max), val(window), val(cluster_labels), val(sort_by), path(assignments)
    path step_script   // bin/heatmap_plot_enriched.R

    output:
    path "${heatmap}.enriched_heatmap.svg", emit: svg
    path "command.sh", emit: command_sh

    script:
    """
    LABELS='${labels.toString().replace("'", "'\\''")}'
    LEGEND='${(legend_titles ?: '').toString().replace("'", "'\\''")}'
    CLUSTER_LABELS='${(cluster_labels ?: '').toString().replace("'", "'\\''")}'
    ${runner(step_script)} "${step_script}" \\
        --norm-mats=${norm_mats} \\
        --targets=${targets} \\
        --window=${window} \\
        --out-prefix=${heatmap} \\
        --names=${sample_ids} \\
        --labels="\$LABELS" \\
        --colors="${colors}" \\
        ${legend_titles ? '--legend-titles="$LEGEND"' : ''} \\
        ${z_min ? "--z-min=${z_min}" : ''} \\
        ${z_max ? "--z-max=${z_max}" : ''} \\
        ${y_max ? "--y-max=${y_max}" : ''} \\
        ${assignments ? "--assignments=${assignments}" : ''} \\
        ${assignments ? '--cluster-labels="$CLUSTER_LABELS"' : ''} \\
        ${assignments && sort_by ? "--sort-by=${sort_by}" : ''} \\
        ${params.heatmap_plot_enriched_args ?: ''}
    """

    stub:
    """
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "${heatmap}.enriched_heatmap.svg"
    """
}
