include { runner; extraArgs } from './functions'
include { shq; optArg } from './deg_functions'

process DEG_HEATMAP {
    tag "${hm.name}"
    container '/project/rephilllab/nf_containers/omics-r_1.0.sif'   // R 4 + DESeq2 etc., built from the omics Docker image R env lock
    cache 'deep'
    publishDir "${params.outdir}/deg", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/deg/${hm.model}/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "heatmap_${hm.name}_command.sh" }

    input:
    // hm: one row of deg_heatmaps with the defaults filled in (name, model, comparison, clust_method, cutree_rows, col_swap,
    // palette_n, fontsize, fontsize_row, width, height, title). Files: the comparison's results table and the model's
    // rlogCounts.rds / vstCounts.rds (both staged flat), and the row's gene set (an empty list = none).
    tuple val(hm), path(results), path(counts), path(gene_set, stageAs: 'gene_set/*')
    path label_rules    // params.deg_label_rules, or an empty list when unset (nothing staged)
    path common         // bin/deg_common.R, or the project's `deg_common` replacement
    val color_down      // params.deg_color_down
    val color_up        // params.deg_heatmap_color_up
    val color_mid       // params.deg_color_mid
    path step_script    // bin/deg_heatmap.R, or the project's `deg_heatmap` script replacement

    output:
    path "${hm.model}/heatmap_${hm.name}.svg", optional: true, emit: svg   // not written (exit 0) when fewer than 2 genes pass
    path "command.sh",                         emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --common="${common}" \\
        --name=${shq(hm.name)} \\
        --model=${shq(hm.model)} \\
        --comparison=${shq(hm.comparison)} \\
        --counts="${counts}" \\
        --results="${results}" \\
        ${optArg('gene-set', gene_set)} \\
        ${optArg('label-rules', label_rules)} \\
        --clust-method=${shq(hm.clust_method)} \\
        --cutree-rows=${shq(hm.cutree_rows)} \\
        --col-swap=${shq(hm.col_swap)} \\
        --palette-n=${hm.palette_n} \\
        --fontsize=${hm.fontsize} \\
        --fontsize-row=${hm.fontsize_row} \\
        --width=${hm.width} \\
        --height=${hm.height} \\
        --title=${shq(hm.title)} \\
        --color-down=${shq(color_down)} \\
        --color-up=${shq(color_up)} \\
        --color-mid=${shq(color_mid)} \\
        --out-dir=. \\
        ${extraArgs('deg_heatmap')}
    """

    stub:
    """
    mkdir -p ${shq(hm.model)}
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > ${shq("${hm.model}/heatmap_${hm.name}.svg")}
    """
}
