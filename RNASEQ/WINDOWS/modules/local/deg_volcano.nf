include { runner; extraArgs } from './functions'
include { shq; rbool } from './deg_functions'

process DEG_VOLCANO {
    tag "${model}"
    container '/project/rephilllab/nf_containers/omics-r_1.0.sif'   // R 4 + DESeq2 etc., built from the omics Docker image R env lock
    cache 'deep'
    publishDir "${params.outdir}/deg", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/deg/${model}/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "volcano_command.sh" }

    input:
    tuple val(model), path(results)   // the model's results_<comparison>.csv files (DEG_RESULTS), staged flat
    path comparisons    // comparisons.tsv (DEG_PREPARE)
    path common         // bin/deg_common.R, or the project's `deg_common` replacement
    val alpha           // params.deg_alpha
    val lfc             // params.deg_lfc
    val sval_threshold  // params.deg_sval_threshold
    val shrink          // 'true' unless params.deg_shrink_type is none: also draw the shrunk-LFC volcano
    val color_down      // params.deg_color_down
    val color_up        // params.deg_color_up
    val color_ns        // params.deg_color_ns
    val label_genes     // params.deg_label_genes
    val xlim            // params.deg_volcano_xlim
    val ylim            // params.deg_volcano_ylim
    val common_xlim     // params.deg_volcano_common_xlim
    val common_ylim     // params.deg_volcano_common_ylim
    val combined_layout // params.deg_volcano_combined_layout: stacked | side
    val pdfs            // params.deg_run_volcano_pdfs (true | false)
    val combined        // params.deg_run_volcano_combined (true | false)
    path step_script    // bin/deg_volcano.R, or the project's `deg_volcano` script replacement

    output:
    path "${model}/*/volcano_*.svg",                      optional: true, emit: individual   // volcano_all_genes_<c>.svg, volcano_shrunk_all_genes_<c>.svg
    path "${model}/volcano_plots_*.pdf",                  optional: true, emit: pdfs         // only with pdfs true
    path "${model}/volcano_combined_*_comparisons_*.svg", optional: true, emit: combined     // only with combined true and >= 2 comparisons with a combined_order
    path "command.sh",                                    emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --common="${common}" \\
        --model=${shq(model)} \\
        --comparisons="${comparisons}" \\
        --alpha=${alpha} \\
        --lfc=${lfc} \\
        --sval-threshold=${sval_threshold} \\
        --shrink=${rbool(shrink)} \\
        --color-down=${shq(color_down)} \\
        --color-up=${shq(color_up)} \\
        --color-ns=${shq(color_ns)} \\
        --label-genes=${shq(label_genes)} \\
        --xlim=${shq(xlim)} \\
        --ylim=${shq(ylim)} \\
        --common-xlim=${shq(common_xlim)} \\
        --common-ylim=${shq(common_ylim)} \\
        --combined-layout=${shq(combined_layout)} \\
        --individual=TRUE \\
        --pdfs=${rbool(pdfs)} \\
        --combined=${rbool(combined)} \\
        --out-dir=. \\
        ${extraArgs('deg_volcano')}
    """

    stub:
    // The comparison names come from the staged results file names. The combined figure is written when >= 2 comparisons of
    // the model have a combined_order (column read by name from the stub comparisons.tsv).
    def comps = [results].flatten().collect { f -> f.name.replaceFirst(/^results_/, '').replaceFirst(/\.csv$/, '') }
    def svg = "echo '<svg xmlns=\"http://www.w3.org/2000/svg\"/>'"
    def cmds = []
    comps.each { c ->
        cmds << "mkdir -p \"${model}/${c}\""
        cmds << "${svg} > \"${model}/${c}/volcano_all_genes_${c}.svg\""
        if (rbool(shrink) == 'TRUE') cmds << "${svg} > \"${model}/${c}/volcano_shrunk_all_genes_${c}.svg\""
    }
    if (rbool(pdfs) == 'TRUE') cmds << "touch \"${model}/volcano_plots_yscaled.pdf\" \"${model}/volcano_plots_xscaled.pdf\" \"${model}/volcano_plots_freescale.pdf\""
    if (rbool(combined) == 'TRUE') {
        cmds << "n=\$(awk -F'\\t' -v m=${shq(model)} 'NR == 1 { for (i = 1; i <= NF; i++) h[\$i] = i; next } \$(h[\"model\"]) == m && \$(h[\"combined_order\"]) != \"\" { n++ } END { print n + 0 }' \"${comparisons}\")"
        cmds << "if [ \"\$n\" -ge 2 ]; then ${svg} > \"${model}/volcano_combined_\${n}_comparisons_${combined_layout}.svg\"; fi"
    }
    """
    mkdir -p "${model}"
    ${cmds.join('\n    ')}
    """
}
