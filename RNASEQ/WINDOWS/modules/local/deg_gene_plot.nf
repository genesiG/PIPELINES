include { runner; extraArgs } from './functions'
include { shq; optArg } from './deg_functions'

process DEG_GENE_PLOT {
    tag "gene_plots"
    container '/project/rephilllab/nf_containers/omics-r_1.0.sif'   // R 4 + DESeq2 etc., built from the omics Docker image R env lock
    cache 'deep'
    publishDir "${params.outdir}/deg", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/deg/gene_plots/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "gene_plots_command.sh" }

    input:
    path counts         // rlogCounts.rds / vstCounts.rds of the first model (DEG_TRANSFORM)
    path metadata       // sample_metadata_edited.rds (DEG_PREPARE)
    path colors         // params.deg_colors, or an empty list when unset (nothing staged)
    path common         // bin/deg_common.R, or the project's `deg_common` replacement
    val genes           // params.deg_gene_plot_genes: comma list
    val x_by            // the metadata column on the x axis and for the colours (params.deg_gene_plot_x, else deg_color_by, else the first design variable)
    path step_script    // bin/deg_gene_plot.R, or the project's `deg_gene_plot` script replacement

    output:
    path "gene_plots/*.svg", emit: svg   // one per gene found in the table (genes absent from it are skipped)
    path "command.sh",       emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --common="${common}" \\
        --counts="${counts}" \\
        --metadata="${metadata}" \\
        --genes=${shq(genes)} \\
        --x-by=${shq(x_by)} \\
        --color-by=${shq(x_by)} \\
        ${optArg('colors', colors)} \\
        --out-dir=. \\
        ${extraArgs('deg_gene_plot')}
    """

    stub:
    def gene_list = genes.toString().split(',').collect { g -> g.trim() }.findAll { g -> g }
    """
    mkdir -p gene_plots
    for g in ${gene_list.collect { g -> shq(g) }.join(' ')}; do
        echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "gene_plots/\$g.svg"
    done
    """
}
