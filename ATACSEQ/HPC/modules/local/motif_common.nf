/*
 * Motifs significant in several comparisons, one region set (bin/motif_common.R; legacy find_common_motifs.Rmd).
 * names: comma-separated comparison names, same order as tables (labels = names with motif_common_label_strip removed).
 * summaries: MOTIF_ENRICHMENT NAME.motif_summary.csv files of this region set, or assets/NO_FILE when there are none.
 * Published to motif/common/REGION/.
 */
include { runner; extraArgs } from './functions'

process MOTIF_COMMON {
    tag "${region}"
    container '/project/rephilllab/nf_containers/omics-r_1.0.sif'
    publishDir path: { "${params.outdir}/motif/common/${region}" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir path: { "${params.outdir}/motif/common/${region}/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${region}_motif_common_command.sh" }

    input:
    tuple val(region), val(names), path(tables), path(summaries)   // tables: NAME_motif_analysis.csv of this region set
    val label_strip                 // motif_common_label_strip: regex removed from the names ('' = none)
    val padj                        // motif_common_padj
    val min_lines                   // motif_common_min_lines
    val plot_font                   // plot_font: Helvetica (font family of the plot)
    val plot_format                 // plot_format: svg | pdf
    path step_script                // bin/motif_common.R, or the project's `motif_common:` replacement

    output:
    path "common_significant_motifs.csv",    emit: table
    path "venn_common_motifs.{svg,pdf}",     optional: true, emit: venn
    path "motif_analysis_summary.csv",       optional: true, emit: summary
    path "command.sh",                       emit: command_sh

    script:
    def table_names   = [tables].flatten().collect { it.name }.join(',')
    def summary_names = [summaries].flatten().findAll { it.name != 'NO_FILE' }.collect { it.name }.join(',')
    """
    ${runner(step_script)} "${step_script}" \\
        --inputs=${table_names} \\
        --labels=${names} \\
        --label-strip='${label_strip ?: ''}' \\
        --padj=${padj} \\
        --min-lines=${min_lines} \\
        --summaries=${summary_names} \\
        --out-prefix=. \\
        --font=${plot_font} \\
        --plot-format=${plot_format} \\
        ${extraArgs('motif_common')}
    """

    stub:
    """
    printf 'motif.id,motif.name,label,n_cell_lines_significant\\nMA0004.1,ARNT,ARNT,2\\n' > common_significant_motifs.csv
    printf '<svg xmlns="http://www.w3.org/2000/svg" width="576pt" height="504pt"><text font-family="Helvetica">venn</text></svg>\\n' > venn_common_motifs.svg
    printf 'comparison,region,total_peaks\\n' > motif_analysis_summary.csv
    """
}
