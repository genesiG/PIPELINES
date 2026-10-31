include { runner; extraArgs } from './functions'
include { shq; rbool; optArg } from './deg_functions'

process DEG_QC {
    tag "qc"
    container '/project/rephilllab/nf_containers/omics-r_1.0.sif'   // R 4 + DESeq2 etc., built from the omics Docker image R env lock
    cache 'deep'
    publishDir "${params.outdir}/deg/qc", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/deg/qc/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "qc_command.sh" }

    input:
    path counts         // counts_edited.rds (DEG_PREPARE)
    path metadata       // sample_metadata_edited.rds (DEG_PREPARE)
    path comparisons    // comparisons.tsv (DEG_PREPARE)
    path colors         // params.deg_colors, or an empty list when unset (nothing staged)
    path common         // bin/deg_common.R, or the project's `deg_common` replacement
    val color_by        // params.deg_color_by
    val barplot         // params.deg_run_barplot (true | false)
    val mds             // params.deg_run_mds
    val correlation     // params.deg_run_correlation
    val mds_input       // params.deg_mds_input: counts | logcpm
    val mds_label_by    // params.deg_mds_label_by
    val mds_legend_title  // params.deg_mds_legend_title
    val legacy_legend_col // params.deg_legacy_mds_legend_col
    val cor_method      // params.deg_cor_method
    path step_script    // bin/deg_qc.R, or the project's `deg_qc` script replacement

    output:
    path "counts_per_index.pdf",   optional: true, emit: counts_pdf   // only with barplot true
    path "counts_per_sample.tsv",  optional: true, emit: counts_tsv   // only with barplot true
    path "pca.pdf",                optional: true, emit: pca          // only with mds true
    path "*_correlation_heatmap.svg", optional: true, emit: correlation  // only with correlation true (<cor_method>_correlation_heatmap.svg)
    path "command.sh",             emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \
        --common="${common}" \
        --counts="${counts}" \
        --metadata="${metadata}" \
        --comparisons="${comparisons}" \
        --color-by=${shq(color_by)} \
        ${optArg('colors', colors)} \
        --barplot=${rbool(barplot)} \
        --mds=${rbool(mds)} \
        --correlation=${rbool(correlation)} \
        --mds-input=${shq(mds_input)} \
        --mds-label-by=${shq(mds_label_by)} \
        --mds-legend-title=${shq(mds_legend_title)} \
        --legacy-mds-legend-col=${shq(legacy_legend_col)} \
        --cor-method=${shq(cor_method)} \
        --out-dir=. \
        ${extraArgs('deg_qc')}
    """

    stub:
    def plots = []
    if (rbool(barplot) == 'TRUE') plots << 'touch counts_per_index.pdf; echo stub > counts_per_sample.tsv'
    if (rbool(mds) == 'TRUE') plots << 'touch pca.pdf'
    if (rbool(correlation) == 'TRUE') plots << "echo '<svg xmlns=\"http://www.w3.org/2000/svg\"/>' > ${cor_method}_correlation_heatmap.svg"
    """
    ${plots ? plots.join('\n    ') : ':'}
    """
}
