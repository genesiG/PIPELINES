/*
 * Peak annotation of the peak sets against the GTF, in R (bin/peak_annotation.R). Output prefix is fixed to "annotation".
 * labels / bed_names: comma-separated, same order as the staged beds. colors: "label=#hex,label2=#hex" (optional).
 */
include { runner } from './functions'

process PEAK_ANNOTATION {
    tag 'annotation'
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    publishDir { "${params.outdir}/annotation" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/annotation/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "annotation_peak_annotation_command.sh" }

    input:
    tuple val(labels), val(bed_names), path(beds)
    path gtf
    val colors              // optional, "label=#hex,label2=#hex"; empty = default colours
    val kegg                // true | false
    val tss_region          // bp each side of the TSS
    path step_script        // bin/peak_annotation.R

    output:
    path "annotation.region_metrics.tsv", emit: region_metrics
    path "annotation.summary_stats.tsv", emit: summary
    path "annotation.overlap_matrix.tsv", emit: overlap
    path "annotation.*.svg", emit: svg
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --bed-files="${bed_names}" \\
        --labels="${labels}" \\
        --gtf="${gtf}" \\
        --out-prefix=annotation \\
        --tss-region=${tss_region} \\
        --kegg=${kegg} \\
        ${colors ? "--colors='${colors}'" : ''} \\
        --cpus=${task.cpus} \\
        ${params.peak_annotation_args ?: ''}
    """

    stub:
    """
    touch annotation.region_metrics.tsv annotation.summary_stats.tsv annotation.overlap_matrix.tsv
    for n in annotation_bar dist_to_tss kegg chromosome_coverage; do
        echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "annotation.\${n}.svg"
    done
    IFS=',' read -ra LB <<< "${labels}"
    for l in "\${LB[@]}"; do
        slug="\$(printf %s "\$l" | tr -c 'A-Za-z0-9_' '_')"
        echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "annotation.annotation_pie_\${slug}.svg"
    done
    """
}
