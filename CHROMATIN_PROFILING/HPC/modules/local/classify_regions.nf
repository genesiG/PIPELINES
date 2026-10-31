/*
 * Classifies one set of named regions into retained / lost / unclassified from the ref and treat BAMs (bin/classify_regions.R).
 * Writes <id>.retained.bed, <id>.lost.bed, <id>.unclassified.bed, <id>.per_region.tsv.gz and <id>.summary.tsv.
 * *_ids / *_names: comma-separated; names are the staged BAM names in the order of the matching ids list.
 * scale_factors: the norm_run scale-factor TSV. The regions file is passed by its staged path.
 */
include { runner; maxFrag } from './functions'

process CLASSIFY_REGIONS {
    tag "${classification}"
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    publishDir { "${params.outdir}/classify/${classification}" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/classify/${classification}/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${classification}_classify_regions_command.sh" }

    input:
    tuple val(classification), path(regions, stageAs: 'regions_in/*'), val(ref_ids), val(ref_names), path(ref_bams), path(ref_bais), val(treat_ids), val(treat_names), path(treat_bams), path(treat_bais), val(norm_run), path(scale_factors), val(lfc_retained), val(lfc_lost), val(fdr)
    path step_script   // bin/classify_regions.R

    output:
    tuple val(classification), path("${classification}.retained.bed"), path("${classification}.lost.bed"), path("${classification}.unclassified.bed"), emit: beds
    tuple val(classification), path("${classification}.per_region.tsv.gz"), emit: per_region
    tuple val(classification), path("${classification}.summary.tsv"), emit: summary
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --regions=${regions} \\
        --ref-bams=${ref_names} \\
        --treat-bams=${treat_names} \\
        --ref-ids=${ref_ids} \\
        --treat-ids=${treat_ids} \\
        --scale-factors=${scale_factors} \\
        --id=${classification} \\
        --lfc-retained=${lfc_retained} \\
        --lfc-lost=${lfc_lost} \\
        --fdr=${fdr} \\
        --paired-end=${params.paired_end.toString().toUpperCase()} \\
        --minq=${params.mapq} \\
        --dedup=${params.remove_duplicates.toString().toUpperCase()} \\
        --max-frag=${maxFrag()} \\
        --norm-run=${norm_run} \\
        --cpus=${task.cpus} \\
        --outdir=. \\
        ${params.classify_regions_args ?: ''}
    """

    stub:
    """
    printf 'chr1\\t100\\t200\\t%s\\t0\\t+\\n' "${classification}_retained" > "${classification}.retained.bed"
    printf 'chr1\\t300\\t400\\t%s\\t0\\t+\\n' "${classification}_lost" > "${classification}.lost.bed"
    : > "${classification}.unclassified.bed"
    printf 'chr\\tstart0\\tend\\tname\\tlogFC\\tlogCPM\\tPValue\\tFDR\\tclass\\nchr1\\t100\\t200\\tstub\\t2.5\\t6.0\\t0.001\\t0.01\\tretained\\n' | gzip > "${classification}.per_region.tsv.gz"
    printf 'class\\tn\\tMb\\nretained\\t1\\t0.0001\\nlost\\t1\\t0.0001\\nunclassified\\t0\\t0\\ntotal\\t2\\t0.0002\\n' > "${classification}.summary.tsv"
    printf 'norm_run\\t%s\\nlfc_retained\\t%s\\nlfc_lost\\t%s\\nfdr\\t%s\\nn_ref\\t%s\\nn_treat\\t%s\\n' "${norm_run}" "${lfc_retained}" "${lfc_lost}" "${fdr}" "${ref_ids.toString().split(',').size()}" "${treat_ids.toString().split(',').size()}" >> "${classification}.summary.tsv"
    """
}
