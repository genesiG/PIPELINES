/*
 * csaw differential-window peak calling, one task per peakset (bin/callpeaks_csaw.R).
 * *_names: comma-separated staged BAM names in the order of the matching *_ids list.
 */
include { runner; boolStr; maxFrag } from './functions'

process CSAW {
    tag "${peakset}"
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    publishDir { "${params.outdir}/peaks/csaw" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/peaks/csaw/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${peakset}_csaw_command.sh" }

    input:
    tuple val(peakset), val(target_ids), val(target_names), path(target_bams), path(target_bais), val(control_ids), val(control_names), path(control_bams), path(control_bais), val(small_filt), val(large_filt), val(merge_bp), val(max_width), val(fdr), val(lfc)
    path step_script   // bin/callpeaks_csaw.R

    output:
    path "${peakset}.peaks.bed", emit: peaks
    path "${peakset}.csaw_results.tsv.gz", emit: results
    path "${peakset}.mds.svg", emit: mds
    path "${peakset}.rds", emit: rds
    path "${peakset}.norm_sensitivity.tsv", emit: sensitivity
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --peakset=${peakset} \\
        --target-bams=${target_names} \\
        --control-bams=${control_names} \\
        --target-ids=${target_ids} \\
        --control-ids=${control_ids} \\
        --window=${params.csaw_window} \\
        --spacing=${params.csaw_spacing} \\
        --bin=${params.csaw_bin} \\
        --small-filt=${small_filt} \\
        --large-filt=${large_filt} \\
        --merge-bp=${merge_bp} \\
        --max-width=${max_width} \\
        --fdr=${fdr} \\
        --lfc=${lfc} \\
        --mapq=${params.mapq} \\
        --max-frag=${maxFrag()} \\
        --norm-bins=${params.csaw_norm_bins} \\
        --large-window=${params.csaw_large_window} \\
        --large-spacing=${params.csaw_large_spacing} \\
        --bg-max-enrich=${params.csaw_bg_max_enrich} \\
        --count-filter=${params.csaw_count_filter} \\
        --lfc-from=${params.csaw_lfc_from} \\
        --disjoint-regions=${boolStr(params.csaw_disjoint_regions)} \\
        --dedup=${params.remove_duplicates.toString().toUpperCase()} \\
        --cpus=${task.cpus} \\
        --out-prefix=${peakset} \\
        ${params.callpeaks_csaw_args ?: ''}
    """

    stub:
    """
    touch "${peakset}.peaks.bed"
    echo stub | gzip > "${peakset}.csaw_results.tsv.gz"
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "${peakset}.mds.svg"
    touch "${peakset}.rds"
    printf 'peakset\\tnorm_bins\\tlfc_from\\tdisjoint_regions\\tcount_filter\\tn_regions_csaw_merge\\tn_overlapping_regions_csaw_merge\\tn_peaks_csaw_merge\\tn_overlapping_peaks_csaw_merge\\n%s\\tstub\\t${params.csaw_lfc_from}\\tfalse\\t10\\t0\\t0\\t0\\t0\\n' "${peakset}" > "${peakset}.norm_sensitivity.tsv"
    """
}
