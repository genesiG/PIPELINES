/*
 * Optional post-hoc csaw filter-cutoff scan, one task per peakset (bin/csaw_cutoff_scan.R).
 * Same inputs as CSAW; never feeds thresholds automatically.
 */
include { runner; boolStr; maxFrag } from './functions'

process CSAW_CUTOFF_SCAN {
    tag "${peakset}"
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    publishDir { "${params.outdir}/peaks/csaw/cutoff_scan" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/peaks/csaw/cutoff_scan/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${peakset}_csaw_cutoff_scan_command.sh" }

    input:
    tuple val(peakset), val(target_ids), val(target_names), path(target_bams), path(target_bais), val(control_ids), val(control_names), path(control_bams), path(control_bais), val(small_filt), val(large_filt), val(merge_bp), val(max_width), val(fdr), val(lfc)
    path step_script   // bin/csaw_cutoff_scan.R

    output:
    path "${peakset}.cutoff_scan.tsv", emit: scan
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
        --small-grid=${params.csaw_scan_small_grid} \\
        --large-grid=${params.csaw_scan_large_grid} \\
        --out-prefix=${peakset} \\
        ${params.csaw_cutoff_scan_args ?: ''}
    """

    stub:
    """
    printf 'peakset\\tsmall_filt\\tlarge_filt\\tn_peaks\\tlfc_from\\tdisjoint_regions\\tcount_filter\\n' > "${peakset}.cutoff_scan.tsv"
    printf '%s\\t1\\t1\\t0\\tbest\\tfalse\\t10\\n' "${peakset}" >> "${peakset}.cutoff_scan.tsv"
    """
}
