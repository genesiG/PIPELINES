include { runner; extraArgs } from './functions'

process QC_PEAK_STATS {
    tag "peak_stats"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    cache 'deep'
    publishDir "${params.outdir}/qc/peaks", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/peaks/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "peak_stats_command.sh" }

    input:
    path peaks          // all peak files (*.bed, *.narrowPeak, *.broadPeak), staged flat; names must be unique
    path step_script    // bin/qc_peak_stats.py, or the project's `qc_peak_stats:` replacement

    output:
    path "peak_stats.tsv", emit: table
    path "command.sh",     emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --peaks-dir . \\
        --out peak_stats.tsv \\
        ${extraArgs('qc_peak_stats')}
    """

    stub:
    """
    printf 'peak_set\\tfile\\tn_peaks\\ttotal_bp\\tmin_width\\tmedian_width\\tmean_width\\tmax_width\\tiqr_width\\tn_chroms\\tpct_autosomal\\nstub1\\tstub1.bed\\t100\\t50000\\t200\\t500.0\\t500.0\\t900\\t100.0\\t22\\t95.00\\n' > peak_stats.tsv
    """
}
