include { runner } from './functions'

process QC_ALIGNMENT_PLOT {
    tag "alignment_plots"
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    cache 'deep'
    publishDir "${params.outdir}/qc/alignment", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/alignment/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "alignment_plot_command.sh" }

    input:
    path table          // alignment_stats.tsv
    path theme          // bin/my_theme.R
    path step_script    // bin/qc_alignment_plot.R, or the project's `qc_alignment_plot:` replacement

    output:
    path "alignment_*.svg", emit: svg
    path "command.sh",      emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --table="${table}" \\
        --theme="${theme}" \\
        ${params.qc_alignment_plot_args ?: ''}
    """

    stub:
    """
    for n in read_composition read_composition_pct pct_spikein n_spikein pct_chrM duplicates final_reads; do
        echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "alignment_\${n}.svg"
    done
    """
}
