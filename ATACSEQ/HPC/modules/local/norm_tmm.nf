/*
 * csaw/edgeR TMM scale factors for one TMM group (bin/norm_tmm.R, written by the statistician).
 * sample_ids: comma-separated sample ids, same order as bams (the staged BAM names are computed here).
 * Not cached: the input is BAMs.
 */
include { runner; extraArgs } from './functions'

process NORM_TMM {
    tag "${group}"
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    publishDir "${params.outdir}/scale_factors", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/scale_factors/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${group}_norm_tmm_command.sh" }

    input:
    tuple val(group), val(sample_ids), path(bams), path(bais)   // bams/bais: final (or tn5) BAMs of the group, staged flat
    val bin_size            // tmm_bin_size
    val mapq                // mapq
    val max_frag            // tmm_max_frag: 'auto' or a positive integer
    val bigwig_units        // bigwig_units: none | cpm | rpkm (sets bigwig_multiplier)
    path step_script        // bin/norm_tmm.R, or the project's `norm_tmm:` replacement

    output:
    tuple val(group), path("${group}.scale_factors.tsv"), emit: sf
    path "${group}.ma_plot.{png,pdf}",                    emit: plots
    path "${group}.max_frag.txt",                         emit: max_frag
    path "command.sh",                                    emit: command_sh

    script:
    def bam_names = [bams].flatten().collect { it.name }.join(',')
    """
    ${runner(step_script)} "${step_script}" \\
        --bams=${bam_names} \\
        --sample-ids=${sample_ids} \\
        --norm-run=${group} \\
        --bin-size=${bin_size} \\
        --mapq=${mapq} \\
        --dedup=FALSE \\
        --max-frag=${max_frag} \\
        --pe=both \\
        --bigwig-units=${bigwig_units} \\
        --out=${group}.scale_factors.tsv \\
        ${extraArgs('norm_tmm')}
    """

    stub:
    """
    printf 'sample_id\\tmethod\\traw_factor\\tbigwig_multiplier\\tedger_norm_factor\\tqc_pass\\tnote\\n' > "${group}.scale_factors.tsv"
    for id in \$(echo "${sample_ids}" | tr ',' ' '); do
        printf '%s\\ttmm\\t1\\t1\\t1\\tTRUE\\tstub\\n' "\$id" >> "${group}.scale_factors.tsv"
    done
    touch "${group}.ma_plot.png" "${group}.ma_plot.pdf"
    echo 500 > "${group}.max_frag.txt"
    """
}
