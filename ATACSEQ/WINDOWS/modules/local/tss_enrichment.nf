include { runner; extraArgs } from './functions'

process TSS_ENRICHMENT {
    tag "${sample_id}"
    container 'quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0'
    // bigWig input: no cache 'deep'
    publishDir "${params.outdir}/qc/tss_enrichment", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/tss_enrichment/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_tss_command.sh" }

    input:
    tuple val(sample_id), path(bigwig)
    path tss_bed
    val flank           // tss_flank: bp on each side of the TSS (computeMatrix -b/-a)
    path step_script    // bin/qc_tss_enrichment.sh, or the project's `qc_tss_enrichment:` replacement

    output:
    path "${sample_id}.tss_enrichment.tsv", emit: tsv
    path "${sample_id}.tss_matrix.gz",      emit: matrix
    path "${sample_id}.tss_profile.svg",    emit: svg
    path "command.sh",                      emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bigwig "${bigwig}" \\
        --tss-bed "${tss_bed}" \\
        --flank ${flank} \\
        --threads ${task.cpus} \\
        ${extraArgs('qc_tss_enrichment')}
    """

    stub:
    """
    printf 'sample_id\\tn_regions\\twindow_up\\twindow_down\\tbin_bp\\ttss_signal\\tflank_mean\\tenrichment\\tprofile_max\\tprofile_min\\tmax_over_min\\n${sample_id}\\t1000\\t${flank}\\t${flank}\\t50\\t2.0\\t1.0\\t2.0\\t2.0\\t0.9\\t2.2\\n' > "${sample_id}.tss_enrichment.tsv"
    touch "${sample_id}.tss_matrix.gz"
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "${sample_id}.tss_profile.svg"
    """
}
