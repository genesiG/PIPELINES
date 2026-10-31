include { runner } from './functions'

process TSS_ENRICHMENT {
    tag "${sample_id}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    // bigWig input: no cache 'deep'
    publishDir "${params.outdir}/qc/tss_enrichment", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/tss_enrichment/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_tss_command.sh" }

    input:
    tuple val(sample_id), path(bigwig)
    path tss_ref        // the GTF (params.gtf, preferred); the TSS BED only when no GTF is configured
    val  ref_kind       // 'gtf' or 'bed'
    path step_script    // bin/qc_tss_enrichment.py, or the project's `qc_tss_enrichment:` replacement

    output:
    path "${sample_id}.tss_enrichment.tsv", emit: tsv
    path "${sample_id}.tss_profile.svg",    emit: svg
    path "command.sh",                      emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bigwig "${bigwig}" \\
        --${ref_kind == 'gtf' ? 'gtf' : 'tss-bed'} "${tss_ref}" \\
        --threads ${task.cpus} \\
        ${params.qc_tss_enrichment_args ?: ''}
    """

    stub:
    """
    printf 'sample_id\\tn_regions\\twindow_up\\twindow_down\\tbin_bp\\ttss_signal\\tflank_mean\\tenrichment\\tprofile_max\\tprofile_min\\tmax_over_min\\n${sample_id}\\t1000\\t2000\\t2000\\t50\\t2.0\\t1.0\\t2.0\\t2.0\\t0.9\\t2.2\\n' > "${sample_id}.tss_enrichment.tsv"
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "${sample_id}.tss_profile.svg"
    """
}
