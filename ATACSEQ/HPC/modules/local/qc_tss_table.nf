include { runner; extraArgs } from './functions'

process QC_TSS_TABLE {
    tag "tss_enrichment"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    cache 'deep'
    publishDir "${params.outdir}/qc/tss_enrichment", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/tss_enrichment/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "tss_table_command.sh" }

    input:
    path tss_files      // all S.tss_enrichment.tsv, staged flat
    path samplesheet    // samplesheet.csv (metadata only: condition, replicate, tmm_group, peak_group)
    path step_script    // bin/qc_tss_table.py, or the project's `qc_tss_table:` replacement

    output:
    path "tss_enrichment.tsv", emit: table
    path "command.sh",         emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --input-dir . \\
        --samplesheet "${samplesheet}" \\
        --out tss_enrichment.tsv \\
        ${extraArgs('qc_tss_table')}
    """

    stub:
    """
    printf 'sample_id\\tcondition\\treplicate\\ttmm_group\\tpeak_group\\tn_regions\\twindow_up\\twindow_down\\tbin_bp\\ttss_signal\\tflank_mean\\tenrichment\\tprofile_max\\tprofile_min\\tmax_over_min\\nstub1\\tcond1\\t1\\tgrp1\\tpg1\\t1000\\t2000\\t2000\\t50\\t2.0\\t1.0\\t2.0\\t2.0\\t0.9\\t2.2\\n' > tss_enrichment.tsv
    """
}
