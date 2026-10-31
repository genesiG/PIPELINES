include { runner; extraArgs } from './functions'

process QC_FRIP_TABLE {
    tag "frip"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    cache 'deep'
    publishDir "${params.outdir}/qc/frip", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/frip/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "frip_table_command.sh" }

    input:
    path frip_files     // all S.<peak_source>.frip.tsv, staged flat
    path samplesheet    // samplesheet.csv (metadata only: condition, replicate, tmm_group, peak_group)
    path step_script    // bin/qc_frip_table.py, or the project's `qc_frip_table:` replacement

    output:
    path "frip.tsv",   emit: table
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --input-dir . \\
        --samplesheet "${samplesheet}" \\
        --out frip.tsv \\
        ${extraArgs('qc_frip_table')}
    """

    stub:
    """
    printf 'sample_id\\tcondition\\treplicate\\ttmm_group\\tpeak_group\\tpeak_source\\tpeaks_file\\tn_peaks\\ttotal_mapped_alignments\\treads_in_peaks\\tfrip\\nstub1\\tcond1\\t1\\tgrp1\\tpg1\\tmacs3\\tpeaks.bed\\t100\\t1000\\t200\\t0.2\\n' > frip.tsv
    """
}
