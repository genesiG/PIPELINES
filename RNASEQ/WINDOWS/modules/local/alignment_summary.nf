include { runner; extraArgs } from './functions'

process ALIGNMENT_SUMMARY {
    tag "alignment_summary"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    cache 'deep'
    publishDir "${params.outdir}/qc/hisat2", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/hisat2/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "alignment_summary_command.sh" }

    input:
    path '*'           // all S.hisat2_summary.txt files (ALIGN), collected and staged flat
    val min_pct        // params.alignment_min_pct: samples below this overall rate are flagged (warning only)
    path step_script   // bin/alignment_summary.py, or the project's `alignment_summary` script replacement

    output:
    path "alignment_summary.tsv", emit: table
    path "command.sh",            emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --logs-dir . \\
        --min-pct ${min_pct} \\
        --out alignment_summary.tsv \\
        ${extraArgs('alignment_summary')}
    """

    stub:
    """
    printf 'sample_id\\ttotal_pairs\\taligned_concordantly_0\\taligned_concordantly_0_pct\\taligned_concordantly_1\\taligned_concordantly_1_pct\\taligned_concordantly_multi\\taligned_concordantly_multi_pct\\taligned_discordantly_1\\taligned_discordantly_1_pct\\toverall_alignment_rate\\tbelow_min_pct\\n' > alignment_summary.tsv
    printf 'stub1\\t10000\\t1000\\t10.00\\t8000\\t80.00\\t1000\\t10.00\\tNA\\tNA\\t90.00\\tfalse\\n' >> alignment_summary.tsv
    """
}
