include { runner; extraArgs } from './functions'

process QC_STRANDEDNESS_TABLE {
    tag "strandedness"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    cache 'deep'
    publishDir "${params.outdir}/qc/strandedness", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/strandedness/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "strandedness_table_command.sh" }

    input:
    path '*'           // all S.infer_experiment.txt files (INFER_STRANDEDNESS), collected and staged flat
    val strandedness   // configured strandedness (params.strandedness)
    val threshold      // params.strand_call_threshold
    path step_script   // bin/qc_strandedness_table.py, or the project's `qc_strandedness_table` script replacement

    output:
    path "strandedness.tsv", emit: table
    path "command.sh",       emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --dir . \\
        --strandedness ${strandedness} \\
        --threshold ${threshold} \\
        --out strandedness.tsv \\
        ${extraArgs('qc_strandedness_table')}
    """

    stub:
    """
    printf 'sample_id\\tfraction_forward\\tfraction_reverse\\tfraction_undetermined\\tcall\\tconfigured\\tmismatch\\n' > strandedness.tsv
    printf 'stub1\\t0.5000\\t0.4900\\t0.0100\\tunstranded\\tunstranded\\tfalse\\n' >> strandedness.tsv
    """
}
