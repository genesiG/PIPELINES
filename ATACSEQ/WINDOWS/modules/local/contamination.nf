include { runner; extraArgs } from './functions'

process CONTAMINATION {
    tag "${sample_id}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    publishDir "${params.outdir}/qc/contamination", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/contamination/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_contamination_command.sh" }

    input:
    tuple val(sample_id), path(bam)   // S.aligned.bam (ALIGN output, unsorted)
    path index_files                  // all *.bt2 / *.bt2l files of the contaminant index, staged flat
    val index_base                    // bare file-name prefix of the contaminant index
    path step_script                  // bin/contamination.sh, or the project's `contamination:` replacement

    output:
    path "${sample_id}.contamination.tsv",         emit: table
    path "${sample_id}.contamination_summary.txt", emit: summary
    path "command.sh",                             emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --index-base "${index_base}" \\
        --threads ${task.cpus} \\
        ${extraArgs('contamination')}
    """

    stub:
    """
    printf 'organism\\taligned_reads\\tpct_of_unmapped\\ncontaminant_a\\t100\\t1.0000\\n' > "${sample_id}.contamination.tsv"
    echo "stub contamination summary" > "${sample_id}.contamination_summary.txt"
    """
}
