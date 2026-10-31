include { runner; extraArgs } from './functions'

process COUNTS_MERGE {
    tag "${project_name}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    cache 'deep'
    publishDir "${params.outdir}/counts", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/counts/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "counts_merge_command.sh" }

    input:
    path '*'           // every per-sample S_counts.txt, S_featurecounts.txt and S_featurecounts.txt.summary, collected and staged flat
    val sample_ids     // comma list of the non-excluded sample ids, samplesheet order
    val project_name   // params.project_name: names the counts matrix
    path step_script   // bin/counts_merge.py, or the project's `counts_merge` script replacement

    output:
    path "${project_name}_counts.txt", emit: counts
    path "gene_length.txt",            emit: gene_length
    path "featurecounts_summary.tsv",  emit: summary
    path "command.sh",                 emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --dir . \\
        --samples "${sample_ids}" \\
        --project-name "${project_name}" \\
        --out-counts "${project_name}_counts.txt" \\
        ${extraArgs('counts_merge')}
    """

    stub:
    """
    printf 'Gene\\tstub1\\nstub_gene\\t5\\n' > "${project_name}_counts.txt"
    printf 'Geneid\\tStrand\\tLength\\nstub_gene\\t+\\t100\\n' > gene_length.txt
    printf 'Status\\tstub1\\nAssigned\\t5\\n' > featurecounts_summary.tsv
    """
}
