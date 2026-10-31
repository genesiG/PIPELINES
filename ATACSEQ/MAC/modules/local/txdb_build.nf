/*
 * TxDb (gene and transcript models) from the annotation GTF, saved as SQLite (bin/txdb_build.R).
 * Shared by DIFF_ACCESS (promoters) and DIFF_ANNOTATE (peak annotation), so both use one annotation.
 */
include { runner; extraArgs } from './functions'

process TXDB_BUILD {
    container '/project/rephilllab/nf_containers/omics-r_1.0.sif'
    publishDir "${params.outdir}/annotation/txdb", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/annotation/txdb/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "txdb_build_command.sh" }

    input:
    path gtf                // annotation GTF (T2T/hs1 NCBI RefSeq curated GTF, from the project config)
    val species             // species: organism of the TxDb ('Homo sapiens')
    path step_script        // bin/txdb_build.R, or the project's `txdb_build:` replacement

    output:
    path "txdb.sqlite",     emit: txdb
    path "command.sh",      emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --gtf=${gtf} \\
        --out=txdb.sqlite \\
        --species="${species}" \\
        ${extraArgs('txdb_build')}
    """

    stub:
    """
    touch txdb.sqlite
    """
}
