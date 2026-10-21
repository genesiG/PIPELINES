/*
 * Discovers the constant vector sequence between each library's barcode and
 * its guide (bin/discover_anchor.py — greedy consensus extension, no prior
 * knowledge of the anchor's sequence or length needed), from one
 * representative sample per library. GUIDE_ANCHOR_TRIM then uses it to
 * position every sample's guide correctly regardless of any variable-length
 * spacing in the amplicon, instead of assuming a fixed trim5 offset.
 *
 * This library's MAGeCK-format guides are passed too (--library), so the
 * consensus anchor is refined against them: it must end where this
 * library's guides begin, not inside them (bin/discover_anchor.py's header
 * explains the G + 19-nt case this fixes). The refinement table is
 * published next to the discovery trace.
 *
 * Only runs when params.trim5 is unset (main.nf) — a project can still set
 * trim5 explicitly to bypass discovery entirely and use the old
 * fixed-position behavior.
 */
include { runner } from './functions'

process GUIDE_ANCHOR_DISCOVERY {
    tag "${lib_name}"
    container 'quay.io/biocontainers/cutadapt:5.2--py39hbcbf7aa_0'   // python3 stdlib only; reused to avoid an extra pull
    // command.sh excluded here, published separately below -- see
    // guide_representation.nf for why (no pattern would otherwise also match it).
    publishDir "${params.outdir}/counts/qc/anchor_discovery", mode: params.publish_mode, pattern: "*_anchor_discovery.csv"
    publishDir "${params.outdir}/counts/qc/anchor_discovery", mode: params.publish_mode, pattern: "*_anchor_refinement.csv"
    publishDir "${params.outdir}/counts/qc/anchor_discovery/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${lib_name}_discover_command.sh" }

    input:
    tuple val(lib_name), path(sample_fastq), path(mageck_library)   // one representative (first) sample of this library + its guides
    path discover_script   // bin/discover_anchor.py, or the project's `discover_anchor:` replacement

    output:
    tuple val(lib_name), path("${params.project_name}_${lib_name}_anchor.txt"), emit: anchor
    path "${params.project_name}_${lib_name}_anchor_discovery.csv",              emit: qc
    path "${params.project_name}_${lib_name}_anchor_refinement.csv",             emit: refinement
    path "command.sh",                                                          emit: command_sh

    script:
    """
    ${runner(discover_script)} "${discover_script}" \\
        --fastq "${sample_fastq}" \\
        --project-name "${params.project_name}" \\
        --library-name "${lib_name}" \\
        --library "${mageck_library}" \\
        --outdir . \\
        ${params.discover_anchor_args ?: ''}
    """

    stub:
    """
    echo "NNNNNNNN" > "${params.project_name}_${lib_name}_anchor.txt"
    touch "${params.project_name}_${lib_name}_anchor_discovery.csv"
    touch "${params.project_name}_${lib_name}_anchor_refinement.csv"
    """
}
