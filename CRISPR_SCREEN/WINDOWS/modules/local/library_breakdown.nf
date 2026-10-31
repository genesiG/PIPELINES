/*
 * How each sample's reads split across the run's libraries
 * (bin/library_breakdown.py): per sample, reads matching each library it was
 * counted against, reads matching any library (a guide shared by several
 * libraries counted once), and reads matching none. Runs when the run has
 * more than one library; most informative when every sample is counted
 * against every library (no samplesheet `library` column).
 *
 * Reads only MAGECK_COUNT's outputs (count tables + <prefix>.sample_reads.tsv)
 * and the MAGeCK-format libraries, collected across every library as
 * parallel lists keyed by library name (main.nf builds them from one joined
 * channel, so entry i of every list belongs to the same library). The tables
 * are published next to the guide-representation tables and are picked up by
 * GUIDE_REPRESENTATION_COMPARE_REPORT.
 */
include { runner } from './functions'

process LIBRARY_BREAKDOWN {
    tag "${params.project_name}"
    // Small inputs downstream of a deterministic tool (see guide_representation.nf).
    cache 'deep'
    publishDir "${params.outdir}/counts/qc/tables", mode: params.publish_mode, pattern: "*_library_breakdown*.csv"
    publishDir "${params.outdir}/counts/qc/tables/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${params.project_name}_library_breakdown_command.sh" }

    input:
    tuple val(lib_names), path(count_tables, stageAs: 'counts/*'), path(sample_reads, stageAs: 'reads/*'), path(mageck_libraries, stageAs: 'libs/*')
    val sample_names   // samplesheet order, comma-joined
    path breakdown_script   // bin/library_breakdown.py, or the project's `library_breakdown:` replacement

    output:
    path "${params.project_name}_library_breakdown.csv",      emit: table
    path "${params.project_name}_library_breakdown_wide.csv", emit: wide
    path "command.sh",                                        emit: command_sh

    script:
    """
    ${runner(breakdown_script)} "${breakdown_script}" \\
        --libraries="${[lib_names].flatten().join(',')}" \\
        --count-tables="${[count_tables].flatten().join(' ')}" \\
        --sample-reads="${[sample_reads].flatten().join(' ')}" \\
        --mageck-libraries="${[mageck_libraries].flatten().join(' ')}" \\
        --samples="${sample_names}" \\
        --project="${params.project_name}" \\
        --outdir=. \\
        ${params.library_breakdown_args ?: ''}
    """

    stub:
    """
    printf 'sample,library,reads,mapped_reads,pct_of_sample_reads,pct_of_matched_reads,status\\n' > "${params.project_name}_library_breakdown.csv"
    printf 'sample,reads\\n' > "${params.project_name}_library_breakdown_wide.csv"
    """
}
