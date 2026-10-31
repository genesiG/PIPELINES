/*
 * Fails fast (seconds, not after mageck count) if contrasts.csv can't be
 * resolved: a pattern matching no sample, or treatment and control that only
 * exist in different libraries' count tables. Each contrast is resolved to
 * the one count table holding both sides; MAGECK_TEST then runs against that
 * table. Exits 3 on any problem (a data problem — not retried).
 */
include { runner } from './functions'

process VALIDATE_CONTRASTS {
    tag "${params.project_name}"
    // command.sh excluded here, published separately below -- see
    // guide_representation.nf for why (no pattern would otherwise also match it).
    publishDir "${params.outdir}/mageck_test", mode: params.publish_mode, pattern: "validated_contrasts.csv"
    publishDir "${params.outdir}/mageck_test/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${params.project_name}_validate_contrasts_command.sh" }

    input:
    path count_tables    // every library's count table
    path contrasts_csv
    path helper_script   // bin/validate_contrasts.py, or the project's `validate_contrasts:` replacement

    output:
    path "validated_contrasts.csv", emit: contrasts
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(helper_script)} "${helper_script}" \\
        --count-table ${count_tables} \\
        --contrasts "${contrasts_csv}" \\
        --output validated_contrasts.csv \\
        ${params.validate_contrasts_args ?: ''}
    """

    stub:   // runs the real (stdlib-only) validation, so -stub-run exercises it too
    """
    ${runner(helper_script)} "${helper_script}" \\
        --count-table ${count_tables} \\
        --contrasts "${contrasts_csv}" \\
        --output validated_contrasts.csv \\
        ${params.validate_contrasts_args ?: ''}
    """
}
