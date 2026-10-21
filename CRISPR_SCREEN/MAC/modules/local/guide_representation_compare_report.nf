/*
 * One HTML report comparing every library in this run side by side: the
 * distributional plots (guide representation, ECDF, Lorenz, density, ...)
 * faceted one panel per pool (3 columns), plus the headline metrics
 * (Gini, 90/10 ratio, % undetected, mapping rate) as combined bar charts.
 * Only runs when the samplesheet's `library` column produced more than one
 * library — with a single library, this would just repeat that library's
 * own GUIDE_REPRESENTATION_REPORT.
 *
 * Reads only GUIDE_REPRESENTATION_WRANGLE's tables (collected across every
 * library) and rebuilds its own comparison plots from them — it never reruns
 * `mageck count` or reads another library's count table directly.
 *
 * Each comparison plot is also written as SVG + PNG under pools_comparison/,
 * published to counts/qc/plots/pools_comparison.
 */
include { runner } from './functions'

process GUIDE_REPRESENTATION_COMPARE_REPORT {
    tag "${params.project_name}"
    // Content-based caching: a re-run upstream step that writes byte-identical
    // output doesn't force this step to re-run (see guide_representation.nf).
    cache 'deep'
    // command.sh excluded here, published separately below -- see
    // guide_representation.nf for why (no pattern would otherwise also match it).
    publishDir "${params.outdir}/counts/qc/report", mode: params.publish_mode, pattern: "${params.project_name}_pools_comparison_qc_report.html"
    publishDir "${params.outdir}/counts/qc/plots/pools_comparison", mode: params.publish_mode, pattern: "pools_comparison/*", saveAs: { fn -> file(fn).name }
    publishDir "${params.outdir}/counts/qc/report/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${params.project_name}_pools_comparison_command.sh" }

    input:
    path tables   // every library's tables/*.csv from GUIDE_REPRESENTATION_WRANGLE, collected
    path report_script   // bin/guide_representation_compare_report.R, or the project's `guide_representation_compare_report:` replacement
    path theme

    output:
    path "${params.project_name}_pools_comparison_qc_report.html", emit: report
    path "pools_comparison/*", emit: plots
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(report_script)} "${report_script}" \\
        --tables-dir=. \\
        --project="${params.project_name}" \\
        --theme="${theme}" \\
        --outdir=. \\
        ${params.guide_representation_compare_report_args ?: ''}
    """

    stub:
    """
    touch "${params.project_name}_pools_comparison_qc_report.html"
    mkdir -p pools_comparison
    touch "pools_comparison/${params.project_name}_pools_stub.svg" "pools_comparison/${params.project_name}_pools_stub.png"
    """
}
