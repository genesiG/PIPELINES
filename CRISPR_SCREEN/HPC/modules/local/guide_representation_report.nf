/*
 * Renders one self-contained HTML QC report per library, assembling what
 * GUIDE_REPRESENTATION_WRANGLE (tables) and GUIDE_REPRESENTATION_PLOT (SVGs)
 * already produced: for each plot, what it measures, how to read it, and
 * what to check if it looks poor — so a wet-lab scientist can triage a run
 * without knowing the underlying statistics by heart.
 *
 * Downstream of both WRANGLE and PLOT (not a third parallel branch): it only
 * reads their outputs, so it never recomputes a table or re-renders a plot.
 *
 * In all-vs-all runs (every sample counted against every library) it also
 * waits for the library breakdown (LIBRARY_BREAKDOWN) and passes its wide
 * table to every per-library report, so a low mapping % can be told apart
 * from reads that belong to another library. Otherwise no breakdown is passed.
 */
include { runner; extraArgs } from './functions'

process GUIDE_REPRESENTATION_REPORT {
    tag "${prefix}"
    container 'rocker/tidyverse:4.3'
    // Content-based caching: a re-run upstream step that writes byte-identical
    // output doesn't force this step to re-run (see guide_representation.nf).
    cache 'deep'
    // command.sh excluded here (via saveAs, not pattern -- publishDir's
    // pattern: does NOT support a closure the way its destination path and
    // saveAs: do; a closure there silently matches nothing, a real bug this
    // session shipped and only caught by re-checking a real run's actual
    // published report file, not just that command.sh appeared), published
    // separately below instead.
    publishDir "${params.outdir}/counts/qc/report", mode: params.publish_mode,
        saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/counts/qc/report/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${prefix}_command.sh" }

    input:
    tuple val(prefix), path(tables), path(plots)   // tables/*.csv (WRANGLE) + plots/*.svg (PLOT), staged flat
    path breakdown_wide   // wide library-breakdown table in all-vs-all runs, or an empty list [] otherwise
    path report_script   // bin/guide_representation_report.R, or the project's `guide_representation_report:` replacement

    output:
    path "${prefix}_qc_report.html", emit: report
    path "command.sh", emit: command_sh

    script:
    // [] stages nothing, so the optional --breakdown argument is only added when a file was staged
    def bd = [breakdown_wide].flatten().findAll { it }
    def breakdown_arg = bd ? "--breakdown=\"${bd[0]}\"" : ''
    """
    ${runner(report_script)} "${report_script}" \\
        --tables-dir=. \\
        --plots-dir=. \\
        --prefix="${prefix}" \\
        --project="${params.project_name}" \\
        --outdir=. ${breakdown_arg} \\
        ${extraArgs('guide_representation_report')}
    """

    stub:
    """
    touch "${prefix}_qc_report.html"
    """
}
