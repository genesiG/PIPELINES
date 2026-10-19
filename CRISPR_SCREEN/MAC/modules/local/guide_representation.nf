/*
 * Guide-representation QC, split into two steps so plots can be restyled
 * and re-run without recomputing anything:
 *
 *   GUIDE_REPRESENTATION_WRANGLE — bin/guide_representation_wrangle.R:
 *     count table (+ MAGeCK library, + mageck count log, + unmapped-reads
 *     file) -> CSV tables (Lorenz/Gini, ECDF input, percentile ratio,
 *     undetected guides, mapping rate, GC vs depth, replicate correlations,
 *     dropout, near-miss/far-miss read classification, QC summary).
 *   GUIDE_REPRESENTATION_PLOT — bin/guide_representation_plot.R:
 *     those tables -> PNGs. Reads nothing else.
 *
 * Both scripts (and the plot theme) are task inputs, replaceable per project
 * from config.yaml (key = script name, e.g. `guide_representation_plot:`). Because the plot step's only
 * inputs are the tables, its script and theme, changing any of those three
 * re-runs just the plot step under -resume.
 */
include { runner; extraArgs } from './functions'

process GUIDE_REPRESENTATION_WRANGLE {
    tag "${count_table.simpleName}"
    // Content-based caching: a re-run upstream step that writes byte-identical
    // output (e.g. mageck count after a script-only change — its count table
    // is deterministic) doesn't force this step to re-run. Nextflow's default
    // cache key uses input path + size + mtime, and a re-run step always
    // writes to a new work dir. Inputs here are small, so hashing is cheap.
    cache 'deep'
    // command.sh (root-level, from afterScript) is excluded here (pattern
    // scopes to the tables/ subdir) and published separately below.
    publishDir "${params.outdir}/counts/qc/tables", mode: params.publish_mode, pattern: "tables/*", saveAs: { fn -> file(fn).name }
    publishDir "${params.outdir}/counts/qc/tables/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${lib_name}_command.sh" }

    input:
    // one tuple per library, so a pool's count table is always paired with
    // its own log, unmapped-reads file, and library (separate channels
    // would pair by arrival order)
    tuple val(lib_name), path(count_table), path(mageck_log), path(unmapped_reads), path(mageck_library)
    path wrangle_script

    output:
    path "tables/*", emit: tables   // <prefix>_*.csv; main.nf recovers <prefix> from the file names
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(wrangle_script)} "${wrangle_script}" \\
        --count-table="${count_table}" \\
        --library="${mageck_library}" \\
        --mageck-log="${mageck_log}" \\
        --unmapped="${unmapped_reads}" \\
        --outdir=tables \\
        ${extraArgs('guide_representation_wrangle')}
    """

    stub:
    def stub_prefix = count_table.name - '.count.txt'
    """
    mkdir -p tables
    touch "tables/${stub_prefix}_counts_long.csv" "tables/${stub_prefix}_qc_summary.csv"
    """
}

process GUIDE_REPRESENTATION_PLOT {
    tag "${prefix}"
    // Content-based caching: a re-run upstream step that writes byte-identical
    // output (e.g. mageck count after a script-only change — its count table
    // is deterministic) doesn't force this step to re-run. Nextflow's default
    // cache key uses input path + size + mtime, and a re-run step always
    // writes to a new work dir. Inputs here are small, so hashing is cheap.
    cache 'deep'
    // command.sh excluded here, published separately below -- see WRANGLE above.
    publishDir "${params.outdir}/counts/qc/plots", mode: params.publish_mode, pattern: "plots/*", saveAs: { fn -> file(fn).name }
    publishDir "${params.outdir}/counts/qc/plots/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${prefix}_command.sh" }

    input:
    tuple val(prefix), path(tables)   // staged flat into the task directory
    path plot_script
    path theme

    output:
    path "plots/*", emit: plots
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(plot_script)} "${plot_script}" \\
        --tables-dir=. \\
        --prefix="${prefix}" \\
        --theme="${theme}" \\
        --outdir=plots \\
        ${extraArgs('guide_representation_plot')}
    """

    stub:
    // Real name, not a generic "_stub" placeholder: main.nf finds this
    // library's prefix from the one plot every run always produces with no
    // suffix collisions (<prefix>_guide_representation.svg), including
    // under -stub-run.
    """
    mkdir -p plots
    touch "plots/${prefix}_guide_representation.svg"
    """
}
