/*
 * Density plot of barcode match position per sample, from the small CSV
 * DEMULTIPLEX_SE/PE writes. Plotting script and theme are task inputs
 * (bin/plot_barcode_position.R and bin/my_theme.R, or the project's
 * `plot_barcode_position:` / `my_theme:` replacements in config.yaml).
 *
 * Output capture is `*`: Nextflow excludes staged inputs from output globs,
 * so this picks up whatever the script writes — a custom script that makes
 * extra plots gets those published too, with no module change.
 */
include { runner; extraArgs } from './functions'

process PLOT_BARCODE_POSITION {
    tag "${params.project_name}"
    container 'rocker/tidyverse:4.3'
    // Content-based caching: a re-run upstream step that writes byte-identical
    // output (e.g. mageck count after a script-only change — its count table
    // is deterministic) doesn't force this step to re-run. Nextflow's default
    // cache key uses input path + size + mtime, and a re-run step always
    // writes to a new work dir. Inputs here are small, so hashing is cheap.
    cache 'deep'
    publishDir "${params.outdir}/demultiplexing/qc/demux_plots", mode: params.publish_mode
    // Nextflow doesn't support a negated publishDir pattern to exclude
    // command.sh from the catch-all "*" above, and this process deliberately
    // keeps that catch-all (a custom plot_barcode_position: script's extra
    // outputs need to flow through unchanged). So instead of a root-level
    // command.sh (conf/base.config's default afterScript), place it under
    // log/ here — "*" matches that directory too and publishes it whole,
    // landing at .../demux_plots/log/command.sh with no separate publishDir
    // or output declaration needed.
    afterScript 'mkdir -p log && cp .command.sh log/command.sh 2>/dev/null || true'

    input:
    path counts_csv
    path plot_script
    path theme

    output:
    path "*", emit: plots

    script:
    """
    ${runner(plot_script)} "${plot_script}" \\
        --counts="${counts_csv}" \\
        --project="${params.project_name}" \\
        --theme="${theme}" \\
        --outdir=. \\
        ${extraArgs('plot_barcode_position')}
    """

    stub:
    """
    touch "${params.project_name}_barcode_position_density.svg"
    """
}
