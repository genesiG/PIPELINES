/*
 * Resolves group:NAME max values of one locus's pyGenomeTracks ini to numbers (bin/tracks_limits.py).
 * Writes SLUG.resolved.ini (the ini PYGENOMETRACKS renders for this locus) and SLUG_limits.tsv (one row per grouped track).
 * ini: the TRACKS_INI ini of the track set; bigwigs: the bigWigs it names, staged beside it.
 */
include { runner } from './functions'

process TRACKS_LIMITS {
    tag "${slug}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    publishDir { "${params.outdir}/tracks/limits" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/tracks/limits/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${slug}_tracks_limits_command.sh" }

    input:
    tuple val(track_set), val(slug), val(region), path(ini), path(bigwigs)
    val limit_stat     // mean_bins | max_bin (params.tracks_limit_stat)
    path step_script   // bin/tracks_limits.py

    output:
    tuple val(track_set), val(slug), val(region), path("${slug}.resolved.ini"), emit: ini
    path "${slug}_limits.tsv", emit: limits
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --ini "${ini}" \\
        --region "${region}" \\
        --out-prefix "${slug}" \\
        --limit-stat ${limit_stat} \\
        ${params.tracks_limits_args ?: ''}
    """

    stub:
    """
    cp "${ini}" "${slug}.resolved.ini"
    printf 'group\\ttrack\\tfile\\tstat\\ttrack_peak\\tlimit\\n' > "${slug}_limits.tsv"
    """
}
