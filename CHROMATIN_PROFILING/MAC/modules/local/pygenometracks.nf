/*
 * pyGenomeTracks figure for one locus of one track set (bin/tracks_pygenometracks.sh).
 * The INI (TRACKS_INI) refers to bigWigs, BEDs and the GTF by relative names: they are staged next to it
 * (BEDs under tracks_bed_<n>/, as the INI expects).
 * The automatic peak_qc set gets a wider label column (0.15 of the plot width instead of pyGenomeTracks' default
 * 0.05) so that titles up to ~25 characters are not clipped.
 * No cache 'deep': bigWig and BED inputs.
 */
include { runner } from './functions'

process PYGENOMETRACKS {
    tag "${slug}"
    container 'quay.io/biocontainers/pygenometracks:3.9--pyhdfd78af_0'
    publishDir { "${params.outdir}/tracks/plots" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/tracks/plots/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${slug}_pygenometracks_command.sh" }

    input:
    tuple val(track_set), val(slug), val(region), path(ini), path(bigwigs), path(beds, stageAs: 'tracks_bed_?/*'), path(log2_bigwigs)
    path gtf            // staged so that the INI's [Genes] section finds it by its relative name
    val width           // cm
    val dpi
    path step_script    // bin/tracks_pygenometracks.sh

    output:
    path "${slug}.svg", emit: svg
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --ini "${ini}" \\
        --region "${region}" \\
        --out-prefix "${slug}" \\
        --width ${width} \\
        --dpi ${dpi} \\
        ${track_set == 'peak_qc' ? '--trackLabelFraction 0.15' : ''} \\
        ${params.tracks_pygenometracks_args ?: ''}
    """

    stub:
    """
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "${slug}.svg"
    """
}
