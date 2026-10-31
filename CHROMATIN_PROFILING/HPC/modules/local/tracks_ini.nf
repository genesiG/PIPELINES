/*
 * pyGenomeTracks configuration for one track set (bin/make_tracks_ini.py). Writes NAME.tracks.ini (relative file names,
 * staged by PYGENOMETRACKS) and NAME.regions.tsv (slug<TAB>region, one line per locus).
 * loci / flanks / samples / bw_names / labels / colors / ... are comma-separated lists in the same order.
 * A locus may carry a LABEL= prefix (the label becomes its slug in regions.tsv; the stub honours it too).
 * max_values: one entry per sample (number | auto | group:NAME), empty = pyGenomeTracks autoscale.
 * bed_colors / bed_labels: one entry per bed file (empty = make_tracks_ini.py defaults).
 * Styling options (number of bins, track heights, font sizes, log2 limits) come from params.track_*.
 * Per-set style overrides (e.g. peak_qc: scale bar, gene label size and height) arrive as extra_args, a string of
 * options placed after the params.track_* options (later argparse options win).
 */
include { runner } from './functions'

process TRACKS_INI {
    tag "${track_set}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    publishDir { "${params.outdir}/tracks/ini" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/tracks/ini/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${track_set}_tracks_ini_command.sh" }

    input:
    tuple val(track_set), val(loci), val(flanks), val(sample_ids), val(bw_names), val(labels), val(colors), val(max_values), val(bed_names), val(bed_colors), val(bed_labels), val(scalebar), val(log2_names), val(log2_for), val(extra_args)
    path gtf
    path step_script   // bin/make_tracks_ini.py

    output:
    tuple val(track_set), path("${track_set}.tracks.ini"), path("${track_set}.regions.tsv"), emit: ini
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --track-set "${track_set}" \\
        --loci "${loci}" \\
        --flanks "${flanks}" \\
        --gtf "${gtf}" \\
        --gtf-name "${gtf.name}" \\
        --samples "${sample_ids}" \\
        --bigwigs "${bw_names}" \\
        --labels "${labels}" \\
        --colors "${colors}" \\
        ${max_values ? "--max-values '${max_values}'" : ''} \\
        ${bed_names ? "--bed-files '${bed_names}'" : ''} \\
        ${bed_colors ? "--bed-colors='${bed_colors}'" : ''} \\
        ${bed_labels ? "--bed-labels='${bed_labels}'" : ''} \\
        ${scalebar ? "--scalebar '${scalebar}'" : ''} \\
        ${log2_names ? "--log2-bigwigs '${log2_names}'" : ''} \\
        ${log2_for ? "--log2-for '${log2_for}'" : ''} \\
        --number-of-bins=${params.track_number_of_bins} \\
        --height=${params.track_height} \\
        --fontsize=${params.track_fontsize} \\
        --bed-height=${params.track_bed_height} \\
        --scalebar-fontsize=${params.track_scalebar_fontsize} \\
        --genes-height=${params.track_genes_height} \\
        --log2-min=${params.track_log2_min} \\
        --log2-max=${params.track_log2_max} \\
        ${extra_args ?: ''} \\
        ${params.make_tracks_ini_args ?: ''}
    """

    stub:
    """
    printf '[Genes]\\nfile = x.gtf\\ntype = genes\\n' > "${track_set}.tracks.ini"
    LOCI=\$(echo "${loci}" | sed -E 's/([0-9]),([0-9]{3})/\\1\\2/g; s/([0-9]),([0-9]{3})/\\1\\2/g')
    IFS=',' read -ra L <<< "\$LOCI"
    IFS=',' read -ra F <<< "${flanks}"
    : > "${track_set}.regions.tsv"
    for i in "\${!L[@]}"; do
        locus="\${L[\$i]}"
        if [ "\${#F[@]}" -eq 1 ]; then f="\${F[0]}"; else f="\${F[\$i]}"; fi
        label=""
        if [[ "\$locus" == *=* ]]; then label="\${locus%%=*}"; locus="\${locus#*=}"; fi
        if [[ "\$locus" == *:* ]]; then
            slug="\$(printf '%s' "\$locus" | tr -d ',' | sed 's/[:-]/_/g')"
        elif [ \$(( f % 1000 )) -eq 0 ]; then
            slug="\${locus}_\$(( f / 1000 ))kb"
        else
            slug="\${locus}_\${f}bp"
        fi
        if [ -n "\$label" ]; then slug="\$label"; fi
        printf '%s\\t%s\\n' "${track_set}.\$slug" "chr1:1-1000000" >> "${track_set}.regions.tsv"
    done
    """
}
