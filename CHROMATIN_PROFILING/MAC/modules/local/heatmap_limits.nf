/*
 * Shared colour limits for one heatmap (bin/heatmap_limits.py), run on the plotted computeMatrix output.
 * sample_ids / antibodies: comma-separated, same order as the sample columns of the matrix.
 * Writes NAME_limits.tsv (sample, antibody, profile_max, limit); main.nf substitutes the limit into z_max / y_max.
 */
include { runner } from './functions'

process HEATMAP_LIMITS {
    tag "${heatmap}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    publishDir { "${params.outdir}/heatmaps/limits" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/heatmaps/limits/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${heatmap}_heatmap_limits_command.sh" }

    input:
    tuple val(heatmap), path(matrix), val(sample_ids), val(antibodies)
    path step_script   // bin/heatmap_limits.py

    output:
    tuple val(heatmap), path("${heatmap}_limits.tsv"), emit: limits
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --matrix "${matrix}" \\
        --out-prefix "${heatmap}" \\
        --samples "${sample_ids}" \\
        --antibodies "${antibodies}" \\
        ${params.heatmap_limits_args ?: ''}
    """

    stub:
    """
    printf 'sample\\tantibody\\tprofile_max\\tlimit\\n' > "${heatmap}_limits.tsv"
    IFS=',' read -ra S <<< "${sample_ids}"
    IFS=',' read -ra A <<< "${antibodies}"
    for i in "\${!S[@]}"; do
        printf '%s\\t%s\\t5\\t10\\n' "\${S[\$i]}" "\${A[\$i]}" >> "${heatmap}_limits.tsv"
    done
    """
}
