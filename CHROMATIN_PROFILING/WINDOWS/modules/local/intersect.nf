/*
 * Overlap of two peak sets (bin/intersect_beds.py): A records that overlap B (output records) or the overlap pieces (output overlap),
 * A records that do not overlap, and a summary row.
 * bed_a / bed_b are staged under a/ and b/ so that the two files never share a name.
 */
include { runner } from './functions'

process INTERSECT {
    tag "${name}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    publishDir { "${params.outdir}/intersect" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/intersect/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${name}_intersect_command.sh" }

    input:
    tuple val(name), val(label_a), val(label_b), val(min_overlap), val(reciprocal), val(out_mode), path(bed_a, stageAs: 'a/*'), path(bed_b, stageAs: 'b/*')
    path step_script   // bin/intersect_beds.py

    output:
    path "${name}.intersect.bed", emit: bed
    path "${name}.only_a.bed", emit: only_a
    path "${name}.intersect_summary.tsv", emit: summary
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --name "${name}" \\
        --bed-a "${bed_a}" \\
        --bed-b "${bed_b}" \\
        --label-a "${label_a}" \\
        --label-b "${label_b}" \\
        ${min_overlap ? "--min-overlap '${min_overlap}'" : ''} \\
        ${reciprocal ? "--reciprocal '${reciprocal}'" : ''} \\
        --output "${out_mode}" \\
        ${params.intersect_beds_args ?: ''}
    """

    stub:
    """
    touch "${name}.intersect.bed" "${name}.only_a.bed"
    printf 'name\\tlabel_a\\tlabel_b\\tn_a\\tn_b\\tn_a_overlapping_b\\tn_a_only\\tpct_a_overlapping_b\\tn_b_overlapping_a\\tmin_overlap\\treciprocal\\toutput\\tn_pieces\\n' > "${name}.intersect_summary.tsv"
    printf '%s\\t%s\\t%s\\t0\\t0\\t0\\t0\\t0\\t0\\tstub\\tstub\\t%s\\t0\\n' "${name}" "${label_a}" "${label_b}" "${out_mode}" >> "${name}.intersect_summary.tsv"
    """
}
