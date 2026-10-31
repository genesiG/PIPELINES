/*
 * deepTools computeMatrix for one heatmap (bin/heatmap_computematrix.sh).
 * regions: staged BED file(s); region_names: comma-separated staged names of those files.
 * bw_names: comma-separated staged bigWig names, same order as sample_ids.
 * COMPUTEMATRIX_CLUSTERED is this same process, included under an alias in main.nf.
 * reference-point, missing-as-zero and average-type-bins come from params.heatmap_*.
 */
include { runner; boolStr } from './functions'

process COMPUTEMATRIX {
    tag "${heatmap}"
    container 'quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0'
    publishDir { "${params.outdir}/heatmaps/matrix" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/heatmaps/matrix/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${out_name}_computematrix_command.sh" }

    input:
    tuple val(heatmap), val(out_name), val(region_names), path(regions), val(sample_ids), val(bw_names), path(bigwigs), val(window), val(bin)
    path step_script   // bin/heatmap_computematrix.sh

    output:
    tuple val(heatmap), path("${out_name}.matrix.gz"), emit: matrix
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --out-prefix "${out_name}" \\
        --regions "${region_names}" \\
        --samples "${sample_ids}" \\
        --bigwigs "${bw_names}" \\
        --window ${window} \\
        --bin-size ${bin} \\
        --threads ${task.cpus} \\
        --reference-point ${params.heatmap_reference_point} \\
        --missing-as-zero ${boolStr(params.heatmap_missing_as_zero)} \\
        --average-type-bins ${params.heatmap_average_type_bins} \\
        ${params.heatmap_computematrix_args ?: ''}
    """

    stub:
    """
    echo stub | gzip > "${out_name}.matrix.gz"
    """
}
