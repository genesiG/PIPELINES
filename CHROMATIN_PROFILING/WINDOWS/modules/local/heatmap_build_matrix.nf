/*
 * Per-base-pair normalised signal matrices for one heatmap, in R (bin/heatmap_build_matrix.R).
 * sample_ids / bw_names / bigwigs: same order. cluster_samples: subset used for clustering (empty = none written).
 */
include { runner } from './functions'

process BUILD_MATRIX {
    tag "${heatmap}"
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    publishDir { "${params.outdir}/heatmaps/matrix" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/heatmaps/matrix/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${heatmap}_build_matrix_command.sh" }

    input:
    tuple val(heatmap), path(regions), val(sample_ids), val(bw_names), path(bigwigs), val(window), val(bin), val(cluster_samples)
    path step_script   // bin/heatmap_build_matrix.R

    output:
    tuple val(heatmap), path("${heatmap}.norm_mats.rds"), path("${heatmap}.targets.rds"), emit: mats
    tuple val(heatmap), path("${heatmap}.cluster_matrix.tsv.gz"), optional: true, emit: cluster_matrix
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --regions-bed=${regions} \\
        --bigwigs=${bw_names} \\
        --names=${sample_ids} \\
        --window=${window} \\
        --bin-size=${bin} \\
        --out-prefix=${heatmap} \\
        --cluster-samples=${cluster_samples} \\
        --cpus=${task.cpus} \\
        ${params.heatmap_build_matrix_args ?: ''}
    """

    stub:
    """
    touch "${heatmap}.norm_mats.rds" "${heatmap}.targets.rds"
    if [ -n "${cluster_samples}" ]; then
        printf '1\\t2\\n3\\t4\\n' | gzip > "${heatmap}.cluster_matrix.tsv.gz"
    fi
    """
}
