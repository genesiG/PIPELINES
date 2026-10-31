/*
 * k-means clustering of heatmap rows (bin/heatmap_cluster.py): one k, or a range with gap statistic.
 * matrix: deepTools computeMatrix output (or a plain numeric TSV); regions: the original region BED (same row order).
 * cluster_samples: comma-separated sample ids used for clustering (empty = all); k_spec: "3" or "2-6".
 */
include { runner } from './functions'

process HEATMAP_CLUSTER {
    tag "${heatmap}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    publishDir { "${params.outdir}/heatmaps/cluster" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/heatmaps/cluster/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${heatmap}_heatmap_cluster_command.sh" }

    input:
    tuple val(heatmap), path(matrix), path(regions), val(cluster_samples), val(k_spec)
    val seed
    path step_script   // bin/heatmap_cluster.py

    output:
    tuple val(heatmap), path("${heatmap}.cluster_assignments.csv"), emit: assignments
    path "${heatmap}.cluster_k.tsv", emit: k_table
    path "${heatmap}.clusters.tsv", emit: clusters
    tuple val(heatmap), path("${heatmap}.cluster*.bed"), emit: beds
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --matrix "${matrix}" \\
        --regions "${regions}" \\
        --samples "${cluster_samples}" \\
        --k "${k_spec}" \\
        --seed ${seed} \\
        --out-prefix "${heatmap}" \\
        --threads ${task.cpus} \\
        ${params.heatmap_cluster_args ?: ''}
    """

    stub:
    """
    kk="${k_spec}"
    kk="\${kk%%-*}"
    for i in 0 1 2 3 4 5; do echo \$(( i % kk + 1 )); done > "${heatmap}.cluster_assignments.csv"
    printf 'heatmap\\tk\\tsilhouette\\tgap\\tgap_se\\tchosen\\n' > "${heatmap}.cluster_k.tsv"
    printf '%s\\t%s\\tNA\\tNA\\tNA\\tTRUE\\n' "${heatmap}" "\$kk" >> "${heatmap}.cluster_k.tsv"
    printf 'heatmap\\tk\\tcluster\\tn_regions\\tmean_signal\\n' > "${heatmap}.clusters.tsv"
    for c in \$(seq 1 "\$kk"); do
        printf '%s\\t%s\\t%s\\t%s\\t%s\\n' "${heatmap}" "\$kk" "\$c" 1 0 >> "${heatmap}.clusters.tsv"
        printf 'chr1\\t0\\t100\\n' > "${heatmap}.cluster\${c}.bed"
    done
    """
}
