include { runner; extraArgs } from './functions'
include { shq; rbool } from './deg_functions'

process DEG_RESULTS {
    tag "${model}"
    container '/project/rephilllab/nf_containers/omics-r_1.0.sif'   // R 4 + DESeq2 etc., built from the omics Docker image R env lock
    cache 'deep'
    publishDir "${params.outdir}/deg", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/deg/${model}/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "results_command.sh" }

    input:
    tuple val(model), path(deseq)   // DESeq.rds of the model (DEG_DESEQ2)
    path comparisons    // comparisons.tsv (DEG_PREPARE)
    path common         // bin/deg_common.R, or the project's `deg_common` replacement
    val alpha           // params.deg_alpha
    val lfc             // params.deg_lfc
    val lfc_threshold_in_results  // params.deg_lfc_threshold_in_results (true | false)
    val shrink_type     // params.deg_shrink_type: ashr | apeglm | none
    val sval_threshold  // params.deg_sval_threshold
    val legacy_summary_counts  // params.deg_legacy_summary_counts (true | false)
    path step_script    // bin/deg_results.R, or the project's `deg_results` script replacement

    output:
    tuple val(model), path("${model}/*/results_*.csv"), emit: results   // one per comparison of the model, written to <model>/<comparison>/
    path "${model}/analysis_summary.csv",  emit: summary_csv
    path "${model}/analysis_summary.txt",  emit: summary_txt
    path "${model}/normalized_counts.csv", emit: normalized
    path "command.sh",                     emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --common="${common}" \\
        --deseq="${deseq}" \\
        --comparisons="${comparisons}" \\
        --model=${shq(model)} \\
        --alpha=${alpha} \\
        --lfc=${lfc} \\
        --lfc-threshold-in-results=${rbool(lfc_threshold_in_results)} \\
        --shrink-type=${shq(shrink_type)} \\
        --sval-threshold=${sval_threshold} \\
        --legacy-summary-counts=${rbool(legacy_summary_counts)} \\
        --out-dir=. \\
        ${extraArgs('deg_results')}
    """

    stub:
    // one results table per comparison of this model; the comparison names are read from the (stub) comparisons.tsv by column name
    """
    mkdir -p "${model}"
    awk -F'\\t' -v m=${shq(model)} 'NR == 1 { for (i = 1; i <= NF; i++) h[\$i] = i; next } \$(h["model"]) == m { print \$(h["comparison"]) }' "${comparisons}" |
    while read -r c; do
        mkdir -p "${model}/\$c"
        printf 'Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj\\nstub_gene,10,1,0.1,5,0.01,0.02\\n' > "${model}/\$c/results_\$c.csv"
    done
    printf 'comparison,up,down\\n' > "${model}/analysis_summary.csv"
    echo "stub summary" > "${model}/analysis_summary.txt"
    printf 'Gene,s1\\nstub_gene,1\\n' > "${model}/normalized_counts.csv"
    """
}
