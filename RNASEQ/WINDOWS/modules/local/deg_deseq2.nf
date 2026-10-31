include { runner; extraArgs } from './functions'
include { shq } from './deg_functions'

process DEG_DESEQ2 {
    tag "${model}"
    container '/project/rephilllab/nf_containers/omics-r_1.0.sif'   // R 4 + DESeq2 etc., built from the omics Docker image R env lock
    cache 'deep'
    publishDir "${params.outdir}/deg", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/deg/${model}/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "deseq2_command.sh" }

    input:
    tuple val(model), val(design)   // one model of deg_comparisons: its folder name and its design formula
    path counts         // counts_edited.rds (DEG_PREPARE)
    path metadata       // sample_metadata_edited.rds (DEG_PREPARE)
    path common         // bin/deg_common.R, or the project's `deg_common` replacement
    val min_count       // params.deg_min_count
    val min_samples     // params.deg_min_samples
    val fit_type        // params.deg_fit_type
    path step_script    // bin/deg_deseq2.R, or the project's `deg_deseq2` script replacement

    output:
    tuple val(model), path("${model}/DESeq.rds"),        emit: deseq     // fitted model (DEG_RESULTS)
    tuple val(model), path("${model}/DESeqDataSet.rds"), emit: dataset   // data set with the design (DEG_TRANSFORM)
    path "${model}/design.txt",           emit: design
    path "${model}/dispersion_plots.svg", emit: dispersion
    path "${model}/prefilter.tsv",        emit: prefilter
    path "command.sh",                    emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --common="${common}" \\
        --counts="${counts}" \\
        --metadata="${metadata}" \\
        --model=${shq(model)} \\
        --design=${shq(design)} \\
        --min-count=${min_count} \\
        --min-samples=${min_samples} \\
        --fit-type=${shq(fit_type)} \\
        --out-dir=. \\
        ${extraArgs('deg_deseq2')}
    """

    stub:
    """
    mkdir -p "${model}"
    echo stub > "${model}/DESeq.rds"
    echo stub > "${model}/DESeqDataSet.rds"
    echo ${shq(design)} > "${model}/design.txt"
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "${model}/dispersion_plots.svg"
    printf 'gene\\tkept\\n' > "${model}/prefilter.tsv"
    """
}
