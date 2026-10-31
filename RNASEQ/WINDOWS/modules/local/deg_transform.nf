include { runner; extraArgs } from './functions'
include { shq; rbool } from './deg_functions'

process DEG_TRANSFORM {
    tag "${model}"
    container '/project/rephilllab/nf_containers/omics-r_1.0.sif'   // R 4 + DESeq2 etc., built from the omics Docker image R env lock
    cache 'deep'
    publishDir "${params.outdir}/deg", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/deg/${model}/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "transform_command.sh" }

    input:
    tuple val(model), path(dataset)   // DESeqDataSet.rds of the model (DEG_DESEQ2)
    path common         // bin/deg_common.R, or the project's `deg_common` replacement
    val method          // params.deg_transform_method: rlog | vst
    val blind           // params.deg_transform_blind (true | false)
    path step_script    // bin/deg_transform.R, or the project's `deg_transform` script replacement

    output:
    tuple val(model), path("${model}/${method}Counts.rds"), emit: counts   // data.frame with a Gene column (heatmaps, gene plots)
    path "${model}/${method}Data.rds", emit: data
    path "command.sh",                 emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --common="${common}" \\
        --dataset="${dataset}" \\
        --model=${shq(model)} \\
        --method=${shq(method)} \\
        --blind=${rbool(blind)} \\
        --fit-type=${shq(params.deg_fit_type)} \\
        --out-dir=. \\
        ${extraArgs('deg_transform')}
    """

    stub:
    """
    mkdir -p "${model}"
    echo stub > "${model}/${method}Data.rds"
    echo stub > "${model}/${method}Counts.rds"
    """
}
