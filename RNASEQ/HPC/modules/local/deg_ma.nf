include { runner; extraArgs } from './functions'
include { shq } from './deg_functions'

process DEG_MA {
    tag "${model}"
    container '/project/rephilllab/nf_containers/omics-r_1.0.sif'   // R 4 + DESeq2 etc., built from the omics Docker image R env lock
    cache 'deep'
    publishDir "${params.outdir}/deg", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/deg/${model}/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "ma_command.sh" }

    input:
    tuple val(model), path(results)   // the model's results_<comparison>.csv files (DEG_RESULTS), staged flat
    path comparisons    // comparisons.tsv (DEG_PREPARE)
    path common         // bin/deg_common.R, or the project's `deg_common` replacement
    val alpha           // params.deg_alpha
    val lfc             // params.deg_lfc
    val color_down      // params.deg_color_down
    val color_up        // params.deg_color_up
    val color_ns        // params.deg_color_ns
    path step_script    // bin/deg_ma.R, or the project's `deg_ma` script replacement

    output:
    path "${model}/*/ma_*.svg", emit: svg
    path "command.sh",          emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --common="${common}" \\
        --model=${shq(model)} \\
        --comparisons="${comparisons}" \\
        --alpha=${alpha} \\
        --lfc=${lfc} \\
        --color-down=${shq(color_down)} \\
        --color-up=${shq(color_up)} \\
        --color-ns=${shq(color_ns)} \\
        --out-dir=. \\
        ${extraArgs('deg_ma')}
    """

    stub:
    def comps = [results].flatten().collect { f -> f.name.replaceFirst(/^results_/, '').replaceFirst(/\.csv$/, '') }
    """
    for c in ${comps.collect { c -> shq(c) }.join(' ')}; do
        mkdir -p "${model}/\$c"
        echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "${model}/\$c/ma_\$c.svg"
    done
    """
}
