/*
 * CPM "scale factors" of one norm_run (bin/norm_cpm.py): library size only, every factor 1. Standard library only.
 * sample_ids: comma-separated, run order. Output <norm_run>.cpm.scale_factors.tsv (7-column schema).
 */
include { runner } from './functions'

process NORM_CPM {
    tag "${norm_run}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    publishDir { "${params.outdir}/scale_factors" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn.endsWith('.scale_factors.tsv') ? "methods/${fn}" : fn }
    publishDir { "${params.outdir}/scale_factors/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${norm_run}_cpm_command.sh" }

    input:
    tuple val(norm_run), val(sample_ids)
    path step_script   // bin/norm_cpm.py

    output:
    tuple val(norm_run), val('cpm'), path("${norm_run}.cpm.scale_factors.tsv"), emit: sf
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \
        --ids ${sample_ids} \
        --run ${norm_run} \
        --out ${norm_run}.cpm.scale_factors.tsv \
        ${params.norm_cpm_args ?: ''}
    """

    stub:
    """
    printf 'sample_id\tmethod\traw_factor\tbigwig_multiplier\tedger_norm_factor\tqc_pass\tnote\n' > "${norm_run}.cpm.scale_factors.tsv"
    for id in \$(echo "${sample_ids}" | tr ',' ' '); do
        printf '%s\tcpm\t1\t1\t1\tTRUE\tlibrary size only; bigWig denominator = final BAM mapped reads (deepTools); edgeR lib.size = csaw total\n' "\$id" >> "${norm_run}.cpm.scale_factors.tsv"
    done
    """
}
