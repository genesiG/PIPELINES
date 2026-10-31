/*
 * SEACR of one sample against its control. LAB POLICY: SEACR has no broad mode. relaxed = threshold between the knee and the
 * peak of the total-signal curve (more, wider peaks), stringent = the peak of the curve. The mode is per sample (main.nf):
 * seacr_mode = auto gives relaxed for the marks of broad_marks and stringent for the others; relaxed / stringent force it.
 */
include { runner } from './functions'

process SEACR_CALL {
    tag "${sample_id}"
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'   // SEACR 1.3 + R + bedtools (the seacr biocontainer lacks bedtools)
    publishDir { "${params.outdir}/peaks/seacr" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/peaks/seacr/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_seacr_call_command.sh" }

    input:
    tuple val(sample_id), path(target_bg), path(control_bg), val(mode)   // mode: relaxed | stringent
    path step_script        // bin/seacr_call.sh, or the project's replacement

    output:
    tuple val(sample_id), path("${sample_id}.seacr.${mode}.bed"), emit: bed
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --target-bg "${target_bg}" \\
        --control-bg "${control_bg}" \\
        --mode ${mode} \\
        --norm norm \\
        ${params.seacr_call_args ?: ''}
    """

    stub:
    """
    touch "${sample_id}.seacr.${mode}.bed"
    """
}
