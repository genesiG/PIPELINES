include { runner; extraArgs } from './functions'

/*
 * Union of the MACS3 peak files of several peak groups (peak_merges.csv row: name, peak_groups).
 * peak_files: the peak files of this merge, staged flat; their names are passed to the script.
 */
process PEAK_MERGE {
    tag "${name}"
    container 'quay.io/biocontainers/bedops:2.4.42--hd6d6fdc_1'
    publishDir "${params.outdir}/peaks/merged", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/peaks/merged/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${name}_peak_merge_command.sh" }

    input:
    tuple val(name), path(peak_files)   // the peak files of this merge (MACS3 outputs), staged flat
    path step_script                    // bin/peak_merge.sh, or the project's `peak_merge:` replacement

    output:
    tuple val(name), path("${name}.merged.bed"), emit: bed
    path "command.sh",                           emit: command_sh

    script:
    def pnames = [peak_files].flatten().collect { it.name }.join(',')
    """
    ${runner(step_script)} "${step_script}" \\
        --name "${name}" \\
        --peaks "${pnames}" \\
        ${extraArgs('peak_merge')}
    """

    stub:
    """
    printf 'chr1\\t100\\t200\\n' > "${name}.merged.bed"
    """
}
