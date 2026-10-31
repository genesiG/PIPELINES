/*
 * MACS3 callpeak of one group. The mode is per group: `broad` ('true' | 'false') comes from main.nf, which decides it from the
 * group's antibody (broad_marks, macs3_broad = auto) or forces it (macs3_broad = true | false); broad adds --broad
 * --broad-cutoff (macs3_broad_cutoff) and gives <group>_peaks.broadPeak, narrow gives <group>_peaks.narrowPeak.
 */
include { runner } from './functions'

process MACS3 {
    tag "${group}"
    container 'quay.io/biocontainers/macs3:3.0.4--py312h71493bf_0'
    publishDir { "${params.outdir}/peaks/macs3" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/peaks/macs3/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${group}_macs3_command.sh" }

    input:
    tuple val(group), path(target_bams), path(control_bams), val(broad)
    val genome_size
    val qvalue
    val broad_cutoff
    path step_script        // bin/callpeaks_macs3.sh, or the project's `callpeaks_macs3:` replacement

    output:
    tuple val(group), path("${group}_peaks.{narrow,broad}Peak"), emit: peaks
    path "${group}_peaks.xls",   emit: xls
    path "${group}.macs3.log",   emit: log
    path "${group}_summits.bed", optional: true, emit: summits
    path "${group}_peaks.gappedPeak", optional: true, emit: gapped
    path "command.sh", emit: command_sh

    script:
    def tnames = [target_bams].flatten().collect { it.name }.join(',')
    def cnames = [control_bams].flatten().collect { it.name }.join(',')
    """
    ${runner(step_script)} "${step_script}" \\
        --group "${group}" \\
        --targets "${tnames}" \\
        --controls "${cnames}" \\
        --genome-size ${genome_size} \\
        --qvalue ${qvalue} \\
        --broad ${broad} \\
        --broad-cutoff ${broad_cutoff} \\
        ${params.callpeaks_macs3_args ?: ''}
    """

    stub:
    def ext = broad.toString() == 'true' ? 'broadPeak' : 'narrowPeak'
    """
    touch "${group}_peaks.${ext}" "${group}_peaks.xls" "${group}.macs3.log"
    """
}
