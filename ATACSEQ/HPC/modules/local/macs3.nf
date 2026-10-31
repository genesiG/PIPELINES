/*
 * MACS3 callpeak for one peak group: replicates POOLED, NO control (legacy step_4c_macs2.py).
 * Output prefix = legacy out_name: {group}.qval{q}[.minlen{l}][.maxgap{g}][.fe{fe}][.broad]
 * (a part is left out when its macs3_* param is null/off; fe only when > 1). The helper reads params.macs3_*
 * so the names used in the outputs and the stub match the prefix passed to callpeaks_macs3.sh (--prefix).
 * Numeric params are printed as given (0.05, 250, 2.0).
 */
include { runner; extraArgs } from './functions'

// true when a macs3 param is set (not null, not blank)
def macs3Set(v) {
    return v != null && v.toString().trim() != ''
}

def macs3Prefix(group) {
    def p = "${group}.qval${params.macs3_qvalue}"
    if (macs3Set(params.macs3_min_length)) p += ".minlen${params.macs3_min_length}"
    if (macs3Set(params.macs3_max_gap)) p += ".maxgap${params.macs3_max_gap}"
    if (macs3Set(params.macs3_fe_cutoff) && (params.macs3_fe_cutoff.toString() as BigDecimal) > 1) p += ".fe${params.macs3_fe_cutoff}"
    if (String.valueOf(params.macs3_broad).toLowerCase() == 'true') p += '.broad'
    return p
}

process MACS3 {
    tag "${group}"
    container 'quay.io/biocontainers/macs3:3.0.4--py312h71493bf_0'
    publishDir "${params.outdir}/peaks/macs3", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/peaks/macs3/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${group}_macs3_command.sh" }

    input:
    tuple val(group), path(target_bams)     // final (or tn5) BAMs of the peak group, staged flat; no control
    val genome_size                         // effective_genome_size (macs3 -g)
    val qvalue                              // macs3_qvalue
    val min_length                          // macs3_min_length (null = off)
    val max_gap                             // macs3_max_gap (null = off)
    val fe_cutoff                           // macs3_fe_cutoff (<= 1 = off)
    val call_summits                        // macs3_call_summits: true | false
    val broad                               // macs3_broad: true | false
    val broad_cutoff                        // macs3_broad_cutoff
    path step_script                        // bin/callpeaks_macs3.sh, or the project's `callpeaks_macs3:` replacement

    output:
    tuple val(group), path("${macs3Prefix(group)}_peaks.{narrow,broad}Peak"), emit: peaks
    path "${macs3Prefix(group)}_peaks.xls",         emit: xls
    path "${macs3Prefix(group)}.macs3.log",         emit: log
    path "${macs3Prefix(group)}_summits.bed",       optional: true, emit: summits
    path "${macs3Prefix(group)}_peaks.gappedPeak",  optional: true, emit: gapped
    path "command.sh",                              emit: command_sh

    script:
    def prefix = macs3Prefix(group)
    def tnames = [target_bams].flatten().collect { it.name }.join(',')
    """
    ${runner(step_script)} "${step_script}" \\
        --group "${group}" \\
        --prefix "${prefix}" \\
        --targets "${tnames}" \\
        --genome-size ${genome_size} \\
        --qvalue ${qvalue} \\
        --min-length "${min_length ?: ''}" \\
        --max-gap "${max_gap ?: ''}" \\
        --fe-cutoff "${fe_cutoff ?: ''}" \\
        --call-summits ${call_summits} \\
        --broad ${broad} \\
        --broad-cutoff "${broad_cutoff ?: ''}" \\
        ${params.macs3_args ?: ''} \\
        ${extraArgs('callpeaks_macs3')}
    """

    stub:
    def prefix = macs3Prefix(group)
    def is_broad = String.valueOf(params.macs3_broad).toLowerCase() == 'true'
    def ext = is_broad ? 'broadPeak' : 'narrowPeak'
    def extra_file = ''
    if (is_broad) {
        extra_file = "${prefix}_peaks.gappedPeak"
    } else if (String.valueOf(params.macs3_call_summits).toLowerCase() == 'true') {
        extra_file = "${prefix}_summits.bed"
    }
    """
    touch "${prefix}_peaks.${ext}" "${prefix}_peaks.xls" "${prefix}.macs3.log"
    if [ -n "${extra_file}" ]; then touch "${extra_file}"; fi
    """
}
