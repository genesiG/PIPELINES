include { runner; extraArgs } from './functions'

process FRIP {
    tag "${sample_id}.${peak_source}"
    container 'quay.io/biocontainers/bedtools:2.31.1--h13024bc_3'
    // BAM input: no cache 'deep'
    publishDir "${params.outdir}/qc/frip", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/frip/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}.${peak_source}_frip_command.sh" }

    input:
    // peaks_bed = the sample's OWN peak set; final_flagstat = S.final.flagstat.txt of THIS sample (the script exits 3 if it belongs to another sample)
    tuple val(sample_id), path(bam), path(bai), path(peaks_bed), val(peak_source), path(final_flagstat)
    path step_script   // bin/qc_frip.sh, or the project's `qc_frip:` replacement

    output:
    path "${sample_id}.${peak_source}.frip.tsv", emit: tsv
    path "command.sh",                           emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --peaks "${peaks_bed}" \\
        --peak-source "${peak_source}" \\
        --flagstat "${final_flagstat}" \\
        ${extraArgs('qc_frip')}
    """

    stub:
    """
    printf 'sample_id\\tpeak_source\\tpeaks_file\\tn_peaks\\ttotal_mapped_alignments\\treads_in_peaks\\tfrip\\n${sample_id}\\t${peak_source}\\tpeaks.bed\\t100\\t1000\\t200\\t0.200000\\n' > "${sample_id}.${peak_source}.frip.tsv"
    """
}
