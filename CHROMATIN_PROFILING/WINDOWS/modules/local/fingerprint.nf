include { runner } from './functions'

process FINGERPRINT {
    tag "${group}"
    container 'quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0'
    // BAM inputs: no cache 'deep'
    publishDir "${params.outdir}/qc/fingerprint", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/fingerprint/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${group}_fingerprint_command.sh" }

    input:
    tuple val(group), val(labels), path(bams), path(bais)   // labels: comma-separated, same order as bams
    path step_script                                         // bin/qc_fingerprint.sh, or the project's `qc_fingerprint:` replacement

    output:
    path "fingerprint_${group}.svg", emit: svg
    path "fingerprint_${group}.tsv", emit: tsv
    path "command.sh",               emit: command_sh

    script:
    def bam_list = (bams instanceof Collection ? bams : [bams]).collect { it.name }.join(',')
    """
    ${runner(step_script)} "${step_script}" \\
        --group "${group}" \\
        --bams "${bam_list}" \\
        --labels "${labels}" \\
        --threads ${task.cpus} \\
        ${params.qc_fingerprint_args ?: ''}
    """

    stub:
    """
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "fingerprint_${group}.svg"
    printf '#plotFingerprint --outRawCounts\\nchr\\tstart\\tend\\tstub1\\nchr1\\t0\\t500\\t1\\n' > "fingerprint_${group}.tsv"
    """
}
