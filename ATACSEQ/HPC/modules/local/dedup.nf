include { runner; extraArgs } from './functions'

process DEDUP {
    tag "${sample_id}"
    container 'quay.io/biocontainers/sambamba:1.0.1--he614052_4'
    // The dedup BAM is an intermediate read by BLACKLIST_FILTER; only the sambamba log is published.
    publishDir "${params.outdir}/bam/stats", mode: params.publish_mode, pattern: "*.sambamba_markdup.log"
    publishDir "${params.outdir}/bam/stats/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_dedup_command.sh" }

    input:
    tuple val(sample_id), path(bam), path(bai)
    val remove          // true = remove duplicates (sambamba -r); false = mark only
    path step_script    // bin/dedup.sh, or the project's `dedup:` replacement

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), emit: bam
    path "${sample_id}.sambamba_markdup.log",             emit: log
    path "command.sh",                                    emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --remove ${remove} \\
        --threads ${task.cpus} \\
        ${extraArgs('dedup')}
    """

    stub:
    """
    touch "${sample_id}.dedup.bam"
    echo "stub sambamba markdup log" > "${sample_id}.sambamba_markdup.log"
    """
}
