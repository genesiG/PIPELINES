include { runner; extraArgs } from './functions'

process TRIM {
    tag "${sample_id}"
    container 'quay.io/biocontainers/trimmomatic:0.40--hdfd78af_0'
    publishDir "${params.intermediate_dir}/trimmed", mode: params.publish_mode, pattern: "*.trim.fastq.gz"
    publishDir "${params.outdir}/qc/trimmomatic", mode: params.publish_mode, pattern: "*.trimmomatic.log"
    publishDir "${params.outdir}/qc/trimmomatic/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_trim_command.sh" }

    input:
    tuple val(sample_id), path(r1), path(r2)
    val adapter_set     // NexteraPE-PE (params.trim_adapter_set)
    val illuminaclip    // ILLUMINACLIP parameters after the adapter file (params.trim_illuminaclip)
    val steps           // extra Trimmomatic steps, space separated (params.trim_steps)
    val minlen
    path step_script    // bin/trim.sh, or the project's `trim:` replacement

    output:
    tuple val(sample_id), path("${sample_id}.R1.trim.fastq.gz"), path("${sample_id}.R2.trim.fastq.gz"), emit: reads
    path "${sample_id}.trimmomatic.log", emit: log
    path "command.sh",                   emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --r1 "${r1}" \\
        --r2 "${r2}" \\
        --adapter-set "${adapter_set}" \\
        --illuminaclip "${illuminaclip}" \\
        --steps "${steps ?: ''}" \\
        --minlen ${minlen} \\
        --threads ${task.cpus} \\
        ${extraArgs('trim')}
    """

    stub:
    """
    touch "${sample_id}.R1.trim.fastq.gz" "${sample_id}.R2.trim.fastq.gz"
    echo "stub trimmomatic log" > "${sample_id}.trimmomatic.log"
    """
}
