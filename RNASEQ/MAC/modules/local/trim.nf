include { runner; extraArgs } from './functions'

process TRIM {
    tag "${sample_id}"
    container 'quay.io/biocontainers/trimmomatic:0.40--hdfd78af_0'
    publishDir "${params.intermediate_dir}/trimmed", mode: params.publish_mode, pattern: "*.trim.fastq.gz", enabled: params.publish_intermediates.toString() == 'true'
    publishDir "${params.outdir}/qc/trimmomatic", mode: params.publish_mode, pattern: "*.trimmomatic.log"
    publishDir "${params.outdir}/qc/trimmomatic/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_trim_command.sh" }

    input:
    tuple val(sample_id), path(r1), path(r2)
    val adapter_set     // TruSeq3-PE-2 (params.trim_adapter_set)
    val illuminaclip    // ILLUMINACLIP parameters after the adapter file (params.trim_illuminaclip)
    val steps           // extra Trimmomatic steps, space separated (params.trim_steps)
    val minlen          // params.trim_minlen
    path step_script    // bin/trim.sh, or the project's `trim` script replacement

    output:
    tuple val(sample_id), path("${sample_id}.R1.trim.fastq.gz"), path("${sample_id}.R2.trim.fastq.gz"), emit: reads
    path "${sample_id}.trimmomatic.log", emit: log
    path "command.sh",                   emit: command_sh

    script:
    def java_gb = Math.max(1, (task.memory.toGiga() * 0.75) as int)   // Java heap from the task's memory
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
        --java-gb ${java_gb} \\
        ${extraArgs('trim')}
    """

    stub:
    """
    printf '@stub/1\\nA\\n+\\nI\\n' | gzip -c > "${sample_id}.R1.trim.fastq.gz"
    printf '@stub/2\\nA\\n+\\nI\\n' | gzip -c > "${sample_id}.R2.trim.fastq.gz"
    echo "stub trimmomatic log" > "${sample_id}.trimmomatic.log"
    """
}
