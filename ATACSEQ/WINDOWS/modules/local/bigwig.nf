include { runner; extraArgs } from './functions'

process BIGWIG {
    tag "${sample_id}"
    container 'quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0'
    publishDir { "${params.outdir}/bigwig" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/bigwig/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_${method}_bigwig_command.sh" }

    input:
    tuple val(sample_id), path(bam), path(bai), val(multiplier), val(method)
    val bin_size
    val units               // none | cpm | rpkm, one value for all samples
    val extend_reads        // true | false
    val center_reads        // true | false (bamCoverage --centerReads)
    val effective_genome_size
    path step_script        // bin/bigwig.sh, or the project's `bigwig:` replacement

    output:
    tuple val(sample_id), path("${sample_id}.${method.toString().toLowerCase()}.${units.toString().toLowerCase()}.bw"), emit: bw
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --multiplier "${multiplier}" \\
        --method "${method}" \\
        --bin-size ${bin_size} \\
        --units ${units} \\
        --extend-reads ${extend_reads} \\
        --center-reads ${center_reads} \\
        --effective-genome-size "${effective_genome_size ?: ''}" \\
        --threads ${task.cpus} \\
        ${extraArgs('bigwig')}
    """

    stub:
    """
    touch "${sample_id}.${method.toString().toLowerCase()}.${units.toString().toLowerCase()}.bw"
    """
}
