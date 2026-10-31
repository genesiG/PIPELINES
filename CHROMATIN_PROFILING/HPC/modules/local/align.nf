include { runner } from './functions'

process ALIGN {
    tag "${sample_id}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    publishDir "${params.intermediate_dir}/aligned", mode: params.publish_mode, pattern: "*.aligned.bam"
    publishDir "${params.outdir}/qc/bowtie2", mode: params.publish_mode, pattern: "*.bowtie2.log"
    publishDir "${params.outdir}/qc/bowtie2/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_align_command.sh" }

    input:
    tuple val(sample_id), path(r1), path(r2)
    path index_files   // all *.bt2 / *.bt2l files of the index prefix, staged flat
    val index_base     // bare file-name prefix of the index
    val mode           // local | end_to_end
    val minins
    val maxins
    path step_script   // bin/align.sh, or the project's `align:` replacement

    output:
    tuple val(sample_id), path("${sample_id}.aligned.bam"), emit: bam
    path "${sample_id}.bowtie2.log", emit: log
    path "command.sh",               emit: command_sh

    script:
    // -I/-X are optional: emitted only when set (null or '' = bowtie2 defaults, as in the manuscript)
    def insArgs = ((minins != null && minins.toString() != '') ? "--minins ${minins} " : '') + ((maxins != null && maxins.toString() != '') ? "--maxins ${maxins} " : '')
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --r1 "${r1}" \\
        --r2 "${r2}" \\
        --index-base "${index_base}" \\
        --mode ${mode} \\
        ${insArgs}\\
        --threads ${task.cpus} \\
        ${params.align_args ?: ''}
    """

    stub:
    """
    touch "${sample_id}.aligned.bam"
    echo "stub bowtie2 log" > "${sample_id}.bowtie2.log"
    """
}
