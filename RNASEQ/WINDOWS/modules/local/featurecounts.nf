include { runner; extraArgs } from './functions'

process FEATURECOUNTS {
    tag "${sample_id}"
    container 'quay.io/biocontainers/subread:2.1.1--h577a1d6_0'
    publishDir "${params.outdir}/counts/per_sample", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/counts/per_sample/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_featurecounts_command.sh" }

    input:
    tuple val(sample_id), path(bam), path(bai)
    path gtf
    val strand_code     // featureCounts -s: 0 unstranded, 1 forward, 2 reverse (main.nf maps params.strandedness)
    val feature_type    // params.featurecounts_feature_type (-t)
    val attribute       // params.featurecounts_attribute (-g)
    path step_script    // bin/featurecounts.sh, or the project's `featurecounts` script replacement

    output:
    path "${sample_id}_counts.txt",                 emit: counts
    path "${sample_id}_featurecounts.txt",          emit: table
    path "${sample_id}_featurecounts.txt.summary",  emit: summary
    path "command.sh",                              emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --gtf "${gtf}" \\
        --strand ${strand_code} \\
        --feature-type "${feature_type}" \\
        --attribute "${attribute}" \\
        --threads ${task.cpus} \\
        ${extraArgs('featurecounts')}
    """

    stub:
    """
    printf '# Program:featureCounts v2.1.1\\n# Command:stub\\nGeneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\t${sample_id}.final.bam\\ngene1\\tchr1\\t1\\t100\\t+\\t100\\t5\\n' > "${sample_id}_featurecounts.txt"
    printf 'Geneid\\t${sample_id}.final.bam\\ngene1\\t5\\n' > "${sample_id}_counts.txt"
    printf 'Status\\t${sample_id}.final.bam\\nAssigned\\t5\\n' > "${sample_id}_featurecounts.txt.summary"
    """
}
