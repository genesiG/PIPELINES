include { runner; extraArgs } from './functions'

process INFER_STRANDEDNESS {
    tag "${sample_id}"
    container 'quay.io/biocontainers/rseqc:5.0.6--pyhdfd78af_0'
    publishDir "${params.outdir}/qc/strandedness", mode: params.publish_mode, pattern: "*.infer_experiment.txt"
    publishDir "${params.outdir}/qc/strandedness/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}_infer_strandedness_command.sh" }

    input:
    tuple val(sample_id), path(bam), path(bai)
    path bed           // genes.bed12 (GENEPRED_TO_BED12)
    val reads          // params.infer_experiment_reads: reads sampled by infer_experiment.py -s
    path step_script   // bin/infer_strandedness.sh, or the project's `infer_strandedness` script replacement

    output:
    path "${sample_id}.infer_experiment.txt", emit: txt
    path "command.sh",                        emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --bed "${bed}" \\
        --reads ${reads} \\
        ${extraArgs('infer_strandedness')}
    """

    stub:
    """
    cat > "${sample_id}.infer_experiment.txt" <<'EOF'
This is PairEnd Data
Fraction of reads failed to determine: 0.0100
Fraction of reads explained by "1++,1--,2+-,2-+": 0.9000
Fraction of reads explained by "1+-,1-+,2++,2--": 0.0900
EOF
    """
}
