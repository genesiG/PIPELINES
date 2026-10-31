// Gene model for the strandedness check (RSeQC needs BED12). Two processes, because a process has one container:
// GTF -> genePred (UCSC gtfToGenePred), then genePred -> BED12 (UCSC genePredToBed). Both run once (value channels).
include { runner; extraArgs } from './functions'

process GTF_TO_GENEPRED {
    tag "gtf_to_genepred"
    container 'quay.io/biocontainers/ucsc-gtftogenepred:482--h0b57e2e_0'
    publishDir "${params.outdir}/qc/strandedness/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "gtf_to_genepred_command.sh" }

    input:
    path gtf
    path step_script   // bin/gtf_to_genepred.sh, or the project's `gtf_to_genepred` script replacement

    output:
    path "genes.genepred", emit: genepred
    path "command.sh",     emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --gtf "${gtf}" \\
        --out genes.genepred \\
        ${extraArgs('gtf_to_genepred')}
    """

    stub:
    """
    printf 'stub_gene\\tchr1\\t+\\t0\\t100\\t0\\t0\\t1\\t0,\\t100,\\n' > genes.genepred
    """
}

process GENEPRED_TO_BED12 {
    tag "genepred_to_bed12"
    container 'quay.io/biocontainers/ucsc-genepredtobed:482--h0b57e2e_0'
    publishDir "${params.outdir}/qc/strandedness", mode: params.publish_mode, pattern: "genes.bed12"
    publishDir "${params.outdir}/qc/strandedness/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "genepred_to_bed12_command.sh" }

    input:
    path genepred
    path step_script   // bin/genepred_to_bed12.sh, or the project's `genepred_to_bed12` script replacement

    output:
    path "genes.bed12", emit: bed12
    path "command.sh",  emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --genepred "${genepred}" \\
        --out genes.bed12 \\
        ${extraArgs('genepred_to_bed12')}
    """

    stub:
    """
    printf 'chr1\\t0\\t100\\tstub_gene\\t0\\t+\\t0\\t0\\t0,0,0\\t1\\t100,\\t0,\\n' > genes.bed12
    """
}
