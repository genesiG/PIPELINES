include { runner; extraArgs } from './functions'

process MULTIQC {
    tag "multiqc"
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/qc/multiqc", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/multiqc/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "multiqc_command.sh" }

    input:
    path '*'            // collected QC files, staged flat: FastQC zips (raw and trimmed), Trimmomatic logs, HISAT2 summaries,
                        // raw and filtered flagstat/idxstats, Picard metrics, featureCounts summaries, RSeQC infer_experiment
                        // files; MultiQC auto-detects
    path step_script    // bin/multiqc.sh, or the project's `multiqc:` replacement

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_report_data", emit: data   // MultiQC derives the data dir name from --filename
    path "command.sh",          emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" ${extraArgs('multiqc')}
    """

    stub:
    """
    touch multiqc_report.html
    mkdir -p multiqc_report_data
    """
}
