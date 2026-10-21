/*
 * Barcode demultiplexing: one multiplexed FASTQ (or R1/R2 pair) -> one FASTQ
 * per sample, plus the barcode-position CSV behind the demux_plots density
 * plot. The work is in bin/demultiplex.sh (cutadapt, then
 * bin/barcode_position_qc.py on cutadapt's --info-file inside the same task,
 * so the multi-GB info file never crosses a process boundary).
 *
 * No `cpus` directive: the LSF request is set only via clusterOptions in
 * conf/base.config (cpus + clusterOptions gives a duplicate span[] that LSF
 * rejects), and the script runs cutadapt with --cores 0.
 */
include { runner } from './functions'

process DEMULTIPLEX_SE {
    tag "${params.project_name}"
    container 'quay.io/biocontainers/cutadapt:5.2--py39hbcbf7aa_0'
    // Per-sample FASTQs land in original_data/fastq/ (not results/ or
    // results/intermediate_files/) — by lab convention that's where every pipeline
    // expects to find per-sample FASTQs, whether they arrived there via a
    // demux step (this pipeline) or were placed there directly (pipelines
    // with no demux step).
    publishDir "${params.original_data_dir}/fastq", mode: params.publish_mode, pattern: "*.fastq.gz"
    publishDir "${params.outdir}/demultiplexing/log", mode: params.publish_mode, pattern: "cutadapt.log", saveAs: { "${params.project_name}_step1_cutadapt.log" }
    publishDir "${params.outdir}/demultiplexing/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${params.project_name}_step1_command.sh" }
    publishDir "${params.outdir}/demultiplexing/qc/demux_plots", mode: params.publish_mode, pattern: "*_barcode_position_counts.csv"

    input:
    path fastq_r1
    val  adapter_args
    val  barcode_len
    path demux_script    // bin/demultiplex.sh, or the project's `demultiplex:` replacement
    path qc_script       // bin/barcode_position_qc.py, or the project's `barcode_position_qc:` replacement

    output:
    path "*.fastq.gz",                  emit: fastqs
    path "cutadapt.log",                 emit: log
    path "*_barcode_position_counts.csv", emit: position_csv
    path "command.sh",                    emit: command_sh

    script:
    """
    ${runner(demux_script)} "${demux_script}" \\
        --r1 "${fastq_r1}" \\
        --adapters "${adapter_args}" \\
        --overlap ${barcode_len} \\
        --position-script "${qc_script}" \\
        --project "${params.project_name}" \\
        ${params.demultiplex_args ?: ''}
    """

    stub:
    """
    for n in \$(echo '${adapter_args}' | grep -o -- '-g [^=]*' | cut -c4-); do echo stub > "\${n}.fastq.gz"; done
    echo stub > cutadapt.log
    touch "${params.project_name}_barcode_position_counts.csv"
    """
}

process DEMULTIPLEX_PE {
    tag "${params.project_name}"
    container 'quay.io/biocontainers/cutadapt:5.2--py39hbcbf7aa_0'
    publishDir "${params.original_data_dir}/fastq", mode: params.publish_mode, pattern: "*.fastq.gz"
    publishDir "${params.outdir}/demultiplexing/log", mode: params.publish_mode, pattern: "cutadapt.log", saveAs: { "${params.project_name}_step1_cutadapt.log" }
    publishDir "${params.outdir}/demultiplexing/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${params.project_name}_step1_command.sh" }
    publishDir "${params.outdir}/demultiplexing/qc/demux_plots", mode: params.publish_mode, pattern: "*_barcode_position_counts.csv"

    input:
    path fastq_r1
    path fastq_r2
    val  adapter_args
    val  barcode_len
    path demux_script    // bin/demultiplex.sh, or the project's `demultiplex:` replacement
    path qc_script       // bin/barcode_position_qc.py, or the project's `barcode_position_qc:` replacement

    output:
    path "*.fastq.gz",                  emit: fastqs
    path "cutadapt.log",                 emit: log
    path "*_barcode_position_counts.csv", emit: position_csv
    path "command.sh",                    emit: command_sh

    script:
    """
    ${runner(demux_script)} "${demux_script}" \\
        --r1 "${fastq_r1}" \\
        --r2 "${fastq_r2}" \\
        --adapters "${adapter_args}" \\
        --overlap ${barcode_len} \\
        --position-script "${qc_script}" \\
        --project "${params.project_name}" \\
        ${params.demultiplex_args ?: ''}
    """

    stub:
    """
    for n in \$(echo '${adapter_args}' | grep -o -- '-g [^=]*' | cut -c4-); do echo stub > "\${n}_R1.fastq.gz"; echo stub > "\${n}_R2.fastq.gz"; done
    echo stub > cutadapt.log
    touch "${params.project_name}_barcode_position_counts.csv"
    """
}
