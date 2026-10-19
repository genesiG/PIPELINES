/*
 * Trims every sample of one library through its (discovered) anchor
 * sequence (bin/guide_anchor_trim.sh — cutadapt, real alignment, tolerant
 * of a few mismatches or an indel in the anchor itself), so the guide sits
 * at position 0 of the output regardless of any variable-length spacing in
 * the amplicon. MAGECK_COUNT then runs with --trim-5 0 on these outputs
 * instead of assuming a fixed offset from the original demultiplexed reads.
 *
 * A read where the anchor can't be found is passed through UNCHANGED (no
 * --discard-untrimmed in the script), not dropped — see the script's own
 * header for why that keeps the total read count, and read classification,
 * fully accounted for.
 */
include { runner; extraArgs } from './functions'

process GUIDE_ANCHOR_TRIM {
    tag "${lib_name}"
    // command.sh excluded here, published separately below -- see
    // guide_representation.nf for why (no pattern would otherwise also match it).
    // One subfolder per library: with every sample counted against every
    // library, each library trims its own copy of the same <sample>_trimmed.fastq.gz.
    publishDir "${params.intermediate_dir}/anchor_trimmed_fastq", mode: params.publish_mode, pattern: "*_trimmed.fastq.gz", saveAs: { fn -> "${lib_name}/${fn}" }
    publishDir "${params.outdir}/counts/qc/anchor_discovery", mode: params.publish_mode, pattern: "anchor_trim_stats.csv", saveAs: { "${params.project_name}_${lib_name}_anchor_trim_stats.csv" }
    publishDir "${params.outdir}/counts/qc/anchor_discovery/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${lib_name}_trim_command.sh" }

    input:
    path fastqs   // every sample's FASTQ, staged by name; only this library's samples are used
    tuple val(lib_name), path(anchor_file), val(fastq_args), val(sample_label_str)
          // fastq_args: ordered List<String> of FASTQ file names matching sample_label_str's order
    path trim_script   // bin/guide_anchor_trim.sh, or the project's `guide_anchor_trim:` replacement

    output:
    tuple val(lib_name), path("*_trimmed.fastq.gz"), emit: trimmed_fastqs
    path "anchor_trim_stats.csv",                      emit: stats
    // cutadapt's own per-sample logs (its exact command line and version),
    // read by METHODS; optional, since a replacement script may not write them
    tuple val(lib_name), path("*_anchor_trim.log"),   emit: logs, optional: true
    path "command.sh",                                  emit: command_sh

    script:
    """
    ${runner(trim_script)} "${trim_script}" \\
        --anchor-file "${anchor_file}" \\
        --fastq "${fastq_args.join(' ')}" \\
        --labels "${sample_label_str}" \\
        --outdir . \\
        ${extraArgs('guide_anchor_trim')}
    """

    stub:
    """
    echo "sample,reads,anchor_found_reads,pct_anchor_found" > anchor_trim_stats.csv
    for n in \$(echo '${sample_label_str}' | tr ',' ' '); do
        echo stub > "\${n}_trimmed.fastq.gz"
        echo "\${n},0,0,0.00" >> anchor_trim_stats.csv
        printf 'This is cutadapt 0.0 with Python 0.0\\nCommand line parameters: -g %s -e 0.1 --overlap 6 --cores 0 -o %s_trimmed.fastq.gz %s.fastq.gz\\n' \\
            "\$(tr -d '[:space:]' < '${anchor_file}')" "\${n}" "\${n}" > "\${n}_anchor_trim.log"
    done
    """
}
