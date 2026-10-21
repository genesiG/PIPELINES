/*
 * Runs `mageck count` (bin/mageck_count.sh) for one library over the samples
 * assigned to it. main.nf decides which samples those are: the samples whose
 * samplesheet `library` value names this library, or — when samplesheet.csv
 * has no `library` column — every sample (the legacy step_2_count.py behaviour).
 *
 * fastq_args is the ordered filename list built in main.nf directly from
 * samplesheet.csv (not from channel emission order), so --fastq and
 * --labels always stay positionally correct.
 *
 * Also emits <prefix>.unmapped.txt (bin/mageck_count.sh always passes
 * --unmapped-to-file): every distinct sequence mageck couldn't match to a
 * guide, one column per sample, feeding GUIDE_REPRESENTATION_WRANGLE's
 * near-miss/far-miss classification.
 *
 * drop_zero (true only when every sample is counted against every library —
 * no samplesheet `library` column and more than one library) lets
 * bin/mageck_count.sh leave a sample with 0 reads in this library out of its
 * count table instead of failing: in that mode a sample holding none of one
 * library's guides is expected, not a mistake. Every run also writes
 * <prefix>.sample_reads.tsv (each sample's input/mapped reads and whether it
 * was counted), which LIBRARY_BREAKDOWN collects across libraries.
 *
 * trim5 is a task input, not read from params directly, because its right
 * value depends on whether main.nf ran GUIDE_ANCHOR_DISCOVERY/TRIM first
 * (see main.nf's autoAnchor branch): 0 when it did (the input fastqs are
 * already positioned at the guide), params.trim5 unchanged when a project
 * set it explicitly to bypass anchor discovery.
 */
include { runner } from './functions'

process MAGECK_COUNT {
    tag "${lib_name}"
    // command.sh (from conf/base.config's afterScript) also matches pattern
    // "*" below — route it, like the .log file, into log/ with a unique
    // per-library name; a second publishDir keyed on the exact filename
    // would double-publish it (once via each matching rule).
    publishDir "${params.outdir}/counts", mode: params.publish_mode, pattern: "*",
        saveAs: { fn -> fn.endsWith('.log') ? "log/${fn}" : fn == 'command.sh' ? "log/${prefix(lib_name)}_command.sh" : fn }

    input:
    tuple val(lib_name), path(mageck_library), val(fastq_args), val(sample_label_str), path(fastqs)
          // fastq_args: ordered List<String> of FASTQ file names matching sample_label_str's order
          // fastqs: staged by name — every sample FASTQ (fixed trim5), or only THIS library's
          // GUIDE_ANCHOR_TRIM output (anchor mode; each library trims its own copies, which
          // share file names across libraries, so they must not be pooled) — main.nf decides
    val drop_zero // true only in every-sample-vs-every-library mode: leave 0-mapped samples out of this table
    val trim5     // 0 when main.nf ran GUIDE_ANCHOR_DISCOVERY/TRIM first (already positioned at the guide);
                  // params.trim5 unchanged when a project set it explicitly, bypassing anchor discovery
    path step_script   // bin/mageck_count.sh, or the project's `mageck_count:` replacement

    output:
    tuple val(lib_name),
          path("${prefix(lib_name)}.count.txt"),
          path("${prefix(lib_name)}.count_normalized.txt"),
          path("${prefix(lib_name)}_mageck_count.log"),
          path("${prefix(lib_name)}.unmapped.txt"),             emit: results
    tuple val(lib_name), path("${prefix(lib_name)}.sample_reads.tsv"), emit: sample_reads
    path "${prefix(lib_name)}*",                                emit: all_outputs
    path "command.sh",                                          emit: command_sh

    script:
    // Exit 3 from the script = samples with 0 mapped reads (a data problem,
    // not retried — see conf/base.config).
    """
    ${runner(step_script)} "${step_script}" \\
        --library   "${mageck_library}" \\
        --fastq     "${fastq_args.join(' ')}" \\
        --labels    "${sample_label_str}" \\
        --prefix    "${prefix(lib_name)}" \\
        --sgrna-len ${params.sgrna_len} \\
        --trim5     ${trim5} \\
        ${drop_zero ? '--drop-zero-samples' : ''} \\
        ${params.mageck_count_args ?: ''}
    """

    stub:
    def stub_prefix = prefix(lib_name)
    """
    printf 'sgRNA\\tGene\\t%s\\n' "\$(echo '${sample_label_str}' | tr ',' '\\t')" > "${stub_prefix}.count.txt"
    cp "${stub_prefix}.count.txt" "${stub_prefix}.count_normalized.txt"
    echo stub > "${stub_prefix}_mageck_count.log"
    touch "${stub_prefix}.unmapped.txt"
    printf 'sample\\tfastq\\treads\\tmapped\\tstatus\\n' > "${stub_prefix}.sample_reads.tsv"
    for n in \$(echo '${sample_label_str}' | tr ',' ' '); do printf '%s\\tstub.fastq.gz\\t100\\t50\\tcounted\\n' "\$n" >> "${stub_prefix}.sample_reads.tsv"; done
    """
}

def prefix(lib_name) {
    "${params.project_name}_${lib_name}"
}
