/*
 * Replicate concordance of libraries of the same samplesheet group (bin/qc_replicate_correlation.R), from the shared
 * BIN_COUNTS table of repcor_bin_size (no BAM). Descriptive values, no pass/fail threshold. Gated by run_replicate_correlation.
 * sample_ids: comma list of the active libraries.
 */
include { runner } from './functions'

process REPLICATE_CORRELATION {
    tag "replicate_correlation"
    cache 'deep'
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    publishDir "${params.outdir}/qc/replicate_correlation", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/replicate_correlation/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "replicate_correlation_command.sh" }

    input:
    tuple val(bin_size), path(counts), path(libsizes)
    path samplesheet
    val sample_ids
    path step_script    // bin/qc_replicate_correlation.R, or the project's `qc_replicate_correlation:` replacement

    output:
    path "replicate_correlation.tsv", emit: table
    path "replicate_correlation_pairs_all.tsv", emit: pairs_all
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --counts=${counts} \\
        --libsizes=${libsizes} \\
        --samplesheet=${samplesheet} \\
        --sample-ids=${sample_ids} \\
        --prior=${params.repcor_prior} \\
        --signal-frac=${params.repcor_signal_frac} \\
        --bin-size=${bin_size} \\
        --out-prefix=replicate_correlation \\
        --outdir=. \\
        ${params.qc_replicate_correlation_args ?: ''}
    """

    stub:
    """
    printf 'group\\tantibody\\tcondition\\tsample_1\\tsample_2\\tstatus\\tn_bins_all\\tpearson_all\\tspearman_all\\tn_bins_signal\\tmedian_pair_count_signal\\tpearson_signal\\trobust_sd_log2ratio_signal\\tnonrep_pearson_all_median\\tnonrep_n\\tbin_size\\tprior\\tsignal_frac\\n' > replicate_correlation.tsv
    printf 'stub_group\\tstub\\tstub\\tstub_1\\tstub_2\\tok\\t3\\t1\\t1\\t1\\t10\\t1\\t0\\tNA\\t0\\t${bin_size}\\t${params.repcor_prior}\\t${params.repcor_signal_frac}\\n' >> replicate_correlation.tsv
    printf 'sample_1\\tsample_2\\tantibody\\tsame_group\\tpearson_all\\tspearman_all\\tpearson_signal\\n' > replicate_correlation_pairs_all.tsv
    printf 'stub_1\\tstub_2\\tstub\\tTRUE\\t1\\t1\\t1\\n' >> replicate_correlation_pairs_all.tsv
    """
}
