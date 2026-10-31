/*
 * Normalisation QC tables for one norm_run (bin/norm_qc.R): MA pairs (replicate and treatment), one reference line per
 * method, numeric cloud summary. Reads the shared BIN_COUNTS table of tmm_bin_size (no BAM) and the run's wide scale-factor
 * table <norm_run>.scale_factors_all.tsv (NORM_COMBINE). ref_condition is resolved in main.nf (norm_runs.csv ref_condition,
 * else the classifications' ref_group condition, else the first library's condition).
 */
include { runner } from './functions'

process NORM_QC {
    tag "${norm_run}"
    cache 'deep'
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    publishDir { "${params.outdir}/norm_qc" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/norm_qc/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${norm_run}_norm_qc_command.sh" }

    input:
    tuple val(norm_run), val(sample_ids), val(ref_condition), val(primary), path(scale_factors_all), path(counts), path(libsizes)
    path samplesheet
    path step_script   // bin/norm_qc.R

    output:
    path "${norm_run}.norm_qc.tsv", emit: tsv
    path "${norm_run}.norm_qc_pairs.tsv", emit: pairs
    path "${norm_run}.norm_qc_libsizes.tsv", emit: libsizes
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --counts=${counts} \\
        --libsizes=${libsizes} \\
        --sample-ids=${sample_ids} \\
        --samplesheet=${samplesheet} \\
        --scale-factors-all=${scale_factors_all} \\
        --primary=${primary} \\
        --ref-condition="${ref_condition}" \\
        --min-count=${params.norm_qc_min_count} \\
        --prior=${params.norm_qc_prior} \\
        --out-prefix=${norm_run} \\
        ${params.norm_qc_args ?: ''}
    """

    stub:
    """
    printf 'chr\\tstart\\tend\\n' > "${norm_run}.norm_qc.tsv"
    printf 'pair_id\\tpair_type\\tref\\tother\\tref_group\\tother_group\\tmethod\\tis_primary\\tstatus\\tline\\tbigwig_offset\\tline_se\\tlog2_factor_ratio_check\\tn_bins\\tmin_count\\tprior\\tA_peak\\tM_peak\\tM_centroid\\tcore_lo\\tcore_hi\\tline_minus_centroid\\n' > "${norm_run}.norm_qc_pairs.tsv"
    printf 'sample_id\\tcsaw_total\\tbam_mapped_reads\\tbam_fragments\\n' > "${norm_run}.norm_qc_libsizes.tsv"
    for id in \$(echo "${sample_ids}" | tr ',' ' '); do
        printf '%s\\t60\\t120\\t60\\n' "\$id" >> "${norm_run}.norm_qc_libsizes.tsv"
    done
    """
}
