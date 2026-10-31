/*
 * csaw-style TMM scale factors for one norm_run (bin/norm_tmm.R, counts mode): the factors are computed from the shared
 * BIN_COUNTS table (counts + libsizes of the run's bin size), so no BAM is read.
 * sample_ids: comma-separated, run order (a subset of the counts columns). role: 'primary' (non-finite factors: exit 3) or
 * 'diagnostic' (--fail-soft=TRUE: an UNAVAILABLE table instead). bin_size is the key of the BIN_COUNTS output the task gets.
 */
include { runner } from './functions'

process TMM {
    tag "${norm_run}"
    cache 'deep'
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    publishDir { "${params.outdir}/scale_factors" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn.endsWith('.scale_factors.tsv') ? "methods/${fn}" : fn }
    publishDir { "${params.outdir}/scale_factors/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${norm_run}_tmm_command.sh" }

    input:
    tuple val(bin_size), val(norm_run), val(sample_ids), val(role), path(counts), path(libsizes)
    path step_script   // bin/norm_tmm.R

    output:
    tuple val(norm_run), val('tmm'), path("${norm_run}.tmm.scale_factors.tsv"), emit: sf
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \
        --counts=${counts} \
        --libsizes=${libsizes} \
        --sample-ids=${sample_ids} \
        --norm-run=${norm_run} \
        --bin-size=${bin_size} \
        --out=${norm_run}.tmm.scale_factors.tsv \
        ${role == 'diagnostic' ? '--fail-soft=TRUE' : ''} \
        ${params.norm_tmm_args ?: ''}
    """

    stub:
    """
    printf 'sample_id\tmethod\traw_factor\tbigwig_multiplier\tedger_norm_factor\tqc_pass\tnote\n' > "${norm_run}.tmm.scale_factors.tsv"
    for id in \$(echo "${sample_ids}" | tr ',' ' '); do
        printf '%s\ttmm\t1\t1\t1\tTRUE\tstub\n' "\$id" >> "${norm_run}.tmm.scale_factors.tsv"
    done
    """
}
