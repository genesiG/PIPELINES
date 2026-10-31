/*
 * Combines the per-method scale-factor tables of one norm_run (bin/norm_combine.py). Publishes
 * <norm_run>.scale_factors.tsv (byte copy of the PRIMARY method's table, the path and content every consumer knows) and
 * <norm_run>.scale_factors_all.tsv (all methods side by side).
 * method_tsvs: the staged <norm_run>.<method>.scale_factors.tsv files; the method is read from the file name.
 * expected: comma list of the methods requested for the run; unavailable: list of 'method=reason' for methods that were
 * not launched (written as UNAVAILABLE rows).
 */
include { runner } from './functions'

process NORM_COMBINE {
    tag "${norm_run}"
    cache 'deep'
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    publishDir { "${params.outdir}/scale_factors" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/scale_factors/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${norm_run}_norm_combine_command.sh" }

    input:
    tuple val(norm_run), val(primary), val(expected), val(unavailable), path(method_tsvs)
    path step_script   // bin/norm_combine.py

    output:
    tuple val(norm_run), path("${norm_run}.scale_factors.tsv"), emit: primary
    tuple val(norm_run), path("${norm_run}.scale_factors_all.tsv"), emit: wide
    path "command.sh", emit: command_sh

    script:
    def pairs = [method_tsvs].flatten().collect { f -> "${f.name.substring(norm_run.length() + 1).replace('.scale_factors.tsv', '')}=${f.name}".toString() }.join(' ')
    def unavail = [unavailable].flatten().findAll { u -> u }.collect { u -> "'${u}'".toString() }.join(' ')
    """
    ${runner(step_script)} "${step_script}" \
        --run ${norm_run} \
        --primary ${primary} \
        --methods-tsv ${pairs} \
        --expected-methods ${expected} \
        ${unavail ? "--unavailable ${unavail}" : ''} \
        --out-prefix ${norm_run} \
        ${params.norm_combine_args ?: ''}
    """

    stub:
    """
    cp "${norm_run}.${primary}.scale_factors.tsv" "${norm_run}.scale_factors.tsv"
    methods="chipseqspikeinfree tmm cpm spikein"
    { printf 'sample_id\tprimary_method'
      for m in \$methods; do printf '\t%s_raw_factor\t%s_bigwig_multiplier\t%s_edger_norm_factor\t%s_qc_pass\t%s_status\t%s_note' \$m \$m \$m \$m \$m \$m; done; printf '\n'
      tail -n +2 "${norm_run}.scale_factors.tsv" | cut -f1 | while read -r id; do
          printf '%s\t%s' "\$id" "${primary}"
          for m in \$methods; do
              if [ -e "${norm_run}.\$m.scale_factors.tsv" ]; then printf '\t1\t1\t1\tTRUE\tok\tstub'; else printf '\tNA\tNA\tNA\tNA\tnot_run\tnot requested'; fi
          done; printf '\n'
      done; } > "${norm_run}.scale_factors_all.tsv"
    """
}
