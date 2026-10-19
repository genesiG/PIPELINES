/*
 * Runs one `mageck test` (bin/mageck_test.sh) per (treatment, control)
 * contrast, fanned out as independent tasks, each against the count table
 * VALIDATE_CONTRASTS resolved it to (the library whose samples include both sides).
 */
include { runner; extraArgs } from './functions'

process MAGECK_TEST {
    tag "${name}"
    // command.sh excluded here (via saveAs, not pattern -- publishDir's
    // pattern: does NOT support a closure the way its destination path and
    // saveAs: do; a closure there silently matches nothing — see
    // guide_representation_report.nf), published separately below instead.
    publishDir { "${params.outdir}/mageck_test/${name}" }, mode: params.publish_mode,
        saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/mageck_test/${name}/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${name}_command.sh" }

    input:
    tuple val(table_name), val(name), val(treat_cols), val(control_cols), path(count_table)
          // name: <treatment>_vs_<control>, file-name safe (validate_contrasts.py)
    path step_script   // bin/mageck_test.sh, or the project's `mageck_test:` replacement

    output:
    path "${name}*", emit: results
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --count-table "${count_table}" \\
        --name "${name}" \\
        --treatment "${treat_cols}" \\
        --control "${control_cols}" \\
        ${extraArgs('mageck_test')}
    """

    stub:
    """
    echo "stub: ${count_table} | ${treat_cols} vs ${control_cols}" > "${name}.gene_summary.txt"
    """
}
