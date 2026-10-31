/*
 * Spike-in scale factors for one norm_run (bin/norm_spikein.py).
 * target_idxstats / spike_idxstats: the staged idxstats files of the run's samples, in sample_ids order
 * (targets S.final.idxstats.txt, spike-in S.spikein.idxstats.txt, so names never clash). The script counts
 * reads from these text files, so it needs no samtools.
 * role: 'primary' (S = 0 stops the run, exit 3) or 'diagnostic' (--fail-soft: an UNAVAILABLE table instead).
 * Output <norm_run>.spikein.scale_factors.tsv, published to scale_factors/methods/.
 */
include { runner } from './functions'

process SPIKEIN {
    tag "${norm_run}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    publishDir { "${params.outdir}/scale_factors" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn.endsWith('.scale_factors.tsv') ? "methods/${fn}" : fn }
    publishDir { "${params.outdir}/scale_factors/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${norm_run}_spikein_command.sh" }

    input:
    tuple val(norm_run), val(sample_ids), val(role), path(target_idxstats), path(spike_idxstats)
    val spike_prefix
    path step_script   // bin/norm_spikein.py

    output:
    tuple val(norm_run), val('spikein'), path("${norm_run}.spikein.scale_factors.tsv"), emit: sf
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --target-idxstats ${target_idxstats} \\
        --spike-idxstats ${spike_idxstats} \\
        --sample-ids ${sample_ids.replace(',', ' ')} \\
        --norm-run ${norm_run} \\
        --out ${norm_run}.spikein.scale_factors.tsv \\
        --spike-prefix "${spike_prefix}" \\
        ${role == 'diagnostic' ? '--fail-soft' : ''} \\
        ${params.norm_spikein_args ?: ''}
    """

    stub:
    """
    printf 'sample_id\\tmethod\\traw_factor\\tbigwig_multiplier\\tedger_norm_factor\\tqc_pass\\tnote\\n' > "${norm_run}.spikein.scale_factors.tsv"
    for id in \$(echo "${sample_ids}" | tr ',' ' '); do
        printf '%s\\tspikein\\t1\\t1\\t1\\tTRUE\\tstub\\n' "\$id" >> "${norm_run}.spikein.scale_factors.tsv"
    done
    """
}
