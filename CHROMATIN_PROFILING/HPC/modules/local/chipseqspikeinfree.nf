/*
 * ChIPseqSpikeInFree scale factors for one norm_run (bin/norm_chipseqspikeinfree.R).
 * bam_names: comma-separated staged BAM file names in sample_ids order.
 * role: 'primary' (the run's method: failures stop the run, exit 3) or 'diagnostic' (--fail-soft=TRUE: a package failure
 * writes an UNAVAILABLE table instead). The method table is <norm_run>.chipseqspikeinfree.scale_factors.tsv, published to
 * scale_factors/methods/; the MA svg stays in scale_factors/ (it is absent when a diagnostic run is UNAVAILABLE).
 */
include { runner } from './functions'

process CHIPSEQSPIKEINFREE {
    tag "${norm_run}"
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    publishDir { "${params.outdir}/scale_factors" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn.endsWith('.scale_factors.tsv') ? "methods/${fn}" : fn }
    publishDir { "${params.outdir}/scale_factors/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${norm_run}_chipseqspikeinfree_command.sh" }

    input:
    tuple val(norm_run), val(sample_ids), val(bam_names), path(bams), path(bais), val(role), val(bin_size), val(cutoff), val(max_turns), val(min_first_turn)
    path chrom_sizes
    path step_script   // bin/norm_chipseqspikeinfree.R

    output:
    tuple val(norm_run), val('chipseqspikeinfree'), path("${norm_run}.chipseqspikeinfree.scale_factors.tsv"), emit: sf
    path "${norm_run}.chipseqspikeinfree.svg", optional: true, emit: plot
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \
        --bams=${bam_names} \
        --sample-ids=${sample_ids} \
        --norm-run=${norm_run} \
        --bin-size=${bin_size} \
        --cutoff=${cutoff} \
        --max-turns=${max_turns} \
        --chrom-sizes=${chrom_sizes} \
        --cpus=${task.cpus} \
        --out-prefix=${norm_run} \
        --sf-out=${norm_run}.chipseqspikeinfree.scale_factors.tsv \
        ${role == 'diagnostic' ? '--fail-soft=TRUE' : ''} \
        ${min_first_turn ? "--min-first-turn=${min_first_turn}" : ''} \
        ${params.norm_chipseqspikeinfree_args ?: ''}
    """

    stub:
    """
    printf 'sample_id\tmethod\traw_factor\tbigwig_multiplier\tedger_norm_factor\tqc_pass\tnote\n' > "${norm_run}.chipseqspikeinfree.scale_factors.tsv"
    for id in \$(echo "${sample_ids}" | tr ',' ' '); do
        printf '%s\tchipseqspikeinfree\t1\t1\t1\tTRUE\tstub\n' "\$id" >> "${norm_run}.chipseqspikeinfree.scale_factors.tsv"
    done
    echo '<svg xmlns="http://www.w3.org/2000/svg"/>' > "${norm_run}.chipseqspikeinfree.svg"
    """
}
