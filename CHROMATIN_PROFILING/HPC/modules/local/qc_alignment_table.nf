include { runner } from './functions'

process QC_ALIGNMENT_TABLE {
    tag "alignment_stats"
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    cache 'deep'
    publishDir "${params.outdir}/qc/alignment", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/alignment/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "alignment_table_command.sh" }

    input:
    path stats          // all S.raw.flagstat/idxstats, S.final.flagstat/idxstats and S.markdup_metrics.txt, staged flat (raw files absent for pre-aligned BAM samples)
    path samplesheet    // samplesheet.csv
    val spike_prefix    // params.spike_prefix
    val mito_chrom      // params.mito_chrom (not in the contract's input list; added so chrM is not hard-coded)
    path step_script    // bin/qc_alignment_table.R, or the project's `qc_alignment_table:` replacement

    output:
    path "alignment_stats.tsv", emit: table
    path "command.sh",          emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --stats-dir=. \\
        --samplesheet="${samplesheet}" \\
        --spike-prefix="${spike_prefix}" \\
        --mito-chrom="${mito_chrom}" \\
        ${params.qc_alignment_table_args ?: ''}
    """

    stub:
    """
    printf 'sample_id\\tantibody\\tcondition\\treplicate\\tgroup\\traw_total_reads\\traw_mapped_reads\\tpct_mapped\\tn_chrM\\tpct_chrM\\tn_spike_in\\tpct_spike_in\\tn_target_mapped\\tpct_unmapped\\tfinal_total_reads\\tfinal_mapped_reads\\tfinal_n_chrM\\tfinal_n_spike_in\\tpct_final_of_raw_mapped\\tdup_pairs_examined\\tdup_pairs_duplicate\\tdup_pairs_optical\\tpct_duplication\\n' > alignment_stats.tsv
    printf 'stub1\\tAb1\\tcond1\\t1\\tgrp1\\t1000\\t900\\t90\\t10\\t1.1\\t20\\t2.2\\t870\\t10\\t700\\t700\\t0\\t0\\t77.8\\t400\\t40\\t1\\t10\\n' >> alignment_stats.tsv
    """
}
