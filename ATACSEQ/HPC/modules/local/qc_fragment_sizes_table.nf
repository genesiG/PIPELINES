include { runner; extraArgs } from './functions'

process QC_FRAGMENT_SIZES_TABLE {
    tag "fragment_sizes"
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    // BAM inputs: no cache 'deep'
    publishDir "${params.outdir}/qc/fragment_sizes", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/fragment_sizes/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "fragment_sizes_table_command.sh" }

    input:
    path bams           // all S.final.bam (or S.tn5.bam) and their .bai, staged flat
    path samplesheet    // samplesheet.csv
    path step_script    // bin/qc_fragment_sizes_table.R, or the project's `qc_fragment_sizes_table:` replacement

    output:
    path "fragment_sizes.tsv",             emit: summary
    path "fragment_size_distribution.tsv", emit: distribution
    path "command.sh",                     emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --bam-dir=. \\
        --samplesheet="${samplesheet}" \\
        ${extraArgs('qc_fragment_sizes_table')}
    """

    stub:
    """
    printf 'sample_id\\tcondition\\treplicate\\ttmm_group\\tpeak_group\\tn_fragments\\tmean_size\\tsd_size\\tmedian_size\\tmode_size\\tmax_size\\tpct_lt_120\\tpct_120_300\\tpct_gt_300\\tn_over_max\\tfrac_nfr\\tfrac_mono\\tfrac_di\\n' > fragment_sizes.tsv
    printf 'stub1\\tcond1\\t1\\tgrp1\\tpg1\\t1000\\t180\\t60\\t170\\t160\\t900\\t10\\t70\\t20\\t0\\t0.1000\\t0.6000\\t0.1000\\n' >> fragment_sizes.tsv
    printf 'sample_id\\tsize\\tcount\\tfraction\\nstub1\\t160\\t500\\t0.5\\nstub1\\t200\\t500\\t0.5\\n' > fragment_size_distribution.tsv
    """
}
