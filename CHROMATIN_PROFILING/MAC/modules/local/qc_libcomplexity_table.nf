include { runner } from './functions'

process QC_LIBCOMPLEXITY_TABLE {
    tag "library_complexity"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    cache 'deep'
    publishDir "${params.outdir}/qc/library_complexity", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/library_complexity/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "library_complexity_table_command.sh" }

    input:
    path results        // all S.complexity.tsv, S.elc.txt and S.markdup_metrics.txt, staged flat
    path samplesheet    // samplesheet.csv
    path step_script    // bin/qc_library_complexity_table.py, or the project's `qc_library_complexity_table:` replacement

    output:
    path "library_complexity.tsv", emit: table
    path "command.sh",             emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --input-dir . \\
        --samplesheet "${samplesheet}" \\
        --out library_complexity.tsv \\
        ${params.qc_library_complexity_table_args ?: ''}
    """

    stub:
    """
    printf 'sample_id\\tantibody\\tcondition\\treplicate\\tgroup\\ttotal_fragments\\tdistinct_fragments\\tonce\\ttwice\\tNRF\\tPBC1\\tPBC2\\telc_read_pairs_examined\\telc_pct_duplication\\telc_estimated_library_size\\tmarkdup_read_pairs_examined\\tmarkdup_pct_duplication\\tmarkdup_estimated_library_size\\n' > library_complexity.tsv
    printf 'stub1\\tAb1\\tcond1\\t1\\tgrp1\\t1000\\t900\\t850\\t40\\t0.9\\t0.944\\t21.25\\t1000\\t10\\t10000\\t1000\\t10\\t10000\\n' >> library_complexity.tsv
    """
}
