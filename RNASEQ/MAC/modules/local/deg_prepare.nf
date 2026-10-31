include { runner; extraArgs } from './functions'
include { shq } from './deg_functions'

process DEG_PREPARE {
    tag "prepare"
    container '/project/rephilllab/nf_containers/omics-r_1.0.sif'   // R 4 + DESeq2 etc., built from the omics Docker image R env lock
    cache 'deep'
    publishDir "${params.outdir}/deg/prepared", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/deg/prepared/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "prepare_command.sh" }

    input:
    path counts, stageAs: 'counts/*'            // merged counts table (Gene + one column per sample): COUNTS_MERGE output or params.deg_counts
    path metadata, stageAs: 'metadata/*'        // sample metadata: params.deg_metadata, or the samplesheet
    path comparisons, stageAs: 'comparisons/*'  // params.deg_comparisons (CSV)
    path common                                 // bin/deg_common.R, or the project's `deg_common` replacement
    val sample_col      // params.deg_sample_col: the metadata column holding the sample names
    val strip_suffix    // params.deg_strip_suffix: regex removed from the count column names ('' = none)
    val reference_levels  // params.deg_reference_levels: "factor=level,factor=level"
    path step_script    // bin/deg_prepare.R, or the project's `deg_prepare` script replacement

    output:
    path "counts_edited.rds",          emit: counts
    path "counts_edited.tsv",          emit: counts_tsv
    path "sample_metadata_edited.rds", emit: metadata
    path "sample_metadata_edited.csv", emit: metadata_csv
    path "comparisons.tsv",            emit: comparisons
    path "prepare_summary.txt",        emit: summary
    path "command.sh",                 emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \
        --common="${common}" \
        --counts="${counts}" \
        --metadata="${metadata}" \
        --comparisons="${comparisons}" \
        --sample-col=${shq(sample_col)} \
        --strip-suffix=${shq(strip_suffix)} \
        --reference-levels=${shq(reference_levels)} \
        --drop-excluded=TRUE \
        --out-dir=. \
        ${extraArgs('deg_prepare')}
    """

    stub:
    """
    echo stub > counts_edited.rds
    echo stub > counts_edited.tsv
    echo stub > sample_metadata_edited.rds
    echo stub > sample_metadata_edited.csv
    sed 's/,/\t/g' "${comparisons}" > comparisons.tsv
    echo "stub prepare summary" > prepare_summary.txt
    """
}
