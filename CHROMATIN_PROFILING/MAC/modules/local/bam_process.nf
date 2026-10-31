include { runner } from './functions'

// MARKDUP + DEDUP_SPLIT + BLACKLIST_FILTER + BAM_FINALIZE fused into one task so that the markdup BAM, the
// target name-sorted BAM and the blacklist-filtered name-sorted BAM (each as big as the filtered BAM) never
// outlive the task. Each step script is run exactly as the former process ran it (same runner, same
// arguments). Each step's command line is written to <sample>_<step>_command.sh and run from there, so the
// per-step command logs keep their former published names. A transient file is deleted only after the step
// that consumes it has succeeded.
process BAM_PROCESS {
    tag "${sample_id}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.3.sif'   // picard + samtools + bedtools (the former MARKDUP, DEDUP_SPLIT, BLACKLIST_FILTER and BAM_FINALIZE images)
    publishDir "${params.outdir}/bam/spikein", mode: params.publish_mode, pattern: "*.spikein.bam*"
    publishDir "${params.outdir}/bam/stats", mode: params.publish_mode, pattern: "*.dedup_counts.tsv"
    publishDir "${params.outdir}/bam/stats", mode: params.publish_mode, pattern: "*.spikein.idxstats.txt"
    publishDir "${params.outdir}/bam/stats", mode: params.publish_mode, pattern: "*.flagstat.txt"
    publishDir "${params.outdir}/bam/stats", mode: params.publish_mode, pattern: "*.final.idxstats.txt"
    publishDir "${params.outdir}/bam/stats", mode: params.publish_mode, pattern: "*.orphan_check.txt"
    publishDir "${params.outdir}/bam/final", mode: params.publish_mode, pattern: "*.final.bam*"
    publishDir "${params.outdir}/qc/markdup", mode: params.publish_mode, pattern: "*.markdup_metrics.txt"
    publishDir "${params.outdir}/qc/markdup/log", mode: params.publish_mode, pattern: "*_markdup_command.sh"
    publishDir "${params.outdir}/bam/stats/log", mode: params.publish_mode, pattern: "*_dedup_split_command.sh"
    publishDir "${params.outdir}/bam/final/log", mode: params.publish_mode, pattern: "*_blacklist_filter_command.sh"
    publishDir "${params.outdir}/bam/final/log", mode: params.publish_mode, pattern: "*_bam_finalize_command.sh"

    input:
    tuple val(sample_id), path(bam), path(bai), val(is_control)
    val remove_duplicates
    val spikein_dedup
    val spike_prefix
    val mito_chrom
    path blacklist
    path markdup_script          // bin/markdup.sh, or the project's `markdup:` replacement
    path dedup_split_script      // bin/dedup_split.sh, or the project's `dedup_split:` replacement
    path blacklist_filter_script // bin/blacklist_filter.sh, or the project's `blacklist_filter:` replacement
    path bam_finalize_script     // bin/bam_finalize.sh, or the project's `bam_finalize:` replacement

    output:
    tuple val(sample_id), path("${sample_id}.final.bam"), path("${sample_id}.final.bam.bai"), emit: bam
    tuple val(sample_id), path("${sample_id}.spikein.bam"), path("${sample_id}.spikein.bam.bai"), emit: spikein
    path "${sample_id}.spikein.idxstats.txt", emit: spike_idxstats
    path "${sample_id}.dedup_counts.tsv", emit: counts
    path "${sample_id}.markdup_metrics.txt", emit: metrics
    path "${sample_id}.final.flagstat.txt", emit: flagstat
    path "${sample_id}.final.idxstats.txt", emit: idxstats
    path "${sample_id}.orphan_check.txt", emit: orphan_check
    path "${sample_id}_markdup_command.sh", emit: markdup_command_sh
    path "${sample_id}_dedup_split_command.sh", emit: dedup_split_command_sh
    path "${sample_id}_blacklist_filter_command.sh", emit: blacklist_filter_command_sh
    path "${sample_id}_bam_finalize_command.sh", emit: bam_finalize_command_sh

    script:
    """
    cat > "${sample_id}_markdup_command.sh" <<'CNR_STEP_EOF'
    #!/bin/bash -euo pipefail
    ${runner(markdup_script)} "${markdup_script}" \\
        --sample "${sample_id}" \\
        --bam "${bam}" \\
        --mem-gb ${task.memory.toGiga()} \\
        --threads ${task.cpus} \\
        ${params.markdup_args ?: ''}
    CNR_STEP_EOF
    bash -euo pipefail "${sample_id}_markdup_command.sh"

    cat > "${sample_id}_dedup_split_command.sh" <<'CNR_STEP_EOF'
    #!/bin/bash -euo pipefail
    ${runner(dedup_split_script)} "${dedup_split_script}" \\
        --sample "${sample_id}" \\
        --bam "${sample_id}.markdup.bam" \\
        --is-control ${is_control} \\
        --remove-duplicates ${remove_duplicates} \\
        --spikein-dedup ${spikein_dedup} \\
        --spike-prefix "${spike_prefix}" \\
        --mito-chrom "${mito_chrom}" \\
        --threads ${task.cpus} \\
        ${params.dedup_split_args ?: ''}
    CNR_STEP_EOF
    bash -euo pipefail "${sample_id}_dedup_split_command.sh"
    # the markdup BAM is no longer needed once the split succeeded
    rm -f "${sample_id}.markdup.bam" "${sample_id}.markdup.bam.bai"

    cat > "${sample_id}_blacklist_filter_command.sh" <<'CNR_STEP_EOF'
    #!/bin/bash -euo pipefail
    ${runner(blacklist_filter_script)} "${blacklist_filter_script}" \\
        --sample "${sample_id}" \\
        --bam "${sample_id}.target.namesorted.bam" \\
        --blacklist "${blacklist}" \\
        ${params.blacklist_filter_args ?: ''}
    CNR_STEP_EOF
    bash -euo pipefail "${sample_id}_blacklist_filter_command.sh"
    # the target name-sorted BAM is no longer needed once the blacklist filter succeeded
    rm -f "${sample_id}.target.namesorted.bam"

    cat > "${sample_id}_bam_finalize_command.sh" <<'CNR_STEP_EOF'
    #!/bin/bash -euo pipefail
    ${runner(bam_finalize_script)} "${bam_finalize_script}" \\
        --sample "${sample_id}" \\
        --bam "${sample_id}.noblacklist.namesorted.bam" \\
        --threads ${task.cpus} \\
        ${params.bam_finalize_args ?: ''}
    CNR_STEP_EOF
    bash -euo pipefail "${sample_id}_bam_finalize_command.sh"
    # the name-sorted BAM is no longer needed once the final BAM exists
    rm -f "${sample_id}.noblacklist.namesorted.bam"
    """

    stub:
    """
    touch "${sample_id}.final.bam" "${sample_id}.final.bam.bai"
    touch "${sample_id}.spikein.bam" "${sample_id}.spikein.bam.bai"
    printf 'ecoli_chr1\\t1000\\t0\\t0\\n*\\t0\\t0\\t0\\n' > "${sample_id}.spikein.idxstats.txt"
    printf 'sample_id\\tstage\\treads\\n${sample_id}\\ttarget_before_dedup\\t0\\n${sample_id}\\ttarget_after_dedup\\t0\\n${sample_id}\\tspikein_before_dedup\\t0\\n${sample_id}\\tspikein_after_dedup\\t0\\n' > "${sample_id}.dedup_counts.tsv"
    printf '## METRICS CLASS\\tpicard.sam.DuplicationMetrics\\nLIBRARY\\tUNPAIRED_READS_EXAMINED\\tREAD_PAIRS_EXAMINED\\n' > "${sample_id}.markdup_metrics.txt"
    echo "0 + 0 in total (QC-passed reads + QC-failed reads)" > "${sample_id}.final.flagstat.txt"
    printf 'stub_chr\\t1000\\t0\\t0\\n' > "${sample_id}.final.idxstats.txt"
    echo OK > "${sample_id}.orphan_check.txt"
    touch "${sample_id}_markdup_command.sh" "${sample_id}_dedup_split_command.sh"
    touch "${sample_id}_blacklist_filter_command.sh" "${sample_id}_bam_finalize_command.sh"
    """
}
