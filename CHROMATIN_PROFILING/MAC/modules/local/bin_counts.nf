/*
 * One csaw bin-count pass over all active libraries (bin/bin_counts.R), one task per distinct bin size. Feeds TMM, NORM_QC,
 * NORM_QC_PLOT, REPLICATE_CORRELATION and REPLICATE_CORRELATION_PLOT (each subsets its own columns).
 * bam_names / sample_ids: comma lists in the same order as the staged BAMs; idxstats: the final idxstats of the same
 * libraries in sample_ids order. max.frag follows maxFrag() (csaw_max_frag, bowtie2_maxins, 500).
 * No cache 'deep': BAM inputs.
 */
include { runner; maxFrag } from './functions'

process BIN_COUNTS {
    tag "${bin_size}"
    container '/project/rephilllab/nf_containers/cutandrun-r_1.3.sif'
    publishDir "${params.intermediate_dir}/bin_counts", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }, enabled: params.publish_intermediates.toString() == 'true'
    publishDir "${params.intermediate_dir}/bin_counts/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "bins_${bin_size}_command.sh" }, enabled: params.publish_intermediates.toString() == 'true'

    input:
    tuple val(bin_size), val(sample_ids), val(bam_names), path(bams), path(bais)
    path idxstats        // S.final.idxstats.txt of every library, sample_ids order
    path blacklist       // [] when the project has none
    path step_script     // bin/bin_counts.R

    output:
    tuple val(bin_size), path("bins_${bin_size}.counts.tsv.gz"), path("bins_${bin_size}.libsizes.tsv"), emit: counts
    path "command.sh", emit: command_sh

    script:
    def idx = (idxstats instanceof Collection ? idxstats : [idxstats]).collect { f -> f.name }.join(',')
    """
    ${runner(step_script)} "${step_script}" \
        --bams=${bam_names} \
        --sample-ids=${sample_ids} \
        --idxstats=${idx} \
        --bin-size=${bin_size} \
        --mapq=${params.mapq} \
        --dedup=${params.remove_duplicates.toString().toUpperCase()} \
        --max-frag=${maxFrag()} \
        ${blacklist ? "--blacklist=${blacklist}" : ''} \
        --cpus=${task.cpus} \
        --out-prefix=bins_${bin_size} \
        ${params.bin_counts_args ?: ''}
    """

    stub:
    """
    ids="${sample_ids}"
    { printf 'chr\tstart\tend\tblacklisted'; for id in \$(echo "\$ids" | tr ',' ' '); do printf '\t%s' "\$id"; done; printf '\n'
      for r in 1 2 3; do
          printf 'chr1\t%s\t%s\tFALSE' "\$(( (r - 1) * ${bin_size} + 1 ))" "\$(( r * ${bin_size} ))"
          for id in \$(echo "\$ids" | tr ',' ' '); do printf '\t%s' "\$(( 10 * r ))"; done; printf '\n'
      done; } | gzip > "bins_${bin_size}.counts.tsv.gz"
    printf 'sample_id\tcsaw_total\tbam_mapped_reads\tbam_fragments\tbin_size\tmapq\tmax_frag\tdedup\n' > "bins_${bin_size}.libsizes.tsv"
    for id in \$(echo "\$ids" | tr ',' ' '); do
        printf '%s\t60\t120\t60\t${bin_size}\t${params.mapq}\t${maxFrag()}\t${params.remove_duplicates.toString().toUpperCase()}\n' "\$id" >> "bins_${bin_size}.libsizes.tsv"
    done
    """
}
