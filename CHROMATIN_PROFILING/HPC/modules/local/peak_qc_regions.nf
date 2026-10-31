/*
 * Automatic peak-QC regions of one run (bin/peak_qc_regions.py). Chooses the regions from the peak sets and writes
 * regions.tsv (NN_label<TAB>chr:start-end), regions.detail.tsv, tracks.loci.txt, params.json and peak_qc_regions.log.
 * Published to results/tracks/peak_qc/ (peak_qc_regions.log included); the command script is saved as tracks/peak_qc/log/peak_qc_regions_command.sh.
 * The TSV inputs (peaksets_tsv, classes_tsv, bigwigs_tsv) hold the STAGED base names of the files listed in peak_files
 * (every bed, classification table and bigWig is staged in the work directory under its own base name, so the
 * columns reference these names). classes_tsv is [] when no classification is used.
 * Optional params (see nextflow.config, peak_qc_*): loci, negative, max_regions, exclude_bed, shape_marks, nopeak_igg_max_fold;
 * spike_prefix and mito_chrom are passed when set.
 */
include { runner } from './functions'

process PEAK_QC_REGIONS {
    tag "peak_qc"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    publishDir { "${params.outdir}/tracks/peak_qc" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir { "${params.outdir}/tracks/peak_qc/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "peak_qc_regions_command.sh" }

    input:
    path peaksets_tsv
    path classes_tsv    // [] when no classification is used
    path bigwigs_tsv
    path peak_files   // all peak BEDs, classification tables and bigWigs, staged by base name
    path gtf
    path chrom_sizes
    path blacklist    // [] when the run has none
    path exclude_beds    // [] unless peak_qc_exclude_bed is set
    val exclude_names    // comma list of the exclude_beds file names, '' when none
    val shape_marks    // comma list or '' (script default)
    val nopeak_igg_max_fold    // number or '': --nopeak-igg-max-fold; '' when off or when no IgG bigWig is staged
    path step_script   // bin/peak_qc_regions.py

    output:
    path "regions.tsv", emit: regions
    path "regions.detail.tsv", emit: detail
    path "tracks.loci.txt", emit: loci
    path "params.json", emit: params_json
    path "peak_qc_regions.log", emit: log
    path "command.sh", emit: command_sh

    script:
    def loci = (params.peak_qc_loci ?: '').toString().trim()
    def negative = (params.peak_qc_negative ?: 'auto').toString().trim()
    def max_regions = params.peak_qc_max_regions ?: 15
    def spike = (params.spike_prefix ?: '').toString().trim()
    def mito = (params.mito_chrom ?: '').toString().trim()
    """
    ${runner(step_script)} "${step_script}" \\
        --peaksets "${peaksets_tsv}" \\
        ${classes_tsv ? "--classes '${classes_tsv}'" : ''} \\
        --bigwigs "${bigwigs_tsv}" \\
        --gtf "${gtf}" \\
        --chrom-sizes "${chrom_sizes}" \\
        ${blacklist ? "--blacklist '${blacklist}'" : ''} \\
        ${exclude_names ? "--exclude-bed '${exclude_names}'" : ''} \\
        ${loci ? "--loci '${loci}'" : ''} \\
        --negative '${negative}' \\
        ${shape_marks ? "--shape-marks '${shape_marks}'" : ''} \\
        ${nopeak_igg_max_fold ? "--nopeak-igg-max-fold '${nopeak_igg_max_fold}'" : ''} \\
        ${spike ? "--spike-prefix '${spike}'" : ''} \\
        ${mito ? "--mito '${mito}'" : ''} \\
        --max-regions=${max_regions} \\
        --outdir . \\
        ${params.peak_qc_regions_args ?: ''}
    """

    stub:
    """
    printf '01_stub_a\\tchr1:1000-200000\\n02_stub_b\\tchr1:300000-500000\\n' > regions.tsv
    printf 'name\\tcategory\\tsource\\tanchor\\tanchor_name\\tkey_values\\trank_tried\\tn_rejected\\tstatus\\n' > regions.detail.tsv
    printf 'chr1:1000-200000,chr1:300000-500000\\n' > tracks.loci.txt
    printf '{}\\n' > params.json
    printf 'stub\\n' > peak_qc_regions.log
    printf '#!/bin/bash\\n' > command.sh
    """
}
