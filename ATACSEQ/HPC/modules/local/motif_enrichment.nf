/*
 * Binned motif enrichment (monaLisa, binomial test, composition-controlled) for one analysis (bin/motif_enrichment.R).
 * mode diff:    name = comparison, region = region key (promoters | peaks | all_sites | ...), table = differential table or annotation table,
 *               target_label = region label.
 * mode gene:    name = motif_targets name, table = all_sites table of the universe comparison, target = gene symbol.
 * mode regions: name = motif_targets name, table = all_sites table of the universe comparison, target_bed = BED of the target regions.
 * target_bed is in the tuple (per analysis): assets/NO_FILE unless mode is regions. pwm_rds is '' (JASPAR2024 from the package) or a
 * cached PWMList path.
 * Published to motif/NAME/REGION/ (diff), motif/NAME/ (gene, regions). A skipped analysis writes NAME.skipped.txt only, so it emits
 * nothing on the table and summary tuples (optional: a missing file emits nothing; MOTIF_COMMON groups the emitted ones by region).
 */
include { runner; extraArgs } from './functions'

process MOTIF_ENRICHMENT {
    tag "${name}"
    container '/project/rephilllab/nf_containers/omics-r_1.0.sif'
    cpus params.motif_cpus
    publishDir path: { region ? "${params.outdir}/motif/${name}/${region}" : "${params.outdir}/motif/${name}" }, mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir path: { (region ? "${params.outdir}/motif/${name}/${region}" : "${params.outdir}/motif/${name}") + "/log" }, mode: params.publish_mode, pattern: "command.sh", saveAs: { "${name}_motif_enrichment_command.sh" }

    input:
    tuple val(name), val(region), val(mode), path(table), val(target), val(target_label), val(min_peaks), path(target_bed)
                                    // target_label: diff = region label (Promoters); gene/regions = bin label of the target
                                    // target_bed: regions mode: BED of the target regions; otherwise assets/NO_FILE
    val species                     // species
    val genome_build                // genome_build
    val bsgenome                    // bsgenome: BSgenome package (BSgenome.Hsapiens.UCSC.hs1)
    val pwm_rds                     // motif_pwm_rds: cached PWMList path, or '' for JASPAR2024 from the package
    val min_score                   // motif_min_score
    val test                        // motif_test
    val max_kmer                    // motif_max_kmer
    val background                  // motif_background
    val seed                        // motif_seed
    val top                         // motif_top
    val padj_primary                // motif_padj_primary
    val padj_relaxed                // motif_padj_relaxed
    val lfc                         // motif_lfc (diff mode)
    val gene_column                 // motif_gene_column (gene mode)
    val match                       // motif_match: exact | substring (gene mode)
    val plot_font                   // plot_font: Helvetica (font family of every plot)
    val plot_format                 // plot_format: svg | pdf
    path step_script                // bin/motif_enrichment.R, or the project's `motif_enrichment:` replacement

    output:
    tuple val(name), val(region), path("${name}_motif_analysis.csv"), optional: true, emit: table
    path "${name}_motif_heatmap_padj.{svg,pdf}",      optional: true, emit: heatmap
    path "${name}_motif_similarity_padj.{svg,pdf}",   optional: true, emit: similarity
    path "${name}.gc_fraction.{svg,pdf}",             optional: true, emit: gc
    path "${name}.dinucleotide_frequency.{svg,pdf}",  optional: true, emit: dinuc
    tuple val(name), val(region), path("${name}.motif_summary.csv"), optional: true, emit: summary
    path "${name}.skipped.txt",                       optional: true, emit: skipped
    path "command.sh",                                emit: command_sh

    script:
    def mode_args = mode == 'diff' ?
        "--region=${region} --region-label=\"${target_label}\" --lfc=${lfc} --min-bin-peaks=${min_peaks}" :
        "--target=${mode == 'regions' ? target_bed : target} --target-label=\"${target_label}\" --min-peaks=${min_peaks}" +
        (mode == 'gene' ? " --gene-column=${gene_column} --match=${match}" : '')
    """
    ${runner(step_script)} "${step_script}" \\
        --mode=${mode} \\
        --name=${name} \\
        --table=${table} \\
        --species="${species}" \\
        --genome-build=${genome_build} \\
        --bsgenome=${bsgenome} \\
        --pwm-rds=${pwm_rds ?: ''} \\
        --min-score=${min_score} \\
        --test=${test} \\
        --max-kmer=${max_kmer} \\
        --background=${background} \\
        --seed=${seed} \\
        --top=${top} \\
        --padj-primary=${padj_primary} \\
        --padj-relaxed=${padj_relaxed} \\
        --font=${plot_font} \\
        --plot-format=${plot_format} \\
        --cpus=${task.cpus} \\
        ${mode_args} \\
        ${extraArgs('motif_enrichment')}
    """

    stub:
    """
    printf 'motif.id,motif.name,log2enr_down,log2enr_up,padj_down,padj_up\\n' > "${name}_motif_analysis.csv"
    printf 'MA0001.1,ARNT,0.1,-0.2,0.5,0.4\\n' >> "${name}_motif_analysis.csv"
    for f in "${name}_motif_heatmap_padj" "${name}_motif_similarity_padj" "${name}.gc_fraction" "${name}.dinucleotide_frequency"; do
        printf '<svg xmlns="http://www.w3.org/2000/svg" width="288pt" height="288pt"><text font-family="Helvetica">%s</text></svg>\\n' "\$f" > "\$f.svg"
    done
    printf 'comparison,region,total_peaks,n_peaks_up,n_peaks_down,n_motifs_tested,n_sig_padj_0.01,n_sig_padj_0.1,n_motifs_plotted\\n' > "${name}.motif_summary.csv"
    printf '${name},${region},100,50,50,1,0,0,1\\n' >> "${name}.motif_summary.csv"
    """
}
