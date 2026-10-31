/*
 * Differential accessibility (csaw + edgeR) for ONE comparison over peaks and/or promoters (bin/diff_accessibility.R).
 * Peaks and promoters of one comparison run in one task so that the 10 kb bins are counted once.
 * ctrl_peaks / case_peaks are staged under ctrl_peaks/ and case_peaks/; txdb is assets/NO_FILE when promoters are not requested.
 * regions and fixed_dispersion are per comparison (diff_regions minus the tables given in diff_tables; the fixed_dispersion column
 * of comparisons.csv, else diff_fixed_dispersion), so they travel in the tuple. The promoters table is emitted only when the task
 * computed promoters (optional on the tuple: a missing file emits nothing).
 */
include { runner; extraArgs } from './functions'

process DIFF_ACCESS {
    tag "${name}"
    container '/project/rephilllab/nf_containers/omics-r_1.0.sif'
    cpus params.diff_cpus
    publishDir "${params.outdir}/differential", mode: params.publish_mode, pattern: "*.{csv,tsv,pdf,svg}", saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/differential/bed", mode: params.publish_mode, pattern: "*.bed"
    publishDir "${params.outdir}/differential/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${name}_diff_access_command.sh" }

    input:
    tuple val(name), val(sample_ids), path(bams), path(bais), val(groups), val(blocks), path(ctrl_peaks, stageAs: 'ctrl_peaks/*'), path(case_peaks, stageAs: 'case_peaks/*'), val(regions), val(fixed_dispersion)
                            // regions: peaks | peaks,promoters (this comparison's region sets still to compute)
                            // fixed_dispersion: this comparison's dispersion for the LRT branch (comparisons.csv, else diff_fixed_dispersion)
    path txdb               // TxDb (TXDB_BUILD), or assets/NO_FILE when promoters are not requested
    val genome_build        // genome_build
    val style               // UCSC
    val minq                // diff_minq
    val pe                  // diff_pe: both | none
    val max_frag            // diff_max_frag
    val dedup               // diff_dedup as resolved by main.nf: true | false | auto
    val bin_size            // diff_bin_size
    val lfc                 // diff_lfc
    val fdr                 // diff_fdr
    val table               // diff_table: best | combined
    val qlf_min_libs        // diff_qlf_min_libs
    val promoter_up         // diff_promoter_up
    val promoter_down       // diff_promoter_down
    val promoter_min_fc     // diff_promoter_min_fc
    val plot_font           // plot_font: Helvetica (font family of every plot)
    val plot_format         // plot_format: svg | pdf
    path step_script        // bin/diff_accessibility.R, or the project's `diff_accessibility:` replacement

    output:
    tuple val(name), path("${name}.csv"),                            emit: peaks
    tuple val(name), path("${name}_promoters.csv"), optional: true, emit: promoters
    path "${name}.{norm_factors.tsv,summary.tsv,replicates.csv}", emit: tables
    path "${name}.{diag_*.svg,diagnostics.pdf}", optional: true,     emit: plots
    path "${name}*.bed", optional: true,                             emit: beds
    path "command.sh",                                               emit: command_sh

    script:
    def bam_names  = [bams].flatten().collect { it.name }.join(',')
    def ctrl_names = [ctrl_peaks].flatten().collect { "ctrl_peaks/${it.name.tokenize('/')[-1]}" }.join(',')
    def case_names = [case_peaks].flatten().collect { "case_peaks/${it.name.tokenize('/')[-1]}" }.join(',')
    """
    ${runner(step_script)} "${step_script}" \\
        --name=${name} \\
        --bams=${bam_names} \\
        --sample-ids=${sample_ids} \\
        --groups=${groups} \\
        --blocks=${blocks} \\
        --control-peaks=${ctrl_names} \\
        --case-peaks=${case_names} \\
        --regions=${regions} \\
        --txdb=${txdb} \\
        --genome-build=${genome_build} \\
        --style=${style} \\
        --minq=${minq} \\
        --pe=${pe} \\
        --max-frag=${max_frag} \\
        --dedup=${dedup} \\
        --bin-size=${bin_size} \\
        --lfc=${lfc} \\
        --fdr=${fdr} \\
        --table=${table} \\
        --qlf-min-libs=${qlf_min_libs} \\
        --fixed-dispersion=${fixed_dispersion} \\
        --promoter-up=${promoter_up} \\
        --promoter-down=${promoter_down} \\
        --promoter-min-fc=${promoter_min_fc} \\
        --font=${plot_font} \\
        --plot-format=${plot_format} \\
        --cpus=${task.cpus} \\
        ${extraArgs('diff_accessibility')}
    """

    stub:
    """
    printf '"",seqnames,start,end,width,strand,PValue,FDR,logFC,score,name\\n' > "${name}.csv"
    printf '"1",chr1,1000,1500,500,*,0.001,0.01,1.2,1.2,region1\\n' >> "${name}.csv"
    printf '"2",chr2,2000,2400,400,*,0.02,0.1,-0.9,-0.9,region2\\n' >> "${name}.csv"
    printf 'chr1\\t1000\\t1500\\tregion1\\t1.2\\t*\\n' > "${name}.bed"
    printf 'chr2\\t2000\\t2400\\tregion2\\t-0.9\\t*\\n' > "${name}_downregulated.bed"
    printf 'chr1\\t1000\\t1500\\tregion1\\t1.2\\t*\\n' > "${name}_upregulated.bed"
    if [[ "${regions}" == *promoters* ]]; then
        printf '"",seqnames,start,end,width,strand,gene_id,PValue,FDR,logFC,score,name\\n' > "${name}_promoters.csv"
        printf '"1",chr1,1000,3200,2200,+,GENE1,0.002,0.02,1.5,1.5,GENE1\\n' >> "${name}_promoters.csv"
        printf '"2",chr2,5000,7200,2200,-,GENE2,0.03,0.2,-0.7,-0.7,GENE2\\n' >> "${name}_promoters.csv"
        printf 'chr1\\t1000\\t3200\\tGENE1\\t1.5\\t+\\n' > "${name}_promoters.bed"
        printf 'chr1\\t1000\\t3200\\tGENE1\\t1.5\\t+\\n' > "${name}_promoters_up.bed"
        printf 'chr2\\t5000\\t7200\\tGENE2\\t-0.7\\t-\\n' > "${name}_promoters_down.bed"
    fi
    printf 'sample_id\\tnorm_factor\\tlib_size\\tgroup\\tblock\\n' > "${name}.norm_factors.tsv"
    for id in \$(echo "${sample_ids}" | tr ',' ' '); do
        printf '%s\\t1\\t1000000\\tNA\\tNA\\n' "\$id" >> "${name}.norm_factors.tsv"
    done
    printf 'region_set\\tmetric\\tvalue\\n' > "${name}.summary.tsv"
    printf 'peaks\\tn_regions\\t2\\n' >> "${name}.summary.tsv"
    printf 'peaks\\tn_increased\\t1\\n' >> "${name}.summary.tsv"
    printf 'peaks\\tn_decreased\\t1\\n' >> "${name}.summary.tsv"
    printf 'Comparison,Case,Control,Case_Replicates,Control_Replicates,Total_Samples\\n' > "${name}.replicates.csv"
    printf '${name},Condition,Control,1,1,2\\n' >> "${name}.replicates.csv"
    for d in diag_ma diag_bcv diag_mds diag_qldisp; do
        printf '<svg xmlns="http://www.w3.org/2000/svg" width="288pt" height="144pt"><text font-family="Helvetica">%s</text></svg>\\n' "\$d" > "${name}.\$d.svg"
    done
    """
}
