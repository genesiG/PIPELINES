/*
 * ChIPseeker annotation of one differential table: all / promoter / non-promoter / distal regions (bin/diff_annotate.R).
 * txdb: TXDB_BUILD output (the same GTF as the promoter table). bsgenome and org_db are package names from params.
 */
include { runner; extraArgs } from './functions'

process DIFF_ANNOTATE {
    tag "${name}"
    container '/project/rephilllab/nf_containers/omics-r_1.0.sif'
    publishDir "${params.outdir}/annotation", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/annotation/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${name}_diff_annotate_command.sh" }

    input:
    tuple val(name), path(table)   // differential table: NAME.csv or NAME_promoters.csv (DIFF_ACCESS)
    path txdb                      // TxDb (TXDB_BUILD)
    val species                    // species
    val genome_build               // genome_build (hs1)
    val bsgenome                   // bsgenome: BSgenome package (BSgenome.Hsapiens.UCSC.hs1)
    val org_db                     // org_db: org.Hs.eg.db
    val tss_region                 // annotation_tss_region: -2000,200
    val fdr                        // annotation_fdr: blank = all regions
    val plot_font                  // plot_font: Helvetica (font family of the plot)
    val plot_format                // plot_format: svg | pdf
    path step_script               // bin/diff_annotate.R, or the project's `diff_annotate:` replacement

    output:
    tuple val(name), path("${name}_*_sites.csv"), emit: tables
    path "${name}_*_sites.bed",                   emit: beds
    path "${name}_peak_distribution.{svg,pdf}", optional: true, emit: pdf
    path "${name}.annotation_summary.tsv",        emit: summary
    path "command.sh",                            emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --name=${name} \\
        --table=${table} \\
        --txdb=${txdb} \\
        --species="${species}" \\
        --genome-build=${genome_build} \\
        --bsgenome=${bsgenome} \\
        --org-db=${org_db} \\
        --tss-region=${tss_region} \\
        --fdr=${fdr ?: ''} \\
        --font=${plot_font} \\
        --plot-format=${plot_format} \\
        --cpus=${task.cpus} \\
        ${extraArgs('diff_annotate')}
    """

    stub:
    """
    printf '"",seqnames,start,end,width,strand,PValue,FDR,logFC,annotation,geneId\\n' > "${name}_all_sites.csv"
    printf '"1",chr1,1000,1500,500,*,0.001,0.01,1.2,Promoter (<=1kb),GENE1\\n' >> "${name}_all_sites.csv"
    printf '"2",chr2,2000,2400,400,*,0.02,0.1,-0.9,Distal Intergenic,GENE2\\n' >> "${name}_all_sites.csv"
    printf '"",seqnames,start,end,width,strand,PValue,FDR,logFC,annotation,geneId\\n' > "${name}_promoter_sites.csv"
    printf '"1",chr1,1000,1500,500,*,0.001,0.01,1.2,Promoter (<=1kb),GENE1\\n' >> "${name}_promoter_sites.csv"
    printf '"",seqnames,start,end,width,strand,PValue,FDR,logFC,annotation,geneId\\n' > "${name}_non-promoter_sites.csv"
    printf '"2",chr2,2000,2400,400,*,0.02,0.1,-0.9,Distal Intergenic,GENE2\\n' >> "${name}_non-promoter_sites.csv"
    printf '"",seqnames,start,end,width,strand,PValue,FDR,logFC,annotation,geneId\\n' > "${name}_distal_sites.csv"
    printf '"2",chr2,2000,2400,400,*,0.02,0.1,-0.9,Distal Intergenic,GENE2\\n' >> "${name}_distal_sites.csv"
    printf 'chr1\\t1000\\t1500\\tregion1\\t1.2\\t*\\n' > "${name}_all_sites.bed"
    printf 'chr1\\t1000\\t1500\\tregion1\\t1.2\\t*\\n' > "${name}_promoter_sites.bed"
    printf 'chr2\\t2000\\t2400\\tregion2\\t-0.9\\t*\\n' > "${name}_non-promoter_sites.bed"
    printf 'chr2\\t2000\\t2400\\tregion2\\t-0.9\\t*\\n' > "${name}_distal_sites.bed"
    printf '<svg xmlns="http://www.w3.org/2000/svg" width="360pt" height="360pt"><text font-family="Helvetica">%s</text></svg>\\n' "${name}" > "${name}_peak_distribution.svg"
    printf 'set\\tn_regions\\nall_sites\\t2\\npromoter_sites\\t1\\nnon-promoter_sites\\t1\\ndistal_sites\\t1\\n' > "${name}.annotation_summary.tsv"
    """
}
