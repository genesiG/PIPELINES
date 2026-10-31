include { runner; extraArgs } from './functions'

process QC_CHROM_FRACTION {
    tag "chrom_fraction"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    cache 'deep'
    publishDir "${params.outdir}/qc/chrom_fraction", mode: params.publish_mode, saveAs: { fn -> fn == 'command.sh' ? null : fn }
    publishDir "${params.outdir}/qc/chrom_fraction/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "chrom_fraction_command.sh" }

    input:
    path '*'           // all S.raw.idxstats.txt (BAM_FILTER), collected and staged flat
    val chroms         // params.remove_chroms (comma list)
    path step_script   // bin/qc_chrom_fraction.py, or the project's `qc_chrom_fraction:` replacement

    output:
    path "chrom_fraction.tsv", emit: table
    path "command.sh",         emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --input-dir . \\
        --chroms=${chroms ?: ''} \\
        --out chrom_fraction.tsv \\
        ${extraArgs('qc_chrom_fraction')}
    """

    stub:
    def names = (chroms ?: '').tokenize(',').collect { it.trim() }.findAll { it }
    def header = (['sample_id', 'raw_mapped_reads'] + names.collectMany { [it + '_reads', it + '_pct'] }).join('\t')
    def row = (['stub1', '1000'] + names.collectMany { ['0', '0.0000'] }).join('\t')
    """
    printf '%s\\n' '${header}' > chrom_fraction.tsv
    printf '%s\\n' '${row}' >> chrom_fraction.tsv
    """
}
