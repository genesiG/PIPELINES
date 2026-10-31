include { runner; extraArgs } from './functions'

// Output prefix of one ALIGN task: S.chunkNNN (NNN = chunk index, 3 digits) when chunked, else S.
def alignPrefix(sample_id, chunk) {
    return chunk >= 0 ? "${sample_id}.chunk${chunk.toString().padLeft(3, '0')}".toString() : "${sample_id}".toString()
}

process ALIGN {
    tag "${alignPrefix(sample_id, chunk)}"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    publishDir "${params.intermediate_dir}/aligned", mode: params.publish_mode, pattern: "*.aligned.bam", saveAs: { fn -> chunk < 0 ? fn : null }, enabled: params.publish_intermediates.toString() == 'true'
    publishDir "${params.outdir}/qc/bowtie2", mode: params.publish_mode, pattern: "*.bowtie2.log", saveAs: { fn -> chunk < 0 ? fn : null }
    publishDir "${params.outdir}/qc/bowtie2/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${sample_id}${chunk >= 0 ? '.chunk' + chunk.toString().padLeft(3, '0') : ''}_align_command.sh" }

    input:
    tuple val(sample_id), val(chunk), val(skip), val(upto), path(r1), path(r2)   // chunk -1 = not chunked (skip/upto unused)
    path index_files   // all *.bt2 / *.bt2l files of the index prefix, staged flat
    val index_base     // bare file-name prefix of the index
    val mode           // local | end_to_end
    val minins         // null = bowtie2 default (no -I)
    val maxins
    val bowtie2_args   // params.bowtie2_args, passed to bowtie2 after `--`
    path step_script   // bin/align.sh, or the project's `align:` replacement

    output:
    tuple val(sample_id), val(chunk), path("${alignPrefix(sample_id, chunk)}.aligned.bam"), path("${alignPrefix(sample_id, chunk)}.bowtie2.log"), emit: aligned
    path "command.sh", emit: command_sh

    script:
    def prefix = alignPrefix(sample_id, chunk)
    def minins_arg = (minins == null || minins.toString().trim() in ['', 'null']) ? '' : minins.toString().trim()
    def chunk_args = chunk >= 0 ? "--chunk-skip ${skip} --chunk-upto ${upto}" : ''
    """
    ${runner(step_script)} "${step_script}" \\
        --sample "${prefix}" \\
        --r1 "${r1}" \\
        --r2 "${r2}" \\
        --index-base "${index_base}" \\
        --mode "${mode}" \\
        --minins "${minins_arg}" \\
        --maxins ${maxins} \\
        --threads ${task.cpus} \\
        ${chunk_args} \\
        -- ${bowtie2_args ?: ''} ${extraArgs('align')}
    """

    stub:
    def prefix = alignPrefix(sample_id, chunk)
    """
    touch "${prefix}.aligned.bam"
    cat > "${prefix}.bowtie2.log" <<'EOF'
10000 reads; of these:
  10000 (100.00%) were paired; of these:
    1000 (10.00%) aligned concordantly 0 times
    8000 (80.00%) aligned concordantly exactly 1 time
    1000 (10.00%) aligned concordantly >1 times
90.00% overall alignment rate
EOF
    """
}
