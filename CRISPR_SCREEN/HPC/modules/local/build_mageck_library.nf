/*
 * Converts a library pool CSV (gene,gRNA_ID,sequence) into the tab-separated
 * <sgRNA_id  sequence  gene> format `mageck count --list-seq` expects.
 * Ported from step_2_count.py's _build_mageck_library().
 */
include { runner; extraArgs } from './functions'

process BUILD_MAGECK_LIBRARY {
    tag "${lib_name}"
    container 'quay.io/biocontainers/cutadapt:5.2--py39hbcbf7aa_0'   // any py3 image; reused to avoid an extra pull
    // command.sh excluded here, published separately below -- see fastqc.nf
    // for why (no pattern here would otherwise also match it).
    publishDir "${params.intermediate_dir}/mageck_libraries", mode: params.publish_mode, pattern: "*_mageck_library.txt"
    publishDir "${params.intermediate_dir}/mageck_libraries/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "${lib_name}_command.sh" }

    input:
    tuple val(lib_name), path(library_csv)   // lib_name = file name without .csv; samplesheet `library` values refer to it
    path helper_script   // bin/build_mageck_library.py, or the project's `build_mageck_library:` replacement

    output:
    tuple val(lib_name), path("${lib_name}_mageck_library.txt"), emit: library
    path "command.sh", emit: command_sh

    script:
    """
    ${runner(helper_script)} "${helper_script}" \\
        --input "${library_csv}" \\
        --output "${lib_name}_mageck_library.txt" \\
        ${extraArgs('build_mageck_library')}
    """

    stub:
    """
    printf 'sgRNA_id\\tsequence\\tgene\\n' > "${lib_name}_mageck_library.txt"
    """
}
