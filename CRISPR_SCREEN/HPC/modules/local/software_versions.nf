/*
 * Records the exact version of each tool inside one container image
 * (bin/software_versions.sh: cutadapt, MAGeCK, FastQC, MultiQC, Python, R
 * and the R packages the pipeline's R scripts load — whichever the image
 * has). main.nf runs it once per distinct image the modules declare; METHODS
 * builds METHODS.md's software table from the results.
 *
 * Kept out of the steps themselves on purpose: probing versions inside
 * DEMULTIPLEX or MAGECK_COUNT would change their scripts and re-run them
 * (hours) in every project, for a few seconds of `--version` calls.
 *
 * executor 'local': runs inside the Nextflow head job (itself an LSF job on
 * a compute node, never the login node), so these few-second tasks don't
 * each wait in the LSF queue.
 */
include { runner } from './functions'

process SOFTWARE_VERSIONS {
    tag "${image}"
    container { image }
    executor 'local'
    publishDir "${params.metadata_dir}/pipeline_info/log", mode: params.publish_mode, pattern: "command.sh",
        saveAs: { "software_versions_${safeName(image)}_command.sh" }

    input:
    val  image             // a container image, e.g. quay.io/biocontainers/mageck:0.5.9.5--py39h2de1943_8
    path versions_script   // bin/software_versions.sh, or the project's `software_versions:` replacement

    output:
    path "${safeName(image)}.versions.tsv", emit: versions
    path "command.sh",                       emit: command_sh

    script:
    """
    ${runner(versions_script)} "${versions_script}" \\
        --image "${image}" \\
        --output "${safeName(image)}.versions.tsv" \\
        ${params.software_versions_args ?: ''}
    """

    stub:
    """
    printf 'image\\ttool\\tversion\\n' > "${safeName(image)}.versions.tsv"
    """
}

def safeName(image) {
    image.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')
}
