/*
 * Records the exact version of each tool inside the pipeline's Docker image
 * (bin/software_versions.sh: cutadapt, MAGeCK, FastQC, MultiQC, Python, R
 * and the R packages the pipeline's R scripts load). Docker version: every
 * step runs in the one image the Nextflow head process itself runs in
 * (docker/Dockerfile), so main.nf passes exactly one image name here and
 * this runs once. METHODS builds METHODS.md's software table from the results.
 *
 * Kept out of the steps themselves on purpose: probing versions inside
 * DEMULTIPLEX or MAGECK_COUNT would change their scripts and re-run them
 * in every project, for a few seconds of `--version` calls.
 */
include { runner; extraArgs } from './functions'

process SOFTWARE_VERSIONS {
    tag "${image}"
    executor 'local'
    publishDir "${params.metadata_dir}/pipeline_info/log", mode: params.publish_mode, pattern: "command.sh",
        saveAs: { "software_versions_${safeName(image)}_command.sh" }

    input:
    val  image             // the pipeline's Docker image, e.g. crisprscreen:1.0
    path versions_script   // bin/software_versions.sh, or the project's `software_versions:` replacement

    output:
    path "${safeName(image)}.versions.tsv", emit: versions
    path "command.sh",                       emit: command_sh

    script:
    """
    ${runner(versions_script)} "${versions_script}" \\
        --image "${image}" \\
        --output "${safeName(image)}.versions.tsv" \\
        ${extraArgs('software_versions')}
    """

    stub:
    """
    printf 'image\\ttool\\tversion\\n' > "${safeName(image)}.versions.tsv"
    """
}

def safeName(image) {
    image.toString().replaceAll(/[^A-Za-z0-9._-]+/, '_')
}
