/*
 * Writes METHODS.md (bin/write_methods.py): the Materials and Methods for
 * this run — only the steps it performed, with the parameters and software
 * versions they actually used — into the project directory (the launch
 * directory) and into results/, so a run with its own outdir (e.g. a second
 * config.yaml in the same project) keeps its own copy.
 *
 * Runs last: its `n_finished` input is the count of every other step's final
 * output (reports, plots), and its evidence comes from the remaining steps.
 * So METHODS.md is only (re)written by a run that succeeded; a failed run
 * leaves the previous one, which names the run it describes, in place.
 *
 * Inputs (all small):
 *   run_info   main.nf runInfo(): the run's resolved configuration
 *   evidence   what the steps recorded (tool logs, anchors, libraries,
 *              per-sample reads, contrasts, software versions), staged
 *              together — main.nf guarantees unique file names
 *   discover_script  the discover_anchor script this run used, whose
 *              argparse defaults write_methods.py reads (so METHODS.md
 *              quotes the parameters that actually applied)
 *
 * executor 'local': a few seconds of Python inside the Nextflow head job
 * (itself an LSF job on a compute node), not worth an LSF queue wait.
 * mode 'copy' regardless of publish_mode: a symlink into the work
 * directory would break once that is cleaned.
 */
include { runner; extraArgs } from './functions'

process METHODS {
    tag "${params.project_name}"
    container 'quay.io/biocontainers/cutadapt:5.2--py39hbcbf7aa_0'   // python3 stdlib only; reused to avoid an extra pull
    executor 'local'
    publishDir "${workflow.launchDir}", mode: 'copy', pattern: "METHODS.md"
    publishDir "${params.outdir}", mode: 'copy', pattern: "METHODS.md"
    publishDir "${params.metadata_dir}/pipeline_info", mode: 'copy', pattern: "{run_info.json,software_versions.tsv}"
    publishDir "${params.metadata_dir}/pipeline_info/log", mode: 'copy', pattern: "command.sh", saveAs: { "methods_command.sh" }

    input:
    path run_info, stageAs: 'inputs/run_info.json'
    path evidence, stageAs: 'evidence/*'
    val  n_finished        // only so this runs after every other step
    path discover_script, stageAs: 'inputs/*'
    path methods_script    // bin/write_methods.py, or the project's `write_methods:` replacement

    output:
    path "METHODS.md",            emit: methods
    path "run_info.json",         emit: run_info
    path "software_versions.tsv", emit: versions
    path "command.sh",            emit: command_sh

    script:
    """
    ${runner(methods_script)} "${methods_script}" \\
        --run-info inputs/run_info.json \\
        --evidence-dir evidence \\
        --discover-script "${discover_script}" \\
        --output METHODS.md \\
        --run-info-out run_info.json \\
        --versions-out software_versions.tsv \\
        ${extraArgs('write_methods')}
    """

    stub:   // runs the real (stdlib-only) writer on the stub outputs, so -stub-run exercises it too
    """
    ${runner(methods_script)} "${methods_script}" \\
        --run-info inputs/run_info.json \\
        --evidence-dir evidence \\
        --discover-script "${discover_script}" \\
        --output METHODS.md \\
        --run-info-out run_info.json \\
        --versions-out software_versions.tsv \\
        ${extraArgs('write_methods')}
    """
}
