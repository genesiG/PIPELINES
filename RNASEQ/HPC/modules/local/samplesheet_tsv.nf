include { runner; extraArgs } from './functions'

process SAMPLESHEET_TSV {
    tag "samplesheet"
    container '/project/rephilllab/nf_containers/cutandrun-py_1.2.sif'
    cache 'deep'
    publishDir "${params.outdir}", mode: params.publish_mode, pattern: "samplesheet.tsv"
    publishDir "${params.outdir}/log", mode: params.publish_mode, pattern: "command.sh", saveAs: { "samplesheet_tsv_command.sh" }

    input:
    path samplesheet   // samplesheet.csv (params.samplesheet)
    val keep_ids       // comma list of the non-excluded sample ids
    path step_script   // bin/samplesheet_tsv.py, or the project's `samplesheet_tsv` script replacement

    output:
    path "samplesheet.tsv", emit: tsv
    path "command.sh",      emit: command_sh

    script:
    """
    ${runner(step_script)} "${step_script}" \\
        --samplesheet "${samplesheet}" \\
        --keep "${keep_ids}" \\
        --out samplesheet.tsv \\
        ${extraArgs('samplesheet_tsv')}
    """

    stub:
    """
    printf 'sample_id\\tcondition\\nstub1\\tA\\n' > samplesheet.tsv
    """
}
