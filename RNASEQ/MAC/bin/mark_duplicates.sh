#!/bin/bash
# mark_duplicates.sh -- Picard MarkDuplicates on the filtered BAM, or a plain copy when marking is off.
# Called by MARK_DUPLICATES (modules/local/mark_duplicates.nf):
#   mark_duplicates.sh --sample S --bam S.filtered.bam --mark true --remove false --pixel-distance 2500 --java-gb G --threads N
# Arguments:
#   --sample          sample id
#   --bam             coordinate-sorted BAM (BAM_FILTER output); its index must sit next to it (BAM + ".bai") when --mark false
#   --mark            true: run MarkDuplicates; false: no marking, the BAM and its index are copied as the final BAM
#   --remove          true or false: REMOVE_DUPLICATES (config remove_duplicates, default false); used when --mark true
#   --pixel-distance  OPTICAL_DUPLICATE_PIXEL_DISTANCE (config optical_duplicate_pixel_distance, default 2500)
#   --java-gb         picard -Xmx in GB (the module passes task.memory * 0.75)
#   --threads         accepted for the common interface; Picard MarkDuplicates is single-threaded
#   other arguments   appended to the MarkDuplicates command
# Command (legacy step_3_process_bam.py block E; the legacy batch shows the same arguments):
#   picard -Xmx<G>g MarkDuplicates INPUT=S.filtered.bam OUTPUT=S.final.bam METRICS_FILE=S.picard_metrics.txt \
#     OPTICAL_DUPLICATE_PIXEL_DISTANCE=PX REMOVE_DUPLICATES=R TMP_DIR=<work dir>/picard_tmp
#   picard -Xmx<G>g BuildBamIndex INPUT=S.final.bam OUTPUT=S.final.bam.bai
# Documented deviation: the Picard image has no samtools, so the index is written with Picard BuildBamIndex
#   (legacy: samtools index; same .bai format).
# Legacy IgG branch (always remove duplicates for IgG samples) is not ported: RNA-seq libraries have no IgG controls.
# Exit 3 on an invalid --mark, --remove, --pixel-distance or --java-gb value, an empty input BAM, or an empty output.
# Outputs: S.final.bam, S.final.bam.bai, S.picard_metrics.txt (--mark true only; not written when --mark false)
set -euo pipefail

sample="" bam="" mark="" remove="" pixel_distance="" java_gb="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)          sample="$2"; shift 2 ;;
        --bam)             bam="$2"; shift 2 ;;
        --mark)            mark="$2"; shift 2 ;;
        --remove)          remove="$2"; shift 2 ;;
        --pixel-distance)  pixel_distance="$2"; shift 2 ;;
        --java-gb)         java_gb="$2"; shift 2 ;;
        --threads)         threads="$2"; shift 2 ;;
        *)                 extra+=("$1"); shift ;;
    esac
done
for v in sample bam mark remove pixel_distance java_gb; do
    [[ -n "${!v}" ]] || { echo "mark_duplicates.sh: missing --${v//_/-}" >&2; exit 3; }
done

case "${mark}" in
    true|false) ;;
    *) echo "mark_duplicates.sh: --mark must be true or false (got '${mark}')" >&2; exit 3 ;;
esac
case "${remove}" in
    true|false) ;;
    *) echo "mark_duplicates.sh: --remove must be true or false (got '${remove}')" >&2; exit 3 ;;
esac
[[ "${pixel_distance}" =~ ^[0-9]+$ ]] || { echo "mark_duplicates.sh: --pixel-distance must be an integer (got '${pixel_distance}')" >&2; exit 3; }
[[ "${java_gb}" =~ ^[0-9]+$ ]] || { echo "mark_duplicates.sh: --java-gb must be an integer (got '${java_gb}')" >&2; exit 3; }
[[ -s "${bam}" ]] || { echo "mark_duplicates.sh: empty or missing input BAM ${bam}" >&2; exit 3; }

if [[ "${mark}" == true ]]; then
    tmp_dir="${PWD}/picard_tmp"
    mkdir -p "${tmp_dir}"
    picard -Xmx"${java_gb}"g MarkDuplicates \
        INPUT="${bam}" \
        OUTPUT="${sample}.final.bam" \
        METRICS_FILE="${sample}.picard_metrics.txt" \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE="${pixel_distance}" \
        REMOVE_DUPLICATES="${remove}" \
        TMP_DIR="${tmp_dir}" \
        ${extra[@]+"${extra[@]}"}
    rm -rf "${tmp_dir}"
    picard -Xmx"${java_gb}"g BuildBamIndex INPUT="${sample}.final.bam" OUTPUT="${sample}.final.bam.bai"
else
    [[ -s "${bam}.bai" ]] || { echo "mark_duplicates.sh: index ${bam}.bai is missing or empty" >&2; exit 3; }
    cp "${bam}" "${sample}.final.bam"
    cp "${bam}.bai" "${sample}.final.bam.bai"
fi

[[ -s "${sample}.final.bam" ]] || { echo "mark_duplicates.sh: empty output ${sample}.final.bam" >&2; exit 3; }
[[ -s "${sample}.final.bam.bai" ]] || { echo "mark_duplicates.sh: empty output ${sample}.final.bam.bai" >&2; exit 3; }
