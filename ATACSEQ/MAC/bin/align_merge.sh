#!/bin/bash
# align_merge.sh -- join the chunk alignments of one sample (ALIGN with align_chunk_pairs > 0) into S.aligned.bam and
#   S.bowtie2.log.
# Called by ALIGN_MERGE (modules/local/align_merge.nf):
#   align_merge.sh --sample S --bams S.chunk000.aligned.bam,S.chunk001.aligned.bam \
#       --logs S.chunk000.bowtie2.log,S.chunk001.bowtie2.log --log-merger merge_bowtie2_logs.py --threads N
# Arguments:
#   --sample      sample id (output prefix)
#   --bams        comma list of the chunk BAMs, in chunk order; kept as given (not re-sorted)
#   --logs        comma list of the chunk bowtie2 logs, same order and same count as --bams
#   --log-merger  path of merge_bowtie2_logs.py
#   --threads     accepted for the module's convention; samtools cat does not use it
# BAM: one chunk = copied to S.aligned.bam; several = `samtools cat -o S.aligned.bam` in the order given. The header
#   of the output is the header of the FIRST chunk (samtools cat), so every chunk must carry the same @SQ lines.
#   The records of the chunks, concatenated in order, are the records of one unchunked run (see bin/align.sh).
# Log: merge_bowtie2_logs.py sums the paired-end bowtie2 summaries of the chunks into S.bowtie2.log (also for one
#   chunk, so that the output has the summary only, without the [WARNING] and timing lines).
# Outputs: S.aligned.bam, S.bowtie2.log (the latter is also written to stderr, as align.sh does)
# Exit 3 on missing arguments, a BAM or log that is missing or empty, or a different number of BAMs and logs.
#   merge_bowtie2_logs.py exits 3 when a log cannot be parsed.
set -euo pipefail

sample="" bams="" logs="" log_merger="" threads=1
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)     sample="$2"; shift 2 ;;
        --bams)       bams="$2"; shift 2 ;;
        --logs)       logs="$2"; shift 2 ;;
        --log-merger) log_merger="$2"; shift 2 ;;
        --threads)    threads="$2"; shift 2 ;;
        *)            shift ;;
    esac
done
[[ -n "${sample}" && -n "${bams}" && -n "${logs}" && -n "${log_merger}" ]] \
    || { echo "align_merge.sh: missing --sample, --bams, --logs or --log-merger" >&2; exit 3; }

IFS=',' read -r -a bam_list <<< "${bams}"
IFS=',' read -r -a log_list <<< "${logs}"
(( ${#bam_list[@]} == ${#log_list[@]} )) \
    || { echo "align_merge.sh: ${#bam_list[@]} BAMs but ${#log_list[@]} logs" >&2; exit 3; }
for f in "${bam_list[@]}" "${log_list[@]}"; do
    [[ -s "${f}" ]] || { echo "align_merge.sh: missing or empty chunk file: '${f}'" >&2; exit 3; }
done

if (( ${#bam_list[@]} == 1 )); then
    cp "${bam_list[0]}" "${sample}.aligned.bam"
else
    samtools cat -o "${sample}.aligned.bam" "${bam_list[@]}"
fi
python3 "${log_merger}" --out "${sample}.bowtie2.log" "${log_list[@]}"

cat "${sample}.bowtie2.log" >&2
