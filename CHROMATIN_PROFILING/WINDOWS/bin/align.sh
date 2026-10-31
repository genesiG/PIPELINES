#!/bin/bash
# align.sh -- bowtie2 paired-end alignment piped straight into an (unsorted) BAM.
# Called by ALIGN (modules/local/align.nf):
#   align.sh --sample S --r1 R1 --r2 R2 --index-base PREFIX --mode local --threads N
# Arguments:
#   --sample      sample id
#   --r1, --r2    FASTQ files
#   --index-base  bare file-name prefix of the bowtie2 index (index files are staged in the work dir)
#   --mode        local (--local --very-sensitive-local) | end_to_end (--end-to-end --very-sensitive)
#   --minins/--maxins  optional bowtie2 -I / -X; omitted when empty, so bowtie2 defaults (-I 0, -X 500) apply, as in the manuscript
#   --threads     bowtie2 -p and samtools -@ (task.cpus)
#   other arguments are appended to bowtie2.
# Fixed: --no-mixed --no-discordant --dovetail -q (legacy: -q --dovetail; no --no-unal).
# Never writes a SAM: bowtie2 stdout | samtools view -b.
# Outputs: S.aligned.bam (unsorted), S.bowtie2.log (bowtie2 stderr summary)
set -euo pipefail

sample="" r1="" r2="" index_base="" mode="" minins="" maxins="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)     sample="$2"; shift 2 ;;
        --r1)         r1="$2"; shift 2 ;;
        --r2)         r2="$2"; shift 2 ;;
        --index-base) index_base="$2"; shift 2 ;;
        --mode)       mode="$2"; shift 2 ;;
        --minins)     minins="$2"; shift 2 ;;
        --maxins)     maxins="$2"; shift 2 ;;
        --threads)    threads="$2"; shift 2 ;;
        *)            extra+=("$1"); shift ;;
    esac
done

for v in sample r1 r2 index_base mode; do
    [[ -n "${!v}" ]] || { echo "align.sh: missing --${v//_/-}" >&2; exit 3; }
done
case "${mode}" in
    local)      mode_args=(--local --very-sensitive-local) ;;
    end_to_end) mode_args=(--end-to-end --very-sensitive) ;;
    *) echo "align.sh: --mode must be local or end_to_end (got '${mode}')" >&2; exit 3 ;;
esac
ls "${index_base}"*.bt2* >/dev/null 2>&1 || { echo "align.sh: no bowtie2 index files for prefix '${index_base}'" >&2; exit 3; }
# -I/-X are optional: "null" or "NULL" (from the Nextflow config) counts as empty.
if [[ "${minins}" == "null" || "${minins}" == "NULL" ]]; then minins=""; fi
if [[ "${maxins}" == "null" || "${maxins}" == "NULL" ]]; then maxins=""; fi
ins_args=()
if [[ -n "${minins}" ]]; then ins_args+=(-I "${minins}"); fi
if [[ -n "${maxins}" ]]; then ins_args+=(-X "${maxins}"); fi

bowtie2 -p "${threads}" "${mode_args[@]}" \
    -q --no-mixed --no-discordant --dovetail \
    ${ins_args[@]+"${ins_args[@]}"} \
    -x "${index_base}" -1 "${r1}" -2 "${r2}" \
    ${extra[@]+"${extra[@]}"} \
    2> "${sample}.bowtie2.log" \
  | samtools view -b -@ "${threads}" -o "${sample}.aligned.bam" -

cat "${sample}.bowtie2.log" >&2
