#!/bin/bash
# align.sh -- bowtie2 paired-end alignment piped straight into an (unsorted) BAM.
# Called by ALIGN (modules/local/align.nf):
#   align.sh --sample S --r1 R1 --r2 R2 --index-base PREFIX --mode local --minins 10 --maxins 2000 --threads N
# Arguments:
#   --sample      sample id (output prefix; ALIGN passes S.chunkNNN for chunks)
#   --r1, --r2    FASTQ files
#   --index-base  bare file-name prefix of the bowtie2 index (index files are staged in the work dir)
#   --mode        local (--local --very-sensitive-local) | end_to_end (--end-to-end --very-sensitive)
#   --minins      bowtie2 -I; default 10 (config bowtie2_minins); empty or 'null' = not given (no -I)
#   --maxins      bowtie2 -X (config bowtie2_maxins, default 2000)
#   --threads     bowtie2 -p and samtools -@ (task.cpus)
#   --chunk-skip N, --chunk-upto M
#                 optional, give both or neither: bowtie2 gets -s N -u M (skip the first N pairs, align the next M).
#                 Exactly one of them, or a non-integer or negative value, is an error (exit 3).
#   --            everything after it goes to bowtie2 unchanged (config bowtie2_args, default empty; the legacy
#                 toggle is '-N 1 -t'), followed by the config align_args
# Fixed flags (CUT&RUN parameters plus -X, user decision 2026-10-28): -q --no-mixed --no-discordant --dovetail
#   -I (from --minins) -X (from --maxins). The legacy -N 1 and -t are no longer fixed; they can only come back
#   through the -- arguments. No --no-unal: unaligned records stay in the BAM for the QC steps.
# --reorder: records come out in input order whatever -p is. Alignments are unchanged, but without it the order
#   depends on thread timing, and the order decides which of two tied duplicates sambamba markdup keeps
#   (1M-pair test, 2026-10-28: same aligned records, 0.2-0.5% of final-BAM records differed between two runs).
# Chunking: the chunk BAMs, concatenated in chunk order, equal the records of one unchunked run. bowtie2 seeds its
#   pseudo-random generator per read from the read name, sequence, qualities and --seed, so a read aligns the same
#   way whatever else is in the input; a chunk with -s N -u M therefore gives the same records as those pairs in the
#   full run. ALIGN_MERGE (bin/align_merge.sh) joins the chunks; no temporary FASTQ is written.
# Never writes a SAM: bowtie2 stdout | samtools view -b.
# Outputs: S.aligned.bam (unsorted, bowtie2 order, includes unmapped records), S.bowtie2.log (bowtie2 stderr summary)
# Exit 3 on missing arguments, an unknown --mode, no bowtie2 index files for --index-base, or a bad --chunk-skip /
#   --chunk-upto (only one given, not an integer >= 0).
set -euo pipefail

sample="" r1="" r2="" index_base="" mode="" minins="10" maxins="" threads=1 chunk_skip="" chunk_upto="" extra=()
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
        --chunk-skip) chunk_skip="$2"; shift 2 ;;
        --chunk-upto) chunk_upto="$2"; shift 2 ;;
        --)           shift; extra+=("$@"); break ;;
        *)            extra+=("$1"); shift ;;
    esac
done

for v in sample r1 r2 index_base mode maxins; do
    [[ -n "${!v}" ]] || { echo "align.sh: missing --${v//_/-}" >&2; exit 3; }
done
case "${mode}" in
    local)      mode_args=(--local --very-sensitive-local) ;;
    end_to_end) mode_args=(--end-to-end --very-sensitive) ;;
    *) echo "align.sh: --mode must be local or end_to_end (got '${mode}')" >&2; exit 3 ;;
esac
ins_args=(-X "${maxins}")
if [[ -n "${minins}" && "${minins}" != "null" ]]; then
    ins_args+=(-I "${minins}")
fi
chunk_args=()
if [[ -n "${chunk_skip}" || -n "${chunk_upto}" ]]; then
    [[ -n "${chunk_skip}" && -n "${chunk_upto}" ]] || { echo "align.sh: --chunk-skip and --chunk-upto must be given together" >&2; exit 3; }
    for v in chunk_skip chunk_upto; do
        [[ "${!v}" =~ ^[0-9]+$ ]] || { echo "align.sh: --${v//_/-} must be an integer >= 0 (got '${!v}')" >&2; exit 3; }
    done
    chunk_args=(-s "${chunk_skip}" -u "${chunk_upto}")
fi
ls "${index_base}"*.bt2* >/dev/null 2>&1 || { echo "align.sh: no bowtie2 index files for prefix '${index_base}'" >&2; exit 3; }

bowtie2 -p "${threads}" "${mode_args[@]}" \
    -q --no-mixed --no-discordant --dovetail \
    --reorder \
    "${ins_args[@]}" \
    ${chunk_args[@]+"${chunk_args[@]}"} \
    -x "${index_base}" -1 "${r1}" -2 "${r2}" \
    ${extra[@]+"${extra[@]}"} \
    2> "${sample}.bowtie2.log" \
  | samtools view -b -@ "${threads}" -o "${sample}.aligned.bam" -

cat "${sample}.bowtie2.log" >&2
