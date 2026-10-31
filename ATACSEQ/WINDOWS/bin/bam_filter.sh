#!/bin/bash
# bam_filter.sh -- sort, stats of the unfiltered alignment, MAPQ + flag filter, optional contig removal.
# Called by BAM_FILTER (modules/local/bam_filter.nf):
#   bam_filter.sh --sample S --bam S.aligned.bam --mapq 30 --remove-chroms chrM --threads N
# Arguments:
#   --sample          sample id (output prefix)
#   --bam             unsorted aligned BAM (ALIGN output)
#   --mapq            samtools view -q (config mapq, default 30)
#   --remove-chroms   comma-separated contig names to drop (config remove_chroms, default chrM);
#                     empty = keep every contig (legacy)
#   --threads         samtools -@ (task.cpus)
#   other arguments   appended to the filtering `samtools view`.
# Steps: coordinate sort + index; flagstat and idxstats of that unfiltered sorted BAM (the "raw" stats);
#   then samtools view -b -q MAPQ -f 3 -F 2828, restricted with -L to the kept contigs when --remove-chroms is set.
#   The -L BED ("chrom 0 length") is built from the raw idxstats: every contig except the removed ones and '*'.
#   A proper pair lies on one contig, so both mates are removed together.
# Flags copied from legacy step_3_process_bam.py (paired-end):
#   "-f 1 -f 2 -F 4 -F 8 -F 256 -F 512 -F 1024 -F 2048". samtools ORs repeated options, so
#   -f 1 -f 2 = -f 3 (paired + proper pair); -F 4,8,256,512,2048 = -F 2828.
#   The legacy -F 1024 is omitted: no duplicates are marked at this stage (DEDUP follows), so it was a no-op here.
# MAPQ: legacy wrote the literal `-q mapq`, which samtools reads with atoi (0), so legacy kept every MAPQ.
#   `--mapq 0` reproduces that; the default 30 is the fix.
# Legacy kept chrM; the default here removes it (`--remove-chroms ''` reproduces that).
# Warns on stderr for a removed contig that is not in the BAM. Exit 3 when the input BAM is missing or empty,
# or when no contig is left after the removal.
# Outputs: S.filtered.bam, S.filtered.bam.bai, S.raw.flagstat.txt, S.raw.idxstats.txt
set -euo pipefail

sample="" bam="" mapq="" remove_chroms="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)         sample="$2"; shift 2 ;;
        --bam)            bam="$2"; shift 2 ;;
        --mapq)           mapq="$2"; shift 2 ;;
        --remove-chroms)  remove_chroms="$2"; shift 2 ;;
        --threads)        threads="$2"; shift 2 ;;
        *)                extra+=("$1"); shift ;;
    esac
done
for v in sample bam mapq; do
    [[ -n "${!v}" ]] || { echo "bam_filter.sh: missing --${v}" >&2; exit 3; }
done
[[ -s "${bam}" ]] || { echo "bam_filter.sh: empty or missing input BAM ${bam}" >&2; exit 3; }

sorted="${sample}.sorted_unfiltered.tmp.bam"
kept_bed="${sample}.kept_contigs.tmp.bed"
samtools sort -m "${ATAC_SORT_MEM:-1G}" -@ "${threads}" -o "${sorted}" "${bam}"
samtools index -@ "${threads}" "${sorted}"
samtools flagstat -@ "${threads}" "${sorted}" > "${sample}.raw.flagstat.txt"
samtools idxstats "${sorted}" > "${sample}.raw.idxstats.txt"

region_args=()
remove_chroms="${remove_chroms// /}"
if [[ -n "${remove_chroms}" ]]; then
    IFS=',' read -ra drop <<< "${remove_chroms}"
    for c in "${drop[@]}"; do
        [[ -n "${c}" ]] || continue
        awk -F'\t' -v c="${c}" '$1 == c { found = 1 } END { exit found ? 0 : 1 }' "${sample}.raw.idxstats.txt" \
            || echo "WARNING: contig '${c}' (--remove-chroms) is not in ${bam}" >&2
    done
    awk -F'\t' -v drop="${remove_chroms}" '
        BEGIN { n = split(drop, a, ","); for (i = 1; i <= n; i++) removed[a[i]] = 1 }
        $1 != "*" && !($1 in removed) { print $1 "\t0\t" $2 }' "${sample}.raw.idxstats.txt" > "${kept_bed}"
    [[ -s "${kept_bed}" ]] || { echo "bam_filter.sh: no contig left in ${bam} after removing '${remove_chroms}'" >&2; exit 3; }
    region_args=(-L "${kept_bed}")
fi

samtools view -b -@ "${threads}" -q "${mapq}" -f 3 -F 2828 \
    ${region_args[@]+"${region_args[@]}"} ${extra[@]+"${extra[@]}"} \
    -o "${sample}.filtered.bam" "${sorted}"
samtools index -@ "${threads}" "${sample}.filtered.bam"
rm -f "${sorted}" "${sorted}.bai" "${kept_bed}"
