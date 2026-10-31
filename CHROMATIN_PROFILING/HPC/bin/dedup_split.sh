#!/bin/bash
# dedup_split.sh -- split the duplicate-marked BAM into spike-in and target (human) reads.
# Called by DEDUP_SPLIT (modules/local/dedup_split.nf):
#   dedup_split.sh --sample S --bam S.markdup.bam --is-control true|false --remove-duplicates true|false
#       --spikein-dedup true|false --spike-prefix ecoli_chr --mito-chrom chrM --threads N
# The split happens BEFORE any human filter (blacklist), as in legacy step_3.
# Contig lists come from `samtools idxstats`:
#   spike contigs  = names starting with --spike-prefix
#   target contigs = contigs with >0 mapped reads (legacy whitelist), excluding '*', --mito-chrom and spike contigs
# Spike BAM: coordinate-sorted; `-F 1024` when --spikein-dedup true, otherwise duplicates kept.
# Target BAM: `-F 1024` when --remove-duplicates true OR --is-control true; NAME-sorted (samtools sort -n).
# Outputs: S.spikein.bam, S.spikein.bam.bai, S.target.namesorted.bam,
#   S.spikein.idxstats.txt (samtools idxstats of the spike-in BAM; input of norm_spikein.py),
#   S.dedup_counts.tsv (sample_id, stage, reads: target_before_dedup, target_after_dedup,
#   spikein_before_dedup, spikein_after_dedup; alignment records)
# Exit 3 when the spike prefix is empty or no target reads remain.
set -euo pipefail

sample="" bam="" is_control="false" rmdup="true" spike_dedup="true" spike_prefix="" mito="chrM" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)            sample="$2"; shift 2 ;;
        --bam)               bam="$2"; shift 2 ;;
        --is-control)        is_control="$2"; shift 2 ;;
        --remove-duplicates) rmdup="$2"; shift 2 ;;
        --spikein-dedup)     spike_dedup="$2"; shift 2 ;;
        --spike-prefix)      spike_prefix="$2"; shift 2 ;;
        --mito-chrom)        mito="$2"; shift 2 ;;
        --threads)           threads="$2"; shift 2 ;;
        *)                   extra+=("$1"); shift ;;
    esac
done
[[ -n "${sample}" && -s "${bam}" ]] || { echo "dedup_split.sh: missing --sample or empty --bam" >&2; exit 3; }
[[ -n "${spike_prefix}" ]] || { echo "dedup_split.sh: --spike-prefix is empty" >&2; exit 3; }

samtools idxstats "${bam}" > idxstats.tmp
mapfile -t spike_contigs < <(awk -v p="${spike_prefix}" 'index($1,p)==1 {print $1}' idxstats.tmp)
mapfile -t target_contigs < <(awk -v p="${spike_prefix}" -v m="${mito}" \
    '$1!="*" && $1!=m && index($1,p)!=1 && $3>0 {print $1}' idxstats.tmp)

[[ ${#target_contigs[@]} -gt 0 ]] || { echo "dedup_split.sh: no target contigs with reads in ${bam}" >&2; exit 3; }
[[ ${#spike_contigs[@]} -gt 0 ]] || echo "WARNING: no contigs with prefix '${spike_prefix}'; spike-in BAM will be header-only" >&2

# ---- (1) spike-in BAM (coordinate-sorted because the input is) ----
spike_flag=()
[[ "${spike_dedup}" == "true" ]] && spike_flag=(-F 1024)
if [[ ${#spike_contigs[@]} -gt 0 ]]; then
    samtools view -b -h -@ "${threads}" ${spike_flag[@]+"${spike_flag[@]}"} "${bam}" "${spike_contigs[@]}" \
        -o "${sample}.spikein.bam"
    spike_before=$(samtools view -c -@ "${threads}" "${bam}" "${spike_contigs[@]}")
else
    samtools view -b -H "${bam}" -o "${sample}.spikein.bam"
    spike_before=0
fi
samtools index -@ "${threads}" "${sample}.spikein.bam"
samtools idxstats "${sample}.spikein.bam" > "${sample}.spikein.idxstats.txt"
spike_after=$(samtools view -c -@ "${threads}" "${sample}.spikein.bam")

# ---- (2) target BAM, name-sorted ----
target_flag=()
if [[ "${rmdup}" == "true" || "${is_control}" == "true" ]]; then target_flag=(-F 1024); fi
target_before=$(samtools view -c -@ "${threads}" "${bam}" "${target_contigs[@]}")
samtools view -b -h -@ "${threads}" ${target_flag[@]+"${target_flag[@]}"} ${extra[@]+"${extra[@]}"} \
    "${bam}" "${target_contigs[@]}" \
  | samtools sort -n -m 1G -@ "${threads}" -o "${sample}.target.namesorted.bam" -
target_after=$(samtools view -c -@ "${threads}" "${sample}.target.namesorted.bam")

# ---- (3) counts ----
{
    printf 'sample_id\tstage\treads\n'
    printf '%s\ttarget_before_dedup\t%s\n' "${sample}" "${target_before}"
    printf '%s\ttarget_after_dedup\t%s\n'  "${sample}" "${target_after}"
    printf '%s\tspikein_before_dedup\t%s\n' "${sample}" "${spike_before}"
    printf '%s\tspikein_after_dedup\t%s\n'  "${sample}" "${spike_after}"
} > "${sample}.dedup_counts.tsv"
rm -f idxstats.tmp
