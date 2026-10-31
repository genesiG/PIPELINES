#!/bin/bash
# bam_filter.sh -- raw statistics of the ALIGN BAM, MAPQ + flag filter, optional contig removal, sort and index.
# Called by BAM_FILTER (modules/local/bam_filter.nf):
#   bam_filter.sh --sample S --bam S.aligned.bam --mapq 30 --remove-chroms "" --threads N [extra samtools view args...]
# Arguments:
#   --sample          sample id (output prefix)
#   --bam             coordinate-sorted, indexed BAM (ALIGN output; its .bai must sit next to it)
#   --mapq            samtools view -q (config mapq, default 30)
#   --remove-chroms   comma-separated contig names to drop (config remove_chroms, default empty = keep every contig)
#   --threads         samtools -@ (task.cpus)
#   other arguments   appended to the filtering `samtools view`
# Steps:
#   1. S.raw.flagstat.txt and S.raw.idxstats.txt from the ALIGN BAM, before any filter. If the .bai is missing, the BAM
#      is copied into the work dir and indexed there (the copy is removed at the end).
#   2. samtools view -@ N -h -b -q MAPQ -f 1 -f 2 -F 4 -F 8 -F 256 -F 512 -F 1024 -F 2048 [-L kept.bed] IN -o S.filtered.unsorted.tmp.bam
#      The -L BED ("contig 0 length") is built from the raw idxstats: every contig except the removed ones and '*'.
#      A proper pair lies on one contig, so both mates are removed together.
#   3. samtools sort -m 1G -@ N -o S.filtered.bam; samtools index; S.filtered.flagstat.txt and S.filtered.idxstats.txt.
# Flags, written exactly as the legacy step_3_process_bam.py block A command (-f 1 -f 2 -F 4 -F 8 -F 256 -F 512 -F 1024 -F 2048):
#   samtools 1.23 (sam_view.c lines 1017 and 1028: settings.flag_on |= and settings.flag_off |=) ORs repeated -f and -F,
#   so the filter equals -f 3 -F 3852: paired and proper pair; not unmapped, mate unmapped, secondary, QC fail,
#   duplicate (1024) or supplementary. 1024 is kept here, unlike the ATACSEQ port, to be identical to the legacy filter.
# Legacy differences: legacy took its pre-filter flagstat/idxstats after the MAPQ filter (the file it wrote as qc.bam);
#   here they come from the ALIGN BAM (spec). Legacy dropped chrM with a region list (KEEP_CHROMS); here the -L BED is
#   built from the raw idxstats.
# Warns on stderr for a removed contig that is not in the BAM. Exit 3 when the input BAM is missing or empty, or when
#   no contig is left after the removal.
# Outputs: S.filtered.bam, S.filtered.bam.bai, S.filtered.flagstat.txt, S.filtered.idxstats.txt,
#          S.raw.flagstat.txt, S.raw.idxstats.txt
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

in_bam="${bam}"
copy_bam=""
if [[ ! -s "${bam}.bai" ]]; then
    copy_bam="${sample}.input_copy.tmp.bam"
    cp "${bam}" "${copy_bam}"
    samtools index -@ "${threads}" "${copy_bam}"
    in_bam="${copy_bam}"
fi

unsorted="${sample}.filtered.unsorted.tmp.bam"
kept_bed="${sample}.kept_contigs.tmp.bed"
samtools flagstat -@ "${threads}" "${in_bam}" > "${sample}.raw.flagstat.txt"
samtools idxstats "${in_bam}" > "${sample}.raw.idxstats.txt"

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

samtools view -@ "${threads}" -h -b -q "${mapq}" -f 1 -f 2 -F 4 -F 8 -F 256 -F 512 -F 1024 -F 2048 \
    ${region_args[@]+"${region_args[@]}"} ${extra[@]+"${extra[@]}"} \
    "${in_bam}" -o "${unsorted}"

samtools sort -m 1G -@ "${threads}" -o "${sample}.filtered.bam" "${unsorted}"
samtools index -@ "${threads}" "${sample}.filtered.bam"
samtools flagstat -@ "${threads}" "${sample}.filtered.bam" > "${sample}.filtered.flagstat.txt"
samtools idxstats "${sample}.filtered.bam" > "${sample}.filtered.idxstats.txt"

rm -f "${unsorted}" "${kept_bed}"
if [[ -n "${copy_bam}" ]]; then
    rm -f "${copy_bam}" "${copy_bam}.bai"
fi
