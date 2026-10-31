#!/bin/bash
# bam_finalize.sh -- coordinate-sort the filtered target BAM, index, stats, orphan check.
# Called by BAM_FINALIZE (modules/local/bam_finalize.nf):
#   bam_finalize.sh --sample S --bam S.noblacklist.namesorted.bam --threads N
# Arguments: --sample, --bam (name-sorted), --threads; other arguments go to `samtools sort`.
# Outputs: S.final.bam, S.final.bam.bai, S.final.flagstat.txt, S.final.idxstats.txt,
#   S.orphan_check.txt (single line OK; if any read name does not occur exactly twice among
#   paired reads the offenders are written there and the script exits 3)
set -euo pipefail

sample="" bam="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)  sample="$2"; shift 2 ;;
        --bam)     bam="$2"; shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        *)         extra+=("$1"); shift ;;
    esac
done
[[ -n "${sample}" && -s "${bam}" ]] || { echo "bam_finalize.sh: missing --sample or empty --bam" >&2; exit 3; }

samtools sort -m "${ATAC_SORT_MEM:-1G}" -@ "${threads}" ${extra[@]+"${extra[@]}"} -o "${sample}.final.bam" "${bam}"
samtools index -@ "${threads}" "${sample}.final.bam"
samtools flagstat -@ "${threads}" "${sample}.final.bam" > "${sample}.final.flagstat.txt"
samtools idxstats "${sample}.final.bam" > "${sample}.final.idxstats.txt"

# samtools collate groups reads by name through 128 temporary files in the work dir, so memory stays bounded. (The
# image's sort is BusyBox: it held every read name in memory, 14 GB for 211 M pairs, and has no -S/-T.)
samtools collate -O -u -n 128 -@ "${threads}" "${sample}.final.bam" "${sample}.collate_tmp" \
    | samtools view -f 1 - | cut -f1 | uniq -c | awk '$1!=2' > "${sample}.orphan_check.txt"
if [[ -s "${sample}.orphan_check.txt" ]]; then
    echo "bam_finalize.sh: orphan reads (read name count != 2) in ${sample}.final.bam; see ${sample}.orphan_check.txt" >&2
    head -5 "${sample}.orphan_check.txt" >&2
    exit 3
fi
echo OK > "${sample}.orphan_check.txt"
