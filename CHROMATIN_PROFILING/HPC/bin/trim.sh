#!/bin/bash
# trim.sh -- Trimmomatic paired-end adapter trimming.
# Called by TRIM (modules/local/trim.nf):
#   trim.sh --r1 R1.fq.gz --r2 R2.fq.gz --sample S --adapter-set TruSeq3-PE-2 --minlen 25 --threads N
# Arguments:
#   --r1, --r2        input paired FASTQ files
#   --sample          sample id (output prefix)
#   --adapter-set     adapter fasta basename: TruSeq3-PE-2 (CUT&RUN) or NexteraPE-PE (CUT&Tag)
#   --minlen          Trimmomatic MINLEN
#   --threads         Trimmomatic -threads (task.cpus)
#   other arguments   appended to the trimmomatic step list (after MINLEN)
# Legacy parameters (step_1_trimmomatic.py, non-ATAC): ILLUMINACLIP:<fa>:2:30:10:8:True MINLEN:25
# (no LEADING/TRAILING/SLIDINGWINDOW for CUT&RUN/CUT&Tag). Unpaired outputs are discarded (/dev/null).
# The adapter fasta is located inside the container.
# Outputs: S.R1.trim.fastq.gz, S.R2.trim.fastq.gz, S.trimmomatic.log (trimmomatic stderr)
set -euo pipefail

r1="" r2="" sample="" adapter_set="" minlen="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --r1)          r1="$2"; shift 2 ;;
        --r2)          r2="$2"; shift 2 ;;
        --sample)      sample="$2"; shift 2 ;;
        --adapter-set) adapter_set="$2"; shift 2 ;;
        --minlen)      minlen="$2"; shift 2 ;;
        --threads)     threads="$2"; shift 2 ;;
        *)             extra+=("$1"); shift ;;
    esac
done

for v in r1 r2 sample adapter_set minlen; do
    [[ -n "${!v}" ]] || { echo "trim.sh: missing --${v//_/-}" >&2; exit 3; }
done

tm="$(command -v trimmomatic)" || { echo "trim.sh: trimmomatic not found" >&2; exit 3; }
roots=("$(dirname "$(readlink -f "${tm}")")/..")
[[ -n "${CONDA_PREFIX:-}" ]] && roots+=("${CONDA_PREFIX}/share")
adapter_fa=""
for root in "${roots[@]}"; do
    adapter_fa="$(find "${root}" -name "${adapter_set}.fa" 2>/dev/null | head -1 || true)"
    [[ -n "${adapter_fa}" ]] && break
done
[[ -n "${adapter_fa}" ]] || { echo "trim.sh: adapter file ${adapter_set}.fa not found in the container" >&2; exit 3; }
echo "Adapter file: ${adapter_fa}" >&2

trimmomatic PE -threads "${threads}" \
    "${r1}" "${r2}" \
    "${sample}.R1.trim.fastq.gz" /dev/null \
    "${sample}.R2.trim.fastq.gz" /dev/null \
    "ILLUMINACLIP:${adapter_fa}:2:30:10:8:True" \
    "MINLEN:${minlen}" ${extra[@]+"${extra[@]}"} \
    2> "${sample}.trimmomatic.log"

cat "${sample}.trimmomatic.log" >&2
