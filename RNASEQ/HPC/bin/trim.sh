#!/bin/bash
# trim.sh -- Trimmomatic paired-end adapter trimming (RNA-seq settings, all set from config.yaml).
# Called by TRIM (modules/local/trim.nf):
#   trim.sh --r1 R1.fastq.gz --r2 R2.fastq.gz --sample S --adapter-set TruSeq3-PE-2 \
#       --illuminaclip 2:30:10:8:True --steps "" --minlen 25 --threads N --java-gb G
# Arguments:
#   --r1, --r2        input paired FASTQ files
#   --sample          sample id (output prefix)
#   --adapter-set     adapter fasta basename, e.g. TruSeq3-PE-2 (located inside the container at run time)
#   --illuminaclip    ILLUMINACLIP parameters after the adapter file, seedMismatches:palindromeClipThreshold:
#                     simpleClipThreshold[:minAdapterLength:keepBothReads] (config trim_illuminaclip, default 2:30:10:8:True)
#   --steps           further Trimmomatic steps, space separated; each word becomes one argument (config trim_steps,
#                     default empty: legacy had none); may be empty
#   --minlen          Trimmomatic MINLEN (config trim_minlen, default 25)
#   --threads         Trimmomatic -threads (task.cpus)
#   --java-gb         optional; exports _JAVA_OPTIONS=-Xmx<G>g (the module passes task.memory * 0.75; legacy used -Xmx28g)
#   other arguments   appended to the trimmomatic command (after MINLEN)
# Legacy (step_1_trimmomatic.py, paired-end, non-ATAC branch): trimmomatic PE -threads 8 R1 R2 o1 o1u o2 o2u
#   ILLUMINACLIP:<adapters>/TruSeq3-PE-2.fa:2:30:10:8:True MINLEN:25, with trimmomatic 0.40 in the chipseq conda env.
#   The legacy command has no -phred flag; none is added here (Trimmomatic detects phred33).
# Unpaired outputs are discarded (/dev/null); legacy wrote them to disk. They are not used by any later step, so there is
#   no effect on counts.
# The adapter fasta is located inside the container (no hard-coded build string).
# Exit 3 when an input FASTQ is missing or empty, the adapter file is not found, or a required option is missing.
# Outputs: S.R1.trim.fastq.gz, S.R2.trim.fastq.gz, S.trimmomatic.log (trimmomatic stderr)
set -euo pipefail

r1="" r2="" sample="" adapter_set="" illuminaclip="" steps="" minlen="" threads=1 java_gb="" extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --r1)           r1="$2"; shift 2 ;;
        --r2)           r2="$2"; shift 2 ;;
        --sample)       sample="$2"; shift 2 ;;
        --adapter-set)  adapter_set="$2"; shift 2 ;;
        --illuminaclip) illuminaclip="$2"; shift 2 ;;
        --steps)        steps="$2"; shift 2 ;;
        --minlen)       minlen="$2"; shift 2 ;;
        --threads)      threads="$2"; shift 2 ;;
        --java-gb)      java_gb="$2"; shift 2 ;;
        *)              extra+=("$1"); shift ;;
    esac
done

for v in r1 r2 sample adapter_set illuminaclip minlen; do
    [[ -n "${!v}" ]] || { echo "trim.sh: missing --${v//_/-}" >&2; exit 3; }
done
for v in r1 r2; do
    [[ -s "${!v}" ]] || { echo "trim.sh: empty or missing input FASTQ ${!v}" >&2; exit 3; }
done

if [[ -n "${java_gb}" ]]; then
    [[ "${java_gb}" =~ ^[0-9]+$ ]] || { echo "trim.sh: --java-gb must be an integer (got '${java_gb}')" >&2; exit 3; }
    export _JAVA_OPTIONS="-Xmx${java_gb}g"
fi

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

# Split the step string on whitespace: one Trimmomatic argument per word.
read -ra step_list <<< "${steps}"

trimmomatic PE -threads "${threads}" \
    "${r1}" "${r2}" \
    "${sample}.R1.trim.fastq.gz" /dev/null \
    "${sample}.R2.trim.fastq.gz" /dev/null \
    "ILLUMINACLIP:${adapter_fa}:${illuminaclip}" \
    ${step_list[@]+"${step_list[@]}"} \
    "MINLEN:${minlen}" ${extra[@]+"${extra[@]}"} \
    2> "${sample}.trimmomatic.log"

cat "${sample}.trimmomatic.log" >&2
