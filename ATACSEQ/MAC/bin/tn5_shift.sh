#!/bin/bash
# tn5_shift.sh -- Tn5 insertion-site shift (+4 on the forward strand, -5 on the reverse) with deepTools
# alignmentSieve --ATACshift, then coordinate sort and index.
# Called by TN5_SHIFT (modules/local/tn5_shift.nf):
#   tn5_shift.sh --sample S --bam S.final.bam --threads N
# Arguments:
#   --sample    sample id
#   --bam       coordinate-sorted final BAM (BAM_FINALIZE output)
#   --threads   alignmentSieve -p (task.cpus)
#   other arguments are appended to alignmentSieve.
# OFF by default (config tn5_shift: false). When switched on, the downstream steps (fragment sizes, TMM, bigWig,
# MACS3, FRiP) read S.tn5.bam instead of S.final.bam. How --ATACshift treats the proper-pair flags (the -f 3
# filter of BAM_FILTER is applied before the shift) is under statistician review; not yet verified.
# alignmentSieve writes in input order. The sort and index use pysam, because the deepTools image may have
# no samtools binary; pysam is a deepTools dependency.
# Outputs: S.tn5.bam, S.tn5.bam.bai (the unsorted intermediate S.tn5.unsorted.tmp.bam is removed)
# Exit 3 when the input BAM is missing or empty, or alignmentSieve leaves no reads.
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
[[ -n "${sample}" && -s "${bam}" ]] || { echo "tn5_shift.sh: missing --sample or empty --bam" >&2; exit 3; }

unsorted="${sample}.tn5.unsorted.tmp.bam"
alignmentSieve --ATACshift -b "${bam}" -o "${unsorted}" -p "${threads}" ${extra[@]+"${extra[@]}"}

python3 - "${unsorted}" "${sample}.tn5.bam" <<'PYEOF'
import sys

import pysam

tmp, out = sys.argv[1], sys.argv[2]
pysam.sort("-o", out, tmp)
pysam.index(out)
n = int(pysam.view("-c", out))
if n == 0:
    sys.stderr.write("tn5_shift.sh: alignmentSieve --ATACshift left no reads in %s\n" % tmp)
    sys.exit(3)
PYEOF

rm -f "${unsorted}"
