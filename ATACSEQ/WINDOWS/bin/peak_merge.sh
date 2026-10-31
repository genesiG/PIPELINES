#!/bin/bash
# peak_merge.sh - union of several peak sets into one BED with bedops --merge (legacy merge_peaks.sh).
# Called by PEAK_MERGE (modules/local/peak_merge.nf):
#   peak_merge.sh --name high_gd2 --peaks a.narrowPeak,b.narrowPeak,c.narrowPeak
# Arguments:
#   --name    name of the merged set; the output is <name>.merged.bed
#   --peaks   the peak files of the set, comma or space separated; they must sit in the task directory
#             (the module stages them there under their own names).
#             Input order does not matter: bedops --merge needs sort-bed order, and MACS3 writes peaks in
#             BAM-header contig order (not sort-bed order), so the script sorts all listed files together
#             with sort-bed into a temp file before the merge.
# Output: <name>.merged.bed (bedops --merge output: chrom, start, end; one row per merged interval)
# Exit 3: a listed file is missing, or the list of files is empty.
#         A listed file of size 0 (no peaks in that group) is merged as empty input, with a warning.
# Exit 2: missing --name.
set -euo pipefail

name="" peaks=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --name)  name="$2"; shift 2 ;;
        --peaks) peaks="$2"; shift 2 ;;
        *)       echo "peak_merge.sh: unknown argument $1" >&2; exit 2 ;;
    esac
done
[[ -n "$name" ]] || { echo "peak_merge.sh: --name is required" >&2; exit 2; }

read -r -a files <<< "${peaks//,/ }"
if [[ ${#files[@]} -eq 0 ]]; then
    echo "peak_merge.sh: no peak files listed for merged set '${name}'" >&2; exit 3
fi
for f in "${files[@]}"; do
    [[ -e "$f" ]] || { echo "peak_merge.sh: peak file '${f}' of merged set '${name}' is missing" >&2; exit 3; }
    [[ -s "$f" ]] || echo "peak_merge.sh: WARNING peak file '${f}' of merged set '${name}' is empty" >&2
done

# sort-bed reads all listed files and writes one sorted stream. A temp file is used instead of
# process substitution, which does not behave reliably inside the Nextflow task shell.
sort-bed "${files[@]}" > sorted.tmp.bed
bedops --merge sorted.tmp.bed > "${name}.merged.bed"
rm -f sorted.tmp.bed
