#!/bin/bash
# seacr_consensus.sh - replicate consensus of SEACR peak BEDs for one group.
# Called by SEACR_CONSENSUS (modules/local/seacr_consensus.nf):
#   seacr_consensus.sh --group G --beds a.bed,b.bed,c.bed --min-replicates 2 --reciprocal-overlap 0.5
# Arguments:
#   --group               group name
#   --beds                the explicit replicate BEDs, comma or space separated (no globbing, so
#                         threshold-style outputs can never be picked up; legacy glob fixed)
#   --min-replicates      a peak is kept if it is supported by at least this many replicates
#                         (itself included). If the group has fewer BEDs, lowered to the number of
#                         BEDs with a warning.
#   --reciprocal-overlap  minimum reciprocal overlap fraction (bedtools intersect -f F -r) for a
#                         peak of another replicate to count as support
# Method: for every replicate BED, each peak counts the other replicates holding a peak with
# reciprocal overlap >= F; peaks with support >= min are kept; kept peaks of all replicates are
# concatenated, sorted and merged (bedtools merge).
# Output: <group>.seacr_consensus.bed (BED3, merged, sorted; may be empty)
set -euo pipefail

group="" beds="" minrep=2 ro=0.5
while [[ $# -gt 0 ]]; do
    case "$1" in
        --group)              group="$2"; shift 2 ;;
        --beds)               beds="$2"; shift 2 ;;
        --min-replicates)     minrep="$2"; shift 2 ;;
        --reciprocal-overlap) ro="$2"; shift 2 ;;
        *) echo "seacr_consensus.sh: unknown argument $1" >&2; exit 2 ;;
    esac
done
[[ -n "$group" && -n "$beds" ]] || { echo "seacr_consensus.sh: --group and --beds are required" >&2; exit 2; }
read -r -a files <<< "${beds//,/ }"
n=${#files[@]}
for f in "${files[@]}"; do [[ -e "$f" ]] || { echo "seacr_consensus.sh: missing input $f" >&2; exit 3; }; done
if ! [[ "$minrep" =~ ^[0-9]+$ ]] || [[ "$minrep" -lt 1 ]]; then echo "seacr_consensus.sh: --min-replicates must be an integer >= 1" >&2; exit 2; fi
if [[ "$minrep" -gt "$n" ]]; then
    echo "seacr_consensus.sh: WARNING group '${group}' has ${n} replicate(s) < min_replicates ${minrep}; using ${n}" >&2
    minrep=$n
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
: > "$tmp/kept.bed"
for ((i=0; i<n; i++)); do
    LC_ALL=C sort -k1,1 -k2,2n "${files[$i]}" | cut -f1-3 > "$tmp/rep$i.bed"
done
for ((i=0; i<n; i++)); do
    awk 'BEGIN{OFS="\t"}{print $0, 1}' "$tmp/rep$i.bed" > "$tmp/cnt.cur"
    for ((j=0; j<n; j++)); do
        [[ $i -eq $j ]] && continue
        bedtools intersect -a "$tmp/rep$i.bed" -b "$tmp/rep$j.bed" -f "$ro" -r -c > "$tmp/ov.tsv"
        paste "$tmp/cnt.cur" <(cut -f4 "$tmp/ov.tsv") | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4+($5>0?1:0)}' > "$tmp/cnt.new"
        mv "$tmp/cnt.new" "$tmp/cnt.cur"
    done
    awk -v m="$minrep" 'BEGIN{OFS="\t"} $4>=m {print $1,$2,$3}' "$tmp/cnt.cur" >> "$tmp/kept.bed"
done

LC_ALL=C sort -k1,1 -k2,2n "$tmp/kept.bed" | bedtools merge -i - > "${group}.seacr_consensus.bed"
echo "seacr_consensus.sh: ${group}: ${n} replicate(s), min ${minrep}, reciprocal overlap ${ro}, $(wc -l < "${group}.seacr_consensus.bed") consensus peaks" >&2
