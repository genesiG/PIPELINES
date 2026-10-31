#!/bin/bash
# heatmap_computematrix.sh - deepTools computeMatrix for one heatmap (reference-point, centre).
# Port of the computeMatrix step of the legacy heatmap scripts (reference-point, centre).
# Called by process COMPUTEMATRIX (also aliased COMPUTEMATRIX_CLUSTERED):
#   heatmap_computematrix.sh --out-prefix NAME --regions a.bed[,b.bed...] --samples id1,id2 \
#       --bigwigs f1.bw,f2.bw --window W --bin-size B --threads N [--region-labels l1,l2] [extra args]
# Arguments:
#   --out-prefix     NAME (the heatmap id); output NAME.matrix.gz
#   --regions        comma list of region BED files; one group per file, group names come from the file names
#   --samples        comma list of sample ids, one per bigWig; passed to computeMatrix --samplesLabel
#   --bigwigs        comma list of bigWig files, same order as --samples
#   --window         W, bp on each side of the reference point (computeMatrix -b W -a W)
#   --bin-size       B, bp per bin (computeMatrix -bs B)
#   --threads        N, computeMatrix -p N
#   --region-labels  accepted and ignored (computeMatrix names the groups from the file names)
#   --reference-point   center, TSS or TES, default center (computeMatrix --referencePoint)
#   --missing-as-zero   true or false, default true; false omits --missingDataAsZero
#   --average-type-bins X, default mean (computeMatrix --averageTypeBins)
#   any other args   passed on to computeMatrix, after the fixed arguments
# Runs: computeMatrix reference-point --referencePoint <reference-point> -b W -a W -bs B [--missingDataAsZero]
#       --averageTypeBins <average-type-bins> -p N -S <bigwigs> -R <regions> --samplesLabel <ids> -o NAME.matrix.gz
# Output: NAME.matrix.gz
# Exit 3: a regions file is missing or has no records (only comments or blank lines); a bigWig is missing;
#         the number of samples and bigWigs differ; W or B is not a positive integer; W % B != 0;
#         --reference-point is not center, TSS or TES; --missing-as-zero is neither true nor false.
# Exit 2: a required argument is missing.
set -euo pipefail

prefix="" regions="" samples="" bigwigs="" window="" bin="" threads=""
reference_point="center" missing_as_zero="true" average_type_bins="mean"
extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --out-prefix)    prefix="$2"; shift 2 ;;
        --regions)       regions="$2"; shift 2 ;;
        --samples)       samples="$2"; shift 2 ;;
        --bigwigs)       bigwigs="$2"; shift 2 ;;
        --window)        window="$2"; shift 2 ;;
        --bin-size)      bin="$2"; shift 2 ;;
        --threads)       threads="$2"; shift 2 ;;
        --region-labels) shift 2 ;;
        --reference-point)   reference_point="$2"; shift 2 ;;
        --missing-as-zero)   missing_as_zero="$2"; shift 2 ;;
        --average-type-bins) average_type_bins="$2"; shift 2 ;;
        *)               extra+=("$1"); shift ;;
    esac
done
[[ -n "$prefix" && -n "$regions" && -n "$samples" && -n "$bigwigs" && -n "$window" && -n "$bin" && -n "$threads" ]] \
    || { echo "heatmap_computematrix.sh: --out-prefix --regions --samples --bigwigs --window --bin-size --threads are required" >&2; exit 2; }

fail() { echo "heatmap_computematrix.sh: ERROR: $*" >&2; exit 3; }

[[ "$window" =~ ^[0-9]+$ ]] && (( window > 0 )) || fail "--window must be a positive integer (got '$window')"
[[ "$bin" =~ ^[0-9]+$ ]] && (( bin > 0 )) || fail "--bin-size must be a positive integer (got '$bin')"
(( window % bin == 0 )) || fail "window $window is not a multiple of bin size $bin"
case "$reference_point" in
    center|TSS|TES) ;;
    *) fail "--reference-point must be center, TSS or TES (got '$reference_point')" ;;
esac
case "$missing_as_zero" in
    true|false) ;;
    *) fail "--missing-as-zero must be true or false (got '$missing_as_zero')" ;;
esac

IFS=, read -ra REGION_FILES <<< "$regions"
for r in "${REGION_FILES[@]}"; do
    [[ -f "$r" ]] || fail "regions file not found: $r"
    n=$(awk '!/^#/ && NF > 0 {n++} END {print n + 0}' "$r")
    (( n > 0 )) || fail "regions file has no records (only comments or blank lines): $r"
done

IFS=, read -ra SAMPLE_IDS <<< "$samples"
IFS=, read -ra BIGWIGS <<< "$bigwigs"
(( ${#SAMPLE_IDS[@]} == ${#BIGWIGS[@]} )) \
    || fail "--samples has ${#SAMPLE_IDS[@]} ids but --bigwigs has ${#BIGWIGS[@]} files"
for b in "${BIGWIGS[@]}"; do
    [[ -f "$b" ]] || fail "bigWig not found: $b"
done

cmd=(computeMatrix reference-point --referencePoint "$reference_point"
     -b "$window" -a "$window" -bs "$bin"
     --averageTypeBins "$average_type_bins" -p "$threads"
     -S "${BIGWIGS[@]}" -R "${REGION_FILES[@]}"
     --samplesLabel "${SAMPLE_IDS[@]}"
     -o "${prefix}.matrix.gz")
if [[ "$missing_as_zero" == true ]]; then
    cmd+=(--missingDataAsZero)
fi
if (( ${#extra[@]} )); then
    cmd+=("${extra[@]}")
fi
"${cmd[@]}"
