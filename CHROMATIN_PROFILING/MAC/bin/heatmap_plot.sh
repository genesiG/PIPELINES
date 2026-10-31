#!/bin/bash
# heatmap_plot.sh - deepTools plotHeatmap for one matrix, no clustering.
# Port of the plotHeatmap step of the legacy heatmap scripts; label cleanup as in the legacy EED version.
# Called by process PLOTHEATMAP (cache 'deep'):
#   heatmap_plot.sh --matrix M --out-prefix NAME --labels "a,b" --colors "#hex,#hex" \
#       [--z-min "0,0"] [--z-max "25,25"] [--y-max "25,25"] [--region-labels "A,B"] \
#       [--sort-using INDEX] [--ref-label "peak center"] [--x-label "distance (bp)"] \
#       [--y-label "signal"] [--height 12] [--width 3.5] [--threads N] [extra plotHeatmap args]
# Arguments:
#   --matrix         deepTools matrix (.gz) from COMPUTEMATRIX
#   --out-prefix     NAME (the heatmap id)
#   --labels         sample labels, comma separated; one per sample (n = number of labels).
#                    One pair of surrounding quotes is stripped; the two characters backslash+n become a newline.
#   --colors         one hex colour per sample; each gives --colorList 'white,<hue>'.
#                    A single colour is used for every sample.
#   --z-min          optional; one value or one per sample (--zMin). Empty: if --z-max is given, zero per sample.
#   --z-max          optional; one value or one per sample (--zMax)
#   --y-max          optional; one value or one per sample (--yMax). Independent of --z-max.
#   --region-labels  optional, comma separated, one per region group (--regionsLabel); same label cleanup
#   --sort-using     optional 1-based sample index (--sortUsingSamples)
#   --ref-label      reference point label (--refPointLabel), default 'peak center'
#   --x-label        x axis label (-x), default 'distance (bp)'
#   --y-label        y axis label (-y), default 'signal'
#   --height         heatmap height in cm (--heatmapHeight), default 12
#   --width          heatmap width in cm (--heatmapWidth), default 3.5
#   --plot-type      plotHeatmap --plotType, default se
#   --average-type-summary  plotHeatmap --averageTypeSummaryPlot, default mean
#   --legend-location       plotHeatmap --legendLocation, default none
#   --threads        optional; accepted for interface symmetry; plotHeatmap is single-threaded
#   any other args   passed on to plotHeatmap
# Runs: plotHeatmap -m M -o NAME.heatmap.svg --outFileSortedRegions NAME.sorted_regions.bed --plotType <plot-type>
#       --averageTypeSummaryPlot <average-type-summary> --whatToShow 'plot, heatmap and colorbar' --heatmapHeight H
#       --heatmapWidth W --colorList white,<hue> ... --samplesLabel ... [--zMin] [--zMax] [--yMax] [--regionsLabel]
#       [--sortUsingSamples] --legendLocation <legend-location> --plotFileFormat svg -x X -y Y --refPointLabel R
# This script never clusters (no --kmeans, --hclust or --silhouette).
# Outputs: NAME.heatmap.svg NAME.sorted_regions.bed
# Exit 3: matrix missing or empty; a list has neither 1 nor n entries (n = number of labels); sort index out of range.
# Exit 2: a required argument is missing.
set -euo pipefail

fail() { echo "heatmap_plot.sh: ERROR: $*" >&2; exit 3; }

matrix="" name="" labels="" colors="" zmin="" zmax="" ymax="" rlabels="" sort_using=""
ref_label="peak center" xlab="distance (bp)" ylab="signal" height=12 width=3.5 threads=""
plot_type="se" average_type_summary="mean" legend_location="none"
extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --matrix)        matrix="$2"; shift 2 ;;
        --out-prefix)    name="$2"; shift 2 ;;
        --labels)        labels="$2"; shift 2 ;;
        --colors)        colors="$2"; shift 2 ;;
        --z-min)         zmin="$2"; shift 2 ;;
        --z-max)         zmax="$2"; shift 2 ;;
        --y-max)         ymax="$2"; shift 2 ;;
        --region-labels) rlabels="$2"; shift 2 ;;
        --sort-using)    sort_using="$2"; shift 2 ;;
        --ref-label)     ref_label="$2"; shift 2 ;;
        --x-label)       xlab="$2"; shift 2 ;;
        --y-label)       ylab="$2"; shift 2 ;;
        --height)        height="$2"; shift 2 ;;
        --width)         width="$2"; shift 2 ;;
        --threads)       threads="$2"; shift 2 ;;
        --plot-type)             plot_type="$2"; shift 2 ;;
        --average-type-summary)  average_type_summary="$2"; shift 2 ;;
        --legend-location)       legend_location="$2"; shift 2 ;;
        *)               extra+=("$1"); shift ;;
    esac
done
[[ -n "$matrix" && -n "$name" && -n "$labels" && -n "$colors" ]] \
    || { echo "heatmap_plot.sh: --matrix --out-prefix --labels --colors are required" >&2; exit 2; }
[[ -s "$matrix" ]] || fail "matrix missing or empty: $matrix"

# Expand a per-sample list: one entry is repeated for every sample; otherwise it needs n entries.
# Sets EXPANDED.
expand_list() { # $1 option name, $2 n samples, $3 comma list
    local -a vals=()
    local i
    IFS=, read -ra vals <<< "$3"
    EXPANDED=()
    if (( ${#vals[@]} == 1 )); then
        for ((i = 0; i < $2; i++)); do EXPANDED+=("${vals[0]}"); done
    elif (( ${#vals[@]} == $2 )); then
        EXPANDED=("${vals[@]}")
    else
        fail "$1 has ${#vals[@]} entries; expected 1 or $2 (one per sample)"
    fi
}

# Label cleanup in place: strip one pair of surrounding quotes, then literal backslash+n -> newline.
clean_labels() { # $1 array name
    local -n arr="$1"
    local i l q
    for i in "${!arr[@]}"; do
        l="${arr[$i]}"
        q="${l:0:1}"
        if (( ${#l} >= 2 )) && [[ "$q" == "'" || "$q" == '"' ]] && [[ "${l: -1}" == "$q" ]]; then
            l="${l:1:${#l}-2}"
        fi
        l=${l//'\n'/$'\n'}
        arr[$i]="$l"
    done
}

IFS=, read -ra LABELS <<< "$labels"
n=${#LABELS[@]}
(( n > 0 )) || fail "--labels is empty"
clean_labels LABELS
# a literal \n in --ref-label becomes a line break, as in the labels
ref_label=${ref_label//'\n'/$'\n'}

expand_list --colors "$n" "$colors"
COLORS=("${EXPANDED[@]}")

ZMIN=() ZMAX=() YMAX=() RLABELS=()
if [[ -n "$zmin" ]]; then
    expand_list --z-min "$n" "$zmin"; ZMIN=("${EXPANDED[@]}")
elif [[ -n "$zmax" ]]; then
    for ((i = 0; i < n; i++)); do ZMIN+=(0); done
fi
if [[ -n "$zmax" ]]; then
    expand_list --z-max "$n" "$zmax"; ZMAX=("${EXPANDED[@]}")
fi
if [[ -n "$ymax" ]]; then
    expand_list --y-max "$n" "$ymax"; YMAX=("${EXPANDED[@]}")
fi
if [[ -n "$rlabels" ]]; then
    IFS=, read -ra RLABELS <<< "$rlabels"
    clean_labels RLABELS
fi
if [[ -n "$sort_using" ]]; then
    [[ "$sort_using" =~ ^[0-9]+$ ]] && (( sort_using >= 1 && sort_using <= n )) \
        || fail "--sort-using must be a sample index 1..$n (got '$sort_using')"
fi

cmd=(plotHeatmap -m "$matrix" -o "${name}.heatmap.svg"
     --outFileSortedRegions "${name}.sorted_regions.bed"
     --plotType "$plot_type" --averageTypeSummaryPlot "$average_type_summary"
     --whatToShow 'plot, heatmap and colorbar'
     --heatmapHeight "$height" --heatmapWidth "$width"
     --colorList)
for c in "${COLORS[@]}"; do
    cmd+=("white,${c}")
done
cmd+=(--samplesLabel "${LABELS[@]}")
if (( ${#ZMIN[@]} )); then cmd+=(--zMin "${ZMIN[@]}"); fi
if (( ${#ZMAX[@]} )); then cmd+=(--zMax "${ZMAX[@]}"); fi
if (( ${#YMAX[@]} )); then cmd+=(--yMax "${YMAX[@]}"); fi
if (( ${#RLABELS[@]} )); then cmd+=(--regionsLabel "${RLABELS[@]}"); fi
if [[ -n "$sort_using" ]]; then cmd+=(--sortUsingSamples "$sort_using"); fi
cmd+=(--legendLocation "$legend_location" --plotFileFormat svg
      -x "$xlab" -y "$ylab" --refPointLabel "$ref_label")
if (( ${#extra[@]} )); then
    cmd+=("${extra[@]}")
fi
"${cmd[@]}"
