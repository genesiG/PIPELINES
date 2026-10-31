#!/bin/bash
# bigwig.sh - one scaled bigWig per sample with deepTools bamCoverage (legacy step_5_generate_bigwig.py settings).
# Called by BIGWIG (modules/local/bigwig.nf):
#   bigwig.sh --sample S --bam S.final.bam --multiplier M --method NAME --bin-size 10 \
#             --units none|cpm|rpkm --extend-reads true|false --center-reads true|false \
#             --effective-genome-size G --threads T [bamCoverage args...]
# Arguments:
#   --sample                 sample id (used verbatim in the output name)
#   --bam                    coordinate-sorted BAM (S.final.bam, or S.tn5.bam when tn5_shift is on); its .bai must sit next to it
#   --multiplier             bamCoverage --scaleFactor, from the scale-factor table (column bigwig_multiplier).
#                            Empty/NA/NaN/non-numeric/<=0 -> exit 3 (never defaulted to 1).
#   --method                 normalisation name (e.g. tmm); lower-cased into the file name
#   --bin-size               bamCoverage --binSize, bp (default 10, legacy)
#   --units                  none | cpm | rpkm -> bamCoverage --normalizeUsing None|CPM|RPKM (default none, legacy;
#                            ONE value for every sample). Scaling then comes entirely from --multiplier.
#   --extend-reads           true: --extendReads (config default true: paired-end fragment coverage, as CUT&RUN).
#                            false = legacy: each mate counted over its own aligned blocks.
#   --center-reads           true: --centerReads (config default false). Legacy passed it without --extendReads,
#                            where deepTools 3.5.6 ignores it for paired-end reads (countReadsPerBin.get_fragment_from_read
#                            returns the aligned blocks before the centring code). With --extendReads it replaces each
#                            fragment by a read-length window at the fragment midpoint.
#   --effective-genome-size  passed to bamCoverage --effectiveGenomeSize when given
#   --threads                bamCoverage -p
# Fixed: --outFileFormat bigwig --skipNonCoveredRegions (legacy).
# Any other arguments (config.yaml `bigwig_args`) are passed to bamCoverage.
# Output: <sample>.<method>.<units>.bw, method and units lower-cased (e.g. S.tmm.none.bw)
#
# Legacy (step_5_generate_bigwig.py, AV_ATACSEQ_2026): units None, --binSize 10 --skipNonCoveredRegions
# --centerReads, --scaleFactor = 1/SF (SF = nf * lib.size / 1e6, computed in norm_tmm.R as bigwig_multiplier
# when bigwig_units is none: 1e6 / (nf * lib.size)).
set -euo pipefail

sample="" bam="" mult="" method="" bin=10 units="none" ext="false" center="true" egs="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)                sample="$2"; shift 2 ;;
        --bam)                   bam="$2"; shift 2 ;;
        --multiplier)            mult="$2"; shift 2 ;;
        --method)                method="$2"; shift 2 ;;
        --bin-size)              bin="$2"; shift 2 ;;
        --units)                 units="$2"; shift 2 ;;
        --extend-reads)          ext="$2"; shift 2 ;;
        --center-reads)          center="$2"; shift 2 ;;
        --effective-genome-size) egs="$2"; shift 2 ;;
        --threads)               threads="$2"; shift 2 ;;
        *)                       extra+=("$1"); shift ;;
    esac
done

[[ -n "$sample" && -n "$bam" && -n "$method" ]] || { echo "bigwig.sh: --sample, --bam and --method are required" >&2; exit 2; }

# multiplier must be a finite number > 0 (legacy silently used 1)
if ! [[ "$mult" =~ ^\+?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$ ]] || ! awk -v m="$mult" 'BEGIN{exit !(m+0 > 0)}'; then
    echo "bigwig.sh: sample '${sample}' has missing/invalid bigwig multiplier '${mult}' (empty, NA, NaN, non-numeric or <= 0); check the scale-factor table" >&2
    exit 3
fi

units_lc=$(echo "$units" | tr '[:upper:]' '[:lower:]')
case "$units_lc" in
    rpkm) norm="RPKM" ;;
    cpm)  norm="CPM" ;;
    none) norm="None" ;;
    *) echo "bigwig.sh: --units must be none, cpm or rpkm (got '${units}')" >&2; exit 2 ;;
esac
method_lc=$(echo "$method" | tr '[:upper:]' '[:lower:]')

args=(-b "$bam" -o "${sample}.${method_lc}.${units_lc}.bw" --outFileFormat bigwig
      --binSize "$bin" --skipNonCoveredRegions
      --normalizeUsing "$norm" --scaleFactor "$mult" -p "$threads")
if [[ "$(echo "$ext" | tr '[:upper:]' '[:lower:]')" == "true" ]]; then args+=(--extendReads); fi
if [[ "$(echo "$center" | tr '[:upper:]' '[:lower:]')" == "true" ]]; then args+=(--centerReads); fi
if [[ -n "$egs" ]]; then args+=(--effectiveGenomeSize "$egs"); fi

echo "bamCoverage ${args[*]} ${extra[*]-}" >&2
bamCoverage "${args[@]}" ${extra[@]+"${extra[@]}"}
