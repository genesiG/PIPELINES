#!/bin/bash
# bigwig.sh - one scaled, normalised bigWig per sample with deepTools bamCoverage.
# Called by BIGWIG (modules/local/bigwig.nf):
#   bigwig.sh --sample S --bam S.final.bam --multiplier M --method NAME --bin-size N \
#             --units rpkm|cpm|none --extend-reads true|false \
#             --effective-genome-size G --threads T [bamCoverage args...]
# Arguments:
#   --sample                 sample id (used verbatim in the output name)
#   --bam                    coordinate-sorted BAM, its .bai must sit next to it
#   --multiplier             bamCoverage --scaleFactor, from the scale-factor table
#                            (column bigwig_multiplier). Empty/NA/NaN/non-numeric/<=0 -> exit 3.
#   --method                 normalisation name; lower-cased into the file name
#   --bin-size               bin size in bp
#   --units                  rpkm | cpm | none  -> bamCoverage --normalizeUsing RPKM|CPM|None.
#                            ONE value for every sample, IgG included.
#   --extend-reads           true: --extendReads (paired-end fragments; legacy used
#                            --centerReads alone = read pileups, fixed here)
#   --effective-genome-size  passed to bamCoverage --effectiveGenomeSize when given
#   --threads                bamCoverage -p
# Any other arguments (config.yaml `bigwig_args`) are passed to bamCoverage.
# Output: <sample>.<method>.<units>.bw
#
# Legacy modes (step_5_generate_bigwig.py) -> units + multiplier. The multiplier now
# always comes from the scale-factor table (bigwig_multiplier):
#   SPIKEIN       bamCoverage --normalizeUsing RPKM --scaleFactor bamCov_SF  => units=rpkm, multiplier=bamCov_SF
#   SPIKEIN_CPM   --normalizeUsing CPM  --scaleFactor bamCov_SF              => units=cpm,  multiplier=bamCov_SF
#   SPIKEIN_RRPM  no depth normalisation, --scaleFactor 1e6/spike-in reads   => units=none, multiplier=1e6/spike reads
#   CHIPSEQSPIKEINFREE / TMM  RPKM with --scaleFactor 1/SF                   => units=rpkm, multiplier=1/SF
#   CPM / RPKM    built-in, scaleFactor 1                                    => units=cpm|rpkm, multiplier=1
# (legacy IgG always used CPM or RPKM with factor 1; now IgG gets the same units and its own multiplier.)
set -euo pipefail

sample="" bam="" mult="" method="" bin=10 units="rpkm" ext="true" egs="" threads=1 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)                sample="$2"; shift 2 ;;
        --bam)                   bam="$2"; shift 2 ;;
        --multiplier)            mult="$2"; shift 2 ;;
        --method)                method="$2"; shift 2 ;;
        --bin-size)              bin="$2"; shift 2 ;;
        --units)                 units="$2"; shift 2 ;;
        --extend-reads)          ext="$2"; shift 2 ;;
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
    *) echo "bigwig.sh: --units must be rpkm, cpm or none (got '${units}')" >&2; exit 2 ;;
esac
method_lc=$(echo "$method" | tr '[:upper:]' '[:lower:]')

args=(-b "$bam" -o "${sample}.${method_lc}.${units_lc}.bw" --outFileFormat bigwig
      --binSize "$bin" --skipNonCoveredRegions
      --normalizeUsing "$norm" --scaleFactor "$mult" -p "$threads")
if [[ "$(echo "$ext" | tr '[:upper:]' '[:lower:]')" == "true" ]]; then args+=(--extendReads); fi
if [[ -n "$egs" ]]; then args+=(--effectiveGenomeSize "$egs"); fi

echo "bamCoverage ${args[*]} ${extra[*]-}" >&2
bamCoverage "${args[@]}" ${extra[@]+"${extra[@]}"}
