#!/bin/bash
# bigwig_log2ratio.sh - log2(target / control) bigWig from two step-BIGWIG tracks
# (deepTools bigwigCompare --operation log2; port of CBX2 step_5c).
# Called by BIGWIG_LOG2RATIO (modules/local/bigwig_log2ratio.nf):
#   bigwig_log2ratio.sh --sample S --target S.x.bw --control C.x.bw --control-id C \
#                       --bin-size N --threads T [--pseudocount auto|NUMBER]
# Arguments:
#   --sample       target sample id
#   --target       target bigWig
#   --control      control bigWig
#   --control-id   control sample id (output name)
#   --bin-size     output bin in bp (bigwigCompare cannot resolve below the source bins)
#   --threads      bigwigCompare -p
#   --pseudocount  'auto' (default) or a number. auto = 5 x the smallest non-zero bin value of
#                  the noisier of the two tracks (= about 5 reads per bin; in the legacy script this
#                  was computed from read depth, here it is estimated from the bigWigs themselves
#                  with pyBigWig over the first ~2e6 intervals); falls back to 1 with a warning.
# Any other arguments are passed to bigwigCompare. Bins empty in both files are skipped.
# Output: <sample>.vs_<control-id>.log2ratio.bw
set -euo pipefail

sample="" target="" control="" cid="" bin=1000 threads=1 pc="auto" extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)      sample="$2"; shift 2 ;;
        --target)      target="$2"; shift 2 ;;
        --control)     control="$2"; shift 2 ;;
        --control-id)  cid="$2"; shift 2 ;;
        --bin-size)    bin="$2"; shift 2 ;;
        --threads)     threads="$2"; shift 2 ;;
        --pseudocount) pc="$2"; shift 2 ;;
        *)             extra+=("$1"); shift ;;
    esac
done
[[ -n "$sample" && -n "$target" && -n "$control" && -n "$cid" ]] || { echo "bigwig_log2ratio.sh: --sample --target --control --control-id are required" >&2; exit 2; }

if [[ "$pc" == "auto" ]]; then
    pc=$(python3 - "$target" "$control" <<'PY' || true
import sys
try:
    import pyBigWig
except Exception:
    sys.exit(0)
vals = []
for path in sys.argv[1:3]:
    bw = pyBigWig.open(path)
    lo, seen = None, 0
    for chrom in bw.chroms():
        ivs = bw.intervals(chrom) or []
        for _, _, v in ivs:
            if v > 0 and (lo is None or v < lo):
                lo = v
        seen += len(ivs)
        if seen > 2000000:
            break
    if lo is not None:
        vals.append(lo)
if vals:
    print("%.6g" % (5 * max(vals)))
PY
)
    if [[ -z "$pc" ]]; then
        echo "bigwig_log2ratio.sh: WARNING could not estimate pseudocount, using 1" >&2
        pc=1
    fi
fi
echo "pseudocount: ${pc}" >&2

bigwigCompare -b1 "$target" -b2 "$control" --operation log2 --pseudocount "$pc" \
    --binSize "$bin" --skipZeroOverZero -p "$threads" --outFileFormat bigwig \
    -o "${sample}.vs_${cid}.log2ratio.bw" ${extra[@]+"${extra[@]}"}
