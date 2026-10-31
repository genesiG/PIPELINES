#!/bin/bash
# tracks_pygenometracks.sh - render one locus of a tracks.ini with pyGenomeTracks.
# Called by process PYGENOMETRACKS (one call per locus line of NAME.regions.tsv):
#   tracks_pygenometracks.sh --ini F.ini --region chr:s-e --out-prefix SLUG \
#       --width 15 --dpi 600 --format svg [extra pyGenomeTracks args]
# Arguments:
#   --ini         tracks.ini from TRACKS_INI (NAME.tracks.ini)
#   --region      chr:start-end to draw
#   --out-prefix  SLUG (the locus slug from NAME.regions.tsv, e.g. NAME.GENEA_125kb)
#   --width       figure width in cm (pyGenomeTracks --width), default 15
#   --dpi         pyGenomeTracks --dpi, default 600
#   --format      output format, default svg
#   any other args passed on to pyGenomeTracks
# Runs: pyGenomeTracks --tracks <ini> --region <region> --outFileName SLUG.<format> --width W --dpi D
# Output: SLUG.<format> (default SLUG.svg)
# Exit 3: the ini file is missing.
# Exit 2: a required argument is missing.
set -euo pipefail

ini="" region="" slug="" width=15 dpi=600 format=svg
extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --ini)        ini="$2"; shift 2 ;;
        --region)     region="$2"; shift 2 ;;
        --out-prefix) slug="$2"; shift 2 ;;
        --width)      width="$2"; shift 2 ;;
        --dpi)        dpi="$2"; shift 2 ;;
        --format)     format="$2"; shift 2 ;;
        *)            extra+=("$1"); shift ;;
    esac
done
[[ -n "$ini" && -n "$region" && -n "$slug" ]] \
    || { echo "tracks_pygenometracks.sh: --ini --region --out-prefix are required" >&2; exit 2; }
[[ -f "$ini" ]] || { echo "tracks_pygenometracks.sh: ERROR: tracks ini not found: $ini" >&2; exit 3; }

cmd=(pyGenomeTracks --tracks "$ini" --region "$region"
     --outFileName "${slug}.${format}" --width "$width" --dpi "$dpi")
if (( ${#extra[@]} )); then
    cmd+=("${extra[@]}")
fi
"${cmd[@]}"
