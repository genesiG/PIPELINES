#!/bin/bash
# callpeaks_macs3.sh - MACS3 callpeak for one peak group: replicates POOLED, NO control (legacy step_4c_macs2.py).
# Called by MACS3 (modules/local/macs3.nf):
#   callpeaks_macs3.sh --group G --prefix G.qval0.05.minlen250.maxgap500.fe2.0 --targets a.bam,b.bam \
#       --genome-size N --qvalue Q --min-length L --max-gap G --fe-cutoff F --call-summits true|false \
#       --broad true|false --broad-cutoff C
# Arguments:
#   --group         peak group name (messages only)
#   --prefix        output prefix = macs3 -n. Computed by the module: {group}.qval{q}[.minlen{l}][.maxgap{g}][.fe{fe}][.broad]
#   --targets       treatment BAMs of the group, comma or space separated. Every BAM must be non-empty (exit 3).
#   --genome-size   effective genome size (macs3 -g)
#   --qvalue        macs3 -q
#   --min-length    --min-length; empty = option left out (legacy: off when null)
#   --max-gap       --max-gap; empty = option left out (legacy: off when null)
#   --fe-cutoff     --fe-cutoff, only when the value is > 1 (legacy: option left out when <= 1)
#   --call-summits  true: --call-summits (writes <prefix>_summits.bed, narrow peaks). Legacy default false.
#   --broad         true: --broad (writes <prefix>_peaks.broadPeak and <prefix>_peaks.gappedPeak). Legacy default false.
#   --broad-cutoff  --broad-cutoff, used only with --broad; empty = macs3 default
# Any other arguments are passed to macs3 callpeak. They come from two places, both placed after the options above:
#   `macs3_args` (nextflow.config, put on the command line by modules/local/macs3.nf) and
#   `callpeaks_macs3_args` (config.yaml, via extraArgs() in the same module).
# Fixed: -f BAMPE --keep-dup all --outdir . ; no -c (no control). The legacy script used the same options.
# The --fe-cutoff option must be checked with `macs3 callpeak --help` in the container.
# Outputs: <prefix>_peaks.narrowPeak | <prefix>_peaks.broadPeak, <prefix>_peaks.xls,
#   <prefix>_summits.bed (only with --call-summits, narrow peaks), <prefix>_peaks.gappedPeak (only with --broad),
#   <prefix>.macs3.log (the macs3 command line, then macs3 stdout/stderr)
# Exit 3: no target BAM, or a target BAM missing or empty. Exit 2: missing required argument.
set -euo pipefail

group="" prefix="" targets="" gs="" q="0.05" minlen="" maxgap="" fe="" summits="false" broad="false" bc="" extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --group)        group="$2"; shift 2 ;;
        --prefix)       prefix="$2"; shift 2 ;;
        --targets)      targets="$2"; shift 2 ;;
        --genome-size)  gs="$2"; shift 2 ;;
        --qvalue)       q="$2"; shift 2 ;;
        --min-length)   minlen="$2"; shift 2 ;;
        --max-gap)      maxgap="$2"; shift 2 ;;
        --fe-cutoff)    fe="$2"; shift 2 ;;
        --call-summits) summits="$2"; shift 2 ;;
        --broad)        broad="$2"; shift 2 ;;
        --broad-cutoff) bc="$2"; shift 2 ;;
        *)              extra+=("$1"); shift ;;
    esac
done
[[ -n "$group" && -n "$prefix" && -n "$gs" && -n "$q" ]] || { echo "callpeaks_macs3.sh: --group, --prefix, --genome-size and --qvalue are required" >&2; exit 2; }

read -r -a tbams <<< "${targets//,/ }"
if [[ ${#tbams[@]} -eq 0 ]]; then
    echo "callpeaks_macs3.sh: peak group '${group}' has no target BAMs" >&2; exit 3
fi
for b in "${tbams[@]}"; do
    [[ -s "$b" ]] || { echo "callpeaks_macs3.sh: peak group '${group}': target BAM '${b}' is missing or empty" >&2; exit 3; }
done

args=(callpeak -f BAMPE -g "$gs" --keep-dup all -q "$q" -n "$prefix" --outdir .)
if [[ -n "$minlen" ]]; then args+=(--min-length "$minlen"); fi
if [[ -n "$maxgap" ]]; then args+=(--max-gap "$maxgap"); fi
if [[ -n "$fe" ]] && awk -v f="$fe" 'BEGIN{exit !(f+0 > 1)}'; then args+=(--fe-cutoff "$fe"); fi
if [[ "$(echo "$summits" | tr '[:upper:]' '[:lower:]')" == "true" ]]; then args+=(--call-summits); fi
if [[ "$(echo "$broad" | tr '[:upper:]' '[:lower:]')" == "true" ]]; then
    args+=(--broad)
    if [[ -n "$bc" ]]; then args+=(--broad-cutoff "$bc"); fi
fi
args+=(-t "${tbams[@]}")

echo "macs3 ${args[*]} ${extra[*]-}" > "${prefix}.macs3.log"
macs3 "${args[@]}" ${extra[@]+"${extra[@]}"} 2>&1 | tee -a "${prefix}.macs3.log"
