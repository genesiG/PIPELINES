#!/bin/bash
# callpeaks_macs3.sh - MACS3 callpeak for one group (replicates POOLED, as in legacy step_6a).
# Called by MACS3 (modules/local/macs3.nf):
#   callpeaks_macs3.sh --group G --targets a.bam,b.bam --controls c.bam,d.bam \
#       --genome-size N --qvalue Q --broad true|false --broad-cutoff C
# Arguments:
#   --group         group name; used as macs3 -n
#   --targets       treatment BAMs, comma or space separated
#   --controls      control BAMs (samplesheet controls of the same condition), comma or space
#                   separated. Empty -> exit 3 (legacy glob never matched and MACS ran control-free).
#   --genome-size   effective genome size (-g)
#   --qvalue        -q
#   --broad         true: --broad --broad-cutoff C
#   --broad-cutoff  broad cutoff
# Any other arguments (config.yaml `callpeaks_macs3_args`) are passed to macs3 callpeak.
# Fixed: -f BAMPE --keep-dup all.
# Outputs: <group>_peaks.narrowPeak | <group>_peaks.broadPeak, <group>_peaks.xls,
#          (narrow: <group>_summits.bed; broad: <group>_peaks.gappedPeak), <group>.macs3.log
set -euo pipefail

group="" targets="" controls="" gs="" q="0.01" broad="false" bc="0.1" extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --group)        group="$2"; shift 2 ;;
        --targets)      targets="$2"; shift 2 ;;
        --controls)     controls="$2"; shift 2 ;;
        --genome-size)  gs="$2"; shift 2 ;;
        --qvalue)       q="$2"; shift 2 ;;
        --broad)        broad="$2"; shift 2 ;;
        --broad-cutoff) bc="$2"; shift 2 ;;
        *)              extra+=("$1"); shift ;;
    esac
done
[[ -n "$group" && -n "$gs" ]] || { echo "callpeaks_macs3.sh: --group and --genome-size are required" >&2; exit 2; }

read -r -a tbams <<< "${targets//,/ }"
read -r -a cbams <<< "${controls//,/ }"
if [[ ${#tbams[@]} -eq 0 ]]; then
    echo "callpeaks_macs3.sh: group '${group}' has no target BAMs" >&2; exit 3
fi
if [[ ${#cbams[@]} -eq 0 ]]; then
    echo "callpeaks_macs3.sh: group '${group}' has no control BAMs (samplesheet control_id); MACS3 is not run without a control" >&2
    exit 3
fi

args=(callpeak -f BAMPE -g "$gs" --keep-dup all -q "$q" -n "$group" --outdir .)
if [[ "$(echo "$broad" | tr '[:upper:]' '[:lower:]')" == "true" ]]; then args+=(--broad --broad-cutoff "$bc"); fi
args+=(-t "${tbams[@]}" -c "${cbams[@]}")

echo "macs3 ${args[*]} ${extra[*]-}" > "${group}.macs3.log"
macs3 "${args[@]}" ${extra[@]+"${extra[@]}"} 2>&1 | tee -a "${group}.macs3.log"
