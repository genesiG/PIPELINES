#!/bin/bash
# seacr_call.sh - SEACR peak calling of one target against its matched control bedGraph.
# Called by SEACR_CALL (modules/local/seacr_call.nf):
#   seacr_call.sh --sample S --target-bg S.fragments.bedgraph --control-bg C.fragments.bedgraph \
#                 --mode relaxed|stringent --norm non|norm
# Arguments:
#   --sample      sample id
#   --target-bg   target fragments bedGraph (sorted)
#   --control-bg  control (IgG) bedGraph; required, exit 3 if empty (every target needs a control)
#   --mode        relaxed | stringent
#   --norm        norm (default; legacy step_6b) | non
# Legacy invocation: SEACR_1.3.sh target.bg control.bg norm relaxed prefix
# Output: <sample>.seacr.<mode>.bed (SEACR's <prefix>.<mode>.bed renamed; if SEACR writes an empty
# file a header-less empty file is kept and a warning printed).
set -euo pipefail

sample="" tbg="" cbg="" mode="relaxed" norm="norm" extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample)     sample="$2"; shift 2 ;;
        --target-bg)  tbg="$2"; shift 2 ;;
        --control-bg) cbg="$2"; shift 2 ;;
        --mode)       mode="$2"; shift 2 ;;
        --norm)       norm="$2"; shift 2 ;;
        *)            extra+=("$1"); shift ;;
    esac
done
[[ -n "$sample" && -n "$tbg" ]] || { echo "seacr_call.sh: --sample and --target-bg are required" >&2; exit 2; }
if [[ -z "$cbg" ]]; then
    echo "seacr_call.sh: sample '${sample}' has no control bedGraph (control_id in samplesheet)" >&2; exit 3
fi
case "$mode" in relaxed|stringent) ;; *) echo "seacr_call.sh: --mode must be relaxed or stringent" >&2; exit 2 ;; esac
case "$norm" in non|norm) ;; *) echo "seacr_call.sh: --norm must be non or norm" >&2; exit 2 ;; esac

seacr=$(command -v SEACR_1.3.sh || command -v SEACR || true)
[[ -n "$seacr" ]] || { echo "seacr_call.sh: SEACR_1.3.sh / SEACR not found in PATH" >&2; exit 2; }

bash "$seacr" "$tbg" "$cbg" "$norm" "$mode" "${sample}" ${extra[@]+"${extra[@]}"}

out="${sample}.${mode}.bed"
[[ -e "$out" ]] || { echo "seacr_call.sh: SEACR did not write ${out}" >&2; exit 3; }
mv "$out" "${sample}.seacr.${mode}.bed"
if [[ ! -s "${sample}.seacr.${mode}.bed" ]]; then
    echo "seacr_call.sh: WARNING SEACR output for ${sample} is empty" >&2
fi
