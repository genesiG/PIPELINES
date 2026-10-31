#!/bin/bash
# mageck_test.sh — `mageck test` for one treatment-vs-control contrast.
# Called by MAGECK_TEST (modules/local/mageck_test.nf):
#   mageck_test.sh --count-table T.count.txt --name TREAT_vs_CTRL \
#                  --treatment "t1,t2" --control "c1,c2" [mageck test args...]
# --treatment/--control are count-table column names (validate_contrasts.py resolves them).
# Any other arguments (config.yaml `mageck_test_args`) are passed to mageck test.
# Must write: NAME* (gene_summary, sgrna_summary, ...)
set -euo pipefail

count_table="" name="" treatment="" control="" extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --count-table) count_table="$2"; shift 2 ;;
        --name)        name="$2"; shift 2 ;;
        --treatment)   treatment="$2"; shift 2 ;;
        --control)     control="$2"; shift 2 ;;
        *)             extra+=("$1"); shift ;;
    esac
done

mageck test \
    -k "${count_table}" \
    -n "${name}" \
    --norm-method median \
    --adjust-method fdr \
    -c "${control}" \
    -t "${treatment}" \
    ${extra[@]+"${extra[@]}"}
