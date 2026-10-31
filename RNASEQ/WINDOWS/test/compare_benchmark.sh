#!/bin/bash
# compare_benchmark.sh — checks a finished run against the RNA-seq benchmark.
#
#   bash test/compare_benchmark.sh PROJECT_DIR BENCHMARK_DIR
#     PROJECT_DIR    the project folder of the run to check (its results/ folder is compared)
#     BENCHMARK_DIR  the benchmark: benchmark/sub or benchmark/full (read-only; e.g. /d/PIPELINES/RNASEQ/benchmark/sub)
#
# (a) every results/counts/*_counts.txt of BENCHMARK_DIR must exist in PROJECT_DIR with the same bytes (cmp);
# (b) every line of BENCHMARK_DIR/benchmark_text_md5.txt must check OK with md5sum -c, run from PROJECT_DIR
#     (the list names the text results as results/...). A line that is missing or different counts as failed.
# Prints one line per counts table and a one-line PASS/FAIL summary with the number of failed lines.
# Exit 0 = PASS, 1 = FAIL, 2 = usage error. Nothing is written to either folder.
set -uo pipefail

if [[ $# -ne 2 ]]; then
    echo "usage: bash test/compare_benchmark.sh PROJECT_DIR BENCHMARK_DIR" >&2
    exit 2
fi
project="$(cd "$1" 2>/dev/null && pwd)" || { echo "not a folder: $1" >&2; exit 2; }
bench="$(cd "$2" 2>/dev/null && pwd)" || { echo "not a folder: $2" >&2; exit 2; }
md5_list="${bench}/benchmark_text_md5.txt"
[[ -f "${md5_list}" ]] || { echo "missing ${md5_list}" >&2; exit 2; }

# ---- (a) counts tables -----------------------------------------------------------
n_counts=0
n_counts_bad=0
for ref in "${bench}"/results/counts/*_counts.txt; do
    [[ -f "${ref}" ]] || continue
    n_counts=$(( n_counts + 1 ))
    name="$(basename "${ref}")"
    if cmp -s "${ref}" "${project}/results/counts/${name}"; then
        echo "counts  OK      ${name}"
    else
        echo "counts  FAILED  ${name} (missing in the project or different)"
        n_counts_bad=$(( n_counts_bad + 1 ))
    fi
done
if (( n_counts == 0 )); then
    echo "counts  FAILED  no results/counts/*_counts.txt in ${bench}"
    n_counts_bad=1
fi

# ---- (b) text results, md5 ---------------------------------------------------------
md5_total=$(grep -c . "${md5_list}")
md5_out="$(cd "${project}" && md5sum -c "${md5_list}" 2>&1)"
md5_ok=$(printf '%s\n' "${md5_out}" | grep -c ': OK$')
md5_failed=$(( md5_total - md5_ok ))
printf '%s\n' "${md5_out}" | grep -v ': OK$' | grep -v 'computed checksum' | sed 's/^/md5     /' | head -20

# ---- summary -------------------------------------------------------------------------
failed_lines=$(( n_counts_bad + md5_failed ))
counts_ok=$(( n_counts > 0 ? n_counts - n_counts_bad : 0 ))
if (( failed_lines == 0 )); then
    echo "PASS: 0 failed lines (counts ${counts_ok} of ${n_counts} identical; md5 ${md5_ok} of ${md5_total} lines OK)"
    exit 0
fi
echo "FAIL: ${failed_lines} failed line(s) (counts ${counts_ok} of ${n_counts} identical; md5 ${md5_failed} of ${md5_total} lines failed)"
exit 1
