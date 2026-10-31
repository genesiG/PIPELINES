#!/bin/bash
# mageck_count.sh — `mageck count` of one library's samples.
# Called by MAGECK_COUNT (modules/local/mageck_count.nf):
#   mageck_count.sh --library lib_mageck_library.txt --fastq "s1.fastq.gz s2.fastq.gz" \
#                   --labels "s1,s2" --prefix OUT --sgrna-len 20 --trim5 22 \
#                   [--drop-zero-samples] [mageck count args...]
# --fastq and --labels are in the same (samplesheet) order.
# Any other arguments (config.yaml `mageck_count_args`) are passed to mageck count.
# Must write: OUT.count.txt, OUT.count_normalized.txt, OUT_mageck_count.log, OUT.unmapped.txt,
#             OUT.sample_reads.tsv (sample, fastq, reads, mapped, status — one row per
#             sample passed in; status = counted | dropped_zero_mapped)
# A sample with 0 mapped reads crashes mageck count (see below). By default
# that exits 3 (a data problem; the pipeline doesn't retry it), naming the
# samples. With --drop-zero-samples (main.nf passes it only when every sample
# is counted against every library, where a sample holding none of one
# library's guides is expected), those samples are left out of this
# library's count table instead: mageck count re-runs without them, the first
# attempt's log is kept as OUT_mageck_count_attempt1.log, and they are listed
# in OUT.sample_reads.tsv as dropped_zero_mapped. If no sample maps at all,
# it still exits 3.
#
# --unmapped-to-file: mageck count's own guide matching is always exact
# (verified via `mageck count --help` — there's no mismatch-tolerance flag;
# --sgrna-len only tells mageck how to format this file, since the library
# file itself determines the actual matching length). OUT.unmapped.txt is
# one row per DISTINCT unmapped sequence across every sample passed here (in
# --labels order), tab-separated: sequence, sequence again, then one count
# column per sample. guide_representation_wrangle.R reads it to classify
# unmapped reads as "near-miss" (1 substitution from a real guide — lost to
# sequencing/PCR error) vs "far-miss" (unrelated — contamination, wrong
# library, adapter reads).
set -euo pipefail

library="" fastqs="" labels="" prefix="" sgrna_len="" trim5="" drop_zero=0 extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --library)   library="$2"; shift 2 ;;
        --fastq)     fastqs="$2"; shift 2 ;;
        --labels)    labels="$2"; shift 2 ;;
        --prefix)    prefix="$2"; shift 2 ;;
        --sgrna-len) sgrna_len="$2"; shift 2 ;;
        --trim5)     trim5="$2"; shift 2 ;;
        --drop-zero-samples) drop_zero=1; shift ;;
        *)           extra+=("$1"); shift ;;
    esac
done
read -ra fastq_files <<< "${fastqs}"
IFS=',' read -ra sample_labels <<< "${labels}"
orig_labels="${labels}"

run_mageck() {   # $1 = log file; uses the current fastq_files / labels
    local rc=0
    mageck count \
        --list-seq      "${library}" \
        --fastq         "${fastq_files[@]}" \
        --sample-label  "${labels}" \
        --output-prefix "${prefix}" \
        --sgrna-len     "${sgrna_len}" \
        --trim-5        "${trim5}" \
        --norm-method   median \
        --unmapped-to-file \
        ${extra[@]+"${extra[@]}"} \
        2>&1 | tee "$1" || rc=$?
    return "${rc}"
}

# Per-FASTQ "Total:" / "Mapped:" lines of a mageck count log, as
# "<fastq>\t<reads>\t<mapped>" (mageck logs them per file before it
# normalizes, so they're present even when it then crashes).
log_reads() {
    awk '/Parsing FASTQ file/ { f = $NF; gsub(/[.][.][.]/, "", f) }
         /INFO .*: Total: [0-9]/ { t = $NF; gsub(/[.]/, "", t) }
         /INFO .*: Mapped: [0-9]/ { m = $NF; gsub(/[.]/, "", m); print f "\t" t "\t" m }' "$1"
}

# mageck count normalizes by each sample's total mapped reads (even with
# --norm-method median), so a single sample with 0 reads mapping to the
# library crashes it with ZeroDivisionError. Name those samples plainly.
log="${prefix}_mageck_count.log"
rc=0
run_mageck "${log}" || rc=$?
dropped=()
first_log=""
if [[ "${rc}" -ne 0 ]]; then
    zero=$(log_reads "${log}" | awk -F'\t' '$3 == "0" { print $1 }')
    if [[ -z "${zero}" ]]; then
        exit "${rc}"
    fi
    if [[ "${drop_zero}" -eq 0 ]]; then
        echo "ERROR: mageck count failed: these samples had 0 reads mapping to the library:" >&2
        echo "${zero}" | sed 's/^/    /' >&2
        echo "Remove them from samplesheet.csv (or check they belong to this library / trim5 / sgrna_len) and re-run." >&2
        exit 3
    fi
    # Every-sample-vs-every-library mode: leave the zero-mapped samples out of
    # this library's table and count the rest again.
    keep_fq=() keep_lab=()
    for i in "${!fastq_files[@]}"; do
        if grep -qxF "${fastq_files[$i]}" <<< "${zero}"; then
            dropped+=("${sample_labels[$i]}")
        else
            keep_fq+=("${fastq_files[$i]}"); keep_lab+=("${sample_labels[$i]}")
        fi
    done
    if [[ "${#dropped[@]}" -eq 0 ]]; then
        exit "${rc}"   # the failure wasn't one of this run's samples mapping 0 reads
    fi
    if [[ "${#keep_fq[@]}" -eq 0 ]]; then
        echo "ERROR: no sample has any read matching this library (${library}) — every sample mapped 0 reads." >&2
        echo "Check that this library file belongs to this screen (and trim5 / sgrna_len, if set)." >&2
        exit 3
    fi
    echo "WARNING: these samples have 0 reads matching ${library} and are left out of its count table:" >&2
    printf '    %s\n' "${dropped[@]}" >&2
    first_log="${prefix}_mageck_count_attempt1.log"
    mv "${log}" "${first_log}"
    fastq_files=("${keep_fq[@]}")
    sample_labels=("${keep_lab[@]}")
    labels=$(IFS=','; echo "${keep_lab[*]}")
    run_mageck "${log}"
fi

# OUT.sample_reads.tsv: every sample passed in, with its input and mapped reads.
{
    printf 'sample\tfastq\treads\tmapped\tstatus\n'
    log_reads "${log}" | while IFS=$'\t' read -r fq t m; do
        for i in "${!fastq_files[@]}"; do
            if [[ "${fastq_files[$i]}" == "${fq}" ]]; then
                printf '%s\t%s\t%s\t%s\tcounted\n' "${sample_labels[$i]}" "${fq}" "${t}" "${m}"
            fi
        done
    done
    if [[ -n "${first_log}" ]]; then
        IFS=' ' read -ra all_fq <<< "${fastqs}"
        IFS=',' read -ra all_lab <<< "${orig_labels}"
        log_reads "${first_log}" | while IFS=$'\t' read -r fq t m; do
            for i in "${!all_fq[@]}"; do
                if [[ "${all_fq[$i]}" == "${fq}" ]] && printf '%s\n' "${dropped[@]}" | grep -qxF "${all_lab[$i]}"; then
                    printf '%s\t%s\t%s\t%s\tdropped_zero_mapped\n' "${all_lab[$i]}" "${fq}" "${t}" "${m}"
                fi
            done
        done
    fi
} > "${prefix}.sample_reads.tsv"
