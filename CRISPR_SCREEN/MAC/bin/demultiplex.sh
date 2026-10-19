#!/bin/bash
# demultiplex.sh — split a multiplexed FASTQ (or R1/R2 pair) into one FASTQ per
# sample by inline 5' barcode (cutadapt demux mode), then tabulate where in
# the read each barcode matched (barcode_position_qc.py).
#
# Same cutadapt flags as the pre-Nextflow Scripts/step_1_demultiplex.py
# (-e 0, --no-indels, --discard-untrimmed, --overlap = barcode length), so
# counts stay comparable. --cores 0 uses every core this step may use: the LSF
# allocation on the cluster (conf/base.config gives the task no `cpus`
# directive), the run's CPU budget in the Docker versions (conf/base.config
# reserves all of it for this step).
#
# Called by DEMULTIPLEX_SE / DEMULTIPLEX_PE (modules/local/demultiplex.nf):
#   demultiplex.sh --r1 R1.fastq.gz [--r2 R2.fastq.gz] --adapters "-g S1=ACGT -g S2=..." \
#                  --overlap N --position-script barcode_position_qc.py --project NAME [cutadapt args...]
# Any other arguments (config.yaml `demultiplex_args`) are passed to cutadapt.
# Must write: <sample>.fastq.gz (paired: <sample>_R1/_R2.fastq.gz), cutadapt.log,
#             <project>_barcode_position_counts.csv
set -euo pipefail

r1="" r2="" adapters="" overlap="" position_script="" project="" extra=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --r1)              r1="$2"; shift 2 ;;
        --r2)              r2="$2"; shift 2 ;;
        --adapters)        adapters="$2"; shift 2 ;;
        --overlap)         overlap="$2"; shift 2 ;;
        --position-script) position_script="$2"; shift 2 ;;
        --project)         project="$2"; shift 2 ;;
        *)                 extra+=("$1"); shift ;;
    esac
done
read -ra adapter_args <<< "${adapters}"

if [[ -n "${r2}" ]]; then
    outputs=(--output '{name}_R1.fastq.gz' --paired-output '{name}_R2.fastq.gz' "${r1}" "${r2}")
else
    outputs=(--output '{name}.fastq.gz' "${r1}")
fi

# The info file (one line per read: ~100+ GB for a full lane — 180 GB were
# written per full-lane run on the cluster) is never written to disk:
# cutadapt writes it into a named pipe that barcode_position_qc reads as it
# goes, keeping only per-position counts. Same numbers as writing it out and
# reading it back, without the temporary disk space and its write time
# (cluster storage quota; SSD wear on a laptop).
# For paired-end reads the info file describes R1, where the -g barcodes are
# matched, which is what the position plot shows.
case "${position_script,,}" in
    *.r)  position_runner=Rscript ;;
    *.sh) position_runner=bash ;;
    *)    position_runner=python3 ;;
esac
# The pipe lives in local temporary space, not the work folder: that may be a
# Windows or Mac folder shared into Docker, or a cluster file system, where
# named pipes are not always supported. A directory of its own per task: on a
# cluster node, several demultiplexing jobs may share /tmp.
fifo_dir="$(mktemp -d "${TMPDIR:-/tmp}/crisprscreen_demultiplex.XXXXXX")"
info="${fifo_dir}/info.tsv"
trap 'rm -rf "${fifo_dir}"' EXIT
mkfifo "${info}"
"${position_runner}" "${position_script}" \
    --info-file "${info}" \
    --project-name "${project}" \
    --outdir . &
qc_pid=$!
# If the reader dies before opening the pipe (e.g. a broken replacement
# script), cutadapt would wait for it forever: open and close the pipe's read
# end then, so cutadapt fails at once instead.
( while kill -0 "${qc_pid}" 2>/dev/null; do sleep 2; done
  exec 3<"${info}"; exec 3<&- ) &
watch_pid=$!

rc=0
cutadapt \
    "${adapter_args[@]}" \
    -e 0 \
    --overlap "${overlap}" \
    --discard-untrimmed \
    --no-indels \
    --cores 0 \
    --info-file "${info}" \
    ${extra[@]+"${extra[@]}"} \
    "${outputs[@]}" \
    2>&1 | tee cutadapt.log || rc=$?
kill "${watch_pid}" 2>/dev/null || true
if [[ "${rc}" -ne 0 ]]; then
    # If cutadapt failed before opening the pipe, the reader is still waiting for it.
    kill "${qc_pid}" 2>/dev/null || true
    wait "${qc_pid}" 2>/dev/null || true
    exit "${rc}"
fi
wait "${qc_pid}"
