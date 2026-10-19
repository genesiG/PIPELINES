#!/bin/bash
# run_in_container.sh — runs the pipeline INSIDE the pipeline's Docker image.
# Started by run_crisprscreen (never by hand); the Docker counterpart of the
# cluster version's LSF head job. Mounts (set up by run_crisprscreen):
#   /project   the user's project folder (config.yaml, CSVs, original_data/, results/ ...)
#   /pipeline  this pipeline folder, read-only
#   /inputs/fastq_r1, /inputs/fastq_r2, /inputs/fastq_dir
#              the raw FASTQ file(s) or folder named in config.yaml, read-only; only
#              mounted when they are outside the project folder (they are not moved).
#              Their paths in the container come in CRISPR_FASTQ_R1, CRISPR_FASTQ_R2
#              and CRISPR_FASTQ_DIR; the Windows paths in CRISPR_FASTQ_*_HOST (banner only).
#
#   run_in_container.sh CONFIG_FILE TIMESTAMP
#
# Writes (all inside the project folder, so they outlive the container):
#   metadata/nf_logs/nextflow_<TIMESTAMP>.log   everything Nextflow printed
#   metadata/nf_logs/last_run.status           running | succeeded | failed (exit N)
#   metadata/failed_tasks/<TIMESTAMP>/<step>/  for a failed run: each failed
#       step's .command.sh (what ran), .command.err/.command.out/.command.log
#       (what it said) and .exitcode, copied out of Nextflow's work folder.
set -uo pipefail

config="${1:-config.yaml}"
ts="${2:-$(date +%Y%m%d_%H%M%S)}"
cd /project

log_dir="metadata/nf_logs"
mkdir -p "${log_dir}"
log="${log_dir}/nextflow_${ts}.log"
status_file="${log_dir}/last_run.status"

# ---- Resources ---------------------------------------------------------------
# The run's budget was chosen by run_crisprscreen and is enforced as this
# container's hard limits (docker run --memory / --cpuset-cpus); it is passed
# in as CRISPR_BUDGET_CPUS / CRISPR_BUDGET_MEMORY_MB. Without them (container
# started some other way), the container's own limits are read instead.
cpus=$(nproc)                                       # honours --cpuset-cpus
mem_mb=$(( $(awk '/^MemTotal:/ { print $2 }' /proc/meminfo) / 1024 ))
if [[ -r /sys/fs/cgroup/memory.max ]]; then         # docker run --memory (cgroup v2)
    lim=$(cat /sys/fs/cgroup/memory.max)
    if [[ "${lim}" =~ ^[0-9]+$ ]] && (( lim / 1048576 < mem_mb )); then mem_mb=$(( lim / 1048576 )); fi
fi
if [[ "${CRISPR_BUDGET_CPUS:-}" =~ ^[0-9]+$ ]] && (( CRISPR_BUDGET_CPUS >= 1 && CRISPR_BUDGET_CPUS < cpus )); then
    cpus=${CRISPR_BUDGET_CPUS}
fi
if [[ "${CRISPR_BUDGET_MEMORY_MB:-}" =~ ^[0-9]+$ ]] && (( CRISPR_BUDGET_MEMORY_MB < mem_mb )); then
    mem_mb=${CRISPR_BUDGET_MEMORY_MB}
fi
# Nextflow itself (Java heap capped at 768 MB below, plus the JVM's own
# overhead) lives inside the same budget: the steps get the rest.
nf_reserve_mb=1024
steps_mb=$(( mem_mb - nf_reserve_mb ))
(( steps_mb < 1024 )) && steps_mb=1024
export CRISPR_MAX_CPUS="${cpus}"
export CRISPR_MAX_MEMORY_MB="${steps_mb}"

# ---- Nextflow settings ------------------------------------------------------
# work/ (inside the project folder) holds every step's working files; it is
# what lets a re-run resume instead of starting over. After a successful run
# the folder itself and Nextflow's .nextflow/ cache are deleted unless
# keep_work: true (cleanup_work below); a failed or stopped run keeps both.
# `run_crisprscreen clean` deletes them on demand.
export NXF_WORK="${CRISPR_WORK_DIR:-/project/work}"
export NXF_OPTS="-Xms128m -Xmx768m"
export NXF_ANSI_LOG=false
export NXF_OFFLINE=true
mkdir -p "${NXF_WORK}"

# ---- Raw FASTQ files outside the project folder -------------------------------
# run_crisprscreen mounts them read-only and names them in CRISPR_FASTQ_R1 / _R2 /
# _DIR. They go on the command line, so they override config.yaml's (Windows)
# paths; FASTQs inside the project folder need nothing here (config.yaml is used).
fq_args=()
if [[ -n "${CRISPR_FASTQ_R1:-}" ]]; then fq_args+=(--fastq_r1 "${CRISPR_FASTQ_R1}"); fi
if [[ -n "${CRISPR_FASTQ_R2:-}" ]]; then fq_args+=(--fastq_r2 "${CRISPR_FASTQ_R2}"); fi
if [[ -n "${CRISPR_FASTQ_DIR:-}" ]]; then fq_args+=(--fastq_dir "${CRISPR_FASTQ_DIR}"); fi

start_marker="$(mktemp)"
echo "running" > "${status_file}"

{
    echo "CRISPR_SCREEN (Docker) — run ${ts}"
    echo "  project   : ${CRISPR_PROJECT_HOST_DIR:-/project}"
    echo "  config    : ${config}"
    if [[ -n "${CRISPR_FASTQ_R1_HOST:-}" ]]; then echo "  fastq_r1  : ${CRISPR_FASTQ_R1_HOST}"; fi
    if [[ -n "${CRISPR_FASTQ_R2_HOST:-}" ]]; then echo "  fastq_r2  : ${CRISPR_FASTQ_R2_HOST}"; fi
    if [[ -n "${CRISPR_FASTQ_DIR_HOST:-}" ]]; then echo "  fastq_dir : ${CRISPR_FASTQ_DIR_HOST}"; fi
    echo "  image     : ${CRISPR_IMAGE:-?}"
    echo "  resources : ${cpus} CPUs, $(( mem_mb / 1024 )).$(( mem_mb % 1024 * 10 / 1024 )) GB memory (${CRISPR_BUDGET_REASON:-the container limits}); ${steps_mb} MB of it for the steps"
    echo "  log       : ${log}"
    echo ""
} | tee "${log}"

# In the background + wait, so `run_crisprscreen stop` (docker stop -> SIGTERM)
# reaches Nextflow, which then stops its steps cleanly and stays resumable.
nextflow -log metadata/.nextflow.log run /pipeline -params-file "${config}" -resume ${fq_args[@]+"${fq_args[@]}"} \
    > >(tee -a "${log}") 2>&1 &
nf_pid=$!
trap 'kill -TERM "${nf_pid}" 2>/dev/null' TERM INT
rc=0
wait "${nf_pid}" || rc=$?
while kill -0 "${nf_pid}" 2>/dev/null; do   # signalled: wait for Nextflow to finish shutting down
    rc=0; wait "${nf_pid}" || rc=$?
done
sleep 1   # let tee flush the last lines

# Deletes the work folder and Nextflow's .nextflow/ cache (the project's resume
# state, useless without work/) after a successful run (keep_work: false, the
# default). Only on success: a failed or stopped run keeps both so a re-run
# resumes. Safe because every publishDir copies (publish_mode 'copy'), so
# results/ never points into work/. keep_work: true disables it.
cleanup_work() {
    if grep -qiE '^keep_work:[[:space:]]*[^[:alnum:]]?(true|yes)([^[:alnum:]]|$)' "${config}" 2>/dev/null; then
        return 0
    fi
    # results/ only stays valid without work/ when steps publish copies (the default).
    if grep -qiE '^publish_mode:[[:space:]]*[^[:alnum:]]?(symlink|rellink)' "${config}" 2>/dev/null; then
        echo "Kept the work folder: publish_mode links results/ into it." | tee -a "${log}"
        return 0
    fi
    local kb gb what=""
    if [[ -n "${NXF_WORK}" && -d "${NXF_WORK}" ]]; then
        case "${NXF_WORK}" in
            /|/project) ;;
            *)
                kb=$(du -sk "${NXF_WORK}" 2>/dev/null | cut -f1)
                gb=$(awk -v k="${kb:-0}" 'BEGIN { printf "%.1f", k / 1048576 }')
                find "${NXF_WORK}" -mindepth 1 -delete 2>/dev/null || true
                case "${NXF_WORK}" in
                    /project/work)            # default: remove the (now empty) folder too
                        rmdir "${NXF_WORK}" 2>/dev/null || true
                        what="the work folder (${gb} GB)" ;;
                    *)                        # custom CRISPR_WORK_DIR: only emptied
                        what="the contents of ${NXF_WORK} (${gb} GB)" ;;
                esac
                ;;
        esac
    fi
    # Nextflow's resume cache and its logs, in the project root (launch dir).
    if [[ -e /project/.nextflow ]]; then
        rm -rf /project/.nextflow
        if [[ -n "${what}" ]]; then what="${what} and Nextflow's .nextflow/ cache"; else what="Nextflow's .nextflow/ cache"; fi
    fi
    rm -f /project/.nextflow.log*
    if [[ -n "${what}" ]]; then
        echo "Deleted ${what} to free space; set keep_work: true in ${config} to keep them for resuming." | tee -a "${log}"
    fi
}

if [[ "${rc}" -eq 0 ]]; then
    echo "succeeded" > "${status_file}"
    cleanup_work
else
    echo "failed (exit ${rc})" > "${status_file}"
    # Copy each failed step's own files out of work/, so the reason is in the
    # project folder next to the log, under the step's name.
    out="metadata/failed_tasks/${ts}"
    n=0
    while IFS= read -r ec; do
        code=$(cat "${ec}" 2>/dev/null || echo "?")
        [[ "${code}" == "0" ]] && continue
        dir=$(dirname "${ec}")
        name=$(sed -n 's/^# NEXTFLOW TASK: //p' "${dir}/.command.run" 2>/dev/null | head -1)
        name=$(printf '%s' "${name:-$(basename "${dir}")}" | tr -c 'A-Za-z0-9._()-' '_')
        mkdir -p "${out}/${name}"
        for f in .command.sh .command.err .command.out .command.log .exitcode; do
            [[ -f "${dir}/${f}" ]] && cp "${dir}/${f}" "${out}/${name}/${f#.}"
        done
        echo "work folder: work/${dir#${NXF_WORK}/}" > "${out}/${name}/workdir.txt"
        n=$((n + 1))
    done < <(find "${NXF_WORK}" -mindepth 3 -maxdepth 3 -name .exitcode -newer "${start_marker}" 2>/dev/null)
    {
        echo ""
        if (( n > 0 )); then
            echo "The files of ${n} failed step(s) were copied to ${out}/"
            echo "(command.sh = what ran; command.err / command.log = what it said)."
        fi
    } | tee -a "${log}"
fi
rm -f "${start_marker}"
exit "${rc}"
