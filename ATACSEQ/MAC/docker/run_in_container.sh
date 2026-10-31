#!/bin/bash
# run_in_container.sh — runs the pipeline INSIDE the pipeline's Docker image.
# Started by run_atacseq (never by hand). Mounts (set up by run_atacseq):
#   /project   the user's project folder (config.yaml, CSVs, original_data/, results/ ...)
#   /pipeline  this pipeline folder, read-only
#   /nfwork    a Docker volume (per project) that holds Nextflow's work folder and resume
#              cache when ATAC_WORK_VOLUME=1 (work_volume: true, the default). It lives on
#              Docker's own disk, which is much faster than the project folder's bind mount.
#   /mnt/c, /mnt/d ...   the Windows drives (or /Users on a Mac) that config.yaml and the
#              CSVs name, read-only; their list comes in ATAC_MOUNTS_HOST (banner only).
#              A path written as D:\FASTQ\x.fastq.gz is read as /mnt/d/FASTQ/x.fastq.gz
#              (main.nf hostPath()); nothing is moved or copied.
#
#   run_in_container.sh CONFIG_FILE TIMESTAMP
#
# Writes (all inside the project folder, so they outlive the container):
#   metadata/nf_logs/nextflow_<TIMESTAMP>.log   everything Nextflow printed
#   metadata/nf_logs/last_run.status           running | succeeded | failed (exit N)
#   metadata/failed_tasks/<TIMESTAMP>/<step>/  for a failed run: each failed
#       step's .command.sh (what ran), .command.err/.command.out/.command.log
#       (what it said) and .exitcode, copied out of Nextflow's work folder.
# Nextflow's own trace, report and timeline are named with a launch timestamp by
# nextflow.config (metadata/pipeline_info/), so a resumed run does not overwrite them.
# The work folder and Nextflow's resume cache are NOT in the project folder when
# ATAC_WORK_VOLUME=1: they are in the /nfwork volume (see Nextflow settings below).
# With ATAC_WORK_VOLUME=0 (work_volume: false) they are metadata/work and metadata/.nextflow.
set -uo pipefail

config="${1:-config.yaml}"
ts="${2:-$(date +%Y%m%d_%H%M%S)}"
cd /project

log_dir="metadata/nf_logs"
mkdir -p "${log_dir}"
log="${log_dir}/nextflow_${ts}.log"
status_file="${log_dir}/last_run.status"

# ---- Resources ---------------------------------------------------------------
# The run's budget was chosen by run_atacseq and is enforced as this
# container's hard limits (docker run --memory / --cpus); it is passed
# in as ATAC_BUDGET_CPUS / ATAC_BUDGET_MEMORY_MB. Without them (container
# started some other way), the container's own limits are read instead.
cpus=$(nproc)                                       # honours --cpuset-cpus
mem_mb=$(( $(awk '/^MemTotal:/ { print $2 }' /proc/meminfo) / 1024 ))
if [[ -r /sys/fs/cgroup/memory.max ]]; then         # docker run --memory (cgroup v2)
    lim=$(cat /sys/fs/cgroup/memory.max)
    if [[ "${lim}" =~ ^[0-9]+$ ]] && (( lim / 1048576 < mem_mb )); then mem_mb=$(( lim / 1048576 )); fi
fi
if [[ "${ATAC_BUDGET_CPUS:-}" =~ ^[0-9]+$ ]] && (( ATAC_BUDGET_CPUS >= 1 && ATAC_BUDGET_CPUS < cpus )); then
    cpus=${ATAC_BUDGET_CPUS}
fi
if [[ "${ATAC_BUDGET_MEMORY_MB:-}" =~ ^[0-9]+$ ]] && (( ATAC_BUDGET_MEMORY_MB < mem_mb )); then
    mem_mb=${ATAC_BUDGET_MEMORY_MB}
fi
# Nextflow itself (Java heap capped at 768 MB below, plus the JVM's own
# overhead) lives inside the same budget: the steps get the rest.
nf_reserve_mb=1024
steps_mb=$(( mem_mb - nf_reserve_mb ))
(( steps_mb < 1024 )) && steps_mb=1024
export ATAC_MAX_CPUS="${cpus}"
export ATAC_MAX_MEMORY_MB="${steps_mb}"

# ---- Nextflow settings ------------------------------------------------------
# The work folder holds every step's working files; it is what lets a re-run
# resume instead of starting over. Nextflow's resume cache (NXF_CACHE_DIR) holds
# the run history and the record of every task (absolute work paths, file times).
# The two must stay together, so they live in the same place:
#   ATAC_WORK_VOLUME=1 (work_volume: true, the default): /nfwork/work and
#     /nfwork/.nextflow, both in the per-project Docker volume. The volume is on
#     Docker's own disk, not in the project folder (a bind mount there is about ten
#     times slower). Created, kept and deleted together.
#   ATAC_WORK_VOLUME=0 (work_volume: false, or unset): metadata/work and
#     metadata/.nextflow inside the project folder.
# Both are kept after a run by default (keep_work: true); keep_work: false deletes
# them after a successful run (cleanup_work below); a failed or stopped run always
# keeps both. The launcher's clean command deletes them on demand.
# The Nextflow log (-log below) stays in the project folder: metadata/.nextflow.log.
# NXF_WORK is read by Nextflow itself and is also passed as -work-dir below.
# NXF_HOME (Nextflow's own folder, default $HOME/.nextflow) is left at its default:
# it is inside the container, not in /pipeline (read-only) or /project.
work_vol=0
if [[ "${ATAC_WORK_VOLUME:-}" == "1" ]]; then
    work_vol=1
    export NXF_WORK="${ATAC_WORK_DIR:-/nfwork/work}"
    export NXF_CACHE_DIR=/nfwork/.nextflow
else
    export NXF_WORK="${ATAC_WORK_DIR:-/project/metadata/work}"
    export NXF_CACHE_DIR=/project/metadata/.nextflow
fi
export NXF_OPTS="-Xms128m -Xmx768m"
export NXF_ANSI_LOG=false
export NXF_OFFLINE=true
# Nextflow version: the pin of the cluster launcher (26.04.6).
export NXF_VER=26.04.6
mkdir -p "${NXF_WORK}"

# Nextflow profile(s) and extra Nextflow options for test runs only (the launcher
# does not set them): ATAC_NF_PROFILE (default docker, e.g. docker,test) and
# ATAC_NF_EXTRA_ARGS (word-split, e.g. "-preview" or "-stub-run").
nf_profile="${ATAC_NF_PROFILE:-docker}"
nf_extra=()
read -r -a nf_extra <<< "${ATAC_NF_EXTRA_ARGS:-}"

# ---- Publish mode -------------------------------------------------------------
# publish_mode 'auto' (the default in nextflow.config): results/ and original_data/
# get hard links to the work folder's files when the disk allows (a big file is
# then stored once), else copies. An explicit value in config.yaml is used as it
# is (nothing is passed on the command line). Logged in the run's banner below.
pub_cfg=$(grep -iE '^publish_mode:' "${config}" 2>/dev/null | head -1 | tr -d '\r' \
    | sed -E 's/^[^:]*:[[:space:]]*//; s/(^|[[:space:]])#.*$//; s/^["'"'"']//; s/["'"'"'][[:space:]]*$//; s/[[:space:]]+$//')
# With the work volume (ATAC_WORK_VOLUME=1) the hard-link test is skipped: hard links
# and symlinks cannot cross from the volume to the project folder. Explicit links
# are refused (exit 2, before Nextflow starts); auto or empty means copy.
pub_args=()
pub_lc=$(printf '%s' "${pub_cfg}" | tr '[:upper:]' '[:lower:]')
if (( work_vol )) && [[ "${pub_lc}" == "symlink" || "${pub_lc}" == "rellink" || "${pub_lc}" == "link" ]]; then
    {
        echo "ERROR: publish_mode: ${pub_cfg} cannot be used with the work volume: hard links and links cannot cross filesystems, and results/ would point into the volume."
        echo "Set publish_mode: auto or copy in ${config}, or work_volume: false (the work folder then stays in the project folder)."
    } | tee "${log}"
    echo "failed (exit 2)" > "${status_file}"
    exit 2
elif (( work_vol )) && [[ -z "${pub_lc}" || "${pub_lc}" == "auto" ]]; then
    pub_args=(--publish_mode copy)
    pub_note="publish mode: copy (the work folder is a Docker volume; hard links cannot cross filesystems)"
elif [[ -n "${pub_cfg}" && "${pub_cfg}" != "auto" ]]; then
    pub_note="publish mode: ${pub_cfg} (from ${config})"
else
    # Hard-link test: a temp file in the work folder, linked to a temp name in
    # metadata/ (where results/ and original_data/ are written). Both temp files
    # are removed again, whatever the result.
    lt_src=$(mktemp "${NXF_WORK}/.linktest.XXXXXX" 2>/dev/null || true)
    lt_dst="/project/metadata/.linktest.$$"
    hl=0
    if [[ -n "${lt_src}" ]] && ln "${lt_src}" "${lt_dst}" 2>/dev/null \
        && [[ "$(stat -c %h "${lt_src}" 2>/dev/null)" == "2" ]]; then
        hl=1
    fi
    rm -f "${lt_dst}"
    if [[ -n "${lt_src}" ]]; then rm -f "${lt_src}"; fi
    if (( hl )); then
        pub_args=(--publish_mode link)
        pub_note="publish mode: link (hard links from the work folder: big files are not stored twice)"
    else
        pub_args=(--publish_mode copy)
        pub_note="publish mode: copy (hard links are not possible between the work folder and the project folder)"
    fi
fi

start_marker="$(mktemp)"
echo "running" > "${status_file}"

{
    echo "ATACSEQ (Docker) - run ${ts}"
    echo "  project   : ${ATAC_PROJECT_HOST_DIR:-/project}"
    if [[ -n "${ATAC_MOUNTS_HOST:-}" ]]; then
        mounts_show=$(printf '%s' "${ATAC_MOUNTS_HOST}" | tr '\n' ',' | sed 's/,$//; s/,/, /g')
        echo "  inputs    : read-only from ${mounts_show}"
    fi
    echo "  config    : ${config}"
    echo "  image     : ${ATAC_IMAGE:-?}"
    if (( work_vol )); then
        echo "  work folder : Docker volume ${ATAC_WORK_VOLUME_NAME:-?} (/nfwork/work)"
    else
        echo "  work folder : metadata/work (project folder)"
    fi
    echo "  resources : ${cpus} CPUs, $(( mem_mb / 1024 )).$(( mem_mb % 1024 * 10 / 1024 )) GB memory (${ATAC_BUDGET_REASON:-the container limits}); ${steps_mb} MB of it for the steps"
    echo "  log       : ${log}"
    echo "  ${pub_note}"
    echo ""
} | tee "${log}"

# In the background + wait, so `run_atacseq stop` (docker stop -> SIGTERM)
# reaches Nextflow, which then stops its steps cleanly and stays resumable.
nextflow -log metadata/.nextflow.log run /pipeline -profile "${nf_profile}" -params-file "${config}" -resume -work-dir "${NXF_WORK}" \
    ${pub_args[@]+"${pub_args[@]}"} ${nf_extra[@]+"${nf_extra[@]}"} \
    > >(tee -a "${log}") 2>&1 &
nf_pid=$!
trap 'kill -TERM "${nf_pid}" 2>/dev/null' TERM INT
rc=0
wait "${nf_pid}" || rc=$?
while kill -0 "${nf_pid}" 2>/dev/null; do   # signalled: wait for Nextflow to finish shutting down
    rc=0; wait "${nf_pid}" || rc=$?
done
sleep 1   # let tee flush the last lines

# Deletes the work folder and Nextflow's resume cache (metadata/.nextflow) after a
# successful run, but only when config.yaml says keep_work: false. Without that
# line (keep_work: true is the default) both are kept, so a re-run resumes. A
# failed or stopped run always keeps both. Safe with publish_mode copy and link:
# results/ holds its own copies or hard links, never a pointer into work/. Not
# with symlink/rellink, which keep the work folder (guard below).
cleanup_work() {
    if ! grep -qiE '^keep_work:[[:space:]]*[^[:alnum:]]?(false|no)([^[:alnum:]]|$)' "${config}" 2>/dev/null; then
        return 0
    fi
    # results/ would point into the work folder with symlink/rellink: keep it then.
    if grep -qiE '^publish_mode:[[:space:]]*[^[:alnum:]]?(symlink|rellink)' "${config}" 2>/dev/null; then
        echo "Kept the work folder: publish_mode links results/ into it." | tee -a "${log}"
        return 0
    fi
    local kb gb what=""
    if (( work_vol )); then
        # Work volume: the whole volume goes (work/ and .nextflow/ together, so the
        # resume cache never outlives the work files it points to).
        case "${NXF_WORK}" in
            ''|/|/project) return 0 ;;
        esac
        if [[ -d /nfwork ]]; then
            kb=$(du -sk /nfwork 2>/dev/null | cut -f1)
            gb=$(awk -v k="${kb:-0}" 'BEGIN { printf "%.1f", k / 1048576 }')
            find /nfwork -mindepth 1 -delete 2>/dev/null || true
            echo "Deleted the work volume contents (${gb} GB) because keep_work: false in ${config}." | tee -a "${log}"
        fi
        return 0
    fi
    if [[ -n "${NXF_WORK}" && -d "${NXF_WORK}" ]]; then
        case "${NXF_WORK}" in
            /|/project) ;;
            *)
                kb=$(du -sk "${NXF_WORK}" 2>/dev/null | cut -f1)
                gb=$(awk -v k="${kb:-0}" 'BEGIN { printf "%.1f", k / 1048576 }')
                find "${NXF_WORK}" -mindepth 1 -delete 2>/dev/null || true
                case "${NXF_WORK}" in
                    /project/metadata/work)   # default: remove the (now empty) folder too
                        rmdir "${NXF_WORK}" 2>/dev/null || true
                        what="the work folder (${gb} GB)" ;;
                    *)                        # custom ATAC_WORK_DIR: only emptied
                        what="the contents of ${NXF_WORK} (${gb} GB)" ;;
                esac
                ;;
        esac
    fi
    # Nextflow's resume cache (NXF_CACHE_DIR).
    if [[ -e /project/metadata/.nextflow ]]; then
        rm -rf /project/metadata/.nextflow
        if [[ -n "${what}" ]]; then what="${what} and Nextflow's metadata/.nextflow cache"; else what="Nextflow's metadata/.nextflow cache"; fi
    fi
    if [[ -n "${what}" ]]; then
        echo "Deleted ${what} to free space because keep_work: false in ${config}." | tee -a "${log}"
    fi
}

if [[ "${rc}" -eq 0 ]]; then
    echo "succeeded" > "${status_file}"
    cleanup_work
else
    echo "failed (exit ${rc})" > "${status_file}"
    # Copy each failed step's own files out of the work folder (metadata/work/), so the reason is in the
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
        if (( work_vol )); then
            echo "work folder: Docker volume ${ATAC_WORK_VOLUME_NAME:-?}: ${dir#${NXF_WORK}/}" > "${out}/${name}/workdir.txt"
        else
            echo "work folder: ${NXF_WORK#/project/}/${dir#${NXF_WORK}/}" > "${out}/${name}/workdir.txt"
        fi
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
