#!/bin/bash
# run_script_in_container.sh — runs ONE project script inside the pipeline's Docker image.
# Started by `run_crisprscreen script` (never by hand). Mounts (set up by run_crisprscreen):
#   /project   the user's project folder (the script is inside it; the working directory is /project)
#   /pipeline  this pipeline folder, read-only
#
#   run_script_in_container.sh SCRIPT_REL LOG_REL [ARGS...]
#     SCRIPT_REL  the script's path relative to the project folder, forward slashes (e.g. utils/deseq.R)
#     LOG_REL     the log file, relative to the project folder (metadata/script_logs/<name>_<timestamp>.log)
#
# The log gets a header (what ran, with which script md5, arguments, image and budget), then the
# script's output (stdout and stderr), then a footer with the exit code. All of it is also echoed.
set -u

script_rel="${1:-}"
log_rel="${2:-}"
if [[ -z "${script_rel}" || -z "${log_rel}" ]]; then
    echo "usage: run_script_in_container.sh SCRIPT_REL LOG_REL [ARGS...]" >&2
    exit 2
fi
shift 2
cd /project || exit 2

# ---- The script must be a regular file inside the project folder ---------------
case "${script_rel}" in
    /*) echo "refused: ${script_rel} is an absolute path (give it relative to the project folder)" >&2; exit 2 ;;
esac
case "/${script_rel}/" in
    */../*) echo "refused: ${script_rel} contains '..'" >&2; exit 2 ;;
esac
if [[ ! -f "${script_rel}" ]]; then
    echo "not found: ${script_rel} (no file at that path in the project folder)" >&2
    exit 2
fi
real_path=$(realpath -- "${script_rel}" 2>/dev/null || true)
case "${real_path}" in
    /project/*) ;;
    *) echo "refused: ${script_rel} resolves to a place outside the project folder" >&2; exit 2 ;;
esac

# ---- Interpreter by extension ---------------------------------------------------
case "${script_rel,,}" in
    *.r)  interp=(Rscript) ;;
    *.py) interp=(python3) ;;
    *.sh) interp=(bash) ;;
    *)    echo "unsupported extension (use .R, .py or .sh): ${script_rel}" >&2; exit 2 ;;
esac

# ---- Resources (what the script may use; read by scripts as CRISPR_MAX_*) ------
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
export CRISPR_MAX_CPUS="${cpus}"
export CRISPR_MAX_MEMORY_MB="${mem_mb}"

# ---- Header, then the run -------------------------------------------------------
mkdir -p "$(dirname "${log_rel}")"
args_str="(none)"
if (( $# > 0 )); then args_str=$(printf '%q ' "$@"); args_str="${args_str% }"; fi
{
    echo "CRISPR_SCREEN script run"
    echo "  date      : $(date '+%Y-%m-%d %H:%M:%S %z')"
    echo "  project   : ${CRISPR_PROJECT_HOST_DIR:-/project}"
    echo "  script    : ${script_rel}"
    echo "  md5       : $(md5sum "${script_rel}" | cut -d' ' -f1)"
    echo "  args      : ${args_str}"
    echo "  image     : ${CRISPR_IMAGE:-?}"
    echo "  image id  : ${CRISPR_IMAGE_ID:-?}"
    echo "  budget    : ${cpus} CPUs, $(( mem_mb / 1024 )).$(( mem_mb % 1024 * 10 / 1024 )) GB memory (${CRISPR_BUDGET_REASON:-?})"
    echo "  working directory: /project"
    echo ""
    echo "------------------------------------------------------------------------------"
} | tee -a "${log_rel}"

"${interp[@]}" "${script_rel}" "$@" 2>&1 | tee -a "${log_rel}"
rc=${PIPESTATUS[0]}

echo "-- exit code ${rc} --" | tee -a "${log_rel}"
exit "${rc}"
