#!/bin/bash
# run_crisprscreen.sh — CRISPR_SCREEN pipeline launcher, Docker version for Mac.
# (Also works on Linux with Docker Engine.) Written for macOS's built-in
# bash 3.2: no associative arrays, no ${x,,}, no readlink -f, BSD sed/stat.
#
#   run_crisprscreen.sh setup            one-time: check/install Docker Desktop, build the
#                                        pipeline's Docker image, put this command on PATH
#   run_crisprscreen.sh init [FOLDER]    set up a new project folder (default: the current folder)
#   run_crisprscreen.sh [run] [CONFIG]   run the pipeline for the project in the current folder
#                                        (CONFIG defaults to config.yaml)
#   run_crisprscreen.sh status           is it running? how far along? did it fail, and where?
#   run_crisprscreen.sh logs             follow the live log of the running analysis
#   run_crisprscreen.sh stop             stop a running analysis (re-running resumes it)
#   run_crisprscreen.sh clean            delete the work folder and resume cache (metadata/work, metadata/.nextflow) now
#   run_crisprscreen.sh script utils/X.R [ARGS]  run one of the project's own scripts (R, Python or shell) in the pipeline's software environment
#   run_crisprscreen.sh help
#
# Options: --detach (run: start and return, don't follow the log), --yes
# (answer yes to questions: for AI agents / scripts), --pause (wait for Enter
# before closing; used by the double-clickable RUN_ANALYSIS.command),
# --rebuild (setup: rebuild the Docker image even if it is up to date).
#
# What `run` does (the Docker counterpart of the cluster's run_crisprscreen.sh):
#   1. Checks the FASTQ paths named in config.yaml (fastq_r1, fastq_r2 or
#      fastq_dir). They stay where they are: never moved, copied or renamed, and
#      config.yaml is not rewritten. Outside the project folder they are
#      bind-mounted read-only under /inputs/...; on a Mac, the folder must be one
#      that Docker Desktop shares. Any other path in config.yaml must be inside
#      the project folder (it is rewritten as a relative path if needed).
#   2. Makes sure Docker Desktop is running and the image is built and up to
#      date with docker/ (rebuilt automatically after a pipeline update).
#   3. Starts one container that runs Nextflow and every step
#      (docker/run_in_container.sh), with the project folder mounted at
#      /project and this pipeline folder (read-only) at /pipeline, then
#      follows its log. Closing the window does NOT stop the analysis.

set -u -o pipefail

PIPELINE_NAME="crisprscreen"
IMAGE="${CRISPR_IMAGE_OVERRIDE:-rephilllab/crisprscreen:0.1.0}"

# This script's real folder, following symlinks (e.g. one in ~/bin).
src="${BASH_SOURCE[0]}"
while [ -h "${src}" ]; do
    dir="$(cd -P "$(dirname "${src}")" && pwd)"
    src="$(readlink "${src}")"
    case "${src}" in /*) ;; *) src="${dir}/${src}" ;; esac
done
PIPELINE_DIR="$(cd -P "$(dirname "${src}")" && pwd)"

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
PAUSE=0; YES=0; DETACH=0; REBUILD=0
positional=()
pass_rest=0   # after the command `script`, every argument belongs to the script
for a in "$@"; do
    if (( pass_rest )); then positional+=("${a}"); continue; fi
    case "${a}" in
        --pause|-pause)     PAUSE=1 ;;
        --yes|-yes|-y)      YES=1 ;;
        --detach|-detach)   DETACH=1 ;;
        --rebuild|-rebuild) REBUILD=1 ;;
        *)                  positional+=("${a}")
                            if [[ "${a}" == script ]] && (( ${#positional[@]} == 1 )); then pass_rest=1; fi ;;
    esac
done
COMMAND="run"
if (( ${#positional[@]} > 0 )); then
    case "${positional[0]}" in
        setup|init|run|status|logs|stop|clean|script|help|-h|--help)
            COMMAND="${positional[0]}"
            positional=("${positional[@]:1}") ;;
    esac
fi
npos=${#positional[@]}
arg1=""
(( npos > 0 )) && arg1="${positional[0]}"

# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------
if [[ -t 1 ]]; then C_INFO=$'\033[36m'; C_WARN=$'\033[33m'; C_ERR=$'\033[31m'; C_OFF=$'\033[0m'
else C_INFO=""; C_WARN=""; C_ERR=""; C_OFF=""; fi
say()  { printf '%s\n' "$*"; }
info() { printf '%s%s%s\n' "${C_INFO}" "$*" "${C_OFF}"; }
warn() { printf '%sWARNING: %s%s\n' "${C_WARN}" "$*" "${C_OFF}" >&2; }
done_() {
    if (( PAUSE )); then echo; read -r -p "Press Enter to close this window " _ || true; fi
    exit "$1"
}
fail() { echo; printf '%sERROR: %s%s\n' "${C_ERR}" "$*" "${C_OFF}" >&2; done_ 1; }
confirm() {
    (( YES )) && return 0
    [[ -t 0 ]] || return 1
    local ans; read -r -p "$1 [y/N] " ans || return 1
    [[ "${ans}" == y || "${ans}" == Y || "${ans}" == yes || "${ans}" == YES ]]
}
have() { command -v "$1" >/dev/null 2>&1; }
is_mac() { [[ "$(uname -s)" == Darwin ]]; }

# Absolute, symlink-resolved path of an existing file or folder.
abspath() {
    if [[ -d "$1" ]]; then (cd -P "$1" && pwd)
    else printf '%s/%s\n' "$(cd -P "$(dirname "$1")" && pwd)" "$(basename "$1")"; fi
}

# The file a symlink finally points to: one link at a time (no readlink -f on
# macOS). A relative link target is relative to the link's folder. At most 20
# links; a path that is not a link comes back as abspath gives it.
resolve_symlinks() {   # path
    local p="$1" i=0 dir target
    while [[ -h "${p}" ]] && (( i < 20 )); do
        target="$(readlink "${p}")"
        case "${target}" in
            /*) p="${target}" ;;
            *)  dir="$(cd -P "$(dirname "${p}")" && pwd)"; p="${dir}/${target}" ;;
        esac
        i=$(( i + 1 ))
    done
    abspath "${p}"
}

# ---------------------------------------------------------------------------
# Docker
# ---------------------------------------------------------------------------
install_docker() {
    echo
    say "Docker Desktop is not installed. The pipeline runs inside Docker, so it is needed (one time)."
    if ! is_mac; then
        fail "Install Docker Engine (https://docs.docker.com/engine/install/), then run this again."
    fi
    if have brew && confirm "Install Docker Desktop now with Homebrew (brew install --cask docker)?"; then
        info "Installing Docker Desktop..."
        brew install --cask docker || fail "Homebrew could not install Docker Desktop."
        say ""
        say "Docker Desktop is installed. Next:"
        say "  1. Open Docker from Applications and accept its terms (no account needed: 'Skip')."
        say "  2. Run this command again."
        done_ 0
    fi
    local chip="Apple chip"; [[ "$(uname -m)" == x86_64 ]] && chip="Intel chip"
    fail "Please install Docker Desktop for Mac (${chip}) from https://www.docker.com/products/docker-desktop/ ,
open it once from Applications and accept its terms (no account needed: 'Skip'), then run this again."
}

# Makes sure the docker command exists and the Docker engine is running.
ensure_docker() {   # $1 = allow_install (1/0)
    if ! have docker; then
        local d
        for d in /Applications/Docker.app/Contents/Resources/bin "${HOME}/.docker/bin" /usr/local/bin /opt/homebrew/bin; do
            [[ -x "${d}/docker" ]] && export PATH="${d}:${PATH}"
        done
    fi
    if ! have docker; then
        [[ "${1:-0}" == 1 ]] && install_docker
        fail "Docker is not installed. Run:  bash \"${PIPELINE_DIR}/run_crisprscreen.sh\" setup"
    fi
    docker info >/dev/null 2>&1 && return 0

    if is_mac && [[ -d /Applications/Docker.app ]]; then
        info "Starting Docker Desktop (this can take a minute)..."
        open -a Docker
        local i
        for i in $(seq 1 48); do
            sleep 5
            if docker info >/dev/null 2>&1; then say "Docker is running."; return 0; fi
        done
        fail "Docker Desktop did not start within 4 minutes. Open Docker from Applications, wait until it
says 'Engine running' (accept its terms the first time), then run this again."
    fi
    fail "Docker is installed but its engine is not running. Start Docker (Desktop) and try again."
}

# Apple Silicon: the image is linux/amd64, which Docker Desktop runs through
# Rosetta 2 (fast) — make sure Rosetta is there.
ensure_rosetta() {
    is_mac || return 0
    [[ "$(uname -m)" == arm64 ]] || return 0
    arch -x86_64 /usr/bin/true >/dev/null 2>&1 && return 0
    say "This Mac has an Apple chip; the pipeline's tools need Rosetta 2 (Apple's Intel-app translator)."
    if confirm "Install Rosetta 2 now (softwareupdate --install-rosetta)?"; then
        softwareupdate --install-rosetta --agree-to-license || warn "Rosetta 2 could not be installed; Docker will use slower emulation."
    else
        warn "without Rosetta 2, Docker uses much slower emulation."
    fi
}

# SHA-256 of what the image is built from (docker/Dockerfile + docker/env/*.explicit.txt),
# line endings normalized; stored as an image label, so a changed docker/
# folder (e.g. after `git pull`) triggers a rebuild.
docker_context_hash() {
    local f
    {
        for f in "${PIPELINE_DIR}/docker/Dockerfile" "${PIPELINE_DIR}"/docker/env/*.explicit.txt; do
            [[ -f "${f}" ]] || continue
            printf '%s\n' "$(basename "${f}")"
            tr -d '\r' < "${f}"
        done
    } | { shasum -a 256 2>/dev/null || sha256sum; } | cut -d' ' -f1
}

ensure_image() {   # $1 = force (1/0)
    local want have
    want="$(docker_context_hash)"
    have="$(docker image inspect "${IMAGE}" --format '{{ index .Config.Labels "crisprscreen.docker_hash" }}' 2>/dev/null)"
    local exists=$?
    if [[ "${1:-0}" != 1 && ${exists} -eq 0 && "${have}" == "${want}" ]]; then return 0; fi
    if [[ ${exists} -eq 0 && "${1:-0}" != 1 ]]; then
        info "The pipeline's Docker image is out of date (docker/ changed) — rebuilding it."
    else
        info "Building the pipeline's Docker image (${IMAGE})."
        say  "This happens once (and again only after a pipeline update changes docker/)."
        say  "It downloads about 3 GB and takes roughly 10-30 minutes. Keep the computer awake and online."
    fi
    ensure_rosetta
    docker build --platform linux/amd64 --progress=plain --label "crisprscreen.docker_hash=${want}" \
        -t "${IMAGE}" "${PIPELINE_DIR}/docker" \
        || fail "Building the Docker image failed (see the messages above). Check the internet connection and free disk space, then run this again."
    info "Image ${IMAGE} is ready."
}

# Bash scripts must have LF line endings inside the container; repair a copy
# that came through Windows.
repair_line_endings() {
    local f
    for f in "${PIPELINE_DIR}"/bin/*.sh "${PIPELINE_DIR}"/docker/*.sh "${PIPELINE_DIR}"/docker/env/*.sh; do
        [[ -f "${f}" ]] || continue
        if grep -q $'\r' "${f}"; then
            tr -d '\r' < "${f}" > "${f}.tmp" && cat "${f}.tmp" > "${f}" && rm -f "${f}.tmp"
            warn "converted $(basename "${f}") to Unix line endings (needed inside the container)"
        fi
    done
}

# ---------------------------------------------------------------------------
# config.yaml: flat keys only (as on the cluster), so simple line edits do.
# ---------------------------------------------------------------------------
yaml_get() {   # file key
    grep -E "^$2:" "$1" | head -1 \
        | sed -E "s/^$2:[[:space:]]*//; s/[[:space:]]+#.*\$//; s/[[:space:]]*\$//; s/^['\"]//; s/['\"]\$//" | tr -d '\r'
}
yaml_set() {   # file key value — rewrites the whole `key:` line (portable: no sed -i)
    local file="$1" key="$2" value="$3" tmp
    case "${value}" in
        *[:#\ ]*|[\[\{\&\*\!\|\>%@\`]*) value="'${value//\'/\'\'}'" ;;
    esac
    tmp="${file}.tmp.$$"
    awk -v k="${key}" -v v="${value}" 'BEGIN { done = 0 }
        !done && index($0, k ":") == 1 { print k ": " v; done = 1; next } { print }' "${file}" > "${tmp}" \
        && cat "${tmp}" > "${file}" && rm -f "${tmp}"
}
is_unset() { [[ -z "$1" || "$1" == null || "$1" == "~" ]]; }
is_abs()   { [[ "$1" == /* || "$1" == "~"* ]]; }

# A path from config.yaml as an absolute path (relative = relative to the project).
resolve_config_path() {
    local v="$1"
    case "${v}" in
        "~")   v="${HOME}" ;;
        "~/"*) v="${HOME}/${v#\~/}" ;;
        /*)    ;;
        *)     v="${PROJECT_DIR}/${v}" ;;
    esac
    if [[ -e "${v}" ]]; then abspath "${v}"; else printf '%s\n' "${v}"; fi
}

# Path relative to the project folder, or empty if outside it.
project_relative() {
    case "$1" in
        "${PROJECT_DIR}")   printf '.\n' ;;
        "${PROJECT_DIR}"/*) printf '%s\n' "${1#${PROJECT_DIR}/}" ;;
        *)                  printf '' ;;
    esac
}

# On a Mac, Docker Desktop shares only these folders with containers by default.
mac_docker_shared() {   # physical path
    case "$1" in
        /Users/*|/Volumes/*|/private/*|/tmp/*|/var/folders/*) return 0 ;;
    esac
    return 1
}

# FASTQ paths from config.yaml: fastq_r1 and fastq_r2 are files, fastq_dir is a
# folder. They stay where they are (never moved, copied or renamed; config.yaml
# is not rewritten). Sets RAW_PATH to the absolute path as written ("" if the key
# is unset). A file that is a symlink is followed first: the folder mounted is
# the one holding the real file. A path outside the project folder is
# bind-mounted read-only into the container: its docker run arguments are
# appended to input_args.
resolve_raw_fastq() {   # key
    local key="$1" v full real rel mount_dir upper target cpath found real_rel
    v="$(yaml_get "${CONFIG_FILE}" "${key}")"
    RAW_PATH=""
    is_unset "${v}" && return 0
    full="$(resolve_config_path "${v}")"
    if [[ "${key}" == fastq_dir ]]; then
        [[ -d "${full}" ]] || fail "${key} in ${CONFIG_NAME} is not a folder: ${v} (resolved to ${full})"
        shopt -s nullglob
        found=( "${full}"/*.fastq.gz "${full}"/*.fastq )
        shopt -u nullglob
        (( ${#found[@]} > 0 )) || fail "${key} in ${CONFIG_NAME}: ${v} (resolved to ${full}) contains no *.fastq.gz or *.fastq files."
        real="${full}"
        mount_dir="${full}"
    else
        [[ -f "${full}" ]] || fail "${key} in ${CONFIG_NAME} points to a file that doesn't exist: ${v} (resolved to ${full})"
        real="$(resolve_symlinks "${full}")"
        mount_dir="$(dirname "${real}")"
    fi
    RAW_PATH="${full}"
    rel="$(project_relative "${full}")"
    upper="$(printf '%s' "${key}" | tr '[:lower:]' '[:upper:]')"
    if [[ -n "${rel}" && -n "$(project_relative "${real}")" ]]; then
        say "  ${key}: ${rel} (in the project folder)"
        # Not mounted: the container sees the project folder at /project. The path
        # comes from the real file (a symlink may point elsewhere in the project).
        real_rel="$(project_relative "${real}")"
        cpath="/project"
        [[ "${real_rel}" == . ]] || cpath="/project/${real_rel}"
        input_args+=(-e "CRISPR_${upper}=${cpath}" -e "CRISPR_${upper}_HOST=${full}")
        return 0
    fi
    [[ "${real}" != *,* ]] \
        || fail "${key}: ${real} contains a comma, which breaks the container's mount. Rename or move it (no comma in the path), then run again."
    if is_mac && ! mac_docker_shared "${mount_dir}"; then
        fail "${key}: Docker Desktop cannot see ${mount_dir} (${real}).
Open Docker Desktop > Settings > Resources > File sharing, add the folder (or move the file under /Users or /Volumes), then run again."
    fi
    if [[ "${key}" == fastq_dir ]]; then
        target="/inputs/fastq_dir"; cpath="/inputs/fastq_dir"
    else
        target="/inputs/${key}"; cpath="/inputs/${key}/$(basename "${real}")"
    fi
    input_args+=(--mount "type=bind,source=${mount_dir},target=${target},readonly"
                 -e "CRISPR_${upper}=${cpath}" -e "CRISPR_${upper}_HOST=${full}")
    if [[ "${real}" == "${full}" ]]; then
        say "  ${key}: ${full} (stays where it is; mounted read-only)"
    else
        say "  ${key}: ${full} (a symlink to ${real}; the real file stays where it is; mounted read-only)"
    fi
    return 0
}

# Every other path in config.yaml must be inside the project folder (the
# container sees nothing else); written relative.
check_config_paths() {
    local keys="samplesheet library contrasts outdir original_data_dir intermediate_dir metadata_dir"
    local k v rel line
    while IFS= read -r line; do
        k="${line%%:*}"
        case " ${keys} fastq_r1 fastq_r2 fastq_dir " in *" ${k} "*) continue ;; esac
        v="$(yaml_get "${CONFIG_FILE}" "${k}")"
        case "${v}" in *.R|*.r|*.py|*.sh) keys="${keys} ${k}" ;; esac
    done < <(grep -E '^[A-Za-z0-9_]+:' "${CONFIG_FILE}")
    for k in ${keys}; do
        v="$(yaml_get "${CONFIG_FILE}" "${k}")"
        is_unset "${v}" && continue
        if is_abs "${v}"; then
            rel="$(project_relative "$(resolve_config_path "${v}")")"
            [[ -n "${rel}" ]] || fail "${k}: ${v} is outside the project folder, which the pipeline cannot see.
Copy it into ${PROJECT_DIR} and set ${k} to its name there (e.g. ${k}: $(basename "${v}"))"
            say "  ${k}: ${v} -> ${rel} (relative to the project folder)"
            yaml_set "${CONFIG_FILE}" "${k}" "${rel}"
        fi
    done
    v="$(yaml_get "${CONFIG_FILE}" samplesheet)"
    if ! is_unset "${v}" && [[ ! -e "$(resolve_config_path "${v}")" ]]; then
        fail "samplesheet: ${v} not found in the project folder. Create it (see examples/samplesheet_example.csv)."
    fi
}

# ---------------------------------------------------------------------------
# Project identity / container naming
# ---------------------------------------------------------------------------
container_name() {
    local leaf h
    leaf="$(basename "$1" | sed 's/[^A-Za-z0-9_.-]/_/g')"
    h="$(printf '%s' "$1" | cksum | awk '{ printf "%08x", $1 }')"
    printf '%s_%s_%s\n' "${PIPELINE_NAME}" "${leaf}" "${h}"
}
is_running() { [[ "$(docker ps --filter "name=^/$1\$" --format '{{.Names}}' 2>/dev/null)" == "$1" ]]; }

project_from_args() {
    PROJECT_DIR="$(pwd -P)"
    CONFIG_NAME="config.yaml"
    if [[ -n "${arg1}" ]]; then
        if [[ -d "${arg1}" ]]; then PROJECT_DIR="$(abspath "${arg1}")"
        elif [[ -f "${arg1}" ]]; then
            PROJECT_DIR="$(cd -P "$(dirname "${arg1}")" && pwd)"; CONFIG_NAME="$(basename "${arg1}")"
        else fail "${arg1} not found (expected a config .yaml file or a project folder)."; fi
    fi
    CONFIG_FILE="${PROJECT_DIR}/${CONFIG_NAME}"
}

# ---------------------------------------------------------------------------
# Resources: one budget per run, enforced as the container's hard limits.
#
#   memory = 75% of the computer's RAM, but always leaving 3 GB to macOS
#            (a plain percentage starves small laptops: 80% of 8 GB leaves 1.6 GB)
#   CPUs   = all logical cores but 1 (but 2 from 8 cores up), so the desktop
#            stays responsive; heavy CPU use only slows a computer, it can't harm it
#   ... and never more than Docker Desktop itself may use (Settings > Resources).
#   config.yaml max_cpus / max_memory_gb replace the default (up to Docker's limit).
#
# docker run --memory makes the budget a hard ceiling (a step that outgrows it
# is killed inside the container and retried, the Mac never starts swapping);
# --cpuset-cpus makes tools that use "all cores" (cutadapt --cores 0) see
# only the budget's cores.
# ---------------------------------------------------------------------------
requirement() {   # key default: the value as written (integer or decimal)
    local v num_re='^[0-9]+(\.[0-9]+)?$'
    v="$(grep -E "^[[:space:]]*$1[[:space:]]*=" "${PIPELINE_DIR}/docker/requirements.conf" 2>/dev/null | head -1 | sed -E 's/.*=[[:space:]]*([0-9]+(\.[0-9]+)?).*/\1/')"
    [[ "${v}" =~ ${num_re} ]] && echo "${v}" || echo "$2"
}
# key default_gb: the same value in whole MB (GB x 1024, rounded)
requirement_mb() {
    awk -v g="$(requirement "$1" "$2")" 'BEGIN { printf "%d", g * 1024 + 0.5 }'
}

# HOST_MEM_MB HOST_CPUS HOST_FREE_MB
host_resources() {
    HOST_MEM_MB=0; HOST_CPUS=1; HOST_FREE_MB=0
    if is_mac; then
        HOST_MEM_MB=$(( $(sysctl -n hw.memsize 2>/dev/null || echo 0) / 1048576 ))
        HOST_CPUS=$(sysctl -n hw.logicalcpu 2>/dev/null || echo 1)
        # free + inactive + speculative pages: what macOS can hand out without swapping
        HOST_FREE_MB=$(vm_stat 2>/dev/null | awk '
            /page size of/ { ps = $8 }
            /^Pages (free|inactive|speculative):/ { gsub(/\./, "", $NF); n += $NF }
            END { if (ps > 0) printf "%d", n * ps / 1048576; else print 0 }')
    else
        HOST_MEM_MB=$(awk '/^MemTotal:/ { m = $2 } END { printf "%d", m / 1024 }' /proc/meminfo 2>/dev/null || echo 0)
        HOST_CPUS=$(nproc 2>/dev/null || echo 1)
        HOST_FREE_MB=$(awk '/^MemAvailable:/ { m = $2 } END { printf "%d", m / 1024 }' /proc/meminfo 2>/dev/null || echo 0)
    fi
}

# DOCKER_MEM_MB DOCKER_CPUS (what Docker Desktop's VM has; the host itself on Linux)
docker_resources() {
    local n m
    n="$(docker info --format '{{.NCPU}}' 2>/dev/null)"
    m="$(docker info --format '{{.MemTotal}}' 2>/dev/null)"
    [[ "${n}" =~ ^[0-9]+$ ]] && DOCKER_CPUS=${n} || DOCKER_CPUS=0
    [[ "${m}" =~ ^[0-9]+$ ]] && DOCKER_MEM_MB=$(( m / 1048576 )) || DOCKER_MEM_MB=0
}

# POLICY_MEM_MB POLICY_CPUS: the default budget, from the computer alone.
policy_budget() {
    local a b spare=1
    a=$(( HOST_MEM_MB * 3 / 4 )); b=$(( HOST_MEM_MB - 3072 ))
    POLICY_MEM_MB=$(( a < b ? a : b )); (( POLICY_MEM_MB < 0 )) && POLICY_MEM_MB=0
    (( HOST_CPUS >= 8 )) && spare=2
    POLICY_CPUS=$(( HOST_CPUS - spare )); (( POLICY_CPUS < 1 )) && POLICY_CPUS=1
}

gb() { awk -v m="$1" 'BEGIN { printf "%.1f GB", m / 1024 }'; }

# BUDGET_CPUS BUDGET_MEM_MB BUDGET_REASON for this run, or stops with an explanation.
run_budget() {
    host_resources; docker_resources; policy_budget
    local dock_mem=$(( DOCKER_MEM_MB - 512 )) dock_cpu=${DOCKER_CPUS} want_mem want_cpu min_mem min_cpu
    (( dock_mem > 0 )) || dock_mem=${POLICY_MEM_MB}
    (( dock_cpu > 0 )) || dock_cpu=${POLICY_CPUS}
    want_mem=""; want_cpu=""   # no config.yaml (the script command): the default budget
    if [[ -f "${CONFIG_FILE}" ]]; then
        want_mem="$(yaml_get "${CONFIG_FILE}" max_memory_gb)"
        want_cpu="$(yaml_get "${CONFIG_FILE}" max_cpus)"
    fi
    BUDGET_REASON=""
    local num_re='^[0-9]+(\.[0-9]+)?$'
    if [[ "${want_mem}" =~ ${num_re} ]]; then
        BUDGET_MEM_MB=$(awk -v g="${want_mem}" 'BEGIN { printf "%d", g * 1024 + 0.5 }')
        if (( BUDGET_MEM_MB > dock_mem )); then
            warn "max_memory_gb: ${want_mem} is more than Docker Desktop may use; using $(gb ${dock_mem})."
            BUDGET_MEM_MB=${dock_mem}
        fi
        BUDGET_REASON="max_memory_gb: ${want_mem} in ${CONFIG_NAME}"
    else
        # Capped at what the pipeline can use: more only fills Docker's VM with
        # file cache that macOS gets back slowly (docker/requirements.conf).
        local useful=$(requirement_mb max_useful_memory_gb 0)
        if (( useful > 0 && useful < POLICY_MEM_MB && useful < dock_mem )); then
            BUDGET_MEM_MB=${useful}
            BUDGET_REASON="the most this pipeline can use (docker/requirements.conf)"
        elif (( POLICY_MEM_MB <= dock_mem )); then
            BUDGET_MEM_MB=${POLICY_MEM_MB}
            BUDGET_REASON="75% of this computer's $(gb ${HOST_MEM_MB}), leaving at least 3 GB to the system"
        else
            BUDGET_MEM_MB=${dock_mem}
            BUDGET_REASON="all that Docker Desktop may use: $(gb ${DOCKER_MEM_MB}) of this computer's $(gb ${HOST_MEM_MB})"
        fi
    fi
    if [[ "${want_cpu}" =~ ^[0-9]+$ ]] && (( want_cpu >= 1 )); then
        BUDGET_CPUS=$(( want_cpu < dock_cpu ? want_cpu : dock_cpu ))
        BUDGET_REASON="${BUDGET_REASON}; max_cpus: ${want_cpu}"
    else
        BUDGET_CPUS=$(( POLICY_CPUS < dock_cpu ? POLICY_CPUS : dock_cpu ))
    fi

    min_mem=$(requirement_mb min_memory_gb 2.5)
    min_cpu=$(requirement min_cpus 1)
    if (( BUDGET_MEM_MB < min_mem )); then
        if [[ "${want_mem}" =~ ${num_re} ]]; then
            fail "max_memory_gb: ${want_mem} in ${CONFIG_NAME} is below the $(gb ${min_mem}) this pipeline needs. Raise it or remove the line."
        elif (( POLICY_MEM_MB >= min_mem )); then
            fail "Docker Desktop may only use $(gb ${DOCKER_MEM_MB}), and this pipeline needs about $(gb ${min_mem}) plus Docker's own share.
Give Docker more: Docker Desktop > Settings > Resources > Memory (e.g. $(( POLICY_MEM_MB / 1024 )) GB) > Apply & restart."
        else
            fail "This computer has $(gb ${HOST_MEM_MB}) of memory. Keeping 3 GB for the system leaves less than the $(gb ${min_mem})
this pipeline needs. Run it on a computer with at least $(( (min_mem + 3072 + 1023) / 1024 )) GB, or set max_memory_gb in config.yaml
to use more (the computer may then become slow while it runs)."
        fi
    fi
    (( BUDGET_CPUS < min_cpu )) && warn "only ${BUDGET_CPUS} CPU core(s) for the analysis: it will work, but slowly."
    if (( HOST_FREE_MB > 0 && HOST_FREE_MB < BUDGET_MEM_MB )); then
        warn "only $(gb ${HOST_FREE_MB}) of memory is free right now (other programs use the rest), and the analysis may use
up to $(gb ${BUDGET_MEM_MB}). Closing large programs (e.g. a browser with many tabs) keeps the computer responsive."
    fi
    return 0
}

# setup: report what each run will get, and how to give Docker more if it is the limit.
report_resources() {
    host_resources; docker_resources; policy_budget
    local dock_mem=$(( DOCKER_MEM_MB - 512 )) cpus mem target useful
    echo
    say "This computer: ${HOST_CPUS} CPU cores, $(gb ${HOST_MEM_MB}) memory. Docker may use: ${DOCKER_CPUS} cores, $(gb ${DOCKER_MEM_MB})."
    cpus=$(( POLICY_CPUS < DOCKER_CPUS ? POLICY_CPUS : DOCKER_CPUS ))
    target=${POLICY_MEM_MB}
    useful=$(requirement_mb max_useful_memory_gb 0)
    (( useful > 0 && useful < target )) && target=${useful}
    mem=$(( target < dock_mem ? target : dock_mem ))
    say "Each analysis will use at most ${cpus} CPU cores and $(gb ${mem}) of memory (max_cpus / max_memory_gb in a project's config.yaml change this)."
    if is_mac && (( DOCKER_MEM_MB < target + 512 )); then
        say "The pipeline could use up to $(gb ${target}) (never more than 75% of the memory, leaving at least 3 GB to macOS), but Docker"
        say "Desktop stops it at $(gb ${DOCKER_MEM_MB}). To give it more: Docker Desktop > Settings > Resources > Memory"
        say "(e.g. $(( (target + 512 + 1023) / 1024 )) GB) > Apply & restart. Optional — smaller screens run fine as is."
    fi
    if (( mem < $(requirement_mb min_memory_gb 2.5) )); then
        warn "that is less than the $(gb $(requirement_mb min_memory_gb 2.5)) this pipeline needs: runs will refuse to start."
    fi
}

host_os() {
    if is_mac; then printf 'Docker Desktop on macOS %s (%s)\n' "$(sw_vers -productVersion 2>/dev/null)" "$(uname -m)"
    else printf 'Docker on %s (%s)\n' "$( (. /etc/os-release 2>/dev/null && printf '%s' "${PRETTY_NAME}") || uname -sr)" "$(uname -m)"; fi
}

latest_log() { ls -1 "${PROJECT_DIR}"/metadata/nf_logs/nextflow_*.log 2>/dev/null | sort | tail -1; }

# One line if the project still has the work/ or .nextflow/ folders of an older version.
legacy_note() {
    local names=() list
    [[ -d "${PROJECT_DIR}/work" ]] && names+=("work/")
    [[ -d "${PROJECT_DIR}/.nextflow" ]] && names+=(".nextflow/")
    (( ${#names[@]} > 0 )) || return 0
    list="${names[0]}"
    if (( ${#names[@]} > 1 )); then list="${names[0]} and ${names[1]}"; fi
    say "Note: ${list} in the project folder are from an older version and are no longer used; remove them with: run_crisprscreen.sh clean"
}

# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------
cmd_help() {
    cat <<EOF
CRISPR_SCREEN pipeline (Docker version for Mac)

  run_crisprscreen.sh setup            one-time: Docker Desktop + the pipeline's Docker image
  run_crisprscreen.sh init [FOLDER]    set up a new project folder (default: current folder)
  run_crisprscreen.sh                  run the pipeline for the project in the current folder
  run_crisprscreen.sh run CONFIG.yaml  ... with another config file
  run_crisprscreen.sh status           progress of the current/last run (and why it failed)
  run_crisprscreen.sh logs             follow the live log
  run_crisprscreen.sh stop             stop the running analysis (re-running resumes it)
  run_crisprscreen.sh clean            delete the work folder and resume cache (metadata/work, metadata/.nextflow) now
  run_crisprscreen.sh script utils/X.R [ARGS]  run one of the project's own scripts (R, Python or shell) in the pipeline's software environment

Options: --detach  --yes  --pause  --rebuild
Pipeline folder: ${PIPELINE_DIR}
EOF
    done_ 0
}

cmd_setup() {
    info "CRISPR_SCREEN setup (Docker version for Mac)"
    ensure_docker 1
    repair_line_endings
    ensure_image "${REBUILD}"

    report_resources

    echo
    docker run --rm --platform linux/amd64 "${IMAGE}" bash -c 'echo "cutadapt $(cutadapt --version)"; echo "MAGeCK $(mageck -v 2>&1 | tail -1)"; fastqc --version; echo "MultiQC $(multiqc --version | sed "s/.*version //")"; python3 --version; R --version | head -1; nextflow -version 2>&1 | grep -i version'

    # Put run_crisprscreen.sh on PATH for new terminals (zsh is macOS's default shell).
    local rc line="export PATH=\"${PIPELINE_DIR}:\$PATH\"   # CRISPR_SCREEN pipeline (run_crisprscreen.sh)"
    for rc in "${HOME}/.zshrc" "${HOME}/.bash_profile"; do
        if [[ "${rc}" == *.zshrc || -f "${rc}" ]]; then
            if ! grep -qF "${PIPELINE_DIR}:" "${rc}" 2>/dev/null; then
                printf '\n%s\n' "${line}" >> "${rc}"
                say "Added ${PIPELINE_DIR} to your PATH in ${rc}"
            fi
        fi
    done
    chmod +x "${PIPELINE_DIR}/run_crisprscreen.sh" 2>/dev/null || true
    echo
    info "Setup complete. Next: create a project folder with   run_crisprscreen.sh init <folder>"
    say  "(in a NEW terminal window, run_crisprscreen.sh works from any folder)"
    done_ 0
}

# Lowercase md5 of a file (macOS: md5 -q; Linux: md5sum).
file_md5() {
    if have md5; then md5 -q "$1"; else md5sum "$1" | awk '{ print $1 }'; fi | tr 'A-F' 'a-f'
}

cmd_init() {
    local target="${arg1:-$(pwd)}" repo_root f d
    mkdir -p "${target}" || fail "cannot create ${target}"
    target="$(abspath "${target}")"
    repo_root="$(cd -P "${PIPELINE_DIR}/.." && pwd)"
    case "${target}/" in
        "${repo_root}/"*) fail "Don't put a project inside the pipeline's own folder (${repo_root}). Choose another folder, e.g. ~/Documents/CRISPR_screens/my_screen" ;;
    esac
    info "Setting up a new ${PIPELINE_NAME} project in ${target} ..."
    mkdir -p "${target}"/original_data/fastq \
             "${target}"/metadata "${target}"/results "${target}"/utils "${target}"/examples
    if [[ -f "${target}/config.yaml" ]]; then say "  config.yaml already exists — leaving it alone."
    else cp "${PIPELINE_DIR}/examples/config.yaml" "${target}/config.yaml"; say "  Fetched config.yaml"; fi
    for f in "${PIPELINE_DIR}"/examples/*; do
        [[ -e "${target}/examples/$(basename "${f}")" ]] || cp "${f}" "${target}/examples/"
    done
    say "  Fetched examples/ (samplesheet_example.csv, library_example.csv, contrasts_example.csv)"

    # Default barcode-position plot for utils/: copied once, md5 recorded (md5sum format).
    local u_src="${PIPELINE_DIR}/bin/plot_barcode_position.R" u_dst="${target}/utils/plot_barcode_position.R"
    local u_rel="utils/plot_barcode_position.R" u_md5="${target}/utils/.pipeline_defaults.md5" u_sum tmp
    if [[ -f "${u_src}" ]]; then
        if [[ -e "${u_dst}" ]]; then say "  ${u_rel} already exists — leaving it alone."
        else
            cp "${u_src}" "${u_dst}"
            u_sum="$(file_md5 "${u_dst}")"
            tmp="${u_md5}.tmp.$$"
            { [[ -f "${u_md5}" ]] && awk -v rel="${u_rel}" 'substr($0, 35) != rel' "${u_md5}"
              printf '%s  %s\n' "${u_sum}" "${u_rel}"; } > "${tmp}" \
                && cat "${tmp}" > "${u_md5}" && rm -f "${tmp}"
            say "  Fetched ${u_rel} (default barcode-position plot; md5 recorded in utils/.pipeline_defaults.md5)"
        fi
    fi

    # AI-agent settings (.claude/), copied file by file; existing files are never overwritten.
    local c_src="${PIPELINE_DIR}/agent/.claude" c_n=0 rel
    if [[ -d "${c_src}" ]]; then
        while IFS= read -r rel; do
            rel="${rel#./}"
            [[ -n "${rel}" && ! -e "${target}/.claude/${rel}" ]] || continue
            mkdir -p "$(dirname "${target}/.claude/${rel}")"
            cp "${c_src}/${rel}" "${target}/.claude/${rel}"
            c_n=$(( c_n + 1 ))
        done <<< "$(cd "${c_src}" && find . -type f)"
        if (( c_n > 0 )); then say "  Fetched .claude/ (${c_n} files for AI agents; existing files left alone)"; fi
    fi

    # Agent instructions, with this computer's pipeline location filled in.
    for f in AGENTS.md CLAUDE.md; do
        if [[ -f "${target}/${f}" ]]; then say "  ${f} already exists — leaving it alone."; continue; fi
        awk -v p="${PIPELINE_DIR}" '{ gsub(/\{\{PIPELINE_DIR\}\}/, p); print }' "${PIPELINE_DIR}/agent/${f}" > "${target}/${f}"
        say "  Fetched ${f}"
    done

    # Launchers inside the project (rewritten each time: they hold the pipeline's location).
    # The terminal launcher lives in utils/ (the project root keeps only the double-click RUN_ANALYSIS.command).
    local q; q="$(printf '%q' "${PIPELINE_DIR}/run_crisprscreen.sh")"
    mkdir -p "${target}/utils"
    cat > "${target}/utils/run_crisprscreen.sh" <<EOF
#!/bin/bash
# Runs the CRISPR_SCREEN pipeline (Docker version) for THIS project folder.
# Same commands as the pipeline's own launcher: (no argument) = run, status, logs, stop, clean, script utils/X.R [ARGS].
# Generated by 'run_crisprscreen.sh init'; the pipeline lives in:
#   ${PIPELINE_DIR}
cd "\$(dirname "\$0")/.." && exec bash ${q} "\$@"
EOF
    cat > "${target}/RUN_ANALYSIS.command" <<EOF
#!/bin/bash
# Double-click to run (or resume) the CRISPR screen analysis of this folder.
cd "\$(dirname "\$0")" && exec bash ${q} run --pause
EOF
    chmod +x "${target}/utils/run_crisprscreen.sh" "${target}/RUN_ANALYSIS.command"
    # An older init wrote the terminal launcher to the project root: remove it only if it is ours (marker line).
    if [[ -f "${target}/run_crisprscreen.sh" ]] && grep -qF "# Runs the CRISPR_SCREEN pipeline (Docker version) for THIS project folder." "${target}/run_crisprscreen.sh"; then
        rm -f "${target}/run_crisprscreen.sh"
        say "Removed the old run_crisprscreen.sh from the project root (now utils/run_crisprscreen.sh)"
    fi
    say "  Created utils/run_crisprscreen.sh (terminal) and RUN_ANALYSIS.command (double-click)"

    # .gitattributes: LF line endings for the scripts, so a utils/ copy stays
    # byte-identical to the pipeline's default on every OS. Never overwritten.
    if [[ ! -e "${target}/.gitattributes" ]]; then
        cat > "${target}/.gitattributes" <<'EOF'
# Keep scripts byte-identical to the pipeline's copies on every OS: utils/
# replacements are compared with the pipeline's default scripts, and bash
# scripts need LF line endings.
utils/** text eol=lf
*.bat text eol=crlf
*.ps1 text eol=crlf
*.sh text eol=lf
*.R text eol=lf
*.py text eol=lf
EOF
        say "  Fetched .gitattributes (LF line endings for scripts)"
    fi

    if [[ -d "${target}/.git" ]]; then
        say "  Already a git repository — leaving it alone."
    elif have git && git --version >/dev/null 2>&1; then
        [[ -f "${target}/.gitignore" ]] || cp "${PIPELINE_DIR}/examples/project.gitignore" "${target}/.gitignore"
        local paths="config.yaml examples .gitignore .gitattributes AGENTS.md CLAUDE.md RUN_ANALYSIS.command utils"
        if [[ -d "${target}/.claude" ]]; then paths="${paths} .claude"; fi
        ( cd "${target}" && git init -q . \
            && git add ${paths} 2>/dev/null \
            && git -c user.name=crisprscreen -c user.email=crisprscreen@localhost commit -q -m "Initialize ${PIPELINE_NAME} project" ) >/dev/null 2>&1 || true
        say "  Initialized git repository (raw data / results / metadata excluded — see .gitignore; utils/ and .claude/ are included in the first commit)"
    else
        say "  (git not available — skipped version control; optional)"
    fi

    echo
    say "Next steps:"
    say "  1. Put your input tables in this folder (start from the files in examples/):"
    say "       samplesheet.csv, library.csv, and — when you are ready for mageck test — contrasts.csv"
    say "  2. Edit config.yaml: project_name, and fastq_r1 (multiplexed FASTQ) or fastq_dir"
    say "     (folder of already-demultiplexed FASTQs). Your FASTQ files stay where they are: put"
    say "     their full path in config.yaml and keep it stable."
    say "  3. Run:  double-click RUN_ANALYSIS.command   (or type  ./utils/run_crisprscreen.sh  in this folder)"
    done_ 0
}

cmd_run() {
    project_from_args
    [[ -f "${CONFIG_FILE}" ]] || fail "${CONFIG_NAME} not found in ${PROJECT_DIR}. Set the project up first:  run_crisprscreen.sh init"
    local name; name="$(container_name "${PROJECT_DIR}")"
    ensure_docker 0
    is_running "${name}" && fail "this project's analysis is already running. Watch it with:  run_crisprscreen.sh logs   (or stop it: run_crisprscreen.sh stop)"
    repair_line_endings
    run_budget

    info "Checking the FASTQ paths in ${CONFIG_FILE} (they stay where they are)..."
    input_args=()
    resolve_raw_fastq fastq_dir
    if [[ -n "${RAW_PATH}" ]]; then
        local r1; r1="$(yaml_get "${CONFIG_FILE}" fastq_r1)"
        is_unset "${r1}" || say "  fastq_dir is set, so the samples are treated as already demultiplexed (fastq_r1 is ignored)."
    else
        is_unset "$(yaml_get "${CONFIG_FILE}" fastq_r1)" \
            && fail "config.yaml needs fastq_r1 (a multiplexed FASTQ) or fastq_dir (a folder of already-demultiplexed FASTQs)."
        resolve_raw_fastq fastq_r1
        resolve_raw_fastq fastq_r2
    fi
    check_config_paths

    ensure_image 0

    local ts commit="" modified="false" image_id docker_version
    ts="$(date +%Y%m%d_%H%M%S)"
    if have git && commit="$(git -c safe.directory='*' -C "${PIPELINE_DIR}" rev-parse --short HEAD 2>/dev/null)"; then
        [[ -n "$(git -c safe.directory='*' -C "${PIPELINE_DIR}" status --porcelain -- . 2>/dev/null)" ]] && modified="true"
    else
        commit=""
    fi
    image_id="$(docker image inspect "${IMAGE}" --format '{{.Id}}' 2>/dev/null)"
    docker_version="$(docker version --format '{{.Server.Version}}' 2>/dev/null)"

    local user_args=()
    # Linux hosts: write files as the user, not root (Docker Desktop on Mac maps this itself).
    is_mac || user_args=(--user "$(id -u):$(id -g)" -e HOME=/tmp)

    legacy_note

    # Mark the run as started now: until the container rewrites it, `status`
    # (and an agent polling it) would otherwise still see the previous run's outcome.
    mkdir -p "${PROJECT_DIR}/metadata/nf_logs"
    echo "running" > "${PROJECT_DIR}/metadata/nf_logs/last_run.status"

    docker run -d --rm --init --name "${name}" --platform linux/amd64 \
        --label "crisprscreen.project=${PROJECT_DIR}" \
        --memory "${BUDGET_MEM_MB}m" --memory-swap "${BUDGET_MEM_MB}m" --cpuset-cpus "0-$(( BUDGET_CPUS - 1 ))" \
        -e "CRISPR_BUDGET_CPUS=${BUDGET_CPUS}" -e "CRISPR_BUDGET_MEMORY_MB=${BUDGET_MEM_MB}" -e "CRISPR_BUDGET_REASON=${BUDGET_REASON}" \
        --mount "type=bind,source=${PROJECT_DIR},target=/project" \
        --mount "type=bind,source=${PIPELINE_DIR},target=/pipeline,readonly" \
        ${input_args[@]+"${input_args[@]}"} \
        ${user_args[@]+"${user_args[@]}"} \
        -e "CRISPR_IMAGE=${IMAGE}" \
        -e "CRISPR_IMAGE_ID=${image_id}" \
        -e "CRISPR_CONTAINER_ENGINE=docker" \
        -e "CRISPR_CONTAINER_ENGINE_VERSION=${docker_version}" \
        -e "CRISPR_HOST_OS=$(host_os)" \
        -e "CRISPR_PIPELINE_COMMIT=${commit}" \
        -e "CRISPR_PIPELINE_MODIFIED=${modified}" \
        -e "CRISPR_PIPELINE_HOST_DIR=${PIPELINE_DIR}" \
        -e "CRISPR_PROJECT_HOST_DIR=${PROJECT_DIR}" \
        "${IMAGE}" bash /pipeline/docker/run_in_container.sh "${CONFIG_NAME}" "${ts}" >/dev/null \
        || { echo "failed to start" > "${PROJECT_DIR}/metadata/nf_logs/last_run.status"
             fail "Docker could not start the pipeline container (see the message above)."; }

    # Keep the Mac from sleeping (idle sleep) for as long as the container runs.
    if is_mac && have caffeinate; then
        nohup caffeinate -i bash -c "while docker ps -q --filter 'name=^/${name}\$' 2>/dev/null | grep -q .; do sleep 60; done" \
            >/dev/null 2>&1 &
    fi

    echo
    info "The analysis is running inside Docker."
    say  "  Log file : ${PROJECT_DIR}/metadata/nf_logs/nextflow_${ts}.log"
    say  "  Uses at most ${BUDGET_CPUS} CPU cores and $(gb ${BUDGET_MEM_MB}) memory (${BUDGET_REASON})."
    say  "  Progress : run_crisprscreen.sh status      Stop: run_crisprscreen.sh stop"
    say  "  Closing this window does NOT stop the analysis. Keep the Mac plugged in and the lid"
    say  "  open until it finishes (it is kept from idle-sleeping while the analysis runs)."
    (( DETACH )) && done_ 0
    say  "  (Press Ctrl+C to stop watching; the analysis keeps running.)"
    echo
    follow_run "${name}"
}

follow_run() {
    local name="$1" st=""
    trap 'echo; say "Stopped watching; the analysis keeps running (run_crisprscreen.sh status)."; done_ 0' INT
    docker logs -f "${name}" 2>&1
    trap - INT
    sleep 1
    [[ -f "${PROJECT_DIR}/metadata/nf_logs/last_run.status" ]] && st="$(cat "${PROJECT_DIR}/metadata/nf_logs/last_run.status")"
    echo
    if [[ "${st}" == succeeded ]]; then
        info "Pipeline finished successfully."
        say  "  Results : ${PROJECT_DIR}/results/   (start with results/counts/qc/report/*_qc_report.html)"
        say  "  Methods : ${PROJECT_DIR}/METHODS.md"
        done_ 0
    elif [[ "${st}" == running ]] && is_running "${name}"; then
        say "Stopped watching; the analysis is still running."
        done_ 0
    else
        printf '%sPipeline did not finish: %s%s\n' "${C_ERR}" "${st}" "${C_OFF}"
        say "  See:  run_crisprscreen.sh status   (it names the failed step and where its messages were saved)"
        say "  Fix the problem, then run again: finished steps are reused, not redone."
        done_ 1
    fi
}

cmd_status() {
    project_from_args
    legacy_note
    local name running=0 log st trace ft
    name="$(container_name "${PROJECT_DIR}")"
    if have docker && docker info >/dev/null 2>&1 && is_running "${name}"; then running=1; fi
    log="$(latest_log)"
    [[ -n "${log}" ]] || { say "No run yet in ${PROJECT_DIR}."; done_ 0; }
    if (( running )); then info "Status: RUNNING"
    elif [[ -f "${PROJECT_DIR}/metadata/nf_logs/last_run.status" ]]; then
        st="$(cat "${PROJECT_DIR}/metadata/nf_logs/last_run.status")"
        [[ "${st}" == running ]] && st="interrupted (Docker or the computer stopped mid-run; run again to resume)"
        info "Status: last run ${st}"
    fi
    say "Log: ${log}"
    trace="${PROJECT_DIR}/metadata/pipeline_info/execution_trace.txt"
    if [[ -s "${trace}" ]]; then
        echo
        say "Steps (this run):"
        awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
            { n[$c["status"]]++; if ($c["status"] == "FAILED") f[++k] = $c["name"] " (exit " $c["exit"] ")" }
            END { for (s in n) printf "  %-10s %d\n", s, n[s]; for (i = 1; i <= k; i++) print "  FAILED: " f[i] }' "${trace}"
    fi
    ft="${PROJECT_DIR}/metadata/failed_tasks/$(basename "${log}" .log | sed 's/^nextflow_//')"
    [[ -d "${ft}" ]] && { echo; say "Files of the failed step(s): ${ft}"; }
    echo
    say "Last lines of the log:"
    tail -n 25 "${log}" | sed 's/^/  /'
    done_ 0
}

cmd_logs() {
    project_from_args
    local name log; name="$(container_name "${PROJECT_DIR}")"
    ensure_docker 0
    if is_running "${name}"; then say "(Ctrl+C stops watching; the analysis keeps running.)"; follow_run "${name}"; fi
    log="$(latest_log)"
    [[ -n "${log}" ]] || { say "No run yet."; done_ 0; }
    say "Not running. Last log: ${log}"
    tail -n 60 "${log}"
    done_ 0
}

cmd_stop() {
    project_from_args
    local name; name="$(container_name "${PROJECT_DIR}")"
    ensure_docker 0
    is_running "${name}" || { say "Nothing is running for this project."; done_ 0; }
    info "Stopping the analysis (finished steps are kept; running it again resumes)..."
    docker stop -t 60 "${name}" >/dev/null 2>&1
    echo "stopped by the user" > "${PROJECT_DIR}/metadata/nf_logs/last_run.status"
    say "Stopped."
    done_ 0
}

cmd_clean() {
    project_from_args
    local name; name="$(container_name "${PROJECT_DIR}")"
    local dirs=() d
    for d in metadata/work metadata/.nextflow work .nextflow; do
        [[ -d "${PROJECT_DIR}/${d}" ]] && dirs+=("${d}")
    done
    (( ${#dirs[@]} > 0 )) || { say "Nothing to clean (no metadata/work, metadata/.nextflow, work/ or .nextflow/ folder)."; done_ 0; }
    ensure_docker 0
    is_running "${name}" && fail "the analysis is running — stop it first (run_crisprscreen.sh stop)."
    say "This deletes:"
    for d in "${dirs[@]}"; do say "  ${PROJECT_DIR}/${d}"; done
    say "They are Nextflow's working copies of every step and its resume cache."
    say "results/, original_data/ and METHODS.md are NOT touched, but the next run will start from scratch."
    confirm "Delete them?" || { say "Nothing deleted."; done_ 0; }
    for d in "${dirs[@]}"; do rm -rf "${PROJECT_DIR}/${d}"; done
    say "Deleted."
    done_ 0
}

# run_crisprscreen.sh script utils/X.R [ARGS]: one of the project's own scripts, run in the
# pipeline's software environment. The project is the current folder; all ARGS go to the script.
cmd_script() {
    local script_arg="${arg1}" src full rel base stem name ts logrel image_id rc
    local user_args=() script_args=()
    if [[ -z "${script_arg}" ]]; then
        fail "usage: run_crisprscreen.sh script utils/<name>.<R|py|sh> [args...]   (run from the project folder)"
    fi
    PROJECT_DIR="$(pwd -P)"
    CONFIG_NAME="config.yaml"
    CONFIG_FILE="${PROJECT_DIR}/${CONFIG_NAME}"

    # The script: relative to the project folder (or absolute), inside it, with an extension the container runs.
    src="${script_arg}"
    case "${src}" in
        "~")   src="${HOME}" ;;
        "~/"*) src="${HOME}/${src#\~/}" ;;
        /*)    ;;
        *)     src="${PROJECT_DIR}/${src}" ;;
    esac
    [[ -f "${src}" ]] || fail "${script_arg} not found, or not a regular file (looked for ${src})."
    full="$(abspath "${src}")"
    rel="$(project_relative "${full}")"
    [[ -n "${rel}" ]] || fail "The script must be inside the project folder (${PROJECT_DIR})."
    base="$(basename "${full}")"
    case "${base}" in
        *.R|*.r|*.py|*.sh) ;;
        *) fail "${script_arg}: unsupported extension (use .R, .py or .sh)" ;;
    esac
    stem="${base%.*}"
    if (( ${#positional[@]} > 1 )); then script_args=("${positional[@]:1}"); fi

    ensure_docker 0
    name="$(container_name "${PROJECT_DIR}")"
    is_running "${name}" && warn "the project's analysis is running: this script shares the computer's memory with it."
    repair_line_endings
    run_budget
    ensure_image 0

    is_mac || user_args=(--user "$(id -u):$(id -g)" -e HOME=/tmp)
    ts="$(date +%Y%m%d_%H%M%S)"
    logrel="metadata/script_logs/${stem}_${ts}.log"
    mkdir -p "${PROJECT_DIR}/metadata/script_logs"
    image_id="$(docker image inspect "${IMAGE}" --format '{{.Id}}' 2>/dev/null)"

    info "Running ${rel} in the pipeline's software environment (at most ${BUDGET_CPUS} CPU cores and $(gb ${BUDGET_MEM_MB}) memory). Log: ${PROJECT_DIR}/${logrel}"
    docker run --rm --init --platform linux/amd64 \
        --label "crisprscreen.project=${PROJECT_DIR}" \
        --memory "${BUDGET_MEM_MB}m" --memory-swap "${BUDGET_MEM_MB}m" --cpuset-cpus "0-$(( BUDGET_CPUS - 1 ))" \
        -e "CRISPR_BUDGET_CPUS=${BUDGET_CPUS}" -e "CRISPR_BUDGET_MEMORY_MB=${BUDGET_MEM_MB}" -e "CRISPR_BUDGET_REASON=${BUDGET_REASON}" \
        --mount "type=bind,source=${PROJECT_DIR},target=/project" \
        --mount "type=bind,source=${PIPELINE_DIR},target=/pipeline,readonly" \
        -w /project \
        ${user_args[@]+"${user_args[@]}"} \
        -e "CRISPR_IMAGE=${IMAGE}" \
        -e "CRISPR_IMAGE_ID=${image_id}" \
        -e "CRISPR_PROJECT_HOST_DIR=${PROJECT_DIR}" \
        "${IMAGE}" bash /pipeline/docker/run_script_in_container.sh "${rel}" "${logrel}" \
        ${script_args[@]+"${script_args[@]}"}
    rc=$?

    echo
    say "Script finished (exit code ${rc}). Log: ${PROJECT_DIR}/${logrel}"
    done_ "${rc}"
}

case "${COMMAND}" in
    setup)  cmd_setup ;;
    init)   cmd_init ;;
    run)    cmd_run ;;
    status) cmd_status ;;
    logs)   cmd_logs ;;
    stop)   cmd_stop ;;
    clean)  cmd_clean ;;
    script) cmd_script ;;
    *)      cmd_help ;;
esac
