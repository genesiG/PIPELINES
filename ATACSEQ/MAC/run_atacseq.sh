#!/bin/bash
# run_atacseq.sh — ATAC-seq pipeline launcher, Docker version for Mac.
# (Also works on Linux with Docker Engine.) Written for macOS's built-in
# bash 3.2: no associative arrays, no ${x,,}, no readlink -f, BSD sed/stat.
#
#   run_atacseq.sh setup            one-time: check/install Docker Desktop, download the
#                                     pipeline's Docker image, put this command on PATH
#   run_atacseq.sh init [FOLDER]    set up a new project folder (default: the current folder)
#   run_atacseq.sh [run] [CONFIG]   run the pipeline for the project in the current folder
#                                     (CONFIG defaults to config.yaml)
#   run_atacseq.sh status           is it running? how far along? did it fail, and where?
#   run_atacseq.sh logs             follow the live log of the running analysis
#   run_atacseq.sh stop             stop a running analysis (re-running resumes it)
#   run_atacseq.sh clean            delete the work folder, resume cache and work volume (metadata/work, metadata/.nextflow, Docker volume) now
#   run_atacseq.sh script utils/X.R [ARGS]  run one of the project's own scripts (R, Python or shell) in the pipeline's software environment
#   run_atacseq.sh help
#   mem-plan [--total-ram-gb N] [--docker-gb N]   hidden (maintainers, tests): the Docker memory decision for hypothetical inputs
#   scan-paths [CONFIG]                           hidden (maintainers, tests): the folders a run would mount
#   selftest                                      hidden (maintainers, tests): PASS/FAIL checks of the pure helpers (no Docker needed)
#
# Options: --detach (run: start and return, don't follow the log), --yes
# (answer yes to questions: for AI agents / scripts), --pause (wait for Enter
# before closing; used by the double-clickable RUN_ANALYSIS.command).
#
# What `run` does (the Docker counterpart of the cluster's run_atacseq.sh):
#   1. Checks the input files named in config.yaml and samplesheet.csv. The inputs
#      stay where they are: never moved, copied or renamed. An absolute path such as
#      /Users/me/x.fastq.gz or /Volumes/LabNAS/y.bed is visible inside the container at
#      the same path, read-only: the top-level folder holding it is mounted (/Users,
#      /Volumes/<name>, or e.g. /data). Network paths (smb://..., //server/share) are
#      refused: mount the share in Finder and use its /Volumes/<name>/... path. Docker
#      Desktop must be allowed to share those folders (Settings > Resources > File sharing).
#   2. Checks the free disk space on the project's volume (warns; stops only below
#      disk_hard_min_gb in docker/requirements.conf).
#   3. Checks Docker Desktop's memory against this Mac (memory_plan). If Docker is set too
#      low for the analysis it says how to raise it. It never changes Docker's settings.
#   4. Makes sure the pipeline's image is downloaded, then starts one container that runs
#      Nextflow and every step (docker/run_in_container.sh), with the project folder mounted
#      at /project and this pipeline folder (read-only) at /pipeline, then follows its log.
#      Closing the window does NOT stop the analysis.

set -u -o pipefail

PIPELINE_NAME="atacseq"
IMAGE_FILE=""                # set below: docker/image.txt
IMAGE=""                     # read from docker/image.txt by the commands that need it
NUM_RE='^[0-9]+(\.[0-9]+)?$'
# Paths in config.yaml and CSV files (grep -E). An absolute (or ~/) path starts after the
# start of a line, a space, a comma or a quote; a network path is a scheme or a leading //.
TOKEN_RE="(^|[[:space:],'\"])(~/|/)[^[:space:],'\"]*"
NET_RE="(smb|cifs|afp|nfs)://[^[:space:],'\"]*"

# This script's real folder, following symlinks (e.g. one in ~/bin).
src="${BASH_SOURCE[0]}"
while [ -h "${src}" ]; do
    dir="$(cd -P "$(dirname "${src}")" && pwd)"
    src="$(readlink "${src}")"
    case "${src}" in /*) ;; *) src="${dir}/${src}" ;; esac
done
PIPELINE_DIR="$(cd -P "$(dirname "${src}")" && pwd)"
IMAGE_FILE="${PIPELINE_DIR}/docker/image.txt"

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
PAUSE=0; YES=0; DETACH=0
positional=()
pass_rest=0   # after the command `script`, every argument belongs to the script
for a in "$@"; do
    if (( pass_rest )); then positional+=("${a}"); continue; fi
    case "${a}" in
        --pause|-pause)     PAUSE=1 ;;
        --yes|-yes|-y)      YES=1 ;;
        --detach|-detach)   DETACH=1 ;;
        *)                  positional+=("${a}")
                            if [[ "${a}" == script ]] && (( ${#positional[@]} == 1 )); then pass_rest=1; fi ;;
    esac
done
COMMAND="run"
if (( ${#positional[@]} > 0 )); then
    case "${positional[0]}" in
        setup|init|run|status|logs|stop|clean|script|help|mem-plan|scan-paths|selftest|-h|--help)
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

# Size of a file in bytes (GNU stat first; macOS stat has no -c).
file_size() {   # file
    stat -c %s "$1" 2>/dev/null || stat -f %z "$1" 2>/dev/null || echo 0
}
# bytes -> "1.2 GB" (from 1 GB up) or "350 MB"
fmt_bytes() {   # bytes
    awk -v b="$1" 'BEGIN { if (b >= 1073741824) printf "%.1f GB", b / 1073741824; else printf "%.0f MB", b / 1048576 }'
}
# MB -> "x.y GB"
gb() { awk -v m="$1" 'BEGIN { printf "%.1f GB", m / 1024 }'; }
# MB -> whole GB, rounded to the nearest (macOS and Docker report sizes a little off).
round_gb() { awk -v m="$1" 'BEGIN { printf "%d", m / 1024 + 0.5 }'; }
# "~/x" -> "$HOME/x"
expand_tilde() {
    case "$1" in
        "~/"*) printf '%s\n' "${HOME}/${1#\~/}" ;;
        *)     printf '%s\n' "$1" ;;
    esac
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
        fail "Docker is not installed. Run:  bash \"${PIPELINE_DIR}/run_atacseq.sh\" setup"
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

# The image tag: the one line of docker/image.txt (ATAC_IMAGE_OVERRIDE replaces it).
get_image_tag() {
    if [[ -n "${ATAC_IMAGE_OVERRIDE:-}" ]]; then
        IMAGE="$(printf '%s' "${ATAC_IMAGE_OVERRIDE}" | tr -d ' \t\r')"
        return 0
    fi
    [[ -f "${IMAGE_FILE}" ]] || fail "docker/image.txt is missing from the pipeline folder (${PIPELINE_DIR}). It names the pipeline's software image."
    IMAGE="$(tr -d '\r' < "${IMAGE_FILE}" | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//' | grep -v '^#' | grep -v '^$' | head -1)"
    [[ -n "${IMAGE}" ]] || fail "docker/image.txt is empty. It must hold the image name, e.g. ghcr.io/genesig/omics:2026.10.2"
}

# Downloads the pipeline's image (setup always does; run only when it is missing).
pull_image() {
    info "Downloading the pipeline's software image (${IMAGE})."
    say  "This happens once (and again only when docker/image.txt names a new version)."
    say  "It is several GB and can take a while. Keep the Mac plugged in and online."
    if have caffeinate; then caffeinate -i docker pull --platform linux/amd64 "${IMAGE}"
    else docker pull --platform linux/amd64 "${IMAGE}"; fi \
        || fail "Downloading the pipeline's image failed (see the messages above). Check the internet connection and free disk space, then run this again."
    info "Image ${IMAGE} is ready."
}

ensure_image() {
    ensure_rosetta
    docker image inspect "${IMAGE}" >/dev/null 2>&1 || pull_image
}

# Bash scripts must have LF line endings inside the container; repair a copy
# that came through Windows.
repair_line_endings() {
    local f
    for f in "${PIPELINE_DIR}"/bin/*.sh "${PIPELINE_DIR}"/docker/*.sh; do
        [[ -f "${f}" ]] || continue
        if grep -q $'\r' "${f}"; then
            tr -d '\r' < "${f}" > "${f}.tmp" && cat "${f}.tmp" > "${f}" && rm -f "${f}.tmp"
            warn "converted $(basename "${f}") to Unix line endings (needed inside the container)"
        fi
    done
}

# ---------------------------------------------------------------------------
# config.yaml: flat keys only (as on the cluster), so simple line reads do.
# ---------------------------------------------------------------------------
yaml_get() {   # file key
    grep -E "^$2:" "$1" | head -1 \
        | sed -E "s/^$2:[[:space:]]*//; s/[[:space:]]+#.*\$//; s/[[:space:]]*\$//; s/^['\"]//; s/['\"]\$//" | tr -d '\r'
}
is_unset() { [[ -z "$1" || "$1" == null || "$1" == "~" ]]; }

# A path from config.yaml or a CSV as an absolute path (relative = relative to the project).
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
mac_docker_shared() {   # physical path (the shared folder itself, or anything under it)
    case "$1" in
        /Users|/Volumes|/private|/tmp|/var/folders) return 0 ;;
        /Users/*|/Volumes/*|/private/*|/tmp/*|/var/folders/*) return 0 ;;
    esac
    return 1
}

# ---------------------------------------------------------------------------
# Input paths. Inputs are never moved. An absolute path named in config.yaml, in
# the project's CSV files or in a script's arguments is visible inside the
# container at the same path: the top-level folder holding it is mounted
# read-only (/Users, /Volumes/<name>, or /data and the like). Network paths
# cannot be mounted and are refused.
# ---------------------------------------------------------------------------
MOUNT_ROOTS=()     # folders to mount read-only, at the same path (set by scan_inputs / add_root_for)
MOUNT_ARGS=()      # docker run arguments for MOUNT_ROOTS (set by build_mounts)
MOUNTS_TEXT=""     # ATAC_MOUNTS_HOST: lines like /Users=/Users (set by build_mounts)
EXTRA_ENV=()       # optional test settings passed on to the container (add_optional_env)
WARNED="|"         # messages already printed (no repeats)

warn_once() {   # message
    case "${WARNED}" in *"|$1|"*) return 0 ;; esac
    WARNED="${WARNED}$1|"
    warn "$1"
}

# The folder (or file) a mount needs for an absolute path: /Users/x -> /Users,
# /Volumes/NAS/x -> /Volumes/NAS, /data/x -> /data. Prints nothing (and fails)
# when that folder does not exist.
mount_root_of() {   # path
    local p="$1" rest first root
    case "${p}" in
        /Volumes/*) rest="${p#/Volumes/}"; first="${rest%%/*}"; root="/Volumes/${first}" ;;
        /*)         rest="${p#/}";         first="${rest%%/*}"; root="/${first}" ;;
        *) return 1 ;;
    esac
    [[ -n "${first}" && -e "${root}" ]] || return 1
    printf '%s\n' "${root}"
}

add_mount_root() {   # root
    local i
    for (( i = 0; i < ${#MOUNT_ROOTS[@]}; i++ )); do
        [[ "${MOUNT_ROOTS[$i]}" == "$1" ]] && return 0
    done
    MOUNT_ROOTS+=("$1")
}

# Mounts for one absolute path: the folder it is in, and the folder of the real
# file when it is a symlink (its target may live elsewhere).
add_root_for() {   # path
    local p="$1" r real
    if r="$(mount_root_of "${p}")"; then add_mount_root "${r}"; fi
    if [[ -e "${p}" ]]; then
        real="$(resolve_symlinks "${p}")"
        if r="$(mount_root_of "${real}")"; then add_mount_root "${r}"; fi
    fi
}

refuse_network() {   # what, path
    fail "$1 names a network path ($2). Docker Desktop cannot read network shares. Mount the share in Finder (it then appears under /Volumes/<name>), then write the path that starts with /Volumes/<name>/."
}

# One line per path named in a file: "key<TAB>path" (key = the config.yaml key, or "-"
# for a CSV). Comments are removed from config.yaml first. Network paths come out
# too (leading // or a scheme), so the caller can refuse them.
named_paths() {   # file is_config(1/0)
    local file="$1" is_cfg="$2" line key tok
    while IFS= read -r line || [[ -n "${line}" ]]; do
        line="${line%$'\r'}"
        key="-"
        if (( is_cfg )); then
            line="$(printf '%s\n' "${line}" | sed -E 's/(^|[[:space:]])#.*$//')"
            case "${line}" in [A-Za-z_]*:*) key="${line%%:*}" ;; esac
        fi
        while IFS= read -r tok; do
            printf '%s\t%s\n' "${key}" "${tok}"
        done < <(
            printf '%s\n' "${line}" | grep -oE "${NET_RE}"
            printf '%s\n' "${line}" | grep -oE "${TOKEN_RE}" | sed -E "s/^[[:space:],'\"]//"
        )
    done < "${file}"
}

# Every absolute path named in config.yaml and in the project's CSV files (root folder
# only). Refuses network paths, fills MOUNT_ROOTS, and warns about paths that do not exist.
scan_inputs() {
    local files=() f key tok path is_cfg csv i
    MOUNT_ROOTS=()
    [[ -f "${CONFIG_FILE}" ]] && files+=("${CONFIG_FILE}")
    for csv in "${PROJECT_DIR}"/*.csv; do [[ -f "${csv}" ]] && files+=("${csv}"); done
    for (( i = 0; i < ${#files[@]}; i++ )); do
        f="${files[$i]}"
        is_cfg=0; [[ "${f}" == "${CONFIG_FILE}" ]] && is_cfg=1
        while IFS=$'\t' read -r key tok; do
            [[ "${key}" == "-" ]] && key=""
            case "${tok}" in //*|*://*) refuse_network "$(basename "${f}")" "${tok}" ;; esac
            path="$(expand_tilde "${tok}")"
            add_root_for "${path}"
            if [[ "${key}" == genome_index || "${key}" == contamination_index ]]; then
                ls "${path}"*.bt2* >/dev/null 2>&1 \
                    || warn_once "${key}: no ${path##*/}*.bt2* index files found in $(dirname "${path}")"
            elif [[ ! -e "${path}" ]]; then
                warn_once "$(basename "${f}"): ${tok} does not exist on this Mac."
            fi
        done < <(named_paths "${f}" "${is_cfg}")
    done
}

# The reference files of config.yaml that are written relative to the project
# (absolute ones are checked by scan_inputs).
check_config_files() {
    local k v item f
    for k in genome_index contamination_index blacklist chrom_sizes gtf tss_bed; do
        v="$(yaml_get "${CONFIG_FILE}" "${k}")"
        is_unset "${v}" && continue
        while IFS= read -r item; do
            item="$(printf '%s' "${item}" | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//')"
            [[ -n "${item}" ]] || continue
            case "${item}" in /*|"~"*) continue ;; esac
            f="$(resolve_config_path "${item}")"
            if [[ "${k}" == genome_index || "${k}" == contamination_index ]]; then
                ls "${f}"*.bt2* >/dev/null 2>&1 || warn "${k}: no ${item##*/}*.bt2* index files found in $(dirname "${f}")"
            elif [[ ! -e "${f}" ]]; then
                warn "${k}: file not found: ${item} (looked for ${f})"
            fi
        done <<< "$(printf '%s' "${v}" | tr ',' '\n')"
    done
}

# The output folders are written inside the project folder: an absolute path cannot be written.
check_output_dirs() {
    local k v d
    for k in outdir metadata_dir intermediate_dir; do
        v="$(yaml_get "${CONFIG_FILE}" "${k}")"
        is_unset "${v}" && continue
        case "${v}" in
            /*|"~"*)
                case "${k}" in
                    outdir)       d="./results" ;;
                    metadata_dir) d="./metadata" ;;
                    *)            d="./intermediate_files" ;;
                esac
                fail "${k}: ${v} is an absolute path. The results are written into the project folder, so write it relative to the project folder instead, e.g. ${k}: ${d}." ;;
        esac
    done
}

# Every folder the analysis mounts, and the project folder, must be shared with Docker Desktop.
check_mounts_shared() {
    local i r
    for (( i = 0; i < ${#MOUNT_ROOTS[@]}; i++ )); do
        r="${MOUNT_ROOTS[$i]}"
        case "${r}" in *,*) fail "${r} contains a comma, which breaks the container's mount. Rename or move it (no comma in the path), then run again." ;; esac
        if is_mac && ! mac_docker_shared "${r}"; then
            fail "Docker Desktop cannot see ${r}, which an input file is in.
Open Docker Desktop > Settings > Resources > File sharing, add the folder (or move the file under /Users or /Volumes), then run again."
        fi
    done
    case "${PROJECT_DIR}" in *,*) fail "the project folder (${PROJECT_DIR}) contains a comma, which breaks the container's mount. Rename or move it, then run again." ;; esac
    if is_mac && ! mac_docker_shared "${PROJECT_DIR}"; then
        fail "Docker Desktop cannot see the project folder (${PROJECT_DIR}).
Open Docker Desktop > Settings > Resources > File sharing, add the folder (or move the project under /Users or /Volumes), then run again."
    fi
}

# docker run arguments for the mounts (MOUNT_ARGS) and the ATAC_MOUNTS_HOST text (MOUNTS_TEXT).
build_mounts() {
    local i r
    MOUNT_ARGS=(); MOUNTS_TEXT=""
    for (( i = 0; i < ${#MOUNT_ROOTS[@]}; i++ )); do
        r="${MOUNT_ROOTS[$i]}"
        MOUNT_ARGS+=(--mount "type=bind,source=${r},target=${r},readonly")
        MOUNTS_TEXT="${MOUNTS_TEXT}${r}=${r}"$'\n'
    done
    MOUNTS_TEXT="${MOUNTS_TEXT%$'\n'}"
}

samplesheet_path() {
    local v; v="$(yaml_get "${CONFIG_FILE}" samplesheet)"
    is_unset "${v}" && v="samplesheet.csv"
    resolve_config_path "${v}"
}

check_samplesheet() {
    local ss; ss="$(samplesheet_path)"
    [[ -f "${ss}" ]] || fail "the samplesheet (${ss}) is not in the project folder. Copy the template from examples/samplesheet.csv into the project folder, fill it in, and run this again."
}

# The values of the fastq_1 and fastq_2 columns of a CSV, one per line. Quoted
# fields are read simply (no quotes inside a field).
csv_input_values() {   # file
    awk '
    function split_csv(line, out,   n, i, c, f, q) {
        n = 0; f = ""; q = 0
        for (i = 1; i <= length(line); i++) {
            c = substr(line, i, 1)
            if (c == "\"") q = !q
            else if (c == "," && !q) { out[++n] = f; f = "" }
            else f = f c
        }
        out[++n] = f
        return n
    }
    function clean(s) { gsub(/^[ \t"]+|[ \t"]+$/, "", s); return s }
    NR == 1 {
        sub(/\r$/, "")
        nh = split_csv($0, hdr)
        for (j = 1; j <= nh; j++) {
            h = clean(hdr[j])
            if (h == "fastq_1" || h == "fastq_2") want[j] = 1
        }
        next
    }
    {
        sub(/\r$/, "")
        n = split_csv($0, row)
        for (j in want) if (j <= n) { v = clean(row[j]); if (v != "") print v }
    }' "$1"
}

# Free space on the project's volume against the input files named in the samplesheet.
check_disk() {
    local ss v full real total=0 count=0 seen="|" free_kb vol free_b factor res_factor hard_gb need_b work_b
    ss="$(samplesheet_path)"
    while IFS= read -r v; do
        [[ -n "${v}" ]] || continue
        full="$(resolve_config_path "${v}")"
        case "${seen}" in *"|${full}|"*) continue ;; esac
        seen="${seen}${full}|"
        [[ -f "${full}" ]] || continue
        real="$(resolve_symlinks "${full}")"
        total=$(( total + $(file_size "${real}") ))
        count=$(( count + 1 ))
    done < <(csv_input_values "${ss}")
    free_kb="$(df -Pk "${PROJECT_DIR}" 2>/dev/null | awk 'NR == 2 { print $4 }')"
    vol="$(df -Pk "${PROJECT_DIR}" 2>/dev/null | awk 'NR == 2 { print $6 }')"
    free_b=$(( ${free_kb:-0} * 1024 ))
    factor="$(requirement disk_factor 5.6)"
    res_factor="$(requirement disk_results_factor 1.0)"
    hard_gb="$(requirement disk_hard_min_gb 20)"
    # With the work folder in Docker's volume, this drive holds only the results; the work files go to Docker's disk.
    if (( WORK_VOLUME )); then
        need_b="$(awk -v t="${total}" -v f="${res_factor}" 'BEGIN { printf "%.0f", t * f }')"
        work_b="$(awk -v t="${total}" -v f="${factor}" -v r="${res_factor}" 'BEGIN { printf "%.0f", t * (f - r) }')"
    else
        need_b="$(awk -v t="${total}" -v f="${factor}" 'BEGIN { printf "%.0f", t * f }')"
        work_b=0
    fi
    if (( count == 0 )); then
        warn "none of the input files named in the samplesheet were found. Check the fastq_1 / fastq_2 paths in $(basename "${ss}")."
    fi
    if (( WORK_VOLUME )); then
        say "  Disk: $(fmt_bytes "${free_b}") free on ${vol:-project volume}. Input files: $(fmt_bytes "${total}"); the analysis may need about $(fmt_bytes "${need_b}") on this drive (results) and $(fmt_bytes "${work_b}") in Docker's disk (work files)."
    else
        say "  Disk: $(fmt_bytes "${free_b}") free on ${vol:-project volume}. Input files: $(fmt_bytes "${total}"); the analysis may need about $(fmt_bytes "${need_b}")."
    fi
    if awk -v f="${free_b}" -v h="${hard_gb}" 'BEGIN { exit !(f / 1073741824 < h) }'; then
        fail "only $(fmt_bytes "${free_b}") is free on ${vol}. The analysis needs at least ${hard_gb} GB free to start. Free up space (or move the project to a drive with more room), then run this again."
    fi
    if (( free_b < need_b )); then
        if (( WORK_VOLUME )); then
            warn "${vol} has $(fmt_bytes "${free_b}") free, but the results may need about $(fmt_bytes "${need_b}") on this drive. The work files are in the Docker work volume (see Docker's disk below), so this drive mainly needs room for the results. If the disk fills up the run stops and can be resumed after you make room."
        else
            warn "${vol} has $(fmt_bytes "${free_b}") free, but the analysis may need about $(fmt_bytes "${need_b}"). If the disk fills up the run stops and can be resumed after you make room. To free space: run_atacseq.sh clean deletes metadata/work (the working copies; the next run rebuilds them). keep_work: false in config.yaml deletes them by itself after a successful run."
        fi
    fi
    if (( WORK_VOLUME )); then check_docker_disk "${work_b}"; fi
}

# Free space on Docker's disk, where the work volume lives: the filesystem that holds Docker Desktop's
# data folder (its settings name it), else the home folder's filesystem (on a Mac, Docker's disk is the
# Mac's own disk). Not the df inside the VM, which reports the VM's virtual disk. Silent when neither
# can be read. Stops when the work files do not fit; warns when they fit with less than 20% to spare.
check_docker_disk() {   # need_b
    local free_b
    docker info >/dev/null 2>&1 || return 0
    docker image inspect "${IMAGE}" >/dev/null 2>&1 || return 0
    free_b="$(docker_data_free_bytes)"
    [[ -n "${free_b}" ]] || free_b="$(df_avail_bytes "$(df -Pk "${HOME}" 2>/dev/null)")"
    [[ -n "${free_b}" ]] || return 0
    say "  Docker's disk: $(fmt_bytes "${free_b}") free (the work volume lives there)."
    if awk -v f="${free_b}" -v n="$1" 'BEGIN { exit !(f < n) }'; then
        fail "Docker's disk has $(fmt_bytes "${free_b}") free, but the work files need about $(fmt_bytes "$1"). Set work_volume: false in ${CONFIG_NAME} (the work folder then goes into the project folder, on a drive with more room), or free space on that disk, then run this again. run_atacseq.sh clean deletes the work volume."
    fi
    if awk -v f="${free_b}" -v n="$1" 'BEGIN { exit !(f < n * 1.2) }'; then
        warn "Docker's disk has $(fmt_bytes "${free_b}") free, and the work files may need about $(fmt_bytes "$1"), so there is little room to spare. If the disk fills up the run stops and can be resumed after you make room. To free space: run_atacseq.sh clean deletes the work volume (keep_work: false in config.yaml deletes it after a successful run)."
    fi
    return 0
}

# The work folder: the Docker volume (work_volume: true) or the project folder. Stops on a
# publish_mode that cannot be used with the volume. Makes no Docker calls.
check_work_volume() {
    local mode="" prob
    set_work_volume_vars
    [[ -f "${CONFIG_FILE}" ]] && mode="$(yaml_get "${CONFIG_FILE}" publish_mode)"
    prob="$(publish_mode_problem "${mode}" "${WORK_VOLUME}")"
    [[ -z "${prob}" ]] || fail "${prob}"
    if (( WORK_VOLUME )); then
        say "  Work folder: Docker volume ${WORK_VOLUME_NAME} (on Docker's disk, fast; resume survives restarts)"
    else
        say "  Work folder: metadata/work in the project folder (work_volume: false)"
    fi
}

# Everything the run reads from the project, checked before anything is started.
check_project_inputs() {
    local i
    info "Checking the input files named in ${CONFIG_FILE} and samplesheet.csv"
    check_output_dirs
    scan_inputs
    check_mounts_shared
    check_samplesheet
    check_config_files
    check_work_volume
    check_disk
    for (( i = 0; i < ${#MOUNT_ROOTS[@]}; i++ )); do
        say "  ${MOUNT_ROOTS[$i]} is mounted read-only inside the analysis (same path); nothing is copied from it."
    done
    if (( ${#MOUNT_ROOTS[@]} == 0 )); then say "  All inputs are inside the project folder."; fi
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

# ---------------------------------------------------------------------------
# Work volume: Nextflow's work folder lives in a per-project Docker volume (Docker's own
# disk; bind mounts of the project folder are about ten times slower). config.yaml:
# work_volume (default true; false/no/off/0 turn it off).
# ---------------------------------------------------------------------------
# The volume's name for a project folder (deterministic: the same name every run).
work_volume_name() {   # dir
    printf '%s_work\n' "$(container_name "$1")"
}

# Prints nothing when publish_mode can be used with the work volume, else the message.
# Links cannot cross from the volume (Docker's disk) to the project folder, and results/
# would point into the volume.
publish_mode_problem() {   # mode work_volume(1/0)
    if [[ "$2" == 1 ]]; then
        case "$1" in
            symlink|rellink|link)
                printf '%s\n' "publish_mode: $1 cannot be used with work_volume: true. Hard links and links cannot cross from the work volume (Docker's disk) to the project folder, and results/ would point into the volume. Set publish_mode: auto or copy, or work_volume: false." ;;
        esac
    fi
}

# Free bytes from `df -Pk` text: the Available column (4th) of the last line, in KB x 1024.
# Prints nothing when it cannot be read.
df_avail_bytes() {   # df_text
    printf '%s\n' "$1" | awk '
        NF > 0 { last = $0 }
        END {
            n = split(last, f, " ")
            if (n >= 4 && f[4] ~ /^[0-9]+$/) printf "%.0f\n", f[4] * 1024
        }'
}

# Docker Desktop's data folder from its settings JSON text ("DataFolder" or "dataFolder").
# Prints nothing when it is absent or the text is not such a file.
docker_data_folder() {   # json_text
    printf '%s\n' "$1" | sed -n -E 's/.*"(DataFolder|dataFolder)"[[:space:]]*:[[:space:]]*"([^"]*)".*/\2/p' | head -1
}

# Free bytes of the filesystem holding Docker Desktop's data folder (named in its settings
# file). Prints nothing when no settings file names a folder.
docker_data_free_bytes() {
    local f folder
    for f in "${HOME}/Library/Group Containers/group.com.docker/settings-store.json" \
             "${HOME}/Library/Group Containers/group.com.docker/settings.json"; do
        [[ -f "${f}" ]] || continue
        folder="$(docker_data_folder "$(cat "${f}" 2>/dev/null)")"
        [[ -n "${folder}" ]] || continue
        df_avail_bytes "$(df -Pk "$(expand_tilde "${folder}")" 2>/dev/null)"
        return 0
    done
}

# True (status 0) for the config words that mean off.
is_false_word() {   # value
    case "$(printf '%s' "$1" | tr '[:upper:]' '[:lower:]')" in
        false|no|off|0) return 0 ;;
    esac
    return 1
}

# WORK_VOLUME (1/0) and WORK_VOLUME_NAME for this project, from config.yaml (default: on).
set_work_volume_vars() {
    local v=""
    [[ -f "${CONFIG_FILE}" ]] && v="$(yaml_get "${CONFIG_FILE}" work_volume)"
    if is_false_word "${v}"; then WORK_VOLUME=0; else WORK_VOLUME=1; fi
    WORK_VOLUME_NAME="$(work_volume_name "${PROJECT_DIR}")"
}

# KB used by a work volume: du in a tiny container. Prints nothing when that cannot run
# (the pipeline image must be present locally; Docker must be up).
work_volume_kb() {   # volume
    docker image inspect "${IMAGE}" >/dev/null 2>&1 || return 0
    docker run --rm --memory 256m --mount "type=volume,source=$1,target=/w,readonly" "${IMAGE}" du -sk /w 2>/dev/null \
        | awk '$1 ~ /^[0-9]+$/ { print $1; exit }'
}

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

# The project must not be inside this pipeline's own repository folder.
check_not_inside_pipeline() {   # dir
    local repo_root
    repo_root="$(cd -P "${PIPELINE_DIR}/.." && pwd)"
    case "$1/" in
        "${repo_root}/"*) fail "Don't put a project inside the pipeline's own folder (${repo_root}). Choose another folder, e.g. ~/Documents/my_atacseq" ;;
    esac
}

# ---------------------------------------------------------------------------
# Resources: one budget per run, enforced as the container's hard limits.
#
#   memory = 75% of the Mac's RAM, but always leaving 3 GB to macOS
#            (a plain percentage starves small laptops: 80% of 8 GB leaves 1.6 GB)
#   CPUs   = all logical cores but 1 (but 2 from 8 cores up), so the desktop
#            stays responsive; heavy CPU use only slows a Mac, it can't harm it
#   ... and never more than Docker Desktop itself may use (Settings > Resources).
#   The pipeline's own cap (docker/requirements.conf) applies on top.
#   config.yaml max_cpus / max_memory_gb replace the default (up to Docker's limit).
#
# docker run --memory makes the budget a hard ceiling (a step that outgrows it
# is killed inside the container and retried, the Mac never starts swapping);
# --cpuset-cpus makes tools that use "all cores" (cutadapt --cores 0) see
# only the budget's cores.
# ---------------------------------------------------------------------------
requirement() {   # key default: the value as written (integer or decimal)
    local v
    v="$(grep -E "^[[:space:]]*$1[[:space:]]*=" "${PIPELINE_DIR}/docker/requirements.conf" 2>/dev/null | head -1 | sed -E 's/.*=[[:space:]]*([0-9]+(\.[0-9]+)?).*/\1/')"
    [[ "${v}" =~ ${NUM_RE} ]] && echo "${v}" || echo "$2"
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

# DOCKER_MEM_MB DOCKER_CPUS (what Docker Desktop's VM has; 0 when Docker is not running)
docker_resources() {
    local n m
    n="$(docker info --format '{{.NCPU}}' 2>/dev/null)"
    m="$(docker info --format '{{.MemTotal}}' 2>/dev/null)"
    [[ "${n}" =~ ^[0-9]+$ ]] && DOCKER_CPUS=${n} || DOCKER_CPUS=0
    [[ "${m}" =~ ^[0-9]+$ ]] && DOCKER_MEM_MB=$(( m / 1048576 )) || DOCKER_MEM_MB=0
}

# The default memory budget from this Mac alone (MB): 75%, but 3 GB left to macOS.
policy_mem_mb() {   # host_mb
    local a=$(( $1 * 3 / 4 )) b=$(( $1 - 3072 )) m
    m=$(( a < b ? a : b )); (( m < 0 )) && m=0
    echo "${m}"
}

# POLICY_CPUS: the default CPU budget from this Mac alone.
policy_cpus() {
    local spare=1
    (( HOST_CPUS >= 8 )) && spare=2
    POLICY_CPUS=$(( HOST_CPUS - spare )); (( POLICY_CPUS < 1 )) && POLICY_CPUS=1
}

# Docker Desktop's memory setting for analyses, by this Mac's memory (whole GB).
# The same tiers as the Windows launcher (run_atacseq.ps1, Get-WslTierGB).
mac_tier_gb() {   # total memory, whole GB
    if (( $1 < 12 )); then echo 6; elif (( $1 < 24 )); then echo 8; else echo 16; fi
}

# The plain-language note about Docker Desktop's memory. $1 is the sentence that ends it.
docker_note() {
    printf 'Docker Desktop is set to use %s GB. For this analysis it works best with about %s GB: open Docker Desktop > Settings > Resources, set Memory to %s GB, click Apply & restart. %s' \
        "${PLAN_DOCKER_GB}" "${PLAN_TIER_GB}" "${PLAN_TIER_GB}" "$1"
}

# The memory decision for one run. total_mb: this Mac's memory; docker_mb: Docker
# Desktop's memory (0 = Docker not running); want: config.yaml max_memory_gb ("" = not set);
# useful_mb: the pipeline's cap (0 = none). Changes nothing. Sets:
#   PLAN_POLICY_MB      the default budget (75% rule)
#   PLAN_BUDGET_MB      what a run gets, PLAN_REASON says why
#   PLAN_DOCKER_LIMITED 1 when Docker Desktop's memory holds the budget back
#   PLAN_TIER_GB        the setting for analyses (by this Mac's memory)
#   PLAN_DOCKER_GB      Docker Desktop's memory, whole GB
#   PLAN_NOTE           1 when Docker is below the setting and holds the budget back
#   PLAN_OVER           1 when max_memory_gb is more than Docker Desktop gives
#   PLAN_MIN_MB, PLAN_FAIL  the pipeline's minimum, and 1 when the budget is below it
memory_plan() {   # total_mb docker_mb want_gb useful_mb
    local total_mb=$1 docker_mb=$2 want=$3 useful_mb=$4 dock
    PLAN_POLICY_MB=$(policy_mem_mb "${total_mb}")
    dock=$(( docker_mb - 512 ))                      # Docker's VM less its own kernel and services
    (( dock > 0 )) || dock=${PLAN_POLICY_MB}
    PLAN_OVER=0; PLAN_DOCKER_LIMITED=0
    if [[ -n "${want}" ]]; then
        PLAN_REASON="max_memory_gb: ${want} in ${CONFIG_NAME:-config.yaml}"
        PLAN_BUDGET_MB=$(awk -v g="${want}" 'BEGIN { printf "%d", g * 1024 + 0.5 }')
        if (( PLAN_BUDGET_MB > dock )); then PLAN_OVER=1; PLAN_DOCKER_LIMITED=1; PLAN_BUDGET_MB=${dock}; fi
    else
        PLAN_BUDGET_MB=${PLAN_POLICY_MB}
        PLAN_REASON="75% of this Mac's $(gb "${total_mb}"), leaving at least 3 GB to the system"
        if (( useful_mb > 0 && useful_mb < PLAN_BUDGET_MB )); then
            PLAN_BUDGET_MB=${useful_mb}
            PLAN_REASON="the most this pipeline can use (docker/requirements.conf)"
        fi
        if (( dock < PLAN_BUDGET_MB )); then
            PLAN_BUDGET_MB=${dock}; PLAN_DOCKER_LIMITED=1
            PLAN_REASON="all that Docker Desktop may use: $(gb "${docker_mb}") of this Mac's $(gb "${total_mb}")"
        fi
    fi
    PLAN_MIN_MB=$(requirement_mb min_memory_gb 2.5)
    PLAN_TIER_GB=$(mac_tier_gb "$(round_gb "${total_mb}")")
    PLAN_DOCKER_GB=$(round_gb "${docker_mb}")
    PLAN_NOTE=0
    if (( PLAN_DOCKER_LIMITED && docker_mb > 0 && PLAN_DOCKER_GB < PLAN_TIER_GB )); then PLAN_NOTE=1; fi
    PLAN_FAIL=0
    if (( PLAN_BUDGET_MB < PLAN_MIN_MB )); then PLAN_FAIL=1; fi
}

# BUDGET_CPUS BUDGET_MEM_MB BUDGET_REASON for this run, or stops with an explanation.
run_budget() {
    host_resources; docker_resources; policy_cpus
    local want_mem="" want_cpu="" dock_cpu min_cpu
    if [[ -f "${CONFIG_FILE}" ]]; then
        want_mem="$(yaml_get "${CONFIG_FILE}" max_memory_gb)"
        want_cpu="$(yaml_get "${CONFIG_FILE}" max_cpus)"
    fi
    [[ "${want_mem}" =~ ${NUM_RE} ]] || want_mem=""
    memory_plan "${HOST_MEM_MB}" "${DOCKER_MEM_MB}" "${want_mem}" "$(requirement_mb max_useful_memory_gb 0)"
    BUDGET_MEM_MB=${PLAN_BUDGET_MB}; BUDGET_REASON="${PLAN_REASON}"
    if (( PLAN_OVER )); then
        warn "max_memory_gb: ${want_mem} is more than Docker Desktop may use; using $(gb "${PLAN_BUDGET_MB}")."
    fi

    dock_cpu=${DOCKER_CPUS}; (( dock_cpu > 0 )) || dock_cpu=${POLICY_CPUS}
    if [[ "${want_cpu}" =~ ^[0-9]+$ ]] && (( want_cpu >= 1 )); then
        BUDGET_CPUS=$(( want_cpu < dock_cpu ? want_cpu : dock_cpu ))
        BUDGET_REASON="${BUDGET_REASON}; max_cpus: ${want_cpu}"
    else
        BUDGET_CPUS=$(( POLICY_CPUS < dock_cpu ? POLICY_CPUS : dock_cpu ))
    fi

    if (( PLAN_FAIL )); then
        if [[ -n "${want_mem}" ]]; then
            fail "max_memory_gb: ${want_mem} in ${CONFIG_NAME} is below the $(gb "${PLAN_MIN_MB}") this pipeline needs. Raise it or remove the line."
        elif (( PLAN_DOCKER_LIMITED )); then
            fail "$(docker_note "The analysis needs at least $(gb "${PLAN_MIN_MB}"), so it cannot start with Docker's current setting. Make the change, then run this again.")"
        else
            fail "This Mac has $(gb "${HOST_MEM_MB}") of memory. Keeping 3 GB for the system leaves less than the $(gb "${PLAN_MIN_MB}")
this pipeline needs. Run it on a Mac with at least $(( (PLAN_MIN_MB + 3072 + 1023) / 1024 )) GB, or set max_memory_gb in config.yaml
to use more (the Mac may then become slow while it runs)."
        fi
    fi
    if (( PLAN_NOTE )); then say "$(docker_note "(Run continues now with ${PLAN_DOCKER_GB} GB.)")"; fi

    min_cpu=$(requirement min_cpus 1)
    if (( BUDGET_CPUS < min_cpu )); then warn "only ${BUDGET_CPUS} CPU core(s) for the analysis: it will work, but slowly."; fi
    if (( HOST_FREE_MB > 0 && HOST_FREE_MB < BUDGET_MEM_MB )); then
        warn "only $(gb "${HOST_FREE_MB}") of memory is free right now (other programs use the rest), and the analysis may use up to $(gb "${BUDGET_MEM_MB}"). Closing large programs (e.g. a browser with many tabs) keeps the Mac responsive."
    fi
    return 0
}

# setup: report what each run will get, and how to give Docker more if it is the limit.
report_resources() {
    host_resources; docker_resources; policy_cpus
    local cpus=${DOCKER_CPUS}
    (( cpus > 0 )) || cpus=${POLICY_CPUS}
    (( POLICY_CPUS < cpus )) && cpus=${POLICY_CPUS}
    memory_plan "${HOST_MEM_MB}" "${DOCKER_MEM_MB}" "" "$(requirement_mb max_useful_memory_gb 0)"
    echo
    say "This Mac: ${HOST_CPUS} CPU cores, $(gb "${HOST_MEM_MB}") memory. Docker Desktop may use: ${DOCKER_CPUS} cores, $(gb "${DOCKER_MEM_MB}")."
    say "Each analysis will use at most ${cpus} CPU cores and $(gb "${PLAN_BUDGET_MB}") of memory (max_cpus / max_memory_gb in a project's config.yaml change this)."
    if (( PLAN_NOTE )); then say "$(docker_note "(Change it before the first analysis for the best result.)")"; fi
    if (( PLAN_FAIL )); then
        warn "that is less than the $(gb "${PLAN_MIN_MB}") this pipeline needs: runs will refuse to start."
    fi
}

# Test settings from this Mac's environment, passed on when set (see docker/run_in_container.sh).
add_optional_env() {
    EXTRA_ENV=()
    if [[ -n "${ATAC_NF_PROFILE:-}" ]]; then EXTRA_ENV+=(-e "ATAC_NF_PROFILE=${ATAC_NF_PROFILE}"); fi
    if [[ -n "${ATAC_NF_EXTRA_ARGS:-}" ]]; then EXTRA_ENV+=(-e "ATAC_NF_EXTRA_ARGS=${ATAC_NF_EXTRA_ARGS}"); fi
    if [[ -n "${ATAC_WORK_DIR:-}" ]]; then EXTRA_ENV+=(-e "ATAC_WORK_DIR=${ATAC_WORK_DIR}"); fi
    return 0
}

# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------
cmd_help() {
    cat <<EOF
ATAC-seq pipeline (Docker version for Mac)

  run_atacseq.sh setup            one-time: Docker Desktop + the pipeline's software image
  run_atacseq.sh init [FOLDER]    set up a new project folder (default: current folder)
  run_atacseq.sh                  run the pipeline for the project in the current folder
  run_atacseq.sh run CONFIG.yaml  ... with another config file
  run_atacseq.sh status           progress of the current/last run (and why it failed)
  run_atacseq.sh logs             follow the live log
  run_atacseq.sh stop             stop the running analysis (re-running resumes it)
  run_atacseq.sh clean            delete the work folder, resume cache and work volume (metadata/work, metadata/.nextflow, Docker volume) now
  run_atacseq.sh script utils/X.R [ARGS]  run one of the project's own scripts (R, Python or shell) in the pipeline's software environment

Options: --detach  --yes  --pause
Pipeline folder: ${PIPELINE_DIR}
EOF
    done_ 0
}

cmd_setup() {
    info "ATAC-seq setup (Docker version for Mac)"
    ensure_docker 1
    repair_line_endings
    ensure_rosetta
    pull_image

    report_resources

    echo
    docker run --rm --platform linux/amd64 "${IMAGE}" bash -c 'echo nextflow: $(nextflow -version 2>&1 | grep -i version | head -1); echo R: $(R --version 2>&1 | head -1); echo python3: $(python3 --version 2>&1); echo samtools: $(samtools --version 2>&1 | head -1); echo bowtie2: $(bowtie2 --version 2>&1 | head -1)'

    # Put run_atacseq.sh on PATH for new terminals (zsh is macOS's default shell).
    local rc line="export PATH=\"${PIPELINE_DIR}:\$PATH\"   # ATAC-seq pipeline (run_atacseq.sh)"
    for rc in "${HOME}/.zshrc" "${HOME}/.bash_profile"; do
        if [[ "${rc}" == *.zshrc || -f "${rc}" ]]; then
            if ! grep -qF "${PIPELINE_DIR}:" "${rc}" 2>/dev/null; then
                printf '\n%s\n' "${line}" >> "${rc}"
                say "Added ${PIPELINE_DIR} to your PATH in ${rc}"
            fi
        fi
    done
    chmod +x "${PIPELINE_DIR}/run_atacseq.sh" 2>/dev/null || true
    echo
    info "Setup complete. Next: create a project folder with   run_atacseq.sh init <folder>"
    say  "(in a NEW terminal window, run_atacseq.sh works from any folder)"
    done_ 0
}

cmd_init() {
    local target="${arg1:-$(pwd)}" f csv_n=0 src_agent="${PIPELINE_DIR}/agent" c_src rel c_n=0
    mkdir -p "${target}" || fail "cannot create ${target}"
    target="$(abspath "${target}")"
    check_not_inside_pipeline "${target}"
    info "Setting up a new ${PIPELINE_NAME} project in ${target} ..."
    mkdir -p "${target}/original_data/fastq" "${target}/metadata" "${target}/results" "${target}/utils" "${target}/examples"
    if [[ -f "${target}/config.yaml" ]]; then say "  config.yaml already exists — leaving it alone."
    else cp "${PIPELINE_DIR}/examples/config.yaml" "${target}/config.yaml"; say "  Fetched config.yaml"; fi

    # CSV templates: copied into examples/ (never overwritten). Copy the ones you need to the project folder and fill them in.
    for f in "${PIPELINE_DIR}"/examples/*.csv; do
        [[ -f "${f}" ]] || continue
        if [[ ! -e "${target}/examples/$(basename "${f}")" ]]; then
            cp "${f}" "${target}/examples/"; csv_n=$(( csv_n + 1 ))
        fi
    done
    say "  Fetched examples/ (${csv_n} CSV templates, e.g. samplesheet.csv)"

    # AI agent instructions and agents, when this pipeline folder has them (agent/ holds them).
    for f in AGENTS.md CLAUDE.md; do
        [[ -f "${src_agent}/${f}" ]] || continue
        if [[ -e "${target}/${f}" ]]; then say "  ${f} already exists — leaving it alone."; continue; fi
        awk -v p="${PIPELINE_DIR}" '{ gsub(/\{\{PIPELINE_DIR\}\}/, p); print }' "${src_agent}/${f}" > "${target}/${f}"
        say "  Fetched ${f}"
    done
    c_src="${src_agent}/.claude/agents"
    if [[ -d "${c_src}" ]]; then
        while IFS= read -r rel; do
            rel="${rel#./}"
            [[ -n "${rel}" && ! -e "${target}/.claude/agents/${rel}" ]] || continue
            mkdir -p "$(dirname "${target}/.claude/agents/${rel}")"
            cp "${c_src}/${rel}" "${target}/.claude/agents/${rel}"
            c_n=$(( c_n + 1 ))
        done <<< "$(cd "${c_src}" && find . -type f)"
        if (( c_n > 0 )); then say "  Fetched .claude/agents/ (${c_n} files for AI agents; existing files left alone)"; fi
    fi

    # Launchers inside the project (rewritten each time: they hold the pipeline's location).
    # The terminal launcher lives in utils/ (the project root keeps only the double-click RUN_ANALYSIS.command).
    local q; q="$(printf '%q' "${PIPELINE_DIR}/run_atacseq.sh")"
    cat > "${target}/utils/run_atacseq.sh" <<EOF
#!/bin/bash
# Runs the ATAC-seq pipeline (Docker version) for THIS project folder.
# Same commands as the pipeline's own launcher: (no argument) = run, status, logs, stop, clean, script utils/X.R [ARGS].
# Generated by 'run_atacseq.sh init'; the pipeline lives in:
#   ${PIPELINE_DIR}
cd "\$(dirname "\$0")/.." && exec bash ${q} "\$@"
EOF
    cat > "${target}/RUN_ANALYSIS.command" <<EOF
#!/bin/bash
# Double-click to run (or resume) the ATAC-seq analysis of this folder.
cd "\$(dirname "\$0")" && exec bash ${q} run --pause
EOF
    chmod +x "${target}/utils/run_atacseq.sh" "${target}/RUN_ANALYSIS.command"
    say "  Created utils/run_atacseq.sh (terminal) and RUN_ANALYSIS.command (double-click)"

    # .gitattributes: LF line endings for the scripts, so a utils/ copy stays byte-identical
    # to the pipeline's default on every OS. Never overwritten.
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
    elif have git; then
        [[ -f "${target}/.gitignore" ]] || cp "${PIPELINE_DIR}/examples/project.gitignore" "${target}/.gitignore"
        local paths="" p
        for p in config.yaml examples .gitignore .gitattributes AGENTS.md CLAUDE.md RUN_ANALYSIS.command utils .claude; do
            if [[ -e "${target}/${p}" ]]; then paths="${paths} ${p}"; fi
        done
        ( cd "${target}" && git init -q . \
            && git add ${paths} 2>/dev/null \
            && git -c user.name=atacseq -c user.email=atacseq@localhost commit -q -m "Initialize ${PIPELINE_NAME} project" ) >/dev/null 2>&1 || true
        say "  Initialized git repository (raw data / results / metadata excluded - see .gitignore)"
    else
        say "  (git not available — skipped version control; optional)"
    fi

    echo
    say "Next steps:"
    say "  1. Copy the input tables you need from examples/ into this folder and fill them in:"
    say "       samplesheet.csv (one row per sample) and, if you merge peaks, peak_merges.csv (see config.yaml)."
    say "  2. Edit config.yaml: project_name, samplesheet, and the reference files (genome_index, contamination_index, blacklist, chrom_sizes, gtf, tss_bed)."
    say "     Write Mac paths in full, e.g. /Users/you/data/x.fastq.gz, or relative to this folder, e.g. original_data/x.fastq.gz."
    say "     Docker Desktop shares /Users, /Volumes, /private and /tmp by default; another folder must be added under"
    say "     Settings > Resources > File sharing."
    say "     Network paths (smb://..., //server/share/...) are not supported: mount the share in Finder and use its /Volumes/<name>/... path."
    say "     Your input files stay where they are: they are never moved or copied."
    say "  3. Run:  double-click RUN_ANALYSIS.command   (or type  ./utils/run_atacseq.sh  in this folder)"
    done_ 0
}

# keep_work: false with the work volume: after a successful run the container has emptied
# the volume, so it is removed here. Quiet when Docker is down (clean removes it later).
drop_work_volume() {
    local keep=""
    [[ -f "${CONFIG_FILE}" ]] || return 0
    set_work_volume_vars
    (( WORK_VOLUME )) || return 0
    keep="$(yaml_get "${CONFIG_FILE}" keep_work)"
    is_false_word "${keep}" || return 0
    docker info >/dev/null 2>&1 || return 0
    docker volume inspect "${WORK_VOLUME_NAME}" >/dev/null 2>&1 || return 0
    docker volume rm "${WORK_VOLUME_NAME}" >/dev/null 2>&1 && say "  Work volume removed (keep_work: false)."
    return 0
}

# Follows the live log of a running container, then reports its outcome.
follow_run() {
    local name="$1" st="" statusf="${PROJECT_DIR}/metadata/nf_logs/last_run.status"
    trap 'echo; say "Stopped watching; the analysis keeps running (run_atacseq.sh status)."; done_ 0' INT
    docker logs -f "${name}" 2>&1
    trap - INT
    sleep 1
    [[ -f "${statusf}" ]] && st="$(cat "${statusf}")"
    echo
    if [[ "${st}" == succeeded ]]; then
        info "Pipeline finished successfully."
        say  "  Results : ${PROJECT_DIR}/results/"
        drop_work_volume
        done_ 0
    elif [[ "${st}" == running ]] && is_running "${name}"; then
        say "Stopped watching; the analysis is still running."
        done_ 0
    else
        printf '%sPipeline did not finish: %s%s\n' "${C_ERR}" "${st}" "${C_OFF}"
        say "  See:  run_atacseq.sh status   (it names the failed step and where its messages were saved)"
        say "  Fix the problem, then run again: finished steps are reused, not redone."
        done_ 1
    fi
}

# An earlier run's work folder (metadata/work) holds the resume state, which the work volume
# does not have: with work_volume: true the run starts from scratch, so ask first.
check_old_work_folder() {
    local old="${PROJECT_DIR}/metadata/work"
    (( WORK_VOLUME )) || return 0
    [[ -d "${old}" && -n "$(ls -A "${old}" 2>/dev/null)" ]] || return 0
    if docker volume inspect "${WORK_VOLUME_NAME}" >/dev/null 2>&1; then
        warn_once "metadata/work (an earlier work folder) is not used: the work volume holds the work now. run_atacseq.sh clean deletes it."
        return 0
    fi
    warn "an earlier work folder (metadata/work) was found. With work_volume: true the run starts from scratch, because the resume state is in that folder."
    say  "  To resume that run instead, set work_volume: false in ${CONFIG_NAME} and run again."
    say  "  To continue with the work volume, answer yes. The old folder stays until you run clean, which deletes it."
    confirm "Start from scratch with the work volume?" \
        || fail "nothing was started. Set work_volume: false in ${CONFIG_NAME} to resume the earlier run, or answer yes to start from scratch."
}

cmd_run() {
    project_from_args
    [[ -f "${CONFIG_FILE}" ]] || fail "${CONFIG_NAME} not found in ${PROJECT_DIR}. Set the project up first:  run_atacseq.sh init"
    check_not_inside_pipeline "${PROJECT_DIR}"
    local name ts image_id work_args
    name="$(container_name "${PROJECT_DIR}")"
    ensure_docker 0
    is_running "${name}" && fail "this project's analysis is already running. Watch it with:  run_atacseq.sh logs   (or stop it: run_atacseq.sh stop)"
    set_work_volume_vars
    check_old_work_folder
    repair_line_endings
    check_project_inputs
    run_budget
    ensure_image
    build_mounts
    add_optional_env

    ts="$(date +%Y%m%d_%H%M%S)"
    image_id="$(docker image inspect "${IMAGE}" --format '{{.Id}}' 2>/dev/null)"

    # The work volume (created when new; creating an existing one is harmless), or the old setting.
    work_args=()
    if (( WORK_VOLUME )); then
        docker volume create --label "atacseq.project=${PROJECT_DIR}" "${WORK_VOLUME_NAME}" >/dev/null \
            || fail "Docker could not create the work volume ${WORK_VOLUME_NAME} (see the message above)."
        work_args=(--mount "type=volume,source=${WORK_VOLUME_NAME},target=/nfwork" -e "ATAC_WORK_VOLUME=1" -e "ATAC_WORK_VOLUME_NAME=${WORK_VOLUME_NAME}")
    else
        work_args=(-e "ATAC_WORK_VOLUME=0")
    fi

    # Mark the run as started now: until the container rewrites it, `status`
    # (and an agent polling it) would otherwise still see the previous run's outcome.
    mkdir -p "${PROJECT_DIR}/metadata/nf_logs"
    echo "running" > "${PROJECT_DIR}/metadata/nf_logs/last_run.status"

    docker run -d --rm --init --name "${name}" --platform linux/amd64 \
        --label "atacseq.project=${PROJECT_DIR}" \
        --memory "${BUDGET_MEM_MB}m" --memory-swap "${BUDGET_MEM_MB}m" --cpuset-cpus "0-$(( BUDGET_CPUS - 1 ))" \
        -e "ATAC_BUDGET_CPUS=${BUDGET_CPUS}" -e "ATAC_BUDGET_MEMORY_MB=${BUDGET_MEM_MB}" -e "ATAC_BUDGET_REASON=${BUDGET_REASON}" \
        --mount "type=bind,source=${PROJECT_DIR},target=/project" \
        --mount "type=bind,source=${PIPELINE_DIR},target=/pipeline,readonly" \
        ${MOUNT_ARGS[@]+"${MOUNT_ARGS[@]}"} \
        ${work_args[@]+"${work_args[@]}"} \
        -e "ATAC_HDD=0" \
        -e "ATAC_MOUNTS_HOST=${MOUNTS_TEXT}" \
        -e "ATAC_IMAGE=${IMAGE}" \
        -e "ATAC_IMAGE_ID=${image_id}" \
        -e "ATAC_PROJECT_HOST_DIR=${PROJECT_DIR}" \
        ${EXTRA_ENV[@]+"${EXTRA_ENV[@]}"} \
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
    say  "  Uses at most ${BUDGET_CPUS} CPU cores and $(gb "${BUDGET_MEM_MB}") memory (${BUDGET_REASON})."
    say  "  Progress : run_atacseq.sh status      Stop: run_atacseq.sh stop"
    say  "  Closing this window does NOT stop the analysis. Keep the Mac plugged in and the lid"
    say  "  open until it finishes (it is kept from idle-sleeping while the analysis runs)."
    (( DETACH )) && done_ 0
    say  "  (Press Ctrl+C to stop watching; the analysis keeps running.)"
    echo
    follow_run "${name}"
}

latest_log() { ls -1 "${PROJECT_DIR}"/metadata/nf_logs/nextflow_*.log 2>/dev/null | sort | tail -1; }

# Where the work is kept and how big the work volume is. Never fails: Docker may be down.
status_work_volume() {
    local kb
    set_work_volume_vars
    if (( ! WORK_VOLUME )); then say "Work folder: metadata/work (work_volume: false)"; return 0; fi
    if ! have docker || ! docker info >/dev/null 2>&1; then say "Work volume: not checked (Docker is not running)"; return 0; fi
    if ! docker volume inspect "${WORK_VOLUME_NAME}" >/dev/null 2>&1; then say "Work volume: none yet"; return 0; fi
    if ! docker image inspect "${IMAGE}" >/dev/null 2>&1; then
        say "Work volume: ${WORK_VOLUME_NAME} (size not checked: the pipeline image is not downloaded)"; return 0
    fi
    kb="$(work_volume_kb "${WORK_VOLUME_NAME}")"
    if [[ -n "${kb}" ]]; then say "Work volume: ${WORK_VOLUME_NAME}, size $(fmt_bytes $(( kb * 1024 )))"
    else say "Work volume: ${WORK_VOLUME_NAME} (size not checked)"; fi
    return 0
}

cmd_status() {
    project_from_args
    local name running=0 log st trace ft statusf
    name="$(container_name "${PROJECT_DIR}")"
    statusf="${PROJECT_DIR}/metadata/nf_logs/last_run.status"
    if have docker && docker info >/dev/null 2>&1 && is_running "${name}"; then running=1; fi
    log="$(latest_log)"
    [[ -n "${log}" ]] || { say "No run yet in ${PROJECT_DIR}."; done_ 0; }
    if (( running )); then info "Status: RUNNING"
    elif [[ -f "${statusf}" ]]; then
        st="$(cat "${statusf}")"
        [[ "${st}" == running ]] && st="interrupted (Docker or the Mac stopped mid-run; run again to resume)"
        info "Status: last run ${st}"
    fi
    say "Log: ${log}"
    status_work_volume
    # The trace is named by the launch timestamp (nextflow.config): the newest one is this project's latest run.
    trace="$(ls -1 "${PROJECT_DIR}"/metadata/pipeline_info/execution_trace*.txt 2>/dev/null | sort | tail -1)"
    if [[ -s "${trace}" ]]; then
        echo
        say "Steps (this run):"
        awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next }
            { n[$c["status"]]++; if ($c["status"] == "FAILED") f[++k] = $c["name"] " (exit " $c["exit"] ")" }
            END { for (s in n) printf "  %-10s %d\n", s, n[s]; for (i = 1; i <= k; i++) print "  FAILED: " f[i] }' "${trace}"
    fi
    ft="${PROJECT_DIR}/metadata/failed_tasks/$(basename "${log}" .log | sed 's/^nextflow_//')"
    if [[ -d "${ft}" ]]; then echo; say "Files of the failed step(s): ${ft}"; fi
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
    local name d sz kb=0 vol_kb=0 dirs=() docker_up=0 vol_exists=0
    name="$(container_name "${PROJECT_DIR}")"
    WORK_VOLUME_NAME="$(work_volume_name "${PROJECT_DIR}")"
    # metadata/ holds the work folder and the resume cache.
    for d in metadata/work metadata/.nextflow; do
        [[ -e "${PROJECT_DIR}/${d}" ]] && dirs+=("${d}")
    done
    # The work volume can only be looked at with Docker running; without it the volume is skipped.
    if have docker && docker info >/dev/null 2>&1; then docker_up=1; fi
    if (( docker_up )) && docker volume inspect "${WORK_VOLUME_NAME}" >/dev/null 2>&1; then vol_exists=1; fi
    if (( ${#dirs[@]} == 0 && vol_exists == 0 )); then
        if (( docker_up == 0 )) && have docker; then
            say "Nothing to clean in the project folder. The work volume was skipped: Docker is not running (run clean again with Docker running to check it)."
        else
            say "Nothing to clean (no metadata/work or metadata/.nextflow folder, and no work volume)."
        fi
        done_ 0
    fi
    is_running "${name}" && fail "the analysis is running: stop it first (run_atacseq.sh stop)."
    for d in ${dirs[@]+"${dirs[@]}"}; do
        sz="$(du -sk "${PROJECT_DIR}/${d}" 2>/dev/null | awk '{ print $1 }')"
        kb=$(( kb + ${sz:-0} ))
    done
    if (( vol_exists )); then
        vol_kb="$(work_volume_kb "${WORK_VOLUME_NAME}")"
        kb=$(( kb + ${vol_kb:-0} ))
    fi
    if (( docker_up == 0 )) && have docker; then
        say "Docker is not running: the work volume is not checked and is skipped."
    fi
    say "This deletes:"
    for d in ${dirs[@]+"${dirs[@]}"}; do say "  ${PROJECT_DIR}/${d}"; done
    if (( vol_exists )); then say "  Docker volume ${WORK_VOLUME_NAME} (the work folder when work_volume is true)"; fi
    say "They are Nextflow's working copies of every step (about $(gb $(( kb / 1024 )))) and its resume cache."
    say "results/, original_data/ and your config and CSV files are NOT touched, but the next run will start from scratch."
    confirm "Delete them?" || { say "Nothing deleted."; done_ 0; }
    if (( ${#dirs[@]} > 0 )) && docker image inspect "${IMAGE}" >/dev/null 2>&1; then
        # Deleted from inside Linux, as on Windows: the work folder holds Linux symlinks.
        docker run --rm --platform linux/amd64 --mount "type=bind,source=${PROJECT_DIR},target=/project" "${IMAGE}" \
            rm -rf /project/metadata/work /project/metadata/.nextflow >/dev/null 2>&1 || true
    fi
    for d in ${dirs[@]+"${dirs[@]}"}; do rm -rf "${PROJECT_DIR}/${d}"; done
    if (( vol_exists )); then
        docker volume rm "${WORK_VOLUME_NAME}" >/dev/null 2>&1 \
            || warn "could not delete the Docker volume ${WORK_VOLUME_NAME}: run clean again, or remove it in Docker Desktop (Volumes)."
    fi
    say "Deleted. About $(gb $(( kb / 1024 ))) of disk space is free again."
    done_ 0
}

# run_atacseq.sh script utils/X.R [ARGS]: one of the project's own scripts, run in the
# pipeline's software environment. Foreground only; the project is the current folder.
cmd_script() {
    local script_arg="${arg1}" src full rel base stem name ts logrel image_id rc a i
    local script_args=()
    if [[ -z "${script_arg}" ]]; then
        fail "usage: run_atacseq.sh script utils/<name>.<R|py|sh> [args...]   (run from the project folder)"
    fi
    PROJECT_DIR="$(pwd -P)"
    CONFIG_NAME="config.yaml"
    CONFIG_FILE="${PROJECT_DIR}/${CONFIG_NAME}"

    # The script: relative to the project folder (or absolute), inside it, with an extension the container runs.
    src="$(resolve_config_path "${script_arg}")"
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
    is_running "${name}" && warn "the project's analysis is running: this script shares the Mac's memory with it."
    repair_line_endings
    scan_inputs
    # Folders named in the script's arguments are mounted too.
    for (( i = 0; i < ${#script_args[@]}; i++ )); do
        a="${script_args[$i]}"
        case "${a}" in
            //*|*://*) refuse_network "A script argument" "${a}" ;;
            /*|"~/"*) add_root_for "$(expand_tilde "${a}")" ;;
        esac
    done
    check_mounts_shared
    run_budget
    ensure_image
    build_mounts

    ts="$(date +%Y%m%d_%H%M%S)"
    logrel="metadata/script_logs/${stem}_${ts}.log"
    mkdir -p "${PROJECT_DIR}/metadata/script_logs"
    image_id="$(docker image inspect "${IMAGE}" --format '{{.Id}}' 2>/dev/null)"

    info "Running ${rel} in the pipeline's software environment (at most ${BUDGET_CPUS} CPU cores and $(gb "${BUDGET_MEM_MB}") memory). Log: ${PROJECT_DIR}/${logrel}"
    docker run --rm --init --platform linux/amd64 \
        --label "atacseq.project=${PROJECT_DIR}" \
        --memory "${BUDGET_MEM_MB}m" --memory-swap "${BUDGET_MEM_MB}m" --cpuset-cpus "0-$(( BUDGET_CPUS - 1 ))" \
        -e "ATAC_BUDGET_CPUS=${BUDGET_CPUS}" -e "ATAC_BUDGET_MEMORY_MB=${BUDGET_MEM_MB}" -e "ATAC_BUDGET_REASON=${BUDGET_REASON}" \
        --mount "type=bind,source=${PROJECT_DIR},target=/project" \
        --mount "type=bind,source=${PIPELINE_DIR},target=/pipeline,readonly" \
        ${MOUNT_ARGS[@]+"${MOUNT_ARGS[@]}"} \
        -e "ATAC_MOUNTS_HOST=${MOUNTS_TEXT}" \
        -e "ATAC_IMAGE=${IMAGE}" \
        -e "ATAC_IMAGE_ID=${image_id}" \
        -e "ATAC_PROJECT_HOST_DIR=${PROJECT_DIR}" \
        -w /project \
        "${IMAGE}" bash /pipeline/docker/run_script_in_container.sh "${rel}" "${logrel}" \
        ${script_args[@]+"${script_args[@]}"}
    rc=$?

    echo
    say "Script finished (exit code ${rc}). Log: ${PROJECT_DIR}/${logrel}"
    done_ "${rc}"
}

# Hidden (maintainers, tests): the Docker memory decision for hypothetical inputs. Changes nothing.
cmd_mem_plan() {
    local i=0 a v total_mb="" docker_mb=""
    while (( i < ${#positional[@]} )); do
        a="${positional[$i]}"
        case "${a}" in
            --total-ram-gb|--docker-gb)
                v="${positional[$(( i + 1 ))]:-}"
                [[ "${v}" =~ ${NUM_RE} ]] || fail "${a} needs a number of GB, e.g. ${a} 16"
                if [[ "${a}" == --total-ram-gb ]]; then
                    total_mb=$(awk -v g="${v}" 'BEGIN { printf "%d", g * 1024 + 0.5 }')
                else
                    docker_mb=$(awk -v g="${v}" 'BEGIN { printf "%d", g * 1024 + 0.5 }')
                fi
                i=$(( i + 1 ))
                ;;
            *) fail "mem-plan: unknown argument ${a} (use --total-ram-gb N and --docker-gb N)" ;;
        esac
        i=$(( i + 1 ))
    done
    if [[ -z "${total_mb}" ]]; then
        host_resources; total_mb=${HOST_MEM_MB}
        say "(no --total-ram-gb given: this Mac's memory is used)"
    fi
    if [[ -z "${docker_mb}" ]]; then
        docker_resources; docker_mb=${DOCKER_MEM_MB}
        if (( docker_mb == 0 )); then
            docker_mb=$(( total_mb / 2 ))
            say "(Docker is not running: its default, half of the memory, is used)"
        fi
    fi
    memory_plan "${total_mb}" "${docker_mb}" "" "$(requirement_mb max_useful_memory_gb 0)"
    say "Plan for a Mac with $(gb "${total_mb}") of memory and Docker Desktop set to $(gb "${docker_mb}") (nothing is read from or changed in Docker):"
    say "  Setting for analyses : ${PLAN_TIER_GB} GB (by this Mac's memory)"
    say "  Docker Desktop now   : ${PLAN_DOCKER_GB} GB"
    say "  Run budget           : $(gb "${PLAN_BUDGET_MB}") (${PLAN_REASON})"
    if (( PLAN_NOTE && PLAN_FAIL )); then
        say "  Note                 : $(docker_note "The analysis needs at least $(gb "${PLAN_MIN_MB}"), so it cannot start with Docker's current setting.")"
    elif (( PLAN_NOTE )); then
        say "  Note                 : $(docker_note "(Run continues now with ${PLAN_DOCKER_GB} GB.)")"
    else
        say "  Note                 : none"
    fi
    if (( PLAN_FAIL )); then
        if (( PLAN_DOCKER_LIMITED )); then say "  Result               : the run would stop (Docker's memory is below what the analysis needs)"
        else say "  Result               : the run would stop (this Mac has less memory than the pipeline needs)"; fi
    else
        say "  Result               : the run continues with $(gb "${PLAN_BUDGET_MB}")"
    fi
    done_ 0
}

# Hidden (maintainers, tests): the folders a run would mount, from config.yaml and the
# project's CSV files. Starts nothing.
cmd_scan_paths() {
    project_from_args
    [[ -f "${CONFIG_FILE}" ]] || fail "${CONFIG_NAME} not found in ${PROJECT_DIR}."
    say "Scanning ${CONFIG_NAME} and the CSV files in ${PROJECT_DIR}"
    check_output_dirs
    scan_inputs
    local i r
    if (( ${#MOUNT_ROOTS[@]} == 0 )); then say "  No absolute input folder is named: only the project folder is visible inside the analysis."; fi
    for (( i = 0; i < ${#MOUNT_ROOTS[@]}; i++ )); do
        r="${MOUNT_ROOTS[$i]}"
        say "  would mount ${r} read-only at ${r}"
        if ! mac_docker_shared "${r}"; then
            warn "${r} is not shared by Docker Desktop by default: add it under Settings > Resources > File sharing, or the run stops."
        fi
    done
    say "  would mount ${PROJECT_DIR} read-write at /project"
    check_config_files
    done_ 0
}

# One PASS/FAIL line for selftest: name expected actual.
selftest_check() {
    if [[ "$2" == "$3" ]]; then say "PASS  $1"
    else say "FAIL  $1 (expected: '$2', got: '$3')"; ST_FAILS=$(( ST_FAILS + 1 )); fi
    ST_TOTAL=$(( ST_TOTAL + 1 ))
}

# Hidden (maintainers): PASS/FAIL checks of the pure helpers (work volume, df, Docker settings).
# Needs no Docker. Exit 1 when a check fails.
cmd_selftest() {
    local re='^[a-zA-Z0-9][a-zA-Z0-9_.-]*$' v ok m got
    ST_FAILS=0; ST_TOTAL=0

    v="$(work_volume_name "/Users/me/ATAC-seq project")"
    selftest_check "volume name is the same every time" "${v}" "$(work_volume_name "/Users/me/ATAC-seq project")"
    if [[ "${v}" =~ $re ]]; then ok=yes; else ok=no; fi
    selftest_check "volume name is a valid Docker name" "yes" "${ok}"
    case "${v}" in *_work) ok=yes ;; *) ok=no ;; esac
    selftest_check "volume name ends with _work" "yes" "${ok}"

    for m in symlink rellink link; do
        got="$(publish_mode_problem "${m}" 1)"
        if [[ -n "${got}" ]]; then ok=1; else ok=0; fi
        selftest_check "publish_mode ${m} with the work volume gives a message" "1" "${ok}"
    done
    for m in auto copy; do
        selftest_check "publish_mode ${m} with the work volume gives no message" "" "$(publish_mode_problem "${m}" 1)"
    done
    selftest_check "publish_mode symlink without the work volume gives no message" "" "$(publish_mode_problem symlink 0)"

    selftest_check "df: Available of the last line x 1024" "1005568000000" \
        "$(df_avail_bytes $'Filesystem 1024-blocks Used Available Capacity Mounted on\noverlay 1055762868 20000000 982000000 3% /')"
    selftest_check "df: garbage gives nothing" "" "$(df_avail_bytes "garbage text, not a df table")"
    selftest_check "df: empty text gives nothing" "" "$(df_avail_bytes "")"

    selftest_check "Docker settings: DataFolder" "/Volumes/Data/docker" \
        "$(docker_data_folder '{"DataFolder": "/Volumes/Data/docker", "cpus": 4}')"
    selftest_check "Docker settings: dataFolder" "/Users/me/dd" \
        "$(docker_data_folder '{"cpus":4, "dataFolder":"/Users/me/dd"}')"
    selftest_check "Docker settings: invalid JSON gives nothing" "" "$(docker_data_folder 'not json at all {')"
    selftest_check "Docker settings: no folder key gives nothing" "" "$(docker_data_folder '{"cpus": 4}')"

    say "${ST_TOTAL} checks, ${ST_FAILS} failed."
    if (( ST_FAILS > 0 )); then done_ 1; fi
    done_ 0
}

# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------
case "${COMMAND}" in
    setup|run|script|clean|status) get_image_tag ;;
esac
case "${COMMAND}" in
    setup)      cmd_setup ;;
    init)       cmd_init ;;
    run)        cmd_run ;;
    status)     cmd_status ;;
    logs)       cmd_logs ;;
    stop)       cmd_stop ;;
    clean)      cmd_clean ;;
    script)     cmd_script ;;
    mem-plan)   cmd_mem_plan ;;
    scan-paths) cmd_scan_paths ;;
    selftest)   cmd_selftest ;;
    *)          cmd_help ;;
esac
