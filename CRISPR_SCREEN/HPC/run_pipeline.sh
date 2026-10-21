#!/bin/bash
# run_pipeline.sh — this pipeline's init/run implementation.
#
# Labmates reach this via a pipeline-specific name on their PATH
# (/project/rephilllab/bin/run_crisprscreen.sh -> this script) — one shared
# PATH entry added once, ever, regardless of how many lab pipelines exist.
# Each pipeline gets its own distinctly-named symlink there, so the command
# name itself (unambiguous to a human or an AI agent driving it) says which
# pipeline is running — no config file needs to be read to find out.
#
#   run_crisprscreen.sh init             — set up a new project directory here
#   run_crisprscreen.sh [config.yaml]    — run the pipeline
#   run_crisprscreen.sh script [--cpus N] [--mem GB] [--time HH:MM] utils/X.R [ARGS]
#                                        — run one of the project's own scripts
#                                          (.R, .py or .sh) in the pipeline's software
#                                          environment as an LSF job
#
# `script` submits one LSF job that runs the script under Apptainer, from the
# project folder. The flags (defaults 1 core, 8 GB, 4:00) are read only between
# `script` and the script path (an optional `--` may follow them); everything after
# the script path is passed to the script verbatim. The script must be inside the
# project folder. Its log goes to metadata/script_logs/ in the project.
#
# `init` only scaffolds directories, fetches config.yaml + examples/, the
# default barcode-position plot script into utils/ (its md5 goes in
# utils/.pipeline_defaults.md5) and agent/.claude/ into .claude/, then
# git-inits the project and commits them — it never runs anything. The plain
# (no-argument) form is the one that submits real compute jobs.
#
# What the run form does:
#   1. Checks the FASTQ paths named in config.yaml (fastq_r1/fastq_r2, or
#      fastq_dir). The raw FASTQs stay where they are: never moved, copied or
#      renamed, and config.yaml is not rewritten. The paths are resolved to
#      absolute paths and handed to Nextflow.
#   2. Makes sure a shared Nextflow launcher exists (downloads it once, to
#      shared lab storage, if missing — every labmate then runs the same
#      pinned version).
#   3. Submits the actual pipeline run as an LSF job, since the Nextflow
#      head process needs more memory/time than the login node allows and
#      must not run there (see AGENTS.md rule 1 in any project built on
#      this pipeline).
#
# The head job itself just orchestrates — it submits/polls the real compute
# jobs (cutadapt, mageck, etc.) via LSF, each in its own container. After a
# successful run it deletes the Nextflow work folder and the .nextflow/ cache
# unless keep_work: true in config.yaml (a failed or stopped run keeps both so
# a re-run resumes).

set -euo pipefail

PIPELINE_NAME="crisprscreen"

# Resolve through symlinks: this script is normally invoked as
# /project/rephilllab/bin/run_crisprscreen.sh, a symlink into this repo.
# `dirname "${BASH_SOURCE[0]}"` alone would resolve to the symlink's own
# directory (the shared bin/), not this repo — readlink -f follows the
# link to find where main.nf/.nextflow-bin/examples/ actually live.
SCRIPT_REAL_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
PIPELINE_DIR="$(cd "$(dirname "${SCRIPT_REAL_PATH}")" && pwd)"

# ---------------------------------------------------------------------------
# init — scaffold a new project directory, fetch the starter files, git init.
# Never runs the pipeline.
# ---------------------------------------------------------------------------
if [[ "${1:-}" == "init" ]]; then
    echo "Setting up a new ${PIPELINE_NAME} project in $(pwd)..."
    mkdir -p original_data/fastq metadata results utils

    if [[ -f config.yaml ]]; then
        echo "  config.yaml already exists — leaving it alone."
    else
        cp "${PIPELINE_DIR}/examples/config.yaml" ./config.yaml
        echo "  Fetched config.yaml"
    fi

    mkdir -p examples && cp -rn "${PIPELINE_DIR}/examples/." ./examples/ 2>/dev/null || true   # contents only: re-running init must not nest examples/examples/
    echo "  Fetched examples/ (samplesheet_example.csv, library_example.csv, contrasts_example.csv)"

    # Default barcode-position plot script. Never overwrites an existing copy;
    # the md5 of what was fetched is recorded so later edits can be detected.
    utils_src="${PIPELINE_DIR}/bin/plot_barcode_position.R"
    utils_dst="utils/plot_barcode_position.R"
    if [[ -f "${utils_src}" ]]; then
        if [[ -e "${utils_dst}" ]]; then
            echo "  ${utils_dst} already exists — leaving it alone."
        else
            cp "${utils_src}" "${utils_dst}"
            utils_md5="$(md5sum "${utils_dst}" | awk '{print $1}')"
            md5_file="utils/.pipeline_defaults.md5"
            md5_tmp="utils/.pipeline_defaults.md5.tmp.$$"
            touch "${md5_file}"
            awk -v p="${utils_dst}" '{ if (substr($0, 35) != p) print }' "${md5_file}" > "${md5_tmp}"
            printf '%s  %s\n' "${utils_md5}" "${utils_dst}" >> "${md5_tmp}"
            cat "${md5_tmp}" > "${md5_file}"; rm -f "${md5_tmp}"
            echo "  Fetched ${utils_dst} (default barcode-position plot; md5 recorded in utils/.pipeline_defaults.md5)"
        fi
    fi

    # Agent settings for Claude Code (.claude/), copied file by file so that
    # existing files are never overwritten.
    claude_src="${PIPELINE_DIR}/agent/.claude"
    if [[ -d "${claude_src}" ]]; then
        n=0
        claude_list="$(mktemp)"
        (cd "${claude_src}" && find . -type f) > "${claude_list}"
        while IFS= read -r rel; do
            rel="${rel#./}"
            mkdir -p ".claude/$(dirname "${rel}")"
            if [[ ! -e ".claude/${rel}" ]]; then
                cp "${claude_src}/${rel}" ".claude/${rel}"
                n=$((n+1))
            fi
        done < "${claude_list}"
        rm -f "${claude_list}"
        if (( n > 0 )); then
            echo "  Fetched .claude/ (${n} files for AI agents; existing files left alone)"
        fi
    fi

    # Instructions for AI agents working in this project (how to run, monitor,
    # customize; what never to touch). Never overwrites an existing copy.
    for f in AGENTS.md CLAUDE.md; do
        if [[ -f "${f}" ]]; then
            echo "  ${f} already exists — leaving it alone."
        else
            cp "${PIPELINE_DIR}/agent/${f}" "./${f}"
            echo "  Fetched ${f}"
        fi
    done

    # .gitattributes: LF line endings for the scripts, so a utils/ copy stays
    # byte-identical to the pipeline's default on every OS. Never overwritten.
    if [[ ! -e .gitattributes ]]; then
        cat > .gitattributes <<'EOF'
# Keep scripts byte-identical to the pipeline's copies on every OS: utils/
# replacements are compared with the pipeline's default scripts, and bash
# scripts need LF line endings.
utils/** text eol=lf
*.sh text eol=lf
*.R text eol=lf
*.py text eol=lf
EOF
        echo "  Fetched .gitattributes (LF line endings for scripts)"
    fi

    if [[ -d .git ]]; then
        echo "  Already a git repository — leaving it alone."
    else
        git init -q .
        if [[ ! -f .gitignore ]]; then
            cp "${PIPELINE_DIR}/examples/project.gitignore" ./.gitignore
        fi
        git_paths=( config.yaml examples .gitignore .gitattributes AGENTS.md CLAUDE.md utils )
        if [[ -d .claude ]]; then
            git_paths+=( .claude )
        fi
        git add "${git_paths[@]}" 2>/dev/null || true
        git commit -q -m "Initialize ${PIPELINE_NAME} project" || true
        echo "  Initialized git repository (raw data / results / work dirs excluded — see .gitignore; utils/ and .claude/ are included in the first commit)"
    fi

    echo ""
    echo "Next steps:"
    echo "  1. Copy the example CSVs you need from examples/ into this directory"
    echo "     (samplesheet.csv, library.csv, and — if you're ready to run mageck"
    echo "     test — contrasts.csv) and fill in your real data."
    echo "  2. Edit config.yaml: project_name, and fastq_r1 (multiplexed FASTQ) or"
    echo "     fastq_dir (folder of already-demultiplexed FASTQs). Your FASTQ files stay"
    echo "     where they are: put their full path in config.yaml and keep it stable."
    echo "  3. Run:  run_${PIPELINE_NAME}.sh"
    exit 0
fi

# ---------------------------------------------------------------------------
# script — run one of the project's own scripts in the pipeline's software
# environment, as one LSF job. Nothing runs on the login node: only bsub is
# called here. Needs no config.yaml.
# ---------------------------------------------------------------------------
if [[ "${1:-}" == "script" ]]; then
    shift
    PROJECT_DIR="$(pwd)"
    SCRIPT_CPUS=1
    SCRIPT_MEM=8
    SCRIPT_TIME="4:00"

    # Flags are read only between `script` and the script path.
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --cpus|--mem|--time)
                if [[ $# -lt 2 ]]; then
                    echo "ERROR: $1 needs a value." >&2
                    exit 1
                fi
                case "$1" in
                    --cpus) SCRIPT_CPUS="$2" ;;
                    --mem)  SCRIPT_MEM="$2" ;;
                    *)      SCRIPT_TIME="$2" ;;
                esac
                shift 2 ;;
            --) shift; break ;;
            *) break ;;
        esac
    done

    if [[ ! "${SCRIPT_CPUS}" =~ ^[0-9]+$ ]] || (( 10#${SCRIPT_CPUS} < 1 )); then
        echo "ERROR: --cpus must be an integer >= 1 (got '${SCRIPT_CPUS}')." >&2
        exit 1
    fi
    if [[ ! "${SCRIPT_MEM}" =~ ^[0-9]+$ ]] || (( 10#${SCRIPT_MEM} < 1 )); then
        echo "ERROR: --mem must be an integer number of GB >= 1 (got '${SCRIPT_MEM}')." >&2
        exit 1
    fi
    if [[ ! "${SCRIPT_TIME}" =~ ^[0-9]+:[0-9]{2}$ ]]; then
        echo "ERROR: --time must be HH:MM (got '${SCRIPT_TIME}')." >&2
        exit 1
    fi

    if [[ $# -eq 0 ]]; then
        echo "ERROR: usage: run_${PIPELINE_NAME}.sh script [--cpus N] [--mem GB] [--time HH:MM] utils/<name>.<R|py|sh> [args...]" >&2
        exit 1
    fi
    SCRIPT_ARG="$1"
    shift

    # ---- The script: inside the project folder, a regular file, known type --
    if [[ "${SCRIPT_ARG}" == /* ]]; then
        SCRIPT_ABS="${SCRIPT_ARG}"
    else
        SCRIPT_ABS="${PROJECT_DIR}/${SCRIPT_ARG}"
    fi
    if [[ ! -e "${SCRIPT_ABS}" ]]; then
        echo "ERROR: script not found: ${SCRIPT_ARG} (looked for ${SCRIPT_ABS})" >&2
        exit 1
    fi
    if [[ ! -f "${SCRIPT_ABS}" ]]; then
        echo "ERROR: ${SCRIPT_ARG} is not a regular file." >&2
        exit 1
    fi
    SCRIPT_REAL="$(readlink -f "${SCRIPT_ABS}")"
    PROJECT_REAL="$(readlink -f "${PROJECT_DIR}")"
    if [[ "${SCRIPT_REAL}" != "${PROJECT_REAL}/"* ]]; then
        echo "ERROR: the script must be inside the project folder" >&2
        exit 1
    fi
    SCRIPT_REL="${SCRIPT_REAL#"${PROJECT_REAL}/"}"
    SCRIPT_BASE="$(basename "${SCRIPT_REL}")"
    SCRIPT_STEM="${SCRIPT_BASE%.*}"

    # ---- Runner and container: the same images as the pipeline's own steps --
    case "${SCRIPT_BASE}" in
        *.R|*.r) RUNNER="Rscript"; IMAGE="rocker/tidyverse:4.3" ;;
        *.py)    RUNNER="python3"; IMAGE="quay.io/biocontainers/cutadapt:5.2--py39hbcbf7aa_0" ;;
        *.sh)    RUNNER="bash";    IMAGE="rocker/tidyverse:4.3" ;;
        *)
            echo "ERROR: unsupported extension (use .R, .py or .sh)" >&2
            exit 1 ;;
    esac

    # Nextflow's Apptainer cache (apptainer.cacheDir in conf/pmacs.config) stores
    # each image as <image with / and : replaced by ->.img. If it is not there,
    # the job pulls it once from docker:// into the same cache.
    CACHE_DIR="/project/rephilllab/nf_containers"
    IMG_STEM="$(printf '%s' "${IMAGE}" | tr '/:' '--')"
    IMG_FILE=""
    if [[ -f "${CACHE_DIR}/${IMG_STEM}.img" ]]; then
        IMG_FILE="${CACHE_DIR}/${IMG_STEM}.img"
    else
        shopt -s nullglob
        cached=( "${CACHE_DIR}"/*"${IMG_STEM}"* )
        shopt -u nullglob
        for cand in "${cached[@]}"; do
            if [[ -f "${cand}" ]]; then
                IMG_FILE="${cand}"
                break
            fi
        done
    fi
    if [[ -n "${IMG_FILE}" ]]; then
        IMG_REF="${IMG_FILE}"
        IMG_NOTE="${IMG_FILE}"
    else
        IMG_REF="docker://${IMAGE}"
        IMG_NOTE="not cached in ${CACHE_DIR}: pulling docker://${IMAGE} once into APPTAINER_CACHEDIR"
    fi

    # ---- Arguments for the script: quoted for the job, shown in the log ------
    ARGS_Q=""
    ARGS_LOG="(none)"
    if [[ $# -gt 0 ]]; then
        ARGS_Q="$(printf '%q ' "$@")"
        ARGS_LOG="${ARGS_Q% }"
    fi

    MEM_MB=$(( 10#${SCRIPT_MEM} * 1000 ))
    TS="$(date +%Y%m%d_%H%M%S)"
    LOG_DIR="${PROJECT_DIR}/metadata/script_logs"
    LOG="${LOG_DIR}/${SCRIPT_STEM}_${TS}.log"
    mkdir -p "${LOG_DIR}"

    echo "Submitting script job to LSF (queue: rhel9)..."
    echo "  Script   : ${PROJECT_DIR}/${SCRIPT_REL}"
    echo "  Resources: ${SCRIPT_CPUS} cores, ${SCRIPT_MEM} GB, ${SCRIPT_TIME} h"
    echo "  Image    : ${IMG_NOTE}"

    bsub_out="$(bsub <<EOF
#BSUB -J "${PIPELINE_NAME}_script_${SCRIPT_STEM}_${TS}"
#BSUB -q rhel9
#BSUB -oo "${LOG%.log}.lsf.out"
#BSUB -eo "${LOG%.log}.lsf.err"
#BSUB -n ${SCRIPT_CPUS}
#BSUB -M ${MEM_MB}
#BSUB -R "rusage[mem=${MEM_MB}]"
#BSUB -W ${SCRIPT_TIME}

cd "${PROJECT_DIR}"
mkdir -p metadata/script_logs
export APPTAINER_CACHEDIR="${CACHE_DIR}"
export CRISPR_MAX_CPUS=${SCRIPT_CPUS}
export CRISPR_MAX_MEMORY_MB=${MEM_MB}
{
echo "CRISPR_SCREEN script run (cluster)"
echo "date      : \$(date)"
echo "project   : ${PROJECT_DIR}"
echo "script    : ${SCRIPT_REL}"
echo "md5       : \$(md5sum "${SCRIPT_REL}" | awk '{print \$1}')"
echo "args      : ${ARGS_LOG}"
echo "container : ${IMAGE}"
echo "image file: ${IMG_NOTE}"
echo "resources : ${SCRIPT_CPUS} cores, ${SCRIPT_MEM} GB, ${SCRIPT_TIME} h"
echo "LSF job   : \$LSB_JOBID"
echo "host      : \$(hostname)"
echo ""
echo "------------------------------------------------------------"
} >> "${LOG}"
apptainer exec -B /usr/lib64/ossl-modules --bind "${PROJECT_DIR}" --pwd "${PROJECT_DIR}" "${IMG_REF}" ${RUNNER} "${SCRIPT_REL}" ${ARGS_Q} >> "${LOG}" 2>&1
script_rc=\$?
printf '\n-- exit code %s --\n' "\$script_rc" >> "${LOG}"
exit \$script_rc
EOF
)"
    echo "${bsub_out}"
    job_id="$(printf '%s\n' "${bsub_out}" | sed -n 's/^Job <\([0-9]*\)>.*/\1/p' | head -1)"
    if [[ -z "${job_id}" ]]; then
        echo "ERROR: could not read the job id from bsub's output (see above)." >&2
        exit 1
    fi

    echo ""
    echo "Submitted job ${job_id}."
    echo "  Log      : ${LOG}"
    echo "  LSF logs : ${LOG%.log}.lsf.out / ${LOG%.log}.lsf.err"
    echo "Follow it with:"
    echo "  tail -f ${LOG}"
    echo "  bjobs ${job_id}"
    echo "  bpeek ${job_id}"
    exit 0
fi

# ---------------------------------------------------------------------------
# Run — everything below here submits real compute jobs.
# ---------------------------------------------------------------------------
CONFIG_FILE="${1:-config.yaml}"
PROJECT_DIR="$(pwd)"
NXF_BIN_DIR="${PIPELINE_DIR}/.nextflow-bin"
NXF_BIN="${NXF_BIN_DIR}/nextflow"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "ERROR: ${CONFIG_FILE} not found in $(pwd)." >&2
    echo "Run 'run_${PIPELINE_NAME}.sh init' first." >&2
    exit 1
fi

# ---- Minimal flat-YAML get/set for config.yaml's top-level keys -----------
# config.yaml is intentionally kept flat (no nested keys) specifically so
# this works without a YAML-parsing dependency — every labmate's account
# has bash, nothing else required.
yaml_get() {
    local key="$1"
    grep -E "^${key}:" "${CONFIG_FILE}" | head -1 \
        | sed -E "s/^${key}:[[:space:]]*//; s/[[:space:]]*(#.*)?\$//; s/^['\"]//; s/['\"]\$//"
}

PROJECT_NAME="$(yaml_get "project_name")"
PROJECT_NAME="${PROJECT_NAME:-unnamed}"

# keep_work: true keeps the Nextflow work folder after a successful run (for
# -resume); anything else, or no such key, deletes it (see the job script below).
KEEP_WORK="$(yaml_get "keep_work" | tr 'A-Z' 'a-z' || true)"
case "${KEEP_WORK}" in true|yes) KEEP_WORK=true ;; *) KEEP_WORK=false ;; esac
# results/ only stays valid without work/ when steps publish copies (the default).
case "$(yaml_get "publish_mode" | tr -d "\"'" || true)" in symlink|rellink) KEEP_WORK=true ;; esac
# publish_mode: auto is the Docker versions' value; on the cluster it means copy
# (results/ stays standalone, as with the default). Any other value passes nothing.
EXTRA_ARGS=""
if [[ "$(yaml_get "publish_mode" | tr -d "\"'" || true)" == "auto" ]]; then
    EXTRA_ARGS="--publish_mode copy "
fi

# ---- Raw FASTQ paths: check them where they are, never move or copy them ---
# The raw FASTQs stay where the user keeps them, and config.yaml is only read,
# never rewritten. Each path is checked, made absolute (leading ~ expanded,
# relative paths taken from the project folder) and passed to Nextflow on its
# command line, which overrides the same keys in -params-file.
# config_path KEY: the value of a path key, or empty when it is unset.
config_path() {
    local value
    value="$(yaml_get "$1" || true)"
    case "${value}" in
        ""|null|"~") value="" ;;
    esac
    printf '%s' "${value}"
}

# resolve_input_path KEY file|dir: checks the path, prints it and sets
# RESOLVED_ABS (left empty when the key is unset).
resolve_input_path() {
    local key="$1" kind="$2" raw abs files=()
    RESOLVED_ABS=""
    raw="$(config_path "${key}")"
    if [[ -z "${raw}" ]]; then
        return 0   # not set (e.g. fastq_r2 for single-end) — nothing to check
    fi
    abs="${raw}"
    if [[ "${abs}" == "~/"* ]]; then
        abs="${HOME}/${abs:2}"
    fi
    if [[ "${abs}" != /* ]]; then
        abs="${PROJECT_DIR}/${abs}"
    fi
    if [[ "${kind}" == "dir" ]]; then
        if [[ ! -d "${abs}" ]]; then
            echo "ERROR: ${key} in ${CONFIG_FILE} points to a folder that doesn't exist: ${raw} (looked for ${abs})" >&2
            exit 1
        fi
        shopt -s nullglob
        files=( "${abs}"/*.fastq.gz "${abs}"/*.fastq )
        shopt -u nullglob
        if (( ${#files[@]} == 0 )); then
            echo "ERROR: ${key} ${abs} contains no *.fastq.gz or *.fastq files." >&2
            exit 1
        fi
    elif [[ ! -f "${abs}" ]]; then
        echo "ERROR: ${key} in ${CONFIG_FILE} points to a file that doesn't exist: ${raw} (looked for ${abs})" >&2
        exit 1
    fi
    # The folder is resolved, the name kept as given: a symlink stays a symlink.
    RESOLVED_ABS="$(cd "$(dirname "${abs}")" && pwd)/$(basename "${abs}")"
    echo "  ${key}: ${RESOLVED_ABS} (stays where it is)"
}

# ---- Which FASTQ input the run uses -----------------------------------------
# fastq_dir set means the samples are already demultiplexed: demultiplexing is
# skipped and fastq_r1/fastq_r2 are ignored (not checked, not passed), whatever
# they say. Otherwise fastq_r1 is required and fastq_r2 is checked if set.
# FASTQ_ARGS holds the absolute paths as shell-quoted Nextflow parameters.
FASTQ_ARGS=""
resolve_input_path "fastq_dir" "dir"
if [[ -n "${RESOLVED_ABS}" ]]; then
    if [[ -n "$(config_path "fastq_r1")" ]]; then
        echo "  fastq_dir is set, so the samples are treated as already demultiplexed (fastq_r1 is ignored)."
    fi
    FASTQ_ARGS+="$(printf -- '--fastq_dir %q ' "${RESOLVED_ABS}")"
else
    resolve_input_path "fastq_r1" "file"
    if [[ -z "${RESOLVED_ABS}" ]]; then
        echo "ERROR: config.yaml needs fastq_r1 (a multiplexed FASTQ) or fastq_dir (a folder of already-demultiplexed FASTQs)." >&2
        exit 1
    fi
    FASTQ_ARGS+="$(printf -- '--fastq_r1 %q ' "${RESOLVED_ABS}")"
    resolve_input_path "fastq_r2" "file"
    if [[ -n "${RESOLVED_ABS}" ]]; then
        FASTQ_ARGS+="$(printf -- '--fastq_r2 %q ' "${RESOLVED_ABS}")"
    fi
fi

# ---- One-time shared Nextflow install --------------------------------------
if [[ ! -x "${NXF_BIN}" ]]; then
    echo "Nextflow launcher not found — downloading once to ${NXF_BIN_DIR} (shared by all labmates)..."
    mkdir -p "${NXF_BIN_DIR}"
    curl -fsSL -o "${NXF_BIN}" https://get.nextflow.io
    chmod +x "${NXF_BIN}"
fi

# ---- Where LSF should write logs for the head job --------------------------
# Lives under metadata/ (not the project root) — this is Nextflow's own
# run bookkeeping, not a scientific input or output. .nextflow/ (Nextflow's
# resume cache/history) is the one thing that can't move here: it always
# sits at the directory Nextflow was launched from, and moving the launch
# directory itself would break every relative path in config.yaml. After a
# successful run with keep_work: false it is deleted together with the work
# folder (see the job script below).
LOG_DIR="${PROJECT_DIR}/metadata/nf_logs"
mkdir -p "${LOG_DIR}"
TS="$(date +%Y%m%d_%H%M%S)"

# Shared, group-writable scratch for Nextflow's intermediate work files —
# keeps large per-task staging off /home and off node-local /tmp (which
# would break -resume once a different node picks up the next task).
export NXF_WORK="/project/rephilllab/nf_work/${USER}/$(basename "${PROJECT_DIR}")"
mkdir -p "${NXF_WORK}"

echo "Submitting pipeline head job to LSF (queue: rhel9)..."
echo "  Pipeline : ${PIPELINE_DIR}"
echo "  Config   : ${PROJECT_DIR}/${CONFIG_FILE}"
echo "  Work dir : ${NXF_WORK}"
echo "  Keep work: ${KEEP_WORK}"
echo "  Logs     : ${LOG_DIR}/nextflow_${TS}.out / .err"

# FASTQ_ARGS is already shell-quoted (printf %q), so the job script gets it as-is.
bsub <<EOF
#BSUB -J "${PIPELINE_NAME}_${PROJECT_NAME}_${TS}"
#BSUB -q rhel9
#BSUB -oo "${LOG_DIR}/nextflow_${TS}.out"
#BSUB -eo "${LOG_DIR}/nextflow_${TS}.err"
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "rusage[mem=8000]"
#BSUB -W 72:00

cd "${PROJECT_DIR}"
mkdir -p metadata
export NXF_WORK="${NXF_WORK}"
export NXF_ANSI_LOG=false
"${NXF_BIN}" -log metadata/.nextflow.log run "${PIPELINE_DIR}" -profile pmacs -params-file "${CONFIG_FILE}" -resume ${FASTQ_ARGS}${EXTRA_ARGS}
nf_rc=\$?
# Success + keep_work false: free the disk (publish_mode is 'copy', so results/ does not point into work/)
# The guard needs two components below nf_work (.../nf_work/<user>/<project>), so it can never match nf_work or the user folder.
if [ "\$nf_rc" -eq 0 ] && [ "${KEEP_WORK}" != "true" ]; then
  freed=""
  case "${NXF_WORK}" in
    /project/rephilllab/nf_work/?*/?*)
      if [ -d "${NXF_WORK}" ]; then
        kb=\$(du -sk "${NXF_WORK}" 2>/dev/null | cut -f1)
        gb=\$(awk -v k="\${kb:-0}" 'BEGIN{printf "%.1f", k/1048576}')
        find "${NXF_WORK}" -mindepth 1 -delete 2>/dev/null || true
        rmdir "${NXF_WORK}" 2>/dev/null || true
        freed="the work folder (\${gb} GB)"
      fi ;;
  esac
  # Nextflow's resume cache, at the project root (the directory Nextflow was launched from).
  if [ -e "${PROJECT_DIR}/.nextflow" ]; then
    rm -rf "${PROJECT_DIR}/.nextflow"
    if [ -n "\$freed" ]; then freed="\$freed and "; fi
    freed="\${freed}Nextflow's .nextflow/ cache"
  fi
  rm -f "${PROJECT_DIR}/.nextflow.log"*
  if [ -n "\$freed" ]; then
    echo "Deleted \$freed to free space; set keep_work: true in ${CONFIG_FILE} to keep them for resuming."
  fi
fi
exit \$nf_rc
EOF

echo ""
echo "Submitted. Monitor with:  bjobs -u \$USER"
echo "Follow progress with:     tail -f ${LOG_DIR}/nextflow_${TS}.out"
