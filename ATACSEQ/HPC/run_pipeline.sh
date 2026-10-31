#!/bin/bash
# run_pipeline.sh — this pipeline's init/run implementation.
#
# Labmates reach this via a pipeline-specific name on their PATH
# (/project/rephilllab/bin/run_atacseq.sh -> this script) — one shared
# PATH entry added once, ever, regardless of how many lab pipelines exist.
# Each pipeline gets its own distinctly-named symlink there, so the command
# name itself (unambiguous to a human or an AI agent driving it) says which
# pipeline is running — no config file needs to be read to find out.
#
#   run_atacseq.sh init             — set up a new project directory here
#   run_atacseq.sh [config.yaml]    — run the pipeline
#
# `init` only scaffolds directories, fetches config.yaml + examples/, and
# git-inits the project — it never runs anything. The plain (no-argument)
# form is the one that submits real compute jobs.
#
# What the run form does:
#   1. If config.yaml's fastq_dir is set, moves that directory's FASTQs into
#      original_data/fastq/ (idempotent — skipped if already staged) and
#      rewrites fastq_dir to point there. If fastq_dir is not set, nothing is
#      staged: samplesheet.csv carries explicit fastq_1/fastq_2 paths,
#      which are used exactly as given.
#   2. Makes sure a shared Nextflow launcher exists (downloads it once, to
#      shared lab storage, if missing — every labmate then runs the same
#      pinned version).
#   3. Submits the actual pipeline run as an LSF job, since the Nextflow
#      head process needs more memory/time than the login node allows and
#      must not run there (see the login-node rules in the project's
#      AGENTS.md, fetched by init).
#
# The head job itself just orchestrates — it submits/polls the real compute
# jobs (bowtie2, samtools, MACS3, etc.) via LSF, each in its own container.

set -euo pipefail

PIPELINE_NAME="atacseq"

# Resolve through symlinks: this script is normally invoked as
# /project/rephilllab/bin/run_atacseq.sh, a symlink into this repo.
# `dirname "${BASH_SOURCE[0]}"` alone would resolve to the symlink's own
# directory (the shared bin/), not this repo — readlink -f follows the
# link to find where main.nf/.nextflow-bin/examples/ actually live.
SCRIPT_REAL_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
PIPELINE_DIR="$(cd "$(dirname "${SCRIPT_REAL_PATH}")" && pwd)"

# Nextflow is pinned: every project on the login node runs this exact version.
# NXF_VER makes the launcher use (and, if needed, fetch) that version even if a
# different binary was downloaded earlier.
NXF_VERSION_PIN="26.04.6"
export NXF_VER="${NXF_VERSION_PIN}"

# ---------------------------------------------------------------------------
# init — scaffold a new project directory, fetch the starter files, git init.
# Never runs the pipeline.
# ---------------------------------------------------------------------------
if [[ "${1:-}" == "init" ]]; then
    echo "Setting up a new ${PIPELINE_NAME} project in $(pwd)..."
    mkdir -p original_data/fastq metadata intermediate_files results utils

    # config.yaml and input CSV templates at the project root (never
    # overwritten; a template missing from examples/ is skipped, not an error).
    for f in config.yaml samplesheet.csv peak_merges.csv; do
        if [[ -f "${f}" ]]; then
            echo "  ${f} already exists — leaving it alone."
        elif [[ -f "${PIPELINE_DIR}/examples/${f}" ]]; then
            cp -n "${PIPELINE_DIR}/examples/${f}" "./${f}"
            echo "  Fetched ${f}"
        else
            echo "  ${f} not found in ${PIPELINE_DIR}/examples — skipped."
        fi
    done

    cp -rn "${PIPELINE_DIR}/examples" ./examples 2>/dev/null || true
    echo "  Fetched examples/ (config.yaml, samplesheet.csv, peak_merges.csv, ready-made projects)"

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

    if [[ -d .git ]]; then
        echo "  Already a git repository — leaving it alone."
    else
        git init -q .
        if [[ ! -f .gitignore ]]; then
            cp "${PIPELINE_DIR}/examples/project.gitignore" ./.gitignore
        fi
        git add config.yaml samplesheet.csv examples .gitignore AGENTS.md CLAUDE.md 2>/dev/null || true
        if [[ -f peak_merges.csv ]]; then
            git add peak_merges.csv 2>/dev/null || true
        fi
        git commit -q -m "Initialize ${PIPELINE_NAME} project" || true
        echo "  Initialized git repository (raw data / results / work dirs excluded — see .gitignore)"
    fi

    echo ""
    echo "Next steps:"
    echo "  1. Edit samplesheet.csv: one row per library, with absolute fastq_1/fastq_2"
    echo "     paths, sample_id, condition, replicate, tmm_group (blank = not normalised),"
    echo "     tmm_group_legacy, peak_group (blank = no peak calling) and exclude (true = skip)."
    echo "  2. Edit peak_merges.csv (optional; one row per merged peak set: name,peak_groups)."
    echo "  3. Edit config.yaml: project_name, genome_index, blacklist, tss_bed,"
    echo "     effective_genome_size and tmm_group_column (tmm_group | tmm_group_legacy)."
    echo "  4. Run:  run_${PIPELINE_NAME}.sh"
    exit 0
fi

# ---------------------------------------------------------------------------
# Run — everything below here submits real compute jobs.
# ---------------------------------------------------------------------------
CONFIG_FILE="${1:-config.yaml}"
PROJECT_DIR="$(pwd)"
NXF_BIN_DIR="${PIPELINE_DIR}/.nextflow-bin"
NXF_BIN="${NXF_BIN_DIR}/nextflow"
ORIGINAL_DATA_DIR="original_data"

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
yaml_set() {
    local key="$1" value="$2"
    sed -i -E "s|^${key}:.*|${key}: ${value}|" "${CONFIG_FILE}"
}

PROJECT_NAME="$(yaml_get "project_name")"
PROJECT_NAME="${PROJECT_NAME:-unnamed}"

# ---- Optional: stage a FASTQ directory, idempotently ------------------------
# samplesheet.csv carries explicit fastq_1/fastq_2 paths, which are used as
# given. Only when fastq_dir is set in config.yaml: by convention FASTQs live
# in original_data/fastq/, so unless fastq_dir already is that directory, its
# *.fastq.gz / *.fastq files are moved there first and config.yaml is
# rewritten to point at it.
FASTQ_DEST="${ORIGINAL_DATA_DIR}/fastq"
stage_fastq_dir() {
    local current
    current="$(yaml_get "fastq_dir")"
    if [[ -z "${current}" || "${current}" == "null" ]]; then
        return 1   # not in fastq_dir mode
    fi
    if [[ ! -d "${current}" ]]; then
        echo "ERROR: fastq_dir in ${CONFIG_FILE} is not a directory: ${current}" >&2
        exit 1
    fi
    mkdir -p "${FASTQ_DEST}"
    if [[ "$(cd "${current}" && pwd -P)" == "$(cd "${FASTQ_DEST}" && pwd -P)" ]]; then
        echo "  fastq_dir: already ${FASTQ_DEST}"
        return 0
    fi

    local files=()
    shopt -s nullglob
    files=( "${current}"/*.fastq.gz "${current}"/*.fastq )
    shopt -u nullglob
    if (( ${#files[@]} == 0 )); then
        echo "ERROR: fastq_dir ${current} contains no *.fastq.gz or *.fastq files." >&2
        exit 1
    fi
    local f clash=()
    for f in "${files[@]}"; do
        [[ -e "${FASTQ_DEST}/$(basename "${f}")" ]] && clash+=( "$(basename "${f}")" )
    done
    if (( ${#clash[@]} > 0 )); then
        echo "ERROR: these files already exist in ${FASTQ_DEST}/ — not overwriting them:" >&2
        printf '    %s\n' "${clash[@]}" >&2
        exit 1
    fi
    echo "  Moving ${#files[@]} FASTQ file(s) from ${current} into ${FASTQ_DEST}/"
    for f in "${files[@]}"; do
        if ! mv "${f}" "${FASTQ_DEST}/" 2>/dev/null; then
            echo "  [WARN] mv failed for ${f} (read-only source, or cross-device?) — copying instead" >&2
            cp "${f}" "${FASTQ_DEST}/"
        fi
    done
    yaml_set "fastq_dir" "${FASTQ_DEST}"
}

if ! stage_fastq_dir; then
    echo "fastq_dir not set: samplesheet.csv fastq_1/fastq_2 paths are used as given."
fi

# ---- One-time shared Nextflow install --------------------------------------
if [[ ! -x "${NXF_BIN}" ]]; then
    echo "Nextflow launcher not found — downloading once to ${NXF_BIN_DIR} (shared by all labmates)..."
    mkdir -p "${NXF_BIN_DIR}"
    curl -fsSL -o "${NXF_BIN}" "https://github.com/nextflow-io/nextflow/releases/download/v26.04.6/nextflow"
    chmod +x "${NXF_BIN}"
fi

# ---- Where LSF should write logs for the head job --------------------------
# Lives under metadata/ (not the project root) — this is Nextflow's own
# run bookkeeping, not a scientific input or output. .nextflow/ (Nextflow's
# resume cache/history) is the one thing that can't move here: it always
# sits at the directory Nextflow was launched from, and moving the launch
# directory itself would break every relative path in config.yaml.
LOG_DIR="${PROJECT_DIR}/metadata/nf_logs"
mkdir -p "${LOG_DIR}"
TS="$(date +%Y%m%d_%H%M%S)"

# Shared, group-writable scratch for Nextflow's intermediate work files —
# keeps large per-task staging off /home and off node-local /tmp (which
# would break -resume once a different node picks up the next task).
export NXF_WORK="/project/rephilllab/nf_work/${USER}/$(basename "${PROJECT_DIR}")"
mkdir -p "${NXF_WORK}"

# LSF dispatch watchdog (lib/lsf_watchdog.sh) runs in the head job unless
# config.yaml sets `lsf_watchdog: false`. Read here, expanded into the heredoc.
LSF_WATCHDOG_ON=1
if [[ "$(yaml_get "lsf_watchdog")" == "false" ]]; then
    LSF_WATCHDOG_ON=0
fi

echo "Submitting pipeline head job to LSF (queue: rhel9)..."
echo "  Pipeline : ${PIPELINE_DIR}"
echo "  Config   : ${PROJECT_DIR}/${CONFIG_FILE}"
echo "  Work dir : ${NXF_WORK}"
echo "  Logs     : ${LOG_DIR}/nextflow_${TS}.out / .err"

bsub <<EOF
#BSUB -J "${PIPELINE_NAME}_${PROJECT_NAME}_${TS}"
#BSUB -q rhel9
#BSUB -oo "${LOG_DIR}/nextflow_${TS}.out"
#BSUB -eo "${LOG_DIR}/nextflow_${TS}.err"
#BSUB -n 1
#BSUB -M 8000
#BSUB -R "rusage[mem=8000]"
#BSUB -W 168:00

cd "${PROJECT_DIR}"
mkdir -p metadata
export NXF_WORK="${NXF_WORK}"
export NXF_VER="${NXF_VERSION_PIN}"
export NXF_ANSI_LOG=false
"${NXF_BIN}" -log metadata/.nextflow.log run "${PIPELINE_DIR}" -profile pmacs -params-file "${CONFIG_FILE}" -resume &
NXF_PID=\$!
# Stops the run if an LSF job is aborted at dispatch (Nextflow would hang on it).
if [[ "${LSF_WATCHDOG_ON}" == "1" ]]; then
    "${PIPELINE_DIR}/lib/lsf_watchdog.sh" metadata/.nextflow.log "\${NXF_PID}" &
    WATCHDOG_PID=\$!
fi
wait "\${NXF_PID}"
NXF_STATUS=\$?
[[ -n "\${WATCHDOG_PID:-}" ]] && kill "\${WATCHDOG_PID}" 2>/dev/null
exit "\${NXF_STATUS}"
EOF

echo ""
echo "Submitted. Monitor with:  bjobs -u \$USER"
echo "Follow progress with:     tail -f ${LOG_DIR}/nextflow_${TS}.out"
