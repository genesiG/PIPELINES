#!/bin/bash
# make_locks.sh - (re)generate docker/omics/env/<name>.lock.txt from <name>.yml (developers only).
#
# Each <name>.yml is a human-readable spec (channels conda-forge, bioconda; strict
# priority; no defaults). This script solves it in a pinned micromamba container
# (linux-64) and writes the exact solution (build, URL, md5) as a conda explicit
# lock, which the Dockerfile installs verbatim.
#
#   bash docker/omics/env/make_locks.sh              lock every env/*.yml
#   bash docker/omics/env/make_locks.sh r qc         lock only these
#   bash docker/omics/env/make_locks.sh --verify [names]
#                                                    install each existing lock in a
#                                                    fresh container (must succeed)
#
# Needs Docker + bash (Git Bash on Windows, or Linux/macOS).
set -euo pipefail
IMAGE="mambaorg/micromamba:2.0.5"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MSYS_NO_PATHCONV=1   # Git Bash on Windows: don't rewrite /env etc.

mode=lock
if [[ "${1:-}" == "--verify" ]]; then mode=verify; shift; fi
names=("$@")
if [[ ${#names[@]} -eq 0 ]]; then
    for f in "${here}"/*.yml; do names+=("$(basename "${f}" .yml)"); done
fi

hostdir="${here}"
command -v cygpath >/dev/null 2>&1 && hostdir="$(cygpath -w "${here}")"

# Run a bash snippet in the micromamba container; env/ is mounted at /env.
# A named volume caches downloaded packages between runs and retries.
run() {
    docker run --rm --platform linux/amd64 -u root \
        -e CONDA_SUBDIR=linux-64 -e MAMBA_ROOT_PREFIX=/opt/micromamba \
        -v "${hostdir}:/env" -v omics-locks-pkgs:/opt/micromamba/pkgs \
        --entrypoint /bin/bash "${IMAGE}" -c "$1"
}

retry() {   # retry <label> <command...>: docker/micromamba occasionally abort ("netlink descriptor")
    local label="$1" attempt; shift
    for attempt in 1 2 3 4 5; do
        "$@" && return 0
        echo "   ${label} failed (attempt ${attempt} of 5), retrying" >&2; sleep 10
    done
    return 1
}

solve() {   # solve <name>
    local name="$1" out="${here}/${name}.lock.txt.part"
    # Solve only (dry run; nothing is installed, so no post-link scripts run). In the JSON
    # each LINK record lists "md5" before "url": print url#md5. The solver's order (dependencies
    # first) is kept, NOT sorted: an explicit install runs in file order, and the post-link
    # scripts of Bioconductor data packages need yq/curl/R linked before them.
    run "micromamba create -n x --strict-channel-priority --dry-run --json -f /env/${name}.yml | awk -F'\"' '/\"md5\":/{m=\$4} /\"url\":/{print \$4 \"#\" m}'" > "${out}"
    [[ -s "${out}" ]]
}

for name in "${names[@]}"; do
    if [[ "${mode}" == verify ]]; then
        echo "== verify ${name}"
        retry "verify ${name}" run "micromamba create -y -q -p /tmp/t -f /env/${name}.lock.txt && echo OK: ${name} installs"
        continue
    fi
    echo "== lock ${name}"
    [[ -f "${here}/${name}.yml" ]] || { echo "no ${name}.yml" >&2; exit 1; }
    retry "solve ${name}" solve "${name}"
    {
        echo "# Pinned package list for the '${name}' environment of the omics Docker image."
        echo "# Solved from ${name}.yml (conda-forge, bioconda; strict priority; linux-64)"
        echo "# by docker/omics/env/make_locks.sh in ${IMAGE}."
        echo "@EXPLICIT"
        cat "${here}/${name}.lock.txt.part"
    } > "${here}/${name}.lock.txt"
    rm -f "${here}/${name}.lock.txt.part"
    echo "   $(grep -c '^http' "${here}/${name}.lock.txt") packages -> ${name}.lock.txt"
done
