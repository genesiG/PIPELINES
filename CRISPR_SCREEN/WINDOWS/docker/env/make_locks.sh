#!/bin/bash
# make_locks.sh — (re)generate docker/env/<tool>.explicit.txt (developers only).
#
# The lock files pin every package of each tool's environment (exact build,
# URL and md5), so the image is the same on every computer, whenever it is
# built. They are taken from the very BioContainers images the cluster
# version runs (images.txt: tool -> image), by reading each image's
# conda-meta/ records: the Docker image then holds the same package builds as
# the cluster's per-step containers, not merely the same tool version.
#
#   bash docker/env/make_locks.sh            (from WINDOWS/ or MAC/; needs Docker + Git Bash/bash)
#
# Then copy docker/ to the other OS folder (full repository:
# tools/sync_os_folders.sh) and commit; the next run_crisprscreen notices
# docker/ changed and rebuilds.
# To move a tool to a new version, point images.txt at the new BioContainers
# tag (the same one the cluster's module uses).
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
tmp="$(mktemp -d)"
cid=""
trap 'rm -rf "${tmp}"; [[ -n "${cid}" ]] && docker rm "${cid}" >/dev/null 2>&1 || true' EXIT
export MSYS_NO_PATHCONV=1   # Git Bash on Windows: don't rewrite /usr/local/...

while read -r tool image; do
    [[ -z "${tool}" || "${tool}" == \#* ]] && continue
    echo "== ${tool}: ${image}"
    docker pull -q --platform linux/amd64 "${image}" >/dev/null
    cid="$(docker create --platform linux/amd64 "${image}")"
    rm -rf "${tmp}/meta"
    # The local side in Windows form under Git Bash (path conversion is off above).
    dest="${tmp}/meta"
    command -v cygpath >/dev/null 2>&1 && dest="$(cygpath -w "${dest}")"
    docker cp "${cid}:/usr/local/conda-meta" "${dest}" >/dev/null
    docker rm "${cid}" >/dev/null; cid=""
    {
        echo "# Pinned package list for the ${tool} environment of the CRISPR_SCREEN Docker image."
        echo "# Same package builds as ${image} (the cluster version's container),"
        echo "# read from its conda-meta/ by docker/env/make_locks.sh."
        echo "@EXPLICIT"
        for f in "${tmp}"/meta/*.json; do
            # Top-level fields of a conda-meta record; nested entries carry no url/md5.
            url="$(grep -o '"url": *"[^"]*"' "${f}" | head -1 | sed 's/.*"url": *"//; s/"$//')"
            md5="$(grep -o '"md5": *"[^"]*"' "${f}" | head -1 | sed 's/.*"md5": *"//; s/"$//')"
            fn="$(grep -o '"fn": *"[^"]*"' "${f}" | head -1 | sed 's/.*"fn": *"//; s/"$//')"
            subdir="$(grep -o '"subdir": *"[^"]*"' "${f}" | head -1 | sed 's/.*"subdir": *"//; s/"$//')"
            # A package the image built from a local channel: use its public copy on bioconda.
            case "${url}" in
                http*) ;;
                *) url="https://conda.anaconda.org/bioconda/${subdir:-noarch}/${fn}" ;;
            esac
            [[ -n "${md5}" ]] && echo "${url}#${md5}" || echo "${url}"
        done | sort
    } > "${here}/${tool}.explicit.txt"
    echo "   $(grep -c '^http' "${here}/${tool}.explicit.txt") packages -> ${tool}.explicit.txt"
done < "${here}/images.txt"
