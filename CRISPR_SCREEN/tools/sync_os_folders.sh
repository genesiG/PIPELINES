#!/bin/bash
# sync_os_folders.sh — keep the Docker versions in WINDOWS/ and MAC/ identical
# where they must be (maintainers only).
#
# The two folders differ ONLY in: the launcher (run_crisprscreen.ps1/.bat vs
# run_crisprscreen.sh), README.md, AGENTS.md, CLAUDE.md, agent/ (the project
# AGENTS.md that `init` copies) and examples/config.yaml (path style).
# Everything else — main.nf, nextflow.config, conf/, modules/, bin/, docker/,
# examples/*.csv, examples/project.gitignore, test/make_test_project.py — is shared.
# test/make_test_project.py also exists in HPC/test/ and must stay identical to WINDOWS/test/
# (--check verifies it; HPC/test/ is not synced by this script).
#
#   bash tools/sync_os_folders.sh --check          list shared files and agent kit copies that differ (exit 1 if any)
#   bash tools/sync_os_folders.sh WINDOWS MAC      copy shared files WINDOWS -> MAC
#   bash tools/sync_os_folders.sh MAC WINDOWS      copy shared files MAC -> WINDOWS
#
# Run it after editing shared files in one folder, then commit both. Changes
# that should also reach the cluster version go into HPC/ by hand: main.nf,
# nextflow.config, conf/ and modules/ deliberately differ there (see
# WINDOWS/AGENTS.md "Changing the pipeline"). bin/ must match HPC/bin/ except
# write_methods.py; --check also reports any other difference there.
#
# The agent kit (subagent definitions) lives in <v>/.claude/agents/ and <v>/agent/.claude/agents/
# for each of HPC, WINDOWS and MAC. All six copies must be identical. --check verifies this;
# nothing is synced. HPC/agent/.claude/agents is the canonical copy: after an edit, copy it by hand,
# e.g. `cp -R HPC/agent/.claude/agents/. WINDOWS/agent/.claude/agents/`, and the same to
# WINDOWS/.claude/agents/ and to MAC's two places. There is no root .claude/ any more.
#
# Each version folder (HPC/, WINDOWS/, MAC/) must stay self-contained: labmates clone only one
# of them with git sparse checkout. This script is for the full repository only (maintainers).
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
shared=(main.nf nextflow.config conf modules bin docker examples/contrasts_example.csv
        examples/library_example.csv examples/samplesheet_example.csv examples/project.gitignore
        test/make_test_project.py)

if [[ "${1:-}" == "--check" ]]; then
    rc=0
    for p in "${shared[@]}"; do
        if ! diff -rq --strip-trailing-cr "${root}/WINDOWS/${p}" "${root}/MAC/${p}" >/dev/null 2>&1; then
            diff -rq --strip-trailing-cr "${root}/WINDOWS/${p}" "${root}/MAC/${p}" || true
            rc=1
        fi
    done
    (( rc == 0 )) && echo "WINDOWS/ and MAC/ shared files are identical."
    if ! diff -q --strip-trailing-cr "${root}/HPC/test/make_test_project.py" "${root}/WINDOWS/test/make_test_project.py" >/dev/null 2>&1; then
        echo "HPC/test/make_test_project.py and WINDOWS/test/make_test_project.py differ"; rc=1
    else
        echo "HPC/test/make_test_project.py and WINDOWS/test/make_test_project.py are identical."
    fi
    if ! cmp -s "${root}/HPC/ANALYSIS.md" "${root}/WINDOWS/ANALYSIS.md" || ! cmp -s "${root}/HPC/ANALYSIS.md" "${root}/MAC/ANALYSIS.md"; then
        echo "HPC/ANALYSIS.md, WINDOWS/ANALYSIS.md and MAC/ANALYSIS.md differ"; rc=1
    else
        echo "HPC/ANALYSIS.md, WINDOWS/ANALYSIS.md and MAC/ANALYSIS.md are identical."
    fi
    bin_rc=0
    for n in $( (ls "${root}/HPC/bin"; ls "${root}/WINDOWS/bin") | sort -u ); do
        [[ "${n}" == write_methods.py ]] && continue
        if ! diff -q --strip-trailing-cr "${root}/HPC/bin/${n}" "${root}/WINDOWS/bin/${n}" >/dev/null 2>&1; then
            echo "HPC/bin/${n} and WINDOWS/bin/${n} differ (only write_methods.py may)"; bin_rc=1; rc=1
        fi
    done
    (( bin_rc == 0 )) && echo "HPC/bin/ and the Docker bin/ match (except write_methods.py)."
    kit_rc=0
    for pair in "HPC/.claude/agents:HPC/agent/.claude/agents" \
                "WINDOWS/.claude/agents:WINDOWS/agent/.claude/agents" \
                "MAC/.claude/agents:MAC/agent/.claude/agents" \
                "HPC/agent/.claude/agents:WINDOWS/agent/.claude/agents" \
                "HPC/agent/.claude/agents:MAC/agent/.claude/agents"; do
        a="${root}/${pair%%:*}" b="${root}/${pair#*:}"
        if [[ ! -d "${a}" || ! -d "${b}" ]]; then
            [[ -d "${a}" ]] || echo "missing directory: ${pair%%:*}"
            [[ -d "${b}" ]] || echo "missing directory: ${pair#*:}"
            kit_rc=1
        elif ! diff -rq --strip-trailing-cr "${a}" "${b}"; then
            kit_rc=1
        fi
    done
    if (( kit_rc == 0 )); then
        echo "agent kits (<v>/.claude/agents and <v>/agent/.claude/agents, HPC, WINDOWS, MAC) are identical."
    else
        rc=1
    fi
    exit "${rc}"
fi

from="${1:?usage: sync_os_folders.sh --check | FROM TO}"; to="${2:?usage: sync_os_folders.sh FROM TO}"
for p in "${shared[@]}"; do
    src="${root}/${from}/${p}" dst="${root}/${to}/${p}"
    mkdir -p "$(dirname "${dst}")"
    if [[ -d "${src}" ]]; then rm -rf "${dst}"; cp -R "${src}" "${dst}"; else cp "${src}" "${dst}"; fi
    echo "synced ${p}"
done
