#!/bin/bash
# setup_environment.sh — one-time PATH setup for the Rephill lab's shared
# pipelines. Adds a single line to ~/.bashrc that sources
# /project/rephilllab/.rephilllab.profile (which puts every lab pipeline's
# launcher, e.g. run_crisprscreen.sh, on PATH). Safe to re-run: it's a no-op
# if the line is already there.
#
# Usage (run once, from anywhere):
#   bash /project/rephilllab/setup_environment.sh
#
# This script only ever appends one guarded line to your ~/.bashrc. It never
# overwrites or removes anything already there.
set -euo pipefail

PROFILE="/project/rephilllab/.rephilllab.profile"
LINE="[ -f ${PROFILE} ] && . ${PROFILE}"
BASHRC="${HOME}/.bashrc"

if [[ ! -f "${PROFILE}" ]]; then
    echo "ERROR: ${PROFILE} not found. Nothing to set up — contact the lab." >&2
    exit 1
fi

touch "${BASHRC}"

if grep -qxF "${LINE}" "${BASHRC}"; then
    echo "Already set up: ${BASHRC} already sources ${PROFILE}."
else
    printf '\n# Rephill lab shared shell profile (pipeline launchers, etc.)\n%s\n' "${LINE}" >> "${BASHRC}"
    echo "Added to ${BASHRC}:"
    echo "  ${LINE}"
fi

echo ""
echo "To use it in THIS terminal right now, run:"
echo "  source ${BASHRC}"
echo "Every new terminal from now on will pick it up automatically."
echo ""
echo "Verify with:"
echo "  tail -3 ${BASHRC}"
echo "  which run_crisprscreen.sh"
