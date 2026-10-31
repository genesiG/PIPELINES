#!/bin/bash
# lsf_watchdog.sh — stops a Nextflow run whose LSF job was aborted at dispatch.
#
# Usage:  lsf_watchdog.sh <nextflow_log> <nextflow_pid> [interval_s=300] [grace_s=600]
#
# Why: on PMACS LSF a job can be aborted at dispatch: it runs 0 s, never writes
# .command.begin or .exitcode, and vanishes from `bjobs`. Nextflow 26.04 then
# keeps that task in SUBMITTED forever (nextflow-io/nextflow#1604) and the run
# hangs for hours.
#
# What it does, while <nextflow_pid> is alive, every interval_s:
#   1. Reads every submission line in <nextflow_log>:
#        GridTaskHandler - [LSF] submitted process <NAME> > jobId: <N>; workDir: <DIR>
#   2. For each job whose work dir has neither .command.begin nor .exitcode,
#      asks `bjobs -noheader -o "stat" <N>`:
#        - "is not found", EXIT, or DONE (without .command.begin) is BAD;
#        - PEND, RUN, SSUSP, PSUSP or USUSP is OK;
#        - any other answer, or a bjobs failure that is not "not found", is
#          UNKNOWN and that job is skipped this cycle (its timer is kept).
#      A job that has .command.begin or .exitcode is OK.
#   3. A BAD job starts a per-job timer. Its timer is cleared as soon as the
#      job is OK. If a job stays BAD for more than grace_s, the task is stuck:
#      a block is printed to stderr (with the last 5 lines of <DIR>/.command.log
#      if present), SIGTERM is sent to <nextflow_pid> so Nextflow shuts down
#      cleanly and marks the run failed, the watchdog waits up to 120 s, then
#      sends SIGKILL if Nextflow is still alive, and exits 0.
#
# The watchdog exits 0 as soon as <nextflow_pid> is gone (normal end of run).
# It never dies on a transient error: no `set -e` and no `set -u`; a failed
# bjobs call or a missing log only skips that cycle.
#
# Pure bash + grep/sed/awk. State lives in an associative array (bad_since).
# Started by run_pipeline.sh in the head job; see README.

set -o pipefail

usage() {
    echo "usage: $(basename "$0") <nextflow_log> <nextflow_pid> [interval_s=300] [grace_s=600]" >&2
}

if (( $# < 2 )); then
    usage
    exit 2
fi

LOG="$1"
NF_PID="$2"
INTERVAL="${3:-300}"
GRACE="${4:-600}"

for v in "${NF_PID}" "${INTERVAL}" "${GRACE}"; do
    if [[ ! "${v}" =~ ^[0-9]+$ ]]; then
        echo "RNASEQ WATCHDOG: numeric argument expected, got '${v}'." >&2
        usage
        exit 2
    fi
done

declare -A bad_since   # LSF job id -> epoch second it was first seen BAD

# One "jobId<TAB>name<TAB>workDir" line per submission in the Nextflow log.
submissions() {
    grep -F '[LSF] submitted process' "${LOG}" 2>/dev/null \
        | sed -n 's/.*\[LSF\] submitted process \(.*\) > jobId: \([0-9][0-9]*\); workDir: \(.*[^[:space:]]\)[[:space:]]*$/\2\t\1\t\3/p'
}

# Prints BAD, OK or UNKNOWN for one LSF job that has no start/exit marker.
lsf_state() {
    local jid="$1" out rc status
    out="$(bjobs -noheader -o "stat" "${jid}" 2>&1)"
    rc=$?
    if [[ "${out}" == *"is not found"* ]]; then
        echo BAD
        return
    fi
    if (( rc != 0 )); then
        echo UNKNOWN
        return
    fi
    status="$(printf '%s\n' "${out}" | awk 'NF { print $1; exit }')"
    case "${status}" in
        EXIT|DONE)                    echo BAD ;;
        PEND|RUN|SSUSP|PSUSP|USUSP)   echo OK ;;
        *)                            echo UNKNOWN ;;
    esac
}

# Prints the block, sends SIGTERM to Nextflow, then SIGKILL after 120 s.
report_and_stop() {
    local jid="$1" name="$2" dir="$3" i
    {
        echo "RNASEQ WATCHDOG: LSF job ${jid} for task ${name} ended without starting (work dir ${dir}; LSF log: ${dir}/.command.log). Nextflow would wait for it forever. Stopping the run; relaunch with run_rnaseq.sh — finished steps are reused (-resume)."
        if [[ -f "${dir}/.command.log" ]]; then
            echo "Last 5 lines of ${dir}/.command.log:"
            tail -n 5 "${dir}/.command.log"
        fi
    } >&2
    kill -TERM "${NF_PID}" 2>/dev/null
    for (( i = 0; i < 120; i++ )); do
        kill -0 "${NF_PID}" 2>/dev/null || return 0
        sleep 1
    done
    echo "RNASEQ WATCHDOG: Nextflow (pid ${NF_PID}) still alive 120 s after SIGTERM; sending SIGKILL." >&2
    kill -KILL "${NF_PID}" 2>/dev/null
}

while kill -0 "${NF_PID}" 2>/dev/null; do
    sleep "${INTERVAL}"
    kill -0 "${NF_PID}" 2>/dev/null || break

    while IFS=$'\t' read -r jid name dir; do
        [[ -n "${jid}" ]] || continue
        if [[ -e "${dir}/.command.begin" || -e "${dir}/.exitcode" ]]; then
            unset "bad_since[${jid}]"
            continue
        fi
        case "$(lsf_state "${jid}")" in
            BAD)
                now="$(date +%s)"
                [[ -n "${bad_since[${jid}]:-}" ]] || bad_since[${jid}]="${now}"
                if (( now - ${bad_since[${jid}]} > GRACE )); then
                    report_and_stop "${jid}" "${name}" "${dir}"
                    exit 0
                fi
                ;;
            OK)
                unset "bad_since[${jid}]"
                ;;
            *)
                ;;
        esac
    done < <(submissions)
done

exit 0
