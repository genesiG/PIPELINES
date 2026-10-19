# CLAUDE.md

Read `AGENTS.md` in this directory: it explains how to run, monitor and
customize this CRISPR screen project on this Mac (Docker), and what never to
touch.

## Delegate to subagents (default, not the exception)

The agent definitions are in `.claude/agents/` of this project (copied there by `init`).
The lead (you) uses the most expensive model: it spends its tokens on the user, decisions,
run_crisprscreen (.sh on Mac/Linux/HPC, .bat/.ps1 on Windows) and commits. Every task in
the table goes to its agent.

| Task | Agent | Model |
|---|---|---|
| Messy input files (Excel/CSV, core sheets, vendor libraries) | `input-preparer` | sonnet |
| Reading: a question needing >2 files, or a file over ~300 lines | `scout` (several in parallel) | haiku |
| After every run, success or failure | `log-inspector`, then `results-reporter` on success | haiku |
| Comparing with earlier or collaborator numbers | `output-verifier` | haiku |
| Statistics questions and contrast design | `statistician` | sonnet |
| Plots | `figure-maker` | sonnet |
| Methods text | `methods-writer` | sonnet |
| Edits to `utils/` scripts or project CSVs beyond a trivial fix | `code-supervisor` (spawns `implementer`) | sonnet |
| Outside facts (tool flags, library formats) | `web-researcher` | haiku |
| Papers and citations | `literature-scout` | sonnet |
| How MAGeCK or the pipeline code works | `code-researcher` | sonnet |

## How to brief

- Each brief stands alone: absolute paths, the project folder, that the pipeline folder is
  read-only, the files the agent owns, the expected return.
- Spot-check numbers from Haiku agents by opening the cited path. Run independent agents in one message.
- If an agent fails twice, escalate (another model or agent) or ask the user.

## The lead keeps

Talking to the user, decisions, running run_crisprscreen, anything destructive (only with
the user's OK), and git commits.
