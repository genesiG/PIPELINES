# CLAUDE.md

This folder is the cluster (PMACS/LSF) version of the lab-wide CRISPR_SCREEN
pipeline. Read `AGENTS.md` in this directory first, and only change anything
here when the user explicitly asks to change it for everyone.

- **First-time cluster setup** (your account not set up yet): `SETUP.md`, once.
- **Running an analysis**: `README.md` in this folder. A project gets its own
  `AGENTS.md` from `run_crisprscreen.sh init`; follow that one inside the project.
- **Personal computers** (Windows or Mac) use the Docker versions, which are
  separate folders of the full repository (`WINDOWS/`, `MAC/`), not this one.

## Delegate to subagents (default, not the exception)

The agent definitions are in `.claude/agents/` of this folder (the same kit `init`
copies into projects as `agent/.claude/agents/`). The lead (you) uses the most
expensive model: it spends its tokens on the user, decisions, running
run_crisprscreen and commits. Every task in the table goes to its agent.

| Task | Agent | Model |
|---|---|---|
| Messy input files (Excel/CSV, core sheets, vendor libraries) | `input-preparer` | sonnet |
| Setup help: `SETUP.md` steps, what a PMACS/LSF/Apptainer error means | `scout` (local files), `web-researcher` (outside facts) | haiku |
| Reading: a question needing >2 files, or a file over ~300 lines | `scout` (several in parallel) | haiku |
| After every LSF/Nextflow run, success or failure | `log-inspector`, then `results-reporter` on success | haiku |
| Comparing with earlier or collaborator numbers | `output-verifier` | haiku |
| Statistics questions and contrast design | `statistician` | sonnet |
| Plots | `figure-maker` | sonnet |
| Methods text | `methods-writer` | sonnet |
| Changes to this pipeline (`main.nf`, `modules/`, `bin/`, `conf/`) | `code-supervisor` (spawns `implementer`) | sonnet (implementer: haiku) |
| Changes to a project's own files (`utils/` scripts, CSVs) beyond a trivial fix | `code-supervisor` (spawns `implementer`) | sonnet (implementer: haiku) |
| Outside facts (tool flags, library formats) | `web-researcher` | haiku |
| Papers and citations | `literature-scout` | sonnet |
| How MAGeCK or the pipeline code works | `code-researcher` | sonnet |

## Mandatory triggers (do not do these yourself)

- Preparing messy inputs: `input-preparer`.
- After every run, success or failure: `log-inspector`, then `results-reporter`
  on success.
- Any statistics question: `statistician`.
- Any change to the pipeline: `code-supervisor`.
- Any cluster error you do not recognise: `log-inspector` on the run first, then
  `web-researcher` (outside facts) or `code-researcher` (pipeline code).

## How to brief

- Each brief stands alone: absolute paths (e.g. `/project/rephilllab/CRISPR_SCREEN/HPC/...`),
  the project folder if there is one, that this pipeline folder is read-only
  unless the user asked for a change, the files the agent owns, and the expected return.
- Spot-check numbers from Haiku agents by opening the cited path. Run independent
  agents in one message.
- If an agent fails twice, escalate (another model or agent) or ask the user.

## The lead keeps

Talking to the user, decisions, running run_crisprscreen, anything destructive
(only with the user's OK), and git commits.
