# CLAUDE.md

This folder is the cluster (PMACS/LSF) version of the lab-wide CRISPR_SCREEN
pipeline. Read `AGENTS.md` in this directory first, and only change anything
here when the user explicitly asks to change it for everyone.

- **First-time cluster setup** (your account not set up yet): `SETUP.md`, once.
- **Running an analysis**: `README.md` in this folder. A project gets its own
  `AGENTS.md` from `run_crisprscreen.sh init`; follow that one inside the project.
- **Personal computers** (Windows or Mac) use the Docker versions, which are
  separate folders of the full repository (`WINDOWS/`, `MAC/`), not this one.

## After `init`: switch to the project's instructions

Users usually ask you to set up a project before running it, so the project's `.claude/agents/`,
`CLAUDE.md` and `AGENTS.md` are created mid-session by `run_crisprscreen.sh init`. Immediately
after `init`, read the new project's `CLAUDE.md` and `AGENTS.md` and follow them from then on
(delegation table and triggers included). If Claude Code was started in this pipeline folder,
its `.claude/agents` is the same kit, so the agents are available now. If they are not listed
as available subagent types, either spawn a general-purpose agent whose prompt is the agent
file's body (after the frontmatter) plus the brief, with the model from the file's
frontmatter, or ask the user to restart Claude Code in the project folder (recommended for
every later session, so it loads that project's `CLAUDE.md` and agents).

## Analyses beyond the pipeline

Custom analysis for a project is allowed: a script `utils/<name>.<R|py|sh>` in the project that
writes to a new `results/<name>/`, run with `run_crisprscreen.sh script utils/<name>.R` (an LSF
job; new, not yet tested on the cluster). The project's `AGENTS.md` has the workflow and the
rules. The pipeline's own steps are still only run with the launcher's `run`.

## Ask about the experiment; never assume it

Never assume the biology or how the experiment was done or processed (screen
type, what each sample is, replicates, timepoints, controls, expected hits, how
the sequencing library was made). Ask the user, in short batches, before
preparing inputs, choosing contrasts, interpreting results or writing Methods;
record the answers in the project's `CONTEXT.md`. When the data seem to
contradict what the user said, show them the numbers and ask; never
reinterpret the design yourself. Details: "Ask about the experiment; never
assume it" in `AGENTS.md` of the project folder (copied by `init`).

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
- Any statistics question, and the design of any analysis beyond the pipeline: `statistician`.
- Any change to the pipeline: `code-supervisor`.
- Any cluster error you do not recognise: `log-inspector` on the run first, then
  `web-researcher` (outside facts) or `code-researcher` (pipeline code).

## How to brief

- Each brief stands alone: absolute paths (e.g. `/project/rephilllab/CRISPR_SCREEN/HPC/...`),
  the project folder if there is one, that this pipeline folder is read-only
  unless the user asked for a change, the files the agent owns, the expected return, and the experimental context the
  user confirmed (or the path to `CONTEXT.md`).
- Spot-check numbers from Haiku agents by opening the cited path. Run independent
  agents in one message.
- If an agent fails twice, escalate (another model or agent) or ask the user.

## The lead keeps

Talking to the user, decisions, running run_crisprscreen, anything destructive
(only with the user's OK), and git commits.
