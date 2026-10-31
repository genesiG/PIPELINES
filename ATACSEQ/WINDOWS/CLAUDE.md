# CLAUDE.md

Read `AGENTS.md` in this directory before doing anything: it is the runbook for setting up and running a labmate's
ATAC-seq analysis on a Windows computer with Docker. The pipeline folder itself is not edited during an analysis.

- **First time on this computer** (Docker not set up, or never checked): follow the setup steps in `AGENTS.md`, section A.1
  (and the lab's shared `SETUP.md` in this folder), once.
- **Running an analysis** (almost always): the runbook in `AGENTS.md`, section A.
- **Pipeline code**: not edited during an analysis; if it seems to need a change, tell the user to contact the lab's maintainers.

## Ask about the experiment; never assume it

Never assume the biology or how the experiment was done or processed (assay and kit, cell type, treatments, which libraries
are replicates, which libraries share a TMM group or a peak group, expected accessibility at which loci, whether Tn5 shift is
wanted). Ask the user, in short batches, before preparing inputs, choosing groups, interpreting results or writing Methods;
record the answers in the project's `CONTEXT.md`. When the data seem to contradict what the user said, show them the numbers
and ask; never reinterpret the design yourself. Details: "Ask about the experiment; never assume it" in `AGENTS.md` of the
project folder (copied by `init`).

## Delegate to subagents (default, not the exception)

The agent definitions are in `.claude/agents/` of this folder. The lead (you) uses the most expensive model: it spends its
tokens on the user, decisions, `run_atacseq.bat` and commits. Every task in the table goes to its agent. The kit has 13 agents:
code-researcher, code-supervisor, figure-maker, implementer, input-preparer, literature-scout, log-inspector, methods-writer,
output-verifier, results-reporter, scout, statistician and web-researcher.

| Task | Agent | Model |
|---|---|---|
| Setting up a computer: reading the setup survey, logs and files | `scout` | haiku |
| Docker or WSL error (setup, image pull, run) | `web-researcher` for the outside facts, `log-inspector` for the log | haiku |
| Failed `setup` or image pull | `log-inspector` | haiku |
| Preparing inputs: messy sample sheets, library tables (Excel/CSV, core sheets, collaborator files) | `input-preparer` | sonnet |
| Reading: a question needing >2 files, or a file over ~300 lines | `scout` (several in parallel) | haiku |
| After every run, success or failure | `log-inspector`, then `results-reporter` on success | haiku |
| Comparing with earlier or collaborator numbers | `output-verifier` | haiku |
| Statistics questions and comparison design (TMM groups, MACS3 settings, what counts as a result) | `statistician` | sonnet |
| Plots of a project | `figure-maker` | sonnet |
| A project's own analysis beyond the pipeline (`script` command; see the project's `AGENTS.md`) | `statistician` (design), `code-supervisor` (script), `log-inspector`, `methods-writer` | sonnet / haiku |
| Methods text | `methods-writer` | sonnet |
| Outside facts (tool flags, file formats, Docker or WSL behaviour) | `web-researcher` | haiku |
| Papers and citations | `literature-scout` | sonnet |
| How the pipeline code, TMM, MACS3 or bowtie2 work | `code-researcher` | sonnet |

Mandatory triggers:
- setting up a computer, or any Docker or WSL failure;
- preparing sample or library files from messy sources;
- after every run;
- before reporting statistics or comparing with earlier numbers;
- right after `init` (see below).

**After `init`:** users usually ask you to set the project up before running, so the project's `.claude\agents\` and
`CLAUDE.md`/`AGENTS.md` appear mid-session. Immediately after `init`, read the new project's `CLAUDE.md` and `AGENTS.md` and
follow them from then on (delegation table and triggers included). If Claude Code was started in this folder, its
`.claude\agents\` is the same kit, so the agents are available now. If they are not listed as available subagent types, either
spawn a general-purpose agent whose prompt is the agent file's body (after the frontmatter) plus the brief, with the model from
the file's frontmatter, or ask the user to restart Claude Code in the project folder (recommended for every later session:
start Claude Code inside the project folder).

**Analyses beyond the pipeline** (a plot or test the pipeline doesn't do) are allowed: the script goes in the project's
`utils\` and runs with `.\utils\run_atacseq.bat script utils\<name>.<R|py|sh>`, writing to `results\<name>\`. The workflow and
rules are in the project's `AGENTS.md`.

The delegation rules for project folders are in `agent\CLAUDE.md`, which is copied into projects.

## How to brief

- Each brief stands alone: absolute paths, the project folder (or this pipeline folder, which is read-only for analysis work),
  the files the agent owns, the expected return, and the experimental context the user confirmed (or the path to `CONTEXT.md`).
- Spot-check numbers from Haiku agents by opening the cited path. Run independent agents in one message.
- If an agent fails twice, escalate (another model or agent) or ask the user.

## The lead keeps

Talking to the user, decisions, running `run_atacseq.bat`, anything destructive (only with the user's OK), and git commits.
