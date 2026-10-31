# CLAUDE.md

Read `AGENTS.md` in this directory: it explains how to run, monitor and
customize this CUT&RUN project, and what never to touch.

## Ask about the experiment; never assume it

Never assume the biology or how the experiment was done or processed (assay, antibodies, cell
line, what each sample is, replicates, IgG pairing, spike-in, timepoints, the sequencing library,
which comparisons answer the question). Ask the user, in short batches, before preparing inputs,
choosing peak sets or comparisons, interpreting results or writing Methods; record the answers in
the project's `CONTEXT.md`. When the data seem to contradict what the user said, show them the
numbers and ask; never reinterpret the design yourself. Details: "Ask about the experiment; never
assume it" in `AGENTS.md`.

## Delegate to subagents (default, not the exception)

The agent definitions are meant to be in `.claude/agents/` of this project (copied there by `init`
once the kit is in the pipeline; see "Subagents" in `AGENTS.md`). The lead (you) uses the most
expensive model: it spends its tokens on the user, decisions, running `run_cutandrun.sh` and
commits. Every task in the table goes to its agent.

| Task | Agent | Model |
|---|---|---|
| Messy input files (Excel/CSV sample sheets, `samplesheet.csv`, `norm_runs.csv`, `peaksets.csv`, collaborator or vendor tables) | `input-preparer` | sonnet |
| Reading: a question needing >2 files, or a file over ~300 lines | `scout` (several in parallel) | haiku |
| After every run, success or failure | `log-inspector`, then `results-reporter` on success | haiku |
| Comparing with earlier or collaborator numbers | `output-verifier` | haiku |
| Statistics questions and contrast or peak-set design | `statistician` | sonnet |
| Plots | `figure-maker` | sonnet |
| Methods text | `methods-writer` | sonnet |
| Edits to `utils/` scripts or project CSVs beyond a trivial fix, and new analysis scripts | `code-supervisor` (spawns `implementer`) | sonnet |
| Outside facts (tool flags, library formats) | `web-researcher` | haiku |
| Papers and citations | `literature-scout` | sonnet |
| How csaw, edgeR, ChIPseqSpikeInFree or the pipeline code works | `code-researcher` | sonnet |

If these agents are not listed as available subagent types (the project and its `.claude/agents/`
were created mid-session by `init`), either spawn a general-purpose agent whose prompt is the agent
file's body (after the frontmatter) plus the brief, with the model from the file's frontmatter, or
ask the user to restart Claude Code in the project folder. Starting Claude Code inside the project
folder is recommended for every later session.

## Analyses beyond the pipeline

Custom analysis is allowed, as a script `utils/<name>.<R|py|sh>` that writes to a new
`results/<name>/` folder. This pipeline has no `script` subcommand yet, so the script is not run by
`run_cutandrun.sh`; the rules for submitting it are in `AGENTS.md`. Flow: ask the user (`CONTEXT.md`),
`statistician` designs it, `code-supervisor` writes it, the lead runs it (only with the user's OK),
`log-inspector` reads its logs, `results-reporter`/`figure-maker` review, `methods-writer` writes
`results/<name>/METHODS.md`. Outside the pinned images (conda) only with the user's explicit OK.
Details: "Analyses beyond the pipeline (extra steps)" in `AGENTS.md`.

## How to brief

- Each brief stands alone: absolute paths, the project folder, that the pipeline folder is
  read-only, the files the agent owns, the expected return, and the experimental context the
  user confirmed (or the path to `CONTEXT.md`).
- Spot-check numbers from Haiku agents by opening the cited path. Run independent agents in one message.
- If an agent fails twice, escalate (another model or agent) or ask the user.

## The lead keeps

Talking to the user, decisions, running `run_cutandrun.sh`, anything destructive (only with
the user's OK), and git commits.
