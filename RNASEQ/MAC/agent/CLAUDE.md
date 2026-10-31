# CLAUDE.md

Read `AGENTS.md` in this directory: it explains how to run, monitor and customize this RNA-seq project on this Mac (Docker),
and what never to touch.

## Ask about the experiment; never assume it

Never assume the biology or how the experiment was done or processed (the question, RNA source, library prep and kit,
strandedness, sequencing, replicates, treatments, comparisons, expected genes). Ask the user, in short batches, before
preparing the samplesheet, choosing a strandedness, interpreting results or writing Methods; record the answers in the
project's `CONTEXT.md`. When the data seem to contradict what the user said, show them the numbers and ask; never
reinterpret the design yourself. Details: "Ask about the experiment; never assume it" in `AGENTS.md`.

## Delegate to subagents (default, not the exception)

The agent definitions are in `.claude/agents/` of this project (copied there by `init`). The lead (you) uses the most
expensive model: it spends its tokens on the user, decisions, `utils/run_rnaseq.sh` and commits. Every task in the table
goes to its agent.

| Task | Agent | Model |
|---|---|---|
| Messy input files (Excel/CSV sample sheets, `samplesheet.csv`, collaborator or core tables) | `input-preparer` | sonnet |
| Reading: a question needing >2 files, or a file over ~300 lines | `scout` (several in parallel) | haiku |
| After every run, success or failure | `log-inspector`, then `results-reporter` on success | haiku |
| Comparing with earlier, benchmark or collaborator numbers | `output-verifier` | haiku |
| Statistics questions, counting choices and any design beyond counts | `statistician` | sonnet |
| Plots | `figure-maker` | sonnet |
| Analysis beyond the pipeline (`script` command; "Analyses beyond the pipeline" in `AGENTS.md`) | `statistician` (design), `code-supervisor` (script), `log-inspector` (script log), `methods-writer` | sonnet / haiku |
| Methods text | `methods-writer` | sonnet |
| Edits to `utils/` scripts or the samplesheet beyond a trivial fix | `code-supervisor` (spawns `implementer`) | sonnet |
| Writing or editing code that `code-supervisor` hands over | `implementer` | haiku |
| Outside facts (tool flags, library formats) | `web-researcher` | haiku |
| Papers and citations | `literature-scout` | sonnet |
| How HISAT2, featureCounts, RSeQC or the pipeline code works | `code-researcher` | sonnet |

If the agents are not listed as available subagent types (e.g. the project was created mid-session), either spawn a
general-purpose agent whose prompt is the agent file's body (after the frontmatter) plus the brief, with the model from the
file's frontmatter, or ask the user to restart Claude Code in the project folder (recommended for every later session, so
it loads this `CLAUDE.md` and the agents).

Analyses beyond the pipeline are allowed: a script in `utils/`, run with `bash utils/run_rnaseq.sh script utils/<name>.<R|py|sh>`,
writing to `results/<name>/`. Workflow and rules: "Analyses beyond the pipeline" in `AGENTS.md`.

## How to brief

- Each brief stands alone: absolute paths, the project folder, that the pipeline folder is read-only, the files the agent
  owns, the expected return, and the experimental context the user confirmed (or the path to `CONTEXT.md`).
- Spot-check numbers from Haiku agents by opening the cited path. Run independent agents in one message.
- If an agent fails twice, escalate (another model or agent) or ask the user.

## The lead keeps

Talking to the user, decisions, running `utils/run_rnaseq.sh` (run, status, logs, stop, clean, script), anything destructive
(only with the user's OK), the Docker memory setting (the user sets it in Docker Desktop, only with the user's OK), and git
commits.
