# CLAUDE.md

Read `AGENTS.md` in this directory: it explains how to run, monitor and
customize this CUT&RUN / chromatin profiling project on this Windows computer (Docker), and
what never to touch.

## Ask about the experiment; never assume it

Never assume the biology or how the experiment was done or processed (assay and
protocol, antibodies and the profile each should give, cell line, treatments, replicates,
IgG pairing, spike-in, sequencing, which comparisons answer the question, expected loci).
Ask the user, in short batches, before preparing inputs, choosing peak sets or comparisons,
interpreting results or writing Methods; record the answers in the project's `CONTEXT.md`.
When the data seem to contradict what the user said, show them the numbers and ask; never
reinterpret the design yourself. Details: "Ask about the experiment; never assume it" in
`AGENTS.md`.

## Delegate to subagents (default, not the exception)

The agent definitions are in `.claude\agents\` of this project (copied there by `init`).
The lead (you) uses the most expensive model: it spends its tokens on the user, decisions,
`utils\run_cutandrun.bat` and commits. Every task in the table goes to its agent.

| Task | Agent | Model |
|---|---|---|
| Messy input files (Excel/CSV sample sheets, `samplesheet.csv`, `norm_runs.csv`, `peaksets.csv`, collaborator or vendor tables) | `input-preparer` | sonnet |
| Reading: a question needing >2 files, or a file over ~300 lines | `scout` (several in parallel) | haiku |
| After every run, success or failure | `log-inspector`, then `results-reporter` on success | haiku |
| Comparing with earlier or collaborator numbers | `output-verifier` | haiku |
| Statistics questions and comparison or peak-set design | `statistician` | sonnet |
| Plots | `figure-maker` | sonnet |
| Analysis beyond the pipeline (`script` command; "Analyses beyond the pipeline" in `AGENTS.md`) | `statistician` (design), `code-supervisor` (script), `log-inspector` (script log), `methods-writer` | sonnet / haiku |
| Methods text | `methods-writer` | sonnet |
| Edits to `utils\` scripts or project CSVs beyond a trivial fix | `code-supervisor` (spawns `implementer`) | sonnet |
| Writing or editing code that `code-supervisor` hands over | `implementer` | haiku |
| Outside facts (tool flags, library formats) | `web-researcher` | haiku |
| Papers and citations | `literature-scout` | sonnet |
| How csaw, edgeR, ChIPseqSpikeInFree or the pipeline code works | `code-researcher` | sonnet |

If the agents are not listed as available subagent types (e.g. the project was
created mid-session), either spawn a general-purpose agent whose prompt is the
agent file's body (after the frontmatter) plus the brief, with the model from the
file's frontmatter, or ask the user to restart Claude Code in the project folder
(recommended for every later session, so it loads this `CLAUDE.md` and the agents).

Analyses beyond the pipeline are allowed: a script in `utils\`, run with
`.\utils\run_cutandrun.bat script utils\<name>.<R|py|sh>`, writing to
`results\<name>\`. Workflow and rules: "Analyses beyond the pipeline (extra steps)"
in `AGENTS.md`.

## How to brief

- Each brief stands alone: absolute paths, the project folder, that the pipeline folder is
  read-only, the files the agent owns, the expected return, and the experimental context the
  user confirmed (or the path to `CONTEXT.md`).
- Spot-check numbers from Haiku agents by opening the cited path. Run independent agents in one message.
- If an agent fails twice, escalate (another model or agent) or ask the user.

## The lead keeps

Talking to the user, decisions, running `utils\run_cutandrun.bat` (run, status, logs, stop,
clean, script), anything destructive (only with the user's OK), the Docker memory question
of the launcher (only with the user's OK), and git commits.
