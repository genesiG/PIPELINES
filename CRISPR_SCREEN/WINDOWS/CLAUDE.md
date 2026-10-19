# CLAUDE.md

Read `AGENTS.md` in this directory before doing anything: it is the runbook
for setting up and running a labmate's CRISPR screen analysis on a Windows
computer with Docker, and it says when (rarely) this pipeline folder itself
may be changed.

- **First time on this computer** (Docker not set up, or never checked): follow
  `SETUP.md` in this directory, once.
- **Running an analysis** (almost always): the runbook in `AGENTS.md`, section A.
- **Changing the pipeline**: only when the user explicitly asks for it for the whole
  lab; see `AGENTS.md`, section B.

## Delegate to subagents (default, not the exception)

The agent definitions are in `.claude/agents/` of this folder. The lead (you) uses the
most expensive model: it spends its tokens on the user, decisions, `run_crisprscreen.bat`
and commits. Every task in the table goes to its agent.

| Task | Agent | Model |
|---|---|---|
| Setting up a computer: reading `SETUP.md`'s survey, logs and files | `scout` | haiku |
| Docker or WSL error (setup, image build, run) | `web-researcher` for the outside facts, `log-inspector` for the log | haiku |
| Failed `setup` or image build | `log-inspector` | haiku |
| Preparing inputs: messy sample, library or contrast files (Excel/CSV, core sheets, vendor libraries) | `input-preparer` | sonnet |
| Reading: a question needing >2 files, or a file over ~300 lines | `scout` (several in parallel) | haiku |
| After every run, success or failure | `log-inspector`, then `results-reporter` on success | haiku |
| Comparing with earlier or collaborator numbers | `output-verifier` | haiku |
| Statistics questions and contrast design | `statistician` | sonnet |
| Plots of a project | `figure-maker` | sonnet |
| Methods text | `methods-writer` | sonnet |
| Edits to this pipeline (section B only) | `code-supervisor` (spawns `implementer`) | sonnet |
| Outside facts (tool flags, library formats) | `web-researcher` | haiku |
| Papers and citations | `literature-scout` | sonnet |
| How MAGeCK or the pipeline code works | `code-researcher` | sonnet |

Mandatory triggers:
- setting up a computer, or any Docker or WSL failure;
- preparing sample, library or contrast files from messy sources;
- after every run;
- before reporting statistics or comparing with earlier numbers;
- any change to the pipeline.

The delegation rules for project folders are in `agent\CLAUDE.md`, which is copied into projects.

## How to brief

- Each brief stands alone: absolute paths, the project folder (or this pipeline folder,
  which is read-only for analysis work), the files the agent owns, the expected return.
- Spot-check numbers from Haiku agents by opening the cited path. Run independent agents in one message.
- If an agent fails twice, escalate (another model or agent) or ask the user.

## The lead keeps

Talking to the user, decisions, running run_crisprscreen.bat, anything destructive (only with
the user's OK), and git commits.
