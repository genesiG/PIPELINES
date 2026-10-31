# AGENTS.md — CRISPR_SCREEN pipeline, Docker version for **Mac** (Richard Phillips Lab)

This folder is the lab's pooled-CRISPR-screen pipeline (Nextflow: demultiplex →
FastQC/MultiQC → `mageck count` → guide-representation QC → optional `mageck test`),
packaged to run on a labmate's own Mac inside Docker. (It also runs on Linux
with Docker Engine.) The user is usually a **wet-lab scientist**: do the
technical work for them, explain in plain language, and ask them only for what
you can't know (their samples, their files, their comparisons).

First time on this computer (Docker not set up yet, or never checked): follow
`SETUP.md` in this folder once, then come back here.

Almost always the task is to set up, run or check an analysis. Follow the runbook
below. You never edit anything in this pipeline folder for that: pipeline code is
not edited during an analysis. If the pipeline itself seems to need a change, ask
the maintainers.

If the current folder is a **project** folder (it has `config.yaml`,
`samplesheet.csv`, `results/`), read that project's own `AGENTS.md` instead.

## How to call the launcher

Everything goes through one command, `run_crisprscreen.sh`, in this folder.
Call it with `bash` and its full path (works even if the file lost its
executable bit, and doesn't depend on the user's PATH):

```bash
bash "/Users/<user>/CRISPR_SCREEN/MAC/run_crisprscreen.sh" status
```

| Command | What it does |
|---|---|
| `run_crisprscreen.sh setup --yes` | one-time per Mac: Docker Desktop check/install, Rosetta 2 (Apple chip), build the Docker image (10–30 min the first time), add the command to PATH (`~/.zshrc`) |
| `run_crisprscreen.sh init "<project folder>"` | create a project folder: subfolders, `config.yaml`, `examples/`, `AGENTS.md`, launchers, git |
| `run_crisprscreen.sh run --detach` (in the project folder) | check the FASTQ paths, start the pipeline in Docker, return immediately |
| `run_crisprscreen.sh status` (in the project folder) | running / succeeded / failed, per-step counts, failed step, last log lines |
| `run_crisprscreen.sh logs` | follow the live log (blocks until the run ends — avoid in an agent; use `status`) |
| `run_crisprscreen.sh stop` | stop a run (finished steps are kept; the next run resumes) |
| `run_crisprscreen.sh clean --yes` | delete `metadata/work/` and `metadata/.nextflow/` (resume cache; plus a legacy root `work/` and `.nextflow/` if present) now, to free disk space — the next run starts from scratch, results are kept; without `--yes` it asks first; ask the user first |
| `run_crisprscreen.sh script utils/<name>.<R\|py\|sh> [args...]` (in the project folder) | run the project's own analysis script in the pipeline's software environment (Docker image, same limits as `run`), in the foreground; log in `metadata/script_logs/` — see the project's `AGENTS.md`, "Analyses beyond the pipeline" |

Every project folder also gets its own launcher, `utils/run_crisprscreen.sh`
(same commands, always for that folder; agents run `./utils/run_crisprscreen.sh`
from the project folder), and `RUN_ANALYSIS.command` (double-click for humans:
it waits for Enter at the end, so **agents must not call
`RUN_ANALYSIS.command`**).

Long commands: `setup` (image build) and a run can take a long time. Run
`setup` in the background or with a long timeout, and always start runs with
`--detach`, then poll `status` every few minutes.

## A. Runbook: from "here are my FASTQs" to results

### 1. One-time computer setup

First time on this Mac? Follow `SETUP.md` in this folder (shared by every lab
Docker pipeline: survey the hardware, Docker Desktop + Rosetta, memory, file
sharing, known hiccups and fixes). The steps below are this pipeline's part of it.

1. Check: `docker --version` and `docker info` (Docker Desktop's CLI may live
   in `/Applications/Docker.app/Contents/Resources/bin` or `~/.docker/bin`).
   If both work, go to step 3.
2. `bash ".../MAC/run_crisprscreen.sh" setup --yes`. What can happen:
   - Docker Desktop not installed → with Homebrew it runs
     `brew install --cask docker` (may ask for the user's password — tell
     them); without Homebrew it stops and prints the download link. Then the
     user must **open Docker from Applications once** and accept its terms (no
     account needed — "Skip"). Explain this plainly, wait for them, and re-run
     `setup`.
   - Docker installed but not running → it runs `open -a Docker` and waits.
   - Apple chip (M1/M2/M3/M4): the image is `linux/amd64` (identical on every
     computer), run through Rosetta 2; `setup` installs Rosetta if missing
     (`softwareupdate --install-rosetta`). In Docker Desktop → Settings →
     General, "Use Rosetta for x86_64/amd64 emulation on Apple Silicon" should
     be on (the default); without it everything is several times slower.
3. `setup` builds the image (`rephilllab/crisprscreen:0.1.0`) and prints the
   tool versions. Re-running `setup` is harmless (it rebuilds only if
   `docker/` changed). A run also builds the image if it's missing.
4. Resources: every run gets one budget, enforced by Docker as a hard limit —
   what the pipeline can use (8 GB) but never more than 75% of the RAM, always
   leaving 3 GB to macOS; all CPU cores but 1–2; never more than Docker
   Desktop may use (README.md, "Resources"). `setup`
   reports it and, if Docker Desktop's own memory setting is what holds it
   back, says what to set in Docker Desktop → Settings → Resources → Memory
   (the user does that; optional — this pipeline needs only ~2.5 GB).

### 2. Gather what you need from the user (ask; never guess)

- **A name for the screen** → `project_name` (letters, digits, `_`; used as a prefix
  on every output file) and the project folder's name.
- **Where the FASTQ files are**, and which kind:
  - **multiplexed**: one big file (often `Undetermined_*_R1_001.fastq.gz`)
    holding every sample, split by an inline barcode at the start of the read
    → `fastq_r1` (+ `fastq_r2` and `paired_end: true` if paired) and the
    samplesheet needs `barcode` values;
  - **already demultiplexed**: one file per sample → `fastq_dir`.
- **The samples**: names, and for multiplexed data each sample's barcode
  sequence (as it appears at the start of R1). Sequencing-core sheets often
  hold these; ask for the sheet.
- **The sgRNA library file(s)** (gene, guide ID, guide sequence), and — if
  more than one pool/library — which library each sample was screened with.
- **The comparisons** (optional, for `mageck test`): which samples are
  treatment and which are control. If they don't know yet, run without
  `contrasts` first (counting + QC only) and add it later; re-running then
  only does the new part.

If they give Excel files: ask them to "Save As → CSV UTF-8" (Excel or
Numbers → Export → CSV). Check the result is plain comma-separated text with
normal line endings.

**Then ask about the experiment itself (ask; never assume).** Only the user
knows how it was done; ask even when it seems obvious:

- **The question:** what the screen selects for (proliferation/essentiality,
  drug resistance or sensitivity, a FACS sort on a reporter, survival, ...),
  what a hit would look like (depleted, enriched, or both), and which genes
  they expect to behave how (positive controls, known hits).
- **The screen:** knockout, CRISPRi, CRISPRa, base editing or other; the cell
  line or model; the library (name, version, guides per gene, non-targeting or
  safe-harbour controls); MOI and coverage (cells per guide) at transduction and
  at every passage or sort; the selection (antibiotic, days).
- **Each sample:** what it is (plasmid, pre-selection/day 0, a timepoint, a
  treatment and dose, a sort gate, an unsorted input), which samples are
  biological replicates and which are technical, which samples pair with which,
  and anything that went wrong (low cell numbers, a failed PCR, a contamination,
  a re-sequenced sample).
- **The sequencing library:** how the PCR was done (primers, staggers, inline
  barcodes), single- or paired-end, read length, which run(s) and lanes, whether
  the run also holds other people's samples, and whether the FASTQs are already
  demultiplexed.
- **The comparisons:** which treatment-vs-control comparisons answer their
  question, and in which direction a hit counts.

- Ask in short batches (3-6 questions at a time), in plain language, and say why
  each answer matters. Keep asking until the answers are specific; "the usual"
  or "like last time" is not an answer, so ask what it was.
- **Never infer the design from names** (`D0`, `T14`, `ctrl`, `plasmid`, `R1`,
  `rep2`, `hi`/`lo`) or from what is "typical" for CRISPR screens. Ask what each
  name means.
- Use only what the user confirmed. Where an answer is missing, say what you
  would do instead and wait for an explicit OK; never fill the gap silently,
  and never pick contrasts, controls or replicates yourself.
- **Write the answers down** in `CONTEXT.md` in the project folder (dated, in the
  user's words, with open questions marked as open). Read it at the start of
  every session; ask the user to confirm it is still right before you rely on
  it, and add to it as you learn more. It is tracked by git with the inputs.
- **When the data seem to contradict the context, flag it; don't reinterpret.**
  Examples: a "plasmid" or "day 0" sample with a skewed distribution (high Gini,
  many zero-count guides); a replicate that correlates poorly with its partner;
  essential genes not depleted in a sample described as a late timepoint; a
  sample that maps mostly to another library; a barcode that is never found.
  Show the user the numbers and the file they come from, say what the data
  would be consistent with, and ask. The user's description stands until they
  change it; don't rename samples, drop them or change contrasts on your own.
- Subagents cannot talk to the user: put the confirmed context (or the path to
  `CONTEXT.md`) in every brief, and pass their "Questions for the user" back to
  the user rather than answering them yourself.


### 3. Create the project folder

Put projects **outside** the pipeline's git clone, e.g.
`~/Documents/CRISPR_screens/<screen_name>` (a disk with enough free space:
plan about 2× the size of the raw FASTQ files, plus results). Avoid iCloud-synced folders
(Desktop/Documents with "Optimize Mac Storage" can evict files) for large data
— ask the user if unsure. Use a new folder: if the one you picked already
exists and holds anything, list it and ask the user before using it (`init`
never overwrites files, but it adds its own files and a git repository to
whatever is there). Then:

```bash
bash ".../MAC/run_crisprscreen.sh" init ~/Documents/CRISPR_screens/<screen_name>
```

**Immediately after `init`, read the new project's `CLAUDE.md` and `AGENTS.md` and
follow them from then on** (delegation table and triggers included). Users usually
ask you to set up the project before running it, so the project's `.claude/agents`
and `CLAUDE.md`/`AGENTS.md` are created in the middle of your session. If Claude Code
was started in the pipeline folder, its `.claude/agents` is the same kit, so the agents
are available now. If they are not listed as available subagent types, either spawn a
general-purpose agent whose prompt is the agent file's body (after the frontmatter) plus
the brief, with the model from the file's frontmatter, or ask the user to restart Claude
Code in the project folder (recommended for every later session: start Claude Code inside
the project folder, so it loads that project's `CLAUDE.md` and agents).

### 4. Write the input tables (in the project folder)

**`samplesheet.csv`** — one row per sample:
```
barcode_id,barcode,sample_name
2F1,GCCATCTT,baseline_rep1
2F2,GCTGGATT,treatment_low_rep1
```
- Header names exactly as shown. `barcode_id` is a free label; `barcode` is
  the inline barcode (A/C/G/T, unique per sample — duplicates are rejected).
  Both are only needed for **multiplexed** input; with `fastq_dir`, only
  `sample_name` (+ optional `library`) is needed.
- `sample_name`: unique; use only letters, digits and `_` (it becomes file
  names, count-table columns and contrast patterns). Keep replicate suffixes
  consistent (`_rep1`, `_rep2`).
- With `fastq_dir`, each sample's file must be named `<sample_name>.fastq.gz`
  (or `.fastq`; `<sample_name>_R1.fastq.gz`/`_R2` when paired). If the core's
  files are named differently (`S1_S1_L001_R1_001.fastq.gz`), either use those
  stems as sample names or — **with the user's OK** — rename the files.
- Several libraries: add a `library` column whose values are the library
  **file names without `.csv`**, and set `library: libraries/*.csv` in
  `config.yaml` (library files in `libraries/`). Without that column, every
  sample is counted against every library (useful to find out which pool a
  sample holds; `mageck test` then can't run).

**`library.csv`** — one row per guide: `gene,gRNA_ID,sequence`. `sequence` is
the guide (spacer) only, no PAM, normally 20 nt (`sgrna_len: 20`). Non-targeting
controls: give them all the **same** `gene` value containing `NTC`, `Non-targeting`
or `Control` (e.g. `Control` for `Control-1.1`, `Control-2.2`, ...), not one per
guide: MAGeCK treats every distinct `gene` value as a gene, and the QC then
treats them as controls.

Converting a vendor or collaborator file (check each one: they are rarely
clean):
- `gene` missing: take it from the guide ID (`Aicda.6` -> `Aicda`); keep
  `gRNA_ID` exactly as given. IDs only numbered per gene: number in file order.
- One file holding several pools (a `Pool` column): write one
  `libraries/<pool>.csv` per pool.
- Excel: a header row may sit at the bottom, or be missing; a one-column
  sheet may hold a title (e.g. `Intergenic sequence`) above bare sequences —
  use the title as `gene`. Drop anything that isn't a guide row.
- Afterwards, check the guide count per file against the source, and that
  every sequence is A/C/G/T of the expected length (report any 21-mers etc.
  to the user rather than editing them).

**`contrasts.csv`** (optional) — `treatment,control`, one row per comparison.
Each value is a **case-insensitive regular-expression search** against the
count-table sample names: `treatment_high` matches `treatment_high_rep1` and
`_rep2` (replicates are pooled into one comparison); `Pool_1` also matches
`Pool_10` — anchor with `^Pool_1$`. Then set `contrasts: contrasts.csv`.

Start from `examples/*_example.csv`; save as UTF-8.

### 5. Edit `config.yaml`

- `project_name`
- `fastq_r1` (and `fastq_r2`) **or** `fastq_dir` — Mac paths anywhere, including
  external drives (`/Volumes/MyDrive/run7/Undetermined_S0_L001_R1_001.fastq.gz`)
  and `~/...`. The FASTQ(s) stay where they are: the launcher checks that they
  exist and reads them from that path; it never moves, copies, renames or
  deletes them, and never rewrites these keys. Tell the user where the file is
  read from and that it stays there. Keep the path stable (if the files are
  moved, update `config.yaml`). Docker Desktop shares only `/Users`, `/Volumes`,
  `/private`, `/tmp` and `/var/folders` by default: for a FASTQ elsewhere, add
  its folder in Docker Desktop > Settings > Resources > File sharing. No comma in
  the path. A FASTQ outside the project folder is mounted read-only into Docker
  (its folder, at `/inputs/fastq_r1`, `/inputs/fastq_r2` or `/inputs/fastq_dir`).
  A FASTQ that is a symlink must be given by its real path. Several projects may
  point at the same raw file; it is not duplicated.
- `paired_end`, `contrasts`, `library` (glob for several pools)
- Optional `max_cpus` / `max_memory_gb` if the user needs the Mac to stay
  responsive.
- Every other path (`samplesheet`, `library`, `contrasts`, script overrides)
  must be **inside the project folder**; the launcher refuses paths outside it.
- Leave `trim5` unset (automatic guide positioning) unless the user has a
  reason.

### 6. Check, then run

Check before the first run (the pipeline also checks, but catching it now is
faster): headers exactly right; sample names unique and matching the FASTQ
names (fastq_dir mode); barcodes unique; every `library` value matches a
library file name; contrast patterns match the intended samples.

```bash
cd "<project folder>"
bash utils/run_crisprscreen.sh run --detach
bash utils/run_crisprscreen.sh status        # every few minutes
```

Rough timing: image build 10–30 min (first time only); demultiplexing a full
lane 20–90 min (longer on Apple chips, which emulate amd64); the rest usually
under an hour. The Mac must stay on, plugged in, lid open: the launcher
prevents idle sleep while the analysis runs (`caffeinate`), but closing the lid
or a dead battery pauses it. If Docker restarts, just run again (it resumes).

Nextflow's working files are in `metadata/work/` and its resume cache in
`metadata/.nextflow/` (`metadata/` is git-ignored). They are kept by default
(`keep_work: true` in `config.yaml`), so a re-run resumes and only redoes what
changed (e.g. after adding contrasts, or a plot tweak). With `keep_work: false`
they are deleted after a successful run, and re-running recomputes everything
(about 10–20 min for a full lane); failed or stopped runs always keep them.
By default (`publish_mode: auto`) files in `results/` and `original_data/` are
hard links to the files in `metadata/work/` when the disk allows, so big FASTQs
are not stored twice: never edit a published file in place. A root `work/` or
`.nextflow/` from an older version is legacy and unused (`run`/`status` print a
note; `clean` removes it). To free space use `clean` (ask first). Intermediate
files are in `results/intermediate_files/`.

**Analyses beyond the pipeline.** If the user wants an analysis the pipeline doesn't
do, it is allowed: the script goes in the project's `utils/<name>.<R|py|sh>`, writes
only to a new `results/<name>/`, and is run with `./utils/run_crisprscreen.sh script
utils/<name>.R` (pinned image, same memory limits, log in `metadata/script_logs/`).
The workflow (context from the user, `statistician`, `code-supervisor`, run, `log-inspector`,
reviewers, `methods-writer`), the image's software and the fallback of running on the host
(only after asking the user) are in the project's `AGENTS.md`, "Analyses beyond the
pipeline". Plots: replace the right script in the project's `utils/`, re-run (with
`keep_work` true only the plot step and the report re-run), and review the PNG twins.

### 7. Report results

Status `succeeded` → point the user to (`open <file>` opens it in the browser):
- `results/counts/qc/report/<project>_<library>_qc_report.html` — the QC
  report, one per library; written for non-specialists (start here);
- `results/demultiplexing/qc/multiqc/multiqc_report.html` — read QC;
- `results/counts/*.count.txt` — the counts;
- `results/mageck_test/<treatment>_vs_<control>/*.gene_summary.txt` — hits;
- `METHODS.md` — Materials and Methods of this run, ready to adapt for a paper.

Summarize from `results/counts/qc/tables/*_qc_summary.csv` (`mapping_pct`,
`library_skew_flag`, `undetected_flag`). In an all-vs-all run (more than one library CSV, no
`library` column), a low `mapping_pct` for one library is expected for samples of other
libraries; in `*_library_breakdown_wide.csv` Low mapping means `pct_any_library` < 30, and a
"Mainly <library>" badge is neutral. Otherwise Low mapping means `mapping_pct` < 30, usually a
sample from another library, or a read-structure issue: the read-classification plot tells which (mostly near-miss → read structure /
quality; mostly far-miss → wrong library).

### 8. When a run fails

`run_crisprscreen.sh status` says which step failed. Its own files are copied
to `metadata/failed_tasks/<timestamp>/<step>/` (`command.sh` = what ran;
`command.err` / `command.log` = what it said). Startup errors (bad samplesheet,
missing file, unknown library name, typo in a key) are printed at the end of
`metadata/nf_logs/nextflow_<timestamp>.log`. The project's `AGENTS.md`
("When a run fails") lists the common cases. Fix the input/config, run again:
finished steps are reused. If it looks like a bug in the pipeline itself, tell
the user and stop — don't patch this folder.

## Subagents

The agent definitions are in `.claude/agents/` of this folder, so an agent working here has them.
The same kit is in `agent/.claude/agents/`, which `init` copies into every project. The delegation
table and the mandatory triggers are in `CLAUDE.md` of this folder (and `agent/CLAUDE.md` for
projects). Main ones: `input-preparer` for the user's messy sample and library sheets (A.4), `scout`
for reading, `log-inspector` after runs (A.8), `results-reporter` for A.7, `code-supervisor` for
pipeline changes (only when the maintainers ask), `statistician` for statistics.

## Hard rules

1. **Don't edit this pipeline folder** for a user's analysis. Per-project
   changes go in the project (`config.yaml`, `utils/` script overrides — see
   the project's `AGENTS.md`).
2. **Start the pipeline only with `run_crisprscreen.sh run`; run a project's own
   analysis scripts only with `run_crisprscreen.sh script`.** Don't run `docker run`,
   `nextflow`, `mageck`, `cutadapt`, or the `bin/` scripts yourself.
3. **The raw FASTQs are the user's data.** They stay where the user keeps them
   (`fastq_r1`/`fastq_r2`/`fastq_dir` point to them). Never move, delete,
   modify or overwrite them, and don't change those keys unless the user asks;
   if the files are moved, update `config.yaml`.
4. **Ask before** deleting anything (`results/`, `metadata/work/` and
   `metadata/.nextflow/` via `clean`), installing software, changing Docker Desktop's settings, or
   renaming the user's files.
5. **Don't hand-edit `results/`**; it is regenerated. Intermediate files are
   inside it, in `results/intermediate_files/`.

6. **Never assume the experiment.** Ask the user (step 2), record it in the
   project's `CONTEXT.md`, and flag data that contradict it instead of
   reinterpreting the design.
