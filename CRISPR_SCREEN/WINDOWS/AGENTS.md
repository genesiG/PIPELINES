# AGENTS.md — CRISPR_SCREEN pipeline, Docker version for **Windows** (Richard Phillips Lab)

This folder is the lab's pooled-CRISPR-screen pipeline (Nextflow: demultiplex →
FastQC/MultiQC → `mageck count` → guide-representation QC → optional `mageck test`),
packaged to run on a labmate's own Windows computer inside Docker. The user is
usually a **wet-lab scientist**: do the technical work for them, explain in plain
language, and ask them only for what you can't know (their samples, their files,
their comparisons).

First time on this computer (Docker not set up, or never checked)? Follow
`SETUP.md` in this folder once; it is the one-time computer setup shared by all
the lab's Docker pipelines.

There are two kinds of tasks. Find out which one you have:

- **A. Set up / run / check an analysis** (almost always). Follow the runbook
  below. You never edit anything in this pipeline folder for that.
- **B. Change the pipeline itself** — only when the user explicitly asks to
  change it *for the whole lab*. See "Changing the pipeline" at the end.

If the current folder is a **project** folder (it has `config.yaml`,
`samplesheet.csv`, `results/`), read that project's own `AGENTS.md` instead.

## How to call the launcher

Everything goes through one command, `run_crisprscreen.bat`, in this folder
(full path: `<this folder>\run_crisprscreen.bat`). From the agent's shells:

```powershell
# PowerShell (preferred)
& "D:\path\to\CRISPR_SCREEN\WINDOWS\run_crisprscreen.bat" status
```
```bash
# Git Bash
cmd.exe //c "D:\\path\\to\\CRISPR_SCREEN\\WINDOWS\\run_crisprscreen.bat" status
```

| Command | What it does |
|---|---|
| `run_crisprscreen.bat setup --yes` | one-time per computer: Docker Desktop check/install, build the Docker image (10–30 min the first time), add the command to the user's PATH |
| `run_crisprscreen.bat init "<project folder>"` | create a project folder: subfolders, `config.yaml`, `examples\`, `AGENTS.md`, launchers, git |
| `run_crisprscreen.bat run --detach` (in the project folder) | check the FASTQ paths, start the pipeline in Docker, return immediately |
| `run_crisprscreen.bat status` (in the project folder) | running / succeeded / failed, per-step counts, failed step, last log lines |
| `run_crisprscreen.bat logs` | follow the live log (blocks until the run ends — avoid in an agent; use `status`) |
| `run_crisprscreen.bat stop` | stop a run (finished steps are kept; the next run resumes) |
| `run_crisprscreen.bat clean --yes` | delete `work\` (resume cache) now — rarely needed: it is deleted after a successful run unless `keep_work: true`; needed after failed runs or with `keep_work: true` — ask the user first |

Every project folder also gets its own launcher, `utils\run_crisprscreen.bat`
(same commands, always for that folder; agents run `.\utils\run_crisprscreen.bat`
from the project folder), and `RUN_ANALYSIS.bat` (double-click for humans: it
pauses at the end, so **agents must not call `RUN_ANALYSIS.bat`**).

Long commands: `setup` (image build) and a run can take a long time. Run
`setup` in the background or with a long timeout, and always start runs with
`--detach`, then poll `status` every few minutes.

## A. Runbook: from "here are my FASTQs" to results

### 1. One-time computer setup

First time on this computer? Follow `SETUP.md` in this folder (shared by every
lab Docker pipeline: survey the hardware, Docker Desktop + WSL 2, where
Docker's data and the page file live, memory hogs, known hiccups and fixes).
The steps below are this pipeline's part of it.

1. Check: `docker --version` and `docker info`. If both work, go to step 3.
2. `run_crisprscreen.bat setup --yes`. What can happen:
   - Docker Desktop not installed → it installs WSL 2 + Docker Desktop with
     `winget`. **Windows shows an administrator (UAC) prompt: tell the user to
     approve it.** Then the user must **restart the computer**, open Docker
     Desktop once from the Start menu and accept its terms (no account needed
     — "Skip"). Explain this plainly, then wait for them and re-run `setup`.
   - Docker installed but not running → it starts Docker Desktop and waits.
   - "virtualization is disabled" / WSL errors → hardware virtualization must
     be enabled in the BIOS/UEFI; that's for the user or their IT department.
3. `setup` builds the image (`rephilllab/crisprscreen:0.1.0`) and prints the
   tool versions. Re-running `setup` is harmless (it rebuilds only if
   `docker\` changed). A run also builds the image if it's missing, so a
   forgotten `setup` only costs time.
4. Resources: every run gets one budget, enforced by Docker as a hard limit —
   what the pipeline can use (8 GB) but never more than 75% of the RAM, always
   leaving 3 GB to Windows; all CPU cores but 1–2; never more than Docker
   Desktop may use (README.md, "Resources"). While it runs, the launcher also
   makes Docker drop its file cache when Windows runs low on memory. Docker
   Desktop (WSL 2) may use only half the RAM by default (enough above ~18 GB of
   RAM); `setup` shows this.
   Raising it to the budget edits `%USERPROFILE%\.wslconfig` (keeping a
   `.bak`) and restarts WSL + Docker Desktop, stopping anything else running
   in WSL, so `--yes` never does it. It's optional (this pipeline needs only
   ~2.5 GB): **ask the user**, and only if they agree run
   `run_crisprscreen.bat setup --docker-memory`.

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

If they give Excel files: ask them to "Save As → CSV UTF-8", or, if Excel is
installed on this computer, convert them yourself through Excel's COM object in
PowerShell. Check the result opens as plain comma-separated text.

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
`%USERPROFILE%\Documents\CRISPR_screens\<screen_name>` (a drive with enough free
space: plan about 2× the size of the raw FASTQ files, plus results). Use a new folder: if the one
you picked already exists and holds anything, list it and ask the user before
using it (`init` never overwrites files, but it adds its own files and a git
repository to whatever is there). Then:

```powershell
& "<this folder>\run_crisprscreen.bat" init "$env:USERPROFILE\Documents\CRISPR_screens\<screen_name>"
```

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
  `config.yaml` (library files in `libraries\`). Without that column, every
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
  `libraries\<pool>.csv` per pool.
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

Start from `examples\*_example.csv`; save as UTF-8.

### 5. Edit `config.yaml`

- `project_name`
- `fastq_r1` (and `fastq_r2`) **or** `fastq_dir` — Windows paths are fine as
  they are (`D:\Sequencing\run7\Undetermined_S0_L001_R1_001.fastq.gz`), on a
  local or USB drive. The FASTQ(s) stay where they are: the launcher checks
  that they exist and reads them from that path; it never moves, copies,
  renames or deletes them, and never rewrites these keys. Tell the user where
  the file is read from and that it stays there. Keep the path stable (if the
  files are moved, update `config.yaml`). Refused: UNC paths (`\\server\share`)
  and mapped network drives, because Docker Desktop cannot mount them (copy the
  files to a local or USB drive); a path containing a comma. A FASTQ outside the
  project folder is mounted read-only into Docker (its folder, at
  `/inputs/fastq_r1`, `/inputs/fastq_r2` or `/inputs/fastq_dir`). A FASTQ that is
  a shortcut must be given by its real path. Several projects may point at the
  same raw file; it is not duplicated.
- `paired_end`, `contrasts`, `library` (glob for several pools)
- Optional `max_cpus` / `max_memory_gb` if the user needs the computer to stay
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

```powershell
cd "<project folder>"
& .\utils\run_crisprscreen.bat run --detach
& .\utils\run_crisprscreen.bat status        # every few minutes
```

Rough timing on a laptop: image build 10–30 min (first time only);
demultiplexing a full lane 20–90 min; the rest usually under an hour. On
Windows, a project on a fast internal SSD is much quicker than one on a USB
drive or network share. The computer must stay on and awake (plugged in, sleep
disabled) — if it sleeps, the run pauses; if Docker restarts, just run again
(it resumes).

`work\` is deleted after a successful run unless `keep_work: true` (in
`config.yaml`, default `false`); failed or stopped runs keep it. Re-running a
finished project recomputes everything (about 10–20 min for a full lane); to
resume only what changed (e.g. after adding contrasts), set `keep_work: true`
before the run. Intermediate files are in `results\intermediate_files\`.

### 7. Report results

Status `succeeded` → point the user to:
- `results\counts\qc\report\<project>_<library>_qc_report.html` — the QC
  report, one per library; written for non-specialists (start here);
- `results\demultiplexing\qc\multiqc\multiqc_report.html` — read QC;
- `results\counts\*.count.txt` — the counts;
- `results\mageck_test\<treatment>_vs_<control>\*.gene_summary.txt` — hits;
- `METHODS.md` — Materials and Methods of this run, ready to adapt for a paper.

Summarize from `results\counts\qc\tables\*_qc_summary.csv` (`mapping_pct`,
`library_skew_flag`, `undetected_flag`). Low mapping (well under ~50%) for a
sample usually means it belongs to another library, or a read-structure issue:
the read-classification plot tells which (mostly near-miss → read structure /
quality; mostly far-miss → wrong library).

### 8. When a run fails

`run_crisprscreen.bat status` says which step failed. Its own files are copied
to `metadata\failed_tasks\<timestamp>\<step>\` (`command.sh` = what ran;
`command.err` / `command.log` = what it said). Startup errors (bad samplesheet,
missing file, unknown library name, typo in a key) are printed at the end of
`metadata\nf_logs\nextflow_<timestamp>.log`. The project's `AGENTS.md`
("When a run fails") lists the common cases. Fix the input/config, run again:
finished steps are reused. If it looks like a bug in the pipeline itself, tell
the user and stop — don't patch this folder.

## Subagents

The agent definitions are in `.claude\agents\` of this folder, so an agent working here has them.
The same kit is in `agent\.claude\agents\`, which `init` copies into every project. The delegation
table and the mandatory triggers are in `CLAUDE.md` of this folder (and in `agent\CLAUDE.md` for
projects). Use `input-preparer` for the user's messy sample and library sheets (A.4), `scout` for
reading, `log-inspector` after runs (A.8), `results-reporter` for A.7, `code-supervisor` for pipeline
changes (section B), `statistician` for statistics.

## Hard rules

1. **Don't edit this pipeline folder** for a user's analysis. Per-project
   changes go in the project (`config.yaml`, `utils\` script overrides — see
   the project's `AGENTS.md`).
2. **Only start the pipeline with `run_crisprscreen`.** Don't run `docker run`,
   `nextflow`, `mageck`, `cutadapt`, `Rscript`, or the `bin\` scripts yourself.
3. **The raw FASTQs are the user's data.** They stay where the user keeps them
   (`fastq_r1`/`fastq_r2`/`fastq_dir` point to them). Never move, delete,
   modify or overwrite them, and don't change those keys unless the user asks;
   if the files are moved, update `config.yaml`.
4. **Ask before** deleting anything (`results\`, `work\` via `clean`,
   `.nextflow\`), installing software, editing `.wslconfig`, or renaming the
   user's files.
5. **Don't hand-edit `results\`**; it is regenerated. Intermediate files are
   inside it, in `results\intermediate_files\`.

6. **Never assume the experiment.** Ask the user (step 2), record it in the
   project's `CONTEXT.md`, and flag data that contradict it instead of
   reinterpreting the design.

## B. Changing the pipeline (only when asked to change it for the lab)

- This folder is self-contained: nothing in it may depend on a file outside it, because labmates
  sparse-clone only `WINDOWS/`.
- In the full repository this folder is the Docker port of the cluster pipeline in `HPC/`
  (`HPC/main.nf`, `HPC/bin/`, `HPC/modules/`). Read `HPC/AGENTS.md` for the
  pipeline's conventions and known traps (they all apply here). The resource rules of a local port
  are in "Porting a pipeline to the local Docker versions" below.
- In the full repository, `MAC/` and this folder are identical except for the launcher
  (`run_crisprscreen.ps1`/`.bat` here, `run_crisprscreen.sh` there), the
  README/AGENTS docs, the path style in `examples\config.yaml`, and the
  project-folder launchers that `init` writes. A change to `main.nf`,
  `nextflow.config`, `conf\`, `modules\`, `bin\` or `docker\` must be made in
  **both** folders: edit one, then run `bash tools/sync_os_folders.sh WINDOWS MAC` from the
  repository root (and `--check` to confirm). Usually the cluster version in `HPC/` needs
  it too.
- In the full repository, `bash tools/sync_os_folders.sh --check` verifies that the kit copies
  (`.claude\agents\` here, `agent\.claude\agents\` in every version) are identical.
- `bin\` is identical to the cluster's `HPC/bin/` except `write_methods.py` (Docker/host
  wording in the software paragraph + 3 references).
- Differences from the cluster version: one image for every step
  (`docker\Dockerfile`, no per-module `container`), local executor
  (`conf\docker.config`, resources from the container via
  `docker\run_in_container.sh`), no wall-time limits and memory
  requests set from measured cluster peaks (`conf\base.config`), run provenance (git commit, Docker version,
  image ID, host OS) passed in as `CRISPR_*` environment variables by the
  launcher and recorded by `main.nf` `runInfo()`.
- Tool versions: the `docker\env\*.explicit.txt` lock files are copied from
  the cluster's BioContainers images (`docker\env\images.txt`), so the Docker
  image holds the same package builds. To change a version, point
  `images.txt` at the new BioContainers tag (as the cluster module does), run
  `bash docker/env/make_locks.sh`, then sync. The launcher hashes
  `docker\Dockerfile` + `docker\env\*.explicit.txt` and rebuilds the image on
  everyone's next run when they change.
- Bash scripts must keep LF line endings (`.gitattributes` enforces this; the
  launcher also repairs `bin\*.sh` and `docker\*.sh` if a copy has CRLF).
- Keep `run_crisprscreen.ps1` ASCII-only and compatible with Windows
  PowerShell 5.1 (no `??`, `?:`, `&&`; native-command stderr must not be
  redirected under `$ErrorActionPreference = 'Stop'`).
- Validate with the synthetic test project (`test\README.md` in this folder) before committing.

### Porting a pipeline to the local Docker versions: CPU and memory

The Windows and Mac Docker versions are the reference for running a lab pipeline on a
personal computer. On the cluster each task got its own node; locally every
task shares one budget, so resource handling is the main thing a port must
get right. How this one does it (copy the pattern):

1. **One budget per run, chosen and enforced by the launcher**: what the
   pipeline can use (`docker/requirements.conf` `max_useful_memory_gb`), but
   never more than 75% of the computer's RAM, always leaving 3 GB to the OS (a
   plain percentage starves small machines; the lab's smallest machine is a
   6 GB Mac, so 3 GB is the floor); all logical cores but 1 (2 from 8
   cores up); never more than Docker Desktop's VM. `max_cpus` /
   `max_memory_gb` in `config.yaml` override it. Enforced with
   `docker run --memory/--memory-swap` (a runaway step is killed inside the
   container, the host never swaps) and `--cpuset-cpus` (tools that
   auto-detect cores only see the budget). Don't hand a pipeline more memory
   than it can use: Linux fills Docker's VM with file cache up to the
   container's limit while it streams FASTQs/BAMs, and Windows (WSL 2) or
   macOS gets that memory back only slowly. On a Windows PC with ~16 GB of
   other programs open, a 16 GB VM full of cache froze the machine during an
   image build (2026-10-16; Docker data and the page file on a hard disk).
   The Windows launcher therefore also runs a guard during builds and runs
   that drops the VM's cache (`wsl -d docker-desktop ... drop_caches`) when
   Windows' free memory falls below max(4 GB, 20% of RAM).
2. **Nextflow's local executor packs tasks into that budget**
   (`conf/docker.config`: `executor.cpus/memory` = budget minus ~1 GB for
   Nextflow itself; `process.resourceLimits` caps every request). Fan-out is
   automatic: the cluster's one-job-per-node tasks queue and run side by side
   as far as their requests fit.
3. **Requests = measured peaks, not node sizes.** Set each process's `memory`
   from real cluster runs (Nextflow `execution_trace.txt` `peak_rss`, or the
   LSF "Max Memory" of each job log) with ~2x headroom, and `cpus` from what
   the tool keeps busy (CPU time / wall time). Cluster requests are often
   10x the real peak; copying them only stops steps from running side by side.
   Keep `memory = { X * task.attempt }` with retry on exit 137.
4. **Every tool's thread count comes from `task.cpus`** (`--threads
   ${task.cpus}`), never "all cores" (`--cores 0`, `detectCores()`), unless
   the process reserves the whole budget (as `DEMULTIPLEX_*` does here).
   No script may submit its own jobs (`bsub`/`sbatch`) or size its own worker
   pool: every parallel unit of work is a Nextflow task.
5. **Declare a minimum** in `docker/requirements.conf` (`min_memory_gb`;
   decimals are allowed, e.g. 2.5): the launcher refuses to start below it,
   with an explanation, instead of failing halfway. A step whose real peak is beyond laptops (on the cluster:
   bowtie2 ~35 GB, ChIPseqSpikeInFree ~65 GB, csaw ~50 GB at their maximum;
   a STAR human index ~32 GB) needs a local alternative (prebuilt index,
   larger bins, another tool) before the pipeline is offered locally.
6. **Watch disk, not just memory.** Large intermediates that only feed
   another step (here: cutadapt's `--info-file`, ~100+ GB per full lane)
   should be streamed (named pipe in the container's `/tmp`), not written.

### Git

Commit each coherent change with a message that says what changed, why, and
how it was validated. Shell scripts must stay executable (`git add
--chmod=+x`) and keep LF line endings (`.gitattributes`); `.bat`/`.ps1` keep
CRLF. Don't commit work dirs, `.nextflow*`, project folders or data.
