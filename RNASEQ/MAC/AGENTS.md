# AGENTS.md — RNA-seq pipeline, Docker version for **Mac** (Rephill Lab)

This folder is the lab's RNA-seq pipeline (Nextflow: FASTQ pairs or BAMs to trimming, HISAT2 alignment, filtering,
duplicate marking, the strandedness check, featureCounts, a merged counts table and MultiQC), packaged to run on a
labmate's own Mac inside Docker. It can also run DESeq2 differential expression (`run_deg: true`, off by default; see README.md, "Differential expression"). The user is usually a **wet-lab scientist**: do
the technical work for them, explain in plain language, and ask them only for what you can't know (their samples, their
files, their questions).

First time on this Mac (Docker not set up, or never checked)? Follow `SETUP.md` in this folder once.

Your task is almost always **A. Set up / run / check an analysis**. Follow the runbook below. You never edit anything in this
pipeline folder for that. The pipeline code is not edited during an analysis; if a change to the pipeline seems needed, ask
the maintainers.

If the current folder is a **project** folder (it has `config.yaml`, `samplesheet.csv`, `results/`), read that project's
own `AGENTS.md` instead.

## How to call the launcher

Everything goes through one command, `run_rnaseq.sh`, in this folder (full path: `<this folder>/run_rnaseq.sh`). From the
agent's shells, call it with `bash` and its full path. That works even if the file lost its executable bit, and it doesn't
depend on the user's PATH:

```bash
bash "/Users/<user>/PIPELINES/RNASEQ/MAC/run_rnaseq.sh" status
```

| Command | What it does |
|---|---|
| `run_rnaseq.sh setup --yes` | one-time per Mac: Docker Desktop check/install (with Homebrew, after asking), Rosetta 2 on Apple chips, pull the pipeline's Docker image (`docker/image.txt`), print the tool versions, add the command to PATH (`~/.zshrc`) |
| `run_rnaseq.sh init "<project folder>"` | create a project folder: `original_data/fastq/`, `metadata/`, `results/`, `intermediate_files/`, `utils/`, `examples/`, `config.yaml`, `samplesheet.csv`, launchers, git |
| `run_rnaseq.sh run --detach` (in the project folder) | check the input files and the disk, check Docker's memory, start the pipeline in Docker, return immediately |
| `run_rnaseq.sh status` (in the project folder) | running / succeeded / failed, step counts for this run, the failed step, the last 25 log lines, the Docker work volume's name and size |
| `run_rnaseq.sh logs` | follow the live log (blocks until the run ends — avoid in an agent; use `status`) |
| `run_rnaseq.sh stop` | stop a run (finished steps are kept; the next run resumes) |
| `run_rnaseq.sh clean --yes` | delete Nextflow's work folder and resume cache to free disk: the project's Docker work volume, and `metadata/work/` and `metadata/.nextflow/` if present — the next run starts from scratch, results are kept; without `--yes` it asks first; ask the user first |
| `run_rnaseq.sh script utils/<name>.<R\|py\|sh> [args]` (in the project folder) | run the project's own analysis script in the pipeline's image (foreground; log in `metadata/script_logs/`); see "Analyses beyond the pipeline" below |
| `run_rnaseq.sh help` | the command list |

Options: `--detach` (run: start and return), `--yes` (answer yes to every question: only after the user has agreed to what
is being asked; the questions are the Homebrew install of Docker Desktop, the Rosetta 2 install and the deletion by
`clean`), `--pause` (wait for Enter before closing; used only by `RUN_ANALYSIS.command`).

Every project folder also gets its own launcher, `utils/run_rnaseq.sh` (same commands, always for that folder; agents run
`bash utils/run_rnaseq.sh` from the project folder), and `RUN_ANALYSIS.command` (double-click for humans: it waits for Enter
at the end, so **agents must not call `RUN_ANALYSIS.command`**).

Long commands: `setup` (image pull) and a run can take a long time. Run `setup` in the background or with a long timeout,
and always start runs with `--detach`, then poll `status` every few minutes.

## A. Runbook: from "here are my FASTQs" to results

### 1. One-time Mac setup

The steps below are this pipeline's part of the setup (`SETUP.md` in this folder has the full text).

1. Check: `docker --version` and `docker info`. If both work, go to step 3. (Docker Desktop's command-line tool may be in
   `/Applications/Docker.app/Contents/Resources/bin` or `~/.docker/bin`; the launcher looks there.)
2. `bash "<this folder>/run_rnaseq.sh" setup --yes`. What can happen:
   - Docker Desktop not installed → with Homebrew, the launcher asks and then runs `brew install --cask docker` (it may
     ask for the user's password: tell them). Without Homebrew it stops and prints the download page. Then the user must
     **open Docker from Applications once** and accept its terms (no account needed — "Skip"). Explain this plainly, wait
     for them, and re-run `setup`.
   - Docker installed but not running → it runs `open -a Docker` and waits up to 4 minutes.
   - Apple chip (M1 to M4): the image is `linux/amd64` and runs through Rosetta 2. `setup` offers to install Rosetta
     (`softwareupdate --install-rosetta`); without it Docker uses much slower emulation. In Docker Desktop → Settings →
     General, "Use Rosetta for x86_64/amd64 emulation on Apple Silicon" should be on.
3. `setup` pulls the image (`ghcr.io/genesig/omics`, tag in `docker/image.txt`, always run as `linux/amd64`) and prints the
   versions of the tools. Re-running `setup` is harmless. A run also pulls the image if it is missing, so a forgotten
   `setup` only costs time.
4. Resources: every run gets one budget, enforced by Docker as a hard limit. It is the smallest of: 75% of the Mac's RAM
   (always leaving at least 3 GB to macOS), all logical cores but one (two from 8 cores up), what Docker Desktop may use
   (its Memory setting, less about 0.5 GB for Docker's own system), and the most this pipeline can put to use
   (`max_useful_memory_gb` in `docker/requirements.conf`). `max_cpus` and `max_memory_gb` in `config.yaml` can only lower
   it. The run prints the budget it chose. The pipeline needs at least `min_memory_gb` (`docker/requirements.conf`) and
   refuses to start below that, with an explanation. Both figures are placeholders until a pilot run on a real Mac has
   been measured.
5. **Docker Desktop's memory.** The launcher never changes Docker Desktop's settings. It only recommends a Memory value
   for analyses, by the Mac's RAM: under 12 GB → 6 GB; 12 to under 24 GB → 8 GB; 24 GB or more → 16 GB. The user sets it
   in Docker Desktop → Settings → Resources → Memory, then clicks Apply & restart. Changing it is the user's decision:
   **ask the user first**. If Docker's memory holds the budget back, the run says so in a note (with the setting to
   change); the run stops only when the budget is below the minimum.

### 2. Gather what you need from the user (ask; never guess)

- **A name for the project** → `project_name` (letters, digits, `_`; used as a prefix on the counts table and the other
  project-level outputs) and the project folder's name.
- **Where the reads are**, and whether each sample is:
  - **FASTQ pairs**: `fastq_1` and `fastq_2` (absolute paths, one pair per library), paired-end. If the data are
    single-end, say so before going further: this pipeline reads paired-end data only; or
  - **a BAM**: the `bam` column, with the FASTQ columns left empty. The BAM is converted back to paired FASTQ first and
    then goes through the normal steps.
- **The reference files** (ask for the paths, never guess them): the HISAT2 index **prefix** of the genome
  (`hisat2_index`) and the gene annotation GTF (`gtf`). Mac paths: `/Volumes/<name>/...` for an external drive or a share
  mounted in Finder, `/Users/<name>/...` for the home folder.
- **The samples**: names, which metadata each one has (condition, treatment, cell line, guide, replicate), and which one
  is to be excluded. The sequencing-core sheet often holds these; ask for the sheet.

If they give Excel or Numbers files: ask them to export CSV UTF-8, or convert them yourself inside the pipeline image.
Check the result opens as plain comma-separated text.

**Then ask about the experiment itself (ask; never assume).** Only the user knows how it was done; ask even when it seems
obvious:

- **The question:** what the experiment is meant to show (a treatment's effect, a knockout's effect, a time course), and
  which genes they expect to change or not to change (positive and negative control genes).
- **The RNA and the library:** the source (cell line, tissue, sorted cells), the extraction, the library kit, whether it was
  poly-A selected or total RNA, and whether it is **stranded** (which kit or protocol). The strandedness decides how the
  counts are made.
- **The sequencing:** paired-end, read length, which run(s) and lanes, whether the run also holds other people's samples,
  and whether the FASTQs are demultiplexed.
- **Each sample:** its treatment, dose and timepoint; which samples are biological replicates and which are technical;
  and anything that went wrong (low input, a failed library, a re-sequenced sample).
- **The comparisons:** which treatment-vs-control comparisons answer the question. The pipeline does not make them yet; say
  so, and record the plan in `CONTEXT.md`.

- Ask in short batches (3-6 questions at a time), in plain language, and say why each answer matters. Keep asking until
  the answers are specific; "the usual" or "like last time" is not an answer, so ask what it was.
- **Never infer the design from names** (`A`/`B`, `rep2`, `ctrl`, `DMSO`, `1uM`, `KO`, `hi`/`lo`) or from what is "typical"
  for RNA-seq. Ask what each name means.
- Use only what the user confirmed. Where an answer is missing, say what you would do instead and wait for an explicit OK;
  never fill the gap silently, and never pick the controls, comparisons or replicates yourself.
- **Write the answers down** in `CONTEXT.md` in the project folder (dated, in the user's words, with open questions marked
  as open). Read it at the start of every session; ask the user to confirm it is still right before you rely on it, and add
  to it as you learn more. It is tracked by git with the inputs.
- **When the data seem to contradict the context, flag it; don't reinterpret.** Examples: a strandedness call that differs
  from the configured value (`results/qc/strandedness/strandedness.tsv`, `mismatch`); a low alignment rate
  (`results/qc/hisat2/alignment_summary.tsv`); a sample whose counts look unlike its replicates. Show the user the numbers
  and the file they come from, say what the data would be consistent with, and ask. The user's description stands until
  they change it; don't rename samples, drop them or change comparisons on your own.
- Subagents cannot talk to the user: put the confirmed context (or the path to `CONTEXT.md`) in every brief, and pass their
  "Questions for the user" back to the user rather than answering them yourself.

### 3. Create the project folder

Put projects **outside** the pipeline's git clone, e.g. `~/rnaseq/<project_name>` (a disk with enough free space: see the
disk note in step 6; not an iCloud-synced folder). Use a new folder: if the one you picked already exists and holds
anything, list it and ask the user before using it (`init` never overwrites files, but it adds its own files and a git
repository to whatever is there). Then:

```bash
bash "<this folder>/run_rnaseq.sh" init ~/rnaseq/<project_name>
```

`init` writes the project's `config.yaml` and `samplesheet.csv`, copies the templates into `examples/` (never overwriting
them), and copies the project's `AGENTS.md`/`CLAUDE.md` and the agent kit (`.claude/agents/`) when this folder has them
(`agent/` in this folder). The CSV the project uses is the copy you make or fill in in the project folder itself.

Immediately after `init`, **read the new project's `CLAUDE.md` and `AGENTS.md` and follow them from then on** (their
delegation table and triggers included): users usually ask you to set the project up before running, so the project's
`.claude/agents/` and instructions appear mid-session. If Claude Code was started in this pipeline folder, its
`.claude/agents/` is the same kit, so the agents are available now. If they are not listed as available subagent types,
either spawn a general-purpose agent whose prompt is the agent file's body (after the frontmatter) plus the brief, with
the model from the file's frontmatter, or ask the user to restart Claude Code in the project folder (recommended for every
later session: start Claude Code inside the project folder, so it loads that project's `CLAUDE.md` and agents).

### 4. Write the samplesheet (in the project folder)

The full column rules and the enforced checks are in the project's `AGENTS.md` ("Input files"); the example is
`examples/projects/gg_ezh2i_rnaseq/samplesheet.csv`. Summary:

- **`samplesheet.csv`** (one row per library; named by `samplesheet:` in `config.yaml`):
  `sample_id,fastq_1,fastq_2,bam,exclude` followed by the user's metadata columns. `sample_id` is unique and uses letters,
  digits, `.`, `_`, `-` (it names every output file) and carries no read suffix. Fill `fastq_1` and `fastq_2` **or** `bam`,
  never both and never neither. `exclude` (`true`/`false`) drops a row from the whole run. The metadata columns are copied
  into `results/samplesheet.tsv` and are used by the user's later analysis; they change nothing in the counts.

Start from `examples/projects/gg_ezh2i_rnaseq/samplesheet.csv`; save as UTF-8, no comment lines.

### 5. Edit `config.yaml`

Every key has a default in the pipeline's `nextflow.config`; list only what you change. The key choices:

- `project_name`, `samplesheet`.
- **Reference files**: `hisat2_index` (HISAT2 index prefix) and `gtf` (both required). Write Mac paths in full
  (`/Volumes/<name>/...`, `/Users/<name>/...`) or relative to the project folder. The launcher warns when a file is
  missing, or when `hisat2_index` has no `.ht2` files next to it.
- **Strandedness**: `strandedness` (`unstranded`, `forward` or `reverse`; default `unstranded`, the legacy setting) sets how
  featureCounts counts. Ask the user which library preparation was used; never guess it. The automatic check
  (`strand_call_threshold`, `infer_experiment_reads`) only reports.
- **Input files stay where they are.** The FASTQs and BAMs are never moved, copied, renamed or deleted; the top-level folder
  of each path is mounted read-only into Docker at the same path. Tell the user where the files are read from and that they
  stay there. Keep the paths stable. Refused: network paths (`smb://`, `//server/share`); mount the share in Finder and use
  its `/Volumes/<name>/...` path.
- **`fastq_dir`** is not used by this version: leave it `null`.
- **Output folders must be relative**: `outdir: ./results`, `metadata_dir: ./metadata`, `intermediate_dir:
  ./intermediate_files`. The launcher refuses an absolute path there (such as `/Users/...`).
- **Duplicates**: `mark_duplicates: true` (marked, and counted) and `remove_duplicates: false` (legacy). Setting
  `remove_duplicates: true` changes the counts: ask first. `mark_duplicates: false` removes the duplicate metrics.
- **Trimming and alignment**: `trim: true`, `trim_adapter_set: TruSeq3-PE-2`, `trim_minlen: 25`, `mapq: 30`,
  `hisat2_args: ''`. Changing one of these changes the alignment: ask the user first.
- **Counting**: `featurecounts_feature_type: exon`, `featurecounts_attribute: gene_id`, `featurecounts_args: ''`. Changing a
  counting option changes every count: ask first.
- **chrM**: `remove_chroms: ''` keeps all contigs (the example). Setting `chrM` reproduces the legacy BAMs and changes no
  count on the T2T annotation, but on mouse (mm39) it zeroes the mitochondrial genes: keep it empty there.
- **Run bookkeeping**: `keep_work: true` (default) keeps the work folder and resume cache so a re-run resumes; `work_volume:
  true` (default) keeps them in the project's Docker volume on Docker's disk, `work_volume: false` in `metadata/work/`;
  `publish_mode` is `auto` by default (copies with the work volume; hard links when the disk allows without it).
- **`publish_intermediates`**: `false` by default (the trimmed FASTQs, the FASTQs made from BAMs and the aligned and filtered
  BAMs stay in the work folder). `true` publishes them to `intermediate_files/`; it needs disk space.
- **Resources (optional)**: `max_cpus` and `max_memory_gb` lower the run's budget. Ask the user before raising anything.

### 6. Check, then run

Check before the first run (the launcher also checks, but catching it now is faster): headers exactly right; `sample_id`
unique; FASTQ and BAM paths exist (the launcher names any file it cannot find); `hisat2_index` and `gtf` exist;
`strandedness` is one of the three values, and the user has confirmed it.

```bash
cd "<project folder>"
bash utils/run_rnaseq.sh run --detach
bash utils/run_rnaseq.sh status        # every few minutes
```

The run first checks the inputs and prints the disk check: the free space on Docker's disk (the disk with Docker Desktop's
disk image, which holds the work volume) and on the project's disk (which holds the results), against the total size of the
FASTQ/BAM files in the samplesheet, using the factors in `docker/requirements.conf` (`disk_factor`, `disk_results_factor`).
It warns when the free space is short, and stops below `disk_hard_min_gb` (20 GB) on the project's disk. These factors are
estimates pending the next run. Free space can be given back with `clean` (ask first) or with `keep_work: false`, which
deletes the work files after a successful run.

Poll `status` every few minutes. The Mac must stay on and awake (plugged in, the lid open): if it sleeps, the run pauses; if
Docker restarts, just run again (it resumes).

Nextflow's work folder and its resume cache are kept by default (`keep_work: true`), so a re-run resumes and only redoes what
changed. Failed or stopped runs always keep them. With `work_volume: true` (the default) they are in the project's Docker
volume, on Docker Desktop's own disk; `status` shows its name and size. Files in `results/` are copies of the work files
(`publish_mode: auto`). With `work_volume: false` they are `metadata/work/` and `metadata/.nextflow/` in the project folder.

### 7. Report results

`status` says `succeeded` → point the user to:
- `results/qc/multiqc/multiqc_report.html` — the overview of every QC output (start here);
- `results/counts/<project_name>_counts.txt` — the counts table (one row per gene, one column per sample), with
  `results/counts/gene_length.txt` and `results/counts/featurecounts_summary.tsv`;
- `results/qc/strandedness/strandedness.tsv` — the strandedness call and the `mismatch` column;
- `results/qc/hisat2/alignment_summary.tsv` — the alignment rate per sample and the warning flag;
- `results/qc/picard/` — the duplicate metrics; `results/bam/final/` — the BAMs that were counted;
- `metadata/pipeline_info/execution_report.html` — how long each step took and how much memory it used. A resumed run keeps
  the previous run's trace, report and timeline next to it, named with that run's timestamp.

Summarise from `results/qc/` (the alignment rates, the strandedness table, the duplicate metrics) and the MultiQC report.
Check the strandedness first: a `mismatch` means the configured value and the data disagree, and that is the user's
decision. For statistics or comparisons with earlier numbers, use the `statistician` and `output-verifier` agents (see
"Subagents").

### 8. When a run fails

`bash utils/run_rnaseq.sh status` says which step failed (`FAILED: <step> (exit N)`) and prints the last log lines. The
files to read, in order:

1. `metadata/nf_logs/last_run.status`: `running`, `succeeded`, `failed (exit N)`, `failed to start`, `stopped by the user`,
   or `interrupted` (Docker or the Mac stopped mid-run; run again to resume).
2. `metadata/nf_logs/nextflow_<timestamp>.log`: everything Nextflow printed. Startup errors (bad samplesheet, missing file, a
   typo in a key) are printed at its end.
3. `metadata/failed_tasks/<timestamp>/<step>/`: the failed step's own files (`command.sh` = what ran; `command.err` /
   `command.log` = what it said).
4. `metadata/pipeline_info/execution_trace.txt`: one row per step (status, exit code, `peak_rss` = the memory it used,
   `realtime`). Use it to see whether a step ran out of memory or disk.

Common cases (the project's `AGENTS.md`, "When a run fails", lists more):
- **Exit 137, 140 or 143**: out of memory or killed. These are retried automatically with more memory, up to the budget.
  If the step still fails, check the budget in the run banner and the step's `peak_rss`. Raising the budget is the user's
  decision (`max_memory_gb`, or Docker's Memory setting, which the user changes in Docker Desktop): ask first.
- **Exit 3**: a known data problem, not retried; the message says what to fix. Typical: an empty BAM (ALIGN); a BAM with no
  read pairs (BAM_TO_FASTQ); no contig left after `remove_chroms`; an empty count table; a gene column that differs between
  samples (COUNTS_MERGE).
- **Disk full**: the run stops and can be resumed after space is made. Check the free space on Docker's disk and on the
  project's disk (`status` and the disk check in `run`); `clean` (ask first) deletes the work volume (or the work folder with
  `work_volume: false`).
- **Startup errors**: fix the named column or key in the project folder, then run again.
- **A strandedness WARN** at the end of a run is not a failure: the counts are as configured. Tell the user the numbers.

Fix the input or config, then run again: `run --detach`. Finished steps are reused. If it looks like a bug in the pipeline
itself, tell the user and stop: don't patch this folder.

### 9. Analyses beyond the pipeline

The user may want a table or plot the pipeline doesn't make. That is allowed, through the project's `script` command:
`bash utils/run_rnaseq.sh script utils/<name>.<R|py|sh>` runs the script in the pipeline's image (**Tier 1**), and the log goes
to `metadata/script_logs/`. The script lives in `utils/`, writes only to a new `results/<name>/` folder, and is designed,
written, run and reviewed as described in "Analyses beyond the pipeline" in the project's `AGENTS.md` (`init` copies it
there). The `script` command is foreground only (no `--detach`).

Running `Rscript`, `python` or a shell script **on the host** (**Tier 2**) is only done with the user's OK, with the package
versions recorded next to the output. Never run it without asking the user.

### 10. Differential expression (optional, `run_deg`)

Set `run_deg: true` only when the user asks for differential expression. You need `deg_comparisons` (CSV: `model`, `design`,
`comparison` and one contrast per row) and the sample metadata with the design columns (`deg_metadata`, or the samplesheet's own
columns). The other inputs and the outputs are in README.md, "Differential expression". Choosing the comparisons, the design,
`deg_alpha`, `deg_lfc` or the shrinkage is a scientific decision: ask the user, do not choose it for them. To run only the DEG steps
on a counts table the user already has, set `deg_counts` too.

## Subagents

The agent definitions are in `.claude/agents/` of this folder, so an agent working here has them. The same kit is in
`agent/.claude/agents/`, which `init` copies into every project. The delegation table and the mandatory triggers are in
`CLAUDE.md` of this folder (and in `agent/CLAUDE.md` for projects). Use `input-preparer` for the user's messy sample
sheets (A.4), `scout` for reading, `log-inspector` after runs (A.8), `results-reporter` for A.7, `code-supervisor` for
pipeline changes, `statistician` for statistics and for any design beyond counts.

## Hard rules

1. **Don't edit this pipeline folder** for a user's analysis. Per-project changes go in the project (`config.yaml`,
   `utils/` script overrides, see the project's `AGENTS.md`).
2. **Only start the pipeline with the launcher's `run`**, and run the project's own analysis scripts with its `script`
   command. Don't run `docker run`, `nextflow`, `hisat2`, `samtools`, `featureCounts`, `trimmomatic`, `picard`, or the `bin/`
   scripts yourself. Commands are taken from `run_rnaseq.sh` only.
3. **The raw FASTQs and BAMs are the user's data.** They stay where the user keeps them (the paths in the samplesheet point
   to them). Never move, delete, modify or overwrite them, and don't change those paths unless the user asks; if the files
   are moved, update the paths in the samplesheet and `config.yaml`.
4. **Ask before** deleting anything (`results/`, the Docker work volume, `metadata/work/` and `metadata/.nextflow/` via
   `clean`), installing software (Homebrew, Docker Desktop, Rosetta 2), changing Docker Desktop's memory, changing
   `strandedness`, a duplicate setting, a trimming or counting option, or renaming the user's files.
5. **Don't hand-edit `results/`**; it is regenerated. Published intermediate files are inside `intermediate_files/`.
6. **Never assume the experiment.** Ask the user (step 2), record it in the project's `CONTEXT.md`, and flag data that
   contradict it instead of reinterpreting the design.
7. **Never assume a path.** Reference files and FASTQ locations come from the user.
8. **Never change `strandedness` on your own.** The strandedness check reports what the libraries look like; it never changes
   the counts. Whether the configured value is right is the user's decision.
9. **The Mac's settings are the user's.** Don't change Docker Desktop's memory or other settings by hand, and don't install
   software without asking. The launcher never changes Docker's settings; if it says Docker's memory is too low, tell the
   user what to set and let them set it.
10. **No run without the user's OK.** Starting a run needs the user's explicit OK, after the inputs have been checked with
    them.
