# CRISPR_SCREEN — cluster (HPC) version

A reproducible, containerized Nextflow pipeline for pooled CRISPR screen
analysis on UPenn PMACS (LSF). Demultiplex → FastQC/MultiQC → `mageck count`
→ guide-representation QC (tables, plots, and a per-library HTML report) →
(optional) `mageck test` per contrast.

> **On your own computer instead of the cluster?** Use the Docker versions:
> the `WINDOWS/` or `MAC/` folder of the full repository. Same pipeline,
> same tools and versions, no bioinformatics software to install — only
> Docker Desktop, which their `setup` command installs for you. AI agents:
> start in that folder, at its `AGENTS.md`. The rest of this README
> describes the cluster version; the analysis itself (steps, anchor discovery,
> read classification, outputs) is described in `ANALYSIS.md`, which applies to all three.

This pipeline is **shared across the lab**: the repository is cloned to
`/project/rephilllab/CRISPR_SCREEN`, and this `HPC/` folder is what runs on
the cluster (`/project/rephilllab/bin/run_crisprscreen.sh` is a symlink to
`/project/rephilllab/CRISPR_SCREEN/HPC/run_pipeline.sh`).
It has no project-specific data or sample names baked in — every screen
(yours, a labmate's) runs it the same way, pointed at its own inputs via a
`config.yaml`.

## One-time setup

The cluster's one-time setup is described in [`SETUP.md`](SETUP.md).

Run the lab's setup script once, ever, regardless of how many lab pipelines
or tools ever get added later:
```bash
bash /project/rephilllab/setup_environment.sh
```
It adds a single guarded line to your `~/.bashrc` that sources the lab's
shared shell profile (`/project/rephilllab/.rephilllab.profile`), which puts
every lab pipeline's launcher on your `PATH` — `run_crisprscreen.sh` today,
future `run_rnaseq.sh`/`run_cutandrun.sh`/etc. land in the same shared
`bin/`. Safe to re-run: it's a no-op if you've already run it. It only ever
appends that one line — never overwrites or removes anything already in
your `~/.bashrc`. Then either run what it tells you to (`source ~/.bashrc`)
or open a new terminal.

That's the only thing you ever add to your own `~/.bashrc`, and the only
manual step. Any future lab-wide shell config goes into
`/project/rephilllab/.rephilllab.profile` itself, not your `~/.bashrc` —
everyone picks up an update to it at their next login.

<details>
<summary>What the setup script does, if you'd rather not run someone else's script blindly (reasonable)</summary>

It's two things: appending this exact line to `~/.bashrc`,
```bash
[ -f /project/rephilllab/.rephilllab.profile ] && . /project/rephilllab/.rephilllab.profile
```
guarded by `grep -qxF` so running it twice doesn't duplicate the line — read
the script itself at `/project/rephilllab/setup_environment.sh`, it's short.
</details>

This is PMACS support's recommended pattern for lab-specific config: global
`/etc/profile.d`-style files are for cluster-wide config only, and a lab's
own additions belong under its own `/project` space.

## Quick start (no HPC experience required)

1. Make a directory for your screen anywhere you have write access, e.g.
   `~/my_screen/`, and set it up in one command:
   ```bash
   mkdir ~/my_screen && cd ~/my_screen
   run_crisprscreen.sh init
   ```
   This creates the directory skeleton below, fetches `config.yaml` +
   `examples/` (starter CSVs) + `AGENTS.md`/`CLAUDE.md` (instructions for an
   AI agent helping you with this project), and `git init`s the directory —
   it does **not** run anything. It also copies the default barcode-position
   plot script `bin/plot_barcode_position.R` to `utils/plot_barcode_position.R`
   (only if that file is absent; its md5sum is recorded in
   `utils/.pipeline_defaults.md5`, so a later edit to the copy can be detected)
   and the agent definitions `agent/.claude/` into `.claude/` (only files that
   don't exist yet). `utils/` and `.claude/` go into the initial git commit.
   If an AI agent is helping you, it should read the new project's `CLAUDE.md`
   and `AGENTS.md` right after `init` and follow them from then on. For later
   sessions, start Claude Code inside the project folder so it loads them.
2. Copy the example CSVs you need from `examples/` into the project root
   (next to `config.yaml`) and fill in your real data:
   - **`samplesheet.csv`** — one row per sample: `barcode_id,barcode,sample_name`
     (barcodes only needed when the pipeline demultiplexes; optional `library`
     column — see "Several libraries" below).
   - **`library.csv`** — your sgRNA library: `gene,gRNA_ID,sequence`.
   - **`contrasts.csv`** (optional — only if you're ready to run `mageck test`):
     `treatment,control`, one row per comparison. Each value is matched
     (case-insensitive, substring) against the sample columns of your count
     table — it does not need to be an exact column name.
3. Edit `config.yaml`: set `project_name` and your reads: `fastq_r1` (and
   `fastq_r2` if paired-end) for a multiplexed FASTQ the pipeline
   demultiplexes, or `fastq_dir` for samples that are already demultiplexed
   (see below; when `fastq_dir` is set, `fastq_r1` is ignored). `samplesheet`/`library` already point at the
   files above by default. The FASTQ(s) stay where they are: the pipeline reads
   them from that path and never moves or copies them.
4. Run it:
   ```bash
   run_crisprscreen.sh
   ```
   That's the whole thing. It checks that your raw FASTQ exists (it stays where
   it is), submits itself to LSF, runs every tool in
   its own container (nothing to install), and resumes automatically if
   re-run after a failure.
5. Watch progress:
   ```bash
   bjobs -u $USER
   tail -f metadata/nf_logs/nextflow_<timestamp>.out
   ```
6. Results land in `./results/` (see layout below). Every plot is saved as an
   SVG and as a PNG (same name, ~150 dpi) next to it: per-library plots in
   `results/counts/qc/plots/`, the barcode-position plot in
   `results/demultiplexing/qc/demux_plots/`, and, with more than one library, the
   pool-vs-pool comparison plots in `results/counts/qc/plots/pools_comparison/`
   (`<project>_pools_<name>.svg/.png`). Nextflow's own
   `metadata/pipeline_info/execution_report.html` — per-task runtime, memory,
   and success/failure — is useful for troubleshooting or just seeing what ran.

If anything fails, **just re-run `run_crisprscreen.sh` again** — Nextflow's
`-resume` (already on by default in this wrapper) restarts only the failed
step, not the whole pipeline.

The Nextflow work folder (`work/`, on HPC `/project/rephilllab/nf_work/<user>/<project>/`)
is deleted after a successful run (together with Nextflow's `.nextflow/` cache)
unless `config.yaml` has `keep_work: true` (default `false`). Failed or stopped runs
always keep both, so a re-run resumes.
Re-running a finished project recomputes everything (~10-20 min for a full lane);
to resume only what changed (e.g. after adding contrasts), set `keep_work: true`
before the run. This is safe: `results/` holds copies, not links into the work
folder (`publish_mode` is `copy`; `auto` is also accepted and means copy on the
cluster). Intermediate files are in `results/intermediate_files/`.

## Where the raw FASTQ lives

The raw FASTQ(s) named by `fastq_r1`/`fastq_r2`, and an already-demultiplexed
`fastq_dir` folder, stay where you keep them. The pipeline never moves, copies,
renames or deletes them, and never rewrites those keys in `config.yaml`. It
checks that they exist and stops with a clear message if one is missing. A
relative path is relative to the project folder; an absolute path works too.
Keep the path stable: if you move the files, update `config.yaml`.

The cluster job reads these files under your account, so you need read
permission on them. Apptainer mounts cluster filesystems automatically
(`autoMounts = true`), so a file in another folder you can read is visible to
the job.

## Already-demultiplexed samples

If your samples were already demultiplexed (e.g. they were pooled on a lane
with another screen that did the demultiplexing), set `fastq_dir` in
`config.yaml`:

```yaml
fastq_dir: /path/to/demultiplexed_fastqs   # <sample_name>.fastq.gz (or .fastq); _R1/_R2 if paired_end
```

Whenever `fastq_dir` is set, demultiplexing and its QC are skipped and
`fastq_r1`/`fastq_r2` are ignored, whatever they say; FastQC/MultiQC, counting and
everything after run as usual. The folder stays where it is: the pipeline reads
its `*.fastq.gz` / `*.fastq` files in place and never moves, copies or renames
them, and never rewrites `fastq_dir`. `samplesheet.csv` then only needs
`sample_name` (plus `library` if used).
A sample with no matching FASTQ stops the run at startup, naming it.

## Several libraries (pools)

When samples were screened with different sgRNA libraries, give
`samplesheet.csv` a `library` column naming each sample's library, and point
`library:` at all the library files:

```
samplesheet.csv                          config.yaml
sample_name,library                      library: libraries/*.csv
AM_Pool_1_rep1,AM_Pool_1
AM_Pool_1_rep2,AM_Pool_1                 libraries/AM_Pool_1.csv
AM_Pool_2_rep1,AM_Pool_2                 libraries/AM_Pool_2.csv
```

A `library` value is the library's file name without `.csv`. Each library is
counted separately with **only its own samples**, so it gets its own count
table and its own guide-representation tables and plots — no near-zero
columns from other pools. Replicates and conditions that share a library land
in the same table, so replicate correlation and `mageck test` compare them.
Each contrast is tested within the library holding both of its samples; a
contrast whose treatment and control are in different libraries is rejected
(counts from different libraries aren't comparable). Library files no sample
refers to are skipped with a warning.

With more than one library, `results/counts/qc/report/` also gets one
`<project_name>_pools_comparison_qc_report.html` alongside the per-library
reports: the same plots, faceted one panel per pool (3 columns) so you can
spot a pool that looks different from the rest, plus the headline metrics
(Gini, 90/10 ratio, % undetected, mapping rate) as one bar chart per metric
across every pool. The per-library reports still have the full detail for
one pool (replicate correlation, the exact dropped-out guides, ...).

Without a `library` column, every sample is counted against every library
(identical to before for a single library). Use this when you don't know, or
want to check, which pool(s) each sample holds. In that mode:

- a sample with no read at all for one library is left out of that
  library's count table (with a warning in its mageck log) instead of
  failing the run — it's expected when samples hold different pools.
  It still shows in the breakdown below, as 0%.
- `results/counts/qc/tables/<project_name>_library_breakdown.csv` (one row
  per sample × library) and `..._library_breakdown_wide.csv` (one row per
  sample) give, for every sample, the % of its reads matching each library,
  matching any library (a guide shared by several libraries counted once),
  and matching none, plus its main library. The pools comparison report
  opens with the same breakdown as a stacked bar chart and table.
- each library's `mapping_pct` in `*_qc_summary.csv` counts only that library's
  reads, so a low value for a pool the sample does not hold is expected. Here,
  "Low mapping" means `pct_any_library` < 30, and a sample mainly matching
  another library gets a neutral "Mainly <library>" badge. With a `library`
  column (or one library), Low mapping means `mapping_pct` < 30.
- every sample appears in every pool's QC panels; a sample that holds few
  reads of a pool looks like heavy dropout there, so read those panels
  together with the breakdown.
- `mageck test` contrasts can't be resolved (each sample is in every count
  table); give the samplesheet a `library` column before testing.

The breakdown is also produced with a `library` column (each sample against
its own library only).

## Customizing a step for one project (without changing it for the lab)

Never edit the scripts in `/project/rephilllab/CRISPR_SCREEN/HPC/bin/` to
customize your own analysis — every labmate runs those. Every step of the
pipeline runs one script from `bin/`, and any of them can be replaced for
your project alone: copy it into your project's `utils/`, edit the copy, and
add a `config.yaml` line **named after the script** (its file name without
the extension) giving the copy's path.

```bash
# 1. copy the script into your project
cp /project/rephilllab/CRISPR_SCREEN/HPC/bin/plot_barcode_position.R utils/
# 2. edit utils/plot_barcode_position.R (add a plot, change colors, ...)
# 3. in config.yaml:
#      plot_barcode_position: utils/plot_barcode_position.R
# 4. re-run: that step (and anything downstream) re-executes
run_crisprscreen.sh
```

Note that on the cluster the work folder is deleted after a successful run (see
above), so a re-run only reuses the earlier steps if `keep_work: true` was set
before that run; otherwise the whole pipeline recomputes.

| Script in `bin/` (= `config.yaml` key) | Step                                                                                                                                              |
| ----------------------------------------| ---------------------------------------------------------------------------------------------------------------------------------------------------|
| `demultiplex.sh`                       | cutadapt demultiplexing (replacing it re-runs demultiplexing)                                                                                     |
| `barcode_position_qc.py`               | barcode-position counts (runs inside demultiplexing, so replacing it re-runs demultiplexing)                                                      |
| `cutadapt_stats.py`                    | demux metrics                                                                                                                                     |
| `plot_barcode_position.R`              | barcode-position density plot                                                                                                                     |
| `fastqc.sh`, `multiqc.sh`              | read QC                                                                                                                                           |
| `build_mageck_library.py`              | library CSV → MAGeCK format                                                                                                                       |
| `discover_anchor.py`                   | anchor discovery, one sample per library (only when `trim5` is unset; replacing it re-runs trimming and counting)                                 |
| `guide_anchor_trim.sh`                 | cutadapt trimming of every read through its library's anchor (only when `trim5` is unset; replacing it re-runs counting)                          |
| `mageck_count.sh`                      | `mageck count`, one run per library                                                                                                               |
| `library_breakdown.py`                 | with more than one library: % of each sample's reads matching each library (and none)                                                             |
| `guide_representation_wrangle.R`       | guide-representation QC: count table → CSV tables                                                                                                 |
| `guide_representation_plot.R`          | guide-representation QC: CSV tables → plots                                                                                                       |
| `guide_representation_report.R`        | guide-representation QC: tables + plots → one HTML report per library, with a description, how to interpret it, and troubleshooting for each plot |
| `guide_representation_compare_report.R` | with more than one library: tables from every library → one HTML report comparing all of them side by side (only runs when there is more than one) |
| `validate_contrasts.py`                | contrast validation                                                                                                                               |
| `mageck_test.sh`                       | `mageck test`, one run per contrast                                                                                                               |
| `software_versions.sh`                 | exact tool versions inside each container image (for METHODS.md)                                                                                  |
| `write_methods.py`                     | METHODS.md, the Materials and Methods of the run (see below)                                                                                      |
| `my_theme.R`                           | the look of every plot (must define a ggplot2 theme named `my_theme`)                                                                             |

Each script's header lists the arguments the pipeline passes it and the
files it must write; a replacement must accept the same arguments and write
the same files. A replacement may be in another language: `.py` runs with
`python3`, `.R` with `Rscript`, anything else with `bash`.

To keep the default script but pass it extra arguments, add
`<script name>_args`, e.g. `plot_barcode_position_args: "--bw=1"` or
`mageck_count_args: "--count-n"`. The `.sh` scripts hand extra arguments to
their tool (cutadapt, fastqc, multiqc, mageck). Changing a `<script>_args`
value in `config.yaml` re-runs that step.

Every plotting step is split from its computation: the plot step reads only
the tables the computation step wrote. So changing a plot script, the
theme, or a plot's `_args` re-runs only the plot step (and the report built
from it), with no recounting and no recomputing of statistics, provided the
earlier steps are still cached (`keep_work: true` was set for the previous
run). The plots are drawn in `guide_representation_plot.R` (per-library plots;
`--dpi` sets the PNG resolution, default 150), `guide_representation_compare_report.R`
(the pools comparison plots, drawn inside the report script) and
`plot_barcode_position.R`; `my_theme.R` is the shared look.

Replacements are just files, so they're picked up by your project's git repo
along with `config.yaml`. A custom plotting script can write extra plots —
the step publishes everything it produces. A replacement path that doesn't
exist, or a script path under a key that isn't a script name (a typo), stops
the run at startup with a message.

## Running your own analysis script (beyond the pipeline)

*New, not yet tested on the cluster.* For analysis the pipeline doesn't do, write
a script in your project's `utils/` and run it as one LSF job (never on the login
node):

```bash
run_crisprscreen.sh script [--cpus N] [--mem GB] [--time HH:MM] utils/<name>.<R|py|sh> [args...]
```

- Defaults: 1 core, 8 GB, 4 h. The launcher prints the job ID and the
  `tail -f`, `bjobs` and `bpeek` commands to follow it.
- The script must be inside the project folder. `.R` runs with Rscript, `.py`
  with python3, `.sh` with bash, in the project root, so it can read
  `results/...` and should write to its own `results/<name>/` folder (a
  convention, not enforced; don't write into the pipeline's own result folders).
- Containers: `rocker/tidyverse:4.3` for `.R` and `.sh` (R 4.3 + tidyverse; no
  Bioconductor) and the cutadapt biocontainer for `.py` (Python standard library
  plus cutadapt's dependencies only). A lab image with more tools (Bioconductor,
  samtools, deepTools, ...) is planned.
- Log: `metadata/script_logs/<script>_<YYYYMMDD_HHMMSS>.log`, with a header (date,
  project, script, md5, arguments, image, image ID, CPU/memory budget, working
  directory) and a last line `-- exit code N --`. The LSF output
  (`.lsf.out`, `.lsf.err`) is written next to it.

## Directory convention

The same layout is used by every pipeline built this way (this one, and
future RNA-seq/CUT&RUN/etc. pipelines) — not a CRISPR-screen-specific choice:

```
<project>/
├── .git/                 every project is its own git repo (config.yaml + your
│                          input CSVs are worth version-controlling; raw data,
│                          results, and run bookkeeping are not — see .gitignore)
├── config.yaml
├── samplesheet.csv, library.csv, contrasts.csv   your own small input CSVs — at
│                                                   the project root, not nested
├── original_data/
│   └── fastq/           per-sample FASTQs: this pipeline's demultiplexing
│                         output (the raw FASTQs are not copied here: config.yaml
│                         fastq_r1/fastq_r2/fastq_dir point to where they are kept)
├── metadata/             pipeline-generated run bookkeeping: nf_logs/ (LSF head-job
│                         logs), .nextflow.log, pipeline_info/ (execution report/
│                         trace/timeline/DAG), script_logs/ (logs of `script` runs) — NOT this pipeline's own inputs (see
│                         above); kept here for consistency with pipelines that do
│                         populate metadata/ with their own generated metadata
├── utils/                your project's replacements for pipeline scripts (see "Customizing")
│                          and your own analysis scripts (see "Running your own analysis script")
├── results/              final outputs, one subfolder per pipeline step (below),
│                          plus one subfolder per analysis script you run (results/<name>/)
│   └── intermediate_files/  tool-format conversions (e.g. the MAGeCK library file)
└── examples/             starter CSVs fetched by `init`, for reference
```

One exception, not relocatable: `.nextflow/` (Nextflow's own resume cache +
run history) always lives at the project root, alongside `config.yaml` —
Nextflow ties it to whatever directory it's launched from, and moving the
launch directory itself would break every relative path in `config.yaml`.
Verified directly (not assumed) that `-log` can redirect `.nextflow.log`
into `metadata/` but nothing moves `.nextflow/` independently.

`run_crisprscreen.sh init` creates this whole skeleton, fetches the starter
files, and git-inits the directory in one command — no more copy-pasting
several `cp` calls per new project. It also copies the default plot script
into `utils/` and the agent definitions into `.claude/` (missing files only);
`utils/` and `.claude/` go into the initial git commit.

## What the analysis does

See [`ANALYSIS.md`](ANALYSIS.md): Materials and methods (`METHODS.md`), what each step does, guide position (automatic anchor discovery), read classification (exact / near-miss / far-miss matches), and the output layout.

## Requirements this pipeline assumes about the HPC environment

Verified on PMACS (`pennhpc`, RHEL 9.6, LSF 10.1.0.13) on 2026-09-22:
Apptainer 1.4.5, Java 17, outbound network to quay.io/hub.docker.com.
**Nothing needs to be installed by an admin** — `run_crisprscreen.sh`
bootstraps a shared Nextflow launcher into `.nextflow-bin/` on first use.

If you're porting this to a different cluster, see `conf/pmacs.config` for
the one non-obvious fix required: PMACS compute nodes run OpenSSL in FIPS
mode, which makes Apptainer panic on every container run unless
`/usr/lib64/ossl-modules` is bind-mounted in. Check whether your cluster
needs the same.

## Known v1 limitations

- **Steps 4–5 (CRISPRO scoring, PDB structural coloring)** are not yet
  ported from the legacy `Scripts/step_4_run_crispro.py` /
  `step_5_color_pdb.py`. Neither tool has a public container, so these will
  run via `-profile pmacs,conda` against the lab's existing `crispro` conda
  environment once added — the interfaces in the legacy scripts are the
  starting point.
- **`samplesheet.csv` / `library.csv` schemas are now fixed** (named
  columns: `barcode_id,barcode,sample_name` and `gene,gRNA_ID,sequence`)
  rather than the legacy `config.py`'s configurable column-index scheme.
  If your CSV has different headers, rename them — this is deliberate: a
  fixed schema is much harder for a non-programmer to misconfigure than a
  1-indexed column-position setting.
- **Raw FASTQs are read in place.** Nothing is moved or copied before the LSF
  job is submitted (see `run_pipeline.sh`, reached as `run_crisprscreen.sh`).
  The launcher only checks that the paths exist and resolves them to absolute
  paths, which Nextflow passes to the job; the job reads them from there (see
  "Where the raw FASTQ lives").

## For pipeline developers (not typical labmates)

See `AGENTS.md` in this repo for the rules and conventions an agent (or a
person) should follow when changing the shared pipeline, and `agent/` for
the project-level instructions `init` copies into every project.

- `main.nf` — the DAG. Reads `samplesheet.csv` once at parse time; every
  downstream file name (fastq_args, sample_label_str) is built from that
  parsed order rather than trusted channel-emission order, so `mageck count`
  never silently misaligns a sample.
- `modules/local/*.nf` — one process (or process pair) per tool. Each
  declares its own container; there is no shared mega-image.
- `bin/*` — every step's script (bash wrappers around the tools, Python/R
  helpers), CLI-argument driven. Each module takes its script as a staged
  `path` input (resolved in `main.nf` by `stepScripts()`: the project's
  replacement if `config.yaml` has a key named after the script, else the
  `bin/` default) and runs it via `bash`/`python3`/`Rscript` explicitly — so
  a replacement needn't be executable or on `PATH`, and `-resume` keys the
  step's cache on the script's content. Plotting scripts
  take the theme as `--theme=<file>` (defining `my_theme`); run standalone
  without it, they fall back to `my_theme.R` next to themselves via
  Rscript's `--file=` argument. Don't use `Sys.which()` for this — it
  silently returns `""` for a non-executable file (a real past failure; see
  git history).
- `conf/base.config` — per-process cpu/memory/time, with automatic
  retry-with-escalating-memory on OOM (LSF exit 137) or wall-time kill (140).
- `conf/pmacs.config` — LSF executor + Apptainer + the FIPS bind fix.
- `conf/conda.config` — fallback profile for tools without a container
  (currently unused by the default DAG; reserved for CRISPRO/PyMOL).
- `run_pipeline.sh` — this pipeline's implementation, reached by labmates as
  `run_crisprscreen.sh` (a symlink in the shared `/project/rephilllab/bin/`,
  which is the only thing ever added to `PATH`). `init` scaffolds a project,
  fetches `examples/`, copies the default plot script into `utils/` and the
  agent definitions into `.claude/` (missing files only), and `git init`s it
  (`utils/` and `.claude/` go into the first commit), then exits; the plain form checks
  the raw FASTQ paths and submits the real run. Resolves through the symlink via
  `readlink -f` before taking `dirname` — plain `dirname "$BASH_SOURCE"`
  would resolve to the shared `bin/` directory, not this repo, breaking
  every path built from `PIPELINE_DIR`. A future pipeline's own
  `run_pipeline.sh` gets its own symlink (`run_rnaseq.sh`, etc.) added to
  the same shared `bin/` — no change needed to any existing pipeline or to
  labmates' `PATH`.

To add a new step: write the tool logic as a `bin/` script if it needs more
than a couple of lines of shell, wrap it in a `modules/local/<name>.nf`
process with a pinned container, `include` it in `main.nf`, and wire it into
the `workflow {}` block.
