# RNA-seq pipeline for Mac (Docker)

The Rephill Lab's paired-end RNA-seq pipeline: from FASTQ pairs or BAMs to one counts table per project.
Steps: FastQC on the raw and trimmed reads, Trimmomatic trimming, HISAT2 alignment, samtools filtering, Picard
duplicate marking, an automatic strandedness check (RSeQC), featureCounts, a merged counts table, and MultiQC. It is
the same pipeline as on the cluster (same steps, same outputs; the tools run in one Docker image here, see
"Differences from the cluster"), so you don't install any bioinformatics software yourself.

Optionally (`run_deg: true`) it then runs a DESeq2 differential expression analysis of the counts table; by default the counts table is the end point (see "Differential expression").

> **Status:** the Mac version has not yet been run on a lab Mac. Before you trust its numbers for a dataset, check it
> against the benchmark (section "Checking a run against the benchmark").

> **First time on this Mac?** Follow `SETUP.md` in this folder once (Docker Desktop, Rosetta on an Apple-chip Mac, where
> Docker's data lives, memory). Then come back here.
>
> The agent instructions for this pipeline (`AGENTS.md`, `CLAUDE.md` and the agent definitions in `.claude/agents/`) are
> copied into each new project by `init`. The steps below are complete on their own.

## What you need

- macOS 13 (Ventura) or newer, on an Apple chip or an Intel chip.
- **Memory:** the full RNA-seq data (for example the 12-library benchmark) needs a larger budget than a small test. The
  minimum budget the run needs, and the most it can use, are in `docker/requirements.conf`. Below the minimum the launcher
  refuses to start a run and says so. The budget a Mac gets depends on its RAM and Docker Desktop's Memory setting (see
  "Resources").
- Free disk space: on the disk that holds Docker Desktop's disk image (the Mac's internal disk by default), room for
  the work files (with the default `work_volume: true`; the RNA-seq figure is not measured yet, the launcher computes
  the check before each run; see "Disk space") plus the Docker image (several GB, once); on the disk that holds your
  project, room for the results (the final BAMs are the largest part: the cluster benchmark of 12 libraries, 57 GB of
  FASTQ, gave 32 GB of final BAMs) and at least 20 GB free. The raw FASTQs are not copied.
- Internet the first time (to download Docker and the pipeline's image).
- A reference set: a HISAT2 index of the genome (the index **prefix**) and the gene annotation GTF (featureCounts and the
  strandedness check use it). Any location works as long as you write its path in `config.yaml`; the layout in
  `SETUP.md` ("Where data and projects go") is a suggestion.

## One-time setup (per Mac)

1. Get the pipeline (once), as a sparse clone so that only this folder is downloaded:
   ```
   git clone --filter=blob:none --sparse <the lab's GitHub URL> ~/PIPELINES
   cd ~/PIPELINES
   git sparse-checkout set RNASEQ/MAC
   cd ~/PIPELINES/RNASEQ/MAC
   ```
   Then work in the `RNASEQ/MAC` folder; `git pull` later updates it.
   (The folder name of the clone is whatever the lab's repository is called.)
2. In the `RNASEQ/MAC` folder, run in Terminal:
   ```
   bash run_rnaseq.sh setup
   ```
   - If Docker Desktop isn't installed, it offers to install it with Homebrew (if you have Homebrew; it asks first),
     or tells you where to download it (choose *Apple chip* or *Intel chip*). Then open **Docker** from Applications
     once, accept its terms (no account needed: click *Skip*), and run `setup` again.
   - On an Apple-chip Mac it offers to install Rosetta 2 if needed (the pipeline's image is for Intel chips and runs
     through Rosetta). Answer yes; it asks for your password.
   - It then downloads the pipeline's Docker image (`ghcr.io/genesig/omics`, the tag is in `docker/image.txt`): several
     GB, once. It takes a while; keep the Mac awake and online.
   - It prints the versions of the main tools, and adds `run_rnaseq.sh` to your PATH (in `~/.zshrc`), so new Terminal
     windows can use it from any folder.

## Analyze a dataset

1. **Create a project folder** (one per dataset, *outside* the `RNASEQ/MAC` folder):
   ```
   run_rnaseq.sh init ~/rnaseq/my_dataset
   ```
   It creates the folder skeleton (`original_data/fastq/`, `metadata/`, `results/`, `intermediate_files/`, `utils/`),
   `config.yaml`, `samplesheet.csv`, `examples/` (starter templates), `AGENTS.md` and `CLAUDE.md` with the agent kit,
   `RUN_ANALYSIS.command` (double-click) and `utils/run_rnaseq.sh` (the command used below). When git is installed, it
   also starts a git repository in the folder with a first commit. It doesn't run anything.
2. **Write your samplesheet** `samplesheet.csv` in the project folder: one row per library. Start from
   `examples/projects/gg_ezh2i_rnaseq/samplesheet.csv` (or the template in `examples/`). See "Input table" below.
   Numbers or Excel: *Export as CSV* (UTF-8).
3. **Edit `config.yaml`** (TextEdit or any editor is fine). Set `project_name`, `hisat2_index` (the HISAT2 index prefix,
   e.g. `~/REFERENCES/hisat2/t2t/chm13-2`) and `gtf`. Set `strandedness` (`unstranded`, `forward` or `reverse`) from the
   library preparation: ask the user, it changes the counts. Every other key has a default; list only the ones you
   change. Key summary below.
   Tip: drag a file from Finder into the Terminal window to see its full path.
4. **Run it:** double-click **`RUN_ANALYSIS.command`** in the project folder, or in Terminal in that folder:
   `bash utils/run_rnaseq.sh run --detach`.
   - Before it starts, the launcher checks your inputs (paths, the samplesheet, the output folders, the free disk space
     on the project's disk and on Docker's disk) and checks Docker Desktop's memory against this Mac (see "Resources"
     below). It never changes Docker's settings: if the memory is too low, it says what to set.
   - The analysis runs inside Docker in the background; the window shows its progress. Closing the window does **not**
     stop it. Keep the Mac plugged in with the lid open until it finishes.
   - Check on it any time: `bash utils/run_rnaseq.sh status`.
5. **Results** are in `results/`. Start with `results/qc/multiqc/multiqc_report.html` (the overview of all QC). The counts
   table is `results/counts/<project_name>_counts.txt` (one row per gene, one column per sample). The strandedness check
   is in `results/qc/strandedness/strandedness.tsv` and the alignment rates in `results/qc/hisat2/alignment_summary.tsv`.
   The layout is listed in `ANALYSIS.md`.

**If something fails**, `bash utils/run_rnaseq.sh status` says which step failed and shows the last lines of the log.
The files of the failed step are in `metadata/failed_tasks/<timestamp>/`. Fix the input and run again: the steps that
already finished are reused, not redone.

**Re-running a finished project** only redoes what changed (for example after you change the strandedness, or a trimming
setting), because Nextflow's work folder and its resume cache are kept by default (`keep_work: true`). With the default
`work_volume: true` they are in a Docker volume named after the project, on Docker Desktop's own disk (not in the project
folder; the volume is not visible in Finder). With `work_volume: false` they are `metadata/work/` and `metadata/.nextflow/`
in the project folder. To free disk space, run `bash utils/run_rnaseq.sh clean` (asks first; results are kept) or set
`keep_work: false` in `config.yaml` (deletes the work files after a successful run). The next run then recomputes
everything. A failed or stopped run always keeps them.

**Changing `work_volume`** starts the run from scratch, because the resume state is in the other place: a project whose
`metadata/work/` already has content is asked first (`--yes` answers yes). To resume the earlier run, set
`work_volume: false`. The old `metadata/work/` is left alone (`clean` deletes it).

**Your own analysis** (a table or plot the pipeline doesn't do): put the script in `utils/` (`.R`, `.py` or `.sh`) and run
it with `bash utils/run_rnaseq.sh script utils/my_analysis.R`. It runs in the pipeline's software environment and its log
goes to `metadata/script_logs/`. See "Customizing one project's analysis" below.

## All commands

Run from the project folder (`bash utils/run_rnaseq.sh ...`) or from anywhere after setup (`run_rnaseq.sh ...`, acts on the
current folder):

| Command | |
|---|---|
| `run_rnaseq.sh setup` | one-time: Docker Desktop (offers to install it with Homebrew), Rosetta 2 on Apple chips, the pipeline's Docker image |
| `run_rnaseq.sh init FOLDER` | create a project folder (default: the current folder) |
| `run_rnaseq.sh` (or `run`) | run / resume the analysis of the current folder (`--detach`: start and return, don't follow the log) |
| `run_rnaseq.sh run other_config.yaml` | ... with another config file in the folder |
| `run_rnaseq.sh status` | running? finished? which step failed and why? the work volume's name and size |
| `run_rnaseq.sh logs` | watch the live log (Ctrl+C stops watching, not the analysis) |
| `run_rnaseq.sh stop` | stop the analysis (running again resumes it) |
| `run_rnaseq.sh clean` | delete Nextflow's work folder and resume cache to free disk space: the project's Docker work volume, and `metadata/work/` and `metadata/.nextflow/` if present; asks first (`--yes` skips); the next run starts from scratch, results are kept |
| `run_rnaseq.sh script utils/NAME.R [args]` | run your own script (`.R`, `.py`, `.sh`; inside the project folder) in the pipeline's software environment, in the foreground; log in `metadata/script_logs/` |
| `run_rnaseq.sh help` | the list of commands |

Options: `--yes` answers yes to the launcher's questions (the Homebrew install of Docker Desktop, the Rosetta 2 install,
the deletion by `clean`; for scripts and AI agents); `--pause` waits for Enter before the window closes (used by
`RUN_ANALYSIS.command`). The Mac launcher never changes Docker's settings, so there is no option for that.

## Input table: samplesheet.csv (one row per library)

`samplesheet.csv` is in the project folder root (next to `config.yaml`). Its header is

`sample_id,fastq_1,fastq_2,bam,exclude` followed by any number of metadata columns of your own
(for example `condition,cell_line,treatment`).

| Column | Meaning |
|---|---|
| `sample_id` | Unique library name, used verbatim in every file name. Letters, digits, `.`, `_` and `-` only (`^[A-Za-z0-9._-]+$`). No read suffix (no `.R1`). |
| `fastq_1`, `fastq_2` | Paths of the R1 and R2 FASTQ files (paired-end only). Give both, or give `bam` instead. |
| `bam` | Path of a BAM file. It is converted back to paired FASTQ first and then goes through the normal steps. Give a FASTQ pair or a BAM for each row, not both. |
| `exclude` | `true` drops the row from the whole run; blank or `false` keeps it. Excluded rows stay in the file. |
| metadata columns | Free: copied into `results/samplesheet.tsv` for later analysis. Their names are checked only for being unique. |

Rules the launcher checks: at least one row is kept; each row has a FASTQ pair or a BAM (not both, not neither); every
listed file exists; the header names are unique.

**Comments are not allowed in the CSV**: a `#` line is read as a row and fails the `sample_id` check. Put notes in this
README or in `config.yaml` (which allows `#` comments).

## Config keys (`config.yaml`)

Every key has a default in `nextflow.config`. A project lists only what it changes. `hisat2_index` and `gtf` have no
default. Mac paths are fine (`/Volumes/...`, `/Users/...`); relative paths are relative to the project folder.

**Inputs, outputs and reference**

| Key | Default | Meaning |
|---|---|---|
| `project_name` | `my_rnaseq` | prefix of the counts table (`<project_name>_counts.txt`) and of project-level outputs. Set it. |
| `samplesheet` | `samplesheet.csv` | the library table |
| `fastq_dir` | `null` | not used by this version: name the FASTQs in the samplesheet (they are never moved or copied) |
| `outdir`, `metadata_dir`, `intermediate_dir` | `./results`, `./metadata`, `./intermediate_files` | output folders; **must be relative** to the project (an absolute path such as `/Users/...` is refused) |
| `publish_intermediates` | `false` | Docker versions only. `true` publishes the trimmed FASTQs, the FASTQs made from BAMs, the aligned BAMs and the filtered BAMs to `intermediate_dir`; `false` keeps them only in the work folder (kept while `keep_work: true`), so `intermediate_files/` stays empty |
| `publish_mode` | `auto` | `auto` = hard links from the work folder when the disk allows, else copies (with `work_volume: true`: copies, since hard links cannot cross into the volume); `copy` also accepted |
| `keep_work` | `true` | keep the work folder and resume cache after a successful run (the Docker volume, or `metadata/work` and `metadata/.nextflow` with `work_volume: false`); `false` deletes them after a successful run |
| `work_volume` | `true` | Docker versions only. `true` = Nextflow's work folder and resume cache in a per-project Docker volume on Docker Desktop's disk (fast; not visible in Finder); `false` = `metadata/work` and `metadata/.nextflow` in the project folder. Switching loses the resume (see "Analyze a dataset") |
| `max_cpus`, `max_memory_gb` | `null` | this Mac's budget for the run; see "Resources" |
| `hisat2_index` | `null` | HISAT2 index **prefix**, without `.N.ht2` (for example `/Users/you/REFERENCES/hisat2/t2t/chm13-2`). Required. The launcher stops if `<prefix>.1.ht2` (or `.1.ht2l`) is not there. |
| `gtf` | `null` | gene annotation GTF, used by featureCounts and by the strandedness check (gene model). Required. |

**Strandedness and counting**

| Key | Default | Meaning |
|---|---|---|
| `strandedness` | `unstranded` | the library type given to featureCounts: `unstranded` (`-s 0`, the legacy setting), `forward` (`-s 1`) or `reverse` (`-s 2`). Any other value is an error. The automatic check never changes it. |
| `strand_call_threshold` | `0.8` | the check calls a sample forward or reverse when at least this fraction of the informative reads fits that orientation. Must be in (0.5, 1] |
| `infer_experiment_reads` | `200000` | reads sampled per BAM by RSeQC `infer_experiment.py` (positive integer) |
| `featurecounts_feature_type` | `exon` | GTF feature type counted (featureCounts `-t`) |
| `featurecounts_attribute` | `gene_id` | GTF attribute the counts are summed by (featureCounts `-g`) |
| `featurecounts_args` | `''` | extra featureCounts arguments, appended after the pipeline's own |

**Trimming, alignment and filtering**

| Key | Default | Meaning |
|---|---|---|
| `trim` | `true` | Trimmomatic paired-end trimming. `false`: the reads pass through unchanged |
| `trim_adapter_set` | `TruSeq3-PE-2` | adapter file name (without `.fa`) found in the Trimmomatic tool |
| `trim_illuminaclip` | `'2:30:10:8:True'` | ILLUMINACLIP settings (seed mismatches, palindrome clip threshold, simple clip threshold, minimum adapter length, keep both reads) |
| `trim_steps` | `''` | extra Trimmomatic steps, appended after ILLUMINACLIP (for example `'LEADING:3 TRAILING:3'`). Empty = none (legacy GG) |
| `trim_minlen` | `25` | Trimmomatic `MINLEN` |
| `hisat2_args` | `''` | extra HISAT2 arguments, passed in addition to the pipeline's own (`-p`, `-q`, `--dta`, `--no-discordant`, `--no-mixed`, the index and the summary file). Empty = legacy GG |
| `alignment_min_pct` | `50` | overall HISAT2 alignment rate below which ALIGNMENT_SUMMARY prints a warning. Warns only, never fails |
| `mapq` | `30` | MAPQ filter (samtools `-q`). `0` = no MAPQ filter. Must be a non-negative integer |
| `remove_chroms` | `''` | comma list of contigs to drop before counting. Empty = keep all. Legacy GG removed `chrM`: no effect on the T2T counts, but it zeroes the mitochondrial genes on mm39 (keep this empty for mm39) |

**Duplicates and reports**

| Key | Default | Meaning |
|---|---|---|
| `mark_duplicates` | `true` | Picard MarkDuplicates on the filtered BAM; the marked BAM is the one counted. `false`: the filtered BAM is the final BAM and there are no duplicate metrics |
| `remove_duplicates` | `false` | `true` removes duplicates instead of marking them. Needs `mark_duplicates: true`. Legacy GG: `false` (mark only) |
| `optical_duplicate_pixel_distance` | `2500` | Picard `OPTICAL_DUPLICATE_PIXEL_DISTANCE` |
| `run_fastqc` | `true` | FastQC on the raw and trimmed reads, both mates |
| `run_multiqc` | `true` | the MultiQC report over all QC outputs |

**Differential expression** (off by default): `run_deg` and the `deg_*` keys are described in "Differential expression" below.

**Script replacement**: each pipeline script can be replaced for one project by a key named after it (see "Customizing
one project's analysis"). `run_rnaseq.sh` reads them itself.

## Customizing one project's analysis

Every step runs one script from `bin/`. To change a step for one project only, copy its script from this folder's `bin/`
into the project's `utils/`, edit it, and add a `config.yaml` key named after the script (file name without extension):

```yaml
featurecounts: utils/featurecounts.sh   # replace a step's script for this project only
featurecounts_args: "..."               # or pass extra arguments to the shipped script
```

Each script's header documents its arguments and outputs; a replacement must keep that contract. Changing a script key
re-runs that step and the steps after it. The keys are the base names of the scripts in `bin/` (`bam_to_fastq`, `fastqc`,
`trim`, `align`, `alignment_summary`, `bam_filter`, `mark_duplicates`, `gtf_to_genepred`, `genepred_to_bed12`,
`infer_strandedness`, `qc_strandedness_table`, `featurecounts`, `counts_merge`, `samplesheet_tsv`, `multiqc`, and the `deg_*` scripts of "Differential expression"). Some keys
are both a switch and a script name: `trim` and `mark_duplicates`. A `true` or `false` value is the switch; a path is a
replacement script. Analyses that the pipeline does not do at all are better written as your own script in `utils/`, run
with `run_rnaseq.sh script` (see "Analyze a dataset" above).

## Inputs are never copied

- FASTQs and BAMs stay where they are. The pipeline reads them in place, never moves, copies, renames or deletes them,
  and `config.yaml` and the samplesheet keep pointing to them. Keep the path stable (if you move the files, update the
  samplesheet and `config.yaml`).
- Write **Mac paths** in full (`/Volumes/SeqDrive/FASTQ/run7/x_R1.fastq.gz`, or `~/FASTQ/run7/x_R1.fastq.gz`). Relative paths
  are relative to the project folder (`original_data/fastq/x_R1.fastq.gz`).
- The top-level folder that holds a path named in `config.yaml` or in the samplesheet is mounted **read-only** in the
  analysis, at the same path (`/Users`, or `/Volumes/<name>` for a drive or a share mounted in Finder). Docker Desktop must
  be allowed to share that folder. By default it shares `/Users`, `/Volumes`, `/private`, `/tmp` and `/var/folders`; for
  any other folder add it in Docker Desktop > Settings > Resources > File sharing, otherwise the run stops with a message.
- **Network paths** (`smb://...`, `//server/share/...`) are refused: Docker Desktop cannot read network shares. Mount the
  share in Finder (it then appears under `/Volumes/<name>`) and write the path that starts with `/Volumes/<name>/`.
- The launcher records the paths as you wrote them and as the analysis used them in `metadata/host_paths.tsv`, and logs the
  mapping at the start of the run.

## Resources: how much of the Mac an analysis uses

Each run has one **budget** of memory and CPU cores, and never goes over it. Docker enforces it as a hard limit, so a run
cannot freeze the Mac by running out of memory.

**Memory.** The budget is the smallest of:
- 75% of the Mac's memory, but at least 3 GB left for macOS;
- Docker Desktop's memory minus 0.5 GB (the Docker VM's own share);
- the most this pipeline can use (`max_useful_memory_gb` in `docker/requirements.conf`).

**CPU cores.** All cores but one (all but two from 8 cores up), so the Mac stays responsive.

**The largest step** is ALIGN (HISAT2 with the genome index in memory); its memory request is set per step in
`conf/docker_tiers.config` for the Mac's tier, and it is capped to the budget. The minimum budget the run needs is
`min_memory_gb` in `docker/requirements.conf`. Both figures there are placeholders until a pilot run on a real Mac
confirms them. The run refuses to start, with an explanation, if the budget is below the minimum.

`config.yaml` can lower the budget: `max_memory_gb` (decimals allowed, e.g. `2.5`) and `max_cpus`. Neither can go above
what Docker Desktop allows (the launcher warns). Results do not change, only speed. Each run prints its budget at the start.

**Docker Desktop's memory setting (recommended, not set by the launcher).** Docker Desktop caps the memory all its
containers may use: Docker Desktop > Settings > Resources > Memory. The launcher never changes it. It recommends a value
by the Mac's total memory:

| Mac memory | Docker Desktop memory (recommended) |
|---|---|
| less than 12 GB | 6 GB |
| 12 to less than 24 GB | 8 GB |
| 24 GB or more | 16 GB |

Set Memory to that value and click *Apply & restart*; wait for "Engine running" again. `setup` and each run print the
recommendation when Docker's memory is holding the analysis back. If Docker's memory is below what the analysis needs
(the pipeline's minimum, `min_memory_gb` in `docker/requirements.conf`), the run refuses to start and says what to set.
That figure, like `max_useful_memory_gb`, is a placeholder until a pilot run on a real Mac confirms it.


**Keep memory free.** The pipeline caps its own memory, but other programs share it. Close browsers with many tabs, large
IDE sessions and other analyses before a long run. At the start of a run the launcher warns when the Mac's free memory is
below the run's budget.

## Speed on a desktop

- **Drives.** Put the project folder and the FASTQ files on the internal SSD when you can. With `work_volume: true` the
  work files go to Docker's disk (where Docker keeps its data); with `work_volume: false` they go to the project
  folder. Choose the faster or bigger disk for the work files.
- **Hard disk.** This launcher does not detect hard disks. On a USB or other hard disk, `work_volume: true` with the
  work files on the internal SSD is the faster choice. With `work_volume: false` the disk-heavy steps (FASTQC, TRIM,
  BAM_FILTER, MARK_DUPLICATES) run at most 2 at a time, because several streams at once make a hard disk slower.
- **CPU and memory.** More cores (`max_cpus`) let more samples run at once. The steps use 1 to 8 cores each (FASTQC 1,
  TRIM 2, ALIGN up to 8, MARK_DUPLICATES 2). Each step requests what it measured needs (ALIGN 8.5 GB, TRIM 2 GB,
  MARK_DUPLICATES 8 GB), so a budget of 22 GB lets two alignments run together.
- **Expected time.** 12 libraries of about 4.7 GB gzipped each (56 GB of FASTQ) take about 5 hours
  (measured on a 20-core Windows desktop with 22 GB for Docker, FASTQs and index on its SSD: 4 h 50 min). Trimmomatic
  uses about 1.6 cores whatever it is given, so each trim takes about 30 min here against 12-17 min on the cluster.

## Disk space

- **Plan for the work files on Docker's disk and the results on the project's disk.** With the default
  `work_volume: true`, the work files (the working files of every step, kept so a re-run resumes) are on the disk that
  holds Docker Desktop's disk image, and the results are on the project's disk. The launcher computes both before every
  run from the total size of the FASTQ (or BAM) files named in the samplesheet, using the factors in
  `docker/requirements.conf` (`disk_factor`, `disk_results_factor`, `disk_hard_min_gb`). These factors are estimates
  pending the next run. It warns when the free space is below the figure, and stops before the run only when less than
  `disk_hard_min_gb` (20 GB) is free on the project's disk.
- **Cluster benchmark (12 libraries, 57 GB of FASTQ):** the final BAMs were about 32 GB, and the counts and QC tables are
  small. The work files of the Docker version were not measured yet.
- **To estimate:** in Finder, select the FASTQ folder of your samplesheet and press ⌘I to read its size.
- What uses the space: the work volume on Docker's disk, the trimmed FASTQs, FASTQs made from BAMs, aligned and filtered
  BAMs in `intermediate_files/` (only with `publish_intermediates: true`; the default on Mac is `false`, so the folder
  stays empty), and the results in `results/`. The raw FASTQs are not copied. With the work volume, files in `results/`
  are copies. With `work_volume: false` they are hard links to the work folder when the disk allows (the run log says
  so). Never edit a published file in place; replacing or deleting it is fine.
- To free space: `run_rnaseq.sh clean` (results are kept), or `keep_work: false` (deletes the work files after a
  successful run).
- **Keep the data out of iCloud-synced folders.** Do not keep the FASTQs, the project or the work files in Desktop or
  Documents when *Optimize Mac Storage* is on: macOS may move files to the cloud.
- **The project can be on another disk.** With the work volume, only the results (and `intermediate_files/` with
  `publish_intermediates: true`) are written to the project folder. Docker's disk image is what matters most for speed:
  keep it on the internal SSD (see "Good to know").

## Differential expression

The pipeline can also run the DESeq2 differential expression analysis of the legacy `Scripts/step_6_deg_analysis.Rmd`. It is **off by
default** (`run_deg: false`): with the default nothing changes and the pipeline ends at the counts table. Everything below is
about `run_deg: true`.

### Quick start

```yaml
# config.yaml, on top of the usual keys
run_deg: true
deg_comparisons: deg_comparisons.csv      # required: what to compare (format below)
deg_heatmaps: deg_heatmaps.csv            # optional: heatmaps of the DEGs
deg_colors: deg_colors.csv                # optional: group colours
deg_label_rules: deg_label_rules.csv      # optional: nicer sample names in the heatmaps
deg_reference_levels: treatment=LUC_DMSO  # the reference level of each design factor
```

The sample metadata (the columns the designs use, for example `treatment`) comes from `deg_metadata`, or from `samplesheet.csv` when
that is not set (rows with `exclude: true` are dropped). The count table is this run's `results/counts/<project_name>_counts.txt`.
The analysis starts after COUNTS_MERGE, so the reads are aligned and counted first.

### From an existing counts table

To run only the DEG steps on a counts table you already have, set `deg_counts` (and `run_deg: true`):

```yaml
run_deg: true
deg_counts: /Users/you/RNASEQ/earlier_run/results/counts/my_project_counts.txt   # Gene + one column per sample
deg_metadata: deg_metadata.csv          # or leave it out to use samplesheet.csv (it then needs only sample_id + the metadata columns)
deg_comparisons: deg_comparisons.csv
```

No FASTQ, BAM, `hisat2_index`, `gtf`, QC or MultiQC is needed or checked, and none of the upstream steps runs. `deg_counts` without
`run_deg: true` is a startup error. The column names of the counts table must equal the sample names in the metadata; a suffix such
as the legacy `.R1` is removed with `deg_strip_suffix: '\.R1$'`.

### What runs

| Step | Runs | What it does |
|---|---|---|
| DEG_PREPARE | once | Checks the counts table, the metadata and the comparisons table, aligns the samples, sets the reference levels, writes the cleaned copies. |
| DEG_QC | once (`deg_run_qc`) | Reads per sample (bar plot), MDS plot (`pca.pdf`, it is an MDS plot), sample correlation heatmap. |
| DEG_DESEQ2 | per model | One DESeq2 fit per model: all rows of `deg_comparisons` with the same `model` share it. |
| DEG_RESULTS | per model | A results table for every comparison of the model, shrunk fold changes, an analysis summary, normalised counts. |
| DEG_TRANSFORM | per model that needs it | rlog (or vst) counts for the heatmaps and gene plots. Slow, so only models used by a heatmap (and the first model for gene plots) are transformed. |
| DEG_VOLCANO | per model (`deg_run_volcano`) | Volcano plots: one SVG per comparison (and a shrunk one), three multi-page PDFs, and a combined figure. |
| DEG_MA | per model (`deg_run_ma`, off) | MA plots. |
| DEG_HEATMAP | per row of `deg_heatmaps` | Heatmap of the DEGs of one comparison. |
| DEG_GENE_PLOT | once (`deg_gene_plot_genes`) | One expression plot per gene, from the first model's transformed counts. |

### deg_comparisons.csv (required)

One row per comparison. The header names are exact; the order of the columns is free.

| Column | Required | Meaning |
|---|---|---|
| `model` | yes | Name of the model and of its output folder (`results/deg/<model>/`). Letters, digits, `.`, `_`, `+`, `-`. `prepared`, `qc` and `gene_plots` are reserved. Rows with the same `model` share one DESeq2 fit. |
| `design` | yes | The DESeq2 design formula, for example `~treatment` or `~ treatment + batch`. Every row of one model must have the same design. |
| `comparison` | yes | Name of the comparison (unique within its model; same characters as `model`; `log` is reserved). It names the folder and the files `results_<comparison>.csv`. |
| `factor`, `numerator`, `denominator` | one contrast kind | A contrast between two levels of a design factor, for example `treatment`, `CBX2_1uM`, `LUC_DMSO`: log2 fold change of numerator over denominator. All three together. |
| `coef` | one contrast kind | The name of one model coefficient (as `resultsNames(dds)` lists it), for example `drug_1uM_vs_DMSO` (design `~ drug`). |
| `add`, `subtract` | one contrast kind | A sum of coefficients: names separated by `;`, for example `treatment_B_vs_A;treatment_C_vs_A` (design `~ treatment`, sum of two coefficients). `subtract` (optional) lists the coefficients that are subtracted. |
| `label` | no | Plot title of the comparison. Blank = the comparison name with `_` replaced by spaces. |
| `combined_order`, `combined_title` | no | A number puts the comparison in the combined volcano figure, in that order, with that title. The figure is made when at least two comparisons of a model have a number. |
| `xlim`, `ylim` | no | Volcano axis limits of this comparison as `min,max` (put the cell in quotes: `"-3,3"`). They win over `deg_volcano_xlim` / `deg_volcano_ylim`. |

Each row needs exactly one contrast kind (a row with two kinds, or none, stops the run at startup). `deg_shrink_type: apeglm` works with
`coef` rows only.

```
model,design,comparison,factor,numerator,denominator,coef,add,subtract,label,combined_order,combined_title,xlim,ylim
treatment,~treatment,CBX2_1uM_vs_LUC_DMSO,treatment,CBX2_1uM,LUC_DMSO,,,,CBX2 KO + EZH2i,3,CBX2 KO + EZH2i,,
```

### deg_heatmaps.csv (optional)

One row per heatmap. A heatmap shows the genes that are differentially expressed (not `NS`) in one comparison, from the model's
transformed counts, scaled by row.

| Column | Required | Meaning |
|---|---|---|
| `name` | yes | Unique name; the file is `results/deg/<model>/heatmap_<name>.svg`. |
| `model`, `comparison` | yes | Which results table to use. The pair must exist in `deg_comparisons`. |
| `gene_set` | no | A text file with one gene per line: only DEGs in the set are drawn. |
| `clust_method` | no | `hclust` method; blank = `ward.D2`. |
| `cutree_rows` | no | Number of row clusters to cut; blank = none. |
| `col_swap` | no | `none` (default), `scaled` or `unscaled`: swap the two main column branches (clustering on the scaled or unscaled matrix). |
| `palette_n` | no | Number of colours; blank = 1000. |
| `fontsize`, `fontsize_row`, `width`, `height` | no | Blank = 9, 6.5, 6.5, 4.5 (inches). |
| `title` | no | `{n}` is replaced by the number of genes. Blank = `DEGs: <label> (N = {n})`. |

### deg_colors.csv and deg_label_rules.csv (optional)

`deg_colors.csv` has the columns `group,color` (a group is a level of the colour-by column, a colour a name or a hex code). A group
without a row gets a default colour.

`deg_label_rules.csv` has the columns `pattern,replacement`. The rules are applied in order, as regular expressions, to the sample names
shown in the heatmaps (the same as the legacy chain of `gsub()` calls). A typed `\n` in the replacement is a line break. Quote
cells with leading or trailing spaces, and write an empty replacement as `""`:

```
pattern,replacement
1uM,\nEZH2i
A_," A"
DIPGVI_,""
```

### Parameters

All default to the legacy Rmd (except the colours). The same list with comments is in `examples/config.yaml`. Hex colours need quotes
in YAML (`'#46627A'`): an unquoted `#` starts a comment.

| Key | Default | Meaning |
|---|---|---|
| `run_deg` | `false` | Master switch. |
| `deg_counts` | `null` | Existing counts table: skips all upstream steps (needs `run_deg: true`). |
| `deg_metadata` | `null` | Sample metadata CSV/TSV; `null` = the samplesheet. Must contain the `deg_sample_col` column. |
| `deg_comparisons` | `null` | Comparisons CSV; required with `run_deg`. |
| `deg_heatmaps`, `deg_colors`, `deg_label_rules` | `null` | Optional CSVs (above). |
| `deg_sample_col` | `sample_id` | Metadata column with the sample names. |
| `deg_strip_suffix` | `''` | Regular expression removed from the count column names (and the metadata sample names), for example `'\.R1$'`. |
| `deg_reference_levels` | `''` | `factor=level,factor=level`. Factors not listed keep the alphabetical first level (the run prints which). |
| `deg_min_count`, `deg_min_samples` | `10`, `2` | Keep genes with at least `deg_min_count` reads in at least `deg_min_samples` samples. Positive integers. |
| `deg_fit_type` | `parametric` | DESeq2 dispersion fit: `parametric`, `local` or `mean`. Used by DEG_DESEQ2 and DEG_TRANSFORM. |
| `deg_alpha` | `0.05` | FDR cut-off, between 0 and 1 (not included). |
| `deg_lfc` | `1` | Absolute log2 fold change cut-off, 0 or more. |
| `deg_lfc_threshold_in_results` | `false` | `true` passes `lfcThreshold` to `results()`. |
| `deg_shrink_type` | `ashr` | `ashr`, `apeglm` (coef rows only) or `none`. |
| `deg_sval_threshold` | `0.005` | s-value cut-off of the shrunk effect. |
| `deg_transform_method`, `deg_transform_blind` | `rlog`, `false` | Transformation for heatmaps and gene plots (`rlog` or `vst`). |
| `deg_color_down`, `deg_color_up`, `deg_color_ns`, `deg_color_mid` | `'#46627A'`, `'#8c6bb1'`, `gray`, `white` | Volcano / MA / heatmap colours. |
| `deg_color_by` | `''` | Metadata column that colours the QC plots; blank = the first design variable of the first model. |
| `deg_mds_input` | `logcpm` | `logcpm` (log2 CPM, prior count 2) or `counts` (raw counts as the legacy MDS, dominated by library size). |
| `deg_mds_label_by`, `deg_mds_legend_title` | `''` | Point labels of the MDS plot; legend title (blank = the colour-by column). |
| `deg_legacy_mds_legend_col` | `''` | A metadata column name reproduces the legacy MDS legend. |
| `deg_cor_method` | `spearman` | `pearson`, `spearman` or `kendall`. The heatmap file is `<method>_correlation_heatmap.svg`. |
| `deg_label_genes` | `''` | Comma list of genes named in the volcano plots. |
| `deg_volcano_xlim`, `deg_volcano_ylim` | `''` | `min,max` of every volcano; blank = automatic. |
| `deg_volcano_common_xlim`, `deg_volcano_common_ylim` | `''` | Limits of the combined figure and the scaled PDFs. |
| `deg_volcano_combined_layout` | `stacked` | `stacked` or `side`. |
| `deg_gene_plot_genes` | `''` | Comma list of genes to plot; blank = no gene plots. Genes missing from the table are skipped. |
| `deg_gene_plot_x` | `''` | x axis and colour column of the gene plots; blank = `deg_color_by`, else the first design variable. |
| `deg_legacy_summary_counts` | `false` | `true` reproduces the legacy analysis summary (always 0 up / 0 down). |
| `deg_run_qc`, `deg_run_barplot`, `deg_run_mds`, `deg_run_correlation` | `true` | Switch the QC step and its three plots. |
| `deg_run_volcano`, `deg_run_volcano_pdfs`, `deg_run_volcano_combined` | `true` | Volcano plots, the three PDFs, the combined figure. |
| `deg_run_heatmap` | `true` | Heatmaps (also need `deg_heatmaps`). |
| `deg_run_ma` | `false` | MA plots. |
| `deg_run_gene_plots` | `true` | Gene plots (run only when `deg_gene_plot_genes` is set). |

One more legacy switch of the R scripts is not a config key; pass it with the `_args` key: `deg_qc_args: "--legacy-color-index=TRUE"`
(reproduce the legacy QC colours).

### Custom statistical tests (not shipped)

The pipeline runs per-model Wald tests (design, `coef` and `add` / `subtract` contrasts). An interaction model, or any other custom
statistical test, is not part of the pipeline: write it as a custom script in the project's `utils/` folder. It can read
`results/deg/prepared/counts_edited.rds` and `results/deg/prepared/sample_metadata_edited.rds`.
Run such a script with `run_rnaseq script utils/X.R`.

### Outputs

```
results/deg/
├── prepared/        counts_edited.rds/.tsv, sample_metadata_edited.rds/.csv, comparisons.tsv, prepare_summary.txt
├── qc/              counts_per_index.pdf, counts_per_sample.tsv, pca.pdf (MDS), <method>_correlation_heatmap.svg
├── <model>/         DESeqDataSet.rds, DESeq.rds, design.txt, dispersion_plots.svg, prefilter.tsv,
│                    analysis_summary.csv/.txt, normalized_counts.csv, rlogData.rds, rlogCounts.rds,
│                    volcano_plots_{yscaled,xscaled,freescale}.pdf, volcano_combined_<n>_comparisons_<layout>.svg,
│                    heatmap_<name>.svg, log/ (command.sh of every step, <step>_command.sh)
│   └── <comparison>/  results_<comparison>.csv, volcano_all_genes_<comparison>.svg, volcano_shrunk_all_genes_<comparison>.svg, ma_<comparison>.svg
└── gene_plots/      <gene>.svg
```

`results_<comparison>.csv` has the columns `Gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, shrunkLFC, svalue, effect,
shrunk_effect` (the shrunk columns are missing with `deg_shrink_type: none`); `effect` is `Upregulated`, `Downregulated` or `NS`.

The DEG steps run in the same omics image as the other steps (it contains R with DESeq2).

### Re-running only the plots

The work folder and the resume cache are kept after every successful run (`keep_work: true`, the default). Only `clean` removes them (or `keep_work: false`, which removes them after a successful run)
(`bash utils/run_rnaseq.sh clean`, which keeps the results). To change a figure, edit `config.yaml` or a CSV next to it (for example
`deg_comparisons.csv`, `deg_heatmaps.csv`, or a `deg_*_args` value such as `deg_volcano_args`), then run `bash utils/run_rnaseq.sh` again.
The run resumes from the work folder and reruns only the plot steps that changed; the steps that did not change are reused.

When the work folder is gone (after `clean`, with `keep_work: false`, or on another computer), use the `deg_counts` mode instead:
`run_deg: true` with `deg_counts` set to a counts table. Only the DEG steps run then (see "From an existing counts table" above).

## Strandedness

`strandedness` sets how featureCounts counts: `unstranded` = `-s 0`, `forward` = `-s 1`, `reverse` = `-s 2`. The default
`unstranded` is what the legacy GG counts used.

The automatic check runs RSeQC `infer_experiment.py` on `infer_experiment_reads` reads of each final BAM, against a gene
model made from the GTF. `results/qc/strandedness/strandedness.tsv` reports, per sample, the fraction of informative reads
explained by each strand orientation, the call (`forward`, `reverse` or `unstranded`, by `strand_call_threshold`), the
configured value and a `mismatch` column (`true` when the call differs from the configured value). A mismatch prints a
warning at the end of the run. **The check never changes the counting.** Only the `strandedness` key does. Ask the user
which library preparation was used before you change it.

## Checking a run against the benchmark

The lab keeps a benchmark of this pipeline, outside this folder (ask the user or the lead where the copy is). It has a
quick test (the `sub` set: the first 1,000,000 read pairs of four libraries, plus one library that starts from a BAM;
minutes) and the full GG_EZH2i_RNAseq run (the `full` set: 12 libraries; hours; 57 GB of FASTQ). Use them to check a new
installation or a new image before you trust a dataset:

1. Make a project for the `sub` test with its config and samplesheet, with the paths mapped to the local copies (the
   benchmark's config files carry the cluster's paths; change them to yours).
2. Run it, then compare the project with the benchmark, from this folder in Terminal:
   ```
   bash test/compare_benchmark.sh PROJECT_DIR BENCHMARK_DIR
   ```
   where `PROJECT_DIR` is the project folder and `BENCHMARK_DIR` the benchmark folder for the same test.
3. What must match:
   - the counts table `results/counts/<project_name>_counts.txt` must be **byte-identical** to the benchmark's (the
     benchmark's own counts table already has the legacy `.R1` suffix removed, so compare the files as they are);
   - the other text outputs (gene lengths, per-sample counts, featureCounts summaries, the HISAT2, Trimmomatic, Picard,
     samtools and strandedness QC tables) are listed with their md5 in the benchmark's `benchmark_text_md5.txt`; a
     mismatch in a QC table is compared cell by cell before it is called a difference;
   - Picard metrics and the FastQC and MultiQC files carry dates, paths and versions: compare their numbers, not their md5;
   - the BAMs are not in the benchmark: compare them with `samtools flagstat` and `idxstats`, and `samtools view` piped to
     `md5sum`. HISAT2 output order can depend on the thread count without changing the counts.
4. A difference is explained before it is called a bug: first check the image's tool versions (`setup` output), which may
   differ from the cluster's pins (see `ANALYSIS.md`), and the paths and settings in the config.

## Differences from the cluster version

The steps, the tools and the results are those of the cluster version. Four things differ here:

- **One image holds all tools.** There is no per-step container and no Apptainer; every step runs in the same Docker
  image, on the Mac, through a local executor. The image is an Intel-chip (linux/amd64) image: on an Apple-chip Mac it
  runs through Rosetta 2. The image's tool versions may differ from the cluster's pins (see "Software versions").
- **No cluster scheduler.** The cluster's job submission, its queue, its host exclusion setting and its job watchdog do
  not exist here, so there are no `lsf_*` keys.
- **Intermediates are not published by default.** `publish_intermediates` (default `false`) decides whether the trimmed
  FASTQs, the FASTQs made from BAMs and the aligned and filtered BAMs are copied to `intermediate_files/`.
- **Inputs are mounted, not moved.** The FASTQs and BAMs are read in place, read-only (see "Inputs are never copied").

## Software versions and default settings

The image `ghcr.io/genesig/omics` (tag in `docker/image.txt`) contains the tools. `run_rnaseq.sh setup` prints the main
versions; cite those, and the tool headers in the logs, for a run. The cluster version's pins are:

| Tool | Cluster pin |
|---|---|
| HISAT2 | 2.2.1 |
| samtools | 1.23 |
| Trimmomatic | 0.40 |
| FastQC | 0.12.1 |
| Picard | 2.27.5 (MarkDuplicates) |
| featureCounts (Subread) | 2.1.1 |
| RSeQC (`infer_experiment.py`) | 5.0.6 |
| UCSC gtfToGenePred, genePredToBed | 482 |
| MultiQC | 1.35 |
| Nextflow | 26.04.6 |

The image's versions may differ from these pins; do not claim the pins as the versions of a Docker run.

Defaults that matter for the counts: `unstranded` (the legacy setting), `mapq: 30`, duplicates marked and not removed,
chrM kept (`remove_chroms: ''`), Trimmomatic `TruSeq3-PE-2` with `MINLEN 25`, featureCounts on `exon` features grouped by
`gene_id`.

## Good to know

- **Apple-chip Mac computers** run the pipeline's tools through Rosetta 2. Keep Docker Desktop > Settings > General >
  *Use Rosetta for x86_64/amd64 emulation on Apple Silicon* on (the default). Without it the pipeline still runs, much
  more slowly.
- **Where Docker keeps its data:** Docker Desktop stores its disk image on the Mac's internal disk by default. The work
  volume lives on the same disk, so it needs room for the work files and is best on the internal SSD. If the disk is
  short of room, move it in Docker Desktop > Settings > Resources > Advanced > *Disk image location* (the launcher assumes
  no particular disk).
- **Reproducibility:** every tool runs at a pinned version inside the image, and the image tag is in `docker/image.txt`.
  Every step's `log/` folder holds the exact command that ran (`<sample_id>_<step>_command.sh`). `metadata/pipeline_info/`
  holds the Nextflow execution report, trace and timeline. A resumed run keeps the previous run's trace, report and
  timeline next to it, named with that run's timestamp.
- **Updating the pipeline:** `git pull` in the `RNASEQ/MAC` folder of the clone. If the update changed the image tag, the
  next `setup` or `run` downloads the new image.
- **Example project** in `examples/projects/gg_ezh2i_rnaseq/` (the GG_EZH2i_RNAseq dataset, 12 libraries) uses `/Volumes`
  and `/Users` placeholder paths; copy it and edit its paths before you use it.

For what each step does, which numbers the counts table and the QC tables contain, and how to read them, see `ANALYSIS.md`
in this folder. The analysis is the same in all three versions (cluster, Windows and Mac); only how it is launched differs.
