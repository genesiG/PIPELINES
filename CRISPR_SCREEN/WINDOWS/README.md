# CRISPR_SCREEN for Windows (Docker)

The Phillips Lab's pooled CRISPR screen pipeline — demultiplex → FastQC/MultiQC →
`mageck count` → guide-representation QC (tables, plots, an HTML report per
library) → optional `mageck test` per comparison → a ready-to-adapt
`METHODS.md` — running on **your own Windows computer**. It is the same
pipeline as on the PMACS cluster (same tools, same versions, same outputs), packaged in one Docker image, so you don't install any bioinformatics software yourself.

> **Easiest way:** ask Claude (or another AI agent) to do it for you, e.g.
> *"Clone the lab's CRISPR_SCREEN pipeline from GitHub and analyze my screen:
> the FASTQ is in D:\Sequencing\run7, the samplesheet and library are attached."*
> The agent follows `AGENTS.md` in this folder. Everything below is what it
> does, in case you prefer to do it yourself.
>
> **First time on this computer?** Follow `SETUP.md` in this folder once (Docker Desktop,
> WSL 2, where Docker's data lives, memory). Then come back here.

## What you need

- Windows 10 (22H2) or 11, 64-bit, with hardware virtualization enabled
  (on by default on most computers).
- **6 GB of memory or more** (8 GB or more lets more steps run side by side; the pipeline itself needs about 2.5 GB).
- Free disk space: about **2× the size of your raw FASTQ files** on the project's
  drive (per-sample FASTQs and Nextflow's working folder `metadata\work\`), plus
  space for results and ~8 GB once for the Docker image. The raw FASTQs are not copied.
- Internet the first time (to download Docker and build the image).

## One-time setup (per computer)

1. Get the pipeline (once), as a sparse clone so that only this folder is downloaded:
   ```
   git clone --filter=blob:none --sparse <the lab's GitHub URL>
   cd CRISPR_SCREEN
   git sparse-checkout set WINDOWS
   ```
   Then work in the `WINDOWS` folder; `git pull` later updates it.
2. Open the `CRISPR_SCREEN\WINDOWS` folder and run, in PowerShell or Command Prompt:
   ```
   .\run_crisprscreen.bat setup
   ```
   - If Docker Desktop isn't installed, it offers to install it (Windows asks
     for administrator permission). Then **restart the computer**, open
     **Docker Desktop** from the Start menu once, accept its terms (no account
     needed — click *Skip*), and run `setup` again.
   - It then builds the pipeline's Docker image: ~3 GB download, 10–30 minutes,
     once.
   - It adds `run_crisprscreen` to your PATH, so new terminal windows can use it
     from any folder.

## Analyze a screen

1. **Create a project folder** (one per screen, *outside* the CRISPR_SCREEN folder):
   ```
   run_crisprscreen init C:\Users\<you>\Documents\CRISPR_screens\my_screen
   ```
   It creates the folder skeleton, `config.yaml`, `examples\` (starter CSVs),
   `AGENTS.md`/`CLAUDE.md` (instructions for an AI agent), `RUN_ANALYSIS.bat`
   (double-click) and `utils\run_crisprscreen.bat` (the command, used below).
   It also copies the default barcode-position plot script into `utils\` (only
   if absent) and the agent definitions into `.claude\` (missing files only),
   and both go into the project's first git commit. It doesn't run anything.
2. **Add your tables** to the project folder (copy them from `examples\`;
   edit in Excel and *Save As → CSV UTF-8*):
   - `samplesheet.csv` — one row per sample: `barcode_id,barcode,sample_name`
     (barcodes only when the pipeline demultiplexes; optional `library` column
     when samples come from different library pools).
   - `library.csv` — your sgRNA library: `gene,gRNA_ID,sequence`.
   - `contrasts.csv` (optional — when you're ready for `mageck test`):
     `treatment,control`, one row per comparison.
3. **Edit `config.yaml`** (Notepad is fine): set `project_name`, and either
   - `fastq_r1:` the multiplexed FASTQ (e.g. `D:\Sequencing\run7\Undetermined_S0_L001_R1_001.fastq.gz`), or
   - `fastq_dir:` a folder of already-demultiplexed FASTQs, one per sample
     (`<sample_name>.fastq.gz`).

   The FASTQs **stay where they are**. The analysis reads them from that path
   and never moves, copies, renames or deletes them; `config.yaml` keeps
   pointing to them. Keep that path stable (if you move the files, update
   `config.yaml`). Several projects can point at the same raw file. Use a local
   or USB drive: UNC paths (`\\server\share`) and mapped network drives are
   refused, because Docker Desktop cannot mount them (copy the files to a local
   or USB drive). The path must not contain a comma. Files outside the project
   folder are mounted read-only into Docker. A FASTQ that is a shortcut must be
   given by its real path.
4. **Run it:** double-click **`RUN_ANALYSIS.bat`** in the project folder, or in
   a terminal in that folder: `.\utils\run_crisprscreen.bat`.
   - The analysis runs inside Docker in the background; the window shows its
     progress. Closing the window does **not** stop it — but keep the computer
     on and awake (plugged in, sleep off) until it finishes.
   - Check on it any time: `.\utils\run_crisprscreen.bat status`.
5. **Results** are in `results\` — start with
   `results\counts\qc\report\*_qc_report.html` (open in a browser; every plot
   is explained) and `results\demultiplexing\qc\multiqc\multiqc_report.html`.
   `METHODS.md` describes exactly what was done, for your paper. Every plot is
   saved twice, as SVG and as a PNG with the same name (handy for slides); with
   more than one library, `results\counts\qc\plots\pools_comparison\` also has
   plots comparing the libraries.

**If something fails**, `.\utils\run_crisprscreen.bat status` says which step and
why (details in `metadata\failed_tasks\`). Fix the input and run again: steps
that already finished are reused, not redone.

**Re-running a finished project** only redoes what changed (e.g. after adding
contrasts, or tweaking a plot: just that step and the report), because
Nextflow's working folder `metadata\work\` and its resume cache
`metadata\.nextflow\` are kept by default (`keep_work: true`). To free the disk
space instead, run `.\utils\run_crisprscreen.bat clean` (asks first; results are
kept) or set `keep_work: false` in `config.yaml` (deletes them after a successful
run); the next run then recomputes everything (about 10-20 min for a full lane).
A failed or stopped run always keeps them. If you see a `work\` or `.nextflow\`
folder directly in the project folder, it is from an older version of the
pipeline and is no longer used; `clean` removes it.

**Your own analysis** (a plot or test the pipeline doesn't do): put the script in
`utils\` (`.R`, `.py` or `.sh`) and run it with
`.\utils\run_crisprscreen.bat script utils\my_analysis.R`. It runs in the
pipeline's software environment and its log goes to `metadata\script_logs\`. See
"Analyses beyond the pipeline" in the project's `AGENTS.md`.

## All commands

Run from the project folder (`.\utils\run_crisprscreen.bat ...`) or from anywhere
after setup (`run_crisprscreen ...`, acts on the current folder):

| Command | |
|---|---|
| `run_crisprscreen setup` | one-time: Docker Desktop + the pipeline's Docker image |
| `run_crisprscreen init FOLDER` | create a project folder |
| `run_crisprscreen` (or `run`) | run / resume the analysis of the current folder (`--detach`: don't watch) |
| `run_crisprscreen run other_config.yaml` | ... with another config file |
| `run_crisprscreen status` | running? finished? which step failed and why? |
| `run_crisprscreen logs` | watch the live log (Ctrl+C stops watching, not the analysis) |
| `run_crisprscreen stop` | stop the analysis (running again resumes it) |
| `run_crisprscreen clean` | delete Nextflow's working folder `metadata\work\` and resume cache `metadata\.nextflow\` (and any old root `work\` / `.nextflow\`) to free disk space; asks first (`--yes` skips); the next run starts from scratch, results are kept |
| `run_crisprscreen script utils\NAME.R [args]` | run your own script (`.R`, `.py`, `.sh`; inside the project folder) in the pipeline's software environment; log in `metadata\script_logs\` |

## Resources: how much of the computer an analysis uses

Every run gets one **budget** and can never exceed it — Docker enforces it as
a hard limit, so the analysis cannot freeze Windows by running it out of memory:

| | default | why |
|---|---|---|
| memory | **what the pipeline can use (8 GB), but never more than 75% of the RAM, always leaving 3 GB to Windows** | more would not run faster (below); a plain percentage starves small laptops (80% of 8 GB leaves 1.6 GB) |
| CPU cores | **all but 1** (all but 2 from 8 cores up) | the desktop stays responsive; heavy CPU use only slows a computer, it can't harm it |

…and never more than Docker Desktop itself may use. Examples: 8 GB / 4 cores →
5 GB, 3 cores; 16 GB / 8 cores → 8 GB, 6 cores; 32 GB / 12 cores → 8 GB, 10 cores.
Each run prints its budget at the start.

- **Why not more memory on a big computer:** every step of this pipeline
  needs less than 1 GB. Extra memory would only fill with Docker's file
  cache (the FASTQs it has read), which Windows gets back slowly: on a
  computer whose memory is already busy with other programs, that is what
  makes Windows stop responding. While an analysis (or an image build) runs,
  `run_crisprscreen` also watches Windows' free memory and makes Docker drop
  that cache when it gets low (`docker\requirements.conf`, `max_useful_memory_gb`).
- **Docker Desktop's own limit:** with WSL 2 it may use only half the
  computer's memory by default, which is enough unless the computer has less
  than about 18 GB. `run_crisprscreen setup` then offers to raise it (it writes
  `memory=...` to `C:\Users\<you>\.wslconfig` and restarts Docker Desktop);
  `run_crisprscreen setup --docker-memory` does it without asking.
- **Change it per project** in `config.yaml` (`max_cpus` / `max_memory_gb`):
  see "Exactly what max_cpus and max_memory_gb change" below.
- **Inside the budget**, steps that the cluster ran as separate jobs (each
  sample's FastQC, each library's count, each comparison's test) run side by
  side as far as the budget allows; the rest wait their turn. This pipeline's
  steps all needed less than 1 GB on the cluster, so 2.5 GB is enough, with
  room for only one or two steps at a time; a run refuses to start, with an
  explanation, if the budget is below that.
- **Disk:** plan about 2× the size of your raw FASTQ files (per-sample FASTQs and
  `metadata\work\`), plus space for results. The raw FASTQs are not
  copied. `metadata\work\` and `metadata\.nextflow\` are kept after a run
  (`keep_work: true` in `config.yaml`, the default), so a re-run only redoes what
  changed. Files in `results\` and `original_data\` are hard links to the files in
  `metadata\work\` when the disk allows (the run log says "publish mode: link" or
  "copy"), so big FASTQs are not stored twice; never edit such a file in place
  (replacing or deleting it is fine). To free space: `run_crisprscreen clean`, or
  `keep_work: false` (deletes them after a successful run). Intermediate files
  are in `results\intermediate_files\`. Keep the
  project on an internal SSD:
  on a hard disk drive (HDD), USB drive or network share demultiplexing is
  several times slower, and Windows may respond sluggishly while it runs
  (the launcher warns when the project is on an HDD).

### What protects a small computer

- **The budget is a hard limit.** Docker runs each step with `--memory` set to
  its share of the budget and swap disabled. A step that outgrows it is killed
  inside the container (exit 137) and retried with more memory, up to the
  budget. The computer itself never swaps for the analysis.
- **Windows always keeps at least 3 GB**, whatever the settings say.
- **File cache:** while an analysis or image build runs, the launcher watches
  Windows' free memory and makes Docker drop its file cache when it gets low
  (see "Why not more memory" above).

### Exactly what max_cpus and max_memory_gb change

- **`max_cpus`:** the run is limited to that many cores (Docker `--cpuset-cpus`;
  Nextflow's executor is given the same number). Demultiplexing uses all of
  them, every other step uses 1 core, and steps run side by side up to the
  core count.
- **`max_memory_gb`** (decimals allowed, e.g. `2.5`): Docker `--memory` becomes
  a hard limit of that size. About 1 GB of it is for Nextflow itself, and
  Nextflow packs the steps into the rest. Lower means fewer steps side by side
  and each step is capped: a step that does not fit is retried up to the cap,
  then fails with exit 137. Higher means more steps side by side, never more
  than Docker Desktop allows (the launcher warns when it is capped).
- Results never change, only speed.

## Good to know
- **Where Docker keeps its data:** Docker Desktop stores images in a disk file on
  `C:` by default. If `C:` is nearly full, move it in Docker Desktop →
  Settings → Resources → Advanced → *Disk image location*.
- **Reproducibility:** every tool runs at a pinned version inside the image
  (`docker\Dockerfile`, `docker\env\*.explicit.txt`). `METHODS.md` and
  `metadata\pipeline_info\` record the pipeline's git commit, the image ID,
  Docker's version and every tool's version for each run.
- **Customizing a step** for one project (plot colors, extra MAGeCK
  options, ...) without changing it for the lab: see "Customizing" in the
  project's `AGENTS.md` — copy the script from this folder's `bin\` into the
  project's `utils\` and point `config.yaml` at it. Changing a `<script>_args`
  line in `config.yaml` also re-runs that step.
- **Updating the pipeline:** `git pull` in the CRISPR_SCREEN folder. If the
  update changed the Docker image, the next run rebuilds it automatically.

For how each step works (anchor discovery, read classification, the QC
metrics, output layout), see `ANALYSIS.md` in this folder — the analysis is
identical in all three versions; only how it is launched differs.
