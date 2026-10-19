# CRISPR_SCREEN for Mac (Docker)

The Phillips Lab's pooled CRISPR screen pipeline — demultiplex → FastQC/MultiQC →
`mageck count` → guide-representation QC (tables, plots, an HTML report per
library) → optional `mageck test` per comparison → a ready-to-adapt
`METHODS.md` — running on **your own Mac** (Apple chip or Intel). It is the
same pipeline as on the PMACS cluster (same tools, same versions, same outputs),
packaged in one Docker image, so you don't install any bioinformatics software
yourself.

> **Easiest way:** ask Claude (or another AI agent) to do it for you, e.g.
> *"Clone the lab's CRISPR_SCREEN pipeline from GitHub and analyze my screen:
> the FASTQ is in ~/Sequencing/run7, the samplesheet and library are attached."*
> The agent follows `AGENTS.md` in this folder. Everything below is what it
> does, in case you prefer to do it yourself.
>
> **First time on this Mac?** Read `SETUP.md` in this folder (Docker Desktop,
> Rosetta, memory, file sharing, known hiccups), then come back to the steps below.

## What you need

- macOS 13 (Ventura) or newer.
- **6 GB of memory or more** (8 GB or more lets more steps run side by side; the pipeline itself needs about 2.5 GB).
- Free disk space: about **2× the size of your raw FASTQ files** on the project's
  disk (per-sample FASTQs and, while a run is going, `work/`), plus space for
  results and ~8 GB once for the Docker image. The raw FASTQs are not copied.
- Internet the first time (to download Docker and build the image).

## One-time setup (per Mac)

1. Get the pipeline (once), in Terminal:
   ```bash
   git clone --filter=blob:none --sparse <the lab's GitHub URL> ~/CRISPR_SCREEN
   cd ~/CRISPR_SCREEN
   git sparse-checkout set MAC
   ```
   Only the `MAC` folder is downloaded, and the steps below run from it (`~/CRISPR_SCREEN/MAC`); `git pull` later updates it.
   (or download the ZIP from GitHub and unzip it).
2. Run:
   ```bash
   bash ~/CRISPR_SCREEN/MAC/run_crisprscreen.sh setup
   ```
   - If Docker Desktop isn't installed, it installs it with Homebrew if you
     have it, or tells you where to download it (choose *Apple chip* or *Intel
     chip*). Then **open Docker from Applications** once and accept its terms
     (no account needed — click *Skip*), and run `setup` again.
   - On an Apple-chip Mac it installs Rosetta 2 if needed (the pipeline's
     tools are Intel/amd64 builds, identical on every computer).
   - It then builds the pipeline's Docker image: ~3 GB download, 10–30 minutes,
     once.
   - It adds `run_crisprscreen.sh` to your PATH (`~/.zshrc`), so new Terminal
     windows can use it from any folder.

## Analyze a screen

1. **Create a project folder** (one per screen, *outside* the CRISPR_SCREEN folder):
   ```bash
   run_crisprscreen.sh init ~/Documents/CRISPR_screens/my_screen
   ```
   It creates the folder skeleton, `config.yaml`, `examples/` (starter CSVs),
   `AGENTS.md`/`CLAUDE.md` (instructions for an AI agent), `RUN_ANALYSIS.command`
   (double-click) and `utils/run_crisprscreen.sh` (the command, used below).
   It also copies the default barcode-position plot script into `utils/` (only
   if absent) and the agent definitions into `.claude/` (missing files only),
   and both go into the project's first git commit. It doesn't run anything.
2. **Add your tables** to the project folder (copy them from `examples/`;
   edit in Excel/Numbers and export as *CSV UTF-8*):
   - `samplesheet.csv` — one row per sample: `barcode_id,barcode,sample_name`
     (barcodes only when the pipeline demultiplexes; optional `library` column
     when samples come from different library pools).
   - `library.csv` — your sgRNA library: `gene,gRNA_ID,sequence`.
   - `contrasts.csv` (optional — when you're ready for `mageck test`):
     `treatment,control`, one row per comparison.
3. **Edit `config.yaml`** (TextEdit is fine — keep it plain text): set
   `project_name`, and either
   - `fastq_r1:` the multiplexed FASTQ (e.g. `/Volumes/SeqDrive/run7/Undetermined_S0_L001_R1_001.fastq.gz`), or
   - `fastq_dir:` a folder of already-demultiplexed FASTQs, one per sample
     (`<sample_name>.fastq.gz`).

   Tip: drag a file from Finder into the Terminal window to see its path.
   The FASTQs **stay where they are**. The analysis reads them from that path
   and never moves, copies, renames or deletes them; `config.yaml` keeps
   pointing to them. Keep that path stable (if you move the files, update
   `config.yaml`). Several projects can point at the same raw file. Docker
   Desktop shares only `/Users`, `/Volumes`, `/private`, `/tmp` and
   `/var/folders` by default; for a FASTQ elsewhere, add its folder in Docker
   Desktop > Settings > Resources > File sharing. The path must not contain a
   comma. A FASTQ outside the project folder is mounted read-only into Docker.
   A FASTQ that is a symlink must be given by its real path.
4. **Run it:** double-click **`RUN_ANALYSIS.command`** in the project folder,
   or in Terminal in that folder: `./utils/run_crisprscreen.sh`.
   - The analysis runs inside Docker in the background; the window shows its
     progress. Closing the window does **not** stop it. Keep the Mac plugged
     in with the lid open until it finishes (it's kept from idle-sleeping
     while the analysis runs).
   - Check on it any time: `./utils/run_crisprscreen.sh status`.
5. **Results** are in `results/` — start with
   `results/counts/qc/report/*_qc_report.html` (open in a browser; every plot
   is explained) and `results/demultiplexing/qc/multiqc/multiqc_report.html`.
   `METHODS.md` describes exactly what was done, for your paper.

**If something fails**, `./utils/run_crisprscreen.sh status` says which step and
why (details in `metadata/failed_tasks/`). Fix the input and run again: steps
that already finished are reused, not redone.

**Re-running a finished project** recomputes everything (about 10-20 min for a
full lane), because `work/` (and `.nextflow/`) is deleted after a successful run. To reuse the
finished steps, e.g. after adding contrasts, set `keep_work: true` in
`config.yaml` before the run.

## All commands

Run from the project folder (`./utils/run_crisprscreen.sh ...`) or from anywhere
after setup (`run_crisprscreen.sh ...`, acts on the current folder):

| Command | |
|---|---|
| `run_crisprscreen.sh setup` | one-time: Docker Desktop + the pipeline's Docker image |
| `run_crisprscreen.sh init FOLDER` | create a project folder |
| `run_crisprscreen.sh` (or `run`) | run / resume the analysis of the current folder (`--detach`: don't watch) |
| `run_crisprscreen.sh run other_config.yaml` | ... with another config file |
| `run_crisprscreen.sh status` | running? finished? which step failed and why? |
| `run_crisprscreen.sh logs` | watch the live log (Ctrl+C stops watching, not the analysis) |
| `run_crisprscreen.sh stop` | stop the analysis (running again resumes it) |
| `run_crisprscreen.sh clean` | delete `work/` now (rarely needed: it is deleted after a successful run (with `.nextflow/`) unless `keep_work: true`; use after a failed run) |

## Resources: how much of the Mac an analysis uses

Every run gets one **budget** and can never exceed it — Docker enforces it as
a hard limit, so the analysis cannot freeze the Mac by running it out of memory:

| | default | why |
|---|---|---|
| memory | **what the pipeline can use (8 GB), but never more than 75% of the RAM, always leaving 3 GB to macOS** | more would not run faster (below); a plain percentage starves small laptops (80% of 8 GB leaves 1.6 GB) |
| CPU cores | **all but 1** (all but 2 from 8 cores up) | the desktop stays responsive; heavy CPU use only slows a computer, it can't harm it |

…and never more than Docker Desktop itself may use. Examples: 6 GB / 8 cores →
3 GB, 6 cores (on a 6 GB Mac, Docker Desktop's own memory setting may hold it
lower; then that is the limit); 8 GB / 8 cores → 5 GB, 6 cores; 16 GB / 10 cores
→ 8 GB, 8 cores; 32 GB / 12 cores → 8 GB, 10 cores.
Each run prints its budget at the start.

- **Why not more memory on a big Mac:** every step of this pipeline needs
  less than 1 GB. Extra memory would only fill with Docker's file cache (the
  FASTQs it has read), which macOS gets back slowly
  (`docker/requirements.conf`, `max_useful_memory_gb`).

- **Docker Desktop's own limit:** Docker Desktop → Settings → Resources →
  Memory caps everything. `run_crisprscreen.sh setup` shows it and, if it is
  what holds the analysis back, says what to set it to (Apply & restart).
- **Change it per project** in `config.yaml` (`max_cpus` / `max_memory_gb`):
  see "Exactly what max_cpus and max_memory_gb change" below.
- **Inside the budget**, steps that the cluster ran as separate jobs (each
  sample's FastQC, each library's count, each comparison's test) run side by
  side as far as the budget allows; the rest wait their turn. This pipeline's
  steps all needed less than 1 GB on the cluster, so 2.5 GB is enough, with
  room for only one or two steps at a time; a run refuses to start, with an
  explanation, if the budget is below that.
- **Disk:** plan about 2× the size of your raw FASTQ files (per-sample FASTQs and
  `work/` while a run is going), plus space for results. The raw FASTQs are not
  copied. `work/` and `.nextflow/` are deleted after a successful run unless `keep_work: true` (in
  `config.yaml`, default `false`); failed or stopped runs keep it, so a re-run
  resumes. Intermediate files are in `results/intermediate_files/`.

### What protects a small computer

- **The budget is a hard limit.** Docker runs each step with `--memory` set to
  its share of the budget and swap disabled. A step that outgrows it is killed
  inside the container (exit 137) and retried with more memory, up to the
  budget. The Mac itself never swaps for the analysis.
- **macOS always keeps at least 3 GB**, whatever the settings say.

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

- **Apple chip Macs** run the pipeline's tools through Rosetta 2. Make sure
  Docker Desktop → Settings → General → *Use Rosetta for x86_64/amd64
  emulation on Apple Silicon* is on (the default) — it is several times faster
  than the alternative.
- **Where to keep projects:** an internal disk is fastest. Avoid iCloud-synced
  folders (Desktop/Documents with *Optimize Mac Storage*) for big data — macOS
  may offload files to the cloud.
- **Reproducibility:** every tool runs at a pinned version inside the image
  (`docker/Dockerfile`, `docker/env/*.explicit.txt`). `METHODS.md` and
  `metadata/pipeline_info/` record the pipeline's git commit, the image ID,
  Docker's version and every tool's version for each run.
- **Customizing a step** for one project (plot colors, extra MAGeCK
  options, ...) without changing it for the lab: see "Customizing" in the
  project's `AGENTS.md` — copy the script from this folder's `bin/` into the
  project's `utils/` and point `config.yaml` at it.
- **Updating the pipeline:** `git pull` in the CRISPR_SCREEN folder. If the
  update changed the Docker image, the next run rebuilds it automatically.
- **"permission denied"** running `./run_crisprscreen.sh` in the pipeline
  folder: use `bash run_crisprscreen.sh ...` (or `chmod +x run_crisprscreen.sh`).

For how each step works (anchor discovery, read classification, the QC
metrics, output layout), see `ANALYSIS.md` in this folder — the analysis is
identical in all three versions; only how it is launched differs.
