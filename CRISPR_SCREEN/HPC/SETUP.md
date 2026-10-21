# One-time setup on the PMACS cluster (HPC version)

Audience: an AI agent helping a lab member set up and run the CRISPR_SCREEN
pipeline on PMACS (LSF). The lab member is a wet-lab scientist, not a
programmer. Everything below comes from the files in this folder (`HPC/`).
Anything that lives outside this folder, or that the repository does not
state, is listed at the end under "Not verified from the repository".

On the cluster the repository is cloned in full to `/project/rephilllab/CRISPR_SCREEN`
(lab convention, see `AGENTS.md`). A sparse clone of `HPC/` alone should also
work, because `run_pipeline.sh` only uses paths inside `HPC/`, but that was not
tested.

## 1. Lab member (normal user): one-time setup

Do this once, on a login node:

```bash
bash /project/rephilllab/setup_environment.sh
```

This script lives in the lab's shared `/project/rephilllab/` space, not in this
repository. It appends one guarded line to `~/.bashrc`:

```bash
[ -f /project/rephilllab/.rephilllab.profile ] && . /project/rephilllab/.rephilllab.profile
```

The line is only added if it is not already there, so re-running the script does
nothing. The script never overwrites or removes anything else in `~/.bashrc`.
The shared profile puts the lab's pipeline launchers on `PATH`. Then load it:

```bash
source ~/.bashrc      # or open a new terminal
which run_crisprscreen.sh
```

`which` should print a path inside `/project/rephilllab/bin/`, the shared folder
that holds `run_crisprscreen.sh` (a symlink, see section 2). If it prints
nothing, the profile is not being loaded: check `~/.bashrc` for the line above.

Labmates never edit `/project/rephilllab/CRISPR_SCREEN` and never add anything
else to their `~/.bashrc`. Lab-wide shell configuration belongs in
`/project/rephilllab/.rephilllab.profile`.

## 2. Lab administrator or maintainer: shared cluster setup (once for the lab)

Only needed once per lab, by someone with write access to `/project/rephilllab/`.
Labmates do not repeat these steps.

1. Clone the repository to the shared location:
   ```bash
   git clone <repository URL> /project/rephilllab/CRISPR_SCREEN
   ```
   The repository does not state its clone URL (see "Not verified").
2. Point the shared launcher at this folder's pipeline script (once; re-run it
   only if the pipeline is moved):
   ```bash
   ln -sfn /project/rephilllab/CRISPR_SCREEN/HPC/run_pipeline.sh /project/rephilllab/bin/run_crisprscreen.sh
   ```
   `run_pipeline.sh` resolves this symlink with `readlink -f`, so the pipeline
   finds its own files in `HPC/`. Future lab pipelines add their own symlink to
   the same `/project/rephilllab/bin/`.
3. The shared Nextflow launcher is not installed by an administrator. On the
   first run, `run_pipeline.sh` downloads it with
   `curl -fsSL -o HPC/.nextflow-bin/nextflow https://get.nextflow.io` (only if
   that file is missing or not executable), so it sits in the shared clone and
   every labmate uses the same pinned version. `.nextflow-bin/` is git-ignored.
   The first run therefore needs outbound network access.
4. Container cache: `/project/rephilllab/nf_containers`. Apptainer caches the
   per-tool container images there (`cacheDir` in `conf/pmacs.config`), so
   they are pulled once and reused by everyone. Image files (`*.sif`, and
   `*.img` as named in `AGENTS.md`) are git-ignored and never versioned.
5. Scratch: `/project/rephilllab/nf_work/`. `run_pipeline.sh` creates
   `nf_work/<user>/<project>/` for each run (`mkdir -p`) and uses it as
   `NXF_WORK`. Nothing else needs to be set up. Intermediate files stay off
   `/home` and off node-local `/tmp`.
6. Shared profile: `/project/rephilllab/.rephilllab.profile`. This is where
   lab-wide shell configuration goes, and it is what puts `/project/rephilllab/bin`
   on every labmate's `PATH` (see section 1). Its contents are not in this
   repository.

Requirements recorded in `README.md`: the pipeline was verified on PMACS
(`pennhpc`, RHEL 9.6, LSF 10.1.0.13) on 2026-09-22 with Apptainer 1.4.5, Java 17
and outbound network access to quay.io and hub.docker.com. Nothing needs to be
installed by an administrator. If the cluster has changed since then, re-check
these. If the cluster is not PMACS, `conf/pmacs.config` documents a FIPS-mode
fix that may also be needed.

## 3. First project

```bash
mkdir ~/my_screen && cd ~/my_screen
run_crisprscreen.sh init
```

`init` creates the project skeleton, fetches `config.yaml`, `examples/`,
`AGENTS.md` and `CLAUDE.md`, `utils/plot_barcode_position.R` and the agent
settings in `.claude/`, then runs `git init` and makes the first commit. It does
not run the pipeline. It never overwrites an existing `config.yaml`.

Then:

1. Copy the example CSVs you need from `examples/` into the project root and
   fill them in: `samplesheet.csv`, `library.csv`, and `contrasts.csv` only if
   you will run `mageck test`. `README.md` in this folder gives the columns.
2. Edit `config.yaml`: set `project_name`, and either `fastq_r1` (a multiplexed
   FASTQ the pipeline demultiplexes, plus `fastq_r2` if paired-end) or
   `fastq_dir` (a folder of already-demultiplexed FASTQs). The FASTQ files stay
   where they are. Use an absolute path, or one relative to the project folder.
3. Run the pipeline (no arguments):
   ```bash
   run_crisprscreen.sh
   ```

For a project's details, read the `README.md` and `AGENTS.md` in this folder
(`HPC/`). `init` also copies an `AGENTS.md` into the project folder for an AI
agent working there; follow that one for project work. To check the pipeline
on synthetic data before using real data, see `test/README.md` in this folder.

## 4. Hand-off to the user after a run is submitted

- `run_crisprscreen.sh` (no arguments) checks the FASTQ paths and submits a
  Nextflow head job to LSF, queue `rhel9` (`bsub` in `run_pipeline.sh`: 1 core,
  8 GB, 72 h). The head job submits the individual tool steps as their own LSF
  jobs, each in its own container.
- Monitor with:
  ```bash
  bjobs -u $USER
  tail -f metadata/nf_logs/nextflow_<timestamp>.out
  ```
  The launcher prints the exact log file names when it submits. The Nextflow
  execution report is `metadata/pipeline_info/execution_report.html`.
- If a run fails, run `run_crisprscreen.sh` again. Nextflow resumes and only
  the failed steps re-run.
- Results are in `results/`. QC reports are in `results/counts/qc/report/`, and
  `METHODS.md` (the Materials and Methods of the run) is written to the project
  folder and to `results/`.
- The work folder (`/project/rephilllab/nf_work/<user>/<project>/`, not in the
  project folder) is deleted after a successful run unless `config.yaml` sets
  `keep_work: true` (default `false`). Failed or stopped runs keep it.
  `publish_mode` is `copy` (`auto` also means copy), so `results/` holds real copies.
- To run a project's own analysis script as an LSF job (never on the login node):
  `run_crisprscreen.sh script utils/<name>.R` (new, not yet tested on the cluster).
  Its logs are in `metadata/script_logs/`.
- Never run Nextflow, analysis scripts, or other memory-heavy commands on the
  login node. It has a memory cap of about 2 GB (`AGENTS.md`). Use the LSF
  job submitted by the launcher, or a test submitted with `bsub`.

## 5. Not verified from the repository

These are not stated in the repository, or the repository is not available
from this folder. Do not assume them; check them on the cluster.

- The exact contents of `/project/rephilllab/setup_environment.sh` and of
  `/project/rephilllab/.rephilllab.profile`. This file only describes them.
  The exact `PATH` line is unknown.
- The output of `which run_crisprscreen.sh` on a labmate's account (described
  as "in `/project/rephilllab/bin/`", but not seen).
- Permissions and group ownership of `/project/rephilllab/bin/`,
  `/project/rephilllab/nf_containers/` and `/project/rephilllab/nf_work/`
  (who can write, and whether labmates share a group).
- How `/project/rephilllab/nf_containers/` was populated, and whether the first
  run pulls images from quay.io or whether they are already there.
- The container file suffix. `AGENTS.md` says `<image>.img` and `.gitignore`
  says `*.sif`. The actual suffix on disk was not checked.
- Whether `python3` exists on login nodes and on compute nodes. The test
  generator needs it (see `test/README.md`). The pipeline's own steps run
  Python inside their containers.
- Whether `git` (used by `init`) and `curl` (used to fetch Nextflow) are
  available on the login node. `init` makes its first commit with
  `git commit ... || true`, so a commit that fails (for example with no git
  user name or email configured) is silent. Check `git log` in the new project.
- Whether outbound network access to `get.nextflow.io` is allowed for labmates.
  The repository only records access to quay.io and hub.docker.com.
- The clone URL of the repository (the repository has no remote configured).
- Whether `/project/rephilllab/CRISPR_SCREEN` is currently a full or sparse
  clone, and whether it has been cloned at all.
- `AGENTS.md` mentions `/project/rephilllab/NEXTFLOW_PORTING_GUIDE.md`. That
  file is outside this repository, so its content was not checked.
- The verification of 2026-09-22 (Apptainer 1.4.5, Java 17, LSF 10.1.0.13)
  is what the repository records. The current versions were not checked.
