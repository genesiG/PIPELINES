# CRISPR_SCREEN

The Richard Phillips Lab's pipeline for pooled CRISPR screens (Nextflow): demultiplex
the sequencing run by inline barcode → read QC (FastQC/MultiQC) → find each
guide in its read and count it (`mageck count`) → guide-representation QC with
an HTML report per library → optional hit calling (`mageck test`) → a
ready-to-adapt `METHODS.md`.

Pick the folder for the computer you will run it on; each is complete on its
own and runs the same steps with the same tool versions:

| Folder | Runs on | First time (once) | Then read |
|---|---|---|---|
| [`HPC/`](HPC/README.md) | the PMACS cluster (LSF + Apptainer) | [`HPC/SETUP.md`](HPC/SETUP.md) | [`HPC/README.md`](HPC/README.md) |
| [`WINDOWS/`](WINDOWS/README.md) | your own Windows PC, in Docker | [`WINDOWS/SETUP.md`](WINDOWS/SETUP.md) | [`WINDOWS/README.md`](WINDOWS/README.md) |
| [`MAC/`](MAC/README.md) | your own Mac, in Docker | [`MAC/SETUP.md`](MAC/SETUP.md) | [`MAC/README.md`](MAC/README.md) |

Each folder also has `ANALYSIS.md`, which describes the analysis itself (the same
text in all three), and `test/`, a small synthetic screen for checking a change or
a new computer's setup (see `<folder>/test/README.md`).

## Get only the folder you need

The repository holds all three folders, but you only need one. Clone it as a
sparse checkout:

```bash
git clone --filter=blob:none --sparse https://github.com/<owner>/CRISPR_SCREEN.git
cd CRISPR_SCREEN
git sparse-checkout set WINDOWS     # or MAC, or HPC
```

- The root files (`README.md`, `.gitattributes`, `.gitignore`) always come
  along. The other folders, including `tools/`, do not: the clone holds only the
  folder you chose.
- `git pull` later updates that folder only.
- To add another folder later: `git sparse-checkout add MAC`.
- Each folder is self-contained: its own `README.md`, `SETUP.md`, `AGENTS.md`,
  `CLAUDE.md`, subagent definitions in `.claude/agents/`, and `test/`.

**First time on this computer:** read `<folder>/SETUP.md` and do it once. Then
use `<folder>/README.md`.

AI agents (Claude Code and similar): start them inside the folder you chose, so
they read that folder's `CLAUDE.md` and `AGENTS.md`.

In every version you create a project folder with `run_crisprscreen init`,
fill in `samplesheet.csv`, `library.csv` (and optionally `contrasts.csv`) and
`config.yaml`, then run `run_crisprscreen`. An AI agent can do all of this for
you once it is started in the folder you chose.

[`HPC/ANALYSIS.md`](HPC/ANALYSIS.md), [`WINDOWS/ANALYSIS.md`](WINDOWS/ANALYSIS.md) and
[`MAC/ANALYSIS.md`](MAC/ANALYSIS.md) describe the analysis itself in detail (steps,
automatic guide positioning, read classification, outputs, Methods). The three
files are identical, since the analysis is the same in all three versions.

## For maintainers

`tools/` (the maintenance script `tools/sync_os_folders.sh`) is for maintainers
only and is not needed by users. Changes to the pipeline logic usually belong in
all three folders: edit `HPC/`, mirror the change into `WINDOWS/`, then run
`bash tools/sync_os_folders.sh WINDOWS MAC`. `bash tools/sync_os_folders.sh --check`
verifies that the copies match.
