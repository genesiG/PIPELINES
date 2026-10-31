# AGENTS.md — CUT&RUN / chromatin profiling pipeline, Docker version for **Windows** (Richard Phillips Lab)

This folder is the lab's CUT&RUN / chromatin profiling pipeline (Nextflow: trimming and
bowtie2 alignment against a target + spike-in genome, duplicate removal, ChIPseqSpikeInFree
or E. coli spike-in normalisation, bigWig tracks, csaw peak calling, optional region
classification, heatmaps, tracks, annotation, intersections, and
the QC report), packaged to run on a labmate's own Windows computer inside Docker. The user
is usually a **wet-lab scientist**: do the technical work for them, explain in plain
language, and ask them only for what you can't know (their samples, their files, their
comparisons).

First time on this computer (Docker not set up, or never checked)? The lab's shared
Docker setup is `SETUP.md` (one-time, shared by the lab's Docker pipelines). It is not in
this folder yet, so until it is, the steps in A.1 below are this pipeline's part of the
setup.

Almost every task is **A. Set up / run / check an analysis**: follow the runbook
below. You never edit anything in this pipeline folder for that. Changing the pipeline
itself is not part of an analysis. If the user asks for a change to the pipeline for the
whole lab, ask them to contact the maintainers, and do not edit this folder.

If the current folder is a **project** folder (it has `config.yaml`,
`samplesheet.csv`, `results/`), read that project's own `AGENTS.md` instead.

## How to call the launcher

Everything goes through one command, `run_cutandrun.bat` (which starts
`run_cutandrun.ps1`), in this folder (full path: `<this folder>\run_cutandrun.bat`).
From the agent's shells:

```powershell
# PowerShell (preferred)
& "D:\path\to\CHROMATIN_PROFILING\WINDOWS\run_cutandrun.bat" status
```
```bash
# Git Bash
cmd.exe //c "D:\\path\\to\\CHROMATIN_PROFILING\\WINDOWS\\run_cutandrun.bat" status
```

| Command | What it does |
|---|---|
| `run_cutandrun.bat setup --yes` | one-time per computer: Docker Desktop check/install, pull the pipeline's Docker image (`docker\image.txt`), print the tool versions, add the command to the user's PATH |
| `run_cutandrun.bat init "<project folder>"` | create a project folder: `original_data\fastq\`, `metadata\`, `results\`, `utils\`, `examples\` (CSV templates), `config.yaml`, launchers, git |
| `run_cutandrun.bat run --detach` (in the project folder) | check the input files and the disk, set the memory budget, start the pipeline in Docker, return immediately |
| `run_cutandrun.bat status` (in the project folder) | running / succeeded / failed, step counts for this run, the failed step, the last 25 log lines, the Docker work volume's name and size |
| `run_cutandrun.bat logs` | follow the live log (blocks until the run ends — avoid in an agent; use `status`) |
| `run_cutandrun.bat stop` | stop a run (finished steps are kept; the next run resumes) |
| `run_cutandrun.bat clean --yes` | delete Nextflow's work folder and resume cache to free disk: the project's Docker work volume, and `metadata\work\` and `metadata\.nextflow\` if present — the next run starts from scratch, results are kept; without `--yes` it asks first; ask the user first |
| `run_cutandrun.bat script utils\<name>.<R\|py\|sh> [args]` (in the project folder) | run the project's own analysis script in the pipeline's image (foreground; log in `metadata\script_logs\`); see "Analyses beyond the pipeline" below |
| `run_cutandrun.bat help` | the command list |

Options: `--detach` (run: start and return), `--yes` (answer yes to every question: only
after the user has agreed to what is being asked), `--no-wsl-resize` (run: leave Docker's
memory limit alone; see A.1).

Every project folder also gets its own launcher, `utils\run_cutandrun.bat`
(same commands, always for that folder; agents run `.\utils\run_cutandrun.bat`
from the project folder), and `RUN_ANALYSIS.bat` (double-click for humans: it
pauses at the end, so **agents must not call `RUN_ANALYSIS.bat`**).

Long commands: `setup` (image pull) and a run can take a long time. Run
`setup` in the background or with a long timeout, and always start runs with
`--detach`, then poll `status` every few minutes.

## A. Runbook: from "here are my FASTQs" to results

### 1. One-time computer setup

The steps below are this pipeline's part of the setup (see the note at the top about
`SETUP.md`).

1. Check: `docker --version` and `docker info`. If both work, go to step 3.
2. `run_cutandrun.bat setup --yes`. What can happen:
   - Docker Desktop not installed → it installs WSL 2 + Docker Desktop with
     `winget`. **Windows shows an administrator (UAC) prompt: tell the user to
     approve it.** Then the user must **restart the computer**, open Docker
     Desktop once from the Start menu and accept its terms (no account needed
     — "Skip"). Explain this plainly, then wait for them and re-run `setup`.
   - Docker installed but not running → it starts Docker Desktop and waits.
   - "virtualization is disabled" / WSL errors → hardware virtualization must
     be enabled in the BIOS/UEFI; that's for the user or their IT department.
3. `setup` pulls the image (`ghcr.io/genesig/omics`, tag in `docker\image.txt`, always
   run as `linux/amd64`) and prints the versions of Nextflow, R, Python, samtools and
   bowtie2. Re-running `setup` is harmless. A run also pulls the image if it is missing,
   so a forgotten `setup` only costs time.
4. Resources: every run gets one budget, enforced by Docker as a hard limit. It is the
   smaller of: 75% of the computer's RAM (always leaving at least 3 GB to Windows), all
   logical cores but one (two from 8 cores up), what Docker Desktop may use, and the most
   this pipeline can put to use (14 GB, `docker\requirements.conf`). `max_cpus` and
   `max_memory_gb` in `config.yaml` can only lower it. The run prints the budget it chose
   ("Uses at most N CPU cores and X memory"). The pipeline needs at least 5 GB
   (`min_memory_gb`) and refuses to start below that, with an explanation. Both figures
   are placeholders until a pilot run on a real computer has been measured.
5. **Docker's own memory (the WSL 2 limit).** Docker Desktop runs in a WSL 2 virtual
   machine whose memory is set by the `memory=` line in `%USERPROFILE%\.wslconfig`. For
   chromatin runs the launcher sets it by the computer's RAM: under 12 GB → 6 GB, under
   24 GB → 8 GB, otherwise 16 GB. This is changed for an analysis only, and the old value
   is put back when the analysis finishes (or at the next launcher command if the window
   was closed early). The launcher **asks first**; the user must agree. While other
   containers are running it leaves the limit as it is, and the analysis runs with what
   Docker has. `--no-wsl-resize`, or
   `wsl_resize: false` in `config.yaml`, leaves it alone. Editing `.wslconfig` is a change
   to the user's computer: **ask the user first**, and only pass `--yes` after they agreed,
   because `--yes` answers yes to this question too.

### 2. Gather what you need from the user (ask; never guess)

- **A name for the project** → `project_name` (letters, digits, `_`; used as a prefix on
  the QC report and other project-level outputs) and the project folder's name.
- **Where the FASTQ files are**, and whether each sample is:
  - **FASTQ pairs**: `fastq_1` and `fastq_2` (absolute paths, one pair per library),
    with `paired_end: true` (the pipeline runs paired-end reads only; if the data are
    single-end, say so before going further); or
  - **a pre-aligned, final-filtered BAM**: the `bam` column, with the FASTQ columns
    left empty. Such a row skips trimming, alignment and filtering, and it has no spike-in
    BAM, so it cannot be used for spike-in normalisation.
- **The reference files** (ask for the paths, never guess them): the bowtie2 index prefix
  of the combined target + spike-in genome (`genome_index`), `chrom_sizes`, `blacklist`,
  `gtf` (`tss_bed` only when `gtf` is not set), and the spike-in contig prefix (`spike_prefix`, e.g. `ecoli_chr`) and
  mitochondrial contig name (`mito_chrom`). The lab PC keeps references under
  `D:\REFERENCES` (`bowtie2\t2t_ecoli`, `blacklist`, `chr_size`, `genomes\t2t`,
  `annotation`); that is an example layout, not a requirement.
- **The samples**: names, which antibody/mark each one is, which are IgG controls, and
  the sample-to-control pairing. Sequencing-core sheets often hold these; ask for the sheet.

If they give Excel files: ask them to "Save As → CSV UTF-8", or, if Excel is
installed on this computer, convert them yourself through Excel's COM object in
PowerShell. Check the result opens as plain comma-separated text.

**Then ask about the experiment itself (ask; never assume).** Only the user
knows how it was done; ask even when it seems obvious:

- **The question:** what the experiment is meant to show (a histone mark's spread, a
  factor's binding, a change after a treatment or a knockout, persistence of a mark), and
  which genes or loci they expect to behave how (positive and negative control loci).
- **The assay:** CUT&RUN or CUT&Tag, which kit or protocol, and the sequencing (paired-end,
  read length, which run(s) and lanes, whether the run also holds other people's samples).
- **The antibodies:** which mark or protein each one targets, and what each is expected to
  look like (broad domains or sharp peaks). A "broad" mark called as narrow peaks, or the
  reverse, is worth a question.
- **The biology:** cell line or model; the treatments, doses and timepoints; the genotype.
- **Each sample:** its treatment, which samples are biological replicates and which are
  technical, which IgG control belongs to which sample, and anything that went wrong (low
  cell numbers, a failed step, a contamination, a re-sequenced sample).
- **The spike-in:** whether the E. coli spike-in was added on purpose or is carry-over
  from the sample (this decides the normalisation), and whether any sample has no
  spike-in reads.
- **The normalisation:** which `norm_run` each sample belongs to (ChIPseqSpikeInFree,
  E. coli spike-in, TMM or CPM), and whether a global change is expected (a treatment
  that lowers a mark everywhere would not move a ChIPseqSpikeInFree factor).
- **The comparisons:** which treatment-vs-control comparisons answer their question, and
  which regions count as **retained** or **lost** in each.

- Ask in short batches (3-6 questions at a time), in plain language, and say why
  each answer matters. Keep asking until the answers are specific; "the usual"
  or "like last time" is not an answer, so ask what it was.
- **Never infer the design from names** (`IgG`, `K27me3`, `WT`, `rep2`, `D0`, `ctrl`,
  `hi`/`lo`) or from what is "typical" for chromatin work. Ask what each name means.
- Use only what the user confirmed. Where an answer is missing, say what you
  would do instead and wait for an explicit OK; never fill the gap silently,
  and never pick the controls, comparisons or replicates yourself.
- **Write the answers down** in `CONTEXT.md` in the project folder (dated, in the
  user's words, with open questions marked as open). Read it at the start of
  every session; ask the user to confirm it is still right before you rely on
  it, and add to it as you learn more. It is tracked by git with the inputs.
- **When the data seem to contradict the context, flag it; don't reinterpret.**
  Examples: an IgG control with a high FRiP (its "peaks" are real signal, so the
  control is not one); a replicate that correlates poorly with its partner; a "global
  loss" treatment where the ChIPseqSpikeInFree factors do not move; a mark expected to
  be broad that gives narrow peaks; a sample whose spike-in reads are missing. Show the
  user the numbers and the file they come from, say what the data would be consistent
  with, and ask. The user's description stands until they change it; don't rename
  samples, drop them or change comparisons on your own.
- Subagents cannot talk to the user: put the confirmed context (or the path to
  `CONTEXT.md`) in every brief, and pass their "Questions for the user" back to
  the user rather than answering them yourself.

### 3. Create the project folder

Put projects **outside** the pipeline's git clone, e.g.
`%USERPROFILE%\Documents\chromatin_projects\<project_name>` (a drive with enough free
space: see the disk note in step 6; a chromatin dataset needs about 100 GB or more).
Use a new folder: if the one you picked already exists and holds anything, list it and
ask the user before using it (`init` never overwrites files, but it adds its own files and
a git repository to whatever is there). Then:

```powershell
& "<this folder>\run_cutandrun.bat" init "$env:USERPROFILE\Documents\chromatin_projects\<project_name>"
```

`init` writes the project's `config.yaml`, copies the CSV templates into `examples\`
(never overwriting them), and copies the project's `AGENTS.md`/`CLAUDE.md` and the agent
kit (`.claude\agents\`) when this folder has them (`agent\` in this folder). The CSVs the
project uses are copies you make from `examples\` into the project folder itself.

Immediately after `init`, **read the new project's `CLAUDE.md` and `AGENTS.md` and
follow them from then on** (their delegation table and triggers included): users
usually ask you to set the project up before running, so the project's
`.claude\agents\` and instructions appear mid-session. If Claude Code was started
in this pipeline folder, its `.claude\agents\` is the same kit, so the agents are
available now. If they are not listed as available subagent types, either spawn a
general-purpose agent whose prompt is the agent file's body (after the
frontmatter) plus the brief, with the model from the file's frontmatter, or ask
the user to restart Claude Code in the project folder (recommended for every
later session: start Claude Code inside the project folder, so it loads that
project's `CLAUDE.md` and agents).

### 4. Write the input tables (in the project folder)

The full column rules, the enforced checks and the examples are in the project's
`AGENTS.md` ("Input files"); the templates are in `examples\`. Summary of each table:

- **`samplesheet.csv`** (one row per library; named by `samplesheet:` in `config.yaml`):
  `sample_id,fastq_1,fastq_2,bam,antibody,condition,replicate,control_id,group,norm_run,exclude`.
  `sample_id` is unique and uses letters, digits, `_`, `.`, `-` (it names every output
  file). Fill `fastq_1`/`fastq_2` **or** `bam`, never both. `antibody` is the mark or protein
  (`IgG` marks a control library, any case). `control_id` is the `sample_id` of that
  library's IgG control (empty for the IgG rows). `group` is the unit of replication for
  csaw and the name that `peaksets.csv` and `classifications.csv` refer to (empty = the
  condition and antibody joined). `norm_run` names a row of `norm_runs.csv`. `exclude`
  (`true`/`false`) drops a row from the whole run. Rows in `skip_groups` (default
  none) stay in the file but the core steps skip them.
- **`norm_runs.csv`** (one row per normalisation run; `norm_runs:`):
  `norm_run,method,bin_size,cutoff,max_turns[,min_first_turn][,ref_condition]`. `method` is the
  run's PRIMARY method: `chipseqspikeinfree`, `spikein`, `tmm` or `cpm` (library size only, numbers
  blank). The three ChIPseqSpikeInFree numbers must be numeric for that method; `min_first_turn` is
  optional (blank = 0.5); `ref_condition` is optional (the reference condition of the MA plots; blank =
  the condition of the `ref_group` of the classifications that use the run, else the first library's).
  Results go to `results\scale_factors\<norm_run>.scale_factors.tsv` (primary only; every method is in
  `<norm_run>.scale_factors_all.tsv` and `methods/`).
- **`peaksets.csv`** (one row per csaw peak set; `peaksets:`):
  `peakset,target_group,control_group,small_filt,large_filt,merge_bp,max_width,fdr,lfc`.
  `target_group` and `control_group` are `group` values (target versus IgG). A blank
  threshold cell takes the antibody's `csaw_defaults`, then the global `csaw_<column>`;
  `metadata\peaksets.resolved.tsv` shows where each value came from.
- **`classifications.csv`** (used when `run_classify: true`; `classifications:`):
  `classification,regions,ref_group,treat_group,norm_run,lfc_retained,lfc_lost,fdr`.
  `regions` is a peak set name or a BED path. Each region is **lost** (logFC at or below
  `lfc_lost` and FDR below `fdr`), **retained** (logFC at or above `lfc_retained`, or FDR
  at or above `fdr`) or **unclassified**. Blank cells take `classify_lfc_retained` (-1),
  `classify_lfc_lost` (-2), `classify_fdr` (0.05).
- **`heatmaps.csv`** (used when `run_heatmaps: true`; `heatmaps:`): one row per heatmap:
  `heatmap,regions,samples,labels,window_bp,bin_bp,sort_using,cluster_k,z_min,z_max,y_max,colors,engine`
  (optional `cluster_samples,cluster_labels,legend_title,region_labels`). `engine` is
  `deeptools` or `enrichedheatmap` (the latter takes one region only). `window_bp` must be
  a multiple of `bin_bp`. Clustering is seeded (`heatmap_seed`).
- **`tracks.csv`** (used when `run_tracks: true`; `tracks:`):
  `track_set,locus_or_gene,flank,samples,labels,colors,max_value,beds,scalebar`
  (optional `log2ratio,bed_colors,bed_labels,ini_args`). `locus_or_gene` is a gene name or
  `chr:start-end`; an unknown gene is an error. `max_value` is a number, `auto`,
  `group:NAME` or blank (the shared default).
- **`intersects.csv`** (used when `run_intersect: true`; `intersects:`):
  `name,bed_a,bed_b,label_a,label_b,min_overlap,reciprocal[,output]`. `bed_a` and `bed_b` are peak
  set names, class names (`<id>.retained`, `<id>.lost`) or BED paths. `output` = `records` (default) or `overlap`
  (legacy overlap pieces, one per `bed_a`/`bed_b` pair).

Start from `examples\*.csv`; save as UTF-8.

### 5. Edit `config.yaml`

Every key has a default in the pipeline's `nextflow.config`; list only what you change.
The example is `examples\config.yaml`. The key choices:

- `project_name`, `samplesheet`, and the CSV keys above (`norm_runs`, `peaksets`,
  `classifications`, `heatmaps`, `tracks`, `intersects`).
- **Reference files**: `genome_index` (bowtie2 index prefix), `chrom_sizes`, `blacklist`,
  `gtf`, `tss_bed` (fallback for TSS enrichment, only needed when `gtf` is unset), `spike_prefix`, `mito_chrom`, `effective_genome_size` (3100000000 for
  T2T-CHM13). Write Windows paths in full (`D:\REFERENCES\...`) or relative to the
  project folder. The launcher warns when a file is missing, or when `genome_index` has
  no `.bt2` files next to it.
- **Input files stay where they are.** The FASTQs and BAMs are never moved, copied,
  renamed or deleted; drives they name are mounted read-only into Docker
  (`/mnt/<letter>`). Tell the user where the files are read from and that they stay
  there. Keep the paths stable (if the files are moved, update the CSV and
  `config.yaml`). Refused: UNC paths (`\\server\share`) and mapped network drives, because
  Docker Desktop cannot mount them (copy the files to a local or USB drive).
- **Output folders must be relative**: `outdir: ./results`, `metadata_dir: ./metadata`,
  `intermediate_dir: ./intermediate_files`. The launcher refuses a Windows drive path there.
- **Duplicates**: `remove_duplicates: true` and `spikein_dedup: true` (IgG is always
  deduplicated). Setting either to `false` changes every normalisation factor: ask first.
- **Experiment**: `experiment: cutandrun` or `cutandtag` (selects the trimming adapters);
  `paired_end: true`; `trim: true`, `trim_minlen: 25`; `bowtie2_mode: local`; `mapq: 30`.
- **bowtie2 fragment length**: `bowtie2_maxins` and `bowtie2_minins` are `null` by default
  (bowtie2's own defaults, as in the manuscript). Changing either one re-runs alignment and
  everything downstream: ask the user first.
- **Normalisation**: `normalisation` is informational; the `norm_run` column of the
  samplesheet decides. `bigwig_units` (`rpkm`, `cpm` or `none`), `bigwig_bin_size`
  (10 by default), `bigwig_extend_reads: true`.
- **Normalisation methods (2026-10-27).** `method` in `norm_runs.csv` is the run's **primary** method and the only one that scales the bigWigs and the classification (`cpm` is allowed too: library size only; a run with one library needs it). Every run with at least two libraries also computes the methods of `norm_methods` (default `chipseqspikeinfree,tmm,cpm`; `''` = primary only) as diagnostics, and physical spike-in when `norm_spikein: true` (a run with pre-aligned BAMs: `unavailable`, with a warning). They are in `results\scale_factors\<norm_run>.scale_factors_all.tsv` and `results\scale_factors\methods\`; `<norm_run>.scale_factors.tsv` is the primary alone and unchanged. ChIPseqSpikeInFree as a diagnostic takes blank CSV cells from `csif_bin_size` 10000, `csif_cutoff` 1.2, `csif_max_turns` 0.99. Choose the primary with the user before the results are seen and record it in `CONTEXT.md`; do not switch it because another method looks better in the MA plots.
- **Normalisation QC and replicate concordance.** `results\norm_qc\<norm_run>.ma.svg` (+ `.ma.p2.svg`, ...; `norm_qc_panels_per_page` 12) has one MA panel per replicate pair and per treated-vs-reference pair (reference = `ref_condition` of `norm_runs.csv`, else the condition of the `ref_group` of the classifications that use the run, else the first library's); y = log2 CPM of `other` minus `ref` (negative = less signal); the dashed line of each method is where it declares "no change" (log2 of the ratio of the `edger_norm_factor`s, 0 for CPM; dotted = factor from a QC-failed or imputed library); settings `norm_qc_min_count` 20 and `norm_qc_prior` 0.5. A line on the cloud is automatic for TMM and CPM and is not evidence that they are right when the treatment changes the mark globally; between replicates every reliable method should lie on the cloud. `results\qc\replicate_correlation\replicate_correlation.tsv` gives per replicate pair Pearson (log2 CPM, `repcor_prior` 2) and Spearman correlation over `repcor_bin_size` bins (null = `tmm_bin_size`), the same on the top `repcor_signal_frac` 0.10 of bins, the robust spread of the log2 ratio, and the median correlation with the other groups of the same antibody (context). These are descriptive, with no pass/fail thresholds. `run_replicate_correlation` and `run_deeptools_correlation` (the all-sample Spearman heatmap) are both on by default.
- **csaw**: `run_csaw: true`; `csaw_norm_bins` (`background` default, or `all`); the
  threshold defaults `csaw_small_filt`, `csaw_large_filt`, `csaw_merge_bp`, `csaw_max_width`,
  `csaw_fdr`, `csaw_lfc`, and the per-antibody `csaw_defaults`. Changing a threshold is a
  scientific decision: tell the user and get confirmation first.
- Peak callers: `run_csaw`, `run_macs3` and `run_seacr` are all on by default. MACS3 and SEACR need an IgG `control_id` per target: without one they are skipped with a warning, unless `run_macs3` / `run_seacr` were set to `true` explicitly (then a missing control is a startup error). Broad or narrow is decided per group from its antibody against `broad_marks` (default `H3K27me3, H3K27me2, H3K36me3, H3K4me1, H3K79me2, H3K79me3, H3K9me1, H3K9me2, H3K9me3, H4K20me1, H3F3A`: the ENCODE histone ChIP-seq broad list plus H3K27me2 by lab request and H3K9me3, an ENCODE "exception" treated as broad by St. Jude; case-insensitive, histone prefix optional: `K27me3` = `H3K27me3`). `macs3_broad` `auto` | `true` | `false` (default `auto`): broad uses `--broad --broad-cutoff` `macs3_broad_cutoff` (0.1), narrow gives narrowPeak with `-q` `macs3_qvalue`. `seacr_mode` `auto` | `relaxed` | `stringent` (default `auto` = relaxed for broad marks, stringent otherwise); lab policy: SEACR has no broad mode (relaxed = threshold between the knee and the peak of the total-signal curve, stringent = the peak). Groups mixing broad and narrow antibodies are a startup error. The mode per group is logged and written to `metadata\peak_call_modes.tsv`; tell the user which mode each group got. Other optional steps: `run_log2ratio`, `run_csaw_cutoff_scan` (a post-hoc diagnostic that never feeds thresholds); `run_csaw` and `run_report` are on. `call_peaks: false` turns all of them off (and FRiP, peak QC and peak annotation); bigWigs and normalisation are unaffected.
- **Optional steps, all off by default**: `run_classify`, `run_heatmaps`, `run_tracks`, `run_annotation`, `run_intersect`.
  `run_report: true` builds the QC report.
- **Groups to skip**: `skip_groups` (default none).
- **Resources (optional)**: `max_cpus` and `max_memory_gb` lower the run's budget (A.1,
  step 4); use them only if the user needs the computer to stay responsive.
- **Run bookkeeping**: `keep_work: true` (default) keeps the work folder and resume cache so a
  re-run resumes; `work_volume: true` (default) keeps them in the project's Docker volume on
  Docker's disk, `work_volume: false` in `metadata\work\`; `publish_mode` is `auto` by default
  (copies with the work volume; hard links when the disk allows without it; the example
  `config.yaml` sets `copy`).
- **`wsl_resize: false`** leaves Docker's memory limit alone (A.1, step 5).
- Leave `lsf_watchdog` and the other cluster-only keys alone (they do nothing on this
  computer). Leave `fastq_dir` as `null`: this Windows version does not read a FASTQ folder.

### 6. Check, then run

Check before the first run (the launcher also checks, but catching it now is
faster): headers exactly right; `sample_id` unique; FASTQ and BAM paths exist (the
launcher names any file it cannot find); a control named by every `control_id` exists
and is an IgG row; every `norm_run` matches a `norm_runs.csv` row; `group` values match
what `peaksets.csv` and `classifications.csv` use; the comparisons are the ones the user
confirmed.

```powershell
cd "<project folder>"
& .\utils\run_cutandrun.bat run --detach
& .\utils\run_cutandrun.bat status        # every few minutes
```

The run first checks the inputs and prints the disk check: the free space on Docker's disk
(the drive with Docker Desktop's disk image, which holds the work volume) and on the project's
drive (which holds the results), against the total size of the FASTQ/BAM files in the samplesheet
times two rules in `docker\requirements.conf`: `disk_factor=5.6` (work files and results together;
used when the work folder is in the project folder) and `disk_results_factor=1.0` (results only;
used for the project's drive when the work volume holds the work files, where Docker's disk is
checked against about 4.6 times the input). It warns below that on Docker's disk, and stops below
20 GB free on the project's drive (`disk_hard_min_gb`). To estimate by hand, add up the sizes of the
files named in `samplesheet.csv` and multiply by 5.6 (about 4.6 for Docker's disk and 1 for the
project's drive, with the work volume). These multipliers are estimates pending the next run; a
dataset can need about 100 GB or more of working space. Free space can be given back with
`clean` (ask first) or with `keep_work: false`, which deletes the work files after a successful
run (the Docker volume is removed).

Timing on a laptop has not been measured for this port yet (the cluster takes hours for a
full project). Poll `status` every few minutes. Docker's disk image should be on a fast internal
SSD (or a drive with room): the work files are there. The project folder can be on a hard disk
(the launcher warns): with the work volume only the results are
written there (and `intermediate_files\` with `publish_intermediates: true`). The computer must stay on and awake (plugged in, sleep disabled): if it sleeps, the
run pauses; if Docker restarts, just run again (it resumes).

Nextflow's work folder and its resume cache are kept by default (`keep_work: true`), so a
re-run resumes and only redoes what changed (e.g. after changing a comparison, or a plot
setting). Failed or stopped runs always keep them. With `work_volume: true` (the default) they
are in the project's Docker volume, named `<container name>_work` (e.g.
`cutandrun_<project>_<8 hex>_work`), on Docker Desktop's own disk: the work folder is
`/nfwork/work` and the cache `/nfwork/.nextflow` inside the container. The volume is not visible
in File Explorer; `status` shows its name and size. Files in `results\` are copies of the work
files (`publish_mode: auto`). With `work_volume: false` they are `metadata\work\` and
`metadata\.nextflow\` in the project folder (`metadata\` is git-ignored); files in `results\`
and `original_data\` may then be hard links to the files in `metadata\work\` (`publish_mode:
auto`; the run log banner says "publish mode: link ..." or "copy ..."), so big files are not
stored twice. In either case never edit a published file in place (replacing or deleting is
fine). The trimmed FASTQs, unfiltered aligned BAMs, filtered BAMs, SEACR bedGraphs and BIN_COUNTS
tables are in `intermediate_files\` only with `publish_intermediates: true` (the Windows default is
`false`: they stay in the work folder).

### 7. Report results

`status` says `succeeded` → point the user to:
- `results\qc\report\<project_name>_qc_report.html` — the QC report; written for
  non-specialists (start here);
- `results\qc\multiqc\multiqc_report.html` — read and alignment QC;
- `results\scale_factors\<norm_run>.scale_factors.tsv` — the normalisation table
  (`qc_pass`, `note`); it holds the primary method only (`<norm_run>.scale_factors_all.tsv` has every method); check it before interpreting any normalised track or peak set;
- `results\peaks\csaw\<peakset>.peaks.bed` — the peaks; `results\classify\<id>\` — the
  retained, lost and unclassified regions, if run;
- `results\bigwig\` — the tracks; `results\bam\final\` — the final BAMs;
- `metadata\pipeline_info\execution_report.html` — how long each step took and how much
  memory it used. A resumed run keeps the previous run's trace, report and timeline next to it, named with that run's timestamp (execution_trace.<timestamp>.txt).

Summarise from `results\qc\` (the alignment statistics table and the FRiP and peak
statistics) and the QC report. Check the controls first: an IgG with high FRiP or
a poorly correlated replicate is flagged, not explained away. Spike-in percentages come from
the alignment table. For statistics or comparisons with earlier numbers, use the
`statistician` and `output-verifier` agents (see "Subagents").

**Peak QC tracks (new in this version; not yet validated on real bigWigs).** The step `PEAK_QC_REGIONS` (`bin\peak_qc_regions.py`) picks a few regions from the run's peak sets (your own loci, the default Polycomb marker loci, a negative-control locus, the largest, median and densest peaks per mark) and draws one figure per region in `results\tracks\peak_qc\` (the figures need `run_tracks`). Its settings are `run_peak_qc`, `peak_qc_loci`, `peak_qc_negative`, `peak_qc_max_regions`, `peak_qc_exclude_bed`, `peak_qc_nopeak_igg_max_fold` (IgG cap of the peak-free tile, default 1.5) and `peak_qc_shape_marks`; their defaults are in `nextflow.config`. Ask the user for their own loci of interest and do not assume any. Tell the user that this is new, that its signal criteria are untested on real bigWigs (check `regions.detail.tsv`), and that the loci it picks are a suggestion to check, not a result.

### 8. When a run fails

`run_cutandrun.bat status` says which step failed (`FAILED: <step> (exit N)`) and prints
the last log lines. The files to read, in order:

1. `metadata\nf_logs\last_run.status`: `running`, `succeeded`, `failed (exit N)`,
   `failed to start`, `stopped by the user`, or `interrupted` (Docker or the computer
   stopped mid-run; run again to resume).
2. `metadata\nf_logs\nextflow_<timestamp>.log`: everything Nextflow printed. Startup errors
   (bad samplesheet, missing file, unknown `control_id` or `norm_run`, a typo in a key) are
   printed at its end.
3. `metadata\failed_tasks\<timestamp>\<step>\`: the failed step's own files
   (`command.sh` = what ran; `command.err` / `command.log` = what it said).
4. `metadata\pipeline_info\execution_trace.txt`: one row per step (status, exit code,
   `peak_rss` = the memory it used, `realtime`). Use it to see whether a step ran out of memory
   or disk.

Common cases (the project's `AGENTS.md`, "When a run fails", lists more):
- **Exit 137, 140 or 143**: out of memory or killed. These are retried automatically with
  more memory, up to the budget. If the step still fails, check the budget in the run banner
  and the step's `peak_rss`. Raising the budget is the user's decision (`max_memory_gb`, or
  a bigger Docker limit, which means changing `.wslconfig`): ask first.
- **Exit 3**: a known data problem, not retried; the message says what to fix. Typical:
  a scale factor missing or NaN (a sample missing from, or failing, its `norm_run`); a
  control missing for a peak caller that needs one; an empty BAM.
- **Disk full**: the run stops and can be resumed after space is made. Check the free space
  on Docker's disk and on the project's drive (`status` and the disk check in `run`); `clean`
  (ask first) deletes the work volume (or the work folder with `work_volume: false`).
- **Startup errors**: fix the named CSV or key in the project folder, then run again.
- **ChIPseqSpikeInFree aborts**: the pipeline never falls back to an old scale-factor file;
  read the step's log for the package error and the `norm_runs.csv` numbers of that run.

Fix the input or config, then run again: `run --detach`. Finished steps are reused. If it
looks like a bug in the pipeline itself, tell the user and stop: don't patch this folder.

### 9. Analyses beyond the pipeline

The user may want a plot or test the pipeline doesn't do. That is allowed, through
the project's `script` command: `.\utils\run_cutandrun.bat script utils\<name>.<R|py|sh>`
runs the script in the pipeline's image (**Tier 1**), and the log goes to
`metadata\script_logs\`. The script lives in `utils\`, writes only to a new
`results\<name>\` folder, and is designed, written, run and reviewed as described in
"Analyses beyond the pipeline" in the project's `AGENTS.md` (`init` copies it there). The
`script` command is foreground only (no `--detach`), and script arguments cannot contain
double quotes on Windows.

Running `Rscript`, `python` or a shell script **on the host** (**Tier 2**) is only done with
the user's OK, with the package versions recorded next to the output. Never run it without
asking the user. Figures made this way are reviewed by the `figure-maker` agent from PNG
renders before they are reported.

## Subagents

The agent definitions are in `.claude\agents\` of this folder, so an agent working here has them.
The same kit is in `agent\.claude\agents\`, which `init` copies into every project. The delegation
table and the mandatory triggers are in `CLAUDE.md` of this folder (and in `agent\CLAUDE.md` for
projects). Use `input-preparer` for the user's messy sample, library and contrast sheets (A.4),
`scout` for reading, `log-inspector` after runs (A.8), `results-reporter` for A.7,
`statistician` for statistics.

## Hard rules

1. **Don't edit this pipeline folder** for a user's analysis. Per-project
   changes go in the project (`config.yaml`, `utils\` script overrides — see
   the project's `AGENTS.md`).
2. **Only start the pipeline with the launcher's `run`**, and run the project's
   own analysis scripts with its `script` command. Don't run `docker run`,
   `nextflow`, `bowtie2`, `samtools`, `macs3`, `seacr`, `trimmomatic`, or the `bin\` scripts yourself.
   Commands are taken from `run_cutandrun.ps1` only.
3. **The raw FASTQs and BAMs are the user's data.** They stay where the user keeps them
   (the paths in `samplesheet.csv` point to them). Never move, delete, modify or overwrite
   them, and don't change those paths unless the user asks; if the files are moved, update
   the paths in the CSV and `config.yaml`.
4. **Ask before** deleting anything (`results\`, the Docker work volume, `metadata\work\` and
   `metadata\.nextflow\` via `clean`), installing software, editing `.wslconfig` (or answering yes to the launcher's
   memory question), changing a threshold, or renaming the user's files.
5. **Don't hand-edit `results\`**; it is regenerated. Published intermediate files are
   inside `intermediate_files\`.
6. **Never assume the experiment.** Ask the user (step 2), record it in the
   project's `CONTEXT.md`, and flag data that contradict it instead of
   reinterpreting the design.
7. **Never assume a path.** Reference files and FASTQ locations come from the user.
8. **No run without the user's OK.** Starting a run needs the user's explicit OK, after the
   inputs have been checked with them.
