# AGENTS.md — ATAC-seq pipeline, Docker version for **Mac** (Rephill Lab)

This folder is the ATAC-seq pipeline (Nextflow: trimming, bowtie2 alignment, filtering, duplicate removal, blacklist
filtering, QC, TMM scale factors, bigWig tracks, MACS3 peaks and merges, TSS enrichment, FRiP, MultiQC), packaged to
run on a labmate's own Mac inside one Docker image. The user is usually a **wet-lab scientist**: do the
technical work for them, explain in plain language, and ask them only for what you cannot know (their samples, their
files, their experiment).

First time on this computer (Docker not set up, or never checked)? The lab's shared Docker setup is `SETUP.md` in this
folder (one-time, shared by the lab's Docker pipelines). Read it once; the steps in A.1 below are this pipeline's part of the
setup.

The task here is to set up, run and check an analysis (section A, the runbook below). You never edit anything in this
pipeline folder for that. The pipeline code itself (`main.nf`, `nextflow.config`, `conf/`, `modules/`, `bin/`, `docker/`)
is not meant to be edited during an analysis: if a change seems needed, tell the user and ask them to contact the lab's
maintainers.

If the current folder is a **project** folder (it has `config.yaml`, `samplesheet.csv`, `results/`), read that project's
own `AGENTS.md` instead.

## How to call the launcher

Everything goes through one command, `run_atacseq.sh`, in this folder (full path:
`<this folder>/run_atacseq.sh`). From the agent's shells, call it with `bash` and its full path. That works even if the file
lost its executable bit, and it doesn't depend on the user's PATH:

```bash
bash "/Users/<user>/ATACSEQ/MAC/run_atacseq.sh" status
```

| Command | What it does |
|---|---|
| `run_atacseq.sh setup --yes` | one-time per Mac: Docker Desktop check/install (with Homebrew, after asking), Rosetta 2 on Apple chips, pull the pipeline's Docker image (`docker/image.txt`), print the tool versions, add the command to PATH (`~/.zshrc`) |
| `run_atacseq.sh init "<project folder>"` | create a project folder: `original_data/fastq/`, `metadata/`, `results/`, `utils/`, `examples/` (CSV templates), `config.yaml`, launchers, `.gitattributes`, git when installed |
| `run_atacseq.sh run --detach` (in the project folder) | check the input files and the disk, set the memory budget, start the pipeline in Docker, return immediately |
| `run_atacseq.sh status` (in the project folder) | running / succeeded / failed, the failed step, the last log lines, the Docker work volume's name and size |
| `run_atacseq.sh logs` | follow the live log (blocks until you stop watching; the run goes on). Avoid in an agent: use `status` |
| `run_atacseq.sh stop` | stop a run (finished steps are kept; the next run resumes) |
| `run_atacseq.sh clean --yes` | delete Nextflow's work folder and resume cache to free disk: the project's Docker work volume, and `metadata/work/` and `metadata/.nextflow/` if present. The next run starts from scratch; results are kept. Without `--yes` it asks first; ask the user first |
| `run_atacseq.sh script utils/<name>.<R\|py\|sh> [args]` (in the project folder) | run the project's own analysis script in the pipeline's image (foreground; log in `metadata/script_logs/`); see "Analyses beyond the pipeline" below |
| `run_atacseq.sh help` | the command list |

Also in the launcher, for maintainers and tests (not for users): `mem-plan` (the Docker memory plan for hypothetical
inputs, without Docker), `scan-paths` (the folders a run would mount, from `config.yaml` and the CSVs), and `selftest` (the
pure helper checks of the launcher, no Docker needed).

Options: `--detach` (run: start and return), `--yes` or `-y` (answer yes to every question: only after the user has agreed
to what is being asked), `--pause` (wait for Enter before closing; used by the double-click file). `run CONFIG.yaml` runs with
another config file.

Every project folder also gets its own launcher, `utils/run_atacseq.sh` (same commands, always for that folder; agents
run `bash utils/run_atacseq.sh` from the project folder), and `RUN_ANALYSIS.command` (double-click for humans: it waits for
Enter at the end, so **agents must not call `RUN_ANALYSIS.command`**).

Long commands: `setup` (image pull) and a run can take a long time. Run `setup` in the background or with a long timeout,
and always start runs with `--detach`, then poll `status` every few minutes.

## A. Runbook: from "here are my FASTQs" to results

### 1. One-time Mac setup

The steps below are this pipeline's part of the setup (see the note at the top about `SETUP.md`).

1. Check: `docker --version` and `docker info`. If both work, go to step 3. (Docker Desktop's command-line tool may be in
   `/Applications/Docker.app/Contents/Resources/bin` or `~/.docker/bin`; the launcher looks there.)
2. `bash "<this folder>/run_atacseq.sh" setup --yes`. What can happen:
   - Docker Desktop not installed → with Homebrew, the launcher asks and then runs `brew install --cask docker` (it may ask
     for the user's password: tell them). Without Homebrew it stops and prints the download page. Then the user must **open
     Docker from Applications once** and accept its terms (no account needed: "Skip"). Explain this plainly, wait for them,
     and re-run `setup`.
   - Docker installed but not running → it runs `open -a Docker` and waits up to 4 minutes.
   - Apple chip (M1 to M4): the image is `linux/amd64` and runs through Rosetta 2. `setup` offers to install Rosetta
     (`softwareupdate --install-rosetta`); without it Docker uses much slower emulation. In Docker Desktop → Settings → General,
     "Use Rosetta for x86_64/amd64 emulation on Apple Silicon" should be on; without it everything is several times slower.
3. `setup` pulls the image (`ghcr.io/genesig/omics`, tag in `docker/image.txt`, always run as `linux/amd64`) and prints the
   versions of Nextflow, R and the other tools in the image. Re-running `setup` is harmless. A run also pulls the image if it
   is missing, so a forgotten `setup` only costs time.
4. Resources: every run gets one budget, enforced by Docker as a hard limit. It is the smallest of: 75% of the Mac's RAM
   (always leaving at least 3 GB to macOS), all logical cores but one (two from 8 cores up), what Docker Desktop may use (its
   Memory setting, less about 0.5 GB for Docker's own system), and the most this pipeline can put to use (14 GB,
   `docker/requirements.conf` `max_useful_memory_gb`). `max_cpus` and `max_memory_gb` in `config.yaml` can only lower it. The
   run prints the budget it chose ("Uses at most N CPU cores and X memory"). The pipeline needs at least 5 GB
   (`min_memory_gb` in `docker/requirements.conf`) and refuses to start below that, with an explanation. These figures are
   placeholders until a pilot run on a real Mac has been measured.
5. **Docker Desktop's memory.** The launcher never changes Docker Desktop's settings. It only recommends a Memory value for
   analyses, by the Mac's RAM: under 12 GB → 6 GB; 12 to under 24 GB → 8 GB; 24 GB or more → 16 GB. The user sets it in Docker
   Desktop → Settings → Resources → Memory, then clicks Apply & restart. Changing it is the user's decision: **ask the user
   first**. If Docker's memory holds the budget back, the run says so in a note (with the setting to change); the run stops
   only when the budget is below the 5 GB minimum.

### 2. Gather what you need from the user (ask; never guess)

- **A name for the project** → `project_name` (letters, digits, `_`; used as a prefix on project-level outputs) and the
  project folder's name.
- **Where the FASTQ files are.** This pipeline is **paired-end only**: each library needs `fastq_1` and `fastq_2` (absolute
  paths, one pair per library). If the data are single-end, say so before going further. The FASTQs are named in
  `samplesheet.csv`. A library that already has a **final** BAM (filtered, duplicate-removed, blacklist-filtered) can instead
  give that BAM in the `bam` column, with `fastq_1` and `fastq_2` empty; ask the user whether the BAM is final before using it.
  There is no FASTQ folder mode (`fastq_dir` stays `null`).
- **The reference files** (ask for the paths, never guess them): the bowtie2 index prefix of the genome (`genome_index`),
  `blacklist`, `tss_bed` (required when TSS enrichment runs), `effective_genome_size` (MACS3 `-g`, required when
  peak calling runs; 3100000000 for T2T-CHM13), and optionally `chrom_sizes` (recorded only), `contamination_index` (a
  bowtie2 index of contaminant genomes; null skips CONTAMINATION). The lab's Macs keep references under `~/REFERENCES`; that
  is an example layout, not a requirement.
- **The libraries**: names, which libraries are the replicates of one condition (their `tmm_group`), which are pooled for
  one peak call (their `peak_group`), and which to leave out (`exclude`). Sequencing-core sheets often hold these; ask for
  the sheet.

If they give Excel or Numbers files: ask them to export CSV UTF-8 (Excel: *Save As → CSV UTF-8*; Numbers: *File → Export To →
CSV*), or convert them inside the pipeline image. Check the result opens as plain comma-separated text.

**Then ask about the experiment itself (ask; never assume).** Only the user knows how it was done; ask even when it seems
obvious:

- **The question:** what the experiment is meant to show (accessibility at a set of loci, a change after a treatment or a
  knockout, a difference between cell types), and which genes or loci they expect to open or close (positive and negative
  control loci).
- **The assay:** ATAC-seq (Tn5 transposition), which kit or protocol, the input amount, and the sequencing (paired-end, read
  length, which run(s) and lanes, whether the run also holds other people's samples).
- **The biology:** cell type or model; the treatments, doses and timepoints; the genotype.
- **Each library:** its treatment, which libraries are biological replicates and which are technical, which libraries are
  normalised together (`tmm_group`), which are pooled for peak calling (`peak_group`), and anything that went wrong (low cell
  numbers, a failed library, a contamination, a re-sequenced sample).
- **The normalisation:** whether a global change in accessibility is expected (a treatment that opens or closes most
  regions would move the TMM factors); this version normalises each `tmm_group` with TMM and has no spike-in.
- **The peaks:** whether the peaks should be called per group (as configured) or merged across groups (`peak_merges.csv`),
  and what the expected peak widths are (narrow peaks are the default; `macs3_broad` is off).
- **The Tn5 shift:** whether the user wants the Tn5 offset correction (`tn5_shift`, off by default), and why.
- **Only if the downstream analyses are wanted:** the comparisons (case and control values of which samplesheet column, which
  libraries are blocks such as cell lines, which peak groups each side uses, which lines are pooled), the fixed dispersion for
  a comparison without replicates, the motif targets, and the thresholds. These are scientific decisions: ask for each one.

- Ask in short batches (3-6 questions at a time), in plain language, and say why each answer matters. Keep asking until the
  answers are specific; "the usual" or "like last time" is not an answer, so ask what it was.
- **Never infer the design from names** (`WT`, `rep2`, `D0`, `ctrl`, `hi`/`lo`, `high`/`low`) or from what is "typical" for
  ATAC-seq. Ask what each name means.
- Use only what the user confirmed. Where an answer is missing, say what you would do instead and wait for an explicit OK;
  never fill the gap silently, and never pick the replicates, groups or comparisons yourself.
- **Write the answers down** in `CONTEXT.md` in the project folder (dated, in the user's words, with open questions marked as
  open). Read it at the start of every session; ask the user to confirm it is still right before you rely on it, and add to
  it as you learn more. It is tracked by git with the inputs.
- **When the data seem to contradict the context, flag it; don't reinterpret.** Examples: a low TSS enrichment
  (`results/qc/tss_enrichment/`); a high chrM fraction (`results/qc/chrom_fraction/`); a fragment-size distribution without
  the expected nucleosome pattern (`results/qc/fragment_sizes/`); a replicate that differs strongly from its partner; a
  "global loss" treatment where the TMM factors do not move; a library with a low alignment rate
  (`results/qc/alignment/`). Show the user the numbers and the file they come from, say what the data would be consistent
  with, and ask. The user's description stands until they change it; don't rename samples, drop them or change groups on
  your own.
- Subagents cannot talk to the user: put the confirmed context (or the path to `CONTEXT.md`) in every brief, and pass their
  "Questions for the user" back to the user rather than answering them yourself.

### 3. Create the project folder

Put projects **outside** the pipeline's folder, e.g. `~/atac_projects/<project_name>` (on the internal disk if possible; not
in an iCloud-synced Desktop or Documents folder; see the disk note in step 6). Use a new folder: if the one you picked already
exists and holds anything, list it and ask the user before using it (`init` never overwrites files, but it adds its own files
and a git repository to whatever is there). Then:

```bash
bash "<this folder>/run_atacseq.sh" init "$HOME/atac_projects/<project_name>"
```

`init` writes the project's `config.yaml`, copies the CSV templates into `examples/` (never overwriting them), copies the
project's `AGENTS.md`/`CLAUDE.md` and the agent kit (`.claude/agents/`) when this folder has them (`agent/` in this folder),
and rewrites the two launchers. The CSVs the project uses are copies you make from `examples/` into the project folder
itself.

Immediately after `init`, **read the new project's `CLAUDE.md` and `AGENTS.md` and follow them from then on** (their
delegation table and triggers included): users usually ask you to set the project up before running, so the project's
`.claude/agents/` and instructions appear mid-session. If Claude Code was started in this pipeline folder, its
`.claude/agents/` is the same kit, so the agents are available now. If they are not listed as available subagent types,
either spawn a general-purpose agent whose prompt is the agent file's body (after the frontmatter) plus the brief, with the
model from the file's frontmatter, or ask the user to restart Claude Code in the project folder (recommended for every
later session: start Claude Code inside the project folder, so it loads that project's `CLAUDE.md` and agents).

### 4. Write the input tables (in the project folder)

The full column rules, the enforced checks and the examples are in the project's `AGENTS.md` ("Input files"); the templates
are in `examples/`. Summary of each table:

- **`samplesheet.csv`** (one row per library; named by `samplesheet:` in `config.yaml`):
  `sample_id,fastq_1,fastq_2,condition,replicate,tmm_group,tmm_group_legacy,peak_group,exclude,bam`.
  `sample_id` is unique and uses letters, digits, `_`, `.`, `-` (it names every output file). `fastq_1` and `fastq_2` are
  both required, unless the row has a `bam` (a final BAM ending in `.bam`, existing; then both FASTQ fields are empty; never
  FASTQ and a BAM on one row). `condition` and `replicate` are labels for the QC tables only. `tmm_group` is the TMM normalisation group
  (exact name; blank = the library is not normalised: no bigWig, no TSS enrichment; a group needs at least 2 active
  libraries). `tmm_group_legacy` is a second group column, used only when `tmm_group_column: tmm_group_legacy`. `peak_group`
  is the MACS3 peak group (replicates are pooled; blank = no peak calling and no FRiP for that library). `exclude`
  (`true`/`false`) drops a row from the whole run.
- **`peak_merges.csv`** (optional; `peak_merges:` in `config.yaml`; one row per merged BED): `name,peak_groups`, with
  `peak_groups` a quoted comma list of `peak_group` values. It needs `run_macs3`.
- **Downstream CSVs** (optional; read only when `run_differential`, `run_annotation` or `run_motifs` is on, or when set):
  `comparisons.csv` (`name,case,control,group_column,block_column,case_peaks,control_peaks,pooled,fixed_dispersion`; `case` and
  `control` are quoted comma lists of values of `group_column`, default `condition`), `peak_files.csv` (`name,path`, existing
  peak files), `diff_tables.csv` (`name,regions,path`, existing differential tables, not recomputed) and `motif_targets.csv`
  (`name,type,target,comparison,label,min_peaks`; `type` `gene` or `regions`). Each is named in `config.yaml` (`comparisons:`,
  `peak_files:`, `diff_tables:`, `motif_targets:`). The rules and the startup checks are in the README, "Downstream analyses".

Start from `examples/*.csv`; save as UTF-8.

### 5. Edit `config.yaml`

Every key has a default in the pipeline's `nextflow.config`; list only what you change. The example is `examples/config.yaml`.
The key choices:

- `project_name`, `samplesheet`, `peak_merges`.
- **Reference files**: `genome_index` (bowtie2 index prefix), `blacklist`, `tss_bed`, `gtf` (recorded only), `chrom_sizes`
  (recorded only), `effective_genome_size`, `contamination_index`. Write Mac paths in full (`/Volumes/SeqDrive/REFERENCES/...`,
  `~/REFERENCES/...`) or relative to the project folder. The launcher warns when a file is missing or when `genome_index` has no
  `.bt2` files next to it.
- **Input files stay where they are.** The FASTQs are never moved, copied, renamed or deleted; the top-level folder holding
  them (`/Users`, or `/Volumes/<name>`) is mounted read-only into Docker at the same path. Tell the user where the files are
  read from and that they stay there. Keep the paths stable (if the files are moved, update the CSV and `config.yaml`).
  Refused: network paths (`smb://...`, `//server/share`); mount the share in Finder and use its `/Volumes/<name>/...` path
  (copy the files to a local or external drive if Docker cannot read it).
- **Output folders must be relative**: `outdir: ./results`, `metadata_dir: ./metadata`, `intermediate_dir: ./intermediate_files`.
  The launcher refuses an absolute path there (such as `/Users/...`).
- **Duplicates**: `remove_duplicates: true` (sambamba markdup, duplicates removed; `false` marks them and keeps them). Setting it
  to `false` changes the library sizes and the TMM factors: ask first.
- **Alignment**: `bowtie2_mode: local` (or `end_to_end`); `bowtie2_minins: 10` and `bowtie2_maxins: 2000` (bowtie2 `-I` and
  `-X`; changing either re-runs alignment and everything downstream: ask first); `bowtie2_args` (extra arguments, empty);
  `align_chunk_pairs: 10000000` (read pairs per ALIGN task; `0` = one ALIGN task per library; the output is the same);
  `alignment_min_pct: 20` (warning threshold only, never fails).
- **Trimming**: `trim: true`, `trim_adapter_set: 'NexteraPE-PE'`, `trim_illuminaclip: '2:30:10:8:True'`,
  `trim_steps: 'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20'`, `trim_minlen: 25`.
- **Filtering**: `mapq: 30` (minimum MAPQ kept; `0` reproduces the legacy filter); `remove_chroms: chrM` (comma list of
  contigs removed; empty = keep all).
- **Tn5 shift**: `tn5_shift: false` (`true` = deepTools `alignmentSieve --ATACshift` before the QC, TMM, bigWig and MACS3
  steps; ask the user first).
- **TMM normalisation**: `tmm_bin_size: 10000`, `tmm_max_frag: auto` (the largest proper-pair fragment of the group, or a
  positive integer), `tmm_group_column: tmm_group`.
- **bigWig tracks**: `bigwig_bin_size: 10`, `bigwig_units: none` (`none`, `cpm` or `rpkm`), `bigwig_extend_reads: true`
  (fragment coverage), `bigwig_center_reads: false`, `bigwig_args`.
- **Peaks (MACS3)**: `run_macs3: true`; `macs3_qvalue: 0.05`, `macs3_min_length: 250`, `macs3_max_gap: 500`,
  `macs3_fe_cutoff: 2.0` (a value of 1 or lower turns it off), `macs3_call_summits: false`, `macs3_broad: false`,
  `macs3_broad_cutoff: 0.1`, `macs3_args`. Changing a threshold is a scientific decision: tell the user and get confirmation
  first.
- **QC**: `run_tss_enrichment: true` (needs `tss_bed` when some library has a `tmm_group`); `tss_flank: 2000`; `run_fastqc: true`; `run_multiqc: true`.
- **Step switches**: `run_fastqc`, `trim`, `run_tss_enrichment`, `run_macs3`, `run_multiqc`, `tn5_shift`, `remove_duplicates`.
- **Downstream switches** (all off by default): `run_differential` (DIFF_ACCESS), `run_annotation` (DIFF_ANNOTATE, needs `gtf`),
  `run_motifs` (MOTIF_ENRICHMENT and MOTIF_COMMON). Their CSVs: `comparisons`, `peak_files`, `diff_tables`, `motif_targets`. Their
  values (`diff_*`, `annotation_*`, `motif_*`, `species`, `genome_build`, `bsgenome`, `org_db`, `plot_font`, `plot_format`) are in
  the README table. The statistical ones (`diff_lfc`, `diff_fdr`, `diff_table`, `diff_dedup`, `diff_qlf_min_libs`,
  `diff_fixed_dispersion`, `motif_lfc`, `motif_regions`, `annotation_fdr`) are scientific decisions: confirm them first.
- **Outputs and bookkeeping**: `publish_intermediates: false` (`true` also publishes the trimmed FASTQs, aligned and filtered
  BAMs to `intermediate_dir`); `publish_mode: auto` (hard links from the work folder when the disk allows, else copies;
  `copy`, `link`, `symlink`, `rellink` also accepted; the launcher refuses link modes with the work volume); `keep_work: true`
  (keep the work folder and resume cache after a successful run, so the next run resumes; `false` deletes them);
  `work_volume: true` (the work folder and cache in a Docker volume on Docker's disk; `false` keeps them in `metadata/work/`).
- **Resources (optional)**: `max_cpus` and `max_memory_gb` lower the run's budget (A.1, step 4); use them only if the user
  needs the Mac to stay responsive.
- Leave `fastq_dir` as `null`: this Mac version does not read a FASTQ folder.
- Output keys that are not in this version: there is no QC report HTML from this pipeline; the MultiQC report is the
  overview (`results/qc/multiqc/multiqc_report.html`).

### 6. Check, then run

Check before the first run (the launcher also checks, but catching it now is faster): headers exactly right; `sample_id` unique;
FASTQ paths exist (the launcher names any file it cannot find); every `tmm_group` has at least 2 active libraries; the
`peak_group` values in `peak_merges.csv` exist; a `bam` row has no FASTQ and its file ends in `.bam` and exists; for the
downstream analyses, every comparison's `case` and `control` values occur in the samplesheet, no sample is on both sides, and
the peak groups named in the comparisons exist; the groups are the ones the user confirmed.

```bash
cd "<project folder>"
bash utils/run_atacseq.sh run --detach
bash utils/run_atacseq.sh status        # every few minutes
```

The run first checks the inputs and prints the disk check: the free space on Docker's disk (the disk that holds Docker Desktop's
disk image, which holds the work volume) and on the project's disk (which holds the results), against the total size of the
FASTQ files in the samplesheet times two rules in `docker/requirements.conf`: `disk_factor=5.2` (work files and results
together; used when the work folder is in the project folder) and `disk_results_factor=0.4` (results only; used for the
project's disk when the work volume holds the work files, where Docker's disk is checked against about 4.8 times the input).
It warns below that on Docker's disk, and stops below 20 GB free on the project's disk (`disk_hard_min_gb`). To estimate by
hand, add up the sizes of the FASTQ files named in `samplesheet.csv` and multiply by 5.2 (about 4.8 for Docker's disk and 0.4
for the project's disk, with the work volume). **These multipliers were measured on one full-depth run of two libraries; they remain estimates until
the first full ATAC-seq run gives the real sizes**; a dataset can need about 100 GB or more of working space. Free space can
be given back with `clean` (ask first) or with `keep_work: false`, which deletes the work files after a successful run.

Timing on a Mac has not been measured for this port yet. Poll `status` every few minutes. Docker's disk image should be on the
internal SSD (or a disk with room): the work files are there. The project folder can be on an external or hard disk (the
launcher warns): with the work volume only the results are written there (and `intermediate_files/` with
`publish_intermediates: true`). The Mac must stay on and awake (plugged in, lid open): the launcher keeps it from idle sleep
while the run is going, but if it sleeps anyway the run pauses; if Docker restarts, just run again (it resumes).

Nextflow's work folder and its resume cache are kept by default (`keep_work: true`), so a re-run resumes and only redoes what
changed (e.g. after changing a threshold). Failed or stopped runs always keep them. With `work_volume: true` (the default)
they are in the project's Docker volume, on Docker Desktop's own disk: the work folder is `/nfwork/work` and the cache
`/nfwork/.nextflow` inside the container. The volume is not visible in Finder; `status` shows its name and size. Files
in `results/` are copies of the work files (`publish_mode: auto`). With `work_volume: false` they are `metadata/work/` and
`metadata/.nextflow/` in the project folder (`metadata/` is git-ignored). In either case never edit a published file in place
(replacing or deleting is fine). The trimmed FASTQs, aligned BAMs and filtered BAMs are in `intermediate_files/` only with
`publish_intermediates: true` (the default is `false`: they stay in the work folder).

### 7. Report results

`status` says `succeeded` → point the user to:
- `results/qc/multiqc/multiqc_report.html` — the MultiQC overview (start here);
- `results/qc/alignment/alignment_summary.tsv` — reads and alignment rates per library, with the warning flags;
- `results/scale_factors/<tmm_group>.scale_factors.tsv` — the TMM table (`bigwig_multiplier`, `edger_norm_factor`, `qc_pass`,
  `note`); check it before interpreting any normalised track or peak set; `<tmm_group>.max_frag.txt` gives the fragment size
  used;
- `results/peaks/macs3/<group>...narrowPeak` — the peaks per peak group; `results/peaks/merged/<name>.merged.bed` if merges were
  set;
- `results/bigwig/` — the tracks (`<sample>.<method>.<units>.bw`); `results/bam/final/` — the final BAMs;
- with the downstream switches on: `results/differential/` and `results/differential/bed/` (the differential tables and BEDs; the
  up and down BEDs have no FDR filter), `results/annotation/` (the annotation tables) and `results/motif/` (the motif tables; an
  analysis with `<name>.skipped.txt` was not run). Report them only with the limits in `ANALYSIS.md` ("Differential accessibility",
  "Motif enrichment"), and check the MA plot and the factors (`results/differential/<name>.norm_factors.tsv`) first;
- `metadata/pipeline_info/execution_report_<timestamp>.html` — how long each step took and how much memory it used.

Summarise from `results/qc/` (the alignment table, duplicate and chrM counts, fragment sizes, TSS enrichment, FRiP and peak
statistics) and the MultiQC report. Check the quality first: a library with a low alignment rate, a high chrM fraction, a low
TSS enrichment or a low FRiP is flagged, not explained away. For statistics or comparisons with earlier numbers, use the
`statistician` and `output-verifier` agents (see "Subagents").

**Fragment-size QC.** `results/qc/fragment_sizes/` has the fragment-size table and distribution and the plots of
`qc_fragment_sizes_plot.R`. The expected ATAC pattern (short nucleosome-free fragments and nucleosome-sized fragments) is
something to check against the user's expectation; this pipeline does not decide it.

### 8. When a run fails

`bash utils/run_atacseq.sh status` says which step failed (`FAILED: <step> (exit N)`) and prints the last log lines. The files to
read, in order:

1. `metadata/nf_logs/last_run.status`: `running`, `succeeded`, `failed (exit N)`, `failed to start`, `stopped by the user`,
   or `interrupted` (Docker or the Mac stopped mid-run; run again to resume).
2. `metadata/nf_logs/nextflow_<timestamp>.log`: everything Nextflow printed. Startup errors (bad samplesheet, missing file, an
   unknown column, a TMM group with one library, a typo in a key) are printed at its end.
3. `metadata/failed_tasks/<timestamp>/<step>/`: the failed step's own files (`.command.sh` = what ran; `.command.err` /
   `.command.log` = what it said; `.exitcode`).
4. `metadata/pipeline_info/execution_trace_<timestamp>.txt`: one row per step (status, exit code, `peak_rss` = the memory it
   used, `realtime`). Use it to see whether a step ran out of memory or disk. Use the newest file (one per launch).

Common cases (the project's `AGENTS.md`, "When a run fails", lists more):
- **Exit 137, 140 or 143**: out of memory or killed. These are retried automatically with more memory, up to the budget. If
  the step still fails, check the budget in the run banner and the step's `peak_rss`. Raising the budget is the user's
  decision (`max_memory_gb`, or a bigger Docker Desktop memory setting): ask first.
- **Exit 3**: a known data problem, not retried; the message says what to fix. Typical: a bowtie2 index that is not there; an
  empty or missing BAM; unequal or truncated FASTQ pairs (COUNT_PAIRS); a chunk log that cannot be merged (ALIGN_MERGE);
  orphan reads after blacklist filtering (BAM_FINALIZE); no contig left after `remove_chroms`; a scale factor missing or NaN
  (BIGWIG); a MACS3 target BAM missing or empty; a file listed in a peak merge missing; a TMM group with fewer than 2 libraries
  or no proper pairs (NORM_TMM). Downstream: a comparison side with no BAM, a peak file or a table with no rows, promoters without
  a TxDb, or a BSgenome package that is not installed (DIFF_ANNOTATE, MOTIF_ENRICHMENT). Downstream startup errors name the comparison, peak group, table or value at fault.
- **Memory on this Mac:** DIFF_ACCESS asks for 16 GB and MOTIF_ENRICHMENT for 32 GB (cluster estimates); here the request is
  capped to the run's budget (at most 14 GB). A comparison with many bins or a full-table motif analysis can run out of memory
  here. Say so to the user and offer the cluster; do not change `max_memory_gb` or Docker Desktop's memory setting without asking.
- **Other failures** are retried up to twice by the pipeline, then stop the run.
- **Disk full**: the run stops and can be resumed after space is made. Check the free space on Docker's disk and on the
  project's disk (`status` and the disk check in `run`); `clean` (ask first) deletes the work volume (or the work folder with
  `work_volume: false`).
- **Startup errors**: fix the named CSV or key in the project folder, then run again.

Fix the input or config, then run again: `bash utils/run_atacseq.sh run --detach`. Finished steps are reused. If it looks like a
bug in the pipeline itself, tell the user and stop: don't patch this folder, and ask the user to contact the lab's
maintainers.

### 9. Analyses beyond the pipeline

The user may want a plot or test the pipeline doesn't do. That is allowed, through the project's `script` command:
`bash utils/run_atacseq.sh script utils/<name>.<R|py|sh>` runs the script in the pipeline's image (**Tier 1**), and the log goes
to `metadata/script_logs/`. The script lives in `utils/`, writes only to a new `results/<name>/` folder, and is designed,
written, run and reviewed as described in "Analyses beyond the pipeline" in the project's `AGENTS.md` (`init` copies it there).
The `script` command is foreground only (no `--detach`), and script arguments cannot contain double quotes.

Running `Rscript`, `python` or a shell script **on the host** (**Tier 2**) is only done with the user's OK, with the package
versions recorded next to the output. Never run it without asking the user. Figures made this way are reviewed by the
`figure-maker` agent from PNG renders before they are reported.

## Subagents

The agent definitions are in `.claude/agents/` of this folder, so an agent working here has them. The same kit is in
`agent/.claude/agents/`, which `init` copies into every project. The delegation table and the mandatory triggers are in
`CLAUDE.md` of this folder (and in `agent/CLAUDE.md` for projects). Use `input-preparer` for the user's messy sample sheets
(A.2 and A.4), `scout` for reading, `log-inspector` after runs (A.8), `results-reporter` for A.7, `statistician` for statistics.

## Hard rules

1. **Don't edit this pipeline folder** for a user's analysis. Per-project changes go in the project (`config.yaml`, `utils/`
   script overrides: see the project's `AGENTS.md`).
2. **Only start the pipeline with the launcher's `run`**, and run the project's own analysis scripts with its `script`
   command. Don't run `docker run`, `nextflow`, `bowtie2`, `samtools`, `sambamba`, `bedtools`, `macs3`, `trimmomatic`, or the
   `bin/` scripts yourself. Commands are taken from `run_atacseq.sh` only.
3. **The raw FASTQs are the user's data.** They stay where the user keeps them (the paths in `samplesheet.csv` point to them).
   Never move, delete, modify or overwrite them, and don't change those paths unless the user asks; if the files are moved,
   update the paths in the CSV and `config.yaml`.
4. **Ask before** deleting anything (`results/`, the Docker work volume, `metadata/work/` and `metadata/.nextflow/` via `clean`),
   installing software (Docker Desktop, Homebrew packages, Rosetta), changing Docker Desktop's settings, changing a threshold or
   the Tn5 shift, or renaming the user's files.
5. **Don't hand-edit `results/`**; it is regenerated. Published intermediate files are inside `intermediate_files/`.
6. **Never assume the experiment.** Ask the user (step 2), record it in the project's `CONTEXT.md`, and flag data that
   contradict it instead of reinterpreting the design.
7. **Never assume a path.** Reference files and FASTQ locations come from the user.
8. **No run without the user's OK.** Starting a run needs the user's explicit OK, after the inputs have been checked with them.
