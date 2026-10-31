# CUT&RUN / chromatin profiling pipeline for Windows (Docker)

The Richard Phillips Lab's CUT&RUN and CUT&Tag pipeline: trimming, alignment to the
target and spike-in genome, duplicate marking and removal, QC report, spike-in or
ChIPseqSpikeInFree normalisation (with TMM, CPM and spike-in diagnostics), bigWig tracks, csaw, MACS3 and
SEACR peak calling (MACS3 and SEACR need IgG controls), and optional region classification, heatmaps, genome-browser tracks, peak annotation,
BED intersections. It is the same pipeline as on the PMACS
cluster (same steps, same tools, same outputs; the last four BAM steps run as one process here, see "Differences from the cluster"), packaged in one Docker image, so you
don't install any bioinformatics software yourself.

> **First time on this computer?** Follow `SETUP.md` in this folder once (Docker
> Desktop, WSL 2, where Docker's data lives, memory). Then come back here.
>
> The agent instructions for this pipeline (`AGENTS.md`, `CLAUDE.md` and the agent
> definitions in `.claude\agents\`) are copied into each new project by `init`. The
> steps below are complete on their own.

## What you need

- Windows 10 (22H2) or 11, 64-bit, with hardware virtualization enabled (on by
  default on most computers).
- **8 GB of memory or more** (16 GB or more is more comfortable). With 8 GB the
  analysis gets about 5 GB; with less, the launcher refuses to start a run and
  says why.
- Free disk space: on the drive that holds Docker Desktop's disk image (`C:` by default),
  **about 4.6 times the total size of your FASTQ files** for the work files (with the default
  `work_volume: true`; see "Disk space" below; an estimate, pending the next run), plus room for
  the Docker image (several GB, once); on the drive that holds your project, about 1 times the
  FASTQ total for the results, and at least 20 GB free. The raw FASTQs are not copied.
- Internet the first time (to download Docker and the pipeline's image).
- A reference set: a bowtie2 index of the combined target + spike-in genome, the
  chromosome sizes, a blacklist BED, a gene annotation GTF and a TSS BED. The lab's
  PC keeps them in `D:\REFERENCES` (see "Where data and projects go" in `SETUP.md`);
  any other location works as long as you write its path in `config.yaml`.

## One-time setup (per computer)

1. Get the pipeline (once), as a sparse clone so that only this folder is downloaded:
   ```
   git clone --filter=blob:none --sparse <the lab's GitHub URL>
   cd PIPELINES
   git sparse-checkout set CHROMATIN_PROFILING/WINDOWS
   cd CHROMATIN_PROFILING\WINDOWS
   ```
   Then work in the `CHROMATIN_PROFILING\WINDOWS` folder; `git pull` later updates it.
   (The folder name of the clone is whatever the lab's repository is called.)
2. Open the `CHROMATIN_PROFILING\WINDOWS` folder and run, in PowerShell or Command Prompt:
   ```
   .\run_cutandrun.bat setup
   ```
   - If Docker Desktop isn't installed, it offers to install it (Windows asks for
     administrator permission). Then **restart the computer**, open **Docker Desktop**
     from the Start menu once, accept its terms (no account needed: click *Skip*), and
     run `setup` again.
   - It then downloads the pipeline's Docker image (`ghcr.io/genesig/omics`, the tag
     is in `docker\image.txt`): several GB, once. It takes a while; keep the computer
     awake and online.
   - It prints the versions of the main tools, and adds `run_cutandrun` to your PATH,
     so new terminal windows can use it from any folder.

## Analyze a dataset

1. **Create a project folder** (one per dataset, *outside* the `CHROMATIN_PROFILING`
   folder):
   ```
   run_cutandrun init C:\Users\<you>\Documents\chromatin\my_dataset
   ```
   It creates the folder skeleton (`original_data\fastq\`, `metadata\`, `results\`,
   `utils\`), `config.yaml`, `examples\` (starter CSV templates), `RUN_ANALYSIS.bat`
   (double-click) and `utils\run_cutandrun.bat` (the command used below). It also
   writes `.gitattributes` and, when git is installed, starts a git repository in the
   folder with the first commit. It doesn't run anything.
2. **Add your tables** to the project folder (the project root, next to `config.yaml`).
   Copy the templates you need from `examples\`, fill them in (Excel is fine: *Save As
   → CSV UTF-8*):
   - `samplesheet.csv` (required): one row per library. See "Input tables" below.
   - `norm_runs.csv` (required by the default normalisation): one row per normalisation run.
   - `peaksets.csv` (required by the default peak calling, `run_csaw: true`): one row per csaw peak set.
   - `classifications.csv`, `heatmaps.csv`, `tracks.csv`, `intersects.csv`:
     only when you switch on the matching step in `config.yaml`.
3. **Edit `config.yaml`** (Notepad is fine). Set `project_name`, check the table file
   names, and set the reference paths: `genome_index` (the bowtie2 index prefix, e.g.
   `D:\REFERENCES\bowtie2\t2t_ecoli\<prefix>`), `chrom_sizes`, `blacklist`, `gtf` (and
   `tss_bed` only when `gtf` is not set). Keep `effective_genome_size` for your genome and `experiment` (`cutandrun`
   or `cutandtag`, which selects the adapters). Every other key has a default; list
   only the ones you change. Key summary below.
4. **Run it:** double-click **`RUN_ANALYSIS.bat`** in the project folder, or in a
   terminal in that folder: `.\utils\run_cutandrun.bat`.
   - Before it starts, the launcher checks your inputs (paths, the samplesheet, the
     output folders, the free disk space on the project's drive and on Docker's disk) and
     sets this computer's memory for the
     analysis (see "Resources" below; it asks first).
   - The analysis runs inside Docker in the background; the window shows its
     progress. Closing the window does **not** stop it, but keep the computer on and
     awake (plugged in, sleep off) until it finishes.
   - Check on it any time: `.\utils\run_cutandrun.bat status`.
5. **Results** are in `results\`. Start with `results\qc\report\<project_name>_qc_report.html`
   (open in a browser: the alignment, fragment size, library complexity, correlation,
   FRiP, TSS enrichment and peak QC are explained there). The MultiQC summary is in
   `results\qc\multiqc\`. The normalisation factors are in `results\scale_factors\`
   (the primary table per normalisation run, plus `<run>.scale_factors_all.tsv` with every method), the bigWig tracks in `results\bigwig\`, the csaw
   peaks in `results\peaks\csaw\`. The layout is listed in `ANALYSIS.md`.

**If something fails**, `.\utils\run_cutandrun.bat status` says which step failed and
shows the last lines of the log. The files of the failed step are in
`metadata\failed_tasks\<timestamp>\`. Fix the input and run again: the steps that
already finished are reused, not redone.

**Re-running a finished project** only redoes what changed (for example after you
switch on classification, or change a plot setting), because Nextflow's work folder and
its resume cache are kept by default (`keep_work: true`). With the default
`work_volume: true` they are in a Docker volume named after the project, on Docker
Desktop's own disk (not in the project folder; the volume is not visible in File
Explorer). With `work_volume: false` they are `metadata\work\` and `metadata\.nextflow\`
in the project folder, as before. To free disk space, run
`.\utils\run_cutandrun.bat clean` (asks first; results are kept) or set
`keep_work: false` in `config.yaml` (deletes the work files after a successful run).
The next run then recomputes everything. A failed or stopped run always keeps them.

**Changing `work_volume`** starts the run from scratch, because the resume state is in the
other place: a project whose `metadata\work\` already has content is asked first (`--yes`
answers yes). To resume the earlier run, set `work_volume: false`. The old `metadata\work\`
is left alone (`clean` deletes it). Switching back and forth loses the resume each way.

**Your own analysis** (a plot or test the pipeline doesn't do): put the script in
`utils\` (`.R`, `.py` or `.sh`) and run it with
`.\utils\run_cutandrun.bat script utils\my_analysis.R`. It runs in the pipeline's
software environment and its log goes to `metadata\script_logs\`. See "Customizing one
project's analysis" below.

## All commands

Run from the project folder (`.\utils\run_cutandrun.bat ...`) or from anywhere after
setup (`run_cutandrun ...`, acts on the current folder):

| Command | |
|---|---|
| `run_cutandrun setup` | one-time: Docker Desktop (offers to install it) + the pipeline's Docker image |
| `run_cutandrun init FOLDER` | create a project folder (default: the current folder) |
| `run_cutandrun` (or `run`) | run / resume the analysis of the current folder (`--detach`: start and return, don't follow the log) |
| `run_cutandrun run other_config.yaml` | ... with another config file in the folder (or pass the project folder) |
| `run_cutandrun status` | running? finished? which step failed and why? the work volume's name and size |
| `run_cutandrun logs` | watch the live log (Ctrl+C stops watching, not the analysis) |
| `run_cutandrun stop` | stop the analysis (running again resumes it) |
| `run_cutandrun clean` | delete Nextflow's work folder and resume cache to free disk space: the project's Docker work volume, and `metadata\work\` and `metadata\.nextflow\` if present; asks first (`--yes` skips); the next run starts from scratch, results are kept |
| `run_cutandrun script utils\NAME.R [args]` | run your own script (`.R`, `.py`, `.sh`; inside the project folder) in the pipeline's software environment, in the foreground; log in `metadata\script_logs\` |
| `run_cutandrun help` | the list of commands |

Options: `--yes` answers yes to the launcher's questions, including the Docker memory
question (for scripts and AI agents); `--pause` waits for Enter before the window closes
(used by `RUN_ANALYSIS.bat`); `--no-wsl-resize` (with `run`) leaves Docker Desktop's
memory limit as it is.

## Input tables

All tables live in the project folder root (next to `config.yaml`). Paths in them
follow the rules in "Inputs are never copied" below. Lists inside one cell are
comma-separated and quoted (`"IgG_WT_rep1,K27me3_WT_rep1"`).

### samplesheet.csv (one row per library)

`sample_id,fastq_1,fastq_2,bam,antibody,condition,replicate,control_id,group,norm_run,exclude`

| Column | Meaning |
|---|---|
| `sample_id` | Unique; used verbatim in file names; `^[A-Za-z0-9][A-Za-z0-9_.-]*$` |
| `fastq_1`, `fastq_2` | Paths of the R1 and R2 FASTQ files (paired-end). Give these OR `bam` |
| `bam` | A pre-aligned, final-filtered BAM: skips trimming, alignment and filtering; no spike-in BAM, so it cannot be used in a spike-in normalisation run |
| `antibody` | Mark or protein; `IgG` (any case) marks a control |
| `condition`, `replicate` | Labels, e.g. `DMSO`, `A` |
| `control_id` | `sample_id` of the IgG control of this library; empty for IgG rows. Needed for MACS3, SEACR and log2-ratio tracks |
| `group` | The unit of replication; empty = `<condition>_<antibody>`. Peak sets and classifications refer to groups |
| `norm_run` | A `norm_run` of `norm_runs.csv`; required for every active row unless `normalisation: cpm` |
| `exclude` | `true` drops the row from the run |

### norm_runs.csv (one row per normalisation run)

`norm_run,method,bin_size,cutoff,max_turns,min_first_turn,ref_condition`

- `method` = `chipseqspikeinfree` (`bin_size`, `cutoff` and `max_turns` all required;
  the package defaults are cutoff 1.2 and max_turns 0.99), `spikein` (the three
  numbers left blank), `tmm` (`bin_size` optional, default `tmm_bin_size`) or `cpm`
  (library size only, numbers blank; the only method allowed for a run with one library).
  It is the run's **primary** method: only it scales the bigWigs and the classification.
  The other methods of `norm_methods` (and spike-in when `norm_spikein: true`) are
  computed for comparison in `scale_factors\<run>.scale_factors_all.tsv` and the MA plots.
- `min_first_turn` (`chipseqspikeinfree` only; blank = 0.5, the package default) sets
  the package's minFirstTurn and is written into the scale-factor note. Other methods ignore it.
- `ref_condition` (optional): the reference condition of the run's MA plots. Blank = the
  condition of the `ref_group` of the classifications that use the run (if exactly one), else
  the condition of the run's first library. It must be a condition of the run.

### peaksets.csv (one row per csaw peak set)

`peakset,target_group,control_group,small_filt,large_filt,merge_bp,max_width,fdr,lfc`

- `target_group`, `control_group`: values of the samplesheet `group` column.
- Threshold cells may be blank. Each value resolves in this order: the non-blank
  cell of the row; then `csaw_defaults[<antibody>][<column>]`, where the antibody is
  the samplesheet `antibody` of the row's `target_group`; then the global `csaw_<column>`.
- `max_width` blank = no limit.
- Startup fails if a `target_group` mixes antibodies or a resolved value is not numeric.
  `metadata\peaksets.resolved.tsv` lists every resolved value and where it came from.

**What the FDR and the logFC threshold mean.** `fdr` is the false discovery rate of the
test "target enriched over IgG": a region with FDR at or below `fdr` is enriched over
IgG. `lfc` is a descriptive threshold applied after that test: a region is a peak only
if it passes the FDR and the log2 fold change of its representative window is above
`lfc`. Which window that is follows `csaw_lfc_from` (`combined`, the default, takes
the representative window of the Simes combination; `best` takes the window with the
smallest p-value). `lfc` is an effect-size filter, not a tested claim: the FDR does not
control the probability that a peak is more than 2^`lfc` fold enriched. The `best`
window is chosen on the same p-values it is then filtered on, so its fold change is
the more biased of the two. There is no `glmTreat` option; the default background
normalisation is `csaw_norm_bins: background`.

### classifications.csv (`run_classify: true`)

`classification,regions,ref_group,treat_group,norm_run,lfc_retained,lfc_lost,fdr`

- `classification`: an id (same pattern as `sample_id`), prefix of every output file.
- `regions`: a csaw peak set name (needs `run_csaw`) or the path of a BED file (BED3
  or more; a 4th column gives the name, else `region_<i>`).
- `ref_group`, `treat_group`: `group` values, reference and treated; must differ.
- `norm_run`: empty = the `norm_run` shared by all reference and treated samples.
- `lfc_retained`, `lfc_lost`, `fdr`: empty = `classify_lfc_retained` (-1),
  `classify_lfc_lost` (-2), `classify_fdr` (0.05).

Rules, with logFC (treated over reference) and FDR (Benjamini-Hochberg over all regions of the row):

- **lost**: `logFC <= lfc_lost` and `FDR < fdr`;
- **retained**: `round(logFC, 1) >= lfc_retained`, or `FDR >= fdr`. A region without a
  significant change counts as retained: "no significant loss", not "shown to be retained";
- **unclassified**: the rest. A region both lost and retained is lost.

Counts use the samples' BAMs with the csaw read settings. The normalisation factors are
the `edger_norm_factor` column of the scale-factor table; they are not re-estimated from
the regions. A group with one library gives a warning: the dispersion then comes from the
fixed-dispersion fallback, and the FDR is not valid. Outputs go to
`results\classify\<classification>\`: `<id>.retained.bed`, `<id>.lost.bed`,
`<id>.unclassified.bed` (BED6, 0-based starts, score = logFC), `<id>.per_region.tsv.gz`,
`<id>.summary.tsv` and the MA plot `<id>.ma.svg`.

### heatmaps.csv (`run_heatmaps: true`)

A sample that is not in the samplesheet (or is excluded / in `skip_groups`) is dropped from a row with its `labels`, `colors`, `z_min` and other per-sample entries, with one WARN per row (a row left with no sample is skipped); `tracks.csv` works the same way.

`heatmap,regions,samples,labels,window_bp,bin_bp,sort_using,cluster_k,z_min,z_max,y_max,colors,engine`,
then the optional `cluster_samples,cluster_labels,legend_title,region_labels`.

| Column | Meaning |
|---|---|
| `heatmap` | Id; used in every output name |
| `regions` | a csaw peak set name, a class name (`<id>.retained`, `<id>.lost`, `<id>.unclassified`), an intersect name or a BED path. A quoted list gives several region groups, one per entry. `enrichedheatmap` and rows with `cluster_k` take one region only |
| `samples` | Quoted list of `sample_id`s, one column each |
| `labels` | Quoted list of display labels, one per sample; a literal `\n` makes a line break |
| `window_bp`, `bin_bp` | bp on each side of the region centre, and the bin size; `window_bp` must be a multiple of `bin_bp` |
| `sort_using` | empty, or a `sample_id`: rows are sorted by that sample's signal |
| `cluster_k` | empty or `0`: no clustering; an integer >= 2: that many clusters; a range such as `2-6`: the gap statistic chooses k |
| `z_min`, `z_max`, `y_max` | empty, one number, or one per sample. `z` sets the colour scale; `y_max` the profile height. `z_max` and `y_max` may be `auto` (deepTools engine only) |
| `colors` | Quoted list, one hue per sample (`"#8c96c6"`); the colour ramp runs from white to the hue |
| `engine` | `deeptools` (computeMatrix + plotHeatmap) or `enrichedheatmap` (R EnrichedHeatmap) |
| `cluster_samples` | optional quoted subset of `samples` used for clustering; empty = all |
| `cluster_labels` | optional quoted list of k labels; if the count is not k, the labels are `cluster 1` .. `cluster k` |
| `legend_title` | optional EnrichedHeatmap legend title per sample; empty = the sample's antibody |
| `region_labels` | optional quoted list, one label per region of `regions` |

`auto` limits (deepTools engine): for each sample, the largest value of its mean profiles
(over the region groups). The limit of an antibody is 1.05 times the largest of its
samples' values, rounded up to two significant figures, so the samples of one antibody
share one scale. The table is
`results\heatmaps\limits\<heatmap>_limits.tsv`. Clustering is seeded k-means
(`heatmap_seed`); with a k range the chosen k is found by the gap statistic. The k table is
`cluster_k.tsv`, the cluster BEDs are in `results\heatmaps\cluster\`.

### tracks.csv (`run_tracks: true`)

`track_set,locus_or_gene,flank,samples,labels,colors,max_value,beds,scalebar`, then the
optional `log2ratio,bed_colors,bed_labels,ini_args`.

| Column | Meaning |
|---|---|
| `track_set` | Id; one figure per locus, named `<track_set>.<locus>.svg` |
| `locus_or_gene` | Quoted list of gene names (exact `gene_name` in the GTF) and/or regions `chr:start-end`; an unknown gene stops the run |
| `flank` | bp added on each side of a gene; one value or one per locus; ignored for regions |
| `samples`, `labels`, `colors` | As in heatmaps |
| `max_value` | One entry per sample: a number (fixed), `auto` (autoscale), `group:NAME` (the shared limit of the tracks with that name at this locus), or blank = `tracks_max_value_default` (`shared`: one group per antibody) |
| `beds` | Quoted list of csaw peak set names, class names, intersect names or BED paths, drawn as BED tracks |
| `scalebar` | bp, or `auto` (pyGenomeTracks picks the size); empty = no scale bar |
| `log2ratio` | optional `true` or `false`: adds a log2-ratio track after each sample's bigWig (needs `run_log2ratio: true`) |
| `bed_colors`, `bed_labels` | optional, one entry per item of `beds`; empty = the track colour and the file name |
| `ini_args` | optional; extra options for this track set's make_tracks_ini.py, written as in the command line and starting with `--`, e.g. `--genes-height=2 --merge-transcripts false --gene-rows 2 --arrow-interval 6 --bed-after SAMPLE_ID --gap-after SAMPLE_ID=0.5`. Other options: `--genes-style UCSC\|flybase\|tssarrow`, `--arrowhead-included true\|false`; `--gap-after` takes comma-separated `ID=height` entries (spacer height after that sample, default 0.25); `--bed-after` draws the BED tracks after that sample instead of after the last. No spaces or quotes inside a value |

Shared limits: a `shared` or `group:NAME` track takes its maximum from its group at its
locus (1.05 times the largest per-bin mean over the group's bigWigs, rounded up to two
significant figures). `tracks_limit_stat: max_bin` uses the largest value of the bigWig's own bins
instead. Resolved limits are written to `results\tracks\limits\`. The default makes the
conditions and replicates of one mark share one y-range; IgG forms its own group.

### intersects.csv (`run_intersect: true`)

`name,bed_a,bed_b,label_a,label_b,min_overlap,reciprocal[,output]`

- `bed_a`, `bed_b`: a peak set name, a class name or a BED path. Records of `bed_a`
  that overlap `bed_b` are written once each.
- `min_overlap`: fraction 0-1 of `bed_a` that must overlap; empty = 0.05.
- `reciprocal`: `true` also requires the same fraction of `bed_b`.
- `output`: `records` (default, also when empty) writes each `bed_a` record that overlaps `bed_b` once; `overlap` writes the overlapping pieces, one line per `bed_a`/`bed_b` pair (`bedtools intersect` without `-wa -u`), as the legacy step_7b script and the GG_CBX2 paper heatmaps used. Counts in the summary stay record-based; `n_pieces` counts the lines of `<name>.intersect.bed`.

Outputs: `results\intersect\<name>.intersect.bed`, `<name>.only_a.bed`, `<name>.intersect_summary.tsv`.

## Config keys (`config.yaml`)

Every key has a default in `nextflow.config` (in this folder, with a comment per key).
A project lists only what it changes. The reference paths have no default. Windows
paths are fine (`D:\REFERENCES\...`); relative paths are relative to the project folder.

**Inputs, outputs and reference**

| Key | Default | Meaning |
|---|---|---|
| `project_name` | `my_cutandrun` | prefix of the QC report and project-level outputs |
| `samplesheet`, `norm_runs`, `peaksets` | `samplesheet.csv`, `norm_runs.csv`, `peaksets.csv` | table files |
| `classifications`, `heatmaps`, `tracks`, `intersects` | `null` | table files used when their `run_*` switch is on |
| `fastq_dir` | `null` | not used by this version: name the FASTQs in the samplesheet |
| `outdir`, `metadata_dir`, `intermediate_dir` | `./results`, `./metadata`, `./intermediate_files` | output folders; **must be relative** to the project (a Windows drive path is refused) |
| `publish_intermediates` | `false` | Docker versions only. `true` publishes the trimmed FASTQs, the unfiltered aligned BAMs, the filtered BAMs, the SEACR bedGraphs and the BIN_COUNTS tables to `intermediate_dir`; `false` keeps them only in the work folder (kept while `keep_work: true`), so `intermediate_files\` stays empty |
| `publish_mode` | `auto` | `auto` = hard links from the work folder when the disk allows, else copies (with `work_volume: true`: copies, since hard links cannot cross into the volume); `copy` also accepted; `symlink`, `rellink` and `link` need `work_volume: false` (`symlink` and `rellink` keep the work folder) |
| `keep_work` | `true` | keep the work folder and resume cache after a successful run (the Docker volume, or `metadata\work` and `metadata\.nextflow` with `work_volume: false`); `false` deletes them after a successful run |
| `work_volume` | `true` | Docker versions only. `true` = Nextflow's work folder and resume cache in a per-project Docker volume on Docker Desktop's disk (fast; not visible in File Explorer); `false` = `metadata\work` and `metadata\.nextflow` in the project folder, as before. Switching loses the resume (see "Re-running") |
| `experiment` | `cutandrun` | `cutandrun` (TruSeq3-PE-2 adapters) or `cutandtag` (NexteraPE-PE) |
| `genome_index` | `null` | bowtie2 index **prefix** of the combined target + spike-in genome |
| `chrom_sizes`, `blacklist`, `gtf`, `tss_bed` | `null` | chromosome sizes of the combined genome; blacklist BED (pair-aware removal); gene GTF; TSS BED (TSS enrichment; fallback only, required when `gtf` is unset) |
| `spike_prefix` | `ecoli_chr` | contig-name prefix of the spike-in genome |
| `mito_chrom` | `chrM` | mitochondrial contig, dropped from the target BAM |
| `effective_genome_size` | `3100000000` | MACS3 `-g` and bamCoverage (T2T-CHM13) |
| `paired_end` | `true` | must be true (paired-end only) |
| `skip_groups` | `''` | comma list of `group` values kept in the samplesheet but skipped by all core steps |
| `max_cpus`, `max_memory_gb` | `null` | this computer's budget for the run; see "Resources" |

**Trimming, alignment, duplicates, normalisation**

| Key | Default | Meaning |
|---|---|---|
| `trim`, `trim_minlen` | `true`, `25` | trimmomatic before alignment; its MINLEN |
| `bowtie2_mode` | `local` | `local` (`--local --very-sensitive-local`) or `end_to_end` (`--end-to-end --very-sensitive`) |
| `bowtie2_minins`, `bowtie2_maxins` | `null`, `null` | bowtie2 `-I` and `-X`; `null` = not passed (bowtie2 defaults 0 and 500, as in the manuscript) |
| `mapq` | `30` | minimum MAPQ kept |
| `remove_duplicates` | `true` | remove duplicates from target BAMs (IgG always) |
| `spikein_dedup` | `true` | remove duplicates from the spike-in BAM |
| `normalisation` | `chipseqspikeinfree` | `chipseqspikeinfree`, `spikein`, `tmm` or `cpm`; informational except `cpm` (no scale factors, multiplier 1): the samplesheet `norm_run` decides per sample |
| `tmm_bin_size` | `10000` | bin size for TMM, for the normalisation MA plots and (when `repcor_bin_size` is null) the replicate metrics; csaw counts every library once per bin size |
| `norm_methods` | `'chipseqspikeinfree,tmm,cpm'` | methods computed for every `norm_run` besides its primary (`method` of `norm_runs.csv`); `''` = primary only |
| `norm_spikein` | `false` | also compute physical spike-in factors as a diagnostic (a run with pre-aligned BAMs: `unavailable`, with a warning) |
| `csif_bin_size`, `csif_cutoff`, `csif_max_turns` | `10000`, `1.2`, `0.99` | ChIPseqSpikeInFree settings when it is a diagnostic and the `norm_runs.csv` cells are blank |
| `norm_qc_min_count`, `norm_qc_prior`, `norm_qc_panels_per_page` | `20`, `0.5`, `12` | MA plots: minimum fragments summed over the two libraries of a pair; prior count of the log2 CPM (the reference lines do not depend on it); panels per SVG page |
| `run_replicate_correlation`, `repcor_bin_size`, `repcor_prior`, `repcor_signal_frac` | `true`, `null`, `2`, `0.10` | replicate concordance of libraries of the same group (`qc/replicate_correlation/`); bin size (`null` = `tmm_bin_size`); prior count of the Pearson log2 CPM; fraction of bins with the highest mean signal for the signal-bin metrics |
| `run_deeptools_correlation` | `true` | keep the deepTools all-sample Spearman heatmap (`qc/correlation/`) |
| `bigwig_bin_size`, `bigwig_extend_reads`, `bigwig_units` | `10`, `true`, `rpkm` | bamCoverage `--binSize`, `--extendReads`, units (`rpkm`, `cpm`, `none`) |
| `run_log2ratio` | `false` | also write target-vs-control log2-ratio bigWigs |

**Peaks**

| Key | Default | Meaning |
|---|---|---|
| `call_peaks` | `true` | master switch: `false` = no peak calling (csaw, MACS3, SEACR), FRiP, peak QC or peak annotation; bigWigs, normalisation, QC, heatmaps and tracks still run |
| `run_csaw` | `true` | csaw peak calling per `peaksets.csv` row |
| `csaw_window`, `csaw_spacing`, `csaw_bin` | `150`, `50`, `10000` | window width and spacing (bp); background bin size (bp) |
| `csaw_norm_bins` | `background` | TMM on background bins, or `all` |
| `csaw_small_filt`, `csaw_large_filt`, `csaw_merge_bp`, `csaw_max_width`, `csaw_fdr`, `csaw_lfc` | `3`, `1.5`, `50`, `0`, `0.05`, `1` | global fallbacks for blank `peaksets.csv` cells (`csaw_max_width` 0 = no limit) |
| `csaw_defaults` | `K27me3: {small_filt: 3, large_filt: 1.5, merge_bp: 50}`, `K27me2: {small_filt: 2, large_filt: 1, merge_bp: 50}` | per-antibody fallbacks for blank cells; a value in `config.yaml` replaces the whole map |
| `csaw_large_window`, `csaw_large_spacing` | `2000`, `500` | width and spacing (bp) of the wide window scale |
| `csaw_bg_max_enrich` | `1` | background bins: average abundance at most the median plus this (log2) |
| `csaw_count_filter` | `10` | minimum summed count over all libraries for a window |
| `csaw_max_frag` | `null` | csaw `max.frag` (also the counting of the normalisation MA plots, TMM and the replicate metrics); `null` = `bowtie2_maxins` (500 if that is null too) |
| `csaw_lfc_from` | `combined` | `combined` or `best`: which window's fold change the `lfc` rule uses |
| `csaw_disjoint_regions` | `false` | `true` makes merged regions disjoint (windows assigned by midpoint) |
| `run_csaw_cutoff_scan`, `csaw_scan_small_grid`, `csaw_scan_large_grid` | `false`, `'1,1.5,2,3'`, `'0.5,1,1.5,2'` | optional post-hoc cutoff diagnostic; never feeds the thresholds |
| `run_classify`, `classify_lfc_retained`, `classify_lfc_lost`, `classify_fdr` | `false`, `-1`, `-2`, `0.05` | region classification and the defaults for blank `classifications.csv` cells |
| `broad_marks` | `H3K27me3,H3K27me2,H3K36me3,H3K4me1,H3K79me2,H3K79me3,H3K9me1,H3K9me2,H3K9me3,H4K20me1,H3F3A` | antibodies called as broad marks (case-insensitive, with or without the histone prefix: `K27me3` = `H3K27me3`). Sources: the ENCODE histone ChIP-seq broad list (https://www.encodeproject.org/chip-seq/histone/), plus H3K27me2 by lab request and H3K9me3 (an ENCODE "exception", treated as broad by St. Jude). Any other antibody is narrow |
| `run_macs3`, `macs3_broad`, `macs3_qvalue`, `macs3_broad_cutoff` | `null`, `auto`, `0.01`, `0.1` | MACS3 peaks against the samplesheet control. `null` (not set) = on when every target has a `control_id`, otherwise skipped with a warning; an explicit `true` with a missing control stops the run, `false` turns it off. `macs3_broad`: `auto` (broad for `broad_marks`, narrow otherwise), `true` or `false`; broad uses `--broad --broad-cutoff` `macs3_broad_cutoff`, narrow uses `-q` `macs3_qvalue` |
| `run_seacr`, `seacr_mode`, `seacr_min_replicates`, `seacr_reciprocal_overlap` | `null`, `auto`, `2`, `0.5` | SEACR peaks and per-group consensus; skipped like MACS3 without controls. `seacr_mode`: `auto` (`relaxed` for `broad_marks`, `stringent` otherwise), `relaxed` or `stringent`. Lab policy: SEACR has no broad mode; relaxed = threshold between the knee and the peak of the total-signal curve, stringent = the peak. The mode of every group is written to `metadata\peak_call_modes.tsv` |

**Reports, heatmaps, tracks, annotation, intersections**

| Key | Default | Meaning |
|---|---|---|
| `run_report` | `true` | build the HTML QC report |
| `run_heatmaps`, `heatmap_seed` | `false`, `42` | heatmaps of `heatmaps.csv`; seed of the k-means clustering |
| `heatmap_reference_point`, `heatmap_missing_as_zero`, `heatmap_average_type_bins`, `heatmap_plot_type`, `heatmap_average_type_summary`, `heatmap_height`, `heatmap_width`, `heatmap_x_label`, `heatmap_ref_label`, `heatmap_y_label`, `heatmap_legend_location` | `center`, `true`, `mean`, `se`, `mean`, `12`, `3.5`, `distance (bp)`, `peak center`, `null`, `none` | deepTools heatmap styling (`heatmap_y_label` null = normalised units) |
| `run_tracks`, `tracks_width`, `tracks_dpi` | `false`, `15`, `600` | pyGenomeTracks figures of `tracks.csv`; width (cm) and resolution (dpi) |
| `tracks_max_value_default` | `shared` | max value of a blank `tracks.csv` entry: `shared`, `auto` or a number |
| `tracks_limit_stat` | `mean_bins` | statistic of the shared limits: `mean_bins` or `max_bin` |
| `track_number_of_bins`, `track_height`, `track_fontsize`, `track_bed_height`, `track_bed_color`, `track_scalebar_fontsize`, `track_genes_height`, `track_log2_min`, `track_log2_max` | `700`, `1.0`, `14`, `0.1`, `'#8c96c6'`, `16`, `3`, `-6`, `2` | pyGenomeTracks styling; the log2 limits are fixed. `track_fontsize` (accepted; bigWig track titles use pyGenomeTracks' own size, so it no longer changes them) |
| `run_annotation`, `annotation_colors`, `annotation_kegg`, `annotation_tss_region`, `annotation_extra_beds` | `false`, `''`, `true`, `3000`, `''` | ChIPseeker annotation of the peak sets; label colours (`label=#RRGGBB,...`); KEGG dotplot; bp around the TSS; extra regions annotated next to the peak sets |
| `run_intersect` | `false` | BED intersections of `intersects.csv` |

The optional launcher key `wsl_resize: false` (not in the template) turns off the
Docker memory change described under "Resources".

## Customizing one project's analysis

Every step runs one script from `bin\`. To change a step for one project only, copy
its script from this folder's `bin\` into the project's `utils\`, edit it, and add a
`config.yaml` key named after the script (file name without extension):

```yaml
norm_tmm: utils/norm_tmm.R           # replace a step's script for this project only
callpeaks_csaw_args: "..."           # or pass extra arguments to the shipped script
```

Each script's header documents its arguments and outputs; a replacement must keep that
contract. Changing a script key re-runs that step and the steps after it. Keys are the
base names of the scripts in `bin\`, so the base names are unique (`my_theme.R` for R
plots, `my_theme_py.py` for Python). Analyses that the pipeline does not do at all are
better written as your own script in `utils\`, run with `run_cutandrun script` (see
"Analyze a dataset" above).

## Inputs are never copied

- FASTQs, BAMs and reference files stay where they are. The pipeline reads them in
  place, never moves, copies, renames or deletes them, and `config.yaml` keeps pointing
  to them. Keep the path stable (if you move the files, update `config.yaml` and the CSVs).
- Write **Windows paths** in full (`D:\FASTQ\run7\x_R1.fastq.gz`) or relative to the
  project folder (`original_data\fastq\x_R1.fastq.gz`). Forward slashes work too.
- A drive named in `config.yaml` or in the CSVs is mounted **read-only** in the
  analysis (as `/mnt/<letter>`); a drive that is not connected stops the run with a message.
- **Network paths** (`\\server\share\...`) are refused: Docker Desktop cannot read them.
  Map the share to a drive letter (File Explorer: This PC > Map network drive), or copy
  the files to a local or USB drive, and write the path with that letter.
- The launcher records the paths as you wrote them and as the analysis used them in
  `metadata\host_paths.tsv`, and logs the mapping at the start of the run.

## Resources: how much of the computer an analysis uses

Each run has one **budget** of memory and CPU cores, and never goes over it. Docker
enforces it as a hard limit, so a run cannot freeze Windows by running out of memory.

**Memory.** The budget is the smallest of:
- 75% of the computer's memory, but at least 3 GB left for Windows;
- Docker Desktop's memory minus 0.5 GB (the Docker VM's own share);
- 14 GB, the most this pipeline can use (`max_useful_memory_gb` in `docker\requirements.conf`).

**CPU cores.** All cores but one (all but two from 8 cores up), so the desktop stays
responsive.

**Per-step memory requests** (re-measured): BAM_FILTER 7.5 GB, and BAM_PROCESS (MARKDUP,
DEDUP_SPLIT, BLACKLIST_FILTER and BAM_FINALIZE run in one task) 6 GB, each capped to the budget.
The side QC steps (PICARD_ELC, PRESEQ,
LIBCOMPLEXITY_FRAGS, FINGERPRINT, TSS_ENRICHMENT) run one at a time, so the BAM chain is not
starved.

`config.yaml` can lower the budget: `max_memory_gb` (decimals allowed, e.g. `2.5`) and
`max_cpus`. Neither can go above what Docker Desktop allows (the launcher warns).
Results do not change, only speed. Each run prints its budget at the start.

**Docker's memory limit (WSL 2 memory, Windows only).** With WSL 2, Docker Desktop
may use only a part of the computer's memory, set by the `memory=` line of
`%USERPROFILE%\.wslconfig`. For a chromatin run the launcher sets that limit to suit the
computer, by its total memory:

| Computer memory | Docker limit for the run |
|---|---|
| less than 12 GB | 6 GB |
| 12 to less than 24 GB | 8 GB |
| 24 GB or more | 16 GB |

- It **asks first** (`--yes` answers yes). It then edits only the `memory=` line, restarts
  Docker Desktop (about a minute; `wsl --shutdown` stops every WSL distribution, not only
  Docker), and runs the analysis.
- It is done **only when no other container is running**. Otherwise it keeps the current
  limit and says so.
- When the analysis finishes (or fails, or its window is closed), the previous value is put
  back the same way, if nothing else is running. The file `%USERPROFILE%\.wslconfig.run_cutandrun.bak`
  exists while the changed limit is in use. If the window was closed before the restore,
  the next `run_cutandrun` command that finds Docker running and nothing else running puts
  the old value back.
- `--no-wsl-resize` (on `run`) or `wsl_resize: false` in `config.yaml` leaves the limit alone.
  The budget is then computed from the limit Docker already has.

**The steps.** Each step asks for a share of the budget, set per step in
`conf\docker_tiers.config` for the computer's tier (small, mid or large budget). Steps run
side by side as far as the budget allows; the rest wait their turn. A step that is killed
for lack of memory is retried with more memory, never above the budget, up to three
times. The heavy steps (csaw, duplicate marking, ChIPseqSpikeInFree) need a few GB each,
which is why a small budget runs fewer steps at once. The run refuses to start, with an
explanation, if the budget is below the pipeline's minimum (`min_memory_gb` in
`docker\requirements.conf`). That figure, like `max_useful_memory_gb`, is a placeholder
until a pilot run on a real computer confirms it.

**Examples** (the rules above, for a computer with the Docker limit set as in the table):

| Computer | Docker limit | Memory budget | CPU cores |
|---|---|---|---|
| 8 GB, 4 cores | 6 GB | 5 GB | 3 |
| 16 GB, 8 cores | 8 GB | 7.5 GB | 6 |
| 32 GB, 12 cores | 16 GB | 14 GB | 10 |

**Keep memory free.** The pipeline caps its own memory, but other programs share it. Close
browsers with many tabs, large IDE sessions and other analyses before a long run. While a
run is going, the launcher watches Windows' free memory and makes Docker drop its file
cache when it gets low.

## Disk space

- **Plan about 5.6 times the total size of your FASTQ files** in all: about 4.6 times for the work
  files, on the drive that holds Docker Desktop's disk image, and about 1 time for the results,
  on the project's drive (with the default `work_volume: true`). With `work_volume: false` the work
  files are in the project folder, so the project's drive needs about 5.6 times. These are
  estimates from one 16-library project (27.5 GB of FASTQ): the work folder was 210 GB before the
  BAM steps were fused and is expected to be about 127 GB now (about 4.6 times the input), and the
  results about 26 GB (about 1 time). They are pending the next run. The launcher computes the same
  numbers before every run: it adds up the FASTQs named in the samplesheet and prints how much room
  the analysis may need on the project's drive (results) and on Docker's disk (work files). It warns
  when the free space on Docker's disk is below the work-file figure (the run stops only if that disk
  actually fills up, and can then be resumed). Results need room on the project's drive as well; the
  run stops before it starts only when less than 20 GB is free there.
- **To estimate:** in File Explorer, select the FASTQ folder of your samplesheet and read its
  size. A full dataset can need 100 GB or more, so check before you start.
- What uses the space: the work volume on Docker's disk (the working files of every step, kept
  so a re-run resumes; the markdup BAM and the name-sorted BAMs are deleted during the run as
  soon as the next step has used them), the trimmed FASTQs, unfiltered aligned BAMs, filtered
  BAMs, SEACR bedGraphs and BIN_COUNTS tables in `intermediate_files\` (only with
  `publish_intermediates: true`; the default on Windows is `false`, so the folder stays empty),
  and the results in `results\`. The raw FASTQs are not copied.
  With the work volume, files in `results\` are copies. With `work_volume: false` they are
  hard links to the work folder when the disk allows (the run log says so), so big files are
  not stored twice. Never edit a published file in place; replacing or deleting it is fine.
- To free space: `run_cutandrun clean` (results are kept), or `keep_work: false` (deletes
  the work files after a successful run).
- **The project can be on a hard disk.** With the work volume, only the results (and `intermediate_files\`
  with `publish_intermediates: true`) are written to the project folder, so a hard disk or a USB drive
  slows the steps much less than before (the launcher warns when the project is on a hard
  disk). Docker's disk image is what matters most for speed:
  keep it on an internal SSD, or at least on a drive with room (see "Good to know").

## Peak QC tracks

After the peak callers run, the step `PEAK_QC_REGIONS` (`bin/peak_qc_regions.py`) picks a short list of regions per mark from the peak sets, the region classes and the bigWigs, and draws one figure per region (bigWigs, peak sets and classes), so that the peak calls can be checked by eye. The regions are your own loci, the default Polycomb marker loci, a negative-control locus, and the largest, median and densest peaks per mark. It needs a peak caller and `gtf` in `config.yaml` (without `gtf` it is skipped with a warning); the figures also need `run_tracks`. The peak sets of every caller (csaw, MACS3, SEACR) are drawn as BED tracks, labelled `csaw <set>`, `MACS3 <group>`, `SEACR <group>` (colours black, #0072B2, #D55E00).

Keys in `nextflow.config` (defaults):

- `run_peak_qc = true`: master switch; silently skipped when there is no peak set.
- `peak_qc_loci = ''`: your own loci, `[NAME=]SPEC[:FLANK],...` (gene, glob, `chr:start-end`); empty = default Polycomb marker loci.
- `peak_qc_negative = 'auto'`: negative-control locus: `auto`, `none`, `GENE[:FLANK]` or `chr:start-end`.
- `peak_qc_max_regions = 15`: cap on the automatic regions; your own loci are never dropped.
- `peak_qc_exclude_bed = null`: optional extra BED masks for the peak-free high-signal search; the blacklist is always used.
- `peak_qc_nopeak_igg_max_fold = 1.5`: the peak-free high-signal tile excludes tiles whose IgG signal fold is above this; `null` or `''` = off; only applies when an IgG bigWig exists.
- `peak_qc_shape_marks = ''`: marks that get the largest, median and densest regions (at most 3); empty = every mark, the first 3 when there are more.

Outputs go to `results/tracks/peak_qc/`: the figures (track set `peak_qc`), `regions.tsv`, `regions.detail.tsv`, `tracks.loci.txt`, `params.json` and `peak_qc_regions.log`.

Caution: the signal-based criteria have not been tested on real bigWigs. Check the folds in `regions.detail.tsv` and the figures before trusting the regions.

Agents ask the user for their own loci of interest; they do not assume any.

## Software versions and default settings

The image `ghcr.io/genesig/omics` (tag in `docker\image.txt`) contains the tools below.
Versions are those pinned in the image's build (2026.10.2); the setup step prints the
main ones. R and Bioconductor are newer than in the manuscript (R 4.4.3 with
Bioconductor 3.20; the cluster's previous R image had R 4.3.3 and Bioconductor 3.18), so
numbers can differ slightly from results made with the older containers.

| Tool | Version |
|---|---|
| Nextflow (inside the image) | 26.04.6 |
| R (most steps; `Rscript`) | 4.4.3, Bioconductor 3.20 |
| csaw / edgeR / ChIPseeker / EnrichedHeatmap | 1.40.0 / 4.4.0 / 1.42.0 / 1.36.0 |
| ChIPseqSpikeInFree | 1.2.4 (commit f9cc88f; pinned, do not upgrade) |
| bowtie2 / samtools | 2.5.5 / 1.23.1 |
| trimmomatic / Picard | 0.40 / 2.27.5 (MarkDuplicates) |
| bedtools / FastQC / MultiQC | 2.31.1 / 0.12.1 / 1.35 |
| preseq / deepTools / pyGenomeTracks / MACS3 | 2.0.2 / 3.5.6 / 3.9 / 3.0.4 |
| SEACR | 1.3 (runs on R 4.2.3 in the image) |
| Python (`python3`) | from the image's `py` environment, with pyBigWig 0.3.25 |

Defaults that matter for the alignment: the manuscript's bowtie2 options are
`--local --very-sensitive-local --no-discordant --no-mixed --dovetail`. With
`bowtie2_mode: local` (the default) the pipeline passes those options (`--local
--very-sensitive-local`, `--no-discordant`, `--no-mixed`, `--dovetail`) and also `-q`,
which is carried over from the legacy scripts. No `-I` or `-X` is passed unless
`bowtie2_minins` or `bowtie2_maxins` is set, so bowtie2's defaults (0 and 500) apply, as
in the manuscript.

## Differences from the cluster

The steps, tools and results are those of the cluster version (`HPC`). Three things differ here:

- **The BAM steps run as one process.** MARKDUP, DEDUP_SPLIT, BLACKLIST_FILTER and BAM_FINALIZE
  run as one process, BAM_PROCESS, with one task per sample. The same step scripts run with the
  same arguments, the same results are written, and each step's command log keeps its name. The
  markdup BAM, the name-sorted target BAM and the blacklist-filtered BAM are deleted inside the
  task as soon as the next step has used them.
- **Intermediates are not published by default.** `publish_intermediates` (default `false`)
  decides whether the trimmed FASTQs, unfiltered aligned BAMs, filtered BAMs, SEACR bedGraphs and
  BIN_COUNTS tables are copied to `intermediate_files\`. The cluster always publishes them and has
  no such key.
- **One image holds all tools.** There are no per-step containers; every step runs in the same
  Docker image.

## Good to know

- **Where Docker keeps its data:** Docker Desktop stores images in a disk file on `C:` by
  default. The work volume lives on the same disk, so it needs room for about 4.6 times the
  FASTQ total (an estimate) and is best on an internal SSD. If `C:` is short of room or not an SSD, move it
  in Docker Desktop > Settings > Resources > Advanced > *Disk image location* (the launcher
  assumes no particular drive).
- **Reproducibility:** every tool runs at a pinned version inside the image, and the
  image tag is in `docker\image.txt`. Every step's `log\` folder holds the exact command
  that ran (`<sample_id>_<step>_command.sh`). `metadata\pipeline_info\` holds the Nextflow
  execution report, trace and timeline. A resumed run keeps the previous run's trace, report and timeline next to it, named with that run's timestamp (execution_trace.<timestamp>.txt); the start of each run's log names the image.
- **Updating the pipeline:** `git pull` in the `CHROMATIN_PROFILING` folder of the clone.
  If the update changed the image tag, the next `setup` or `run` downloads the new image.
- **Example projects** in `examples\projects\` (for the nucleation, CBX2 and EED datasets)
  point at the lab's legacy paths; copy one and edit its paths before you use it.

For how each step works (the tools, the normalisation assumptions, the csaw model, the
classification and the output layout), see `ANALYSIS.md` in this folder. The analysis is the
same in all three versions, except that the BAM steps run as one process here (see "Differences
from the cluster"); only how it is launched differs.
