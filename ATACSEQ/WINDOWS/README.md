# ATAC-seq pipeline for Windows (Docker)

The Richard Phillips Lab's ATAC-seq pipeline for paired-end data: adapter trimming, bowtie2 alignment in chunks,
MAPQ and flag filter, duplicate removal, blacklist removal, an optional Tn5 shift, fragment-size and TSS QC, TMM
scale factors on large bins (csaw and edgeR), bigWig tracks, MACS3 peaks per peak group, merged peak sets, FRiP and
a MultiQC report. It is the same pipeline as on the lab's cluster (same steps, same tools, same outputs; see
"Differences from the cluster"), packaged in one Docker image, so you don't install any bioinformatics software
yourself.

> **First time on this computer?** Follow `SETUP.md` in this folder once (Docker
> Desktop, WSL 2, where Docker's data lives, memory). Then come back here.
>
> The agent instructions for this pipeline (`AGENTS.md`, `CLAUDE.md` and the agent
> definitions in `.claude\agents\`) are copied into each new project by `init`. The
> steps below are complete on their own.

## What you need

- Windows 10 (22H2) or 11, 64-bit, with hardware virtualization enabled (on by
  default on most computers).
- **8 GB of memory or more, as Windows reports it** (16 GB or more is recommended). The launcher
  gives one analysis a budget of memory, and the pipeline needs at least 5 GB of that budget. On a
  computer with 8 GB the budget is about 5 GB, so the run works but slowly (the largest steps run
  one at a time). A computer that Windows reports with less than 8 GB is refused, with an
  explanation (see "Resources" below).
- Free disk space: on the drive that holds Docker Desktop's disk image (`C:` by default),
  **about 4.8 times the total size of your FASTQ files** for the work files (with the default
  `work_volume: true`; see "Disk space" below; an estimate, pending the first full ATAC-seq
  run), plus room for the Docker image (several GB, once); on the drive that holds your
  project, about 1 times the FASTQ total for the results, and at least 20 GB free. The raw
  FASTQs are not copied.
- Internet the first time (to download Docker and the pipeline's image).
- A reference set: a bowtie2 index of the genome (for the T2T-CHM13 index, about 5 GB on disk),
  the chromosome sizes, a blacklist BED, a gene annotation GTF, a TSS BED, and, if you want
  the contamination screen, a bowtie2 index of the contaminant genomes. The lab's PC keeps
  them in `D:\REFERENCES` (see "Where data and projects go" in `SETUP.md`); any other
  location works as long as you write its path in `config.yaml`.

## One-time setup (per computer)

1. Get the pipeline (once), as a sparse clone so that only this folder is downloaded:
   ```
   git clone --filter=blob:none --sparse <the lab's GitHub URL>
   cd PIPELINES
   git sparse-checkout set ATACSEQ/WINDOWS
   cd ATACSEQ\WINDOWS
   ```
   Then work in the `ATACSEQ\WINDOWS` folder; `git pull` later updates it.
   (The folder name of the clone is whatever the lab's repository is called.)
2. Open the `ATACSEQ\WINDOWS` folder and run, in PowerShell or Command Prompt:
   ```
   .\run_atacseq.bat setup
   ```
   - If Docker Desktop isn't installed, it offers to install it (Windows asks for
     administrator permission). Then **restart the computer**, open **Docker Desktop**
     from the Start menu once, accept its terms (no account needed: click *Skip*), and
     run `setup` again.
   - It then downloads the pipeline's Docker image (`ghcr.io/genesig/omics`, the tag
     is in `docker\image.txt`, currently `2026.10.3`): several GB, once. It takes a while;
     keep the computer awake and online.
   - It prints the versions of the main tools, and adds `run_atacseq` to your PATH,
     so new terminal windows can use it from any folder.

## Analyze a dataset

1. **Create a project folder** (one per dataset, *outside* the `ATACSEQ\WINDOWS`
   folder):
   ```
   run_atacseq init C:\Users\<you>\Documents\atac\my_dataset
   ```
   It creates the folder skeleton (`original_data\fastq\`, `metadata\`, `results\`,
   `utils\`), `config.yaml`, `examples\` (starter CSV templates), `RUN_ANALYSIS.bat`
   (double-click) and `utils\run_atacseq.bat` (the command used below). It also
   writes `.gitattributes` and, when git is installed, starts a git repository in the
   folder with the first commit. It doesn't run anything.
2. **Add your tables** to the project folder (the project root, next to `config.yaml`).
   Copy the templates you need from `examples\`, fill them in (Excel is fine: *Save As
   → CSV UTF-8*):
   - `samplesheet.csv` (required): one row per library. See "Input tables" below.
   - `peak_merges.csv` (optional): one row per merged peak set. Set `peak_merges` in
     `config.yaml` to point at it.
   - `comparisons.csv`, `peak_files.csv`, `diff_tables.csv` and `motif_targets.csv` (optional,
     only for the downstream analyses; see "Downstream analyses").
3. **Edit `config.yaml`** (Notepad is fine). Set `project_name`, check the table file
   names, and set the reference paths: `genome_index` (the bowtie2 index prefix, e.g.
   `D:\REFERENCES\bowtie2\t2t\<prefix>`), `blacklist`, `chrom_sizes`, `gtf`, `tss_bed`,
   `effective_genome_size` (MACS3 `-g`, required when a sample has a `peak_group`) and
   `contamination_index` (optional; `null` skips the contamination screen). Every other key
   has a default; list only the ones you change. Key summary below. The example
   `config.yaml` marks every line that must be filled in with `# set me`.
4. **Run it:** double-click **`RUN_ANALYSIS.bat`** in the project folder, or in a
   terminal in that folder: `.\utils\run_atacseq.bat`.
   - Before it starts, the launcher checks your inputs (paths, the samplesheet, the
     output folders, the free disk space on the project's drive and on Docker's disk) and
     sets this computer's memory for the analysis (see "Resources" below; it asks first).
   - The analysis runs inside Docker in the background; the window shows its
     progress. Closing the window does **not** stop it, but keep the computer on and
     awake (plugged in, sleep off) until it finishes.
   - Check on it any time: `.\utils\run_atacseq.bat status`.
5. **Results** are in `results\`, one folder per step (`results\qc\multiqc\multiqc_report.html`
   is the MultiQC summary; `results\qc\fragment_sizes\` has the fragment-size plots). The
   bigWig tracks are in `results\bigwig\`, the TMM scale factors in `results\scale_factors\`,
   the peaks in `results\peaks\macs3\` (and `results\peaks\merged\` with `peak_merges`), and
   the final BAMs in `results\bam\final\`. The layout is listed in `ANALYSIS.md`.

**If something fails**, `.\utils\run_atacseq.bat status` says which step failed and
shows the last lines of the log. The files of the failed step are in
`metadata\failed_tasks\<timestamp>\`. Fix the input and run again: the steps that
already finished are reused, not redone.

**Re-running a finished project** only redoes what changed (for example after you
change a bigWig setting or a peak threshold), because Nextflow's work folder and
its resume cache are kept by default (`keep_work: true`). With the default
`work_volume: true` they are in a Docker volume named after the project, on Docker
Desktop's own disk (not in the project folder; the volume is not visible in File
Explorer). With `work_volume: false` they are `metadata\work\` and `metadata\.nextflow\`
in the project folder. To free disk space, run
`.\utils\run_atacseq.bat clean` (asks first; results are kept) or set
`keep_work: false` in `config.yaml` (deletes the work files after a successful run).
The next run then recomputes everything. A failed or stopped run always keeps them.

**Re-running only some steps (for example the plots).** The work volume is never deleted by a run: with the default `keep_work: true` it survives successful runs, and only `run_atacseq clean` removes it. After changing a setting or a CSV, run again: Nextflow resumes, reuses every step whose inputs are unchanged and recomputes only the steps downstream of the change (a plot setting reruns just that plot). Do not run `clean` in a real project: it would force the whole analysis to be recomputed.

**Changing `work_volume`** starts the run from scratch, because the resume state is in the
other place: a project whose `metadata\work\` already has content is asked first (`--yes`
answers yes). To resume the earlier run, set `work_volume: false`. The old `metadata\work\`
is left alone (`clean` deletes it). Switching back and forth loses the resume each way.

**Your own analysis** (a plot or test the pipeline doesn't do): put the script in
`utils\` (`.R`, `.py` or `.sh`) and run it with
`.\utils\run_atacseq.bat script utils\my_analysis.R`. It runs in the pipeline's
software environment and its log goes to `metadata\script_logs\`. See "Customizing one
project's analysis" below.

## All commands

Run from the project folder (`.\utils\run_atacseq.bat ...`) or from anywhere after
setup (`run_atacseq ...`, acts on the current folder):

| Command | |
|---|---|
| `run_atacseq setup` | one-time: Docker Desktop (offers to install it) + the pipeline's Docker image |
| `run_atacseq init FOLDER` | create a project folder (default: the current folder) |
| `run_atacseq` (or `run`) | run / resume the analysis of the current folder (`--detach`: start and return, don't follow the log) |
| `run_atacseq run other_config.yaml` | ... with another config file in the folder |
| `run_atacseq status` | running? finished? which step failed and why? the work volume's name and size |
| `run_atacseq logs` | watch the live log (Ctrl+C stops watching, not the analysis) |
| `run_atacseq stop` | stop the analysis (running again resumes it) |
| `run_atacseq clean` | delete Nextflow's work folder and resume cache to free disk space: the project's Docker work volume, and `metadata\work\` and `metadata\.nextflow\` if present; asks first (`--yes` skips); the next run starts from scratch, results are kept |
| `run_atacseq script utils\NAME.R [args]` | run your own script (`.R`, `.py`, `.sh`; inside the project folder) in the pipeline's software environment, in the foreground; log in `metadata\script_logs\` |
| `run_atacseq help` | the list of commands |

Options: `--yes` answers yes to the launcher's questions, including the Docker memory
question (for scripts and AI agents); `--detach` (with `run`) starts the analysis and
returns without following the log; `--pause` waits for Enter before the window closes
(used by `RUN_ANALYSIS.bat`); `--no-wsl-resize` (with `run`) leaves Docker Desktop's
memory limit as it is.

## Input tables

All tables live in the project folder root (next to `config.yaml`). Paths in them
follow the rules in "Inputs are never copied" below. Lists inside one cell are
comma-separated and quoted (`"peak_1,peak_2"`).

### samplesheet.csv (one row per library)

Header (column order fixed):

```
sample_id,fastq_1,fastq_2,condition,replicate,tmm_group,tmm_group_legacy,peak_group,exclude,bam
```

| Column | Required | Meaning |
|---|---|---|
| `sample_id` | yes | Unique library name, used verbatim in file names. Matches `^[A-Za-z0-9][A-Za-z0-9_.-]*$`. |
| `fastq_1`, `fastq_2` | yes, unless `bam` is set | Paths of the R1 and R2 FASTQ files (paired-end). Windows paths (`D:\FASTQ\run7\x_R1.fastq.gz`) or relative to the project folder. Empty on a row with a `bam`. |
| `condition` | no | Label for the QC tables (fragment sizes, TSS, FRiP). |
| `replicate` | no | Label for the QC tables. |
| `tmm_group` | no | TMM group (exact values). Blank = not normalised, no bigWig, no TSS enrichment. Read when `tmm_group_column` is `tmm_group`. Needs at least 2 active samples per group. |
| `tmm_group_legacy` | no | Second TMM group column. Read when `tmm_group_column` is `tmm_group_legacy` (legacy file-name membership). |
| `peak_group` | no | MACS3 peak group (exact values). Replicates are pooled, no control. Blank = no peak calling and no FRiP. |
| `exclude` | no | `true`, `1` or `yes` drops the row from the whole run. Blank, `false`, `0` or `no` keeps it. Other values are a startup error. |
| `bam` | no | A final BAM of the library instead of FASTQ (BAM_PREALIGNED). Then `fastq_1` and `fastq_2` must be empty. The file must end in `.bam` and exist. The BAM is assumed final (filtered, duplicate-removed, blacklist-filtered). Each active row has FASTQ or a BAM, not both. Windows paths, as for the FASTQs. |

A missing optional column is blank for every row. Any other column name is a startup
error. The file may carry a UTF-8 byte order mark (Excel's *CSV UTF-8*). Quoted fields are
supported. Two FASTQ files may not share one file name, and `fastq_1` and `fastq_2` may not
be the same file. The `examples\samplesheet.csv` template has the columns up to `exclude`; a
`bam` column can be added to it, and a sheet without the column is valid.

Startup errors (before any analysis starts) name the problem: unknown or missing columns,
duplicate `sample_id`, a bad id, a missing FASTQ, a TMM group with one sample,
`tmm_group_column` not a column, every row excluded, an invalid value for a switch or a
number, a missing bowtie2 index, a missing `blacklist`, `chrom_sizes`, `gtf` or
`contamination_index` file, a missing `tss_bed` when TSS enrichment runs with TMM groups, a
missing `effective_genome_size` when MACS3 runs with peak groups, a `bam` that is not a `.bam`
file or is missing, a row with both FASTQ and a BAM, and the downstream checks (see
"Downstream analyses").

### peak_merges.csv (optional)

```
name,peak_groups
merged_peaks,"peak_1,peak_2"
```

One merged BED per row: `results\peaks\merged\<name>.merged.bed` (bedops merge of the listed
groups' MACS3 peak files). `name` is unique and matches the id pattern. `peak_groups` is a
quoted comma list; each group must have an active sample. Setting `peak_merges` to a file
needs `run_macs3: true`. Use `peak_merges: null` for none.

## Config keys (`config.yaml`)

Every key has a default in `nextflow.config` (in this folder, with a comment per key).
A project lists only the keys it changes. Windows paths are fine (`D:\REFERENCES\...`).
Relative paths are relative to the project folder.

**Inputs, outputs and reference**

| Key | Default | Meaning |
|---|---|---|
| `project_name` | `my_atacseq` | prefix of the project-level outputs |
| `samplesheet`, `peak_merges` | `samplesheet.csv`, `null` | table files; `peak_merges: null` = no merged peak sets |
| `fastq_dir` | `null` | not used by this version: name the FASTQs in the samplesheet |
| `outdir`, `metadata_dir`, `intermediate_dir` | `./results`, `./metadata`, `./intermediate_files` | output folders; **must be relative** to the project (a Windows drive path is refused) |
| `publish_intermediates` | `false` | Docker versions only. `true` publishes the trimmed FASTQs (`intermediate_files\trimmed\`), the aligned BAMs (`intermediate_files\aligned\`, joined from the chunks when chunked) and the filtered BAMs (`intermediate_files\filtered\`) to `intermediate_dir`; `false` keeps them only in the work folder (kept while `keep_work: true`), so `intermediate_files\` stays empty |
| `publish_mode` | `auto` | `auto` = hard links from the work folder when the disk allows, else copies (with `work_volume: true`: copies, since hard links cannot cross into the volume); `copy` also accepted; `link`, `symlink` and `rellink` need `work_volume: false` (with `work_volume: true` they are refused, since they would point into the volume) |
| `keep_work` | `true` | keep the work folder and resume cache after a successful run (the Docker volume, or `metadata\work` and `metadata\.nextflow` with `work_volume: false`); `false` deletes them after a successful run |
| `work_volume` | `true` | Docker versions only. `true` = Nextflow's work folder and resume cache in a per-project Docker volume on Docker Desktop's disk (fast; not visible in File Explorer); `false` = `metadata\work` and `metadata\.nextflow` in the project folder. Switching loses the resume (see "Re-running") |
| `genome_index` | `null` | bowtie2 index **prefix** of the genome. Required |
| `chrom_sizes` | `null` | chrom.sizes of the genome. Recorded and checked if set; no step uses it |
| `blacklist` | `null` | blacklist BED for pair-aware removal (both mates dropped when either overlaps). `null` or empty = skip BLACKLIST_FILTER |
| `gtf` | `null` | gene annotation GTF. Recorded and checked if set; no step uses it |
| `tss_bed` | `null` | TSS BED for TSS enrichment. Required when TSS enrichment runs and a sample has a `tmm_group` |
| `effective_genome_size` | `null` | MACS3 `-g` (and bamCoverage `--effectiveGenomeSize` when set). Required when MACS3 runs and a sample has a `peak_group` |
| `contamination_index` | `null` | bowtie2 index **prefix** of contaminant genomes. `null` or empty = skip CONTAMINATION |
| `max_cpus`, `max_memory_gb` | `null` | this computer's budget for the run; see "Resources" |

**Steps and switches**

| Key | Default | Meaning |
|---|---|---|
| `run_fastqc` | `true` | FastQC on raw reads and, with `trim`, on trimmed reads |
| `trim` | `true` | Trimmomatic before alignment. `false` aligns the raw FASTQs |
| `run_tss_enrichment` | `true` | TSS enrichment on the bigWigs (needs `tss_bed`) |
| `run_macs3` | `true` | MACS3 per peak group, with peak statistics and FRiP |
| `run_multiqc` | `true` | MultiQC report |
| `tn5_shift` | `false` | `true` applies deepTools `alignmentSieve --ATACshift`; the shifted BAM feeds every later step |
| `remove_duplicates` | `true` | `sambamba markdup -r`. `false` marks duplicates and keeps them |
| `help` | `false` | print the usage and exit |

**Trimming and alignment**

| Key | Default | Meaning |
|---|---|---|
| `trim_adapter_set` | `NexteraPE-PE` | Trimmomatic adapter file name (without `.fa`), located inside the image |
| `trim_illuminaclip` | `'2:30:10:8:True'` | ILLUMINACLIP parameters after the adapter file |
| `trim_steps` | `'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20'` | further Trimmomatic steps, space separated; may be empty |
| `trim_minlen` | `25` | Trimmomatic `MINLEN` |
| `bowtie2_mode` | `local` | `local` (`--local --very-sensitive-local`) or `end_to_end` (`--end-to-end --very-sensitive`) |
| `bowtie2_minins` | `10` | bowtie2 `-I`. `null` = not passed (legacy) |
| `bowtie2_maxins` | `2000` | bowtie2 `-X` |
| `bowtie2_args` | `''` | further bowtie2 arguments, after the fixed flags. Legacy `'-N 1 -t'` |
| `align_chunk_pairs` | `10000000` | read pairs per ALIGN task (chunked alignment, joined by ALIGN_MERGE). Integer; `0` = no chunking (one ALIGN task per library); negative or non-integer is a startup error. See "Alignment chunks" below |
| `alignment_min_pct` | `20` | overall alignment rate below which ALIGNMENT_SUMMARY writes a warning (never fails) |
| `mapq` | `30` | minimum MAPQ kept (`samtools view -q`). `0` reproduces the legacy behaviour |
| `remove_chroms` | `chrM` | comma list of contigs removed from the BAM. `''` keeps all contigs (legacy) |

**QC, TMM, bigWig and peaks**

| Key | Default | Meaning |
|---|---|---|
| `tss_flank` | `2000` | TSS flank, bp on each side of the TSS |
| `tmm_bin_size` | `10000` | TMM bin width, bp (at least 1000) |
| `tmm_max_frag` | `auto` | maximum fragment size for TMM: `auto` (largest proper-pair fragment of the group) or a positive integer |
| `tmm_group_column` | `tmm_group` | samplesheet column that defines the TMM groups (`tmm_group` or `tmm_group_legacy`) |
| `bigwig_bin_size` | `10` | bamCoverage `--binSize` |
| `bigwig_units` | `none` | `none` (bamCoverage `--normalizeUsing None`, legacy), `cpm` or `rpkm`. Also sets how the multiplier is written |
| `bigwig_extend_reads` | `true` | bamCoverage `--extendReads`: paired-end fragment coverage (legacy `false`: per-mate read coverage) |
| `bigwig_center_reads` | `false` | bamCoverage `--centerReads` (legacy `true`, a no-op without `--extendReads`; with it, read-length windows at fragment midpoints) |
| `bigwig_args` | `''` | further bamCoverage arguments |
| `macs3_qvalue` | `0.05` | MACS3 `-q` |
| `macs3_min_length` | `250` | MACS3 `--min-length`. `null` = off |
| `macs3_max_gap` | `500` | MACS3 `--max-gap`. `null` = off |
| `macs3_fe_cutoff` | `2.0` | MACS3 `--fe-cutoff`. A value of 1 or less = off and no `.fe` in the file name |
| `macs3_call_summits` | `false` | MACS3 `--call-summits` (writes summit BEDs) |
| `macs3_broad` | `false` | MACS3 `--broad` (broadPeak and gappedPeak output) |
| `macs3_broad_cutoff` | `0.1` | MACS3 `--broad-cutoff`, used with `macs3_broad` |
| `macs3_args` | `''` | further MACS3 arguments, placed on the command line by the module |

**Downstream analyses** (all off by default; see "Downstream analyses")

| Key | Default | Meaning |
|---|---|---|
| `run_differential` | `false` | DIFF_ACCESS: csaw and edgeR differential accessibility, one task per comparison (needs `comparisons`) |
| `run_annotation` | `false` | DIFF_ANNOTATE: ChIPseeker annotation of each comparison's peaks table (needs `gtf`) |
| `run_motifs` | `false` | MOTIF_ENRICHMENT per comparison and `motif_regions` set, per `motif_targets` row; MOTIF_COMMON per region set |
| `comparisons` | `null` | comparisons CSV (`name,case,control,group_column,block_column,case_peaks,control_peaks,pooled,fixed_dispersion`). Required when a `run_*` switch is on |
| `peak_files` | `null` | CSV `name,path`: existing peak files usable as peak groups in `case_peaks` and `control_peaks` |
| `diff_tables` | `null` | CSV `name,regions,path`: existing differential tables (`peaks` or `promoters`), not recomputed |
| `motif_targets` | `null` | CSV `name,type,target,comparison,label,min_peaks`: target-vs-other motif analyses (`gene` or `regions`); needs `run_annotation` |
| `species` | `Homo sapiens` | `Homo sapiens` or `Mus musculus`; sets the default BSgenome and annotation package |
| `genome_build` | `hs1` | Genome build of the BAMs and the GTF (AM: `mm39`) |
| `bsgenome` | `null` | BSgenome package for the chromosome lengths; `null` = `BSgenome.Hsapiens.UCSC.<genome_build>` (`BSgenome.Mmusculus.UCSC.<build>` for mouse) |
| `org_db` | `null` | Annotation package for gene symbols; `null` = `org.Hs.eg.db` (`org.Mm.eg.db` for mouse) |
| `diff_regions` | `peaks` | Region sets computed by DIFF_ACCESS: `peaks` or `peaks,promoters` (promoters need `gtf`; the AV example uses `peaks,promoters`) |
| `diff_minq` | `10` | csaw `minq` (minimum mapping quality) |
| `diff_pe` | `both` | csaw `pe`: `both` or `none` (older AV and AM: `none`) |
| `diff_max_frag` | `500` | csaw `max.frag` (csaw default; the legacy never set it) |
| `diff_dedup` | `null` | csaw `dedup`: `null` = `false` when `remove_duplicates` is on, else `true`; or `true`, `false`, `auto` (the legacy rule: TRUE for comparisons with at most two BAMs) |
| `diff_bin_size` | `10000` | Window width (bp) of the normalisation factors |
| `diff_lfc` | `0.585` | absolute logFC for the up and down BEDs and the printed counts (older AV and AM: `1`) |
| `diff_fdr` | `0.1` | FDR for the printed counts only (no table is filtered by it) |
| `diff_table` | `best` | Which overlapResults table gives the statistics: `best` or `combined` (AM legacy: `combined`) |
| `diff_qlf_min_libs` | `3` | Quasi-likelihood F-test when a comparison has at least this many libraries; otherwise fixed dispersion and LRT |
| `diff_fixed_dispersion` | `0.05` | LRT dispersion when the F-test is not used (the legacy gd2 chunk used 0.15: set `fixed_dispersion` per comparison) |
| `diff_promoter_up` | `2000` | Promoter window upstream of the TSS (bp) |
| `diff_promoter_down` | `200` | Promoter window downstream of the TSS (bp) |
| `diff_promoter_min_fc` | `4` | Promoters kept when their background abundance is above log2 of this (filterWindowsGlobal) |
| `diff_cpus` | `4` | BiocParallel workers for the read counting (DIFF_ACCESS cpus) |
| `annotation_tss_region` | `-2000,200` | ChIPseeker `tssRegion`: upstream and downstream of the TSS (bp) |
| `annotation_fdr` | `null` | `null` = annotate all regions (AV). A number keeps regions with FDR below it (older AV: `0.1`) |
| `motif_regions` | `promoters` | Region sets analysed per comparison: `promoters`, `peaks`, `all_sites`, `promoter_sites`, `non-promoter_sites`, `distal_sites` (AM: `promoter_sites,non-promoter_sites,all_sites`) |
| `motif_lfc` | `0.585` | Bins at logFC above `lfc` (increased) and below `-lfc` (decreased) (older AV and AM: `1`) |
| `motif_min_score` | `10` | PWM hit threshold (`min.score`) |
| `motif_test` | `binomial` | Test of `calcBinnedMotifEnrR` |
| `motif_max_kmer` | `3` | `maxKmerSize`: composition correction |
| `motif_background` | `otherBins` | Background of `calcBinnedMotifEnrR` |
| `motif_seed` | `123` | `RNGseed` of the serial BiocParallel parameter |
| `motif_top` | `40` | Motifs shown in the heatmaps when more than this pass the selection |
| `motif_padj_primary` | `2.0` | -log10 adjusted P of the primary selection (legacy padj < 0.01) |
| `motif_padj_relaxed` | `1.0` | -log10 adjusted P of the relaxed selection (legacy padj < 0.1), used when two or fewer motifs pass the primary one |
| `motif_min_bin_peaks` | `10` | Diff mode: an analysis with a bin below this many regions is skipped (`NAME.skipped.txt`) |
| `motif_target_min_peaks` | `15` | Target mode: default minimum regions in the target and in the other set (a CSV `min_peaks` overrides it) |
| `motif_pwm_rds` | `null` | Cached PWMList (legacy `rds/pwms.rds`); `null` = JASPAR2024 from the package |
| `motif_gene_column` | `geneId` | Annotation column searched by a gene target |
| `motif_match` | `substring` | Gene target matching: `substring` (legacy `grepl`) or `exact` |
| `motif_cpus` | `1` | `1` = serial run with the legacy seed (the reproducible result); more uses MulticoreParam and is not bit-reproducible |
| `motif_common_padj` | `0.10` | MOTIF_COMMON: a motif is significant in a comparison at this adjusted P |
| `motif_common_min_lines` | `2` | MOTIF_COMMON: a motif is common when significant in at least this many comparisons |
| `motif_common_label_strip` | `null` | Regex removed from the comparison names to make the labels (AV example: `_high_vs_low.*`) |
| `plot_font` | `Helvetica` | Font family of the text in every downstream plot (svg keeps it editable) |
| `plot_format` | `svg` | `svg` or `pdf` (`pdf` = the legacy file format) |

Not in `nextflow.config`, read by the launcher only: `wsl_resize: false` (not in the
template) turns off the Docker memory change described under "Resources".

## Outputs

```
results\qc\fastqc\raw\, results\qc\fastqc\trimmed\   FastQC zip and html per FASTQ; log\ holds each command
results\qc\trimmomatic\                              Trimmomatic logs
results\qc\bowtie2\                                  bowtie2 alignment logs per sample (summed over chunks); log\ also holds the
                                                     COUNT_PAIRS, ALIGN (per chunk) and ALIGN_MERGE command files
results\qc\alignment\                                alignment_summary.tsv
results\qc\chrom_fraction\                           chrom_fraction.tsv (only when remove_chroms is set)
results\bam\stats\                                   S.raw.flagstat/idxstats (before the filter), S.final.flagstat/idxstats,
                                                     S.orphan_check.txt, S.sambamba_markdup.log
results\bam\final\                                   S.final.bam and .bai per library (also the BAM_PREALIGNED outputs of the bam rows)
results\bam\tn5\                                     S.tn5.bam and .bai (only when tn5_shift is true)
results\qc\contamination\                            S.contamination.tsv, S.contamination_summary.txt (only when contamination_index is set)
results\qc\fragment_sizes\                           fragment_sizes.tsv, fragment_size_distribution.tsv, fragment_sizes_density.svg,
                                                     fragment_sizes_boxplot.svg
results\scale_factors\                               <group>.scale_factors.tsv, <group>.ma_plot.png and .pdf, <group>.max_frag.txt
results\bigwig\                                      <sample>.<method>.<units>.bw, for example sample.tmm.none.bw
results\qc\tss_enrichment\                           S.tss_enrichment.tsv, S.tss_matrix.gz, S.tss_profile.svg, tss_enrichment.tsv
results\peaks\macs3\                                 <group>.qval<q>[.minlen<l>][.maxgap<g>][.fe<fe>][.broad]_peaks.narrowPeak (or broadPeak),
                                                     _peaks.xls, optional _summits.bed and _peaks.gappedPeak, .macs3.log
results\peaks\merged\                                <name>.merged.bed (only with peak_merges)
results\qc\peaks\                                    peak_stats.tsv
results\qc\frip\                                     S.macs3.frip.tsv, frip.tsv
results\qc\multiqc\                                  multiqc_report.html, multiqc_report_data\
results\annotation\txdb\                             txdb.sqlite (only when run_annotation or a promoter table is computed)
results\differential\, results\differential\bed\     differential tables, plots and BEDs per comparison (run_differential)
results\annotation\                                  peak annotation tables, BEDs, plots and summary (run_annotation)
results\motif\<name>\<region>\, results\motif\<name>\, results\motif\common\<region>\   motif analyses (run_motifs; see "Downstream analyses")
intermediate_files\trimmed\                          S.R1.trim.fastq.gz, S.R2.trim.fastq.gz (only with publish_intermediates: true)
intermediate_files\aligned\                          S.aligned.bam (unsorted, includes unmapped records; joined from the chunks; only with publish_intermediates: true)
intermediate_files\filtered\                         S.filtered.bam and .bai (only with publish_intermediates: true)
metadata\nf_logs\                                    nextflow_<timestamp>.log (the log of each run), last_run.status
metadata\pipeline_info\                              execution_report_<timestamp>.html, execution_trace_<timestamp>.txt, execution_timeline_<timestamp>.html, pipeline_dag_<timestamp>.html
metadata\host_paths.tsv                              Windows paths as written and as the analysis used them (only when a drive path was translated)
metadata\failed_tasks\<timestamp>\                   for a failed run: each failed step's command and its messages
```

Every step's `log\` subfolder holds `<name>_command.sh`, the exact command that ran. The
BIGWIG log files are named `<sample>_<method>_bigwig_command.sh`. The DEDUP BAM and the
BLACKLIST_FILTER BAM are not published (they are work files). With chunked alignment the
per-chunk ALIGN BAMs and logs are not published either; their commands are
`<sample>.chunkNNN_align_command.sh`.

## Alignment chunks

COUNT_PAIRS counts each library's read pairs, and ALIGN aligns `align_chunk_pairs` pairs per
task (bowtie2 `-s`/`-u`, no temporary FASTQs). ALIGN_MERGE joins the chunk BAMs in order and
sums their bowtie2 logs, so the merged BAM has the same records as one unchunked run
(`align_chunk_pairs: 0`).

On one computer the chunks share the same cores and memory budget as every other step, so
they run as the budget allows, side by side with each other and with the other libraries.
Chunking does not make alignment faster on its own. The default is 10 million pairs per
chunk; a library with fewer pairs than that is one chunk. The outputs do not depend on
`align_chunk_pairs`.

Downstream steps start per group: NORM_TMM and MACS3 start for a group as soon as its own
libraries are aligned, and FRiP and PEAK_MERGE wait only for their own peak groups.

## Customizing one project's analysis

Every step runs one script from `bin\`. To change a step for one project only, copy
its script from this folder's `bin\` into the project's `utils\`, edit it, and add a
`config.yaml` key named after the script (file name without extension):

```yaml
norm_tmm: utils/norm_tmm.R           # replace a step's script for this project only
fastqc_args: "--nogroup"             # or pass extra arguments to the shipped script
```

Each script's header documents its arguments and outputs; a replacement must keep that
contract. Changing a script key re-runs that step and the steps after it. Keys are the
base names of the scripts in `bin\`, so the base names are unique (`my_theme.R` for R
plots, `my_theme_py.py` for Python). A key that names no script, with a `.sh`, `.py` or `.R`
path as its value, is a startup error. Analyses that the pipeline does not do at all are
better written as your own script in `utils\`, run with `run_atacseq script` (see
"Analyze a dataset" above).

The scripts are:

| Script (`bin\` file = config key) | Purpose |
|---|---|
| `align.sh` | ALIGN: bowtie2 paired-end alignment of one chunk (`--chunk-skip`, `--chunk-upto`) piped into `samtools view -b`; never SAM on disk |
| `align_merge.sh` | ALIGN_MERGE: `samtools cat` of the chunk BAMs in chunk order; the chunk logs are summed by `merge_bowtie2_logs.py` |
| `alignment_summary.py` | ALIGNMENT_SUMMARY: alignment rates from the bowtie2 logs; warns below `alignment_min_pct` |
| `bam_filter.sh` | BAM_FILTER: sort, raw stats, MAPQ and flag filter, `remove_chroms` |
| `bam_finalize.sh` | BAM_FINALIZE: coordinate sort, index, final stats, orphan check (exit 3) |
| `bam_prealigned.sh` | BAM_PREALIGNED: a final BAM from the samplesheet `bam` column, copied (sorted only if needed), indexed, with flagstat and idxstats |
| `bigwig.sh` | BIGWIG: bamCoverage with the scale-factor multiplier |
| `blacklist_filter.sh` | BLACKLIST_FILTER: name sort and `bedtools pairtobed -type neither` |
| `callpeaks_macs3.sh` | MACS3: callpeak per peak group, replicates pooled, no control |
| `contamination.sh` | CONTAMINATION: unmapped pairs screened against the contaminant bowtie2 index |
| `count_pairs.sh` | COUNT_PAIRS: read pairs of a library (R1 and R2 must match); sets the number of ALIGN chunks |
| `dedup.sh` | DEDUP: sambamba markdup, with `-r` when `remove_duplicates` is true |
| `diff_accessibility.R` | DIFF_ACCESS: csaw counts and edgeR tests for one comparison (peaks and promoters); the legacy differential chunks |
| `diff_annotate.R` | DIFF_ANNOTATE: ChIPseeker annotation of one differential table, four sets; the legacy peak_annotation chunk |
| `fastqc.sh` | FASTQC_RAW and FASTQC_TRIM: `fastqc --extract` per FASTQ |
| `merge_bowtie2_logs.py` | helper of ALIGN_MERGE: sums the paired-end bowtie2 summaries of the chunk logs |
| `motif_common.R` | MOTIF_COMMON: motifs significant in several comparisons, per region set (find_common_motifs.Rmd) |
| `motif_enrichment.R` | MOTIF_ENRICHMENT: monaLisa binned motif enrichment, diff mode or a gene or regions target (motif_analysis.Rmd) |
| `multiqc.sh` | MULTIQC: one report over the QC files |
| `my_theme.R` | plot theme for the R figures (used by QC_FRAGMENT_SIZES_PLOT) |
| `my_theme_py.py` | Python plot theme; no step uses it at present |
| `norm_tmm.R` | NORM_TMM: TMM factors for one group, MA plots (the statistician's script) |
| `peak_merge.sh` | PEAK_MERGE: sort-bed then `bedops --merge` of the listed peak files |
| `qc_chrom_fraction.py` | QC_CHROM_FRACTION: percent of mapped reads on each removed contig |
| `qc_fragment_sizes_plot.R` | QC_FRAGMENT_SIZES_PLOT: density and boxplot SVGs from the two fragment tables |
| `qc_fragment_sizes_table.R` | QC_FRAGMENT_SIZES_TABLE: fragment-size distribution and summary per BAM |
| `qc_frip.sh` | FRIP: fraction of reads in peaks for one sample against its own peak file |
| `qc_frip_table.py` | QC_FRIP_TABLE: merges the per-sample FRiP files |
| `qc_peak_stats.py` | QC_PEAK_STATS: peak counts, widths and chromosome summaries for every peak file |
| `qc_tss_enrichment.sh` | TSS_ENRICHMENT: computeMatrix and plotProfile around the TSS for one bigWig |
| `qc_tss_table.py` | QC_TSS_TABLE: merges the per-sample TSS files with the samplesheet metadata |
| `tn5_shift.sh` | TN5_SHIFT: alignmentSieve --ATACshift, then sort and index (pysam) |
| `trim.sh` | TRIM: Trimmomatic paired-end trimming with the ATAC settings |
| `txdb_build.R` | TXDB_BUILD: TxDb (gene and transcript models) from the GTF, saved as SQLite |

Argument keys of the downstream steps: `bam_prealigned_args`, `txdb_build_args`, `diff_accessibility_args`,
`diff_annotate_args`, `motif_enrichment_args` and `motif_common_args` (the same rule as the other `<script>_args` keys).

Plot-only changes (colours, labels, sizes): replace the plot script and/or `my_theme.R`.
Changes to what is computed: replace the compute script and keep its output names. Changes
to normalisation factors, groups, thresholds or models are scientific decisions: confirm them
first.

## Inputs are never copied

- FASTQs, BAMs and reference files stay where they are. The pipeline reads them in
  place, never moves, copies, renames or deletes them, and `config.yaml` keeps pointing
  to them. Keep the path stable (if you move the files, update `config.yaml` and the CSVs).
- Write **Windows paths** in full (`D:\FASTQ\run7\x_R1.fastq.gz`) or relative to the
  project folder (`original_data\fastq\x_R1.fastq.gz`). Forward slashes work too.
- A drive named in `config.yaml` or in the CSVs is mounted **read-only** in the
  analysis (as `/mnt/<letter>`, so `D:\x\y` is read as `/mnt/d/x/y`); a drive that is not
  connected stops the run with a message. Nothing is copied from it.
- **Network paths** (`\\server\share\...`) are refused: Docker Desktop cannot read them.
  Map the share to a drive letter (File Explorer: This PC > Map network drive), or copy
  the files to a local or USB drive, and write the path with that letter.
- The launcher records the paths as you wrote them and as the analysis used them in
  `metadata\host_paths.tsv`, and logs the mapping at the start of the run.

## Downstream analyses

Three switches add the differential, annotation and motif steps after the peaks, before MultiQC. All are off by default
(`run_differential`, `run_annotation`, `run_motifs`). With all three off, none of them is run or checked. Each comparison is a
case side against a control side of a samplesheet column (`group_column`, default `condition`), with its own peak groups.
The steps are the same as in the cluster version (`..\HPC\`).

| Switch | Runs | One task per |
|---|---|---|
| `run_differential` | DIFF_ACCESS (csaw, edgeR), and TXDB_BUILD when promoters are computed | comparison that computes a region set |
| `run_annotation` | DIFF_ANNOTATE (ChIPseeker), and TXDB_BUILD | comparison, on its peaks table |
| `run_motifs` | MOTIF_ENRICHMENT (monaLisa), MOTIF_COMMON | comparison and `motif_regions` set; `motif_targets` row; region set with two or more non-pooled comparisons |

The peaks table is always computed by DIFF_ACCESS unless `diff_tables` gives it. `diff_regions` lists the region sets DIFF_ACCESS
computes (`peaks`, or `peaks,promoters`). A comparison that has a table in `diff_tables` is not recomputed for that region set.
Use `diff_tables` to start the annotation or the motifs from tables that already exist. Use `peak_files` to add existing peak
files as peak groups.

Pre-aligned libraries: a samplesheet row with a final BAM in the `bam` column (and no FASTQ) joins the analysis like a FASTQ
library (BAM_PREALIGNED). Its BAM is assumed final. Its row needs a `tmm_group` or `peak_group` to take part in those steps.

### comparisons.csv

```
name,case,control,group_column,block_column,case_peaks,control_peaks,pooled,fixed_dispersion
BT54_high_vs_low,BT54_high,BT54_low,condition,,BT54_high,BT54_low,false,
```

| Column | Meaning |
|---|---|
| `name` | Output prefix, matches `^[A-Za-z0-9][A-Za-z0-9_.-]*$`, unique |
| `case`, `control` | Comma lists (quoted) of values of `group_column`. The active samples with those values form the sides. A blank side is allowed only when the comparison needs no DIFF_ACCESS task |
| `group_column` | A samplesheet column; default `condition` |
| `block_column` | Optional samplesheet column with a value for every sample of the comparison. Adds a covariate: `~ block + grouping` (legacy gd2: the cell line) |
| `case_peaks`, `control_peaks` | Comma lists of peak group names: MACS3 peak groups (needs `run_macs3`) or names in `peak_files.csv`. Default: the distinct `peak_group` values of the side |
| `pooled` | `true` for a comparison that pools several lines. It is analysed like the others and left out of MOTIF_COMMON |
| `fixed_dispersion` | Optional number: the LRT dispersion of this comparison (default `diff_fixed_dispersion`) |

No sample may be on both sides. Each side needs at least one sample when the comparison is computed.

### peak_files.csv and diff_tables.csv

peak_files.csv:

```
name,path
extra,D:\peaks\extra.narrowPeak
```

diff_tables.csv:

```
name,regions,path
BT54_high_vs_low,promoters,D:\results\differential\BT54_high_vs_low_promoters.csv
```

`peak_files`: names unique and not a `peak_group` of the samplesheet; the file exists and ends in `.narrowPeak`, `.broadPeak`
or `.bed`. `diff_tables`: the name is a comparison; `regions` is `peaks` or `promoters`; the file exists; each (name, regions)
pair is listed once.

### motif_targets.csv

```
name,type,target,comparison,label,min_peaks
ST8SIA1_gene,gene,ST8SIA1,gd2_high_vs_low,ST8SIA1 peaks,
ST8SIA1_cCRE,regions,D:\refs\st8sia1_ccres_hs1.bed,gd2_high_vs_low,ST8SIA1 cCRE,
```

`type` is `gene` (the peaks of the `all_sites` annotation table whose `motif_gene_column` contains the symbol) or `regions`
(peaks that overlap the BED, against the other peaks). `comparison` is the universe: its `all_sites` table. It needs
`run_annotation`. The `name` must differ from every comparison name, because both write under `results\motif\`. A regions BED
must already be in the pipeline genome's coordinates (for `hs1`, T2T-CHM13): the liftover of a BED from another build is a
project pre-step, not done by the pipeline. `label` must not contain quotes, backticks, backslashes or dollar signs.

### Startup checks

Downstream values are checked at startup, whenever their file or switch is set, with a message per problem (the run does not
start): an unknown or duplicate name, a column that is not in the samplesheet, a value of `group_column` that no active sample
has, a sample on both sides, a missing peak group or peak file, a missing table or BED, a `diff_regions` other than `peaks` or
`peaks,promoters`, a `motif_regions` set without its table (promoters and peaks need the table of the comparison; the annotation
sets need `run_annotation`), a missing `gtf` when promoters or annotation need the TxDb, and a wrong number, boolean or
enumerated value.

Errors inside a step: exit 3 is a known data problem (for example fewer than one BAM on one side, no peak file, promoters
without a TxDb, a missing BSgenome package, a differential table with no rows). A **skipped** motif analysis is not an error:
when a bin has fewer regions than `motif_min_bin_peaks` (or fewer than `motif_target_min_peaks` in a target analysis), the step
writes only `<name>.skipped.txt`, with the reason and the counts, and the run goes on.

### Outputs

```
results\annotation\txdb\                 txdb.sqlite (TXDB_BUILD); log\
results\differential\                    <name>.csv (peaks), <name>_promoters.csv, <name>.norm_factors.tsv, <name>.summary.tsv,
                                         <name>.replicates.csv, <name>.diag_*.svg, <name>.diagnostics.pdf (plot_format pdf); log\
results\differential\bed\                <name>.bed, <name>_upregulated.bed, <name>_downregulated.bed, <name>_promoters*.bed
results\annotation\                      <name>_all_sites.csv, <name>_promoter_sites.csv, <name>_non-promoter_sites.csv,
                                         <name>_distal_sites.csv and their .bed files, <name>_peak_distribution.svg,
                                         <name>.annotation_summary.tsv; log\
results\motif\<name>\<region>\           <name>_motif_analysis.csv, <name>.motif_summary.csv, heatmaps, similarity, GC and
                                         dinucleotide plots (diff mode); log\
results\motif\<name>\                    the same for a motif_targets row (gene or regions); <name>.skipped.txt when skipped
results\motif\common\<region>\           common_significant_motifs.csv, venn_common_motifs.svg, motif_analysis_summary.csv; log\
```

Every step writes its command to `log\<tag>_<step>_command.sh`.

### Memory, time and the Docker budget

- The cluster's requests for the new steps are DIFF_ACCESS 16 GB, MOTIF_ENRICHMENT 32 GB, DIFF_ANNOTATE 8 GB, TXDB_BUILD 4 GB,
  MOTIF_COMMON 2 GB and BAM_PREALIGNED 2 GB (per attempt, with retries). They are estimates, not measured on this computer. On a
  Docker computer each request is capped to the run's budget (see "Resources"), which is at most 14 GB here
  (`max_useful_memory_gb`). A comparison with many bins needs a computer with a budget of 16 GB or more; the 14 GB cap does not
  reach that, so such a comparison gets the whole budget and can run out of memory. A step killed for lack of memory is retried
  with more memory, never above the budget.
- In the cluster port's tests, a motif analysis of the full 190,000-peak gd2 table ran out of memory at 13 GB in Docker. Run such
  an analysis on the cluster.
- **Status:** the downstream steps are new to this Docker version and have not yet been run end to end here. The scripts were
  checked against the legacy results, and the wiring with `-stub-run`, in the cluster port's validation (in the same image,
  `ghcr.io/genesig/omics:2026.10.3`). The run times and memory peaks on a
  Windows or Mac computer are not known yet.

### What the statistics can and cannot say

These are the legacy method's limits, ported unchanged. State them when you report a result.

- **Fixed dispersion** (fewer libraries than `diff_qlf_min_libs`, for example one library per side): there is no replication, so
  the variance is assumed (`diff_fixed_dispersion`), not estimated. PValue and FDR are then not valid error rates. Rank regions by
  logFC only; do not call them significantly differential.
- **Blocked or pooled comparisons** (legacy gd2): the unit of replication is the cell line, not the library. The additive block
  assumes one common high-versus-low effect in all lines.
- **Tested regions come from the data:** the union of the peaks of both sides, called on pooled replicates. Regions open only in
  the better-powered group are over-represented. FDR is over this region set only.
- **Normalisation:** TMM assumes that most of the genome does not change. It does not correct efficiency bias. If the compared states
  differ globally in accessibility, check the MA plot (`NAME.diag_ma.svg`) and `NAME.norm_factors.tsv` first.
- **Promoters:** the promoter windows of neighbouring genes overlap, so promoter rows are not independent.
- **Motifs:** the bins are cut on logFC without an FDR filter. Each PWM hit counts as an independent trial, so the p-values are
  anti-conservative, and the JASPAR matrices are correlated. Read the ranking and log2enr, not the count of significant motifs.
  `motif.name` is not unique in JASPAR2024: `motif_common` keys on the name as the legacy did, so a few motifs can be missed (a
  warning is printed).
- **Annotation:** a gene is the nearest TSS with no distance limit. An annotation is not evidence that a region regulates that gene.
- **Decisions:** the comparisons, `block_column`, `pooled`, `fixed_dispersion`, the thresholds (`diff_lfc`, `motif_lfc`, `diff_fdr`)
  and the motif targets are scientific decisions. Confirm them with the user, and with the statistician, before a run.

## Resources: how much of the computer an analysis uses

Each run has one **budget** of memory and CPU cores, and never goes over it. Docker
enforces it as a hard limit, so a run cannot freeze Windows by running out of memory.

**Memory.** The budget is the smallest of:
- 75% of the computer's memory, but at least 3 GB left for Windows;
- Docker Desktop's memory minus 0.5 GB (the Docker VM's own share);
- 14 GB, the most this pipeline can use (`max_useful_memory_gb` in `docker\requirements.conf`).

**CPU cores.** All cores but one (all but two from 8 cores up), so the desktop stays
responsive. `max_cpus` in `config.yaml` can lower this.

**Per-step memory requests.** Each step asks for a share of the budget, set per step in
`conf\docker_tiers.config`. The largest are ALIGN (4.5 GB, fixed, for the bowtie2 index it
loads), BAM_FINALIZE (8 GB), BAM_FILTER (7.5 GB), MACS3 (5.5 GB) and CONTAMINATION (5 GB);
each is capped to the budget. Steps run side by side as far as the budget allows; the rest
wait their turn. A step that is killed for lack of memory is retried with more memory, never
above the budget, up to three times.

`config.yaml` can lower the budget: `max_memory_gb` (decimals allowed, e.g. `2.5`) and
`max_cpus`. Neither can go above what Docker Desktop allows (the launcher warns).
Results do not change, only speed. Each run prints its budget at the start.

**Minimum.** The run refuses to start, with an explanation, if the budget is below 5 GB
(`min_memory_gb` in `docker\requirements.conf`; about 1 GB of it is kept for Nextflow itself).
With the Docker limit set as in the table below, that means a computer with **8 GB of memory
or more, as Windows reports it** (the budget is 3 GB less than the computer's memory, or 75%
of it if that is less, so a computer that reports less than 8 GB is refused). A computer with
8 GB gets a 6 GB Docker limit and a budget of about 5 GB, which meets the minimum: the run goes
ahead, slowly. The minimum and `max_useful_memory_gb` are placeholders until a pilot run on a
real computer confirms them.

**Docker's memory limit (WSL 2 memory, Windows only).** With WSL 2, Docker Desktop
may use only a part of the computer's memory, set by the `memory=` line of
`%USERPROFILE%\.wslconfig`. For an ATAC-seq run the launcher sets that limit to suit the
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
  back the same way, if nothing else is running. The file `%USERPROFILE%\.wslconfig.run_atacseq.bak`
  exists while the changed limit is in use. If the window was closed before the restore,
  the next `run_atacseq` command that finds Docker running and nothing else running puts
  the old value back.
- `--no-wsl-resize` (on `run`) or `wsl_resize: false` in `config.yaml` leaves the limit alone.
  The budget is then computed from the limit Docker already has.

**Examples** (the rules above, for a computer with the Docker limit set as in the table):

| Computer | Docker limit | Memory budget | CPU cores |
|---|---|---|---|
| 8 GB, 4 cores | 6 GB | about 5 GB: runs, slowly (the largest steps one at a time) | 3 |
| 12 GB, 8 cores | 8 GB | 7.5 GB | 6 |
| 16 GB, 8 cores | 8 GB | 7.5 GB | 6 |
| 32 GB, 12 cores | 16 GB | 14 GB | 10 |

**Keep memory free.** The pipeline caps its own memory, but other programs share it. Close
browsers with many tabs, large IDE sessions and other analyses before a long run. While a
run is going, the launcher watches Windows' free memory and makes Docker drop its file
cache when it gets low.

## How long a run takes

Measured with Docker Desktop on a 20-thread PC (budget 6 cpus and 8 GB for the tests), full depth, two libraries of
about 45 million pairs each (atac_av_6349, 11.4 GB of FASTQ): ALIGN took 38-48 min per 10-million-pair chunk
(bowtie2 on 4 threads, about 390% cpu; it is the step that decides the run time), TRIM 22-26 min per library,
raw FastQC 21-29 min per file with a copy of the FASTQ in Docker's scratch folder first (1h13 to 1h30 per file without
the copy when the FASTQs were read from the hard disk D: while other programs used it), BAM_FILTER 2.5-4 min, DEDUP 1.5-2 min,
BLACKLIST_FILTER 10-14 min (one thread), BAM_FINALIZE 2.5-3 min, MACS3 2-3 min, BIGWIG 2 min, TSS_ENRICHMENT 7 min.
The nine alignment chunks of the two libraries run one after the other at 4 cpus (the test budget), so ALIGN alone takes about 6 h there; the run was resumed several times, so its total wall time was not measured. Alignment is about 85% of the cpu time, so a run takes about
(number of read pairs / 10 million) x 160 cpu-minutes of bowtie2, divided by the cpus given to ALIGN (the budget's cpus less 2).
For the 39-library AV dataset (about 347 GB of FASTQ, about 2.7 billion pairs) that is about 44,000 cpu-minutes for ALIGN
and about 7,000 for the other steps: about 2 to 2.5 days of wall time on a PC with 20 threads (16 cpus for ALIGN), if the
FASTQs are on a fast disk. This is an estimate scaled from the two measured libraries, not a measured run of the dataset.

## Disk space

- **Plan about 5.2 times the total size of your FASTQ files** in all: about 4.8 times for the work
  files, on the drive that holds Docker Desktop's disk image, and about 0.4 times for the results,
  on the project's drive (with the default `work_volume: true`). With `work_volume: false` the work
  files are in the project folder, so the project's drive needs about 5.2 times. These are
  estimates: the factors are the values kept in `docker\requirements.conf` (`disk_factor=5.2`,
  `disk_results_factor=0.4`) from the lab's chromatin-profiling Docker pipeline, and they stay estimates until the
  first full ATAC-seq run gives the real sizes. The launcher computes the same numbers before
  every run: it adds up the FASTQs named in the samplesheet and prints how much room the analysis
  may need on the project's drive (results) and on Docker's disk (work files). It warns when the
  free space on Docker's disk is below the work-file figure (the run stops only if that disk
  actually fills up, and can then be resumed). Results need room on the project's drive as well;
  the run stops before it starts only when less than 20 GB is free there.
- **To estimate:** in File Explorer, select the FASTQ folder of your samplesheet and read its
  size. A full dataset can need 100 GB or more, so check before you start.
- What uses the space: the work volume on Docker's disk (the working files of every step, kept
  so a re-run resumes), the trimmed FASTQs, aligned BAMs and filtered BAMs in
  `intermediate_files\` (only with `publish_intermediates: true`; the default on Windows is
  `false`, so the folder stays empty), and the results in `results\`. The raw FASTQs are not
  copied.
  With the work volume, files in `results\` are copies. With `work_volume: false` they are
  hard links to the work folder when the disk allows (the run log says so), so big files are
  not stored twice. Never edit a published file in place; replacing or deleting it is fine.
- To free space: `run_atacseq clean` (results are kept), or `keep_work: false` (deletes
  the work files after a successful run).
- **The project can be on a hard disk.** With the work volume, only the results (and
  `intermediate_files\` with `publish_intermediates: true`) are written to the project folder, so
  a hard disk or a USB drive slows the steps much less than before (the launcher warns when the
  project is on a hard disk). Docker's disk image is what matters most for speed: keep it on an
  internal SSD, or at least on a drive with room (see "Good to know").

## Software versions

The image `ghcr.io/genesig/omics` (tag in `docker\image.txt`, `2026.10.3`) contains the tools
below. Versions are those pinned in the image's package locks (`PIPELINES\docker\omics\env\*.lock.txt`);
the setup step prints the main ones. R and Bioconductor are newer than on the cluster (R 4.4.3 here;
the cluster's R image has R 4.3.3), so numbers can differ slightly from results made with the
older containers. See `ANALYSIS.md` for the cluster versions side by side.

| Tool | Version in the image |
|---|---|
| Nextflow (inside the image) | 26.04.6 |
| bowtie2 / samtools | 2.5.5 / 1.23.1 |
| Trimmomatic / FastQC / MultiQC | 0.40 / 0.12.1 / 1.35 |
| sambamba / bedtools / bedops | 1.0.1 / 2.31.1 / 2.4.42 |
| deepTools / MACS3 | 3.5.6 / 3.0.4 |
| R (`Rscript`) / csaw / edgeR / Rsamtools | 4.4.3 / 1.40.0 / 4.4.0 / 2.22.0 |
| ChIPseeker / monaLisa / txdbmaker / JASPAR2024 (downstream steps) | 1.42.0 / 1.12.0 / 1.2.0 / 0.99.6 |
| Python (`python3`) / pysam | 3.12 / 0.24.1 |

## Differences from the cluster

The steps, tools and results are those of the cluster version (`HPC`). Four things differ here:

- **One image holds all tools.** There are no per-step containers and no scheduler: every
  step runs in the same Docker image with Nextflow's local executor, inside one budget of
  memory and cores for the whole run (see "Resources").
- **Intermediates are not published by default.** `publish_intermediates` (default `false`)
  decides whether the trimmed FASTQs, aligned BAMs and filtered BAMs are copied to
  `intermediate_files\`. The cluster always publishes them and has no such key.
- **Inputs are read where they are.** Drives named in the config are mounted read-only;
  network paths are refused (see "Inputs are never copied"). The cluster reads the same
  files from its own file system.
- **Downstream steps use the image's R packages.** The cluster runs them in its own R image
  (`omics-r_1.0.sif`); here they run in `ghcr.io/genesig/omics`, with the versions above, and
  their memory is capped to this computer's budget (see "Downstream analyses").

## Good to know

- **Where Docker keeps its data:** Docker Desktop stores images in a disk file on `C:` by
  default. The work volume lives on the same disk, so it needs room for about 4.8 times the
  FASTQ total (an estimate) and is best on an internal SSD. If `C:` is short of room or not an SSD, move it
  in Docker Desktop > Settings > Resources > Advanced > *Disk image location* (the
  launcher assumes no particular drive).
- **Reproducibility:** every tool runs at a pinned version inside the image, and the
  image tag is in `docker\image.txt`. Every step's `log\` folder holds the exact command
  that ran (`<sample_id>_<step>_command.sh`). `metadata\pipeline_info\` holds the Nextflow
  execution report, trace and timeline. A resumed run keeps the previous run's trace, report and
  timeline next to it, named with that run's timestamp; the start of each run's log names the image.
- **Stopping:** `run_atacseq stop` stops the analysis cleanly; finished steps are kept and the
  next run resumes. Closing the window does not stop it.
- **Updating the pipeline:** `git pull` in the `ATACSEQ` folder of the clone.
  If the update changed the image tag, the next `setup` or `run` downloads the new image.
- **Example projects** in `examples\` (`config.yaml`, `samplesheet.csv`, `peak_merges.csv`,
  and for the downstream steps `comparisons.csv` and `motif_targets.csv`)
  are the starting templates that `init` copies. `examples\projects\` holds the two dataset
  examples (AV_ATACSEQ_2026 and AM_ATACSEQ) with the lab's legacy paths; copy one and edit its
  paths before you use it.

For how each step works (the tools, the normalisation, the filters and the output layout), see
`ANALYSIS.md` in this folder. The analysis is the same as on the cluster, except for the
differences listed above; only how it is launched differs.
