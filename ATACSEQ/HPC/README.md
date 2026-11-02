# ATACSEQ: ATAC-seq pipeline (Rephill Lab)

A Nextflow pipeline for paired-end ATAC-seq on UPenn PMACS (LSF + Apptainer). It replaces the legacy `step_*` Python and bash
pipelines of two projects, AV_ATACSEQ_2026 (newer; its settings are the defaults) and AM_ATACSEQ (older; its settings are
toggles), with one reusable pipeline, one `config.yaml` and two small CSVs per project. It fixes the known bugs of the legacy
scripts on the way; each fix has a config toggle that reproduces the legacy behaviour.

- For agents working in a project: `agent/AGENTS.md`

Paired-end only. Single-end FASTQ is not supported. Not in scope: deepTools heatmaps and profiles, ChIPseqSpikeInFree, and the
legacy ENCODE and diagnostic helpers.

## What it does

Eighteen numbered steps (3a, 3b and 8b are lettered sub-steps).
Steps in brackets run only when their switch or key is set.

```
FASTQ pairs ──> FASTQC_RAW
            └─> [TRIM (trimmomatic)] ─> [FASTQC_TRIM]
                 ─> [COUNT_PAIRS (read pairs; when align_chunk_pairs > 0)]
                 ─> ALIGN (bowtie2 | samtools view -b; one task per chunk of read pairs) ─> ALIGN_MERGE (samtools cat; summed log)
                      ─> ALIGNMENT_SUMMARY (warning only)
                 ─> BAM_FILTER (sort; raw stats; MAPQ and flag filter; [remove_chroms])
                      ─> [QC_CHROM_FRACTION]
                      ─> DEDUP (sambamba markdup; -r when remove_duplicates)
                           ─> [BLACKLIST_FILTER (pair-aware)] ─> BAM_FINALIZE (final BAM, final stats, orphan check)
                 ─> [CONTAMINATION (contaminant index)]
final BAM ─> [TN5_SHIFT (alignmentSieve --ATACshift)] ─> analysis BAM (S.tn5.bam, or S.final.bam)
analysis BAM ─> QC_FRAGMENT_SIZES_TABLE ─> QC_FRAGMENT_SIZES_PLOT
analysis BAM, per tmm_group ─> NORM_TMM (csaw / edgeR TMM) ─> scale-factor table per group
analysis BAM + multiplier, per sample with a tmm_group ─> BIGWIG (bamCoverage) ─> [TSS_ENRICHMENT] ─> QC_TSS_TABLE
analysis BAM, per peak_group (replicates pooled, no control) ─> MACS3 ─> [PEAK_MERGE from peak_merges.csv]
MACS3 peaks ─> QC_PEAK_STATS;  analysis BAM + own peak group ─> FRIP ─> QC_FRIP_TABLE
all QC outputs ─> MULTIQC
```

The order of the duplicate and blacklist steps is the legacy order (markdup, then blacklist). What changed is the blacklist method:
both mates of a pair are removed when either overlaps, and MAPQ and sequences are kept (legacy converted the BAM to BEDPE and
back).

Chunked alignment: COUNT_PAIRS counts each library's read pairs, and ALIGN aligns `align_chunk_pairs` pairs per task (bowtie2
`-s`/`-u`, no temporary FASTQs). ALIGN_MERGE joins the chunk BAMs in order and sums their bowtie2 logs, so the merged BAM has the
same records as one unchunked run (`align_chunk_pairs: 0`). Downstream steps start per group: NORM_TMM and MACS3 start for a group
as soon as its own libraries are aligned, and FRIP and PEAK_MERGE wait only for their own peak groups. The outputs are the same as
before; only the start times change.

Processes: FASTQC (run twice as FASTQC_RAW and FASTQC_TRIM), TRIM, COUNT_PAIRS, ALIGN, ALIGN_MERGE, ALIGNMENT_SUMMARY, BAM_FILTER,
QC_CHROM_FRACTION, DEDUP, BLACKLIST_FILTER, BAM_FINALIZE, CONTAMINATION, TN5_SHIFT, QC_FRAGMENT_SIZES_TABLE, QC_FRAGMENT_SIZES_PLOT,
NORM_TMM, BIGWIG, TSS_ENRICHMENT, QC_TSS_TABLE, MACS3, PEAK_MERGE, QC_PEAK_STATS, FRIP, QC_FRIP_TABLE, MULTIQC.

Key points:

- The TMM group and the peak group are exact samplesheet columns. Nothing is grouped by file name.
- A sample of a TMM or peak group with no final BAM fails the run. It is never dropped silently.
- Samples with a blank `tmm_group` get no bigWig and no TSS enrichment. Samples with a blank `peak_group` get no peak calling and no FRiP.
- Every step runs a script from `bin/` (or the project's replacement, see "Customisation"). Startup validation checks the
  samplesheet, the peak merges and the reference paths before any job is submitted.

## Install and launcher

The pipeline is in `/project/rephilllab/ATACSEQ`. Labmates reach it through one launcher on the shared PATH:
`/project/rephilllab/bin/run_atacseq.sh`, a symlink to `ATACSEQ/run_pipeline.sh`. The symlink is created by the lab lead. Until it
exists, call `/project/rephilllab/ATACSEQ/run_pipeline.sh` directly.

The launcher pins Nextflow 26.04.6 (`NXF_VER`). It downloads the binary once into `.nextflow-bin/` (shared by everyone who runs
the pipeline) if it is missing. It never runs the pipeline on the login node: `run` submits one LSF head job (queue rhel9, 8 GB,
168 h) that runs `nextflow run /project/rephilllab/ATACSEQ -profile pmacs -params-file config.yaml -resume`. The head job runs each
step as its own LSF job in its container. The head job's work directory is `/project/rephilllab/nf_work/$USER/<project directory name>/`.

```bash
mkdir my_atacseq && cd my_atacseq
run_atacseq.sh init          # scaffolds the project, copies config.yaml, samplesheet.csv, peak_merges.csv, examples/, AGENTS.md, CLAUDE.md; git init; runs nothing
# edit samplesheet.csv, peak_merges.csv (optional) and config.yaml (genome_index, blacklist, tss_bed, effective_genome_size, ...)
run_atacseq.sh               # submits the head job; reads ./config.yaml (or: run_atacseq.sh other.yaml)
bjobs -u $USER
tail -f metadata/nf_logs/nextflow_<timestamp>.out
```

Setting `help: true` in config.yaml prints the usage and exits. `lsf_watchdog: false` in config.yaml turns off the watchdog
(`lib/lsf_watchdog.sh`), which stops a run whose LSF job was aborted at dispatch instead of letting Nextflow wait for it forever.

## Project layout

```
<project>/
├── config.yaml              the run's settings (flat key: value; every key in the table below)
├── samplesheet.csv          one row per library
├── peak_merges.csv          optional: one row per merged peak set
├── comparisons.csv, peak_files.csv, diff_tables.csv, motif_targets.csv   optional: the downstream analyses (see below)
├── original_data/fastq/     FASTQs, only if fastq_dir is used (otherwise the samplesheet paths are used as given)
├── metadata/                nf_logs/ (LSF head-job logs), .nextflow.log, pipeline_info/ (report, trace, timeline, DAG)
├── intermediate_files/      trimmed FASTQ, aligned and filtered BAMs
├── utils/                   this project's replacement scripts (see "Customisation")
├── results/<step>/          final outputs, one folder per step
├── examples/                starter files that init fetched, for reference
└── .nextflow/               Nextflow's resume cache (cannot be moved out of the launch directory)
```

## samplesheet.csv

Header (column order fixed):

```
sample_id,fastq_1,fastq_2,condition,replicate,tmm_group,tmm_group_legacy,peak_group,exclude,bam
```

| Column | Required | Meaning |
|---|---|---|
| `sample_id` | yes | Unique library name, used verbatim in file names. Matches `^[A-Za-z0-9][A-Za-z0-9_.-]*$`. |
| `fastq_1`, `fastq_2` | yes | Absolute paths of the R1 and R2 FASTQ. |
| `condition` | no | Label for the QC tables (fragment sizes, TSS, FRiP). |
| `replicate` | no | Label for the QC tables. |
| `tmm_group` | no | TMM group (exact values). Blank = not normalised, no bigWig, no TSS enrichment. Read when `tmm_group_column` is `tmm_group`. Needs at least 2 active samples per group. |
| `tmm_group_legacy` | no | Second TMM group column. Read when `tmm_group_column` is `tmm_group_legacy` (legacy file-name membership). |
| `peak_group` | no | MACS3 peak group (exact values). Replicates are pooled, no control. Blank = no peak calling and no FRiP. |
| `exclude` | no | `true`, `1` or `yes` drops the row from the whole run. Blank, `false`, `0` or `no` keeps it. Other values are a startup error. |
| `bam` | no | A final BAM of the library instead of FASTQ (BAM_PREALIGNED). Then `fastq_1` and `fastq_2` must be empty. The file must end in `.bam` and exist. The BAM is assumed final (filtered, duplicate-removed, blacklist-filtered). Each active row has FASTQ or a BAM, not both. |

A missing optional column is blank for every row. Any other column name is a startup error. The file may carry a UTF-8 byte
order mark. Quoted fields are supported.

Startup errors (before any job) name the problem: unknown or missing columns, duplicate `sample_id`, a bad id, a missing FASTQ, the
same file as R1 and R2, two FASTQs with one file name, a TMM group with one sample, `tmm_group_column` not a column, every row excluded,
an invalid value for a switch or number, a missing bowtie2 index, a missing `blacklist`, `chrom_sizes` or `gtf` file, a missing
`tss_bed` when TSS enrichment runs with TMM groups, a missing `effective_genome_size` when MACS3 runs with peak groups, a `bam` that
is not a `.bam` file or is missing, a row with both FASTQ and a BAM, and the downstream checks (see "Downstream analyses").

## peak_merges.csv

```
name,peak_groups
merged_peaks,"peak_1,peak_2"
```

One merged BED per row: `results/peaks/merged/<name>.merged.bed` (bedops merge of the listed groups' MACS3 peak files). `name` is
unique and matches the id pattern. `peak_groups` is a quoted comma list; each group must have an active sample. Setting `peak_merges`
to a file needs `run_macs3: true`. Use `peak_merges: null` for none.

## config.yaml keys

Every key has a default in `nextflow.config`. A project lists only the keys it changes. The two example datasets set every key.

| Key | Default | Meaning |
|---|---|---|
| `project_name` | `my_atacseq` | Prefix of the project-level outputs; the LSF job name is `atacseq_<project_name>_<timestamp>` |
| `samplesheet` | `samplesheet.csv` | Library table |
| `peak_merges` | `null` | Peak merge table; `null` = none |
| `fastq_dir` | `null` | A directory of FASTQs that the launcher moves into `original_data/fastq/`; `null` = use the samplesheet paths |
| `outdir` | `./results` | Output root |
| `metadata_dir` | `./metadata` | Run bookkeeping (report, trace, timeline and DAG paths) |
| `intermediate_dir` | `./intermediate_files` | Trimmed, aligned and filtered intermediates |
| `publish_mode` | `copy` | `copy`, or `symlink` (same filesystem only) |
| `lsf_exclude_hosts` | `''` | pmacs profile only: comma list of LSF hosts never used (`bsub -R "select[hname!=...]"`), for a node where every container start fails (seen on a compute node, 2026-10-27: "FIPS mode is enabled, but the required OpenSSL backend is unavailable") |
| `genome_index` | `null` | bowtie2 index prefix of the genome. Required. |
| `chrom_sizes` | `null` | chrom.sizes of the genome. Recorded and checked if set; no step uses it |
| `blacklist` | `null` | Blacklist BED for pair-aware removal. `null` or empty = skip BLACKLIST_FILTER |
| `gtf` | `null` | Gene annotation GTF. Recorded and checked if set; no step uses it |
| `tss_bed` | `null` | TSS BED for TSS enrichment. Required when TSS enrichment runs and a sample has a `tmm_group` |
| `effective_genome_size` | `null` | MACS3 `-g`. Required when MACS3 runs and a sample has a `peak_group` |
| `contamination_index` | `null` | bowtie2 index prefix of contaminant genomes. `null` or empty = skip CONTAMINATION |
| `run_fastqc` | `true` | FastQC on raw reads and, with `trim`, on trimmed reads |
| `trim` | `true` | Trimmomatic before alignment. `false` aligns the raw FASTQs |
| `run_tss_enrichment` | `true` | TSS enrichment on the bigWigs (needs `tss_bed`) |
| `run_macs3` | `true` | MACS3 per peak group, with peak statistics and FRiP |
| `run_multiqc` | `true` | MultiQC report |
| `tn5_shift` | `false` | `true` applies deepTools `alignmentSieve --ATACshift`; the shifted BAM feeds every later step |
| `remove_duplicates` | `true` | `sambamba markdup -r`. `false` marks duplicates and keeps them |
| `help` | `false` | Print the usage and exit |
| `trim_adapter_set` | `NexteraPE-PE` | Trimmomatic adapter file name (without `.fa`), located inside the container |
| `trim_illuminaclip` | `'2:30:10:8:True'` | ILLUMINACLIP parameters after the adapter file |
| `trim_steps` | `'LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20'` | Further Trimmomatic steps, space separated; may be empty |
| `trim_minlen` | `25` | Trimmomatic `MINLEN` |
| `bowtie2_mode` | `local` | `local` (`--local --very-sensitive-local`) or `end_to_end` (`--end-to-end --very-sensitive`) |
| `bowtie2_minins` | `10` | bowtie2 `-I`. `null` = not passed (legacy) |
| `bowtie2_maxins` | `2000` | bowtie2 `-X` |
| `bowtie2_args` | `''` | Further bowtie2 arguments, after the fixed flags. Legacy `'-N 1 -t'` |
| `align_chunk_pairs` | `10000000` | Read pairs per ALIGN task (chunked alignment, joined by ALIGN_MERGE). Integer, `0` = no chunking (one ALIGN task per library); negative or non-integer is a startup error |
| `alignment_min_pct` | `20` | Overall alignment rate below which ALIGNMENT_SUMMARY writes a warning (never fails) |
| `mapq` | `30` | Minimum MAPQ kept (`samtools view -q`). `0` reproduces the legacy behaviour |
| `remove_chroms` | `chrM` | Comma list of contigs removed from the BAM. `''` keeps all contigs (legacy) |
| `tss_flank` | `2000` | TSS flank, bp on each side of the TSS |
| `tmm_bin_size` | `10000` | TMM bin width, bp (at least 1000) |
| `tmm_max_frag` | `auto` | Maximum fragment size for TMM: `auto` (largest proper-pair fragment of the group) or a positive integer |
| `tmm_group_column` | `tmm_group` | Samplesheet column that defines the TMM groups (`tmm_group` or `tmm_group_legacy`) |
| `bigwig_bin_size` | `10` | bamCoverage `--binSize` |
| `bigwig_units` | `none` | `none` (bamCoverage `--normalizeUsing None`, legacy), `cpm` or `rpkm`. Also sets how the multiplier is written |
| `bigwig_extend_reads` | `true` | bamCoverage `--extendReads`: paired-end fragment coverage (legacy `false`: per-mate read coverage) |
| `bigwig_center_reads` | `false` | bamCoverage `--centerReads` (legacy `true`, a no-op without `--extendReads`; with it, read-length windows at fragment midpoints) |
| `bigwig_args` | `''` | Further bamCoverage arguments |
| `macs3_qvalue` | `0.05` | MACS3 `-q` |
| `macs3_min_length` | `250` | MACS3 `--min-length`. `null` = off |
| `macs3_max_gap` | `500` | MACS3 `--max-gap`. `null` = off |
| `macs3_fe_cutoff` | `2.0` | MACS3 `--fe-cutoff`. A value of 1 or less = off and no `.fe` in the file name |
| `macs3_call_summits` | `false` | MACS3 `--call-summits` (writes summit BEDs) |
| `macs3_broad` | `false` | MACS3 `--broad` (broadPeak and gappedPeak output) |
| `macs3_broad_cutoff` | `0.1` | MACS3 `--broad-cutoff`, used with `macs3_broad` |
| `macs3_args` | `''` | Further MACS3 arguments, placed on the command line by the module |
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
| `motif_match` | `exact` | Gene target matching: `exact` (whole gene symbol equals the nearest gene) or `substring` (legacy `grepl`; also matches other genes, e.g. FOS: FOSB, FOSL1, FOSL2) |
| `motif_cpus` | `1` | `1` = serial run with the legacy seed (the reproducible result); more uses MulticoreParam and is not bit-reproducible |
| `motif_common_padj` | `0.10` | MOTIF_COMMON: a motif is significant in a comparison at this adjusted P |
| `motif_common_min_lines` | `2` | MOTIF_COMMON: a motif is common when significant in at least this many comparisons |
| `motif_common_label_strip` | `null` | Regex removed from the comparison names to make the labels (AV example: `_high_vs_low.*`) |
| `plot_font` | `Helvetica` | Font family of the text in every downstream plot (svg keeps it editable) |
| `plot_format` | `svg` | `svg` or `pdf` (`pdf` = the legacy file format) |

Not Nextflow params: `lsf_watchdog` (read by the launcher; default on) and the script keys below. Any config.yaml key named after a
script in `bin/` replaces that script; `<script>_args` adds extra arguments to it (see "Customisation").

## Outputs

```
results/qc/fastqc/raw/, results/qc/fastqc/trimmed/   FastQC zip and html per FASTQ; log/ holds each command
results/qc/trimmomatic/                              Trimmomatic logs
results/qc/bowtie2/                                  bowtie2 alignment logs per sample (summed over chunks); log/ also holds the
                                                     COUNT_PAIRS, ALIGN (per chunk) and ALIGN_MERGE command files
results/qc/alignment/                                alignment_summary.tsv
results/qc/chrom_fraction/                           chrom_fraction.tsv (only when remove_chroms is set)
results/bam/stats/                                   S.raw.flagstat/idxstats (before the filter), S.final.flagstat/idxstats,
                                                     S.orphan_check.txt, S.sambamba_markdup.log
results/bam/final/                                   S.final.bam and .bai per library
results/bam/tn5/                                     S.tn5.bam and .bai (only when tn5_shift is true)
results/qc/contamination/                            S.contamination.tsv, S.contamination_summary.txt (only when contamination_index is set)
results/qc/fragment_sizes/                           fragment_sizes.tsv, fragment_size_distribution.tsv, fragment_sizes_density.svg,
                                                     fragment_sizes_boxplot.svg
results/scale_factors/                               <group>.scale_factors.tsv, <group>.ma_plot.png and .pdf, <group>.max_frag.txt
results/bigwig/                                      <sample>.<method>.<units>.bw, for example sample.tmm.none.bw
results/qc/tss_enrichment/                           S.tss_enrichment.tsv, S.tss_matrix.gz, S.tss_profile.svg, tss_enrichment.tsv
results/peaks/macs3/                                 <group>.qval<q>[.minlen<l>][.maxgap<g>][.fe<fe>][.broad]_peaks.narrowPeak (or broadPeak),
                                                     _peaks.xls, optional _summits.bed and _peaks.gappedPeak, .macs3.log
results/peaks/merged/                                <name>.merged.bed (only with peak_merges)
results/qc/peaks/                                    peak_stats.tsv
results/qc/frip/                                     S.macs3.frip.tsv, frip.tsv
results/qc/multiqc/                                  multiqc_report.html, multiqc_report_data/
results/bam/final/                                   (also) BAM_PREALIGNED outputs of the bam rows: S.final.bam, .bai, flagstat, idxstats
results/annotation/txdb/                             txdb.sqlite (only when run_annotation or a promoter table is computed)
results/differential/, results/differential/bed/     differential tables, plots and BEDs per comparison (run_differential)
results/annotation/                                  peak annotation tables, BEDs, plots and summary (run_annotation)
results/motif/<name>/<region>/, results/motif/<name>/, results/motif/common/<region>/   motif analyses (run_motifs)
intermediate_files/trimmed/                          S.R1.trim.fastq.gz, S.R2.trim.fastq.gz
intermediate_files/aligned/                          S.aligned.bam (unsorted, includes unmapped records; joined from the chunks)
intermediate_files/filtered/                         S.filtered.bam and .bai
metadata/pipeline_info/                              execution_report_<timestamp>.html, execution_trace_<timestamp>.txt, execution_timeline_<timestamp>.html, pipeline_dag_<timestamp>.html
metadata/nf_logs/                                    the LSF head-job stdout and stderr of each run
```

Every step's `log/` subfolder holds `<name>_command.sh`, the exact command that ran. The BIGWIG log files are named
`<sample>_<method>_bigwig_command.sh`. The DEDUP BAM and the BLACKLIST_FILTER BAM are not published (they are work files). With chunked
alignment the per-chunk ALIGN BAMs and logs are not published either; their commands are `<sample>.chunkNNN_align_command.sh`.

## Customisation

Every step runs one script from `bin/`. To replace a step's script for one project, copy it into the project's `utils/`, edit the
copy, and add a `config.yaml` key named after the script (file name without extension):

```yaml
qc_frip: utils/qc_frip.sh            # replace a step's script for this project only
fastqc_args: "--nogroup"             # or pass extra arguments to the shipped script
```

Scripts are task inputs, so `-resume` re-runs only the changed step and what depends on it. Each script's header documents its
arguments and output files; a replacement must keep that contract. A key that names no script, with a `.sh`, `.py` or `.R`
path as its value, is a startup error. Base names are unique across `bin/`, so `my_theme.R` and `my_theme_py.py` have different keys.

| Script (`bin/` file = config key) | Purpose |
|---|---|
| `align.sh` | ALIGN: bowtie2 paired-end alignment of one chunk (`--chunk-skip`, `--chunk-upto`) piped into `samtools view -b`; never SAM on disk |
| `align_merge.sh` | ALIGN_MERGE: `samtools cat` of the chunk BAMs in chunk order; the chunk logs are summed by `merge_bowtie2_logs.py` |
| `alignment_summary.py` | ALIGNMENT_SUMMARY: alignment rates from the bowtie2 logs; warns below `alignment_min_pct` |
| `bam_filter.sh` | BAM_FILTER: sort, raw stats, MAPQ and flag filter (`-q`, `-f 3 -F 2828`), `remove_chroms` |
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
| `merge_bowtie2_logs.py` | Helper of ALIGN_MERGE: sums the paired-end bowtie2 summaries of the chunk logs, in bowtie2's own format. Its own key (`merge_bowtie2_logs`), so a project can replace it |
| `motif_common.R` | MOTIF_COMMON: motifs significant in several comparisons, per region set (find_common_motifs.Rmd) |
| `motif_enrichment.R` | MOTIF_ENRICHMENT: monaLisa binned motif enrichment, diff mode or a gene or regions target (motif_analysis.Rmd) |
| `multiqc.sh` | MULTIQC: one report over the QC files |
| `my_theme.R` | Plot theme for the R figures (used by QC_FRAGMENT_SIZES_PLOT) |
| `my_theme_py.py` | Python plot theme; no step uses it at present |
| `norm_tmm.R` | NORM_TMM: csaw and edgeR TMM factors for one group, MA plots (statistician's script) |
| `peak_merge.sh` | PEAK_MERGE: sort-bed then `bedops --merge` of the listed peak files |
| `qc_chrom_fraction.py` | QC_CHROM_FRACTION: percent of mapped reads on each removed contig |
| `qc_fragment_sizes_plot.R` | QC_FRAGMENT_SIZES_PLOT: density and boxplot SVGs from the two fragment tables |
| `qc_fragment_sizes_table.R` | QC_FRAGMENT_SIZES_TABLE: fragment-size distribution and summary per BAM (Rsamtools) |
| `qc_frip.sh` | FRIP: fraction of reads in peaks for one sample against its own peak file |
| `qc_frip_table.py` | QC_FRIP_TABLE: merges the per-sample FRiP files |
| `qc_peak_stats.py` | QC_PEAK_STATS: peak counts, widths and chromosome summaries for every peak file |
| `qc_tss_enrichment.sh` | TSS_ENRICHMENT: computeMatrix and plotProfile around the TSS for one bigWig |
| `qc_tss_table.py` | QC_TSS_TABLE: merges the per-sample TSS files with the samplesheet metadata |
| `tn5_shift.sh` | TN5_SHIFT: alignmentSieve --ATACshift, then sort and index (pysam) |
| `trim.sh` | TRIM: Trimmomatic paired-end trimming with the ATAC settings |
| `txdb_build.R` | TXDB_BUILD: TxDb (gene and transcript models) from the GTF, saved as SQLite |

Argument keys: `<script>_args` (for example `fastqc_args`, `align_args`, `count_pairs_args`, `align_merge_args`, `bam_filter_args`,
`dedup_args`, `blacklist_filter_args`, `bam_finalize_args`, `tn5_shift_args`, `contamination_args`, `trim_args`, `qc_chrom_fraction_args`, `qc_fragment_sizes_table_args`,
`qc_fragment_sizes_plot_args`, `norm_tmm_args`, `bigwig_args`, `qc_tss_enrichment_args`, `qc_tss_table_args`,
`callpeaks_macs3_args`, `peak_merge_args`, `qc_frip_args`, `qc_frip_table_args`, `qc_peak_stats_args`, `alignment_summary_args`,
`multiqc_args`, `bam_prealigned_args`, `txdb_build_args`, `diff_accessibility_args`, `diff_annotate_args`, `motif_enrichment_args`,
`motif_common_args`). `bigwig_args` is also a declared param. `macs3_args` is a declared param that the MACS3 module puts on the command line
directly. `callpeaks_macs3_args` is the script-level key. Both are passed.

Plot-only changes (colours, labels, sizes): replace the plot script and/or `my_theme.R`. Changes to what is computed: replace the
compute script and keep its output names. Changes to normalisation factors, groups, thresholds or models are scientific decisions:
confirm them first.

## Downstream analyses

Three switches add the differential, annotation and motif steps after the peaks, before MultiQC. All are off by default, and then
nothing of them is run or checked. Each comparison is a case side against a control side of a samplesheet column (`group_column`,
default `condition`), with its own peak groups.

| Switch | Runs | One task per |
|---|---|---|
| `run_differential` | DIFF_ACCESS (csaw, edgeR), and TXDB_BUILD when promoters are computed | comparison that computes a region set |
| `run_annotation` | DIFF_ANNOTATE (ChIPseeker), and TXDB_BUILD | comparison, on its peaks table |
| `run_motifs` | MOTIF_ENRICHMENT (monaLisa), MOTIF_COMMON | comparison and `motif_regions` set; `motif_targets` row; region set with two or more non-pooled comparisons |

The peaks table is always computed by DIFF_ACCESS unless `diff_tables` gives it. `diff_regions` lists the region sets DIFF_ACCESS
computes (`peaks`, or `peaks,promoters`). A comparison that has a table in `diff_tables` is not recomputed for that region set.
Use `diff_tables` to start the annotation or the motifs from tables that already exist. Use `peak_files` to add existing peak files as
peak groups.

Pre-aligned libraries: a samplesheet row with a final BAM in the `bam` column (and no FASTQ) joins the analysis like a FASTQ library
(BAM_PREALIGNED). Its BAM is assumed final. Its row needs a `tmm_group` or `peak_group` to take part in those steps.

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
| `case_peaks`, `control_peaks` | Comma lists of peak group names: MACS3 peak groups (needs `run_macs3`) or names in `peak_files`. Default: the distinct `peak_group` values of the side |
| `pooled` | `true` for a comparison that pools several lines. It is analysed like the others and left out of MOTIF_COMMON |
| `fixed_dispersion` | Optional number: the LRT dispersion of this comparison (default `diff_fixed_dispersion`) |

No sample may be on both sides. Each side needs at least one sample when the comparison is computed.

### peak_files.csv and diff_tables.csv

peak_files.csv:

```
name,path
extra,/abs/path/extra.narrowPeak
```

diff_tables.csv:

```
name,regions,path
BT54_high_vs_low,promoters,/abs/path/BT54_high_vs_low_promoters.csv
```

`peak_files`: names unique and not a `peak_group` of the samplesheet; the file exists and ends in `.narrowPeak`, `.broadPeak` or `.bed`.
`diff_tables`: the name is a comparison; `regions` is `peaks` or `promoters`; the file exists; each (name, regions) pair is listed once.

### motif_targets.csv

```
name,type,target,comparison,label,min_peaks
ST8SIA1_gene,gene,ST8SIA1,gd2_high_vs_low,ST8SIA1 peaks,
ST8SIA1_cCRE,regions,/abs/path/st8sia1_ccres_hs1.bed,gd2_high_vs_low,ST8SIA1 cCRE,
```

`type` is `gene` (the peaks of the `all_sites` annotation table whose `gene_column` contains the symbol) or `regions` (peaks that overlap
the BED against the other peaks). `comparison` is the universe: its `all_sites` table. It needs `run_annotation`. The `name` must differ
from every comparison name, because both write under `results/motif/`. A regions BED must already be in the pipeline genome's
coordinates: the legacy liftOver is a project pre-step. `label` must not contain quotes, backticks, backslashes or dollar signs.

### Startup checks

Downstream values are checked at startup, whenever their file or switch is set, with a message per problem: an unknown or duplicate name, a
column that is not in the samplesheet, a value of `group_column` that no active sample has, a sample on both sides, a missing peak group or
peak file, a missing table or BED, a `diff_regions` other than `peaks` or `peaks,promoters`, a `motif_regions` set without its table
(promoters and peaks need the table of the comparison; the annotation sets need `run_annotation`), a missing `gtf` when promoters or
annotation need the TxDb, and a wrong number, boolean or enumerated value.

### Outputs

```
results/annotation/txdb/                 txdb.sqlite (TXDB_BUILD); log/
results/differential/                    <name>.csv (peaks), <name>_promoters.csv, <name>.norm_factors.tsv, <name>.summary.tsv,
                                         <name>.replicates.csv, <name>.diag_*.svg, <name>.diagnostics.pdf (plot_format pdf); log/
results/differential/bed/                <name>.bed, <name>_upregulated.bed, <name>_downregulated.bed, <name>_promoters*.bed
results/annotation/                      <name>_all_sites.csv, <name>_promoter_sites.csv, <name>_non-promoter_sites.csv,
                                         <name>_distal_sites.csv and their .bed files, <name>_peak_distribution.svg,
                                         <name>.annotation_summary.tsv; log/
results/motif/<name>/<region>/           <name>_motif_analysis.csv, <name>.motif_summary.csv, heatmaps, similarity, GC and
                                         dinucleotide plots (diff mode); log/
results/motif/<name>/                    the same for a motif_targets row (gene or regions); <name>.skipped.txt when skipped
results/motif/common/<region>/           common_significant_motifs.csv, venn_common_motifs.svg, motif_analysis_summary.csv; log/
```

Every step writes its command to `log/<tag>_<step>_command.sh`. A skipped motif analysis writes only `<name>.skipped.txt`, with the
reason and the counts.

## Examples

- `examples/config.yaml`: the generic template, copied by `init`. Every line marked `# set me` must be filled in (`genome_index`,
  `blacklist`, `gtf`, `tss_bed`, `effective_genome_size`, `chrom_sizes`, `contamination_index`, `project_name`).
- `examples/samplesheet.csv` and `examples/peak_merges.csv`: placeholders with `/path/to/...` FASTQ paths and a one-row merge.
- `examples/projects/av_atacseq_2026/`: AV_ATACSEQ_2026, 39 libraries (30 included), T2T-CHM13v2.0, nine TMM groups, the peak groups of
  the legacy replicate patterns, and `peak_merges.csv` with the high and low sets. Its README lists the toggles.
- `examples/projects/am_atacseq/`: AM_ATACSEQ, nine mouse (mm39) libraries, one TMM group (`WT_RAG2_Baseline`), peak groups WT, RAG2
  and Baseline. Its config lists the AM legacy values as `# TOGGLE` lines.
- `examples/comparisons.csv` and `examples/motif_targets.csv`: one comparison and one target, matching `examples/samplesheet.csv`.
  `examples/config.yaml` has the downstream keys commented, with the defaults and the legacy values.
- `examples/projects/av_atacseq_2026/`: also `comparisons.csv` (the five per-line comparisons and gd2) and `motif_targets.csv`
  (ST8SIA1). The AM project has `comparisons.csv` (WT and RAG2 against the controls); its downstream switches are off.

Copy an example into a new directory and adapt it. Do not run an example in place: the paths point at the legacy project folders.

## Tests

- Startup validation runs at the start of every `nextflow run`, including `-preview` and `-stub-run`.

## Containers

Each module has one literal `container '...'` string. The lab images are in `/project/rephilllab/nf_containers` (Apptainer cache,
shared); the `pmacs` profile binds `/usr/lib64/ossl-modules` (FIPS).

| Used by | Image |
|---|---|
| FASTQC_RAW, FASTQC_TRIM | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` |
| TRIM | `quay.io/biocontainers/trimmomatic:0.40--hdfd78af_0` |
| COUNT_PAIRS, ALIGN, ALIGN_MERGE, ALIGNMENT_SUMMARY, BLACKLIST_FILTER, CONTAMINATION, QC_CHROM_FRACTION, QC_TSS_TABLE, QC_PEAK_STATS, QC_FRIP_TABLE | `/project/rephilllab/nf_containers/cutandrun-py_1.2.sif` (python, bowtie2 2.5.5, samtools 1.23, bedtools 2.31.1) |
| BAM_FILTER, BAM_FINALIZE, BAM_PREALIGNED | `quay.io/biocontainers/samtools:1.23--h96c455f_0` |
| DEDUP | `quay.io/biocontainers/sambamba:1.0.1--he614052_4` |
| FRIP | `quay.io/biocontainers/bedtools:2.31.1--h13024bc_3` |
| BIGWIG, TSS_ENRICHMENT, TN5_SHIFT | `quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0` |
| MACS3 | `quay.io/biocontainers/macs3:3.0.4--py312h71493bf_0` |
| PEAK_MERGE | `quay.io/biocontainers/bedops:2.4.42--hd6d6fdc_1` |
| MULTIQC | `quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1` |
| NORM_TMM, QC_FRAGMENT_SIZES_TABLE, QC_FRAGMENT_SIZES_PLOT | `/project/rephilllab/nf_containers/cutandrun-r_1.3.sif` (R 4.3.3, csaw 1.36.0, edgeR 4.0.16) |
| TXDB_BUILD, DIFF_ACCESS, DIFF_ANNOTATE, MOTIF_ENRICHMENT, MOTIF_COMMON | `/project/rephilllab/nf_containers/omics-r_1.0.sif` (downstream steps; csaw, edgeR, ChIPseeker, monaLisa) |

The sambamba and bedops images are not in `nf_containers` yet; Nextflow pulls them on first use. The `.sif` files are the lab's
builds shared with the CUT&RUN pipeline.

## Status and limits

- Version 0.1.0.
- Paired-end only.
- The legacy results are not reproduced byte for byte: the port fixes known legacy bugs.
- Not yet run end to end on real data.
