# CUTANDRUN: CUT&RUN / CUT&Tag pipeline (Rephill Lab)

A Nextflow pipeline for paired-end CUT&RUN and CUT&Tag data on a cluster with
LSF and Apptainer (the PMACS cluster at UPenn). It runs one reusable pipeline with one
`config.yaml` and three small CSVs per project.

- Rules for agents working in a project: `agent/AGENTS.md`

## What it does (stages 1a and 2)

```
FASTQ pairs ──> FASTQC (raw)
            └─> TRIM (trimmomatic) ─> ALIGN (bowtie2 | samtools view -b, combined target + E. coli index)
                  ─> BAM_FILTER (MAPQ + flag filter, raw stats)
                  ─> MARKDUP (Picard)
                  ─> DEDUP_SPLIT (spike-in BAM split off; chrM and spike contigs dropped; dedup)
                  ─> BLACKLIST_FILTER (pair-aware) ─> BAM_FINALIZE (final BAM, final stats, orphan check)
pre-aligned BAM rows ─> BAM_PREALIGNED ───────────────────────────────┘   (same final-BAM channel)

final BAMs ─> QC: alignment stats, fragment sizes, library complexity (NRF/PBC1/PBC2, Picard ELC),
              preseq, fingerprint, correlation, FRiP, TSS enrichment, peak stats ─> MULTIQC + HTML REPORT
final BAMs ─> BIN_COUNTS (csaw bin counts, once) ─> REPLICATE_CORRELATION ─> REPLICATE_CORRELATION_PLOT
per norm_run ─> CHIPSEQSPIKEINFREE | TMM | NORM_CPM | SPIKEIN (norm_spikein) ─> NORM_COMBINE ─> <norm_run>.scale_factors.tsv (primary) + .scale_factors_all.tsv
scale_factors_all + BIN_COUNTS ─> NORM_QC ─> NORM_QC_PLOT
scale factors + final BAMs ─> BIGWIG (bamCoverage) [─> BIGWIG_LOG2RATIO]
final BAMs + peaksets.csv ─> CSAW (default) [─> CSAW_CUTOFF_SCAN]
final BAMs + classifications.csv + scale factors ─> CLASSIFY_REGIONS ─> CLASSIFY_REGIONS_PLOT   (optional: run_classify)
final BAMs + controls ─> MACS3 (run_macs3 null = on when every target has a control);  SEACR_BEDGRAPH -> SEACR_CALL -> SEACR_CONSENSUS (run_seacr, same rule)
stage 2 (optional, see "Heatmaps, tracks, annotation, intersections"):
  heatmaps.csv ─> COMPUTEMATRIX | BUILD_MATRIX ─> HEATMAP_CLUSTER (seeded k-means, optional) ─> PLOTHEATMAP | PLOT_ENRICHED
  heatmaps.csv with `auto` z_max/y_max ─> HEATMAP_LIMITS ─> PLOTHEATMAP (deepTools engine)
  tracks.csv ─> TRACKS_INI ─> TRACKS_LIMITS (shared / group:NAME limits) ─> PYGENOMETRACKS;  peak sets ─> PEAK_ANNOTATION;  intersects.csv ─> INTERSECT
```

Key points:

- Target AND spike-in reads are deduplicated by default (IgG always is).
- The primary normalisation is chosen per `norm_run`: ChIPseqSpikeInFree (NUC, EED), E. coli
  spike-in (CBX2), TMM or CPM. Every method writes the same scale-factor table, and
  downstream steps read only the primary's `bigwig_multiplier` and `edger_norm_factor` columns.
  Every run also computes the other methods (`norm_methods`) for the MA diagnostic.
- csaw, MACS3 and SEACR all run by default. csaw thresholds come from `peaksets.csv`; a blank
  cell takes the per-antibody or global default (see "csaw" below). MACS3 and SEACR need a
  `control_id` per target; without one they are skipped with a warning unless `run_macs3` /
  `run_seacr` were set to `true` explicitly. Broad or narrow mode follows the antibody (`broad_marks`).
- Samples in `skip_groups` stay in the samplesheet but are skipped by all core steps.

**Stage 2 (optional, all off by default):** deepTools and EnrichedHeatmap
heatmaps (`run_heatmaps`), pyGenomeTracks tracks (`run_tracks`), peak annotation
(`run_annotation`) and BED intersections (`run_intersect`); see "Heatmaps, tracks,
annotation, intersections" below.

Single-end data is not supported (startup error).

## Quick start

One-time setup: put `/project/rephilllab/bin` on your `PATH` (it holds one
symlink per lab pipeline: `run_cutandrun.sh -> ../CUTANDRUN/run_pipeline.sh`).

```bash
mkdir my_project && cd my_project
run_cutandrun.sh init          # copies config.yaml, samplesheet.csv, norm_runs.csv, peaksets.csv, classifications.csv, heatmaps.csv, tracks.csv, intersects.csv, examples/, AGENTS.md, CLAUDE.md; git init
# edit config.yaml (project_name, experiment, genome_index, chrom_sizes, blacklist, gtf, tss_bed ...)
# edit samplesheet.csv, norm_runs.csv, peaksets.csv
run_cutandrun.sh               # submits one LSF head job; -resume is built in
bjobs -u $USER
tail -f metadata/nf_logs/nextflow_<timestamp>.out
```

`examples/projects/{nucleation,cbx2,eed}/` hold complete input sets that
re-create each manuscript analysis; copy one into a new project and adapt it. The dataset configs' toggles (the alternatives to the chosen values) are comments next to the keys in `examples/projects/<p>/config.yaml`.
The head job runs in `/project/rephilllab/nf_work/$USER/<project dir name>/`.
A watchdog (`lib/lsf_watchdog.sh`) in the head job stops the run with a `CUTANDRUN WATCHDOG:` message if an LSF job is aborted at dispatch (Nextflow would otherwise hang); relaunch with `run_cutandrun.sh`. `lsf_watchdog: false` in config.yaml turns it off.
Nextflow 26.04.6 is pinned and downloaded once into `.nextflow-bin/`.

## Input files

All live at the project root.

### samplesheet.csv (one row per library)

`sample_id,fastq_1,fastq_2,bam,antibody,condition,replicate,control_id,group,norm_run,exclude`

| Column | Meaning |
|---|---|
| `sample_id` | Unique; used verbatim in file names; `^[A-Za-z0-9][A-Za-z0-9_.-]*$` |
| `fastq_1`, `fastq_2` | Absolute R1/R2 paths. Give these OR `bam` |
| `bam` | Pre-aligned, final-filtered BAM; skips trim/align/filter, no spike-in BAM (not usable in a `spikein` norm_run) |
| `antibody` | Mark/antibody; `IgG` (any case) = control |
| `condition`, `replicate` | Labels (e.g. `DMSO`, `A`) |
| `control_id` | `sample_id` of the IgG control; empty for IgG rows. Required for MACS3, SEACR, log2ratio |
| `group` | Unit of replication; empty = `<condition>_<antibody>` |
| `norm_run` | Row of `norm_runs.csv`; required for non-skipped rows unless `normalisation: cpm` |
| `exclude` | `true` drops the row from the run |

### norm_runs.csv (one row per normalisation run)

`norm_run,method,bin_size,cutoff,max_turns` (+ optional `min_first_turn`, `ref_condition`) with `method` = `chipseqspikeinfree`
(all three numbers required), `spikein` (numbers blank), `tmm` (`bin_size`
optional, default `tmm_bin_size`) or `cpm` (library size only; numbers blank).
`method` is the run's **primary** method: only it scales the bigWigs and the classification. The other methods of
`norm_methods` (and spike-in when `norm_spikein: true`) are computed for comparison in `scale_factors/<norm_run>.scale_factors_all.tsv`
and the MA plots. A run with one library needs `cpm` (startup error otherwise).

`ref_condition` (optional) is the reference condition of the run's MA plots (negative M = less signal than the reference). Blank = the condition
of the `ref_group` of the classifications that use the run (if exactly one), else the condition of the run's first library. It must be a
condition of the run. The numbers of `bin_size`, `cutoff` and `max_turns` of a non-ChIPseqSpikeInFree run are ignored; when ChIPseqSpikeInFree is a
diagnostic and the cells are blank it uses `csif_bin_size`, `csif_cutoff` and `csif_max_turns`.

`min_first_turn` (`chipseqspikeinfree` only; numeric; blank = 0.5, the package default) is passed to `norm_chipseqspikeinfree.R` as `--min-first-turn=` and recorded in the scale-factor `note`. A blank cell passes no argument, so the output is unchanged. Other methods ignore it with a warning.

### peaksets.csv (one row per csaw peak set)

`peakset,target_group,control_group,small_filt,large_filt,merge_bp,max_width,fdr,lfc`
(groups are samplesheet `group` values; blank `max_width` = no limit).

Any threshold cell (`small_filt`, `large_filt`, `merge_bp`, `max_width`, `fdr`, `lfc`)
may be blank. Each value resolves in this order:

1. the non-blank cell of the row;
2. `csaw_defaults[<antibody>][<column>]`, where `<antibody>` is the samplesheet
   `antibody` of the row's `target_group` (`csaw_defaults` is a per-antibody map in `config.yaml`);
3. the global `csaw_<column>` (`csaw_small_filt`, ..., `csaw_lfc`).

Startup fails if a `target_group` mixes antibodies, or if a resolved value is not
numeric. `metadata/peaksets.resolved.tsv` lists the resolved value of every threshold
of every row and its source (`row`, `antibody` or `global`).

Startup validation (`prepareRun()` in `main.nf`) fails with a clear message
before any job runs on: missing columns or files, duplicate `sample_id`, unknown
or non-IgG `control_id`, unknown `norm_run`, a peakset group with no samples,
a target group that mixes antibodies, a resolved threshold that is not numeric,
invalid enumerated params.

### csaw: what the FDR and the logFC threshold mean

- `fdr` is the false discovery rate of the test "target enriched over IgG". A
  region with FDR at or below `fdr` is enriched over IgG.
- `lfc` is a descriptive threshold applied after the test. A region is a peak only
  if it passes the FDR and the log2 fold change of its representative window is above
  `lfc`. Which window is set by `csaw_lfc_from`: `combined` (default) takes the
  representative window of the Simes combination (the legacy rule, and the pipeline's
  rule before 2026-10-23); `best` takes the window with the smallest p-value in the
  region (csaw `getBestTest`).
  It is an effect-size filter, not a tested claim: the FDR does not control
  the probability that a peak is more than 2^`lfc`-fold enriched. The `best` window is
  chosen on the same p-values it is then filtered on, so its logFC is the more
  winner's-curse-biased of the two.
- Decision (2026-10-13): there is no `glmTreat` option. The csaw default is
  `csaw_norm_bins: background` (TMM on background bins); `all` uses every bin.

Further csaw keys. Blank cells of `peaksets.csv` take the values below; the other keys
set the csaw calls directly. Defaults keep the current behaviour. The header of
`bin/callpeaks_csaw.R` documents each argument.

| Key | Default | What it does |
|---|---|---|
| `csaw_small_filt`, `csaw_large_filt`, `csaw_merge_bp`, `csaw_max_width`, `csaw_fdr`, `csaw_lfc` | `3`, `1.5`, `50`, `0`, `0.05`, `1` | Global fallbacks for blank `peaksets.csv` cells (`csaw_max_width` 0 = no limit) |
| `csaw_defaults` | `K27me3: {small_filt: 3, large_filt: 1.5, merge_bp: 50}`, `K27me2: {small_filt: 2, large_filt: 1, merge_bp: 50}` | Per-antibody fallbacks for blank cells, keyed by the samplesheet `antibody` of the target group. Setting `csaw_defaults` in config.yaml replaces the whole map |
| `csaw_large_window`, `csaw_large_spacing` | `2000`, `500` | Width and spacing (bp) of the wide window scale of the two-scale design (the small scale is `csaw_window` / `csaw_spacing`) |
| `csaw_bg_max_enrich` | `1` | Background bins (`csaw_norm_bins: background`) have an average abundance of at most the median plus this value (log2 units) |
| `csaw_count_filter` | `10` | Minimum summed count over all libraries for a window to be kept; the background bins keep csaw's default of 10 |
| `csaw_max_frag` | `null` | csaw `max.frag` for CSAW, CSAW_CUTOFF_SCAN, CLASSIFY_REGIONS and BIN_COUNTS (so NORM_QC, TMM and the replicate metrics); `null` = `bowtie2_maxins` (500 if that is null too) |
| `csaw_lfc_from` | `combined` | `combined` or `best`: which window's logFC the `lfc` rule uses (see `lfc` above) |
| `csaw_disjoint_regions` | `false` | `false`: csaw's own merge, regions may overlap (legacy). `true`: merged regions are made disjoint (windows assigned by midpoint) and the statistics are recomputed |

## Region classification (`run_classify`)

Optional (`run_classify: true`, off by default). Each row of `classifications.csv`
counts the reads of a reference group and a treated group in a fixed set of regions
and classes every region as **lost**, **retained** or **unclassified** from the
treated-versus-reference change (edgeR quasi-likelihood, normalised by the scale
factors of the row's `norm_run`).

```
final BAMs + scale factors + regions ─> CLASSIFY_REGIONS ─> CLASSIFY_REGIONS_PLOT (MA plot)
```

### classifications.csv

`classification,regions,ref_group,treat_group,norm_run,lfc_retained,lfc_lost,fdr`

| Column | Meaning |
|---|---|
| `classification` | Id, `^[A-Za-z0-9][A-Za-z0-9_.-]*$`, unique; prefix of every output file and of the named regions |
| `regions` | A csaw peak set name (requires `run_csaw`) or the path of a BED file (BED3 or more; a 4th column gives the name, else `region_<i>`) |
| `ref_group`, `treat_group` | samplesheet `group` values, reference and treated; must differ; each needs at least one active sample |
| `norm_run` | Empty = the `norm_run` shared by all reference and treated samples. Startup error if they differ or if a sample has none |
| `lfc_retained`, `lfc_lost`, `fdr` | Empty = `classify_lfc_retained` (-1), `classify_lfc_lost` (-2), `classify_fdr` (0.05) |

Rules, with `logFC` (treated over reference) and `FDR` (Benjamini-Hochberg over all
regions of the row):

- **lost**: `logFC <= lfc_lost` and `FDR < fdr`;
- **retained**: `round(logFC, 1) >= lfc_retained`, or `FDR >= fdr`. A region without a
  significant change counts as retained: "no significant loss", not "shown to be retained";
- **unclassified**: the rest. A region that is both lost and retained is lost.

Counts use the samples' BAMs with the csaw read settings (`mapq`, `remove_duplicates`,
`csaw_max_frag` or `bowtie2_maxins`, paired-end). The normalisation factors are the
`edger_norm_factor` column of `results/scale_factors/<norm_run>.scale_factors.tsv`; they
are not re-estimated from the regions. A group with one library gives a warning: the
dispersion then comes from the fixed-dispersion fallback of csaw peak calling (one library
per group) or from the replicated group only, and the FDR is not valid. Reads in overlapping
input regions are counted in each region; a warning gives the number of overlapping regions.

Outputs, in `results/classify/<classification>/`:

| File | Content |
|---|---|
| `<id>.retained.bed`, `<id>.lost.bed`, `<id>.unclassified.bed` | BED6, 0-based starts, name, score = logFC. An empty class gives an empty file |
| `<id>.per_region.tsv.gz` | One row per input region: chr, start0, end, name, logFC, logCPM, PValue, FDR, class |
| `<id>.summary.tsv` | Per class: n and Mb; thresholds; n of reference and treated libraries; their ids; `norm_run`; dispersion |
| `<id>.ma.svg` | MA plot (logCPM against logFC) coloured by class, with the threshold lines |

The BED starts are 0-based. The legacy classifier wrote 1-based starts, so its BEDs
were shifted by 1 bp.

**Named regions.** `<id>.retained`, `<id>.lost` and `<id>.unclassified` are accepted
wherever a region or BED name is: `heatmaps.csv` `regions`, `tracks.csv` `beds`,
`intersects.csv` `bed_a` and `bed_b`. A class name that equals a csaw peak set name is
a startup error. An unknown name stops the run with a message that lists the valid
peak set and class names.

## Heatmaps, tracks, annotation, intersections (stage 2)

Each CSV is read only when its switch is on: `heatmaps` with `run_heatmaps`,
`tracks` with `run_tracks`, `intersects` with `run_intersect`; annotation is
switched by `run_annotation`. The CSVs live at the project root and are
validated at startup. Template rows are in `examples/heatmaps.csv`,
`examples/tracks.csv` and `examples/intersects.csv`.

Quoted lists use commas inside double quotes (`"IgG_WT_rep1,K27me3_WT_rep1"`).
`regions` and `beds` take a csaw peak set name (its `<peakset>.peaks.bed`), a class
name of a classification (`<id>.retained`, `<id>.lost`, `<id>.unclassified`; see
"Region classification"), an intersect name (its `<name>.intersect.bed`, see "intersects.csv") or the path of an existing BED file. INTERSECT runs before the heatmaps, tracks and annotation, so they can name an intersect. An intersect name that equals a csaw peak set or class name is a startup error. Intersect rows themselves take peak sets, classes and paths only.

### heatmaps.csv (`run_heatmaps: true`)

A sample that is not in the samplesheet (or is excluded / in `skip_groups`) is dropped from a row with its `labels`, `colors`, `z_min` and other per-sample entries, with one WARN per row (a row left with no sample is skipped); `tracks.csv` works the same way.

`heatmap,regions,samples,labels,window_bp,bin_bp,sort_using,cluster_k,z_min,z_max,y_max,colors,engine`,
then the optional columns `cluster_samples,cluster_labels,legend_title,region_labels`.

| Column | Meaning |
|---|---|
| `heatmap` | Id of the heatmap; used in every output name |
| `regions` | csaw peak set name, class name (`<id>.retained`, ...), intersect name or BED path. A quoted list gives several region groups: one per entry, passed together to computeMatrix in list order. `enrichedheatmap` and any row with `cluster_k` accept one region only (startup error) |
| `samples` | Quoted list of `sample_id`s, one column each |
| `labels` | Quoted list of display labels, one per sample; a literal `\n` makes a line break |
| `window_bp`, `bin_bp` | bp on each side of the region centre, and the bin size; `window_bp` must be a multiple of `bin_bp` |
| `sort_using` | Empty, or a `sample_id`: rows are sorted by that sample's signal |
| `cluster_k` | Empty or `0`: no clustering. An integer >= 2: that many clusters. A range such as `2-6`: the gap statistic chooses k |
| `z_min`, `z_max`, `y_max` | Empty, one number, or one per sample. `z` sets the colour scale; `y_max` sets the profile height and is independent of `z_max`. `z_max` and `y_max` may be `auto` (deepTools engine only; startup error for `enrichedheatmap`), see below |
| `colors` | Quoted list, one hue per sample (`"#8c96c6"`); the colour ramp runs from white to the hue |
| `engine` | `deeptools` (computeMatrix + plotHeatmap) or `enrichedheatmap` (R EnrichedHeatmap) |
| `cluster_samples` | Optional. Quoted subset of `samples` used for clustering; empty = all |
| `cluster_labels` | Optional. Quoted list of k labels; if the count is not k, the labels are `cluster 1` .. `cluster k` |
| `legend_title` | Optional. EnrichedHeatmap legend title per sample; empty = the antibody of that sample in the samplesheet |
| `region_labels` | Optional. Quoted list, one label per region of `regions` (same count); empty = the entries of `regions` |

`auto` limits (deepTools engine). HEATMAP_LIMITS runs on the plotted matrix. For each
sample it takes the mean profile of each region group (per bin column, NaN ignored; the
lines plotHeatmap draws) and the maximum over those profiles, `profile_max`, so no group's
profile runs off the plot. The limit of an antibody is
`round_up_2sf(1.05 x` the largest `profile_max` of that antibody's samples`)`, rounded up
to two significant figures; a sample without a finite profile gets 1. Each sample uses its
limit as `z_max` and `y_max`; `z_min` is 0 unless the row sets it. The antibody comes from
the samplesheet, so the samples of one antibody share one scale. The table is
`results/heatmaps/limits/<heatmap>_limits.tsv` (`sample`, `antibody`, `profile_max`, `limit`).

Heatmap parameters (`heatmap_*`, see the params table) set the computeMatrix and plotHeatmap
options of the deepTools engine; EnrichedHeatmap plots are unchanged.

Clustering is a separate, seeded step (`heatmap_cluster.py`: scikit-learn
k-means++, `n_init` 100, `random_state` = `heatmap_seed`). For a k range, the
chosen k is the smallest k with gap(k) >= gap(k+1) - se(k+1), otherwise the k with
the largest gap (20 uniform reference datasets). Clusters are relabelled 1..k by
decreasing mean signal. The k table (`cluster_k.tsv`) lists silhouette and gap for
each k. The cluster BEDs are written under `results/heatmaps/cluster/`, and the
clustered deepTools plot recomputes its matrix with those BEDs as region groups.
deepTools never clusters by itself (no `--kmeans`).

### tracks.csv (`run_tracks: true`)

`track_set,locus_or_gene,flank,samples,labels,colors,max_value,beds,scalebar`,
then the optional columns `log2ratio`, `bed_colors`, `bed_labels`, `ini_args`.

| Column | Meaning |
|---|---|
| `track_set` | Id of the track set; one figure per locus, named `<track_set>.<locus>.svg` |
| `locus_or_gene` | Quoted list of gene names (exact `gene_name` in the GTF) and/or regions `chr:start-end`; an unknown gene is an error |
| `flank` | bp added on each side of a gene; one value or one per locus; ignored for regions |
| `samples`, `labels`, `colors` | As in heatmaps: sample ids, display labels, one hue per sample |
| `max_value` | One entry per sample (comma list; blank entries allowed): a number (fixed), `auto` (pyGenomeTracks autoscale), `group:NAME` (the shared limit of the tracks with that name at this locus), or blank = `tracks_max_value_default` (`shared`: one group per samplesheet antibody) |
| `beds` | Quoted list of csaw peak set names, class names (`<id>.retained`, ...), intersect names or BED paths, drawn as BED tracks |
| `scalebar` | bp, or `auto` (pyGenomeTracks picks the size); empty = no scale bar |
| `log2ratio` | Optional, `true` or `false`. Adds a log2-ratio track after each sample's bigWig, read from `S.vs_<control_id>.log2ratio.bw`; these exist only with `run_log2ratio: true` |
| `bed_colors` | Optional. Quoted list, one colour per entry of `beds` (used for the fill and the border); empty = `track_bed_color` |
| `bed_labels` | Optional. Quoted list, one title per entry of `beds`; empty = the file name without extension |
| `ini_args` | Optional. Extra options for this track set's make_tracks_ini.py, written as in the command line and starting with `--`, e.g. `--genes-height=2 --merge-transcripts false --gene-rows 2 --arrow-interval 6 --bed-after SAMPLE_ID --gap-after SAMPLE_ID=0.5`. Other options: `--genes-style UCSC\|flybase\|tssarrow`, `--arrowhead-included true\|false`; `--gap-after` takes comma-separated `ID=height` entries (spacer height after that sample, default 0.25); `--bed-after` draws the BED tracks after that sample instead of after the last. No spaces or quotes inside a value |

Shared limits. A `shared` or `group:NAME` track takes its maximum from its group at its
locus: `round_up_2sf(1.05 x` the largest per-bin mean over the group's bigWigs at the plotted
locus, with `track_number_of_bins` bins`)`. TRACKS_LIMITS computes this per locus and writes
`results/tracks/limits/<slug>.resolved.ini` (the ini with every `group:` token replaced by its
number) and `results/tracks/limits/<slug>_limits.tsv` (`group`, `track`, `file`, `stat`, `track_peak`,
`limit`). PYGENOMETRACKS renders only the resolved ini. The default `shared` makes the
conditions and replicates of one mark share one y-range; the IgG forms its own group. Set
`tracks_max_value_default: auto` to get the previous per-track autoscale. log2-ratio tracks
keep their fixed `track_log2_min` and `track_log2_max`; they are never grouped.

`tracks_limit_stat` chooses the statistic of the shared limit. `mean_bins` (default) is the rule above: the largest per-bin mean over `track_number_of_bins` bins. `max_bin` is the CBX2 rule of 2026-10-15: the largest value of the bigWig's own bins over the whole locus window. The 1.05 factor and the round-up to 2 significant figures are the same in both. The `stat` column of `<slug>_limits.tsv` records the mode and `track_peak` the track's value.

The .ini file uses relative file names only. Locus slugs are `<gene>_<flank/1000>kb`
(or `<gene>_<flank>bp` when the flank is not a multiple of 1000) and `chr_s_e` for regions.

### intersects.csv (`run_intersect: true`)

`name,bed_a,bed_b,label_a,label_b,min_overlap,reciprocal[,output]`

- `bed_a`, `bed_b`: peak set name, class name (`<id>.retained`, ...) or BED path (intersect names are not accepted here). Records of `bed_a` that overlap `bed_b`
  are written once each (`bedtools intersect -wa -u`).
- `min_overlap`: fraction 0-1 of `bed_a` that must overlap (`-f`); empty = the script
  default (0.05).
- `reciprocal`: `true` also requires the same fraction of `bed_b` (`-r`).
- `output`: `records` (default, also when empty) writes each `bed_a` record that overlaps `bed_b` once; `overlap` writes the overlapping pieces, one line per `bed_a`/`bed_b` pair (`bedtools intersect` without `-wa -u`), as the legacy step_7b script and the GG_CBX2 paper heatmaps used. Counts in the summary stay record-based; `n_pieces` counts the lines of `<name>.intersect.bed`.

### Peak annotation (`run_annotation: true`)

ChIPseeker annotation against the `gtf`, with the TxDb built inside each run. It
annotates the csaw peak sets (label `<peakset>`) and, when those steps are on, the
MACS3 sets (`macs3_<group>`) and the SEACR consensus sets (`seacr_<group>`). The
labels are the keys of `annotation_colors`:

- `annotation_colors`: `label=#RRGGBB,label2=#RRGGBB`; labels without a colour get palette colours.
- `annotation_kegg`: `true` draws a KEGG dotplot; `false` writes a placeholder page.
- `annotation_tss_region`: bp on each side of the TSS (default 3000).
- `annotation_extra_beds`: comma list of extra regions annotated next to the peak sets: csaw peak set names, class names (`<id>.retained`, ...), intersect names or BED paths. The label is the name (the file stem for a BED path). A csaw peak set that `run_csaw` already annotates is skipped with a warning. Labels must be unique. Default empty.

### Results layout (stage 2)

```
results/heatmaps/{matrix,cluster,plots,limits}/   matrices; cluster assignments, k table, cluster BEDs; heatmap SVGs; auto limits
results/tracks/{ini,plots,limits}/         generated .ini and region tables; figure SVGs; resolved .ini and limit tables
results/annotation/                        annotation tables, plots, KEGG page
results/intersect/                         <name>.intersect.bed, <name>.only_a.bed, <name>.intersect_summary.tsv
results/norm_qc/                           also <norm_run>.ma.svg (NORM_QC_PLOT)
```

## Configuration (`config.yaml`, flat `key: value`)

Every key has a default in `nextflow.config`; a project lists only what it
changes. Reference paths (`genome_index`, `chrom_sizes`, `blacklist`, `gtf`,
`tss_bed`) have no default.

| Key | Default | Meaning |
|---|---|---|
| `project_name` | `my_cutandrun` | Prefix of the QC report and project-level outputs |
| `samplesheet` | `samplesheet.csv` | Library table |
| `norm_runs` | `norm_runs.csv` | Normalisation runs table |
| `peaksets` | `peaksets.csv` | csaw peak sets table |
| `fastq_dir` | `null` | Optional FASTQ folder; `run_cutandrun.sh` moves its FASTQs to `original_data/fastq/` and rewrites this key. Samplesheet paths are used as given either way |
| `outdir` | `./results` | Final outputs |
| `metadata_dir` | `./metadata` | Nextflow report/trace/timeline/DAG |
| `intermediate_dir` | `./intermediate_files` | Trimmed FASTQ, aligned/filtered/markdup BAMs, SEACR bedGraphs |
| `publish_mode` | `copy` | `copy` or `symlink` |
| `experiment` | `cutandrun` | `cutandrun` (TruSeq3-PE-2) or `cutandtag` (NexteraPE-PE) adapters |
| `genome_index` | `null` | bowtie2 index PREFIX of the combined target + spike-in genome |
| `chrom_sizes` | `null` | chrom.sizes of the combined genome |
| `blacklist` | `null` | Blacklist BED (pair-aware removal) |
| `gtf` | `null` | Gene annotation GTF; TSS enrichment anchors on its transcript start sites |
| `tss_bed` | `null` | TSS BED for TSS enrichment; fallback only, used (and required) when `gtf` is unset |
| `spike_prefix` | `ecoli_chr` | Contig-name prefix of the spike-in genome |
| `mito_chrom` | `chrM` | Mitochondrial contig, dropped from the target BAM |
| `effective_genome_size` | `3100000000` | MACS3 `-g` / bamCoverage (T2T-CHM13) |
| `paired_end` | `true` | Must be true (paired-end only) |
| `trim` | `true` | Run trimmomatic before alignment |
| `trim_minlen` | `25` | Trimmomatic MINLEN |
| `bowtie2_mode` | `local` | `local` (`--local --very-sensitive-local`) or `end_to_end` (`--end-to-end --very-sensitive`) |
| `bowtie2_minins` | `null` | optional bowtie2 `-I`; `null` = omitted (bowtie2 default 0, as in the manuscript); 10 was used in 2026-10 test runs |
| `bowtie2_maxins` | `null` | optional bowtie2 `-X`; `null` = omitted (bowtie2 default 500, as in the manuscript); 700 was used in 2026-10 test runs. Also csaw `max.frag` unless `csaw_max_frag` is set (500 if both are null) |
| `mapq` | `30` | Minimum MAPQ kept |
| `remove_duplicates` | `true` | Remove duplicates from target BAMs (IgG always) |
| `spikein_dedup` | `true` | Remove duplicates from the spike-in BAM |
| `normalisation` | `chipseqspikeinfree` | `chipseqspikeinfree`, `spikein`, `tmm` or `cpm`. Informational except `cpm` (no scale-factor tables, multiplier 1); the samplesheet `norm_run` decides per sample |
| `tmm_bin_size` | `10000` | Bin size for TMM; also the bins of NORM_QC and (when `repcor_bin_size` is null) the replicate metrics. BIN_COUNTS counts every library once per bin size |
| `norm_methods` | `'chipseqspikeinfree,tmm,cpm'` | Methods computed for every `norm_run` besides its primary (`method` of `norm_runs.csv`); `''` = primary only. Diagnostics only feed `scale_factors_all.tsv`, NORM_QC and NORM_QC_PLOT |
| `norm_spikein` | `false` | Also compute physical spike-in factors as a diagnostic (a run with pre-aligned BAMs is `unavailable`, with a warning) |
| `csif_bin_size`, `csif_cutoff`, `csif_max_turns` | `10000`, `1.2`, `0.99` | ChIPseqSpikeInFree settings when it is a diagnostic and the `norm_runs.csv` cells are blank |
| `norm_qc_min_count` | `20` | MA plots: minimum fragments summed over the two libraries of a pair |
| `norm_qc_prior` | `0.5` | MA plots: prior count of the log2 CPM (the reference lines do not depend on it) |
| `norm_qc_panels_per_page` | `12` | MA panels per SVG page (`<run>.ma.svg`, `.ma.p2.svg`, ...) |
| `run_replicate_correlation` | `true` | Replicate concordance table and plots (`qc/replicate_correlation/`) |
| `repcor_bin_size` | `null` | Bin size of the replicate metrics; `null` = `tmm_bin_size` |
| `repcor_prior` | `2` | Prior count of the Pearson log2 CPM of the replicate metrics |
| `repcor_signal_frac` | `0.10` | Fraction of bins (highest mean signal) used for the signal-bin metrics |
| `run_deeptools_correlation` | `true` | Keep the deepTools all-sample Spearman heatmap (`qc/correlation/`) |
| `bigwig_bin_size` | `10` | bamCoverage `--binSize` |
| `bigwig_extend_reads` | `true` | bamCoverage `--extendReads` (fragments) |
| `bigwig_units` | `rpkm` | `rpkm`, `cpm` or `none`; same for targets and IgG |
| `run_log2ratio` | `false` | Also write target-vs-control log2ratio bigWigs |
| `call_peaks` | `true` | master switch: `false` = no peak calling (csaw, MACS3, SEACR), FRiP, peak QC or peak annotation; bigWigs, normalisation, QC, heatmaps and tracks still run |
| `run_csaw` | `true` | csaw peak calling per `peaksets.csv` row |
| `csaw_window` | `150` | csaw window width (bp) |
| `csaw_spacing` | `50` | csaw window spacing (bp) |
| `csaw_bin` | `10000` | csaw background bin size (bp) |
| `csaw_norm_bins` | `background` | TMM on `background` bins or `all` |
| `csaw_small_filt` | `3` | Global fallback for blank `small_filt` cells of `peaksets.csv` |
| `csaw_large_filt` | `1.5` | Global fallback for blank `large_filt` cells |
| `csaw_merge_bp` | `50` | Global fallback for blank `merge_bp` cells (bp) |
| `csaw_max_width` | `0` | Global fallback for blank `max_width` cells; `0` = no limit |
| `csaw_fdr` | `0.05` | Global fallback for blank `fdr` cells |
| `csaw_lfc` | `1` | Global fallback for blank `lfc` cells |
| `csaw_defaults` | `K27me3: {small_filt: 3, large_filt: 1.5, merge_bp: 50}`, `K27me2: {small_filt: 2, large_filt: 1, merge_bp: 50}` | Per-antibody fallbacks for blank cells (keys = samplesheet `antibody` of the target group; a config.yaml value replaces the whole map); see "csaw" |
| `csaw_large_window` | `2000` | Width (bp) of the wide window scale |
| `csaw_large_spacing` | `500` | Spacing (bp) of the wide window scale |
| `csaw_bg_max_enrich` | `1` | Background bins: average abundance at most the median plus this value (log2) |
| `csaw_count_filter` | `10` | Minimum summed count over all libraries for a window |
| `csaw_max_frag` | `null` | csaw `max.frag` for CSAW, CSAW_CUTOFF_SCAN, CLASSIFY_REGIONS and BIN_COUNTS (so NORM_QC, TMM and the replicate metrics); `null` = `bowtie2_maxins` (500 if that is null too) |
| `csaw_lfc_from` | `combined` | `combined` (legacy: Simes representative window) or `best` (logFC of the best window) for the `lfc` rule |
| `csaw_disjoint_regions` | `false` | `true` makes csaw's merged regions disjoint (windows assigned by midpoint) |
| `run_csaw_cutoff_scan` | `false` | Post-hoc cutoff diagnostic; never feeds thresholds |
| `csaw_scan_small_grid` | `'1,1.5,2,3'` | Small-filter grid for the scan |
| `csaw_scan_large_grid` | `'0.5,1,1.5,2'` | Large-filter grid for the scan |
| `run_classify` | `false` | Build the region classifications listed in `classifications` |
| `classifications` | `null` | Classifications CSV, see "Region classification"; used when `run_classify: true` |
| `classify_lfc_retained` | `-1` | `lfc_retained` used for blank cells of `classifications.csv` |
| `classify_lfc_lost` | `-2` | `lfc_lost` used for blank cells |
| `classify_fdr` | `0.05` | `fdr` used for blank cells |
| `broad_marks` | `H3K27me3,H3K27me2,H3K36me3,H3K4me1,H3K79me2,H3K79me3,H3K9me1,H3K9me2,H3K9me3,H4K20me1,H3F3A` | Antibodies called as broad marks (comma list; case-insensitive, with or without the histone prefix: `K27me3` = `H3K27me3`). Sources: the ENCODE histone ChIP-seq broad list (https://www.encodeproject.org/chip-seq/histone/), plus H3K27me2 by lab request and H3K9me3 (an ENCODE "exception", treated as broad by St. Jude). Any other antibody is narrow |
| `run_macs3` | `null` | MACS3 peaks with the samplesheet control. `null` (not set) = on when every target has a `control_id`, otherwise skipped with a warning; an explicit `true` with a missing control stops the run, `false` turns it off |
| `macs3_broad` | `auto` | `auto` (broad for `broad_marks`, narrow otherwise), `true` or `false` (force every group) |
| `macs3_qvalue` | `0.01` | `-q` of narrow groups (narrowPeak) |
| `macs3_broad_cutoff` | `0.1` | `--broad --broad-cutoff` of broad groups (the nf-core/cutandrun default) |
| `run_seacr` | `null` | SEACR peaks and per-group consensus; skipped like `run_macs3` when a target has no control |
| `seacr_mode` | `auto` | `auto` (`relaxed` for `broad_marks`, `stringent` otherwise), `relaxed` or `stringent`. Lab policy: SEACR has no broad mode; relaxed = threshold between the knee and the peak of the total-signal curve, stringent = the peak |
| `seacr_min_replicates` | `2` | Replicates required in the consensus |
| `seacr_reciprocal_overlap` | `0.5` | Reciprocal overlap for the consensus |
| `skip_groups` | `''` | Comma list of `group` values kept but skipped by all core steps |
| `run_report` | `true` | Build the HTML QC report |
| `run_heatmaps` | `false` | Build the heatmaps listed in `heatmaps` (stage 2) |
| `heatmaps` | `null` | Heatmaps CSV, see "Heatmaps, tracks, annotation, intersections"; used when `run_heatmaps: true` |
| `heatmap_seed` | `42` | Seed of the heatmap k-means clustering |
| `heatmap_reference_point` | `center` | computeMatrix reference point: `center`, `TSS` or `TES` |
| `heatmap_missing_as_zero` | `true` | computeMatrix `--missingDataAsZero`; `false` omits it |
| `heatmap_average_type_bins` | `mean` | computeMatrix `--averageTypeBins` |
| `heatmap_plot_type` | `se` | Summary plot type of the deepTools heatmap |
| `heatmap_average_type_summary` | `mean` | Summary statistic of the deepTools heatmap |
| `heatmap_height` | `12` | Height of the deepTools heatmap figure |
| `heatmap_width` | `3.5` | Width of the deepTools heatmap figure |
| `heatmap_x_label` | `distance (bp)` | X-axis label of the deepTools heatmap |
| `heatmap_ref_label` | `peak center` | Label of the reference point |
| `heatmap_y_label` | `null` | Profile label; `null` = `normalised <bigwig_units>` |
| `heatmap_legend_location` | `none` | Legend location of the deepTools summary plot |
| `run_tracks` | `false` | Build the pyGenomeTracks figures listed in `tracks` (stage 2) |
| `tracks` | `null` | Tracks CSV; used when `run_tracks: true` |
| `tracks_width` | `15` | pyGenomeTracks figure width (cm) |
| `tracks_dpi` | `600` | pyGenomeTracks figure resolution (dpi) |
| `tracks_max_value_default` | `shared` | Max value of a blank `tracks.csv` entry: `shared` (one group per antibody), `auto`, or a number |
| `tracks_limit_stat` | `mean_bins` | Statistic of the shared (`group:`) track limits: `mean_bins` (largest per-bin mean; EED rule) or `max_bin` (largest value of the bigWig's own bins over the locus; CBX2 rule). The limit is 1.05 x that, rounded up to 2 significant figures, in both |
| `track_number_of_bins` | `700` | Bins of the shared-limit computation, and `number_of_bins` in the ini |
| `track_height` | `1.0` | Height of each sample bigWig track |
| `track_fontsize` | `14` | Font size of the sample titles (accepted; bigWig track titles use pyGenomeTracks' own size, so it no longer changes them) |
| `track_bed_height` | `0.1` | Height of the BED tracks |
| `track_bed_color` | `'#8c96c6'` | Default colour of BED tracks (`bed_colors` overrides it) |
| `track_scalebar_fontsize` | `16` | Font size of the scale bar |
| `track_genes_height` | `3` | Height of the gene track |
| `track_log2_min` | `-6` | Fixed lower limit of log2-ratio tracks (never grouped) |
| `track_log2_max` | `2` | Fixed upper limit of log2-ratio tracks |
| `run_annotation` | `false` | Peak annotation of the peak sets (stage 2) |
| `annotation_colors` | `''` | Colours of the annotation plots: `label=#RRGGBB,label2=#RRGGBB` |
| `annotation_kegg` | `true` | KEGG dotplot; `false` writes a placeholder page |
| `annotation_tss_region` | `3000` | bp on each side of the TSS in the annotation |
| `annotation_extra_beds` | `''` | Comma list of extra regions annotated next to the peak sets: csaw peak set names, class names, intersect names or BED paths (label = the name) |
| `run_intersect` | `false` | BED intersections listed in `intersects` (stage 2) |
| `intersects` | `null` | Intersections CSV; used when `run_intersect: true` |
| `help` | `false` | Print usage and exit |

## Customizing one project's analysis

Every step runs one script from `bin/`. Copy it into the project's `utils/`,
edit it, and add a `config.yaml` key named after the script (file name without
extension):

```yaml
norm_tmm: utils/norm_tmm.R           # replace a step's script for this project only
callpeaks_csaw_args: "..."           # or pass extra arguments to the shipped script
```

Scripts are task inputs, so `-resume` re-runs only the changed step and what is
downstream. Each script's header documents its arguments and output files; a
replacement must keep that contract. Key = file base name, so base names are
unique across `bin/` (`my_theme.R` for R plots, `my_theme_py.py` for Python).

| Group | Scripts (= keys) |
|---|---|
| Preprocess | `fastqc.sh`, `trim.sh`, `align.sh`, `bam_filter.sh`, `markdup.sh`, `dedup_split.sh`, `blacklist_filter.sh`, `bam_finalize.sh`, `bam_prealigned.sh` |
| QC | `qc_alignment_table.R`, `qc_alignment_plot.R`, `qc_fragment_sizes_table.R`, `qc_fragment_sizes_plot.R`, `library_complexity_frags.sh`, `picard_elc.sh`, `qc_library_complexity_table.py`, `qc_library_complexity_plot.py`, `preseq_ccurve.sh`, `qc_preseq_plot.R`, `qc_fingerprint.sh`, `qc_correlation.sh`, `qc_frip.sh`, `qc_frip_table.py`, `qc_tss_enrichment.py`, `qc_tss_table.py`, `qc_peak_stats.py`, `qc_replicate_correlation.R`, `qc_replicate_correlation_plot.R`, `multiqc.sh`, `build_report.py` |
| Normalise | `bin_counts.R`, `norm_chipseqspikeinfree.R`, `norm_spikein.py`, `norm_tmm.R`, `norm_cpm.py`, `norm_combine.py`, `norm_qc.R`, `norm_qc_plot.R` |
| Tracks and peaks | `bigwig.sh`, `bigwig_log2ratio.sh`, `callpeaks_csaw.R`, `csaw_cutoff_scan.R`, `classify_regions.R`, `classify_regions_plot.R`, `callpeaks_macs3.sh`, `seacr_bedgraph.sh`, `seacr_call.sh`, `seacr_consensus.sh` |
| Heatmaps, tracks, annotation, intersections | `heatmap_computematrix.sh`, `heatmap_cluster.py`, `heatmap_plot.sh`, `heatmap_build_matrix.R`, `heatmap_plot_enriched.R`, `heatmap_limits.py`, `make_tracks_ini.py`, `tracks_limits.py`, `tracks_pygenometracks.sh`, `peak_annotation.R`, `intersect_beds.py` |
| Peak QC tracks (see "Peak QC tracks") | `peak_qc_regions.py` |
| Theme | `my_theme.R`, `my_theme_py.py` |

## Containers

Literal `container '...'` strings in `modules/local/*.nf`. Apptainer cache:
`/project/rephilllab/nf_containers`. The `pmacs` profile binds
`/usr/lib64/ossl-modules` (FIPS).

| Used by | Image |
|---|---|
| FASTQC | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` |
| TRIM | `quay.io/biocontainers/trimmomatic:0.40--hdfd78af_0` |
| ALIGN (bowtie2 2.5.5 + samtools) | `/project/rephilllab/nf_containers/cutandrun-py_1.2.sif` |
| samtools steps | `quay.io/biocontainers/samtools:1.23--h96c455f_0` |
| MARKDUP, PICARD_ELC | `quay.io/biocontainers/picard:2.27.5--hdfd78af_0` |
| BLACKLIST_FILTER, FRIP, SEACR_BEDGRAPH, SEACR_CONSENSUS | `quay.io/biocontainers/bedtools:2.31.1--h13024bc_3` |
| FINGERPRINT, CORRELATION, BIGWIG, BIGWIG_LOG2RATIO, COMPUTEMATRIX, COMPUTEMATRIX_CLUSTERED, PLOTHEATMAP | `quay.io/biocontainers/deeptools:3.5.6--pyhdfd78af_0` |
| PYGENOMETRACKS | `quay.io/biocontainers/pygenometracks:3.9--pyhdfd78af_0` |
| MACS3 | `quay.io/biocontainers/macs3:3.0.4--py312h71493bf_0` |
| SEACR_CALL | `quay.io/biocontainers/seacr:1.3--hdfd78af_2` |
| PRESEQ | `quay.io/biocontainers/preseq:2.0.2--gsl1.16_0` |
| MULTIQC | `quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1` |
| R processes (csaw, edgeR, ChIPseqSpikeInFree, plots, BUILD_MATRIX, PLOT_ENRICHED, PEAK_ANNOTATION, NORM_QC_PLOT, CLASSIFY_REGIONS, CLASSIFY_REGIONS_PLOT) | `/project/rephilllab/nf_containers/cutandrun-r_1.3.sif` |
| Python processes (tables, spike-in, report, HEATMAP_CLUSTER, HEATMAP_LIMITS, TRACKS_INI, TRACKS_LIMITS, INTERSECT) | `/project/rephilllab/nf_containers/cutandrun-py_1.2.sif` |

The two custom images are built from the recipes `containers/*.def` and the pinned
`containers/*.env.yml`.

## Software versions

The R image (`cutandrun-r_1.3.sif`) carries R 4.3.3 and Bioconductor 3.18; csaw is
pinned to 1.36.0 in `containers/cutandrun-r.env.yml`, and ChIPseqSpikeInFree to 1.2.4
(commit `f9cc88f`). The Python image is `cutandrun-py_1.2.sif` (`containers/cutandrun-py.env.yml`).
The Docker image used on Windows and Mac is newer (see "Docker versions (Windows, Mac)").
Methods must cite the versions a run actually used, not the manuscript's: read them from
`containers/` and from the step logs in each `results/*/log/` folder.

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

## Outputs

```
results/qc/{fastqc,trimmomatic,bowtie2,markdup,alignment,fragment_sizes,library_complexity,
            fingerprint,correlation,replicate_correlation,frip,peaks,tss_enrichment,multiqc,report}/
results/bam/final/          final BAM + .bai per library
results/bam/spikein/        spike-in BAMs
results/bam/stats/          raw and final flagstat/idxstats, dedup counts
results/scale_factors/      <norm_run>.scale_factors.tsv  (PRIMARY method: sample_id, method, raw_factor, bigwig_multiplier,
                            edger_norm_factor, qc_pass, note), <norm_run>.scale_factors_all.tsv (all methods side by side),
                            methods/<norm_run>.<method>.scale_factors.tsv (one per method)
results/norm_qc/            <norm_run>.norm_qc.tsv, .norm_qc_pairs.tsv (one row per pair and method, with the reference line), .norm_qc_libsizes.tsv, <norm_run>.ma.svg (+ .ma.p2.svg ...)
results/qc/replicate_correlation/  replicate_correlation.tsv, replicate_correlation_pairs_all.tsv, replicate_correlation.svg, replicate_scatter*.svg
metadata/peak_call_modes.tsv  MACS3 / SEACR mode (broad / narrow, relaxed / stringent) per group and why
results/bigwig/             <sample_id>.<method>.<units>.bw
results/peaks/csaw/         <peakset>.peaks.bed, .csaw_results.tsv.gz, .mds.svg, .rds, .norm_sensitivity.tsv (+ cutoff_scan/)
results/peaks/macs3/        <group>_peaks.{narrowPeak,broadPeak}, .xls  (run_macs3; broad for broad_marks, narrow otherwise)
results/peaks/seacr/        <sample_id>.seacr.<mode>.bed, <group>.seacr_consensus.bed  (run_seacr)
results/classify/<id>/      <id>.{retained,lost,unclassified}.bed, <id>.per_region.tsv.gz, <id>.summary.tsv, <id>.ma.svg  (run_classify)
results/heatmaps/{matrix,cluster,plots,limits}/   matrices; cluster assignments, k table, cluster BEDs; heatmap SVGs; limits/<heatmap>_limits.tsv (auto limits)  (run_heatmaps)
results/tracks/{ini,plots,limits}/ generated .ini and region tables; figure SVGs; limits/<slug>.resolved.ini, <slug>_limits.tsv  (run_tracks)
results/annotation/         annotation tables, plots, KEGG page  (run_annotation)
results/intersect/          <name>.intersect.bed, <name>.only_a.bed, <name>.intersect_summary.tsv  (run_intersect)
intermediate_files/{trimmed,aligned,filtered,markdup,seacr}/
metadata/peaksets.resolved.tsv  resolved csaw thresholds of every peak set row and their source (row, antibody or global)
metadata/pipeline_info/     execution report, trace, timeline, DAG (a resumed run keeps the previous run's as execution_trace.<timestamp>.txt)
```

Every step's `log/` subfolder holds `<name>_command.sh`, the exact command that
ran. Start with `results/qc/report/<project_name>_qc_report.html`.

## Docker versions (Windows, Mac)

The same pipeline runs on a personal computer with Docker: see `../WINDOWS/README.md`
(Windows) and `../MAC/README.md` (Mac). Those folders use one shared image
(`docker/image.txt`) instead of the per-step containers above, and the local executor
instead of LSF. Their `main.nf`, `modules/` and `bin/` match this folder's.

## Status and limits

- Version 0.1.0, local git only. Stages 1a and 2 are in scope (stage 2: heatmaps,
  tracks, annotation, intersections, all optional).
- Paired-end only. Pre-aligned `bam` rows must already be final-filtered and
  cannot take part in spike-in normalisation.
- `bowtie2_minins` and `bowtie2_maxins` default to `null`: no `-I`/`-X` is passed, so bowtie2's own defaults (-I 0, -X 500) apply, as in the manuscript. `10` / `700` were used in 2026-10 test runs and are optional toggles.
- Example projects point at the lab's original data paths (replace them with your own)
  and carry inferred values (see each example README).
