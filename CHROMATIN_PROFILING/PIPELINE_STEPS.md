# CUT&RUN / chromatin profiling pipeline — every step, tool and default

Reference for the pipeline in `HPC/`, `WINDOWS/` and `MAC/` (the step scripts in `bin/` are
identical in the three; `main.nf` and `modules/` are the same except that the Docker versions run
steps 5-8 as one process, BAM_PROCESS, see the rows below). Compiled 2026-10-26 from the code (`main.nf`,
`nextflow.config`, `conf/base.config`, `modules/local/*.nf`, `bin/*`), not from older docs.

How to read it:

- **Order** is the order in which `main.nf` wires the steps. Nextflow runs them as soon as their
  inputs exist, so independent steps overlap (e.g. QC runs beside normalisation).
- **Default** = what happens with the pipeline's own defaults (`nextflow.config`). The example
  dataset configs (`examples/projects/*/config.yaml`) override several of them (e.g. bigWig bin 50 bp,
  `csaw_norm_bins: all`, `csaw_max_frag: 500`, heatmaps/tracks on).
- Values in `<angle brackets>` come from the CSV sheets or `config.yaml`. `cpus`/`memory` are the
  cluster requests (`conf/base.config`); the Docker versions cap them to the run's budget (see the end).
- Every step's exact command is published next to its output as `.../log/*_command.sh`.

## Stage 1 — FASTQ to final BAMs

| # | Step | Runs by default? | Tool and command (defaults filled in) | Output |
|---|---|---|---|---|
| 1 | FASTQC — read quality per FASTQ | yes (FASTQ rows) | FastQC 0.12.1: `fastqc -o . --extract --threads 2 <fastq>` | `qc/fastqc/` |
| 2 | TRIM — adapter trimming | yes (`trim: true`) | Trimmomatic 0.40: `trimmomatic PE -threads <cpus=8> R1 R2 S.R1.trim.fastq.gz /dev/null S.R2.trim.fastq.gz /dev/null ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:True MINLEN:25`. Adapter file `NexteraPE-PE.fa` when `experiment: cutandtag`; `MINLEN` = `trim_minlen` (25). No quality trimming (no LEADING/TRAILING/SLIDINGWINDOW). Unpaired reads discarded. | `intermediate_files/trimmed/` (Docker versions: only with `publish_intermediates: true`), `qc/trimmomatic/` |
| 3 | ALIGN | yes | bowtie2 2.5.5 piped to samtools: `bowtie2 -p <cpus> --local --very-sensitive-local -q --no-mixed --no-discordant --dovetail -x <genome_index> -1 R1 -2 R2 \| samtools view -b`. `-I`/`-X` only when `bowtie2_minins`/`bowtie2_maxins` are set (default null = bowtie2's own 0 / 500, the manuscript command). `bowtie2_mode: end_to_end` switches to `--end-to-end --very-sensitive`. | `intermediate_files/aligned/` (Docker versions: the unfiltered aligned BAM, only with `publish_intermediates: true`), `qc/bowtie2/*.bowtie2.log` |
| 4 | BAM_FILTER — sort, raw stats, filter | yes | samtools: `sort -m 1G -@ 4`; `index`; `flagstat`/`idxstats` (raw); `view -b -q 30 -f 3 -F 2828` = MAPQ ≥ `mapq` (30), paired + properly paired, drop unmapped/mate unmapped/secondary/QC-fail/supplementary | `intermediate_files/filtered/` (Docker versions: only with `publish_intermediates: true`), `bam/stats/*.raw.*` |
| 5 | MARKDUP — mark duplicates (HPC: own task; WINDOWS and MAC: inside BAM_PROCESS) | yes | Picard 2.27.5: `picard -Xmx<memory−2 GB> MarkDuplicates OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 REMOVE_DUPLICATES=false ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true` | `intermediate_files/markdup/` (HPC only; in the Docker versions the markdup BAM is deleted in BAM_PROCESS as soon as DEDUP_SPLIT has read it), `qc/markdup/*.markdup_metrics.txt` |
| 6 | DEDUP_SPLIT — split spike-in, remove duplicates (HPC: own task; WINDOWS and MAC: inside BAM_PROCESS) | yes | samtools: contigs starting with `spike_prefix` (`ecoli_chr`) → spike-in BAM (`-F 1024` when `spikein_dedup: true`); other contigs except `mito_chrom` (`chrM`) → target BAM (`-F 1024` when `remove_duplicates: true`, and always for IgG/controls), name-sorted | `bam/spikein/`, `bam/stats/*.dedup_counts.tsv` (in the Docker versions the name-sorted target BAM is deleted in BAM_PROCESS as soon as BLACKLIST_FILTER has read it) |
| 7 | BLACKLIST_FILTER (HPC: own task; WINDOWS and MAC: inside BAM_PROCESS) | yes | bedtools 2.31.1: `bedtools pairtobed -abam <target> -b <blacklist> -type neither` (a pair is dropped if either mate overlaps) | (intermediate; in the Docker versions the blacklist-filtered name-sorted BAM is deleted in BAM_PROCESS as soon as BAM_FINALIZE has read it) |
| 8 | BAM_FINALIZE (HPC: own task; WINDOWS and MAC: inside BAM_PROCESS) | yes | samtools: `sort -m 1G -@ 4`, `index`, `flagstat`, `idxstats`; orphan-mate check (exit 3 if any read has no mate) | `bam/final/*.final.bam(.bai)`, `bam/stats/*.final.*` |
| 9 | BAM_PREALIGNED — samplesheet `bam` rows instead of FASTQs | only if `bam` rows | samtools: copy if coordinate-sorted, else `sort -m 1G`; `index`, `flagstat`, `idxstats` (no filtering: the BAM is used as given) | `bam/final/` |

## QC (on the final or filtered BAMs)

| # | Step | Runs by default? | Tool and command | Output |
|---|---|---|---|---|
| 10 | QC_ALIGNMENT_TABLE | yes | R: read counts at each stage from flagstat/idxstats, spike-in and chrM fractions | `qc/alignment/alignment_stats.tsv` |
| 11 | QC_ALIGNMENT_PLOT | yes | R/ggplot2: 7 plots of the table | `qc/alignment/*.svg` |
| 12 | QC_FRAGMENT_SIZES_TABLE | yes | R/Rsamtools: \|TLEN\| of first mates of proper pairs (primary), up to 2,000 bp | `qc/fragment_sizes/` |
| 13 | QC_FRAGMENT_SIZES_PLOT | yes | R: density + boxplot (x up to 1,000 bp) | `qc/fragment_sizes/*.svg` |
| 14 | LIBCOMPLEXITY_FRAGS | yes (FASTQ samples) | samtools + awk on the filtered (pre-dedup) BAM: fragments from first mates (`-f 66 -F 2308`) → total, distinct, once, twice → NRF, PBC1, PBC2 | `qc/library_complexity/` |
| 15 | PICARD_ELC | yes (FASTQ samples) | Picard `EstimateLibraryComplexity VALIDATION_STRINGENCY=SILENT` on the filtered BAM | `qc/library_complexity/` |
| 16 | QC_LIBCOMPLEXITY_TABLE | yes | Python: joins steps 14, 15 and the MarkDuplicates metrics | `qc/library_complexity/library_complexity.tsv` |
| 17 | QC_LIBCOMPLEXITY_PLOT | yes | matplotlib: NRF, PBC1, PBC2, Picard % duplication | `qc/library_complexity/*.svg` |
| 18 | PRESEQ | yes (FASTQ samples) | preseq 2.0.2: `preseq c_curve -B -P -v -s 1000000` on the filtered BAM (re-run with a smaller step if < 10 points) | `qc/preseq/` |
| 19 | QC_PRESEQ_PLOT | yes | R: complexity curves | `qc/preseq/preseq_ccurve.svg` |
| 20 | FINGERPRINT — one per group | yes | deepTools 3.5.6: `plotFingerprint --extendReads --skipZeros --numberOfSamples 500000 --binSize 500` | `qc/fingerprint/` |
| 21 | CORRELATION | yes (`run_deeptools_correlation: true`) | deepTools: `multiBamSummary bins --binSize 5000 --distanceBetweenBins 0 --extendReads`; `plotCorrelation --corMethod spearman --whatToPlot heatmap --skipZeros --removeOutliers --colorMap RdYlBu --plotNumbers` (clustered all-sample view; its counting differs from 21a, so the numbers differ) | `qc/correlation/` |
| 21a | BIN_COUNTS — one task per distinct bin size | yes | R/csaw: `windowCounts(bin=TRUE, width=<bin size>, filter=0, readParam(minq=<mapq>, pe="both", max.frag=<csaw_max_frag, else bowtie2_maxins, else 500>, dedup=<remove_duplicates>))` over every active library (IgG included; `skip_groups` excluded); all bins kept, with a `blacklisted` flag (any overlap with the blacklist). Bin sizes: `tmm_bin_size` (10000), `repcor_bin_size` (null = `tmm_bin_size`) and the `bin_size` of TMM-primary runs. Also writes per-library totals (csaw total, final-BAM mapped reads and fragments). Feeds 21b, 25 and 26. | `intermediate_files/bin_counts/bins_<size>.counts.tsv.gz`, `.libsizes.tsv` (Docker versions: only with `publish_intermediates: true`; HPC always) |
| 21b | REPLICATE_CORRELATION | yes (`run_replicate_correlation: true`) | R: per pair of libraries of the same `group`, over non-blacklisted bins (`repcor_bin_size`, default 10 kb) with reads: Pearson of log2 CPM (`repcor_prior` 2), Spearman of the counts, the same Pearson on the top `repcor_signal_frac` (0.10) of bins by mean signal, robust SD of the log2 ratio on those bins, and the median Pearson with the other groups of the same antibody (context). No thresholds. Groups with one library get a `single_library` row. | `qc/replicate_correlation/replicate_correlation.tsv`, `replicate_correlation_pairs_all.tsv` |
| 21c | REPLICATE_CORRELATION_PLOT | same | R: one panel per antibody (replicate pairs vs non-replicate pairs, genome-wide and signal bins) and replicate scatter plots (12 per page) | `qc/replicate_correlation/*.svg` |
| 22 | MULTIQC | yes | MultiQC 1.35: `multiqc . --force` over FastQC, Trimmomatic, bowtie2, MarkDuplicates, flagstat/idxstats, preseq | `qc/multiqc/multiqc_report.html` |

## Normalisation (one task per row of `norm_runs.csv` and method; samples pick a `norm_run` in the samplesheet)

`method` of `norm_runs.csv` is the run's **primary** method (chipseqspikeinfree, spikein, tmm or cpm); only it feeds the bigWigs and
CLASSIFY_REGIONS. Every run with at least 2 libraries also computes the methods of `norm_methods` (default
`chipseqspikeinfree,tmm,cpm`) and spike-in when `norm_spikein: true`, as diagnostics for NORM_QC. A run with one library needs primary `cpm`.

| # | Step | Runs by default? | Tool and command | Output |
|---|---|---|---|---|
| 23 | CHIPSEQSPIKEINFREE | every norm run with ≥ 2 libraries (primary, or diagnostic via `norm_methods`) | ChIPseqSpikeInFree 1.2.4 (commit f9cc88f): `ChIPseqSpikeInFree(bams, chromFile=<chrom_sizes>, binSize=<bin_size>, cutoff_QC=<cutoff>, maxLastTurn=<max_turns>[, minFirstTurn=<min_first_turn>], singleEnd=TRUE, ncores=<cpus>)`. As primary the CSV columns `bin_size, cutoff, max_turns` are required; as a diagnostic blank cells use `csif_bin_size` (10000), `csif_cutoff` (1.2), `csif_max_turns` (0.99). Counting workers = `--cpus`. A sample failing QC gets the largest finite factor (imputed). A diagnostic that fails writes `UNAVAILABLE:` rows (exit 0); a primary that fails completely exits 3. `bigwig_multiplier = 1/SF`; `edger_norm_factor = SF / geometric mean(SF)`. Requests: 8 cpus, 8 GB × attempt. | `scale_factors/methods/<run>.chipseqspikeinfree.scale_factors.tsv`, `scale_factors/<run>.chipseqspikeinfree.svg` |
| 24 | SPIKEIN | when `method: spikein` or `norm_spikein: true` | Python: raw = (spike + target reads) / spike reads; `bigwig_multiplier = raw / min(raw)`; `edger_norm_factor` = (min/raw) / geometric mean; fewer than 1,000 spike reads is flagged (`qc_fail`). Reads come from the final target idxstats and the spike-in BAM's idxstats (`bam/stats/*.spikein.idxstats.txt`). As a diagnostic: 0 spike reads → `UNAVAILABLE`; a run with pre-aligned BAMs is not launched (`unavailable`, startup warning). | `scale_factors/methods/<run>.spikein.scale_factors.tsv` |
| 25 | TMM | every norm run with ≥ 2 libraries | R/csaw from the BIN_COUNTS counts (the run's libraries, bins with run sum ≥ 10): `normFactors(weighted=FALSE)` on `lib.size` = csaw total; `bigwig_multiplier = 1/nf`. Bin size: a TMM primary uses its CSV `bin_size` (blank = `tmm_bin_size` 10000), a diagnostic always `tmm_bin_size`. `dedup` = `remove_duplicates`; `max.frag` = `csaw_max_frag`, else `bowtie2_maxins`, else 500. Counts mode gives the same factors as counting the BAMs directly (1e-12). Requests: 1 cpu, 4 GB × attempt. | `scale_factors/methods/<run>.tmm.scale_factors.tsv` |
| 25a | NORM_CPM | every norm run with ≥ 2 libraries (or primary `cpm`) | Python: library size only; every library `raw_factor = bigwig_multiplier = edger_norm_factor = 1` (bigWig = bamCoverage per million final-BAM mapped reads, `--scaleFactor 1`; edgeR `lib.size` = csaw total, `norm.factors` = 1). 1 cpu, 1 GB. | `scale_factors/methods/<run>.cpm.scale_factors.tsv` |
| 25b | NORM_COMBINE | per norm run | Python: waits for every method table of the run; writes the PRIMARY table unchanged (byte copy) and the side-by-side table `<run>.scale_factors_all.tsv` (columns `<method>_{raw_factor,bigwig_multiplier,edger_norm_factor,qc_pass,status,note}`; status ok, qc_fail, imputed, unavailable, not_run). Exit 3 when the primary is unavailable or has a missing/NaN/≤ 0 factor. Diagnostics never delay the primary: bigWigs read the primary table only. | `scale_factors/<run>.scale_factors.tsv`, `scale_factors/<run>.scale_factors_all.tsv` |
| 26 | NORM_QC | when any norm run with ≥ 2 libraries exists | R from the BIN_COUNTS counts at `tmm_bin_size` (no BAMs): for replicate pairs (same `group`) and treated-vs-reference pairs (reference = `ref_condition` of `norm_runs.csv`, else the condition of the `ref_group` of the classifications using the run, else the first library's), M = log2 CPM of `other` − `ref` with prior `norm_qc_prior` (0.5) on bins outside the blacklist with pair sum ≥ `norm_qc_min_count` (20). One row per pair × method with the reference line `log2(f_other / f_ref)` from `edger_norm_factor` (0 for CPM), `bigwig_offset` (bigWig vs edgeR library size), `line_se` (spike-in) and the cloud centre. Requests: 1 cpu, 8 GB × attempt. | `norm_qc/<run>.norm_qc_pairs.tsv`, `.norm_qc_libsizes.tsv`, `.norm_qc.tsv` |
| 27 | NORM_QC_PLOT | same | R: one MA panel per pair, 12 per page (`norm_qc_panels_per_page`); dashed line per method (Okabe-Ito colours; dotted = factor from a QC-failed or imputed library; primary thicker), the line values written in each panel, an `unavailable` method shown in grey in the legend | `norm_qc/<run>.ma.svg`, `<run>.ma.p2.svg`, ... |

## bigWigs and TSS enrichment

| # | Step | Runs by default? | Tool and command | Output |
|---|---|---|---|---|
| 28 | (factor table) | always | each sample's `bigwig_multiplier` from its norm run; IgG/controls without a norm run → 1 (`unscaled`); no norm runs at all → 1 (`cpm`) | — |
| 29 | BIGWIG | yes | deepTools: `bamCoverage --outFileFormat bigwig --binSize 10 --skipNonCoveredRegions --normalizeUsing RPKM --scaleFactor <multiplier> --extendReads --effectiveGenomeSize 3100000000` (`bigwig_bin_size` 10 — the dataset configs use 50; `bigwig_units` rpkm/cpm/none; `bigwig_extend_reads`) | `bigwig/<sample>.<method>.rpkm.bw` |
| 30 | BIGWIG_LOG2RATIO | no (`run_log2ratio: false`) | deepTools `bigwigCompare -b1 <target> -b2 <IgG> --operation log2 --skipZeroOverZero --binSize 10` | `bigwig/*.log2ratio.bw` |
| 31 | TSS_ENRICHMENT | yes | Python (`qc_tss_enrichment.py`, pyBigWig + numpy; same numbers as deepTools `computeMatrix reference-point --referencePoint TSS -b 2000 -a 2000 --binSize 50 --skipZeros --missingDataAsZero`, but TSSs deduplicated by chromosome, reference position and strand; TSS = transcript start (`+`) or end (`-`) from the configured `gtf`, not the `tss_bed` promoter windows, which are only a fallback when no GTF is set); profile plot (SVG); enrichment = max mean signal within ±250 bp ÷ mean of the outer 200 bp flanks | `qc/tss_enrichment/` |
| 32 | QC_TSS_TABLE | yes | Python: one table of all samples | `qc/tss_enrichment/tss_enrichment.tsv` |

## Stage 2 — peaks, classification, heatmaps, tracks, annotation

| # | Step | Runs by default? | Tool and command | Output |
|---|---|---|---|---|
| 33 | CSAW — one per row of `peaksets.csv` | yes (`run_csaw: true`; needs `peaksets.csv`) | csaw/edgeR: `readParam(minq=30, pe="both", max.frag=<csaw_max_frag, else bowtie2_maxins, else 500>, dedup=TRUE)`; small windows `windowCounts(width=150, spacing=50, filter=10)`; large windows width 2000, spacing 500, filter 10; background `windowCounts(bin=TRUE, width=10000)`. TMM factors (`calcNormFactors(method="TMM", doWeighting=FALSE)`) on `csaw_norm_bins`: `background` (default: bins with aveLogCPM ≤ median + 1) or `all` (the manuscript). Windows kept when `filterWindowsGlobal` enrichment > `small_filt` / > `large_filt`. Design `~0 + group` (target vs IgG); with replicates `estimateDisp` + `glmQLFit` + `glmQLFTest`, without replicates dispersion 0.05 + `glmLRT`. `mergeResultsList(equiweight=TRUE, tol=merge_bp[, max.width])`. Peak = combined FDR ≤ `fdr` (0.05) and logFC > `lfc` (1) (`csaw_lfc_from: combined`). Thresholds per peak set: the `peaksets.csv` cell, else `csaw_defaults[antibody]` (K27me3 3 / 1.5 / merge 50; K27me2 2 / 1 / merge 50), else the global `csaw_*` keys (small 3, large 1.5, merge 50, max_width 0 = none). BED score = min(1000, round(−10·log10 FDR)). | `peaks/csaw/<set>.peaks.bed`, `.csaw_results.tsv.gz`, `.rds`, `.mds.svg`, `.norm_sensitivity.tsv` |
| 34 | CSAW_CUTOFF_SCAN | no (`run_csaw_cutoff_scan: false`) | same counting; peak numbers over the grids small `1,1.5,2,3` × large `0.5,1,1.5,2` (diagnostic only) | `peaks/csaw/cutoff_scan/` |
| 35 | MACS3 — one per target group | yes (`run_macs3` null = on when every target has a control; skipped with a warning when a target has no control and the key was not set explicitly) | MACS3 3.0.4: `macs3 callpeak -f BAMPE -g 3100000000 --keep-dup all -t <targets> -c <their IgG>`, broad groups (`broad_marks`, `macs3_broad: auto`) with `--broad --broad-cutoff 0.1` (`macs3_broad_cutoff`), narrow groups with `-q 0.01` (`macs3_qvalue`, narrowPeak). `macs3_broad: true/false` forces every group. Mode per group: `metadata/peak_call_modes.tsv`. | `peaks/macs3/` |
| 36 | SEACR_BEDGRAPH | yes (`run_seacr` null; same skip rule as 35) | bedtools: `genomecov -ibam <bam> -bg -pc` (fragment pileup), sorted | `intermediate_files/seacr/` (Docker versions: only with `publish_intermediates: true`; HPC always) |
| 37 | SEACR_CALL — per target sample | with 36 | SEACR 1.3: `SEACR_1.3.sh <target.bg> <IgG.bg> norm <relaxed or stringent> <sample>` (`seacr_mode: auto` = relaxed for `broad_marks`, stringent otherwise; or forced). Lab policy: SEACR has no broad mode; relaxed = threshold between the knee and the peak of the total-signal curve, stringent = the peak. | `peaks/seacr/<sample>.seacr.<mode>.bed` |
| 38 | SEACR_CONSENSUS — per group | no | bedtools: a replicate region is kept when supported by ≥ `seacr_min_replicates` (2) replicates at reciprocal overlap ≥ 0.5 (`intersect -f 0.5 -r -c`), then `merge` | `peaks/seacr/<group>.seacr_consensus.bed` |
| 39 | CLASSIFY_REGIONS — per row of `classifications.csv` | no (`run_classify: false`) | csaw `regionCounts` (minq 30, pe both, max.frag as in 33, dedup TRUE); edgeR with the norm run's `edger_norm_factor`; design `~group` (reference first), QL F-test (LRT, dispersion 0.05, without replicates); BH FDR over all regions. **Lost** = logFC ≤ −2 and FDR < 0.05; **retained** = round(logFC, 1) ≥ −1 or FDR ≥ 0.05; else **unclassified**; a region meeting both lost and retained is lost. Thresholds: the CSV cell, else `classify_lfc_retained` −1, `classify_lfc_lost` −2, `classify_fdr` 0.05. | `classify/<id>/<id>.{retained,lost,unclassified}.bed`, `.per_region.tsv.gz`, `.summary.tsv` |
| 40 | CLASSIFY_REGIONS_PLOT | with 39 | R: MA plot with the two logFC lines | `classify/<id>/<id>.ma.svg` |
| 41 | FRIP — per sample × peak source | when any peak set exists | bedtools: reads in peaks (`intersect -u -abam`) ÷ mapped reads (flagstat) | `qc/frip/` |
| 42 | QC_FRIP_TABLE | same | Python: one table | `qc/frip/frip.tsv` |
| 43 | QC_PEAK_STATS | same | Python: n peaks, total bp, width summary, n chromosomes, % autosomal per peak file | `qc/peaks/peak_stats.tsv` |
| 44 | INTERSECT — per row of `intersects.csv` | no (`run_intersect: false`) | bedtools: `intersect -a A -b B -wa -u -f <min_overlap=0.05> [-r if reciprocal]` (`output: records`, default; `output: overlap` drops `-wa -u` = legacy overlap pieces); only-A with `-v` | `intersect/` |
| 45 | COMPUTEMATRIX — per row of `heatmaps.csv` (deepTools engine) | no (`run_heatmaps: false`) | deepTools: `computeMatrix reference-point --referencePoint center -b <window_bp> -a <window_bp> -bs <bin_bp> --averageTypeBins mean --missingDataAsZero` over the normalised bigWigs (window and bin come from the CSV; no default) | `heatmaps/matrix/` |
| 46 | BUILD_MATRIX (EnrichedHeatmap engine) | rows with `engine: enrichedheatmap` | R/EnrichedHeatmap: `normalizeToMatrix(include_target=FALSE, mean_mode="w0", smooth=FALSE, w=<bin_bp>, extend=<window_bp>)` around region centres | `heatmaps/matrix/` |
| 47 | HEATMAP_CLUSTER | rows with `cluster_k` | scikit-learn KMeans (`k-means++`, `lloyd`, `random_state = heatmap_seed` 42, `n_init` 100, `max_iter` 500); a k range is chosen by gap statistic / silhouette | `heatmaps/cluster/` |
| 48 | COMPUTEMATRIX_CLUSTERED | with 47 | as 45, with the cluster BEDs as region groups (deepTools' unseeded `--kmeans` is not used) | `heatmaps/matrix/` |
| 49 | HEATMAP_LIMITS | when `z_max` or `y_max` is `auto` | Python: per antibody, limit = 2-significant-figure round-up of 1.05 × the highest mean profile | `heatmaps/limits/` |
| 50 | PLOTHEATMAP | with 45 | deepTools: `plotHeatmap --plotType se --averageTypeSummaryPlot mean --heatmapHeight 12 --heatmapWidth 3.5 --colorList white,<colour> --legendLocation none -x 'distance (bp)' -y 'normalised rpkm' --refPointLabel 'peak center' [--sortUsingSamples <sort_using>]` (dataset configs use the y label `sf x rpkm`) | `heatmaps/plots/*.svg` |
| 51 | PLOT_ENRICHED | EnrichedHeatmap rows | R/EnrichedHeatmap: heatmap split by clusters, sorted by `sort_using` | `heatmaps/plots/` |
| 52 | PEAK_QC_REGIONS | **yes, whenever a peak set exists** (`run_peak_qc: true`; needs `gtf`; csaw, MACS3 and SEACR peak sets all count) | Python: up to 15 regions (`peak_qc_max_regions`): user loci (`peak_qc_loci`; with a Polycomb mark and no user loci, HOXA, HOXD, CDKN2A by default); per mark the widest, median-width and densest-peak (300 kb) windows; each classification's minority class (≤ 50 % of peaks) and cross-mark discordant region; a negative control (`peak_qc_negative: auto` = peak-free housekeeping gene, else a data-driven long gene, signal ≤ 1.5× the genome median in every mark); the peak-free 250 kb tile with the highest signal (tiles with an IgG fold above `peak_qc_nopeak_igg_max_fold`, default 1.5, are excluded when an IgG bigWig exists; null or `''` = off). Window 100/200/400 kb by anchor width; ≤ 10 % blacklisted; ≤ 20 % overlap with earlier windows. | `tracks/peak_qc/regions.tsv`, `regions.detail.tsv`, `params.json` |
| 53 | TRACKS_INI — per track set (`tracks.csv` rows + the automatic `peak_qc` set) | no for `tracks.csv` (`run_tracks: false`); **yes for `peak_qc`** when 52 runs and `run_tracks` is on | Python: pyGenomeTracks `.ini`: bigWig tracks `type = fill`, 700 bins, height 1.0, BED tracks (peaks/classes; for `peak_qc` every caller's peak sets: `csaw <set>` black, `MACS3 <group>` #0072B2, `SEACR <group>` #D55E00), gene track from the GTF; y-range per antibody shared across samples (`tracks_max_value_default: shared`) | `tracks/ini/` |
| 54 | TRACKS_LIMITS — per locus | with 53 | Python: shared y-limit = 2-s.f. round-up of 1.05 × the largest per-bin mean (`tracks_limit_stat: mean_bins`; or `max_bin`) | `tracks/limits/` |
| 55 | PYGENOMETRACKS — per locus | with 53 | pyGenomeTracks 3.9: `pyGenomeTracks --tracks <ini> --region <locus> --width 15 --dpi 600` (`tracks_width`, `tracks_dpi`), SVG | `tracks/plots/<set>.<locus or label>.svg` |
| 56 | PEAK_ANNOTATION | no (`run_annotation: false`) | ChIPseeker: TxDb from the GTF (`makeTxDbFromGFF`), `annotatePeak(tssRegion = ±3000)` (`annotation_tss_region`), promoters ±3 kb, distance to TSS, widths, overlap matrix; KEGG (`annotation_kegg: true`): clusterProfiler `compareCluster(enrichKEGG, organism="hsa", pvalueCutoff=0.05, BH)`; GC/CpG only if a T2T BSgenome is installed (it is not in the omics image → NA) | `annotation/` |

## Report

| # | Step | Runs by default? | What |
|---|---|---|---|
| 57 | REPORT | yes (`run_report: true`) | self-contained HTML: alignment, fragment sizes, library complexity, fingerprint, correlation, replicate concordance, FRiP, peaks, TSS enrichment, normalisation QC → `qc/report/<project>_qc_report.html` |

## What each dataset can change

Everything below is set per dataset (per project): keys in the project's `config.yaml` override the
pipeline defaults, and the CSV sheets carry per-row settings. Nothing in the pipeline code is specific
to one dataset.

### `config.yaml` keys (default in brackets)

| Group | Keys |
|---|---|
| Inputs (CSV file names) | `samplesheet` [samplesheet.csv], `norm_runs` [norm_runs.csv], `peaksets` [peaksets.csv], `classifications`, `heatmaps`, `tracks`, `intersects` [all null = step off] |
| Assay and reference | `experiment` [cutandrun; cutandtag switches the adapter file], `genome_index`, `chrom_sizes`, `blacklist`, `gtf`, `tss_bed` [paths; `tss_bed` is a fallback for TSS enrichment, required only when `gtf` is unset], `spike_prefix` [ecoli_chr], `mito_chrom` [chrM], `effective_genome_size` [3100000000], `paired_end` [true; single-end not supported in stage 1] |
| Trimming and alignment | `trim` [true], `trim_minlen` [25], `bowtie2_mode` [local; or end_to_end], `bowtie2_minins` [null], `bowtie2_maxins` [null], `mapq` [30] |
| Duplicates | `remove_duplicates` [true; IgG always deduplicated], `spikein_dedup` [true] |
| Normalisation | `tmm_bin_size` [10000]; `normalisation` [chipseqspikeinfree] is informational only: the primary method comes from `norm_runs.csv`; `norm_methods` ['chipseqspikeinfree,tmm,cpm'; '' = primary only], `norm_spikein` [false], `csif_bin_size` [10000], `csif_cutoff` [1.2], `csif_max_turns` [0.99] (ChIPseqSpikeInFree as a diagnostic), `norm_qc_min_count` [20], `norm_qc_prior` [0.5], `norm_qc_panels_per_page` [12] |
| Replicate concordance | `run_replicate_correlation` [true], `repcor_bin_size` [null = `tmm_bin_size`], `repcor_prior` [2], `repcor_signal_frac` [0.10], `run_deeptools_correlation` [true] |
| bigWigs | `bigwig_bin_size` [10], `bigwig_extend_reads` [true], `bigwig_units` [rpkm; cpm, none], `run_log2ratio` [false] |
| Peak calling master switch | `call_peaks` [true]: `false` = no peak calling: steps 33-38 (csaw, CSAW_CUTOFF_SCAN, MACS3, SEACR), 41-43 (FRiP, peak stats), 52 (peak QC) and 56 (annotation) are off; everything else still runs (heatmaps, tracks and classifications then need BED paths, not peak set names) |
| csaw peaks | `run_csaw` [true], `csaw_window` [150], `csaw_spacing` [50], `csaw_bin` [10000], `csaw_norm_bins` [background; all], `csaw_large_window` [2000], `csaw_large_spacing` [500], `csaw_bg_max_enrich` [1], `csaw_count_filter` [10], `csaw_max_frag` [null → bowtie2_maxins → 500], `csaw_lfc_from` [combined; best], `csaw_disjoint_regions` [false]; threshold fallbacks `csaw_small_filt` [3], `csaw_large_filt` [1.5], `csaw_merge_bp` [50], `csaw_max_width` [0 = none], `csaw_fdr` [0.05], `csaw_lfc` [1]; per-antibody fallbacks `csaw_defaults` [K27me3 3/1.5/50, K27me2 2/1/50; a block in config.yaml replaces the whole map]; `run_csaw_cutoff_scan` [false], `csaw_scan_small_grid`, `csaw_scan_large_grid` |
| Classification | `run_classify` [false], `classify_lfc_retained` [−1], `classify_lfc_lost` [−2], `classify_fdr` [0.05] (fallbacks for blank CSV cells) |
| Broad / narrow | `broad_marks` ['H3K27me3,H3K27me2,H3K36me3,H3K4me1,H3K79me2,H3K79me3,H3K9me1,H3K9me2,H3K9me3,H4K20me1,H3F3A'; case-insensitive, histone prefix optional] |
| MACS3 | `run_macs3` [null = on when every target has a control; skipped with a warning without controls unless set explicitly], `macs3_broad` [auto; true, false], `macs3_qvalue` [0.01], `macs3_broad_cutoff` [0.1] |
| SEACR | `run_seacr` [null; same skip rule], `seacr_mode` [auto; relaxed, stringent], `seacr_min_replicates` [2], `seacr_reciprocal_overlap` [0.5] |
| Samples | `skip_groups` [`''` = none skipped; comma list of `group` values] |
| Heatmaps | `run_heatmaps` [false], `heatmap_seed` [42], `heatmap_reference_point` [center; TSS, TES], `heatmap_missing_as_zero` [true], `heatmap_average_type_bins` [mean], `heatmap_plot_type` [se], `heatmap_average_type_summary` [mean], `heatmap_height` [12], `heatmap_width` [3.5], `heatmap_x_label`, `heatmap_ref_label` [peak center], `heatmap_y_label` [null = normalised <units>], `heatmap_legend_location` [none] |
| Tracks | `run_tracks` [false], `tracks_width` [15 cm], `tracks_dpi` [600], `tracks_max_value_default` [shared; auto, a number], `tracks_limit_stat` [mean_bins; max_bin], `track_number_of_bins` [700], `track_height` [1.0], `track_fontsize` [14], `track_bed_height` [0.1], `track_bed_color` [#8c96c6], `track_scalebar_fontsize` [16], `track_genes_height` [3], `track_log2_min` [−6], `track_log2_max` [2] |
| Peak QC tracks | `run_peak_qc` [true], `peak_qc_loci` [''], `peak_qc_negative` [auto; none, a locus], `peak_qc_max_regions` [15], `peak_qc_exclude_bed` [null], `peak_qc_nopeak_igg_max_fold` [1.5; null or '' = off], `peak_qc_shape_marks` [''] |
| Annotation and intersects | `run_annotation` [false], `annotation_colors` [''], `annotation_extra_beds` [''], `annotation_kegg` [true], `annotation_tss_region` [3000], `run_intersect` [false] |
| Report and outputs | `run_report` [true], `project_name`, `outdir` [./results], `metadata_dir` [./metadata], `intermediate_dir` [./intermediate_files], `publish_intermediates` [Docker versions only; false], `publish_mode` [cluster copy; Docker auto] |
| This computer (Docker versions only) | `max_cpus`, `max_memory_gb` [unset = automatic budget], `keep_work` [true], `wsl_resize` [Windows only; true] |
| Any step | `<script name>: utils/<file>` replaces that step's script for this project; `<script name>_args: "..."` appends arguments |

### Per-row settings in the CSV sheets

| Sheet | Columns that change the analysis |
|---|---|
| `samplesheet.csv` | `fastq_1`/`fastq_2` or `bam`, `antibody` (picks `csaw_defaults`, track/heatmap grouping), `condition`, `replicate`, `control_id` (the IgG of each target), `group`, `norm_run`, `exclude` |
| `norm_runs.csv` | `method` = the primary method (chipseqspikeinfree, spikein, tmm, cpm), `bin_size`, `cutoff`, `max_turns`, `min_first_turn` (optional), `ref_condition` (optional: reference condition of the MA plots), per normalisation run |
| `peaksets.csv` | `target_group`, `control_group`, `small_filt`, `large_filt`, `merge_bp`, `max_width`, `fdr`, `lfc` (blank → `csaw_defaults[antibody]` → `csaw_*`) |
| `classifications.csv` | `regions` (a peak set or BED), `ref_group`, `treat_group`, `norm_run`, `lfc_retained`, `lfc_lost`, `fdr` |
| `heatmaps.csv` | `regions`, `samples`, `labels`, `window_bp`, `bin_bp`, `sort_using`, `cluster_k`, `z_min`, `z_max`, `y_max` (number or auto), `colors`, `engine` (deeptools or enrichedheatmap), `cluster_samples`, `cluster_labels`, `legend_title`, `region_labels` |
| `tracks.csv` | `locus_or_gene`, `flank`, `samples`, `labels`, `colors`, `max_value` (shared, auto, a number, group:NAME), `beds`, `scalebar`, `bed_colors`, `bed_labels`, `log2ratio`, `ini_args` |
| `intersects.csv` | `bed_a`, `bed_b`, `label_a`, `label_b`, `min_overlap`, `reciprocal`, `output` (optional: `records` default or `overlap`) |

### How the three manuscript datasets use these (example configs in `HPC/examples/projects/`)

| Setting | NUCLEATION | CBX2 | EED |
|---|---|---|---|
| Normalisation (primary; TMM, CPM as diagnostics) | ChIPseqSpikeInFree (bin 10 kb; cutoff 1.2; max_turns 0.86 / 0.989) | ChIPseqSpikeInFree per mark (bins 5–30 kb); spike-in table kept as the alternative (`norm_runs_spikein.csv`) | ChIPseqSpikeInFree (K27me2 50 kb / 0.7 / 0.876; K27me3 15 kb / 1.2 / 0.9999) |
| bigWig bin | 50 bp | 50 bp | 50 bp |
| csaw | `csaw_norm_bins: all`, `csaw_max_frag: 500`; per-set thresholds in `peaksets.csv` | same; 10 peak sets | same; blank thresholds → `csaw_defaults` (K27me3 3/1.5/50, K27me2 2/1/50) |
| Classification | none | 1 µM vs DMSO, K27me3 and K27me2 | EZH2i 500 nM vs DMSO, K27me3 and K27me2 |
| Heatmaps / tracks | on / on | on / on (`tracks_limit_stat: max_bin`) | on / on (incl. the 13 legacy QC regions as a `tracks.csv` row) |
| Peak QC loci | not set: the selector's built-in Polycomb loci (HOXA, HOXD, CDKN2A) and `auto` negative | HOXA, CDKN2A; negative GAPDH | HOXA, CDKN2A; negative GAPDH |
| Annotation / intersects | on / off | on / on | off / off |

### Fixed in the step scripts (change only by replacing the script or with `<script>_args`)

Trimmomatic `ILLUMINACLIP` values `2:30:10:8:True` and no quality trimming; bowtie2
`--very-sensitive-local --no-mixed --no-discordant --dovetail`; the filter flags `-f 3 -F 2828`; Picard
`OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500`; bamCoverage `--skipNonCoveredRegions`; TSS enrichment window
±2 kb, 50 bp bins, ±250 bp centre; fingerprint (500 bp bins, 500,000 samples) and correlation (5 kb
bins, Spearman) settings; MACS3 `-f BAMPE --keep-dup all`; SEACR `norm`; csaw `pe="both"`; the
classification design `~group`; peak QC window tiers, overlap limits and region categories.

## Settings that apply to every step

- **Replacing a step's script:** put a copy in the project's `utils/` and name it in `config.yaml` under
  the script's base name (e.g. `callpeaks_csaw: utils/callpeaks_csaw.R`); extra arguments with
  `<script>_args` (declared in `nextflow.config`; changing one re-runs that step and what depends on it).
- **Samples:** `exclude: true` drops a row; `skip_groups` (default none) keeps the listed `group` values in
  the samplesheet but skips them in every core step.
- **Errors and retries:** exit 137/140/143 (out of memory/time) → retried with memory × attempt;
  exit 3 (bad input) → stop with a message, no retry; anything else retried once.
- **Publishing:** cluster `copy`; Docker versions `auto` (hard link to the work folder, else copy). In the Docker versions the trimmed FASTQs, unfiltered aligned BAMs, filtered BAMs, SEACR bedGraphs and BIN_COUNTS tables go to `intermediate_dir` only when `publish_intermediates` is true (default false; with false they stay in the work folder). The cluster always publishes them.
- **Docker resource budget** (`WINDOWS|MAC/conf/docker.config` + `docker_tiers.config`): every request
  is capped to the run's budget; per-step settings by budget tier (< 7 GB small, 7–12 GB mid, ≥ 12 GB
  large): CSAW 1 cpu (2 on mid/large), 4 GB; MARKDUP 1 cpu, 4/5/10 GB (Picard heap = memory − 2 GB) on HPC; in the Docker versions BAM_PROCESS (MARKDUP + DEDUP_SPLIT + BLACKLIST_FILTER + BAM_FINALIZE) 4 cpus, 6 GB × attempt (Picard heap = memory − 2 GB, so 4 GB at the first attempt);
  CHIPSEQSPIKEINFREE 2 cpus, 3.5 GB; BIN_COUNTS 2 cpus, 7 GB; TMM 3 GB; NORM_QC 4.5 GB; REPLICATE_CORRELATION(_PLOT) 3 GB; NORM_CPM and NORM_COMBINE 1 GB; CLASSIFY_REGIONS 1 cpu, 3 GB; BUILD_MATRIX 1 cpu, 4 GB; ALIGN all cpus, 4.5 GB; PEAK_ANNOTATION 4 GB. Other steps keep the cluster requests, capped — which makes most
  steps reserve the whole budget and run one at a time (seen in the 2026-10-26 pilot; to be tuned from
  the measured peaks).
- **Startup checks** (`prepareRun` in `main.nf`): samplesheet columns and IDs, one input type per row,
  controls exist, groups named in every CSV exist, reference files exist, thresholds are numbers, peak
  QC loci are well-formed; problems are listed together before any
  step starts.
