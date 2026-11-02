# ATAC-seq — what the analysis does

This file describes the analysis itself. It applies to the cluster version (`HPC\`), the Windows
version (`WINDOWS\`) and the Mac version (`MAC\`); only how each version is launched and where its
software comes from differ. Inputs, config keys and commands are in each version's `README.md`.

The downstream steps (pre-aligned BAMs, differential accessibility, peak annotation, motif enrichment and common motifs) are
described after the upstream steps, in their own sections. They are switched off by default; the switches, the input tables and
the status of the Docker versions are in each version's `README.md` ("Downstream analyses").

## Software versions to cite

The pipeline does not write a Methods section yet. Take versions from where the run actually ran:

| | Cluster (`HPC\`) | Docker (`WINDOWS\`, `MAC\`) |
|---|---|---|
| Image | `cutandrun-py_1.2.sif` (python, bowtie2, samtools, bedtools), `cutandrun-r_1.3.sif` (R) + public biocontainers | `ghcr.io/genesig/omics:2026.10.3` (tag in `docker\image.txt`; package locks in `PIPELINES\docker\omics\env\*.lock.txt`) |
| Nextflow | 26.04.6 | 26.04.6 (inside the image) |
| R / csaw / edgeR | 4.3.3 / 1.36.0 / 4.0.16 | 4.4.3 / 1.40.0 / 4.4.0 (Rsamtools 2.22.0; Bioconductor 3.20) |
| bowtie2 / samtools | 2.5.5 / 1.23 | 2.5.5 / 1.23.1 |
| Trimmomatic / FastQC / MultiQC | 0.40 (biocontainer; the legacy scripts used the 0.32 jar) / 0.12.1 / 1.35 | 0.40 / 0.12.1 / 1.35 |
| sambamba / bedtools / bedops | 1.0.1 / 2.31.1 / 2.4.42 | 1.0.1 / 2.31.1 / 2.4.42 |
| deepTools / MACS3 | 3.5.6 / 3.0.4 | 3.5.6 / 3.0.4 |
| Python / pysam | from `cutandrun-py_1.2.sif`; no pysam listed (the TN5 step uses deepTools' pysam) | 3.12 / 0.24.1 |
| Downstream packages: ChIPseeker / monaLisa / txdbmaker / JASPAR2024 | not recorded in this file (the cluster's downstream image is `omics-r_1.0.sif`) | 1.42.0 / 1.12.0 / 1.2.0 / 0.99.6 (from the image's R lock file) |

The legacy analyses that the downstream port was checked against ran csaw 1.42, edgeR 4.6.3 and Bioconductor 3.21 (not swapped
in the port). The Docker image has csaw 1.40 and edgeR 4.4, so the numbers of the two can differ.

The Docker column comes from the image's package locks (`ngs`, `r` and `py` environments). Anything
not in the locks or in the table is not stated here: check the image lock files
(`PIPELINES\docker\omics\env\`) before you report it. The cluster's Bioconductor release is not recorded here.

The R and Bioconductor packages of the Docker image are newer than those of the cluster's R image, so
numbers can differ slightly between the two. Say which version produced the numbers you report.

## What each step does

| Stage | Step | Tool and key settings |
|---|---|---|
| 1 | Read QC (`run_fastqc`) | FastQC 0.12.1 on each raw FASTQ (`fastqc -o . --extract`) |
| 1 | Trimming (`trim`) | Trimmomatic 0.40, paired-end: `ILLUMINACLIP:NexteraPE-PE.fa:2:30:10:8:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:25` (the defaults of `trim_adapter_set`, `trim_illuminaclip`, `trim_steps`, `trim_minlen`) |
| 2 | Pair count (`align_chunk_pairs` > 0) | COUNT_PAIRS: read pairs of each library; R1 and R2 must have the same count, otherwise the run stops |
| 2 | Alignment | bowtie2 `-p N --local --very-sensitive-local -q --no-mixed --no-discordant --dovetail --reorder -I 10 -X 2000 -x <index>`, piped into `samtools view -b`; one task per chunk of `align_chunk_pairs` pairs (`-s`/`-u`), or one per library when `align_chunk_pairs` is 0 |
| 2 | Chunk merge | ALIGN_MERGE: `samtools cat` of the chunk BAMs in chunk order; the bowtie2 logs are summed |
| 2 | Alignment rate | ALIGNMENT_SUMMARY: overall rate from the bowtie2 log; warns below `alignment_min_pct` (20); never fails |
| 3 | Filter | BAM_FILTER: `samtools sort`, raw flagstat and idxstats, then `samtools view -b -q <mapq=30> -f 3 -F 2828` (properly paired; excludes unmapped, mate unmapped, secondary, QC-fail and supplementary records); with `remove_chroms` (chrM) the listed contigs are dropped (`-L`) |
| 3 | Chromosome fraction | QC_CHROM_FRACTION: percent of mapped reads on each removed contig (only when `remove_chroms` is set) |
| 3 | Duplicates | DEDUP: `sambamba markdup -r` (duplicates removed) when `remove_duplicates` is true; marked and kept when false |
| 3 | Blacklist | BLACKLIST_FILTER (when `blacklist` is set): name sort, then `bedtools pairtobed -type neither`: a pair is removed when either mate overlaps the blacklist |
| 3 | Final BAM | BAM_FINALIZE: coordinate sort, index, final flagstat and idxstats, and an orphan check (every read name must occur exactly twice among the paired reads; otherwise the run stops with exit 3) |
| 3 | Contamination (when `contamination_index` is set) | CONTAMINATION: unmapped pairs (`samtools view -u -f 12 -F 256`) converted to FASTQ, then aligned with bowtie2 `--no-unal` to the contaminant index; counts per reference organism |
| 4 | Tn5 shift (`tn5_shift`, off by default) | TN5_SHIFT: deepTools `alignmentSieve --ATACshift` (+4 bp on the forward strand, -5 bp on the reverse), then sort and index (pysam). When on, every later step reads this BAM instead of the final BAM. How the shift treats the proper-pair flags is under statistician review and is not yet verified |
| 4 | Fragment sizes | QC_FRAGMENT_SIZES_TABLE (Rsamtools: fragment-size distribution and summary per BAM) and QC_FRAGMENT_SIZES_PLOT (density and boxplot SVGs) |
| 5 | TMM scale factors | NORM_TMM, one task per `tmm_group` (see "Normalisation"): csaw counts in 10 kb bins, edgeR TMM factors, MA plots |
| 5 | bigWig tracks | BIGWIG, for each sample with a `tmm_group`: deepTools `bamCoverage --binSize 10 --skipNonCoveredRegions --normalizeUsing <None, CPM or RPKM> --scaleFactor <multiplier> --extendReads -p N`, plus `--effectiveGenomeSize` when `effective_genome_size` is set |
| 5 | TSS enrichment | TSS_ENRICHMENT (deepTools `computeMatrix reference-point --referencePoint TSS`, flank `tss_flank` = 2000 bp each side, bin 50 bp, `--skipZeros --missingDataAsZero`, then `plotProfile`), and QC_TSS_TABLE; only for samples with a bigWig |
| 6 | Peaks | MACS3 3.0.4 per `peak_group`: `callpeak -f BAMPE -g <effective_genome_size> --keep-dup all -q 0.05 --min-length 250 --max-gap 500 --fe-cutoff 2.0`, replicates pooled as one target, no control (see "Peak calling") |
| 6 | Merged peak sets | PEAK_MERGE (with `peak_merges`): sort, then `bedops --merge` of the listed groups' narrowPeak files |
| 6 | Peak statistics and FRiP | QC_PEAK_STATS (peak counts, widths, chromosome summaries); FRIP and QC_FRIP_TABLE: fraction of each sample's mapped alignments that overlap at least one peak of its own peak group |
| 7 | QC report | MULTIQC over the FastQC, Trimmomatic, bowtie2, sambamba, flagstat, idxstats and MACS3 files |
| 8 | Pre-aligned BAM (`bam` column) | BAM_PREALIGNED: the BAM is copied (sorted only when its header is not coordinate-sorted), indexed, with flagstat and idxstats; it is taken as final (see "Pre-aligned BAMs") |
| 8 | TxDb (`run_differential` with promoters, or `run_annotation`) | TXDB_BUILD: `txdbmaker::makeTxDbFromGFF` on the `gtf`, saved as one SQLite file (see "Peak annotation") |
| 9 | Differential accessibility (`run_differential`) | DIFF_ACCESS: csaw counts and edgeR tests for one comparison, over peaks and, optionally, promoters (see "Differential accessibility") |
| 10 | Peak annotation (`run_annotation`) | DIFF_ANNOTATE: ChIPseeker `annotatePeak` of each comparison's table, four sets (see "Peak annotation") |
| 11 | Motif enrichment (`run_motifs`) | MOTIF_ENRICHMENT: monaLisa `calcBinnedMotifEnrR` per comparison and region set, and per motif target (see "Motif enrichment") |
| 11 | Common motifs (`run_motifs`) | MOTIF_COMMON: motifs significant in several comparisons, per region set (see "Common motifs") |

The order of the duplicate and blacklist steps is the legacy order (duplicates, then blacklist).
The blacklist removes both mates of a pair when either overlaps, and keeps MAPQ and the sequences
(the legacy script converted the BAM to BEDPE and back, which lost them).

## Filtering and duplicates

- **MAPQ and flags.** `-q 30` keeps alignments with MAPQ 30 or more (`mapq` in `config.yaml`; `0`
  reproduces the legacy behaviour, where the literal `-q mapq` kept every alignment). `-f 3` keeps
  only properly paired reads, and `-F 2828` removes unmapped, mate unmapped, secondary, QC-fail and
  supplementary records.
- **chrM.** The mitochondrial contig (`remove_chroms: chrM`) is dropped in BAM_FILTER, before
  duplicates are marked, so it does not enter the library sizes or the TMM factors.
- **Duplicates.** `sambamba markdup -r` removes duplicates by default (`remove_duplicates: true`).
  With the ordering of the reads fixed by `--reorder`, the choice of which tied duplicate is kept
  does not depend on the thread count or on the chunking.
- **Orphans.** BAM_FINALIZE checks that every read name occurs exactly twice among the paired reads
  (`S.orphan_check.txt`). A failure stops the run, because a damaged pair would bias the counts.

## Alignment in chunks

With `align_chunk_pairs` above 0 (default 10,000,000), each library is counted (COUNT_PAIRS) and
aligned in chunks of that many read pairs. Each chunk is aligned with bowtie2's `-s` (skip) and `-u`
(number of pairs) options, so no temporary FASTQ files are made. ALIGN_MERGE then joins the chunk
BAMs in chunk order. The merged BAM has the same records, in the same order, as one unchunked run
(`align_chunk_pairs: 0`); the outputs do not depend on the chunk size. This was checked on 1-million-pair
subsamples of both datasets (aligned BAM records identical, and identical peaks, bigWigs, QC tables
and scale factors).

## Normalisation (TMM)

NORM_TMM runs once per `tmm_group`, on the final BAMs of the group's libraries (or the Tn5-shifted
BAMs when `tn5_shift` is on):

- **Counts.** csaw counts the fragments of each library in bins of `tmm_bin_size` (10 kb), with
  `minq = mapq`, paired-end, `filter = 10` (bins with fewer than 10 fragments over all libraries are
  dropped), and `dedup = FALSE` (the BAMs are already duplicate-removed).
- **Fragment length cut-off.** `tmm_max_frag: auto` uses the largest proper-pair fragment of the
  group (the true maximum over all its BAMs); an integer is used as given.
- **Factors.** edgeR's TMM normalisation factors (csaw `normFactors`, unweighted; geometric mean 1).
  The library size is the number of fragments counted by csaw, not the sum over the retained bins.
- **Multiplier for bamCoverage.** With `bigwig_units: none` (the default) the multiplier is
  `1e6 / (nf x lib.size)`, passed as `--scaleFactor` without `--normalizeUsing`. With `cpm` or `rpkm`
  it is `1 / nf`, passed with `--normalizeUsing CPM` or `RPKM`.
- **QC flag.** A library with fewer than 100,000 fragments, or fewer than 1,000 retained bins, is
  marked `qc_pass = FALSE` in the scale-factor table.

The TMM method assumes that most bins are not differentially accessible between the libraries of a
group. Where a treatment changes accessibility across most of the genome, TMM removes that global
difference by construction. The factors are therefore for display and for groups whose members are
expected to share most of the genome; they must not be used to claim or deny a global change in
accessibility. The MA plots (`<group>.ma_plot.png`, `.pdf`) show whether the assumption holds: the
red line should pass through the centre of the dense background cloud. The group composition
(`tmm_group_column`) matters for this: exact groups keep, for example, K27M and K27MKO apart, while
the legacy membership column pools them. Changes to factors, groups or thresholds are scientific
decisions to confirm with the statistician first.

## Peak calling

MACS3 runs once per `peak_group`. The fragments of all libraries in the group are pooled as one
target (`-t` with every BAM; `-f BAMPE` uses the fragment from each mate pair), with `--keep-dup all`
because the duplicates were already removed, and no control sample. The output prefix records the
settings: `<group>.qval<q>[.minlen<l>][.maxgap<g>][.fe<fe>][.broad]`, and a part is left out when its
setting is off.

- `macs3_qvalue` (0.05) is the `-q` cut-off; `macs3_min_length` (250) and `macs3_max_gap` (500) are
  the peak length and gap limits; `macs3_fe_cutoff` (2.0) is the fold-enrichment cut-off (a value of 1
  or less turns it off).
- `macs3_call_summits` writes summit BEDs; `macs3_broad` switches to broad peaks with
  `macs3_broad_cutoff` (0.1). The default is narrow peaks.
- `effective_genome_size` is MACS3's `-g`. It is required when a sample has a `peak_group`.

FRiP is the number of each sample's mapped alignments that overlap at least one peak of its own peak
group (`bedtools intersect -u`), divided by its mapped alignments in the final flagstat. Samples with
an empty peak set get a FRiP of 0 with a warning.

## Tn5 shift (off by default)

`tn5_shift: true` shifts each fragment for the Tn5 insertion offset with deepTools `alignmentSieve
--ATACshift` before the fragment sizes, the TMM factors, the bigWigs, MACS3 and FRiP are computed. It is
not in the legacy analysis and is off in both lab datasets. The treatment of the proper-pair flags
under this shift is still under review (`bin/tn5_shift.sh`), so do not switch it on for reported
results without confirming it.

## TSS enrichment

TSS_ENRICHMENT takes each sample's bigWig and the `tss_bed` file, and computes the mean signal in
bins of 50 bp over 2,000 bp each side of every TSS. The enrichment is the highest mean within 250 bp
of the TSS over the mean of the outer 200 bp on each side. The legacy score (`max_over_min`, the
profile maximum over its minimum) is also written. QC_TSS_TABLE merges the per-sample values with the
samplesheet labels.

## Pre-aligned BAMs (`bam` column)

A samplesheet row can give a final BAM instead of FASTQs (`bam` set, `fastq_1` and `fastq_2` empty). BAM_PREALIGNED copies it
to `S.final.bam` (it is sorted first only when its header is not coordinate-sorted), indexes it, and writes its flagstat and
idxstats. The BAM is taken as **final**: filtered, duplicate-removed and blacklist-filtered by the user. No MAPQ filter, duplicate
removal, blacklist removal, orphan check or raw statistics are run on it. From there it is used like any final BAM (TMM, bigWig,
MACS3, FRiP, and the downstream steps). Whether a pre-aligned BAM meets the same filters as the pipeline's own is the user's
question to answer; the pipeline does not check it.

## Differential accessibility (`run_differential`)

Question: which accessibility regions have a different number of fragments between the case libraries and the control libraries
of one comparison. The legacy csaw and edgeR analysis (`differential_accessibility.Rmd`, plus the promoter part of an older copy)
is ported as `diff_accessibility.R`, one task per comparison.

- **Regions.** The peaks of both sides, each side's peak groups (`case_peaks`, `control_peaks`; default: the side's
  `peak_group` values), merged by union. With `diff_regions: peaks,promoters`, also the promoter windows of the TxDb genes:
  `diff_promoter_up` (2000 bp) upstream and `diff_promoter_down` (200 bp) downstream of the TSS. Promoters whose background
  abundance is not above log2(`diff_promoter_min_fc`) = log2(4) are removed (csaw `filterWindowsGlobal`).
- **Counts.** csaw `regionCounts` on paired-end fragments with MAPQ at least `diff_minq` (10) and fragment length at most
  `diff_max_frag` (500). Duplicates are dropped when `diff_dedup` is TRUE; the default (`null`) is FALSE when `remove_duplicates`
  is on, since the BAMs are then already duplicate-free. The 10 kb bins (`diff_bin_size`) are counted once, for the whole genome:
  their library sizes give the normalisation factors (`normFactors`, TMM, unweighted), as in NORM_TMM.
- **Model.** Without a `block_column`: `~ 0 + grouping` (Control, Condition), the test is Condition minus Control (the legacy
  per-cell-line model). With a `block_column` that has at least two values: `~ block + grouping`, and the test is the Condition
  coefficient (the legacy gd2 model, block = cell line).
- **Test.** The quasi-likelihood F-test (`glmQLFit`, robust, and `glmQLFTest`) when the comparison has at least `diff_qlf_min_libs`
  (3) libraries. Otherwise `glmFit` with the dispersion fixed at `diff_fixed_dispersion` (0.05, or the comparison's
  `fixed_dispersion`) and a likelihood-ratio test (`glmLRT`).
- **Statistics per region.** `overlapResults` combines the bin-level tests; the statistics come from its `best` table
  (`diff_table: best`; AM legacy: `combined`). Each region has PValue, FDR (edgeR's Benjamini-Hochberg over this comparison's
  regions only), logFC, and score = logFC.
- **Outputs.** The table, the BEDs of all regions and of the regions with logFC above `diff_lfc` (0.585) or below minus that
  value. The up and down BEDs have **no FDR filter**: they are fold-change sets. The counts at `diff_fdr` (0.1) are printed in
  the summary only; no table is filtered by FDR.
- **Legacy values that differ.** Older AV and AM values are toggles: `diff_pe: none` (older AV and AM: paired data counted as
  single-end), `diff_lfc: 1` (older). The `max.frag` of 500 is csaw's default; the legacy never set it.

**Limits (state them when reporting).** The legacy method, ported unchanged:
- With fewer libraries than `diff_qlf_min_libs` (for example one library per side) there is no replication, so the variance is
  assumed, not estimated. PValue and FDR are then not valid error rates, and only the ranking by logFC is informative.
- A blocked comparison (gd2: the cell line as block) treats the cell line as the unit of replication and assumes one common
  high-versus-low effect in all lines. A block that holds libraries of one condition only carries no information.
- The tested regions come from the data (peaks called on pooled replicates, the union of both sides), so the test uses the same
  reads as the selection.
- TMM corrects background bias only; efficiency bias (a different signal-to-noise ratio between libraries, FRiP or TSS enrichment)
  moves the logFC of all regions together. Check the MA plot (`NAME.diag_ma.svg`) and the factors (`NAME.norm_factors.tsv`).
- Promoter windows of neighbouring genes overlap, so promoter rows are not independent.

**Checked so far (the cluster port, same image).** Against the legacy result files, the region sets were identical for the five
per-line comparisons and the gd2 comparison. The largest relative difference of PValue, FDR and logFC was of the order of 1e-8,
and about 2e-7 on the logFC of gd2 (where |logFC| is about 0.05). The port attributes these to the csaw and edgeR version
difference (1.40 and 4.4 here, 1.42 and 4.6.3 in the legacy); that attribution is not verified, because the versions were not
swapped. Run times and memory on a Docker computer are not measured.

## Peak annotation (`run_annotation`)

`diff_annotate.R` (ChIPseeker) annotates the regions of one differential table, usually a comparison's peaks table, with the
gene models of the TxDb that `txdb_build.R` makes from the `gtf` (the same GTF as the peaks and the promoters). It runs
`annotatePeak` with `annotation_tss_region` (`-2000,200`) and the gene-symbol package `org.Hs.eg.db` (`org.Mm.eg.db` for mouse).
By default (`annotation_fdr` null) **all** regions of the table are annotated; the older AV analysis kept FDR < 0.1 first
(`annotation_fdr: 0.1`).

The four sets: `all_sites`; `promoter_sites` (annotation contains "Promoter"); `non-promoter_sites`; `distal_sites` (annotation
contains "Distal Intergenic"). Each has a CSV and a BED. The logFC histogram (`NAME_peak_distribution.svg`, Helvetica) is drawn
inside its device; the legacy PDF was empty because `hist()` ran outside `pdf()`. The chromosome lengths come from the BSgenome of
the genome build (`BSgenome.Hsapiens.UCSC.hs1` for `hs1`). The legacy's final whole-genome promoter BED chunk never ran and is not
ported.

**Limits.** The annotation is the nearest TSS gene with no distance limit: a distal intergenic region far from a gene still carries
that gene. A region that overlaps a promoter and another feature takes ChIPseeker's priority order. All-region annotation describes
the tested region set, not the differentially accessible regions; use the FDR and logFC columns for that.

**Checked so far.** Peak annotation of the 6349 comparison was identical to the legacy (78,679 annotated rows) in the cluster port.

## Motif enrichment (`run_motifs`)

`motif_enrichment.R` (monaLisa) runs one of three analyses, all with the same selection code:

- **Diff mode (per comparison and `motif_regions` set).** The regions of a differential table (`promoters` = the promoter
  table, `peaks` = the peak table, or an annotation set) are split into two bins: Increased (logFC above `motif_lfc`, 0.585) and
  Decreased (logFC below minus that). Regions are resized to the mean width around their centre. Sequences come from the BSgenome.
  The motifs are the JASPAR2024 PWMs (or a cached `motif_pwm_rds`); a PWM hit needs a score of at least `motif_min_score` (10).
  `calcBinnedMotifEnrR` runs with the binomial test, the other bin as background (`otherBins`) and `maxKmerSize` 3, which matches
  the background to the foreground's GC and di- and trinucleotide composition. Adjusted P values are Benjamini-Hochberg over all
  motifs and both bins. A motif is selected at -log10 adjusted P above `motif_padj_primary` (2, adjusted P < 0.01), or above
  `motif_padj_relaxed` (1, adjusted P < 0.1) when two or fewer pass the primary threshold. Heatmaps show the top `motif_top` (40).
  The analysis is skipped (only `NAME.skipped.txt`) when a bin has fewer than `motif_min_bin_peaks` (10) regions.
- **Gene target (`type: gene` in `motif_targets.csv`).** The regions of the `all_sites` table of the universe comparison whose
  nearest gene matches the target symbol (`motif_gene_column` = `geneId`; `exact` (default, the whole symbol equals the nearest gene) or `substring` (legacy `grepl`)) against
  all other regions. Each bin needs at least `min_peaks` (default `motif_target_min_peaks`, 15).
- **Regions target (`type: regions`).** The regions of that table that overlap a BED of target regions (in the pipeline genome's
  coordinates; a liftover is a project pre-step) against the other regions. Same minimum.

The selection, the heatmaps and the sequence composition plots are the legacy's; `motif_cpus: 1` runs the serial BiocParallel
parameter with `motif_seed` (123), which is repeatable. A value above 1 uses several cores and is not bit-reproducible.

**Limits.** The bins are cut on the estimated logFC only, with no FDR filter, so for an analysis without replication the extremes
of logFC are partly noise. The reference is the other bin, so every statement is relative. monaLisa's binomial test is one-sided,
and the two columns of a two-bin analysis mirror each other. Each PWM hit counts as one independent trial, the hits are clustered,
and the JASPAR matrices are correlated: the p-values are anti-conservative and the number of significant motifs counts redundant
matrices. Read the ranking and log2enr, not the count. In the gene mode, `substring` also matches other genes (FOS: FOSB, FOSL1, FOSL2; MYC: MYCN, MYCL, MYCBP2; SOX2: SOX2-OT, SOX21).

**Checked so far.** The motif tables were identical to the legacy for the five lines (largest relative difference 1.5e-14). The
gene and regions targets were identical to the legacy chunks on a 5,022-row subsample; the full-table runs ran out of memory at
13 GB in Docker and need the cluster. The JASPAR2024 PWMs of the package were not checked against the cached legacy set.

**Genome build.** For `hs1` the sequences and the chromosome lengths come from `BSgenome.Hsapiens.UCSC.hs1`, and for `mm39` from `BSgenome.Mmusculus.UCSC.mm39`; both are in the Docker image.

## Common motifs (`run_motifs`)

`motif_common.R` (the legacy `find_common_motifs.Rmd`) takes the diff-mode tables of one region set, for the comparisons that are
not `pooled`. A motif is significant in a comparison when its adjusted P for the Increased bin is at most `motif_common_padj`
(0.10). A motif is common when it is significant in at least `motif_common_min_lines` (2) comparisons. The labels are the
comparison names with `motif_common_label_strip` removed (AV example: `_high_vs_low.*`). The output is the table of common motifs,
a Venn diagram (written only when a motif is significant in all inputs), and the stacked summaries.

**Limits.** The count of comparisons is a set-overlap descriptor, not a test of common enrichment, and it inherits the limits of the
motif analysis above. The key is `motif.id` (the JASPAR matrix ID); the output has `motif.name` and a `label` (the name, or name (ID) when several matrices share the name). On the five legacy lines the legacy Venn diagram cannot be
reproduced, because no motif is significant in all five. Comparisons that rest on one library per side should not be mixed with
replicated ones without saying so.

## What the pipeline does not do

- The `match_rna_atac_data` step of the legacy scripts (matching RNA-seq and ATAC-seq data) is out of scope and not ported.
- Choose contrasts, thresholds, normalisation or peak settings for you: they come from the CSVs and
  `config.yaml`, which an agent fills in only from what the user confirmed. This includes the comparisons, the `block_column`,
  the `pooled` and `fixed_dispersion` values, and the motif targets.

## Reproducibility

- Every step's exact command is published next to its output (`.../log/*_command.sh`).
- `metadata/nf_logs/` holds the Nextflow logs; `metadata/pipeline_info/` the execution trace and
  report; the Docker versions also write `metadata/host_paths.tsv` when Windows paths were translated.
- To rerun only some steps (for example the plots), change the setting or CSV and run again: the work folder is kept (`keep_work: true`, only the launcher's `clean` command deletes it), so Nextflow resumes and recomputes only the steps downstream of the change.
- Re-running with the same inputs resumes from the work folder (in the Docker versions a per-project
  Docker volume, kept by default, or `metadata/work` with `work_volume: false`; on the cluster under
  `/project/rephilllab/nf_work/<user>/<project>`).
- The scale-factor table records the settings of each group in its `note` column (`lib_size`, `max.frag`,
  `mapq`, `bin size`, `bigwig_units`).

## Output layout

```
results/
  qc/fastqc/raw, qc/fastqc/trimmed     FastQC reports per FASTQ
  qc/trimmomatic/                      Trimmomatic logs
  qc/bowtie2/                          bowtie2 logs per library, and the command logs
  qc/alignment/                        alignment_summary.tsv
  qc/chrom_fraction/                   chrom_fraction.tsv (with remove_chroms)
  qc/contamination/                    contamination table and summary (with contamination_index)
  qc/fragment_sizes/                   fragment size tables and SVG plots
  qc/tss_enrichment/                   TSS matrices, profiles and the TSS table
  qc/peaks/                            peak_stats.tsv
  qc/frip/                             per-sample and merged FRiP tables
  qc/multiqc/                          MultiQC report
  bam/stats/                           flagstat, idxstats, orphan check, sambamba log
  bam/final/                           final BAMs (+ .bai) used by every later step
  bam/tn5/                             Tn5-shifted BAMs (with tn5_shift)
  bigwig/                              <sample>.<method>.<units>.bw
  scale_factors/                       <group>.scale_factors.tsv, MA plots, <group>.max_frag.txt
  peaks/macs3/                         MACS3 peak files, _peaks.xls, summits and gapped peaks (optional), logs
  peaks/merged/                        <name>.merged.bed (with peak_merges)
  annotation/txdb/                     txdb.sqlite (with run_annotation or a promoter table)
  differential/                        <name>.csv, <name>_promoters.csv, norm factors, summary, replicates, diagnostics (run_differential)
  differential/bed/                    <name>.bed, _upregulated, _downregulated, and the promoter BEDs
  annotation/                          <name>_<set>_sites.csv and .bed, peak distribution plot, annotation summary (run_annotation)
  motif/<name>/<region>/               motif analysis tables, summary, heatmaps, similarity, GC and dinucleotide plots (run_motifs)
  motif/<name>/                        the same for a motif target; <name>.skipped.txt when the analysis was skipped
  motif/common/<region>/               common_significant_motifs.csv, Venn diagram, motif_analysis_summary.csv
intermediate_files/                    trimmed FASTQs, aligned and filtered BAMs; in the Docker versions only
                                       with publish_intermediates: true (default false: the folder stays empty)
```
