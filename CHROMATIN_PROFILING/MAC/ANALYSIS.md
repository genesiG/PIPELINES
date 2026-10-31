# CUT&RUN / chromatin profiling — what the analysis does

This file describes the analysis itself. It applies to all three versions of the pipeline (cluster
`HPC/`, Windows `WINDOWS/`, Mac `MAC/`); only how each version is launched and where its software
comes from differ. Inputs, config keys and commands are in each version's `README.md`.

## Software versions to cite

The pipeline does not write a Methods section yet. Take versions from where the run actually ran:

| | Cluster (`HPC/`) | Docker (`WINDOWS/`, `MAC/`) |
|---|---|---|
| Image | `cutandrun-py_1.2.sif`, `cutandrun-r_1.3.sif` + public biocontainers (`HPC/containers/`) | `ghcr.io/genesig/omics:<tag>` (tag in `docker/image.txt`; package locks in `PIPELINES/docker/omics/env/*.lock.txt`) |
| R / Bioconductor | 4.3.3 / 3.18 (csaw 1.36.0, edgeR 4.0.16) | 4.4.3 / 3.20 (csaw 1.40.0, edgeR 4.4.0) |
| ChIPseqSpikeInFree | 1.2.4 (GitHub commit f9cc88f) | 1.2.4 (same commit) |
| bowtie2 / samtools / Picard | 2.5.5 / 1.23 / 2.27.5 | 2.5.5 / 1.23.1 / 2.27.5 |
| Trimmomatic / deepTools / MACS3 / bedtools | 0.40 / 3.5.6 / 3.0.4 / 2.31.1 | same |
| Nextflow | 26.04.6 | 26.04.6 (inside the image) |

The manuscript's peak sets were made with the cluster versions. The Docker versions use newer
R/Bioconductor packages (accepted by the lab); peak and class counts may differ slightly. Say which
versions produced the numbers you report.

## What each step does

| Stage | Step | Tool and key settings |
|---|---|---|
| 1 | Read QC | FastQC on the raw FASTQs |
| 1 | Trimming (`trim`) | Trimmomatic PE, `ILLUMINACLIP:<TruSeq3-PE-2 (CUT&RUN) or NexteraPE-PE (CUT&Tag)>.fa:2:30:10:8:True MINLEN:<trim_minlen=25>` |
| 1 | Alignment | bowtie2 `--local --very-sensitive-local --no-mixed --no-discordant --dovetail` to the combined human + spike-in index; `-I`/`-X` only when `bowtie2_minins`/`bowtie2_maxins` are set (default: bowtie2's own 0 / 500, as in the manuscript) |
| 1 | Filtering | `samtools view -q <mapq=30> -f 3 -F 2828`: properly paired, primary, both mates mapped |
| 1 | Duplicates | Picard MarkDuplicates (marks); duplicates removed from target BAMs when `remove_duplicates` (always for IgG) and from spike-in BAMs when `spikein_dedup` |
| 1 | Spike-in split | reads on contigs starting with `spike_prefix` (E. coli) go to a separate spike-in BAM; `mito_chrom` reads are dropped |
| 1 | Blacklist | `bedtools pairtobed -type neither`: a pair is removed when either mate overlaps the blacklist |
| 1 | Normalisation | per `norm_run`: primary method plus ChIPseqSpikeInFree, TMM, CPM (and spike-in when toggled) for comparison (see below); writes `scale_factors/<norm_run>.scale_factors.tsv` (primary) and `<norm_run>.scale_factors_all.tsv` (all methods) |
| 1 | bigWigs | deepTools bamCoverage `--binSize <bigwig_bin_size> --extendReads --skipNonCoveredRegions --normalizeUsing RPKM`, scaled by the run's bamCoverage-compatible factor; IgG bigWigs unscaled, same units |
| 1 | QC | alignment table, fragment sizes, library complexity (NRF/PBC1/PBC2 from fragments + Picard EstimateLibraryComplexity), preseq, fingerprint, sample correlation, replicate concordance, TSS enrichment (anchored on GTF transcript starts), FRiP, MultiQC, HTML report |
| 2 | Peaks | csaw, MACS3 (with the control) and SEACR (with the control, consensus over replicates), all on by default; broad or narrow per mark (see "Peak calling with MACS3 and SEACR") |
| 2 | Classification | retained / lost / unclassified regions between a reference and a treated group (`classifications.csv`) |
| 2 | Heatmaps | deepTools computeMatrix (reference point = centre, 100 bp bins, missing as zero) + plotHeatmap, or EnrichedHeatmap; optional seeded k-means clusters |
| 2 | Tracks | pyGenomeTracks panels of `tracks.csv` loci, plus the automatic **peak QC** set (below) |
| 2 | Annotation / intersections | ChIPseeker annotation of peak sets; bedtools intersections (`intersects.csv`) |

## Normalisation

Each sample names a `norm_run` (samplesheet); `norm_runs.csv` gives the run's **primary** method (`method`: chipseqspikeinfree,
spikein, tmm or cpm). Only the primary scales the bigWigs and the classification. For comparison, every run with at least two
libraries also computes ChIPseqSpikeInFree, TMM and CPM factors (`norm_methods`), and physical spike-in factors when
`norm_spikein: true`. All of them are in `scale_factors/<run>.scale_factors_all.tsv`; the primary alone stays in
`<run>.scale_factors.tsv`. The methods assume different things:

- **ChIPseqSpikeInFree** (the lab's choice for the manuscript datasets): a factor from the shape of the cumulative bin-count curve
  (`bin_size`, `cutoff`, `max_turns`, `min_first_turn`). It allows a global change of the mark. It is fragile at low depth or poor
  enrichment; libraries failing its QC get the run's largest factor (flagged). Factors have 2 decimals.
- **Spike-in**: factor from the E. coli spike-in reads per library. It allows a global change and needs the same spike-in amount
  per sample and enough spike-in reads (fewer than 1,000 reads are flagged).
- **TMM** (csaw/edgeR, `tmm_bin_size` bins): assumes most bins do not change. It removes a global change by construction.
- **CPM**: library size only (fragments of the final BAM). It assumes the total amount of the mark per library is the same.

**Normalisation QC** (`norm_qc/`) draws MA plots on 10 kb bins outside the blacklist: y is the log2 CPM difference of a library and
its reference (negative = less signal), x the mean log2 CPM. Panels are each pair of replicates of a group (where no biological
change is expected) and each treated library against the reference-condition library with the same replicate label. The dashed
line of each method is where that method declares "no change": log2 of the ratio of the two libraries' `edger_norm_factor`
(0 for CPM). The lines come from the factors, not from the plot. For the spike-in line the tables also give `line_se` (`norm_qc_pairs.tsv`), a Poisson lower bound on the SE (counting noise of the spike-in reads only). A line on the dense cloud means the method keeps the bulk of the
bins unchanged. That is automatic for TMM and CPM and is not evidence that they are right when the treatment changes the mark
globally. Between replicates every reliable method should lie on the cloud. The primary method is chosen with the user before the
results are seen, and recorded in the project's `CONTEXT.md`.

The reference condition of a run is the optional `ref_condition` column of `norm_runs.csv`; when blank, the condition of the
`ref_group` of the classifications that use the run (if exactly one), else the condition of the run's first library. Lines of a
factor from a QC-failed or imputed library are dotted. Panels are split over pages of `norm_qc_panels_per_page` (12): `<run>.ma.svg`,
`<run>.ma.p2.svg`, ... The bigWigs divide by the final-BAM read count and the edgeR steps by the csaw total (fragments up to
`max.frag`); the small difference per pair is the `bigwig_offset` column of `norm_qc_pairs.tsv`.

## Replicate concordance

After the BAMs are final, reads are counted once in 10 kb bins (csaw, the same counts as the normalisation QC and TMM). For every
pair of libraries of the same samplesheet `group`, `qc/replicate_correlation/replicate_correlation.tsv` gives:

- the Pearson correlation of log2 CPM (prior 2) and the Spearman correlation of counts, over non-blacklisted bins with reads;
- the same Pearson on the 10 % of bins with the highest mean signal;
- the robust spread of the replicate log2 ratio on those bins;
- for context, the median correlation of the two libraries with the other groups of the same antibody.

These are descriptive values without pass/fail thresholds. They depend on depth, bin size and the breadth of the mark, so compare
them within a dataset. The deepTools heatmap of all samples (`qc/correlation/`, 5 kb bins, Spearman) is kept
(`run_deeptools_correlation`).

## Peak calling with csaw

Reads are counted (`readParam(minq = mapq, pe = "both", dedup = remove_duplicates, max.frag =
csaw_max_frag, else bowtie2_maxins, else 500)`) into small windows (`csaw_window` 150 bp every
`csaw_spacing` 50 bp) and large windows (`csaw_large_window` 2,000 bp every `csaw_large_spacing`
500 bp), and into `csaw_bin` (10 kb) bins for normalisation (TMM over `csaw_norm_bins`: `background`
bins, the pipeline default, or `all` bins, which reproduces the manuscript). Windows are filtered by
enrichment over the global background (`small_filt` / `large_filt`, per peak set, else
`csaw_defaults[antibody]`, else the global keys). Target vs IgG is tested with edgeR
quasi-likelihood F-tests (likelihood-ratio tests without replicates); small- and large-window
results are merged (`merge_bp`) and peaks are kept at FDR ≤ `csaw_fdr` and logFC > `csaw_lfc`
(`csaw_lfc_from: combined` = the logFC of `combineTests`). Thresholds per peak set are in
`peaksets.csv` (`HPC/examples/projects/*/peaksets.csv` show the manuscript settings).

## Peak calling with MACS3 and SEACR

`run_csaw`, `run_macs3` and `run_seacr` are on by default (`run_macs3` / `run_seacr` are `null` in `nextflow.config`: on when every target has a control). MACS3 and SEACR need an IgG control per target (`control_id`); `call_peaks: false` turns all of them off (and FRiP, peak QC and peak annotation); bigWigs and normalisation are unaffected.
when a target has none and the user did not set `run_macs3` / `run_seacr` explicitly, they are skipped with a warning (an explicit
`true` stops the run instead). csaw is not affected.

Broad or narrow is decided per group from its antibody, against the list `broad_marks` (default `H3K27me3, H3K27me2, H3K36me3,
H3K4me1, H3K79me2, H3K79me3, H3K9me1, H3K9me2, H3K9me3, H4K20me1, H3F3A`: the ENCODE histone ChIP-seq broad-mark list, plus
H3K27me2 by lab request and H3K9me3, an ENCODE "exception" that St. Jude treats as broad). Matching ignores case and the histone
prefix (`K27me3` = `H3K27me3`); any other antibody is narrow.

- **MACS3** (`macs3_broad`: `auto` | `true` | `false`, default `auto`): broad marks use `--broad --broad-cutoff` `macs3_broad_cutoff`
  (0.1, the nf-core/cutandrun default); narrow marks give a narrowPeak file with `-q` `macs3_qvalue`.
- **SEACR** (`seacr_mode`: `auto` | `relaxed` | `stringent`, default `auto`): lab policy, SEACR has no broad mode. `relaxed` takes the
  threshold between the knee and the peak of the total-signal curve (more, wider peaks); `stringent` takes the peak. `auto` = `relaxed`
  for broad marks, `stringent` for narrow ones.

A group whose antibodies differ in broad / narrow status is a startup error (put them in separate groups, or force the mode with
`macs3_broad` / `seacr_mode`). The mode of every group is logged at startup and written to `metadata/peak_call_modes.tsv`.

## Retained / lost classification

For each `classifications.csv` row, reads are counted in the given regions for the reference and
treated groups and tested with an edgeR quasi-likelihood GLM, using the normalisation run's factors.
A region is **lost** when logFC ≤ `lfc_lost` (default −2) and FDR < `fdr` (0.05), **retained** when
logFC rounded to 1 decimal ≥ `lfc_retained` (default −1) or FDR ≥ `fdr`, and **unclassified**
otherwise (a region meeting both rules is lost). "Retained" means no evidence of loss, not proof of
stability. The full list of steps and defaults is in `../PIPELINE_STEPS.md` at the repository level.

## Peak QC tracks

Whenever a peak caller runs (`run_peak_qc: true`, the default; csaw, MACS3 and SEACR peak sets are all drawn, labelled `csaw <set>`, `MACS3 <group>`, `SEACR <group>`, in black, #0072B2 and #D55E00), `PEAK_QC_REGIONS` chooses up to
`peak_qc_max_regions` (15) regions to look at, and the tracks step draws them with every peak set and
class BED: per mark, the widest, the median-width and the densest-peak window of the reference peak
set; for each classification its minority class (≤ 50 % of the peaks) and a cross-mark discordant
region (lost in one mark, retained in another); a negative control (`peak_qc_negative`: a locus, or
`auto` = a peak-free housekeeping or long gene with low signal); the peak-free 250 kb tile with the
highest target signal (needs bigWigs; tiles whose IgG signal fold is above `peak_qc_nopeak_igg_max_fold`, default 1.5, are excluded when an IgG bigWig exists; `null` or `''` turns the cap off); and the user's `peak_qc_loci`. Windows are 100, 200 or 400 kb
depending on the anchor's width. The rules were reconstructed from the 13 QC regions of the EED
analysis (reproduced exactly or within 2.5 kb). These panels are for checking the calls by eye; they
are extremes chosen by rule, not typical regions, and say nothing about sensitivity or false-negative
rates. Outputs: `results/tracks/peak_qc/` (regions, details with folds and skip reasons) and
`results/tracks/plots/peak_qc.NN_<label>.svg`. These panels use an automatic scale bar, 7 pt gene labels (gene height 1.5, so dense clusters such as HOX stay legible) and a wider track-title
column so that long titles are not cut off.

## What the pipeline does not do

- Differential binding genome-wide (only classification of given regions).
- Motif analysis, gene-set enrichment beyond the optional KEGG panel of the annotation.
- Choose contrasts, thresholds or normalisation for you: they come from the CSVs and `config.yaml`,
  which an agent fills in only from what the user confirmed (`CONTEXT.md`).

## Reproducibility

- Every step's exact command is published next to its output (`.../log/*_command.sh`).
- `metadata/nf_logs/` holds the Nextflow logs; `metadata/pipeline_info/` the execution trace and
  report; Docker versions also write `metadata/host_paths.tsv` when Windows paths were translated.
- Re-running with the same inputs resumes from the work folder (in the Docker versions a
  per-project Docker volume, kept by default, or `metadata/work` with `work_volume: false`; on
  the cluster under `/project/rephilllab/nf_work/<user>/<project>`).

## Output layout

```
results/
  bam/final/            final BAMs (+ .bai) used by every later step; bam/spikein/, bam/stats/
  bigwig/               normalised bigWigs (<sample>.<method>.<units>.bw)
  scale_factors/        <norm_run>.scale_factors.tsv (primary), <norm_run>.scale_factors_all.tsv, methods/ (one table per method) + plot
  norm_qc/              MA plots (<run>.ma.svg, .ma.p2.svg, ...) and tables after normalisation
  peaks/{csaw,macs3,seacr}/
  classify/<id>/        <id>.{retained,lost,unclassified}.bed, per_region.tsv.gz, summary, MA plot
  heatmaps/{matrix,cluster,plots}/
  tracks/{ini,limits,plots,peak_qc}/
  annotation/  intersect/
  qc/                   fastqc, trimmomatic, bowtie2, alignment, markdup, fragment_sizes,
                        library_complexity, preseq, fingerprint, correlation, replicate_correlation, tss_enrichment,
                        frip, peaks, multiqc, report
intermediate_files/     trimmed FASTQs, aligned and filtered BAMs (on the cluster also the markdup BAMs),
                        SEACR bedGraphs and BIN_COUNTS tables; in the Docker versions only with
                        publish_intermediates: true (default false: the folder stays empty; when true, copies
                        of the work files, or hard links with work_volume: false)
```
