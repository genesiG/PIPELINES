# EED example (legacy GG_EED_CUTANDRUN)

Files: `config.yaml`, `samplesheet.csv` (one row per library), `norm_runs.csv`, `peaksets.csv`, `classifications.csv`, `heatmaps.csv`, `tracks.csv`. Absolute FASTQ paths; `fastq_dir` null.

- 35 libraries: DMSO, EEDi_100nM/500nM, EPZ_100nM/500nM, MAK_100nM/500nM x IgG/K27me2/K27me3. Normalisation `chipseqspikeinfree`.
- Each condition has ONE IgG (replicate A); every K27me2/K27me3 library (A and B) uses it as control_id.
  DMSO has K27me2/K27me3 A,B and IgG A.
- group = `<condition>_<antibody>` (dose-specific, e.g. EEDi_100nM_K27me3). The legacy metadata GROUP column collapses doses
  (DMSO/EEDi/EPZ/MAK, 4 groups) and has no antibody, so it was not used.
- Legacy metadata/config mark no sample as excluded; exclude = false.
- norm_runs.csv from config.py: K27me2 50000/0.7/0.876 (round 13, 2026-10-19: tuned to the user's MA-cloud
  targets; was 10000/1.2/0.925), K27me3 15000/1.2/0.9999. IgG is unscaled (empty norm_run).
- 2026-10-27 normalisation and peak-calling update: `norm_runs.csv` has the optional `ref_condition` column (reference condition of the MA plots: DMSO in every run); `method` stays the primary (ChIPseqSpikeInFree), TMM and CPM factors are computed as diagnostics (`norm_methods` default), `norm_spikein: false`. MACS3 and SEACR are on (`run_macs3` / `run_seacr: true`; every target has a control), with `macs3_broad` and `seacr_mode` at `auto` (broad / relaxed for the default `broad_marks`, which include K27me3 and K27me2).

## Latest legacy parameters (2026-10-26)

`config.yaml` opens with the block "Latest legacy parameters (2026-10-26)". The value on each line is the latest legacy
setting; the comment after a key lists its alternatives (toggles). Summary:

- bigWigs 50 bp, RPKM, extendReads (2026-10-15).
- csaw normalises over all 10 kb bins (`csaw_norm_bins: all`, legacy normFactors) with max.frag 500 (`csaw_max_frag: 500`).
- bowtie2 with no `-I`/`-X` (`bowtie2_maxins: null`; bowtie2 default -X 500, as in the manuscript; `700` was a 2026-10 test-run alternative).
- Heatmap y label `sf x rpkm`; shared track scales per antibody (`tracks_max_value_default: shared`, `tracks_limit_stat: mean_bins`).

| Key (config.yaml) | Value here | Alternative(s) |
|---|---|---|
| `bigwig_bin_size` | 50 | 10 (before 2026-10-15); 100 (briefly, with cpm) |
| `bigwig_units` | rpkm | cpm (briefly, with 100 bp bins) |
| `csaw_norm_bins` | all | background (pipeline default) |
| `csaw_max_frag` | 500 | null (= `bowtie2_maxins`) |
| `bowtie2_maxins` | null (bowtie2 default -X 500) | 700 (test-run alternative) |
| `csaw_defaults` K27me2 `large_filt` | 1 | 1.5 (acceptable); 2 (rejected) |
| `norm_runs` K27me2 row | round 13, 50000/0.7/0.876 | round 12, 10000/0.7/0.844 (edit the row in norm_runs.csv) |
| `heatmap_y_label` | `sf x rpkm` | null (normalised rpkm, pipeline default) |
| `tracks_max_value_default` | shared | auto (per-track autoscale); a number (legacy fixed maxima 100 K27me3, 32 K27me2) |
| `tracks_limit_stat` | mean_bins | max_bin |

## Values I had to infer
- condition = file-name prefix (e.g. EEDi_100nM) verbatim.
- norm_run by regex search of `targets` (IgG, K27me2, K27me3) on sample_id; each sample matches exactly one run.
- bigwig_bin_size 50 = config.py BIN_SIZE (since 2026-10-15; 10 before, listed as a toggle).

## csaw peaks and classification (new)
- `peaksets.csv`: `DMSO_K27me3` (DMSO_K27me3 vs DMSO_IgG) and `DMSO_K27me2` (DMSO_K27me2 vs DMSO_IgG). All threshold cells are blank,
  so they resolve from `csaw_defaults` in `config.yaml` (K27me3 3 / 1.5 / 50, K27me2 2 / 1 / 50), then the `csaw_*` globals.
  K27me2 `large_filt` 1 is the chosen setting (1.5 acceptable, 2 rejected; legacy `GG_EED_CUTANDRUN/Logs/round14/k27me2_filter_sweep.tsv`).
  The legacy round-9 EED run made no csaw peak set, so these two rows are new. Legacy round 14 (2026-10-22) called them with
  the legacy scripts; the pipeline reproduces those peak sets exactly (K27me3 25,814, K27me2 66,685; checked 2026-10-23)
  with `csaw_norm_bins: all` and `csaw_max_frag: 500` (legacy settings, set in this example's `config.yaml`; the pipeline
  default `background` is the toggle). The example's peak sets are therefore the round-14 ones.
- `run_csaw: true` and `run_classify: true`.
- `classifications.csv`: `EZH2i500_K27me3` and `EZH2i500_K27me2` (ref DMSO of the same mark, treat EPZ_500nM of the same mark).
  norm_run and thresholds are blank (shared norm_run; `classify_*` defaults). Outputs are in `results/classify/<id>/`;
  the heatmap and track rows below use its `<id>.retained` and `<id>.lost` regions.

## Figure examples (heatmaps.csv, tracks.csv; no intersects.csv)

Legacy sources: `GG_EED_CUTANDRUN/Logs/round9/stage_B1.sh` (heatmaps), `stage_A.sh` (CDKN2A track) and
`GG_EED_CUTANDRUN/Logs/round14/stage_own_regions.sh` with `Scripts/utils/genometracks_own_regions_qc.ini` (QC track set).
`config.yaml` sets `heatmaps` and `tracks`; `run_annotation: false` and `run_intersect: false`.
No `intersects.csv` is written because EED has no intersect figure in round 9.

**heatmaps.csv (4 rows, 10 kb, bin 100)**: `ezh2i500_{retained,lost}_{k27me2,k27me3}_10kb`.
- Regions: one region set per heatmap (`EZH2i500_K27me2.retained`, `EZH2i500_K27me2.lost`, resp. K27me3), labels
  `Retained H3K27me2`, `Depleted H3K27me2` (resp. H3K27me3), as the legacy round-14 `own_{retained,lost}_h3k27me*_10kb` heatmaps.
  One set per heatmap keeps the `auto` scale of each set: the 481 retained K27me3 regions reach ~110 under EZH2i, the
  24k lost ones ~40 in DMSO, so a shared two-group heatmap compresses one of them. Replaces the round-9 `retained_*` BED rows.
- Samples: `DMSO`, `EPZ_500nM`, `EEDi_500nM`, `MAK_500nM` (K27me2 then K27me3, replicate A). Labels as in the round-9 call, e.g. `H3K27me2\nEEDi 1`.
- Colours as in round 9 (#8c96c6 x4 then #8c6bb1 x4); z_min 0.
- z_max and y_max are `auto`: the pipeline computes one limit per antibody from the plotted matrix (legacy fixed values were 15 / 30).
- Y-axis label `sf x rpkm` (global `heatmap_y_label` in `config.yaml`).

**tracks.csv (2 rows)**:
- `eed_cdkn2a`: locus `CDKN2A`, flank 125000, 14 bigWigs (A replicates): K27me3 `DMSO`, `EPZ_100nM`, `EPZ_500nM`,
  `EEDi_100nM`, `EEDi_500nM`, `MAK_100nM`, `MAK_500nM` and the same seven for K27me2. Scalebar 50000.
  - Labels are the legacy ini titles (for example `H3K27me3 EZH2i 100nM`).
  - max_value is blank: shared, so each antibody has one y-range at this locus (legacy fixed maxima were 100 for K27me3 and 32 for K27me2;
    set numbers or `auto` in the cell to get the old behaviour).
  - beds: `EZH2i500_K27me3.retained` (blue, `#2a78d6`, "Retained") and `EZH2i500_K27me3.lost` (red, `#e34948`, "Depleted"), via `bed_colors` and `bed_labels`.
- `eed_own_regions_qc`: the 13 regions of `GG_EED_CUTANDRUN/Analysis_Data/genome_tracks/own_regions_qc/regions.tsv` as
  `chr:start-end` (no flank; the region names of the file are not carried into the track names). Nine bigWigs: K27me3 and K27me2
  DMSO A/B and EPZ_500nM A/B, and `DMSO_IgG_A`. Labels and colours from the ini: DMSO #8c6bb1 (K27me3) / #8c96c6 (K27me2),
  EZH2i #b9a6d3 / #c3c8e0, IgG #808080. Six beds: the csaw DMSO peak sets `DMSO_K27me3`, `DMSO_K27me2` (black) and the EZH2i 500 nM
  `retained` (blue `#1f78b4`) and `lost` (red `#e31a1c`) classes of each mark, labelled as in the ini ("H3K27me3 retained", ...).
  Blank `max_value`: shared per antibody, the IgG its own group. Scalebar 50000 (the ini's scalebar has no length; chosen here).

**The user must supply**:
- The pipeline's own bigWigs (sample ids above) and the GTF from `config.yaml`.

## Peak QC tracks (run_peak_qc, on by default)

With a peak caller on (here csaw) the pipeline also picks about 8-15 informative regions automatically (largest, median and densest peak per mark, minority retained/lost classes, cross-mark discordant regions, a negative control, a peak-free high-signal tile) plus the loci of `peak_qc_loci` in `config.yaml` (HOXA and CDKN2A, K27-mark loci; negative control GAPDH) and draws them as the track set `peak_qc` with the bigWigs of the peak sets' groups (and their IgG controls) and all peak sets and retained/lost BEDs as BED tracks.
Outputs: `results/tracks/peak_qc/` (`regions.tsv`, `regions.detail.tsv`, `params.json`) and `results/tracks/plots/peak_qc.NN_<label>.svg`. Set `run_peak_qc: false` to switch it off.
The hand-written `eed_own_regions_qc` row of `tracks.csv` is kept as the legacy reproduction next to the automatic set.
