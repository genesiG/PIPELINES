# CBX2 (EpiCypher) example (legacy GG_EPICYPHER_CBX2)

Files: `config.yaml`, `samplesheet.csv` (one row per library), `norm_runs.csv`, `norm_runs_spikein.csv` (spike-in alternative),
`peaksets.csv`, `classifications.csv`, `heatmaps.csv`, `tracks.csv`, `intersects.csv`.
Absolute FASTQ paths; `fastq_dir` null.

- 32 libraries, 2 conditions (DMSO, 1uM = `1uMcpd`).
- Normalisation: ChIPseqSpikeInFree (`normalisation: chipseqspikeinfree`, legacy `NORMALIZE_USING = CHIPSEQSPIKEINFREE`,
  CBX2 since 2026-10-14). Target bigWigs are CSIF RPKM at 50 bp with extendReads (round 11, 2026-10-15). IgG bigWigs are unscaled.
- Replicate letters: the RP_### numbers are libraries; A,B,... assigned by RP order within condition x antibody
  (e.g. IgG DMSO: RP_025=A, RP_033=B). Antibodies with one library per condition (H2AK119ub, EZH2, GST*) get A.
- control_id: IgG of the same condition and replicate letter (IgG exists with A and B in both conditions).
- AntiCBX2RP (antibody `CBX2`) and GSTCBX2 are very different experimental settings. They keep separate antibody, group and norm_run values and are never pooled or paired as replicates (peak calling, normalisation, replicate concordance or anything else).
- Broad peak settings (`broad_marks`): the default list plus H2AK119ub, EZH2, CBX2 and GSTCBX2 (lab decision 2026-10-27); H3K4me3, GSTControl and GSTCBX7 use narrow settings.
- group = `<condition>_<antibody>`. The legacy metadata GROUP column is the condition only (DMSO / 1uM).
- Legacy metadata/config mark no sample as excluded; exclude = false.

## Latest legacy parameters (2026-10-26) and toggles

`config.yaml` carries the latest legacy parameters; each key whose value differed between legacy runs lists its
alternative in a comment. Summary:

| Key | This project | Alternative (legacy or pipeline default) |
|---|---|---|
| `normalisation` | `chipseqspikeinfree` | `spikein` (spike-in era) with `norm_runs: norm_runs_spikein.csv` |
| `norm_runs` | `norm_runs.csv` (CSIF, legacy `config.py` NORMALIZATION_RUNS) | `norm_runs_spikein.csv` (spike-in table, same run names) |
| `bigwig_bin_size` | `50` | `100` (CPM bigWigs, briefly 2026-10-15) |
| `bigwig_units` | `rpkm` | `cpm` (spike-in CPM variant) |
| `bigwig_extend_reads` | `true` | - |
| `csaw_norm_bins` | `all` (legacy `normFactors(bins)` over all 10 kb bins) | `background` (pipeline default) |
| `csaw_max_frag` | `500` (csaw default) | `null` = `bowtie2_maxins` |
| `bowtie2_maxins` | `null` (bowtie2 default -X 500, as in the manuscript) | `700` (alternative, 2026-10 test runs) |
| `heatmap_y_label` | `'sf x rpkm'` | `null` = `normalised <bigwig_units>` |
| `tracks_max_value_default` | `shared` | `auto` (per-track autoscale) |
| `tracks_limit_stat` | `max_bin` (CBX2 legacy rule) | `mean_bins` (pipeline default, EED rule) |
| `annotation_extra_beds` | K27me3 classes (`K27me3_1uM.retained,.lost,.unclassified`) | empty (off) |

## Values I had to infer
- condition `1uM` from legacy `GROUP_PATTERNS` (file token `1uMcpd`).
- antibody: from legacy `ANTIBODY_PATTERNS` (RbIgG -> IgG, AntiCBX2RP -> CBX2, H3K27me3 -> K27me3, ...). H3K4me3 -> K4me3.
- norm_runs.csv: method `chipseqspikeinfree`; bin_size, cutoff, max_turns and min_first_turn (K27me2 only) copied from
  legacy `config.py` NORMALIZATION_RUNS (checked 2026-10-26). No IgG row: IgG bigWigs are unscaled. The legacy `IgG` entry
  (10000 / 1.2 / 0.999) is not used here.
- 2026-10-27 normalisation and peak-calling update: `norm_runs.csv` has the optional `ref_condition` column (reference condition of the MA plots: DMSO in every run); `method` stays the primary (ChIPseqSpikeInFree), TMM and CPM factors are computed as diagnostics (`norm_methods` default), `norm_spikein: false`. MACS3 and SEACR are on (`run_macs3` / `run_seacr: true`; every target has a control), with `macs3_broad` and `seacr_mode` at `auto` (broad / relaxed for the default `broad_marks`, which include K27me3 and K27me2).
- norm_run matched by regex search of the legacy `targets` on sample_id (CBX2 run target `AntiCBX2`); each sample matches exactly one run.
- peaksets: `submit_csaw.sh` runs `--small-filt 3 --large-filt 1 --merge-bp 500` for config SAMPLES (H3K27me3, H3K27me2, AntiCBX2RP, H2AK119ub)
  x ANTIBODY_TARGETS (DMSO, 1uMcpd) = 8 sets, plus K27me2 DMSO 2/1/50 and K27me3 DMSO 3/1.5/50. Total 10 rows.
  fdr 0.05 (`qvalue`) and lfc 1 are defaults from csaw_parameters.R, not in command line.
  max_width blank (no max width in legacy); blank cells take the `csaw_max_width` fallback (0 = no limit).
- control_group = IgG group of the same condition (legacy globs the condition's RbIgG BAMs).
- peakset names are new.
- Samples without a peakset (K4me3, GST*, EZH2) are only normalised/bigWig'd, as in the legacy run.

## Classifications (`classifications.csv`, `run_classify: true`)

Classifies the regions of two peak sets, 1uM vs DMSO (ref = DMSO, treat = 1uM):
`K27me3_1uM` on `K27me3_DMSO_s3_l1.5_m50` (DMSO_K27me3 vs 1uM_K27me3) and `K27me2_1uM` on `K27me2_DMSO_s2_l1_m50`
(DMSO_K27me2 vs 1uM_K27me2). norm_run and thresholds are blank (shared norm_run, `classify_*` defaults:
lfc_retained -1, lfc_lost -2, fdr 0.05, the same thresholds as round 11).
Outputs go to `results/classify/<id>/`; `<id>.retained`, `<id>.lost` and `<id>.unclassified` are the class names used in
heatmaps, tracks, intersects and annotation. The legacy classification BEDs under `Analysis_Data/` are not used by these CSVs.

## Round-11 figure examples (heatmaps.csv, tracks.csv, intersects.csv)

Legacy sources: `GG_EPICYPHER_CBX2/Logs/round11/stage_rpkm.sh` (heatmaps), `stage_csif.sh` (classification, intersections),
`redraw_cdkn2a_onerep.sh`, `stage_tracks_fix.sh` and `stage_discordant.sh` (tracks), `Scripts/utils/cbx2.ini`, `cdkn2a.ini`.

**heatmaps.csv (1 row: `cbx2_regions_10kb`)**: six region groups in one deepTools matrix (10 kb window, bin 100, no clustering,
`z_max` and `y_max` = `auto`, one limit per antibody). Regions and labels, in order:
`K27me2_1uM.retained` (Retained H3K27me2), `K27me2_1uM.lost` (Depleted H3K27me2), `K27me3_1uM.retained` (Retained H3K27me3),
`K27me3_1uM.lost` (Depleted H3K27me3), `lost_h3k27me3_in_retained_h3k27me2` (Depleted H3K27me3 intersect Retained H3K27me2),
`lost_h3k27me2_in_lost_h3k27me3` (Depleted H3K27me2 intersect Depleted H3K27me3).
The 8 samples are the round-11 `step_8b` samples: `H3K27me2`, `H3K27me3`, `CBX2`, `EZH2` (DMSO / EZH2i, R1 libraries).
The legacy stage plotted each region set as its own figure; the labels are the legacy `--regions_label` values.

**tracks.csv (6 rows)**:
- `cbx2_cdkn2a`: CDKN2A, flank 125000, six one-replicate tracks (R1 of H3K27me3, H3K27me2 and CBX2, DMSO and 1uM;
  `redraw_cdkn2a_onerep.sh`), bed retained H3K27me2.
- `cbx2_cbx2`: CBX2, flank 25000, 12 tracks (R1 and R2 of each mark), beds retained and lost H3K27me2 (`cbx2.ini`).
- `cbx2_discordant_FOXD1`, `cbx2_discordant_ZIC1_ZIC4`, `cbx2_discordant_HOXD13_EVX2`, `cbx2_discordant_FOXG1`: the region
  windows of `stage_discordant.sh` (flank 0); six tracks as in CDKN2A, bed retained H3K27me2.
- Y limits: blank `max_value` = shared per antibody (see the reference table above). Scale bar 50000.
- Bed colours and titles are in `bed_colors` and `bed_labels` (legacy ini colour #8c96c6).

**intersects.csv (2 rows, class names as input, `-f 0.05 -r`, `output: overlap`)**: `lost_h3k27me3_in_retained_h3k27me2`
(bed_a = `K27me3_1uM.lost`, bed_b = `K27me2_1uM.retained`) and `lost_h3k27me2_in_lost_h3k27me3`
(bed_a = `K27me2_1uM.lost`, bed_b = `K27me3_1uM.lost`). Legacy `step_7b_intersect_beds.py --extra_flags '-f 0.05 -r'`; `output: overlap` reproduces its overlap pieces (the paper heatmaps).
These intersect names are the heatmap regions above.

**annotation**: `run_annotation: true` annotates the csaw peak sets and, through `annotation_extra_beds`, the K27me3 classes
(legacy round 9 annotated K27me3 Lost, Retained and Unclassified).

**Inputs**: the pipeline's own bigWigs (sample ids in the CSVs), the GTF from `config.yaml`, and the classifications and
intersections produced by this run (no legacy files are needed).

## Peak QC tracks (run_peak_qc, on by default)

With a peak caller on (here csaw) the pipeline also picks about 8-15 informative regions automatically (largest, median and densest peak per mark, minority retained/lost classes, cross-mark discordant regions, a negative control, a peak-free high-signal tile) plus the loci of `peak_qc_loci` in `config.yaml` (HOXA and CDKN2A, K27-mark loci; negative control GAPDH) and draws them as the track set `peak_qc` with the bigWigs of the peak sets' groups (and their IgG controls) and all peak sets and retained/lost BEDs as BED tracks.
Outputs: `results/tracks/peak_qc/` (`regions.tsv`, `regions.detail.tsv`, `params.json`) and `results/tracks/plots/peak_qc.NN_<label>.svg`. Set `run_peak_qc: false` to switch it off.
