# nucleation example

Files: `config.yaml` (flat pipeline params, references from legacy `Scripts/config.py`), `samplesheet.csv`
(one row per library), `norm_runs.csv` (= `NORMALIZATION_RUNS`), `peaksets.csv` (= round-9 csaw sets),
`heatmaps.csv`, `tracks.csv`.
All FASTQ/BAM paths are absolute; `fastq_dir` is null.

- 16 CUT&RUN libraries (4 conditions x IgG/K27me3 x replicates A,B), norm_runs.csv method `cpm` (plain RPKM bigWigs, see the lab decision section).
- control_id: each K27me3 library points to the IgG of the same condition and replicate letter.
- group = `<condition>_<antibody>` (the legacy metadata GROUP column is the condition only, e.g. K27MKO_EZH2i).
- Legacy metadata marks no sample as excluded; exclude = false everywhere.

## Latest legacy parameters (2026-10-26)

`config.yaml` reproduces the latest NUCLEATION settings. Keys whose value differed between legacy runs carry
their alternatives in a trailing comment. The table lists each toggle; the value in this example is the first one.

| Key | Value here | What it switches / alternatives |
|---|---|---|
| `normalisation` | `chipseqspikeinfree` | Informational: the primary is `cpm` in `norm_runs.csv` (`norm_run` `cpm` gives `bigwig_multiplier` 1, plain RPKM bigWigs; lab decision 2026-10-28). Setting `cpm` here would skip `norm_runs.csv` and the ChIPseqSpikeInFree/TMM diagnostics. Alternatives: `spikein`, `tmm` |
| `bigwig_bin_size` | `50` | bamCoverage `--binSize`. Alternatives: `10` (before 2026-10-15), `100` (CPM, briefly on 2026-10-15) |
| `bigwig_units` | `rpkm` | bigWig units. Alternatives: `cpm` (before 2026-10-15), `none` |
| `bigwig_extend_reads` | `true` | bamCoverage `--extendReads` (fragments) |
| `csaw_norm_bins` | `all` | TMM normalisation over all 10 kb bins (legacy `callpeaks_csaw.R`). Alternative: `background` (pipeline default; TMM on background bins) |
| `csaw_max_frag` | `500` | csaw `max.frag` (legacy csaw default). Alternative: `null` = `bowtie2_maxins` |
| `bowtie2_maxins` | `null` | bowtie2 `-X` optional toggle; `null` = omitted (bowtie2 default 500, as in the manuscript). Alternative: `700` (2026-10 test runs) |
| `heatmap_y_label` | `'rpkm'` | Profile y-label (unscaled here; `'sf x rpkm'` for CBX2 and EED). Alternative: `null` = `normalised <bigwig_units>` |
| `heatmap_legend_location` | `best` | Legend of the deepTools summary plot. Alternative: `none` (pipeline default) |
| `tracks_max_value_default` | `shared` | Blank `tracks.csv` max_value: one y-range per antibody at a locus. Alternatives: `auto` (per-track autoscale), a number |
| `tracks_limit_stat` | `mean_bins` | Shared-limit statistic (current behaviour). Alternative: `max_bin` (CBX2 2026-10-15 rule: max over the locus window of each bigWig's own bins) |
| `run_intersect` | `false` | BED intersections (`intersects` key, not set). Off: the legacy intersections used the dropped K27MKO classification BEDs |
| `run_annotation` | `true` | Peak annotation of the csaw peak sets |

Shared by all three projects: `bigwig_bin_size 50`, `bigwig_units rpkm`, `bigwig_extend_reads true`,
`tracks_max_value_default shared`.

## Lab decision (2026-10-28): plain RPKM bigWigs

- Lab decision for this paper: the NUCLEATION bigWigs are plain RPKM. `norm_runs.csv` `method` is `cpm` for both runs,
  which gives `bigwig_multiplier` 1, so `<sample>.cpm.rpkm.bw` is plain RPKM.
- ChIPseqSpikeInFree is computed as a diagnostic only. `cutoff` and `max_turns` stay in `norm_runs.csv` for that purpose and are not applied to the bigWigs.
- CBX2 and EED stay ChIPseqSpikeInFree x RPKM (`bigwig_units rpkm`).
- The heatmap y-label is `'rpkm'` here (unscaled); `'sf x rpkm'` for CBX2 and EED.

## Values I had to infer

- condition: legacy `GROUP_PATTERNS` labels (`500nM` -> `EZH2i`: K27M_EZH2i, K27MKO_EZH2i) rather than the file-name token `500nM`.
- sample_id strips `_R1`/`_R2` (K27M_500nM_IgG_A).
- norm_run: each sample matched by regex search of the `targets` pattern on sample_id; every library matches exactly one run.
- norm_runs.csv values verbatim from config.py (K27M_K27me3 0.86, K27MKO_K27me3 0.989, IgG runs 0.99; bin 10000, cutoff 1.2): no difference found.
- 2026-10-27 normalisation and peak-calling update: `norm_runs.csv` has the optional `ref_condition` column (reference condition of the MA plots: K27M_K27me3 -> K27M_DMSO, K27MKO_K27me3 -> K27MKO_DMSO (each run holds one genotype, so the reference condition is that genotype's DMSO)); norm_runs.csv `method` is `cpm` (lab decision 2026-10-28, see the section below); ChIPseqSpikeInFree, TMM and CPM factors are computed as diagnostics (`norm_methods` default), `norm_spikein: false`. MACS3 and SEACR are on (`run_macs3` / `run_seacr: true`; every target has a control), with `macs3_broad` and `seacr_mode` at `auto` (broad / relaxed for the default `broad_marks`, which include K27me3 and K27me2).
- peaksets: small/large/merge are from the round-9 `stage_A.sh` command lines (K27MKO default = csaw_parameters.R 2/0.5/500 for config.SAMPLES K27MKO_DMSO and K27MKO_500nM; K27M_DMSO 3/2/100).
  fdr 0.1 (`qvalue`) and lfc 1 are defaults, not in command line (csaw_parameters.R).
  max_width is blank: the legacy script has no max width.  Needs lead decision on how the pipeline reads blank.
- control_group = IgG group of the same condition (legacy script globs all IgG BAMs of the condition).
- csaw peakset names are new (legacy used file names like K27MKO_DMSO_K27me3.w150...bed).

## Round-9 figure examples (heatmaps.csv, tracks.csv)

Legacy source: the round-9 `stage_B.sh` of the nucleation analysis, with the 2026-10-26 settings above.
Keys: `config.yaml` sets `heatmaps`, `tracks` and the `run_*` switches.

**heatmaps.csv (6 rows)**: `K27MKO_K27me3_DMSO_peaks_{10,5,2}kb_k{2,3}`, region `K27MKO_DMSO_K27me3` (csaw peakset),
samples `K27MKO_DMSO_K27me3_A`, `K27MKO_500nM_K27me3_A`, bin 100, colour `#8c6bb1`, z_min `0,0`.
`z_max` and `y_max` are `auto` (per-sample limit from HEATMAP_LIMITS, one `auto` expands to all samples).
- Source: `step_8a_computematrix.py` (3 windows) then `step_8c_plotheatmap_clustering.py --n_clust 2` (k2 rows, all samples clustered)
  and `--n_clust 3 --cluster_using "1"` (k3 rows, `cluster_samples` = first sample).
- Not reproduced: `step_7e_computematrix.py` + `step_7f_plotheatmap.py` (resistant and vulnerable K27MKO classification BEDs in one matrix,
  4 samples). Dropped with the classification BEDs; the heatmap schema takes one region set per heatmap.
- `sort_using` left empty: legacy `--sortUsingSamples "1 2"` cannot be written as one sample_id.

**tracks.csv (1 row)**: `nucleation_tracks`, loci `CBX2 125000`, `CDKN2A 100000 / 125000 / 200000 / 250000 / 400000`
(six figures, the legacy `step_9a_genome_tracks.py --genes ... --flanks ...` call).
- Samples `K27MKO_DMSO_K27me3_A`, `K27MKO_500nM_K27me3_A`, `K27M_DMSO_K27me3_A`; scalebar 50000.
- Limits: `group:KO,group:KO,group:K27M`. The two K27MKO tracks share one y-range; the K27M DMSO track keeps its own
  scale (legacy ini max 190 vs 19 for the K27MKO tracks).
- Beds: the csaw peak sets `K27MKO_DMSO_K27me3` (`#8c6bb1`, "K27MKO DMSO peaks") and `K27M_DMSO_K27me3`
  (`#8c96c6`, "K27M DMSO peaks"), colours and titles from the legacy `Scripts/utils/.ini`.
- Dropped: the K27MKO retained/lost classification BEDs (legacy step_7d K27MKO classification, declared invalid
  on 2026-10-26).
- The pipeline's own bigWigs (sample ids above) and the GTF from `config.yaml`.

**annotation**: `run_annotation: true` on the csaw peak sets; no CSV is needed.
