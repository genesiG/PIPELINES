# AV_ATACSEQ_2026 example (legacy project /path/to/AV_ATACSEQ_2026)

Files: `config.yaml` (flat, read by `run_pipeline.sh`), `samplesheet.csv` (39 libraries, 30 included), `peak_merges.csv`.
Absolute FASTQ paths in the samplesheet; `fastq_dir` is null. Genome: T2T-CHM13v2.0, paired-end only.

Defaults are the newer legacy settings (AV_ATACSEQ_2026). Where the new pipeline fixes a legacy bug, the fix is the default
and a toggle in `config.yaml` reproduces the legacy behaviour (table below).

## How to run

```
mkdir my_av_atacseq && cd my_av_atacseq
cp /project/rephilllab/ATACSEQ/examples/projects/av_atacseq_2026/{config.yaml,samplesheet.csv,peak_merges.csv} .
/project/rephilllab/bin/run_atacseq.sh        # uses config.yaml in the current directory
```

## Toggles (what each reproduces)

| Key (config.yaml) | Value here | Toggle that reproduces legacy | Legacy behaviour (source) |
|---|---|---|---|
| `mapq` | 30 | `mapq: 0` | `samtools view -q mapq` with the literal word, read as 0 (`Scripts/step_3_process_bam.py:140`) |
| `remove_chroms` | chrM | `remove_chroms: ''` | chrM kept (no removal step in legacy) |
| `tmm_group_column` | tmm_group | `tmm_group_legacy` | TMM groups by file-name substring: K27M includes K27MKO, x42c includes HGGx42c (`Scripts/config.py:127`, `scaleFactors.R:119-124`; membership in `normalization/sample_order_*.txt`) |
| `tn5_shift` | false | `tn5_shift: true` | not in legacy; `true` applies deepTools `alignmentSieve --ATACshift` |
| `macs3_fe_cutoff` | 2.0 | `4.0` | legacy `config.py:150` says 2.0, but `merge_peaks.sh` merged the fe4.0 files (see Peak merges) |
| `tmm_max_frag` | auto | `1999` | legacy used one global max.frag of 1999 (`normalization/max_frag_size.txt`), which came from the first sample only (`getPEsizes.R:78`) |

Other changes vs legacy are listed under "Expected differences" below. They have no toggle.

## Samplesheet

- 39 libraries, all rows of `Original_Data_and_Metadata/Metadata/sample_metadata.csv` (BOM and CRLF stripped). The metadata has no ATAC28 row.
- `exclude: true` for 9 libraries: 11284_whole_A, 11656_whole_A, 11656_whole_B, HGGx6c_low/high/whole_A, x6c_low/high/whole_A.
  They are not in `Metadata/paired_samples.txt` and have 9-12 % alignment. They keep their legacy TMM group names
  (11284, 11656, x6c) as in the spec, so a group whose only members are excluded (11284, 11656, x6c) is a validation point.
- 30 included libraries = the 29 lines of `paired_samples.txt` + BT54_high_A.
- **BT54_high_A is INCLUDED.** Legacy dropped it by accident: it is missing from `paired_samples.txt` (only BT54_high_B is listed),
  but it is in the legacy TMM group (`normalization/sample_order_BT54.txt`), and its alignment rate is 89.31 %.
  The legacy fe2.0 BT54_high narrowPeak file was called with BT54_high_B only (`peaks/log/BT54_high.peak_calling.batch`).
  The new BT54_high peak group pools A and B, so its counts are not comparable with the legacy peak counts.
- **TS603_low_A**: the metadata `sample_name` is `TS503_low_A` (typo; `condition` is TS603_low). The FASTQs are `TS603_low_A_R1/R2`
  and the legacy TMM table uses `TS603_low_A_R1`. The samplesheet uses `TS603_low_A`. The legacy `paired_samples.txt` has the same
  typo in the R2 control column (`TS503_low_A_R2`).
- `tmm_group` (exact): BT54, BT88, TS603, TS667, 6349, x42c, HGGx42c, K27M, K27MKO. `tmm_group_legacy`: the membership of the legacy
  `sample_order_<G>.txt` files (K27M contains K27MKO, x42c contains HGGx42c, x6c contains HGGx6c). Each of the 39 libraries is in
  exactly one legacy file.
- `peak_group`: the legacy `REPLICATES` patterns (`Scripts/config.py:135-143`, applied with `re.search` as in `step_4c_macs2.py:149`)
  resolved to explicit groups. Every included library matches exactly one pattern; no excluded library matches any pattern.

### Peak groups (membership, from samplesheet.csv)

| peak_group | sample_id |
|---|---|
| BT54_low | BT54_low_A |
| BT54_high | BT54_high_A, BT54_high_B |
| BT54_whole | BT54_whole_A, BT54_whole_B |
| BT88_whole | BT88_whole_A, BT88_whole_B |
| BT88_high | BT88_high_A |
| TS603_low | TS603_low_A |
| TS603_high | TS603_high_A |
| TS603_whole | TS603_whole_A, TS603_whole_B |
| TS667 | TS667_whole_A, TS667_whole_B |
| 6349_low | 6349_low_A, 6349_low_B |
| 6349_high | 6349_high_A, 6349_high_B |
| HGGx42c_low | HGGx42c_low_A |
| HGGx42c_high | HGGx42c_high_A |
| HGGx42c_whole | HGGx42c_whole_A |
| x42c_low | x42c_low_A |
| x42c_high | x42c_high_A |
| x42c_whole | x42c_whole_A |
| K27M_low | K27M_low_A |
| K27M_high | K27M_high_A |
| K27M_whole | K27M_whole_A |
| K27MKO_low | K27MKO_low_A |
| K27MKO_high | K27MKO_high_A |
| K27MKO_whole | K27MKO_whole_A |

(Legacy `REPLICATES` pattern `^x42c_*` is anchored, so it does not match `HGGx42c_*`. `K27M_*` does not match `K27MKO_*`. `TS667` matches the TS667_whole A and B libraries.)

## Peak merges (peak_merges.csv)

The CSV has no comments, so here: legacy `Scripts/merge_peaks.sh` ran `bedops --merge` over `*high.qval0.05*fe4.0*.narrowPeak`
(8 files: BT54, BT88, 6349, HGGx42c, x42c, K27M, K27MKO, TS603 high) and `*low...` (7 files: BT54, 6349, HGGx42c, x42c, K27M,
K27MKO, TS603 low) into `high_gd2_peaks_fe4.0.bed` and `low_gd2_peaks_fe4.0.bed`. Those merged files are not in `peaks/` now.
The new rows `high_gd2_peaks` and `low_gd2_peaks` use the same groups, but the peaks of this config's `macs3_fe_cutoff` (2.0 by default),
so the merged regions differ from the legacy fe4.0 merge. Set `macs3_fe_cutoff: 4.0` to merge the legacy fe4.0 set.

## Expected differences vs legacy

Results will differ from the legacy run for these reasons (each one, except the first, is a deliberate change with its toggle above):

1. MAPQ 30 instead of the literal `-q mapq` (= 0): fewer low-MAPQ fragments in every BAM, bigWigs, TMM and peak counts.
2. chrM removed (`remove_chroms: chrM`): chrM reads drop out of the BAMs and of TMM/bigWig/peaks. Legacy kept them.
3. Exact TMM groups: K27M and K27MKO, and x42c and HGGx42c, are normalised separately. The legacy scale factors do not apply.
4. Blacklist: pair-aware removal (`bedtools pairtobed -type neither`, both mates removed if either overlaps) instead of the
   BAM -> BEDPE -> BAM round trip. The round trip set MAPQ 255 and dropped sequences.
5. max.frag: the true maximum over each group's BAMs instead of the first sample's (1999).
6. BT54_high_A included (legacy dropped it). The BT54_high peak group therefore pools two libraries (legacy: one).
7. TS603_low_A sample name fixed (metadata typo TS503_low_A).
8. TSS enrichment path fixed (legacy `genomes/t2t.curated...` does not exist).
9. Trimmomatic: biocontainer 0.40 instead of the legacy 0.32 jar. Same adapter file and steps.
10. Tn5 shift is off (legacy did not shift either, so this is no change).

## Downstream analyses

`config.yaml` turns on the three downstream switches (`run_differential`, `run_annotation`, `run_motifs`) with the AV_ATACSEQ_2026
values, and names two CSVs: `comparisons.csv` and `motif_targets.csv`. The rules of those files are in the top README, section
"Downstream analyses"; the older AV and AM values are `# TOGGLE` lines in `config.yaml`.

- `comparisons.csv`: the five per-line comparisons, `BT54_high_vs_low`, `6349_high_vs_low`, `TS603_high_vs_low`,
  `HGGx42c_high_vs_low` and `K27MKO_high_vs_low` (case = the `_high` condition, control = the `_low` condition, each with its own
  peak group), and `gd2_high_vs_low`: every `_high` condition against every `_low` condition, block = `tmm_group` (the cell line,
  so `HGGx42c` holds the x42c libraries and `K27MKO` and `K27M` are separate lines), `pooled` true, `fixed_dispersion` 0.15 (legacy
  gd2 chunk). The gd2 peak lists are the legacy peak files: 6349, BT54, HGGx42c, K27MKO, K27M, TS603 and x42c, high and low. BT88 has
  no legacy peak files, so its high condition is in the gd2 case samples but not in its peak lists.
- `motif_targets.csv`: the ST8SIA1 gene target (peaks of the gd2 `all_sites` table whose geneId contains ST8SIA1) and the ST8SIA1
  cCRE regions. The cCRE BED is `st8sia1_ccres_hs1.bed`, a placeholder: the legacy cCREs were lifted over to hs1 by a project
  pre-step that the pipeline does not run. Put the lifted BED in this directory (or change the path) before a run with run_motifs on.
- Region sets: `diff_regions: peaks,promoters` (the promoter table, for the motif analysis) and `motif_regions: promoters`.

Before a real run, decide the gd2 sample set: the rule "every `_high` against every `_low`" selects 18 active libraries of this
samplesheet (10 high, 8 low), while the legacy gd2 analysis used 12 BAMs of the processed folder.
