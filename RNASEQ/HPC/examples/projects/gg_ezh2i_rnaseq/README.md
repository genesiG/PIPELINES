# GG_EZH2i_RNAseq example (legacy project /path/to/GG_EZH2i_RNAseq, newest RNA-seq run 2026-08)

Files: `config.yaml` (flat, read by `run_pipeline.sh`), `samplesheet.csv` (12 libraries, DIPGVI cell line, sorted by
`sample_id` like the legacy count columns). Absolute FASTQ paths in the samplesheet; `fastq_dir` is not set.
Genome: T2T-CHM13v2.0 (HISAT2 index `chm13-2`, curated RefSeq GTF). Paired-end only.

The FASTQs in `/path/to/GG_EZH2i_RNAseq/Original_Data/fastq_files/` were regenerated from the hg19-aligned BAMs with
GG step_0d (BAM to FASTQ). They are not the original sequencing FASTQs. All 24 files (12 libraries x R1/R2) exist.

## How to run

```
mkdir my_gg_rnaseq && cd my_gg_rnaseq
cp /project/rephilllab/RNASEQ/examples/projects/gg_ezh2i_rnaseq/{config.yaml,samplesheet.csv,deg_comparisons.csv,deg_heatmaps.csv,deg_colors.csv,deg_label_rules.csv,deg_metadata.csv} .
/project/rephilllab/bin/run_rnaseq.sh        # uses config.yaml in the current directory
```

## Samplesheet

- `sample_id` = the legacy column name without `.R1` (for example `DIPGVI_CBX2.1A_1uM`).
- `cell_line` DIPGVI for all rows. `sgRNA`: sgCBX2.1 (CBX2.1A/B), sgCBX2.3 (CBX2.3A/B), sgLUC (LUCA/LUCB).
- `guide_replicate`: the letter A or B after the guide name. `treatment`: 1uM or DMSO (the last part of the name).
- `drug` (1uM or DMSO) and `target` (CBX2 or LUC) are the last two columns, as in the legacy `Metadata/sample_metadata.csv`.
- Metadata columns are carried into `results/samplesheet.tsv`. The differential expression is described below.

## Differential expression

`config.yaml` has `run_deg: true`: after the counts, the pipeline runs the DESeq2 analysis of the legacy
`Scripts/step_6_deg_analysis.Rmd`. Files:

- `deg_comparisons.csv`: model `treatment` (design `~treatment`, the 6 pairwise contrasts in the legacy loop order, three of them in the
  combined volcano figure).
- `deg_heatmaps.csv`: the treatment heatmap (3 row clusters).
- `deg_colors.csv`: LUC_DMSO #9ebcda, LUC_1uM #8c96c6, CBX2_DMSO #8c6bb1, CBX2_1uM #810f7c.
- `deg_label_rules.csv`: the sample label rewrites of the legacy heatmaps.
- `deg_metadata.csv`: the legacy sample metadata with `sample_id` = the legacy name without `.R1`. The legacy `treatment` column is
  `target_drug` (`CBX2_1uM`, `LUC_DMSO`, ...), which is not the samplesheet's `treatment` (1uM or DMSO), so the DEG steps read this file.

Reference levels are `treatment=LUC_DMSO,drug=DMSO,target=LUC`; volcano limits are those of the legacy plots.

`counts_only/` runs only the DEG steps from an existing counts table (`deg_counts`; no reads, no QC, no MultiQC), with its own
`deg_metadata.csv`. Run it from that directory after pointing `deg_counts` at a counts table.

Legacy outputs to compare with (not modified): `/path/to/GG_EZH2i_RNAseq/Analysis_Data/deseq/...` (per design folder:
`results_<comparison>.csv`, `analysis_summary.csv`, the volcano and heatmap SVGs, `DESeq.rds`, `rlogCounts.rds`). The results tables
(`log2FoldChange`, `padj`, effect) are expected to equal the legacy ones. The differences are the behaviour changes of
the DEG steps: the legacy analysis summary showed 0 up and 0 down genes (fixed), the legacy
QC plots coloured some groups with another group's colour and had a mismatched MDS legend (fixed).
The legacy interaction model (`~ target * drug`) is not part of this pipeline (see the main README, section on custom statistical tests).

**Colour note.** The legacy QC plots indexed the colour vector by factor position, so CBX2_1uM and LUC_1uM had each other's colour.
`deg_colors.csv` maps colours by group name. The volcano down / up colours are the new defaults `#46627A` / `#8c6bb1` (the legacy
volcano used `#810f7c` and gray; the legacy heatmaps used `#4d004b` - white - goldenrod); `deg_color_down` / `deg_color_up` change them.

## Differences from legacy (by design)

- Column names have no `.R1` suffix (legacy R1 stem used as the sample name).
- FastQC runs on both raw and trimmed reads (legacy: trimmed R1 only).
- featureCounts runs with `-T` (legacy had no thread option despite 8 cores); `gene_length.txt` comes from one sample
  (legacy wrote it from each job).
- The RSeQC strandedness check is reported in `results/qc/strandedness/strandedness.tsv` and warns on a mismatch with
  `strandedness: unstranded`; it does not change the counts.
- `remove_chroms` is empty here: chrM is kept. Legacy step_3 removed chrM. `remove_chroms: chrM` reproduces the legacy BAMs;
  it has no effect on T2T counts because the curated T2T GTF has no chrM features.
