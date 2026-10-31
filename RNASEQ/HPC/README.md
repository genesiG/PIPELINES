# RNASEQ: paired-end RNA-seq pipeline (Rephill Lab)

From FASTQ pairs or BAMs to one counts table per project: FastQC, Trimmomatic, HISAT2, samtools filtering, Picard duplicate
marking, an automatic strandedness check (RSeQC), featureCounts, a merged counts table and MultiQC. It runs on UPenn PMACS (LSF),
with every step in its own container. Optionally (`run_deg: true`) it then runs a DESeq2 differential expression analysis; by default the counts table is the end point.

It is a port of the lab's GG_EZH2i_RNAseq analysis (the newest RNA-seq project). The behaviour changes from the legacy scripts
are listed in the GG example's README (`examples/projects/gg_ezh2i_rnaseq/README.md`, "Differences from legacy").

## Quick start

```bash
mkdir my_project && cd my_project
run_rnaseq.sh init                 # scaffolds the project and git-inits it; never runs anything
# edit samplesheet.csv: one row per library
# edit config.yaml: project_name, hisat2_index, gtf, strandedness
run_rnaseq.sh                      # submits one LSF head job and returns
```

`run_rnaseq.sh` is a symlink in `/project/rephilllab/bin/` that points at this repository's `run_pipeline.sh`.
`init` copies `config.yaml`, `samplesheet.csv`, `examples/`, `AGENTS.md` and `CLAUDE.md` into the project (it never overwrites an
existing file) and creates `original_data/fastq`, `metadata`, `intermediate_files`, `results` and `utils`.

The head job is queued on `rhel9` (1 CPU, 8 GB, 168 h). It runs the Nextflow head process, which submits every step as its own
LSF job. A run is re-started by running `run_rnaseq.sh` again: `-resume` reuses every step whose inputs and settings are unchanged.
Startup validation of `config.yaml` and `samplesheet.csv` runs before any job and stops with a message if something is wrong.

## samplesheet.csv

One row per library. The header is

```
sample_id,fastq_1,fastq_2,bam,exclude,condition
```

| Column | Required | Meaning |
|---|---|---|
| `sample_id` | yes | Unique library name, used verbatim in every file name. Must match `^[A-Za-z0-9._-]+$` (letters, digits, `.`, `_`, `-`). Duplicates are a startup error. The sample id has no read suffix (no `.R1`). |
| `fastq_1`, `fastq_2` | for FASTQ input | Absolute paths of the R1 and R2 FASTQ. Paired-end only: a row with one of the two and no `bam` is a startup error. |
| `bam` | for BAM input | Absolute path of a BAM. It is converted to paired FASTQ first (BAM_TO_FASTQ). A row needs either a FASTQ pair or a BAM. Both together is a startup error, and so is neither. |
| `exclude` | no | `true` drops the row from the whole run. Blank or `false` keeps it. Excluded rows stay in the file; `results/samplesheet.tsv` holds only the kept rows. |
| any other column | no | Free metadata, for example `condition`, `cell_line`, `sgRNA`, `treatment`. The names are checked only for being unique. The columns are carried into `results/samplesheet.tsv` for later differential expression. |

Rules: at least one row must not be excluded; every listed file must exist; each column header must appear once.

**Comments are not allowed in the CSV**: a `#` line is read as a row and fails the `sample_id` check. Put any notes in this README or in `config.yaml`
(which does allow `#` comments), not in `samplesheet.csv`.

## config.yaml

A flat `key: value` file. Each key has a default in `nextflow.config`; the keys below are the full list. Values are read by
Nextflow as parameters (`-params-file config.yaml`), so the names must be exact. An unknown key stops the run at startup.

### Parameters

| Key | Default | Meaning |
|---|---|---|
| `project_name` | `my_rnaseq` | Prefix of the counts table (`<project_name>_counts.txt`) and of the LSF job name. Set it. |
| `samplesheet` | `samplesheet.csv` | The library table. |
| `fastq_dir` | `null` | A directory of FASTQs. On the first run the launcher moves its `*.fastq.gz` and `*.fastq` files into `original_data/fastq/` and rewrites this key. `null`: samplesheet paths are used as given. |
| `outdir` | `./results` | Where published outputs go. |
| `metadata_dir` | `./metadata` | Run bookkeeping: `pipeline_info/` (trace, report, timeline, DAG) and the head job's logs. |
| `intermediate_dir` | `./intermediate_files` | Trimmed FASTQs, FASTQs made from BAMs, aligned and filtered BAMs. |
| `publish_mode` | `copy` | How files are published into `outdir` (Nextflow `publishDir` mode). |
| `lsf_exclude_hosts` | `''` | Hosts that compute jobs must avoid (LSF). Empty = none. |
| `hisat2_index` | `null` | HISAT2 index **prefix**, without `.N.ht2`, for example `/project/rephilllab/references/hisat2/t2t/chm13-2`. Required. Startup fails if no `<prefix>.1.ht2` (or `.1.ht2l`) exists. |
| `gtf` | `null` | Gene annotation GTF, used by featureCounts and by the strandedness check (gene model). Required; startup fails if the file is missing. |
| `strandedness` | `unstranded` | Library strandedness used for counting: `unstranded` (featureCounts `-s 0`, legacy GG), `forward` (`-s 1`), `reverse` (`-s 2`). Any other value is a startup error. The automatic check does not change this value (see Strandedness). |
| `strand_call_threshold` | `0.8` | The check calls a sample forward or reverse when at least this fraction of the informative reads fits that orientation. Must be in (0.5, 1]. |
| `infer_experiment_reads` | `200000` | Reads sampled per BAM by RSeQC `infer_experiment.py` (`-s`). Positive integer. |
| `run_fastqc` | `true` | FastQC on the raw reads and on the trimmed reads, both mates. Legacy GG ran it on trimmed R1 only. |
| `trim` | `true` | Trimmomatic paired-end trimming. `false`: the reads pass through unchanged. |
| `trim_adapter_set` | `TruSeq3-PE-2` | Adapter file name (without `.fa`) in the Trimmomatic container. |
| `trim_illuminaclip` | `'2:30:10:8:True'` | ILLUMINACLIP settings (seed mismatches, palindrome clip threshold, simple clip threshold, minimum adapter length, keep both reads). |
| `trim_steps` | `''` | Extra Trimmomatic steps, appended after ILLUMINACLIP (for example `'LEADING:3 TRAILING:3'`). Empty = none (legacy GG). |
| `trim_minlen` | `25` | Trimmomatic `MINLEN`. |
| `hisat2_args` | `''` | Extra HISAT2 arguments, passed in addition to the pipeline's own (`-p`, `-q`, `--dta`, `--no-discordant`, `--no-mixed`, the index and the summary file). Empty = legacy GG. |
| `alignment_min_pct` | `50` | Overall HISAT2 alignment rate below which ALIGNMENT_SUMMARY prints a warning. Warns only, never fails. |
| `mapq` | `30` | MAPQ filter (samtools `-q`). `0` = no MAPQ filter. Must be a non-negative integer. |
| `mark_duplicates` | `true` | Picard MarkDuplicates on the filtered BAM; the marked BAM is the one counted. `false`: the filtered BAM is the final BAM and there are no duplicate metrics. |
| `remove_duplicates` | `false` | `true` removes duplicates instead of marking them. Needs `mark_duplicates: true`. Legacy GG: `false` (mark only). |
| `optical_duplicate_pixel_distance` | `2500` | Picard `OPTICAL_DUPLICATE_PIXEL_DISTANCE`. |
| `remove_chroms` | `''` | Comma list of contigs to drop before counting. Empty = keep all. Legacy GG removed `chrM`: no effect on the T2T counts, but it zeroes the mitochondrial genes on mm39 (keep this empty for mm39). A contig that is not in the BAM is a warning; if no contig is left the run stops (exit 3). |
| `featurecounts_feature_type` | `exon` | GTF feature type counted (featureCounts `-t`). |
| `featurecounts_attribute` | `gene_id` | GTF attribute the counts are summed by (featureCounts `-g`). |
| `featurecounts_args` | `''` | Extra featureCounts arguments, appended after the pipeline's own. |
| `run_multiqc` | `true` | MultiQC report over all QC outputs. |
| `help` | `false` | `true` prints the parameter help and stops. |

`run_timestamp` is not a key: the pipeline sets it at launch, and it names the trace, report, timeline and DAG files in
`metadata/pipeline_info/` (`execution_trace_<timestamp>.txt` and so on), so a resumed run does not overwrite the previous one's files.

### Script replacement keys

Each pipeline script can be replaced for one project by a `config.yaml` key named after the script (its file name without the
extension). Its value is the path of the replacement, relative to the project directory. Extra arguments for a script go in
`<script>_args`. These keys are not declared in `nextflow.config`; they are read by the pipeline itself.

| Key | Script it replaces | Extra arguments |
|---|---|---|
| `bam_to_fastq` | `bin/bam_to_fastq.sh` | `bam_to_fastq_args` |
| `fastqc` | `bin/fastqc.sh` | `fastqc_args` |
| `trim` | `bin/trim.sh` | `trim_args` (a boolean `trim` is the switch, and is only read as a path when it is a string) |
| `align` | `bin/align.sh` | `align_args` |
| `alignment_summary` | `bin/alignment_summary.py` | `alignment_summary_args` |
| `bam_filter` | `bin/bam_filter.sh` | `bam_filter_args` |
| `mark_duplicates` | `bin/mark_duplicates.sh` | `mark_duplicates_args` (a boolean `mark_duplicates` is the switch) |
| `gtf_to_genepred` | `bin/gtf_to_genepred.sh` | `gtf_to_genepred_args` |
| `genepred_to_bed12` | `bin/genepred_to_bed12.sh` | `genepred_to_bed12_args` |
| `infer_strandedness` | `bin/infer_strandedness.sh` | `infer_strandedness_args` |
| `qc_strandedness_table` | `bin/qc_strandedness_table.py` | `qc_strandedness_table_args` |
| `featurecounts` | `bin/featurecounts.sh` | `featurecounts_args` (the declared parameter, read by this script) |
| `counts_merge` | `bin/counts_merge.py` | `counts_merge_args` |
| `samplesheet_tsv` | `bin/samplesheet_tsv.py` | `samplesheet_tsv_args` |
| `multiqc` | `bin/multiqc.sh` | `multiqc_args` |
| `deg_common` | `bin/deg_common.R` (helper sourced by every `deg_*` script, passed to each DEG process) | none |
| `deg_prepare` | `bin/deg_prepare.R` | `deg_prepare_args` |
| `deg_qc` | `bin/deg_qc.R` | `deg_qc_args` |
| `deg_deseq2` | `bin/deg_deseq2.R` | `deg_deseq2_args` |
| `deg_results` | `bin/deg_results.R` | `deg_results_args` |
| `deg_transform` | `bin/deg_transform.R` | `deg_transform_args` (the parameter is `deg_transform_method`) |
| `deg_volcano` | `bin/deg_volcano.R` | `deg_volcano_args` |
| `deg_heatmap` | `bin/deg_heatmap.R` | `deg_heatmap_args` (the parameter is `deg_heatmaps`, with an s) |
| `deg_ma` | `bin/deg_ma.R` | `deg_ma_args` |
| `deg_gene_plot` | `bin/deg_gene_plot.R` | `deg_gene_plot_args` |

A key that looks like a script path (ends in `.sh`, `.py` or `.R`) but names no script is a startup error that lists the valid names.

## Process list

Steps in order. Each one is a process in `modules/local/` and a script in `bin/` (the script is the step's entire work).

| Process | Module | Script | Container | What it does |
|---|---|---|---|---|
| BAM_TO_FASTQ | `bam_to_fastq.nf` | `bam_to_fastq.sh` | samtools 1.23 | BAM rows only: `samtools collate` then `samtools fastq` into `S.R1.fastq.gz` and `S.R2.fastq.gz`. |
| FASTQC_RAW, FASTQC_TRIM | `fastqc.nf` | `fastqc.sh` | fastqc 0.12.1 | FastQC on both mates, raw and trimmed. |
| TRIM | `trim.nf` | `trim.sh` | trimmomatic 0.40 | Paired-end trimming; Trimmomatic's unpaired outputs are discarded. |
| ALIGN | `align.nf` | `align.sh` | `rnaseq-align_1.0.sif` (lab image) | HISAT2 2.2.1, piped into `samtools view -b` and `samtools sort`, then `samtools index`. No SAM on disk. |
| ALIGNMENT_SUMMARY | `alignment_summary.nf` | `alignment_summary.py` | cutandrun-py_1.2.sif | Overall alignment rate per sample from the HISAT2 summaries; warns below `alignment_min_pct`. |
| BAM_FILTER | `bam_filter.nf` | `bam_filter.sh` | samtools 1.23 | Raw flagstat and idxstats; MAPQ and flag filter; optional contig removal; sort, index; filtered flagstat and idxstats. |
| MARK_DUPLICATES | `mark_duplicates.nf` | `mark_duplicates.sh` | picard 2.27.5 | Picard MarkDuplicates (or removal), then `picard BuildBamIndex`. With `mark_duplicates: false`, copies the filtered BAM. |
| GTF_TO_GENEPRED | `gtf_to_bed12.nf` | `gtf_to_genepred.sh` | ucsc-gtftogenepred 482 | GTF to genePred (once per run). |
| GENEPRED_TO_BED12 | `gtf_to_bed12.nf` | `genepred_to_bed12.sh` | ucsc-genepredtobed 482 | genePred to BED12, the gene model for RSeQC (once per run). |
| INFER_STRANDEDNESS | `infer_strandedness.nf` | `infer_strandedness.sh` | rseqc 5.0.6 | RSeQC `infer_experiment.py` per final BAM. |
| QC_STRANDEDNESS_TABLE | `qc_strandedness_table.nf` | `qc_strandedness_table.py` | cutandrun-py_1.2.sif | One table: fractions per sample, the call, the configured value and a mismatch flag. Never fails on a mismatch. |
| FEATURECOUNTS | `featurecounts.nf` | `featurecounts.sh` | subread 2.1.1 | featureCounts on the final BAM with the configured strandedness; per-sample count table. |
| COUNTS_MERGE | `counts_merge.nf` | `counts_merge.py` | cutandrun-py_1.2.sif | Merged counts table, gene length table and the featureCounts summary table. |
| SAMPLESHEET_TSV | `samplesheet_tsv.nf` | `samplesheet_tsv.py` | cutandrun-py_1.2.sif | Non-excluded samplesheet rows, all columns except `exclude`, as TSV. |
| MULTIQC | `multiqc.nf` | `multiqc.sh` | multiqc 1.35 | MultiQC over every QC output, with sample names cleaned to `sample_id`. |
| DEG_PREPARE | `deg_prepare.nf` | `deg_prepare.R` | `omics-r_1.0.sif` (see Containers) | `run_deg` only. Checks and aligns counts, metadata and comparisons; sets reference levels. |
| DEG_QC | `deg_qc.nf` | `deg_qc.R` | same | Reads-per-sample bar plot, MDS plot, sample correlation heatmap. |
| DEG_DESEQ2 | `deg_deseq2.nf` | `deg_deseq2.R` | same | One DESeq2 fit per model. |
| DEG_RESULTS | `deg_results.nf` | `deg_results.R` | same | Results tables, shrunk fold changes, analysis summary, normalised counts (per model). |
| DEG_TRANSFORM | `deg_transform.nf` | `deg_transform.R` | same | rlog / vst counts (per model used by a heatmap or the gene plots). |
| DEG_VOLCANO | `deg_volcano.nf` | `deg_volcano.R` | same | Volcano plots (per model). |
| DEG_HEATMAP | `deg_heatmap.nf` | `deg_heatmap.R` | same | One heatmap per row of `deg_heatmaps`. |
| DEG_MA | `deg_ma.nf` | `deg_ma.R` | same | MA plots (per model; off by default). |
| DEG_GENE_PLOT | `deg_gene_plot.nf` | `deg_gene_plot.R` | same | Expression plots of selected genes (once). |

GTF_TO_BED12 is two processes because a Nextflow process has one container. The step runs once per run.

The DEG processes run only with `run_deg: true` (see Differential expression). With `deg_counts` set, none of the processes above them runs.

## Outputs

```
results/
├── counts/
│   ├── <project_name>_counts.txt        merged counts: header Gene + one column per sample (samplesheet order)
│   ├── gene_length.txt                  Geneid, Strand, Length from the first sample's featureCounts table
│   ├── featurecounts_summary.tsv        Status and one column per sample
│   └── per_sample/                      S_counts.txt, S_featurecounts.txt, S_featurecounts.txt.summary
├── bam/final/                           S.final.bam and S.final.bam.bai (the BAM that was counted)
├── qc/
│   ├── fastqc/raw/, fastqc/trimmed/     FastQC reports, both mates
│   ├── trimmomatic/                     Trimmomatic logs
│   ├── hisat2/                          HISAT2 summaries, alignment_summary.tsv
│   ├── samtools/                        raw and filtered flagstat and idxstats
│   ├── picard/                          Picard duplicate metrics (not with mark_duplicates: false)
│   ├── strandedness/                    S.infer_experiment.txt per sample, strandedness.tsv, genes.bed12
│   └── multiqc/                         multiqc_report.html (start here)
├── samplesheet.tsv                      the kept rows of the samplesheet
├── deg/                                 differential expression (only with run_deg: true; layout under Differential expression)
intermediate_files/
├── fastq_from_bam/                      S.R1.fastq.gz, S.R2.fastq.gz (BAM rows)
├── trimmed/                             S.R1.trim.fastq.gz, S.R2.trim.fastq.gz
├── aligned/                             S.aligned.bam and .bai
└── filtered/                            S.filtered.bam and .bai
metadata/
├── nf_logs/                             the head job's stdout and stderr
└── pipeline_info/                       execution trace, report, timeline and DAG (names carry the launch timestamp)
```

Each published directory has a `log/` subfolder with the exact command that ran for each step, as `<sample>_<step>_command.sh` (a step without a sample uses `<step>_command.sh`).
A step that failed leaves its command only in its Nextflow work directory, which the error message names.

## Differential expression

The pipeline can also run the DESeq2 differential expression analysis of the legacy `Scripts/step_6_deg_analysis.Rmd`. It is **off by
default** (`run_deg: false`): with the default nothing changes and the pipeline ends at the counts table. Everything below is
about `run_deg: true`.

### Quick start

```yaml
# config.yaml, on top of the usual keys
run_deg: true
deg_comparisons: deg_comparisons.csv      # required: what to compare (format below)
deg_heatmaps: deg_heatmaps.csv            # optional: heatmaps of the DEGs
deg_colors: deg_colors.csv                # optional: group colours
deg_label_rules: deg_label_rules.csv      # optional: nicer sample names in the heatmaps
deg_reference_levels: treatment=LUC_DMSO  # the reference level of each design factor
```

The sample metadata (the columns the designs use, for example `treatment`) comes from `deg_metadata`, or from `samplesheet.csv` when
that is not set (rows with `exclude: true` are dropped). The count table is this run's `results/counts/<project_name>_counts.txt`.
The analysis starts after COUNTS_MERGE, so the reads are aligned and counted first.

### From an existing counts table

To run only the DEG steps on a counts table you already have, set `deg_counts` (and `run_deg: true`):

```yaml
run_deg: true
deg_counts: /path/to/earlier_run/results/counts/my_project_counts.txt   # Gene + one column per sample
deg_metadata: deg_metadata.csv          # or leave it out to use samplesheet.csv (it then needs only sample_id + the metadata columns)
deg_comparisons: deg_comparisons.csv
```

No FASTQ, BAM, `hisat2_index`, `gtf`, QC or MultiQC is needed or checked, and none of the upstream steps runs. `deg_counts` without
`run_deg: true` is a startup error. The column names of the counts table must equal the sample names in the metadata; a suffix such
as the legacy `.R1` is removed with `deg_strip_suffix: '\.R1$'`.

### What runs

| Step | Runs | What it does |
|---|---|---|
| DEG_PREPARE | once | Checks the counts table, the metadata and the comparisons table, aligns the samples, sets the reference levels, writes the cleaned copies. |
| DEG_QC | once (`deg_run_qc`) | Reads per sample (bar plot), MDS plot (`pca.pdf`, it is an MDS plot), sample correlation heatmap. |
| DEG_DESEQ2 | per model | One DESeq2 fit per model: all rows of `deg_comparisons` with the same `model` share it. |
| DEG_RESULTS | per model | A results table for every comparison of the model, shrunk fold changes, an analysis summary, normalised counts. |
| DEG_TRANSFORM | per model that needs it | rlog (or vst) counts for the heatmaps and gene plots. Slow, so only models used by a heatmap (and the first model for gene plots) are transformed. |
| DEG_VOLCANO | per model (`deg_run_volcano`) | Volcano plots: one SVG per comparison (and a shrunk one), three multi-page PDFs, and a combined figure. |
| DEG_MA | per model (`deg_run_ma`, off) | MA plots. |
| DEG_HEATMAP | per row of `deg_heatmaps` | Heatmap of the DEGs of one comparison. |
| DEG_GENE_PLOT | once (`deg_gene_plot_genes`) | One expression plot per gene, from the first model's transformed counts. |

### deg_comparisons.csv (required)

One row per comparison. The header names are exact; the order of the columns is free.

| Column | Required | Meaning |
|---|---|---|
| `model` | yes | Name of the model and of its output folder (`results/deg/<model>/`). Letters, digits, `.`, `_`, `+`, `-`. `prepared`, `qc` and `gene_plots` are reserved. Rows with the same `model` share one DESeq2 fit. |
| `design` | yes | The DESeq2 design formula, for example `~treatment` or `~ treatment + batch`. Every row of one model must have the same design. |
| `comparison` | yes | Name of the comparison (unique within its model; same characters as `model`; `log` is reserved). It names the folder and the files `results_<comparison>.csv`. |
| `factor`, `numerator`, `denominator` | one contrast kind | A contrast between two levels of a design factor, for example `treatment`, `CBX2_1uM`, `LUC_DMSO`: log2 fold change of numerator over denominator. All three together. |
| `coef` | one contrast kind | The name of one model coefficient (as `resultsNames(dds)` lists it), for example `drug_1uM_vs_DMSO` (design `~ drug`). |
| `add`, `subtract` | one contrast kind | A sum of coefficients: names separated by `;`, for example `treatment_B_vs_A;treatment_C_vs_A` (design `~ treatment`, sum of two coefficients). `subtract` (optional) lists the coefficients that are subtracted. |
| `label` | no | Plot title of the comparison. Blank = the comparison name with `_` replaced by spaces. |
| `combined_order`, `combined_title` | no | A number puts the comparison in the combined volcano figure, in that order, with that title. The figure is made when at least two comparisons of a model have a number. |
| `xlim`, `ylim` | no | Volcano axis limits of this comparison as `min,max` (put the cell in quotes: `"-3,3"`). They win over `deg_volcano_xlim` / `deg_volcano_ylim`. |

Each row needs exactly one contrast kind (a row with two kinds, or none, stops the run at startup). `deg_shrink_type: apeglm` works with
`coef` rows only.

```
model,design,comparison,factor,numerator,denominator,coef,add,subtract,label,combined_order,combined_title,xlim,ylim
treatment,~treatment,CBX2_1uM_vs_LUC_DMSO,treatment,CBX2_1uM,LUC_DMSO,,,,CBX2 KO + EZH2i,3,CBX2 KO + EZH2i,,
```

### deg_heatmaps.csv (optional)

One row per heatmap. A heatmap shows the genes that are differentially expressed (not `NS`) in one comparison, from the model's
transformed counts, scaled by row.

| Column | Required | Meaning |
|---|---|---|
| `name` | yes | Unique name; the file is `results/deg/<model>/heatmap_<name>.svg`. |
| `model`, `comparison` | yes | Which results table to use. The pair must exist in `deg_comparisons`. |
| `gene_set` | no | A text file with one gene per line: only DEGs in the set are drawn. |
| `clust_method` | no | `hclust` method; blank = `ward.D2`. |
| `cutree_rows` | no | Number of row clusters to cut; blank = none. |
| `col_swap` | no | `none` (default), `scaled` or `unscaled`: swap the two main column branches (clustering on the scaled or unscaled matrix). |
| `palette_n` | no | Number of colours; blank = 1000. |
| `fontsize`, `fontsize_row`, `width`, `height` | no | Blank = 9, 6.5, 6.5, 4.5 (inches). |
| `title` | no | `{n}` is replaced by the number of genes. Blank = `DEGs: <label> (N = {n})`. |

### deg_colors.csv and deg_label_rules.csv (optional)

`deg_colors.csv` has the columns `group,color` (a group is a level of the colour-by column, a colour a name or a hex code). A group
without a row gets a default colour.

`deg_label_rules.csv` has the columns `pattern,replacement`. The rules are applied in order, as regular expressions, to the sample names
shown in the heatmaps (the same as the legacy chain of `gsub()` calls). A typed `\n` in the replacement is a line break. Quote
cells with leading or trailing spaces, and write an empty replacement as `""`:

```
pattern,replacement
1uM,\nEZH2i
A_," A"
DIPGVI_,""
```

### Parameters

All default to the legacy Rmd (except the colours). The same list with comments is in `examples/config.yaml`. Hex colours need quotes
in YAML (`'#46627A'`): an unquoted `#` starts a comment.

| Key | Default | Meaning |
|---|---|---|
| `run_deg` | `false` | Master switch. |
| `deg_counts` | `null` | Existing counts table: skips all upstream steps (needs `run_deg: true`). |
| `deg_metadata` | `null` | Sample metadata CSV/TSV; `null` = the samplesheet. Must contain the `deg_sample_col` column. |
| `deg_comparisons` | `null` | Comparisons CSV; required with `run_deg`. |
| `deg_heatmaps`, `deg_colors`, `deg_label_rules` | `null` | Optional CSVs (above). |
| `deg_sample_col` | `sample_id` | Metadata column with the sample names. |
| `deg_strip_suffix` | `''` | Regular expression removed from the count column names (and the metadata sample names), for example `'\.R1$'`. |
| `deg_reference_levels` | `''` | `factor=level,factor=level`. Factors not listed keep the alphabetical first level (the run prints which). |
| `deg_min_count`, `deg_min_samples` | `10`, `2` | Keep genes with at least `deg_min_count` reads in at least `deg_min_samples` samples. Positive integers. |
| `deg_fit_type` | `parametric` | DESeq2 dispersion fit: `parametric`, `local` or `mean`. Used by DEG_DESEQ2 and DEG_TRANSFORM. |
| `deg_alpha` | `0.05` | FDR cut-off, between 0 and 1 (not included). |
| `deg_lfc` | `1` | Absolute log2 fold change cut-off, 0 or more. |
| `deg_lfc_threshold_in_results` | `false` | `true` passes `lfcThreshold` to `results()`. |
| `deg_shrink_type` | `ashr` | `ashr`, `apeglm` (coef rows only) or `none`. |
| `deg_sval_threshold` | `0.005` | s-value cut-off of the shrunk effect. |
| `deg_transform_method`, `deg_transform_blind` | `rlog`, `false` | Transformation for heatmaps and gene plots (`rlog` or `vst`). |
| `deg_color_down`, `deg_color_up`, `deg_color_ns`, `deg_color_mid` | `'#46627A'`, `'#8c6bb1'`, `gray`, `white` | Volcano / MA / heatmap colours (the heatmap's high end is `deg_heatmap_color_up`). |
| `deg_heatmap_color_up` | `'#88419d'` | DEG heatmap high z-score colour. |
| `deg_color_by` | `''` | Metadata column that colours the QC plots; blank = the first design variable of the first model. |
| `deg_mds_input` | `logcpm` | `logcpm` (log2 CPM, prior count 2) or `counts` (raw counts as the legacy MDS, dominated by library size). |
| `deg_mds_label_by`, `deg_mds_legend_title` | `''` | Point labels of the MDS plot; legend title (blank = the colour-by column). |
| `deg_legacy_mds_legend_col` | `''` | A metadata column name reproduces the legacy MDS legend. |
| `deg_cor_method` | `spearman` | `pearson`, `spearman` or `kendall`. The heatmap file is `<method>_correlation_heatmap.svg`. |
| `deg_label_genes` | `''` | Comma list of genes named in the volcano plots. |
| `deg_volcano_xlim`, `deg_volcano_ylim` | `''` | `min,max` of every volcano; blank = automatic. |
| `deg_volcano_common_xlim`, `deg_volcano_common_ylim` | `''` | Limits of the combined figure and the scaled PDFs. |
| `deg_volcano_combined_layout` | `stacked` | `stacked` or `side`. |
| `deg_gene_plot_genes` | `''` | Comma list of genes to plot; blank = no gene plots. Genes missing from the table are skipped. |
| `deg_gene_plot_x` | `''` | x axis and colour column of the gene plots; blank = `deg_color_by`, else the first design variable. |
| `deg_legacy_summary_counts` | `false` | `true` reproduces the legacy analysis summary (always 0 up / 0 down). |
| `deg_run_qc`, `deg_run_barplot`, `deg_run_mds`, `deg_run_correlation` | `true` | Switch the QC step and its three plots. |
| `deg_run_volcano`, `deg_run_volcano_pdfs`, `deg_run_volcano_combined` | `true` | Volcano plots, the three PDFs, the combined figure. |
| `deg_run_heatmap` | `true` | Heatmaps (also need `deg_heatmaps`). |
| `deg_run_ma` | `false` | MA plots. |
| `deg_run_gene_plots` | `true` | Gene plots (run only when `deg_gene_plot_genes` is set). |

One more legacy switch of the R scripts is not a config key; pass it with the `_args` key: `deg_qc_args: "--legacy-color-index=TRUE"`
(reproduce the legacy QC colours).

### Custom statistical tests (not shipped)

The pipeline runs per-model Wald tests (design, `coef` and `add` / `subtract` contrasts). An interaction model, or any other custom
statistical test, is not part of the pipeline: write it as a custom script in the project's `utils/` folder. It can read
`results/deg/prepared/counts_edited.rds` and `results/deg/prepared/sample_metadata_edited.rds`. On the Docker versions such a script
is run with `run_rnaseq script utils/X.R`, as the CUT&RUN Docker launcher does with `script`.

### Outputs

```
results/deg/
├── prepared/        counts_edited.rds/.tsv, sample_metadata_edited.rds/.csv, comparisons.tsv, prepare_summary.txt
├── qc/              counts_per_index.pdf, counts_per_sample.tsv, pca.pdf (MDS), <method>_correlation_heatmap.svg
├── <model>/         DESeqDataSet.rds, DESeq.rds, design.txt, dispersion_plots.svg, prefilter.tsv,
│                    analysis_summary.csv/.txt, normalized_counts.csv, rlogData.rds, rlogCounts.rds,
│                    volcano_plots_{yscaled,xscaled,freescale}.pdf, volcano_combined_<n>_comparisons_<layout>.svg,
│                    heatmap_<name>.svg, log/ (command.sh of every step, <step>_command.sh)
│   └── <comparison>/  results_<comparison>.csv, volcano_all_genes_<comparison>.svg, volcano_shrunk_all_genes_<comparison>.svg, ma_<comparison>.svg
└── gene_plots/      <gene>.svg
```

`results_<comparison>.csv` has the columns `Gene, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj, shrunkLFC, svalue, effect,
shrunk_effect` (the shrunk columns are missing with `deg_shrink_type: none`); `effect` is `Upregulated`, `Downregulated` or `NS`.

The DEG steps need a container with R and DESeq2 (see Containers).

## Strandedness

`strandedness` sets how featureCounts counts: `unstranded` = `-s 0`, `forward` = `-s 1`, `reverse` = `-s 2`. The default
`unstranded` is what the legacy GG counts used.

The automatic check runs RSeQC `infer_experiment.py` on `infer_experiment_reads` reads of each final BAM, against a gene model
made from the GTF. For each sample, `strandedness.tsv` reports:

- `fraction_forward`: the fraction of informative reads explained by `1++,1--,2+-,2-+`;
- `fraction_reverse`: the fraction explained by `1+-,1-+,2++,2--`;
- `fraction_undetermined`: the rest;
- `call`: `forward` if `fraction_forward` >= `strand_call_threshold`, `reverse` if `fraction_reverse` >= the threshold, else `unstranded`;
- `configured`: the `strandedness` value of the run;
- `mismatch`: `true` when `call` differs from `configured`.

A mismatch prints a `WARN` line at the end of the run. **The check never changes the counting.** Only the `strandedness` key does.

## Containers

Images are cached in `/project/rephilllab/nf_containers` (shared by every lab pipeline; the `.img` and `.sif` files are reused, not
pulled again). Tags are pinned.

| Tool | Image | State |
|---|---|---|
| samtools 1.23 | `quay.io/biocontainers/samtools:1.23--h96c455f_0` | cached |
| fastqc 0.12.1 | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` | cached |
| trimmomatic 0.40 | `quay.io/biocontainers/trimmomatic:0.40--hdfd78af_0` | cached |
| picard 2.27.5 | `quay.io/biocontainers/picard:2.27.5--hdfd78af_0` | cached (no samtools inside) |
| multiqc 1.35 | `quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1` | cached |
| pandas scripts | `/project/rephilllab/nf_containers/cutandrun-py_1.2.sif` | cached (shared with CUT&RUN) |
| HISAT2 2.2.1 + samtools 1.23 | `/project/rephilllab/nf_containers/rnaseq-align_1.0.sif` | built from the recipe under `containers/` (`hisat2==2.2.1`, `samtools==1.23`); ALIGN cannot run until the image exists in `nf_containers` |
| R + DESeq2 (DEG_* steps) | `/project/rephilllab/nf_containers/omics-r_1.0.sif` | built from the R env lock of the omics Docker image, so the R and package versions equal Docker 2026.10.3. It must exist in `nf_containers` before the first DEG run |
| subread 2.1.1 | `quay.io/biocontainers/subread:2.1.1--h577a1d6_0` | pulled on first use |
| RSeQC 5.0.6 | `quay.io/biocontainers/rseqc:5.0.6--pyhdfd78af_0` | pulled on first use |
| UCSC gtfToGenePred / genePredToBed 482 | `quay.io/biocontainers/ucsc-gtftogenepred:482--h0b57e2e_0`, `quay.io/biocontainers/ucsc-genepredtobed:482--h0b57e2e_0` | pulled on first use |

The Java heap of Trimmomatic and Picard is set from the task's memory request. All steps run in
containers through Apptainer with `-B /usr/lib64/ossl-modules` (the PMACS FIPS fix).

## Customising a project

To change one step for one project, copy its script into the project's `utils/`, edit the copy, and add the key from the table
above with the copy's path:

```bash
cp /project/rephilllab/RNASEQ/bin/featurecounts.sh utils/
# edit utils/featurecounts.sh, keeping its header's arguments and outputs
# in config.yaml:   featurecounts: utils/featurecounts.sh
run_rnaseq.sh      # only that step and the steps after it re-run
```

Keep the script's header contract: its arguments and every output file name. The replacement must write the same outputs as the
shipped script, or the downstream steps fail.

## Files

```
run_pipeline.sh              launcher: init and run
main.nf, nextflow.config     the pipeline and its parameters
conf/                        base.config (resources, retries), pmacs.config (LSF, Apptainer), test.config
modules/local/               one process per file (functions.nf: runner, extraArgs; deg_functions.nf: shell quoting helpers of the DEG modules)
bin/                         one script per step
lib/lsf_watchdog.sh          stops a run whose LSF job was aborted at dispatch
containers/                  recipe for rnaseq-align_1.0
examples/                    config.yaml and samplesheet.csv templates, project.gitignore, projects/
agent/                       AGENTS.md and CLAUDE.md copied into each project by init
```
