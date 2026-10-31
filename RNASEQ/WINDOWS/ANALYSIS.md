# RNA-seq: what the analysis does

This file describes the analysis itself. It applies to the cluster version of the pipeline and to the Docker versions
(Windows and Mac); only how each version is launched and where its software comes from differ. Inputs, config keys and
commands are in each version's `README.md`. The script of each step is in `bin/`; its header documents its arguments and outputs.

The pipeline counts reads per gene for each library and merges the counts into one table. Optional DEG steps
(`run_deg: true`) then run DESeq2 on that table.

## Software versions to cite

Take versions from where the run actually ran. The Docker image's tool versions may differ from the cluster pins.

| Tool | Cluster version (pinned) | Docker versions |
|---|---|---|
| Image | lab container images (one per tool) | `ghcr.io/genesig/omics:<tag>` (tag in `docker/image.txt`), one image for every tool |
| HISAT2 | 2.2.1 | as printed by `run_rnaseq setup` |
| samtools | 1.23 | same |
| Trimmomatic | 0.40 | same |
| FastQC | 0.12.1 | same |
| Picard | 2.27.5 | same |
| featureCounts (Subread) | 2.1.1 | same |
| RSeQC (`infer_experiment.py`) | 5.0.6 | same |
| UCSC gtfToGenePred, genePredToBed | 482 | same |
| MultiQC | 1.35 | same |
| Nextflow | 26.04.6 | same |
| DESeq2 | see the DEG image | 1.46 (image 2026.10.3) |

## What each step does

| Process | Tool | What it does |
|---|---|---|
| BAM_TO_FASTQ | samtools | BAM input rows only: `collate`, then `fastq` into paired FASTQs. Stops when no read pairs come out. |
| FASTQC_RAW, FASTQC_TRIM | FastQC | Quality reports for both mates, on the input and on the trimmed reads. |
| TRIM | Trimmomatic PE | ILLUMINACLIP (`trim_adapter_set`, `trim_illuminaclip`), optional `trim_steps`, `MINLEN` = `trim_minlen`. |
| ALIGN | HISAT2 | `-p -q --dta --no-discordant --no-mixed`, index `hisat2_index`, plus `hisat2_args`; piped to samtools view and sort. Unmapped reads are kept. |
| ALIGNMENT_SUMMARY | script | Overall alignment rate per sample; warning below `alignment_min_pct`, never fails. |
| BAM_FILTER | samtools | Flagstat/idxstats before and after; keeps MAPQ >= `mapq`, properly paired, primary, mapped, non-duplicate reads; drops `remove_chroms`. |
| MARK_DUPLICATES | Picard | Marks (or with `remove_duplicates: true` removes) duplicates; the result is the BAM that is counted. |
| GTF_TO_GENEPRED, GENEPRED_TO_BED12 | UCSC tools | Gene model (BED12) from the GTF, for the strandedness check. |
| INFER_STRANDEDNESS, QC_STRANDEDNESS_TABLE | RSeQC | Strandedness call per library; a mismatch with `strandedness` is a warning only. |
| FEATURECOUNTS | featureCounts | Per-gene read-pair counts: `-p --countReadPairs -B -C --primary`, `-s` from `strandedness`. |
| COUNTS_MERGE | script | Merged counts table, gene lengths, featureCounts summary. |
| SAMPLESHEET_TSV | script | The kept samplesheet rows. |
| MULTIQC | MultiQC | One report over every QC output. |
| DEG_* | DESeq2 (R) | Optional: prepare, QC (PCA, MDS on logCPM, correlation), DESeq2 per model, results, rlog transform, volcano plots, heatmaps, gene plots. Defaults: see README, "Differential expression". |

## Counting and strandedness

Counts are read pairs; duplicates are marked, not removed, unless `remove_duplicates: true`. `strandedness` (default
`unstranded`) is what counting uses; the RSeQC check never changes it.

## HISAT2 and the thread count

HISAT2 results depend slightly on the number of threads (`-p`): on the 1M-pair benchmark, 2 threads and 8 threads
differ in about 100 of 1,000,000 pairs, which changes a few multi-mapping genes. ALIGN uses the budget's CPUs less 2 (at
most 8). To reproduce the cluster benchmark byte for byte, run with `max_cpus` of 10 or more (8 ALIGN threads).

## Output layout

`results/counts/<project>_counts.txt`, `gene_length.txt`, `featurecounts_summary.tsv`, `per_sample/`; `results/bam/final/`;
`results/qc/{fastqc,trimmomatic,hisat2,samtools,picard,strandedness,multiqc}/`; `results/samplesheet.tsv`;
`results/deg/{prepared,qc,<model>}/` with `log/` folders holding every step's exact command; `metadata/pipeline_info/`.

## Re-running only the plots

The work volume (Docker) or work folder (cluster) is kept after a successful run (`keep_work: true`, the default). Edit
`config.yaml` or the `deg_*.csv` files and run again: Nextflow resumes and reruns only the steps whose inputs changed (for
example only the volcano or heatmap steps). Only `run_rnaseq clean` (or `keep_work: false`) removes the cache. If the cache
is gone, use DEG-only mode: `run_deg: true` and `deg_counts: <counts table>`, which needs no FASTQs or references.

## Reproducibility

Every step's exact command is in its `log/` folder; keep `config.yaml`, the samplesheet and the `deg_*.csv` files with the
results. `metadata/pipeline_info/` holds the execution trace, report and timeline.
