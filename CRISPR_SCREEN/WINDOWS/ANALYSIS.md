# CRISPR_SCREEN — what the analysis does

This file describes the analysis itself. It applies to all three versions of the pipeline (cluster `HPC/`, Windows `WINDOWS/`, Mac `MAC/`); only how each version is launched differs.

## Materials and methods: METHODS.md

Every successful run ends by writing `METHODS.md` in the project directory
(and a copy in `results/`, or your `outdir`): a Materials and Methods section
for the computational analysis, ready to adapt for a paper, with references,
the exact tool commands (Appendix A) and the MD5 of every input table
(Appendix B). It describes **only the steps the run performed** — no
demultiplexing paragraph for already-demultiplexed input, no statistics
paragraph without `contrasts` — and every parameter in it comes from the run
itself, not from the pipeline's defaults:

- the run's resolved configuration (`metadata/pipeline_info/run_info.json`):
  input mode, samples, which samples were counted against which library,
  anchor discovery or a fixed `trim5`, replaced scripts (path + MD5) and
  `<script>_args`;
- what the steps recorded: cutadapt and MAGeCK log their exact command line
  and version, so a replaced script or extra arguments are described as they
  actually ran; the anchors found for each library; the tool versions probed
  inside each container (`metadata/pipeline_info/software_versions.tsv`);
- what happened to the data: samples left out of a library's count table,
  MAGeCK switching from median to total-count normalisation for a library
  (it does this on its own when a sample has too many zero-count guides),
  which samples MAGeCK used to estimate variance in each comparison, offsets
  chosen by `trim5: AUTO`, duplicated guide sequences.

It covers the computational analysis only: describe library cloning, cell
culture, gDNA extraction, PCR and sequencing yourself. Anything the pipeline
cannot know (for example how already-demultiplexed reads were demultiplexed,
or what a replaced script does) is marked **[To complete: …]**.

`METHODS.md` is written by the last step of the run, so a failed run never
overwrites it; its first lines name the run, configuration file and outdir it
describes. A project running several configurations (each with its own
`outdir`) keeps one copy per `outdir`; the copy in the project directory is
the most recent successful run's.

## What each step does

| Step | Tool | Container |
|---|---|---|
| Demultiplex | cutadapt (`-g NAME=BARCODE`, exact match) | `cutadapt:5.2` |
| QC | FastQC per sample, MultiQC summary, cutadapt-log metrics, barcode-position density (see below) | `fastqc:0.12.1`, `multiqc:1.34` |
| Guide position | Automatic anchor discovery + trim (see below), unless `trim5` is set | `cutadapt:5.2` |
| Count | `mageck count` | `mageck:0.5.9.5` |
| Count QC | Guide representation, Lorenz/Gini, dropout, replicate correlation — computed once into tables, plotted as a separate step, then assembled into one HTML report per library | `rocker/tidyverse:4.3` |
| Test (optional) | `mageck test`, one task per contrast | `mageck:0.5.9.5` |
| Methods | `METHODS.md` from the run's configuration, tool logs and probed tool versions (runs last, after all other steps) | each step's own image (versions); `cutadapt:5.2` (Python) |

Container versions are pinned to the exact same builds as the lab's original
`environment.yml`. With `trim5` set explicitly (bypassing automatic anchor
discovery — see below), counts match byte-for-byte with any prior
non-Nextflow run given the same inputs.

## Guide position: automatic anchor discovery

`mageck count` only ever does exact, fixed-position matching: it takes a
fixed number of bases from the start of each (already demultiplexed) read
and checks that against the library. That fixed number has to be right for
every read, which breaks down the moment the amplicon has any
variable-length stretch — a different UMI, spacer, or synthesis batch —
between the barcode and the guide: reads shift out of position and quietly
fail to match, even though the guide itself was read perfectly.

By default (`trim5` left unset), the pipeline avoids this rather than
asking you to get a fixed number right: it finds the constant vector
sequence between each library's barcode and its guide directly from the
data (`GUIDE_ANCHOR_DISCOVERY` — no prior knowledge of that sequence or its
length needed, just ~200k reads from one of the library's own samples),
then searches for that sequence in every read with real, error-tolerant
alignment (`GUIDE_ANCHOR_TRIM` — cutadapt, not a fixed offset) to position
each guide exactly, regardless of any variable-length spacing before or
after the anchor. `mageck count` then runs on the trimmed reads with
`--trim-5 0`. This is robust to any plasmid design, with nothing to
configure per project.

The anchor is refined against the library's own guides. Consensus extension
alone can't tell a constant vector base from a guide base that is the same
in (nearly) every guide: e.g. a library whose 20-mers all start with the U6
`G` (written as `G` + 19 nt) would get that `G` absorbed into the anchor,
shifting every trimmed read by one base so that nothing maps. Discovery
therefore also tries the anchor with 1, 2, ... (up to 8) bases removed from
its 3' end and keeps the version followed by an exact library guide in the
most sampled reads. `<prefix>_anchor_refinement.csv` lists every candidate
and its score. (For example, `TTGTGGAAAGGACGAAACACCG` for a library of plain
20-mers, but `TTGTGGAAAGGACGAAACACC` for one whose 20-mers include the `G`.)

What discovery assumes: in most reads of the library's first sample, the
anchor starts at the very beginning of the (demultiplexed) read. Variable
spacing *after* the anchor, or a different anchor sequence or length per
library, is handled; a variable-length stagger *before* the anchor (e.g.
staggered primers in already-demultiplexed FASTQs) is not, and discovery
then stops with exit 3 — set `trim5` instead.

`results/counts/qc/anchor_discovery/<prefix>_anchor_discovery.csv` is the
full discovery trace (one row per candidate length, with the fraction of
reads agreeing on it) — useful for confirming what was discovered and why,
or troubleshooting a library where discovery fails (exits 3, with a clear
message, if no anchor at least 8bp long is supported by at least 30% of
sampled reads). `..._anchor_trim_stats.csv` reports, per sample, how many
reads the anchor was actually found in.

Set `trim5` in `config.yaml` to a number to bypass discovery entirely and
use that fixed offset instead, exactly like the pipeline's first version —
useful for reproducing a specific prior run byte-for-byte, or a library
genuinely without a usable constant region.

## Read classification: exact / near-miss / far-miss matches

`mageck count`'s own guide matching is always exact — it has no
mismatch-tolerance option (`mageck count --help` has no such flag), so every
read it counts is, by construction, an exact match to a guide. Mapping Rate
in the QC report already **is** the percent of exact matches; it's just
labeled that way explicitly now.

The pipeline goes one step further and classifies the reads mageck
*couldn't* map (`mageck count --unmapped-to-file`'s output) into:
- **near-miss** — exactly one base different from some guide in the
  library. Most likely lost to sequencing or PCR error; the read almost
  certainly came from a real guide.
- **far-miss** — not within one base of any guide. Contamination,
  adapter-only reads, or a sample that doesn't belong to this library.

This is computed per library in `guide_representation_wrangle.R`
(`<prefix>_read_classification.csv`), plotted per library and — with more
than one library — faceted across pools in the comparison report, so a
sample losing a lot of real signal to sequencing error (high near-miss)
looks very different from one that was never screened with this library at
all (high far-miss).

## Output layout

Every step's own `log/` subfolder below holds both its `.log` (what the
tool said) and a `*_command.sh` (the exact command that ran, fully filled
in) — whatever a labmate or an agent needs to reproduce or troubleshoot it,
without needing Nextflow's own `work/` directories.

```
METHODS.md                      Materials and Methods of the latest successful run (see above)
original_data/fastq/            per-sample FASTQs (demultiplexing output)
metadata/
├── nf_logs/                     run logs, one per run
├── .nextflow.log                Nextflow's own run log
└── pipeline_info/                execution report, timeline, trace, DAG — provenance for every run;
                                  run_info.json (the run's resolved configuration) and
                                  software_versions.tsv (tool versions per container image),
                                  both behind METHODS.md
results/
├── intermediate_files/          tool-format conversions, kept with the results (default
│   │                             intermediate_dir: ./results/intermediate_files)
│   ├── mageck_libraries/        MAGeCK-format library file, converted from library.csv
│   └── anchor_trimmed_fastq/<lib>/  per-sample FASTQs after GUIDE_ANCHOR_TRIM, one folder per
│                                 library (only when trim5 is unset — the input mageck count
│                                 actually counts from)
├── demultiplexing/
│   ├── log/                     cutadapt.log + command.sh
│   └── qc/
│       ├── fastqc/
│       ├── multiqc/
│       ├── cutadapt_stats/      JSON + CSV demux metrics (Gini, CV, per-sample counts)
│       └── demux_plots/         barcode-position density plot (ggplot2 SVG) + raw CSV —
│                                 where in the read the barcode was matched, per sample;
│                                 a shared, sharp peak confirms a consistent read structure
├── counts/                      *.count.txt, *.count_normalized.txt, mageck logs,
│   │                             *.unmapped.txt (every distinct sequence mageck couldn't
│   │                             match to a guide, one column per sample),
│   │                             *.sample_reads.tsv (each sample's input/mapped reads, and
│   │                             whether it was counted or left out with 0 reads)
│   └── qc/
│       ├── anchor_discovery/    <prefix>_anchor_discovery.csv (only when trim5 is unset):
│       │                         the discovery trace — every candidate anchor length tried,
│       │                         and what fraction of reads agreed on it;
│       │                         _anchor_refinement.csv: the candidates checked against the
│       │                         library's guides and which one was kept; and
│       │                         _anchor_trim_stats.csv: per sample, how many reads the
│       │                         anchor was found in.
│       ├── tables/              every number behind the QC plots, as CSVs (Gini, Lorenz,
│       │                         percentile ratio, undetected guides, mapping rate,
│       │                         read classification (exact/near-miss/far-miss), GC vs
│       │                         depth, replicate correlations, dropout, qc_summary);
│       │                         with >1 library also <project_name>_library_breakdown.csv
│       │                         / _wide.csv: % of each sample's reads per library
│       ├── plots/               guide representation / Lorenz / ECDF / ... SVGs, drawn
│       │                         only from tables/
│       └── report/              <prefix>_qc_report.html — one self-contained page per
│                                 library: every plot above, with a description, how to
│                                 interpret it, and what to check if it looks poor. Built
│                                 only from tables/ and plots/, so it never re-runs `mageck
│                                 count`. Start here when triaging a run. With more than one
│                                 library, also <project_name>_pools_comparison_qc_report.html:
│                                 the same plots faceted one panel per pool, to spot the one
│                                 that looks different from the rest.
├── mageck_test/                 one subfolder per contrast, <treatment>_vs_<control> with ^/$ dropped (only if contrasts.csv given)
└── METHODS.md                   this outdir's copy of METHODS.md
```
