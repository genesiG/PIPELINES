# AM_ATACSEQ example (legacy project AM_ATACSEQ)

Files: `config.yaml` (all parameters, defaults = AV_ATACSEQ_2026, AM legacy values as TOGGLE lines), `samplesheet.csv` (9 libraries). No `peak_merges.csv`: the legacy AM scripts did not merge peaks (see Differences, item 11), so `peak_merges: null`. FASTQ paths are absolute (`fastq_dir: null`); the ones in `samplesheet.csv` are the legacy cluster paths, so replace them with the paths of your own FASTQ files before a run. The reference paths in `config.yaml` are `D:/references/...` placeholders.

Provenance: the FASTQ paths in `samplesheet.csv` are the lab's cluster paths; replace them with the local paths of your own FASTQ files.

Dataset: mouse (mm39), paired-end ATAC-seq, cell `fPPP` (metadata column `cell`), 9 `experiment=atacseq` rows of `sample_metadata.csv`. Peak groups WT (4 mice), RAG2 (4 mice), Baseline (1 control). All 9 libraries form one TMM group, `WT_RAG2_Baseline`, as in the legacy `TARGETS = "WT|RAG2|Baseline"`.

## Run

1. New folder for the project: copy this example's `config.yaml` and `samplesheet.csv` into it, then edit them (the placeholder paths, the FASTQ paths).
2. In that folder, start the run with the launcher `run_atacseq.bat` (or `run_atacseq.ps1`). It uses the `config.yaml` in the folder.
3. Outputs go to `results/`, intermediates to `intermediate_files/` (only with `publish_intermediates: true`), run reports to `metadata/pipeline_info/`.

Startup validation (samplesheet columns, duplicate ids, missing FASTQs, group references) runs before the pipeline starts.

## Samplesheet: sample_id mapping to the legacy metadata

The legacy FASTQs were renamed by `Scripts/step_0b_rename_fastq_files.py` to `{genotype}_{mouse_id}_R1/_R2.fastq.gz`. The sample_id here is the name without `_R1`. The legacy `paired_samples.txt` name (`WT_110_R1`) is the same library.

| sample_id | legacy metadata sample_id | genotype (`peak_group`) | mouse_id (`replicate`) | FASTQ (in `Original_Data/fastq_files/`) |
|---|---|---|---|---|
| WT_110 | AM35 | WT | 110 | WT_110_R1/R2.fastq.gz |
| WT_111 | AM36 | WT | 111 | WT_111_R1/R2.fastq.gz |
| WT_113 | AM37 | WT | 113 | WT_113_R1/R2.fastq.gz |
| WT_114 | AM38 | WT | 114 | WT_114_R1/R2.fastq.gz |
| RAG2_115 | AM40 | RAG2 | 115 | RAG2_115_R1/R2.fastq.gz |
| RAG2_117 | AM41 | RAG2 | 117 | RAG2_117_R1/R2.fastq.gz |
| RAG2_118 | AM42 | RAG2 | 118 | RAG2_118_R1/R2.fastq.gz |
| RAG2_119 | AM45 | RAG2 | 119 | RAG2_119_R1/R2.fastq.gz |
| Baseline_CTRL | AM46 | Baseline | CTRL | Baseline_CTRL_R1/R2.fastq.gz |

Columns: `condition` = `genotype` and `replicate` = `mouse_id` (the metadata has no condition or replicate column). `tmm_group` = `tmm_group_legacy` = `WT_RAG2_Baseline`. `peak_group` = the metadata column `genotype`. `exclude` blank. All 18 FASTQs were checked with `ls` on 2026-10-26.

## Toggles

| Key | Value here (AV default) | AM legacy value | What it does |
|---|---|---|---|
| `macs3_fe_cutoff` | 2.0 | 1.0 | MACS3 `--fe-cutoff`; a value <= 1 is off and drops `.fe` from the prefix (legacy prefix `WT.qval0.05`) |
| `macs3_call_summits` | false | true | writes `{prefix}_summits.bed` (legacy wrote summits) |
| `macs3_min_length` | 250 | null | MACS3 `--min-length` (legacy: not set) |
| `macs3_max_gap` | 500 | null | MACS3 `--max-gap` (legacy: not set) |
| `mapq` | 30 | 0 (effective) | legacy wrote the literal `-q mapq` into samtools, which reads it as 0, so no MAPQ filter was applied. `mapq: 0` reproduces it |
| `remove_chroms` | chrM | '' (none) | chrM removal (new). `''` reproduces legacy |
| `tmm_group_column` | tmm_group | tmm_group_legacy | exact groups (here identical, one group) |
| `tn5_shift` | false | false | not legacy. `true` = deepTools `alignmentSieve --ATACshift` |
| `contamination_index` | null | not set | the contaminant index contains mm10 (239 mm10 contigs in its source FASTA), so mouse reads would align to it |

## Differences: AM legacy vs the defaults used here

The shared scripts were compared with `diff -u` (AM `Scripts/` vs AV_ATACSEQ_2026 `Scripts/`, 2026-10-26).

1. MACS3 fe cutoff: 1.0 (AM) vs 2.0 (AV). Here 2.0 (fewer, stronger peaks than legacy).
2. MACS3 summits: on (AM) vs off (AV). Here off.
3. MACS3 min length and max gap: not set (AM) vs 250 / 500 (AV). Here 250 / 500.
4. Contaminant index: set in AV config.py only. Not set here (mm10 issue).
5. bigWig: name `{sample}.tmm.bin_size_10.bw` (AM) vs `{sample}.tmm.bw` (AV); AM uses `bamCoverage -p max`. Bin size 10, `--skipNonCoveredRegions`, `--centerReads`, no `--extendReads` and scale 1/SF are the same. Naming and threads only.
6. TMM (`step_4b`, `scaleFactors.R`): AM hard-codes `TARGETS`, AV reads it from config (same groups here). `scaleFactors.R` differs only in package-loading and message lines (one hunk). Formula SF = normfacs × libSizes / 1e6, bin 10000 and EFFICIENCY_BIAS FALSE are identical.
7. `getPEsizes.R`: only package-loading and message lines differ. The max-frag bug at line 78 (keeps the first sample's maximum) is the same in both. The pipeline uses the true maximum over the group's BAMs (`tmm_max_frag: auto`).
8. MACS3 group match (`step_4c`): substring `group in sample_id` (AM) vs `re.search` (AV). Same result for WT, RAG2 and Baseline; here `peak_group` is an exact column.
9. Identical in AM and AV, so no toggle: trimming (`step_1`, NexteraPE-PE, LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20, MINLEN 25), BAM processing (`step_3`: the `-q mapq` literal, duplicate removal `sambamba markdup -r`, blacklist), `step_3d` (TSS), `step_4a`, `step_4qc`. bowtie2 (`step_2`) is also identical in AM and AV (legacy `--local --very-sensitive-local -N 1 -X 2000`); the pipeline now uses the CUT&RUN parameters plus `-X 2000`, with `-N 1` dropped (toggles `bowtie2_args` and `bowtie2_minins`).
10. Dataset or location only: genome, FASTQ naming (`FORWARD/REVERSE_PATTERN`), directory layout, `SCALINGDIR` (Importable_Data vs Analysis_Data).
11. `merge_peaks.sh`: AM = CBX2 nucleation/spreading subtraction (`bedtools subtract`, CBX2_Fig1 paths), not a merge of this dataset's peaks. AV = `bedops --merge` of high/low peaks. So there is no peak merge here.
12. Blacklist: AM removes by `bedtools intersect -v` on a BED round trip (`bedtobam`, which sets MAPQ 255 and drops sequences). The pipeline removes pairs with `pairtobed -type neither` and keeps MAPQ and sequences. This is a fix and has no toggle.

## Expected differences from the legacy outputs

- **MAPQ**: the legacy effective filter was `-q 0` (none). Here `mapq: 30` removes low-MAPQ reads, so fewer reads reach TMM, bigWigs and MACS3. The scale factors and peak counts will not match the legacy numbers.
- **chrM**: chrM reads are removed (`remove_chroms: chrM`), so library sizes drop by the chrM reads. A `qc/chrom_fraction` table is written.
- **MACS3**: fe 2.0 instead of 1.0, min-length 250 and max-gap 500 added, summits off. Peak counts will not match the legacy counts above; the prefix is `<group>.qval0.05.minlen250.maxgap500.fe2.0`.
- **Trimmomatic**: 0.40 biocontainer here; the legacy run used the 0.32 jar (`/opt/software/Trimmomatic/0.32`). Same steps and adapter set.
- **Blacklist**: pair-aware and sequence-preserving (legacy bed round trip).
- **TMM max frag**: the true maximum over each group's BAMs, not the first sample's. It may differ from 2000.
- **Names**: samples use `WT_110` (legacy `WT_110_R1`). bigWigs are named `<sample>.<method>.<units>.bw` (for example `WT_110.tmm.none.bw`), not the legacy `.tmm.bin_size_10.bw`.
- **Unchanged**: duplicate removal, bowtie2 settings, trimming, TMM groups and bin size, bigWig bin size and units (1/SF), no Tn5 shift.
