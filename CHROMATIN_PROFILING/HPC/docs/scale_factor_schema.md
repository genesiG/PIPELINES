# Scale-factor table schema (norm_spikein.py, norm_chipseqspikeinfree.R, norm_tmm.R, norm_cpm.py, norm_combine.py)

One tab-separated table per normalisation run and method, `<name>.<method>.scale_factors.tsv` (`scale_factors/methods/`), in the
7-column schema below. `<name>.scale_factors.tsv` is the PRIMARY method's table (a byte copy made by `norm_combine.py`; path and
content as before). Every downstream step (bigWigs, classification) reads ONLY the primary table and only these columns, so the
factor-direction bug class (multiply vs divide) cannot recur. The other methods are diagnostics (see "All methods side by side").

| column | meaning |
|---|---|
| `sample_id` | samplesheet id (join key) |
| `method` | `spikein`, `chipseqspikeinfree`, `tmm` or `cpm` (exactly one per table) |
| `raw_factor` | the method's native number: spikein `(S+T)/S`; chipseqspikeinfree the package `SF` (2 dp; NA if the package QC failed); tmm the TMM norm.factor `nf`; cpm 1 |
| `bigwig_multiplier` | value for `bamCoverage --scaleFactor`, to be used **with `--normalizeUsing CPM` or `RPKM`** (deepTools multiplies it into the CPM/RPKM value). Not valid with `None`, `BPM` or `RPGC` |
| `edger_norm_factor` | value for edgeR `norm.factors`; effective library size = `lib.size x edger_norm_factor`; geometric mean over the run = 1 |
| `qc_pass` | `TRUE`/`FALSE` (method-specific QC; see each script header) |
| `note` | free text: inputs, imputation, assumptions (no tabs). A diagnostic method that could not be computed (`--fail-soft`: package error, all libraries failing QC, 0 spike-in reads) has NA factors, `qc_pass` FALSE and a note starting with the token `UNAVAILABLE: <reason>`; a primary method never writes such a table (exit 3) |

## Direction (the one rule)

All methods express the factor as a **divisor of the depth-normalised signal**: `signal = CPM / nf`.

* `edger_norm_factor` = `nf` (centred so that its geometric mean is 1, the edgeR convention).
* `bigwig_multiplier` = `1 / nf` on the method's own scale (not centred, so absolute bigWig levels stay comparable to round 9).
* Invariant: `bigwig_multiplier x edger_norm_factor` is the same constant for every
  sample of a run (`1/geometric_mean(nf)`; exactly 1 for TMM).
* MA reference line of a pair (ref -> other), M = log2 CPM(other) - log2 CPM(ref) (the sign of the classification logFC; negative =
  less signal in `other`): `line = log2(edger_norm_factor_other / edger_norm_factor_ref)` = `log2(bigwig_multiplier_ref / bigwig_multiplier_other)`.
  It is the raw log-CPM difference at which that method declares "no change"; a fixed function of the factors (centring cancels, a
  constant multiple of one method's factors changes nothing), not estimated from the cloud. CPM: exactly 0. Spike-in:
  `log2(raw_ref / raw_other)`, with `line_se = sqrt(2/S_ref + 2/S_other) / ln 2` (a Poisson lower bound on the SE: counting noise of the spike-in reads S only). Only replicate pairs
  have a known expected offset (no biological change); TMM and CPM align the bulk of the bins by construction, so a line on the cloud
  is no evidence that they are right when the treatment changes the mark globally.
* `bigwig_offset` (pairs table of NORM_QC): the bigWigs divide by the final-BAM mapped reads B, edgeR by the csaw total T
  (fragments up to `max.frag`), so the bigWigs of a method are equal at `M = line + delta`, `delta = log2(B_other/T_other) - log2(B_ref/T_ref)`.
  NORM_QC warns when `|delta| > 0.01` for any pair.
* IgG/control libraries are normalised in their own run; never apply a target run's factors to IgG.

| method | `raw_factor` | `bigwig_multiplier` | `edger_norm_factor` before centring |
|---|---|---|---|
| spikein | `(S+T)/S` | `raw / min(raw)` (legacy `bamCov_SF`) | `min(raw) / raw` (legacy `SF`) |
| chipseqspikeinfree | `SF` | `1 / SF_used` | `SF_used` |
| tmm | `nf` (legacy table held `nf x lib.size / 1e6`; see `note`) | `1 / nf` | `nf` |
| cpm | `1` | `1` | `1` (library size only: bigWig = bamCoverage `--scaleFactor 1`; edgeR `lib.size` = csaw total, `norm.factors` = 1) |

`SF_used` = package SF; for libraries that fail the package QC (SF = NA) the run's maximum finite SF is imputed (legacy round-9 rule,
`--impute-failed=TRUE`), `qc_pass = FALSE` and the note says so. Those libraries carry no global-scaling information.

## All methods side by side (`<run>.scale_factors_all.tsv`)

Written by `norm_combine.py` (NORM_COMBINE), one row per library in run order: `sample_id`, `primary_method`, then for each method in the
fixed order chipseqspikeinfree, tmm, cpm, spikein the six columns `<m>_raw_factor`, `<m>_bigwig_multiplier`, `<m>_edger_norm_factor`,
`<m>_qc_pass`, `<m>_status`, `<m>_note`. All four blocks are always present. `<m>_status`: `ok`; `qc_fail` (finite factor, qc_pass FALSE);
`imputed` (chipseqspikeinfree only: qc_pass FALSE, raw_factor NA, factor finite); `unavailable` (no factor: note starts `UNAVAILABLE:`, or the
method was not launched for a stated reason, e.g. spike-in for a run with pre-aligned BAMs); `not_run` (not requested, e.g. spike-in with
`norm_spikein: false`; factor cells NA). Library sizes are not in this table (they are in `norm_qc/<run>.norm_qc_libsizes.tsv`).

## Worked examples (real round-9 NUC numbers, `Analysis_Data/normalization/`)

**ChIPseqSpikeInFree, K27M K27me3 (`max_turns` 0.86).** Package SF: EZH2i A 2.19, EZH2i B 1.85, DMSO A 1.00, DMSO B 1.08
(EZH2i has the larger divisor, i.e. lower signal, as expected after EZH2 inhibition). Geometric mean of SF = 1.4463.

| sample | raw_factor | bigwig_multiplier | edger_norm_factor | product |
|---|---|---|---|---|
| K27M_500nM_K27me3_A_R1 | 2.19 | 0.4566 | 1.5142 | 0.6914 |
| K27M_500nM_K27me3_B_R1 | 1.85 | 0.5405 | 1.2791 | 0.6914 |
| K27M_DMSO_K27me3_A_R1 | 1.00 | 1.0000 | 0.6914 | 0.6914 |
| K27M_DMSO_K27me3_B_R1 | 1.08 | 0.9259 | 0.7467 | 0.6914 |

**Spike-in, K27M K27me3** (`K27M_K27me3_spikein_SF.txt`): raw 391.789711 (EZH2i A), 485.080232 (EZH2i B), 1120.007249 (DMSO A),
1159.517940 (DMSO B); min = 391.789711. `bigwig_multiplier` = 1.000000, 1.238114, 2.858695, 2.959542 (identical to the legacy `bamCov_SF`
column); `min/raw` = 1.0000, 0.8077, 0.3498, 0.3379 (legacy `SF`), geometric mean 0.55586, so `edger_norm_factor` = 1.7990, 1.4530, 0.6293,
0.6079 (product 1.7990 for every sample). DMSO has the low spike fraction (high raw), hence the larger multiplier: more H3K27me3 per cell.

**QC failure with imputation, K27M-KO K27me3 at `max_turns` 0.72.** Package output: EZH2i A 2.29, EZH2i B 2.26, DMSO A NA, DMSO B NA
(`fail: complete loss, input or poor enrichment`). Imputed value 2.29 -> `raw_factor` 2.29, 2.26, NA, NA; `bigwig_multiplier` 0.4367,
0.4425, 0.4367, 0.4367; `edger_norm_factor` 1.0033, 0.9902, 1.0033, 1.0033; `qc_pass` TRUE, TRUE, FALSE, FALSE. The DMSO bigWigs then
carry no DMSO-vs-EZH2i scaling. There is no round-9 TMM table
for NUC; the TMM algebra is therefore only demonstrated on synthetic data.

## How to consume

```
bamCoverage ... --normalizeUsing RPKM --scaleFactor <bigwig_multiplier>
DGEList(counts, lib.size = <csaw totals>, norm.factors = <edger_norm_factor>)
```
A missing/NA factor must stop the consumer (exit 3), never default to 1.

## Known approximations (state them if you report these factors)

* ChIPseqSpikeInFree SF is rounded to 2 decimals by the package; its read counting treats every mate as a read (singleEnd default) and
  applies no MAPQ/duplicate filters (the BAMs must already be final).
* spike-in `raw = (S+T)/S` differs from `T/S` by the factor `1 + S/T` (< 0.3 percent for the NUC libraries); kept from legacy step_4c.
* `bigwig_multiplier` assumes deepTools' CPM/RPKM total is proportional to the csaw library total of the same BAM (true up to the small
  fraction of fragments removed by csaw's `max.frag`); the difference per pair is `bigwig_offset` in `norm_qc_pairs.tsv`. `max.frag` is
  `csaw_max_frag`, else `bowtie2_maxins`, else 500 in BIN_COUNTS (NORM_QC, TMM, replicate metrics), CSAW and CLASSIFY_REGIONS.
* The MA lines of CSIF carry the 2-decimal rounding of SF (up to +-log2(1 + 0.005/SF) per library, about +-0.007 at SF = 1).
