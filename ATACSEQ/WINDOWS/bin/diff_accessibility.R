#!/usr/bin/env Rscript
# diff_accessibility.R -- differential chromatin accessibility (csaw + edgeR) for ONE comparison, over peaks and/or promoters.
# Port of the legacy csaw driver: Scripts/differential_accessibility.Rmd chunks "loop" (per cell line) and "loop 2" (gd2
# additive model), Scripts/utils.R (import_and_merge_peaks, the two diagnostic plots) and Scripts/csaw_parameters.R, plus the
# promoter chunks of the OLDER copy (Processing_and_Analysis/Analysis_Data/csaw/differential_accessibility.Rmd, lines ~436-720:
# "load annotations", "get promoter regions", "test differential accessibility over promoters").
# The computation is the legacy code (same csaw/edgeR calls, arguments and order); only the I/O is replaced by the CLI.
# Peaks and promoters of one comparison run in one task so that the 10 kb bins are counted once.
#
# CLI contract (called by the DIFF_ACCESS module; --key=value):
#   --name=NAME                  REQUIRED. Comparison name = output prefix. ^[A-Za-z0-9][A-Za-z0-9_.-]*$
#   --bams=a.bam,b.bam           REQUIRED. Coordinate-sorted, indexed, filtered PAIRED-END BAMs (names relative to cwd are fine).
#   --sample-ids=x,y             REQUIRED. Same order and length as --bams. Used in the norm-factor and replicate tables.
#   --groups=Control,Condition   REQUIRED. Same order. Control or Condition per BAM. Each side needs >= 1 BAM.
#   --blocks=b1,b2               OPTIONAL. Same order. One label per BAM (the legacy gd2 cell line). With >= 2 distinct labels the
#                                model is ~ block + grouping. Omitted, blank, or one label: no blocks.
#   --control-peaks=f1,f2        REQUIRED. Control peak files (.narrowPeak, .broadPeak or .bed; any number). Merged by union.
#   --case-peaks=f1,f2           REQUIRED. Case peak files; merged by union. Tested peak regions = union(control, case).
#   --regions=peaks[,promoters]  OPTIONAL, default peaks. Region sets tested in this task (one or both).
#   --txdb=F                     OPTIONAL. TxDb SQLite (txdb_build.R). REQUIRED when promoters are requested (exit 3 otherwise).
#   --genome-build=BUILD         REQUIRED. Written to the GRanges with genome<- (e.g. hs1).
#   --style=UCSC                 OPTIONAL, default UCSC. Seqname style of the PEAK outputs (mapSeqlevels). Promoter outputs are
#                                not restyled (legacy).
#   --minq=10                    OPTIONAL, default 10. csaw readParam(minq=): minimum mapping quality.
#   --pe=both|none               OPTIONAL, default both. csaw readParam(pe=).
#   --max-frag=500               OPTIONAL, default 500. csaw readParam(max.frag=) (csaw default; the legacy never set it).
#   --dedup=true|false|auto      OPTIONAL, default auto. auto = TRUE when the comparison has <= 2 BAMs, else FALSE (legacy rule).
#   --bin-size=10000             OPTIONAL, default 10000. Width of the bins that give the normalisation factors.
#   --lfc=0.585                  OPTIONAL, default 0.585. |logFC| threshold for the up/down BEDs and the printed counts.
#   --fdr=0.1                    OPTIONAL, default 0.1. FDR threshold for the summary counts only (no table is filtered by it).
#   --table=best|combined        OPTIONAL, default best. Which overlapResults table gives the statistics: merged$best (legacy) or
#                                merged$combined.
#   --qlf-min-libs=3             OPTIONAL, default 3. Quasi-likelihood F-test when the number of libraries >= this, else the LRT.
#   --fixed-dispersion=0.05      OPTIONAL, default 0.05. Dispersion used by glmFit in the LRT branch (legacy per-line 0.05;
#                                the legacy gd2 chunk used 0.15 and is therefore set per comparison, not hard-coded).
#   --promoter-up=2000           OPTIONAL, default 2000. promoters(upstream=) in bp.
#   --promoter-down=200          OPTIONAL, default 200. promoters(downstream=) in bp.
#   --promoter-min-fc=4          OPTIONAL, default 4. Promoters are kept when background abundance > log2(min-fc) (filterWindowsGlobal).
#   --cpus=1                     OPTIONAL, default 1. >1: BiocParallel MulticoreParam(cpus) for the read counting; 1: SerialParam.
#   --case-label=Condition       OPTIONAL. Free text written to the replicates table (Case column).
#   --control-label=Control      OPTIONAL. Free text written to the replicates table (Control column).
#   --font=Helvetica             OPTIONAL, default Helvetica. Font family of every plot (svg: font-family of the text elements).
#   --plot-format=svg|pdf        OPTIONAL, default svg. svg: one file per diagnostic (svglite, editable <text>). pdf: all diagnostics
#                                as pages of one legacy NAME.diagnostics.pdf.
#   --plot-width=8               OPTIONAL, default 8 (inches). Width of each diagnostic (MA grid: width of the whole grid).
#   --plot-height=4              OPTIONAL, default 4 (inches). Height of one panel row (MA grid: height x number of grid rows).
#
# Outputs (cwd, prefix = --name; NAME below). Peak outputs only when regions include peaks; promoter outputs only when they include promoters.
#   NAME.csv                     peak table, write.csv (row names 1..n): "", seqnames, start, end, width, strand, PValue, FDR, logFC, score, name
#   NAME.bed                     all peak regions (rtracklayer::export.bed; makeGRangesFromDataFrame(keep.extra.columns = TRUE))
#   NAME_upregulated.bed         peaks with logFC > lfc (NO FDR filter: these are fold-change sets, not significant sets)
#   NAME_downregulated.bed       peaks with logFC < -lfc (same)
#   NAME_promoters.csv           promoter table, write.csv (row names = gene_id): "", seqnames, start, end, width, strand, gene_id,
#                                PValue, FDR, logFC, score, name (name = gene_id); promoters with NA score removed
#   NAME_promoters.bed           all promoters with a score
#   NAME_promoters_up.bed        promoters with logFC > lfc (no FDR filter)
#   NAME_promoters_down.bed      promoters with logFC < -lfc (no FDR filter)
#   NAME.norm_factors.tsv        sample_id, norm_factor (TMM, weighted = FALSE), lib_size (csaw totals of the bins), group, block
#   NAME.summary.tsv             region_set, metric, value: counts (increased / decreased by lfc; FDR <= fdr; significant x direction)
#                                and logFC quantiles over all regions and over the FDR-significant ones
#   NAME.replicates.csv          Comparison, Case, Control, Case_Replicates, Control_Replicates, Total_Samples
#   NAME.diag_ma.svg             normalisation MA plots, library 1 vs each other library, all in one grid (2 columns); red line =
#                                log2(nf1/nfj); first region set only
#   NAME.diag_bcv.svg            trended BCV plot (>= 3 libraries); promoters: NAME.diag_promoters_bcv.svg
#   NAME.diag_mds.svg            MDS plot (>= 3 libraries); promoters: NAME.diag_promoters_mds.svg
#   NAME.diag_qldisp.svg         QL dispersion plot (QL branch only); promoters: NAME.diag_promoters_qldisp.svg
#   NAME.diag_promoter_filter.svg  promoters only: histogram of the background abundance (log2 CPM) with the log2(promoter-min-fc) cut
#   NAME.diagnostics.pdf         --plot-format=pdf only: the same plots as pages of this one file (legacy file name)
#
# Exit codes: 0 ok; 2 usage error; 3 known data problem: < 1 BAM on one side, a BAM or its .bai missing, promoters requested
# without a usable --txdb, no peak file or an unknown peak suffix, a peak file missing; --blocks confounded with --groups (no
# block holds both a Control and a Condition library: design not full rank); no residual degrees of freedom for the quasi-likelihood
# branch; an empty region set; regions whose seqnames are not contigs of the BAM header (e.g. peaks "1" vs BAM "chr1"); BAMs with
# different headers. All of these are detected before any read is counted.
#
# ---- Question, model, assumptions -----------------------------------------------------------------------------------------
# Question: which accessibility regions (peaks of the union of the case and control peak sets, or promoters of the TxDb genes)
# have a different fragment count between the Case (Condition) and Control libraries? Output per region: PValue, FDR, logFC
# (edgeR), score = logFC (legacy).
# Counts: paired-end fragments with MAPQ >= minq and fragment length <= max.frag, duplicates dropped when dedup is TRUE,
# summed over each region (csaw regionCounts). The 10 kb bins (windowCounts, bin = TRUE) are counted once, for the whole genome,
# and give the library sizes and the TMM normalisation factors (normFactors, weighted = FALSE). The same factors apply to every
# region set; promoters are first filtered by their background abundance (filterWindowsGlobal > log2(promoter-min-fc)).
# Model A (no --blocks): ~ 0 + grouping (Control, Condition); test = Condition - Control. (legacy per-cell-line model)
# Model B (--blocks with >= 2 labels): ~ block + grouping; test = coefficient groupingCondition, i.e. the Condition effect
# additive to the block effect (legacy gd2 model, block = cell line, first level 6349 as baseline).
# Test: quasi-likelihood F-test (glmQLFit robust = TRUE, glmQLFTest) when the number of libraries >= qlf-min-libs; otherwise
# glmFit with dispersion fixed at fixed-dispersion and a likelihood-ratio test (glmLRT). The two branches are not
# interchangeable: the LRT with a fixed dispersion ignores the uncertainty of the dispersion and is anti-conservative for small
# n; it is used only because the legacy did so for 2-BAM comparisons.
# Overlap: overlapResults() combines the bin-level tests into one row per region; the row statistics come from merged$best
# (legacy, --table=best) or merged$combined (--table=combined).
# Assumptions: (1) most of the genome is not differentially accessible between the libraries (the TMM assumption, as in
# norm_tmm.R); (2) the BAMs are the final duplicate-removed BAMs, so dedup (drop the duplicate flag 0x400) normally changes
# nothing; (3) biological replicates are the libraries, and each library is one sample. FDR is the edgeR BH adjustment over
# the tested regions of this comparison only.
#
# ---- Statistical limits of the legacy method, ported unchanged (state them when reporting results) ---------------------------
# * Fixed dispersion (fewer than --qlf-min-libs libraries, e.g. one library per side, as for K27MKO, HGGx42c and TS603 in the
#   AV project): there is no replication, so the variance is ASSUMED (0.05), not estimated. PValue and FDR are then not valid
#   error rates: they are anti-conservative whenever the true between-library dispersion is larger than the fixed value (the
#   legacy gd2 chunk itself used 0.15 for similar data), and they cannot show that a difference is reproducible. Use such tables
#   to rank regions by logFC only; do not call regions "significantly differentially accessible" from them.
# * Pooled / blocked comparisons (legacy gd2, --blocks = cell line): the unit of replication for a claim about "GD2-high versus
#   GD2-low" is the cell line, not the library. The additive block removes each line's baseline but assumes ONE common
#   high-vs-low effect in all lines (no line x condition interaction). With one library per line and condition, the residual
#   degrees of freedom are the interaction itself (a paired analysis over the lines: valid for the lines studied). Extra
#   replicate libraries of one line (e.g. 6349 high A and B) add residual degrees of freedom that measure within-line noise only,
#   so the p-values can be optimistic for a statement about lines in general. A block that holds libraries of one condition only
#   (BT88_high_A in the legacy 12-BAM set) carries no information on the effect (a warning is printed). The effective n is the
#   number of lines that hold both conditions (printed), not the number of BAMs.
# * Tested regions come from the data: MACS3 peaks are called per group on POOLED replicates and the union of both groups is
#   tested. The selection uses the same reads as the test (symmetric between the groups, so the contrast is not biased much, but
#   a group with more replicates or depth has more power to call a peak, and regions open only in the better-powered group are
#   over-represented). The union also fuses overlapping or adjacent peaks into wider regions. FDR is over this region set only.
# * Normalisation: TMM on 10 kb background bins (weighted = FALSE) assumes that most of the genome does not change and corrects
#   composition (background) bias. It does NOT correct efficiency bias: libraries with a different signal-to-noise ratio
#   (FRiP, TSS enrichment) shift the logFC of all peaks together. If the compared states differ globally in accessibility, the
#   assumption fails and the logFC of every region is biased; check NAME.diag_ma (the red line must run through the bulk of
#   the bins) and the TMM factors in NAME.norm_factors.tsv, and compare with a peak-based normalisation before claiming global
#   gains or losses.
# * Promoters: the 2200 bp windows of neighbouring genes overlap (about 23% of the 28710 hs1 RefSeq promoters overlap another
#   promoter window, 155 share identical coordinates), so promoter rows are not independent. overlapResults() gives each promoter
#   the best window among all windows that overlap it (Holm-adjusted for their number), so a promoter row can carry the statistics
#   and logFC of its neighbour's window, and an overlapped cluster appears in several rows. The abundance filter
#   (filterWindowsGlobal, background log CPM) does not use the contrast and is a valid independent filter.
# * lfc is a threshold on the estimated logFC (not a TREAT test) and is applied independently of FDR: the up/down BEDs contain
#   regions that are not significant.
# * Group assignment: the legacy set Control / Condition and the cell-line block by grepl() on the BAM file names (the whole
#   path), so a name that matched nothing became NA or "Other" (legacy: x42c and K27M libraries fell into one "Other" block) and
#   the promoter chunk called every library without "low" a Condition. Here the groups and blocks are explicit CLI arguments.
#
# ---- Where the port deviates from the legacy, and why ------------------------------------------------------------------------
# * Dedup: legacy set param@dedup <- FALSE in the loop when a comparison had > 2 BAMs, and that global stayed FALSE for every
#   later comparison of the same R session. Here the rule is per comparison (auto = TRUE iff <= 2 BAMs).
# * max.frag: the legacy never set it; csaw's default 500 is now passed explicitly (same value).
# * Fixed dispersion: the legacy per-line branch used 0.05 and the gd2 branch 0.15. The gd2 branch is unreachable with 12
#   libraries; the value is an argument.
# * Stale my.regions (not ported): the legacy gd2 PROMOTER chunk counted reads over the my.regions left over from the previous
#   loop (the peak regions of the last comparison), not over the promoters. Promoters use the same code for every design here.
# * Legacy gd2 chunk: the NA-score filter and the seqnames filter applied before the table was overwritten were no-ops; dropped.
# * The legacy gd2 BEDs were named _up.bed / _down.bed; the names are now _upregulated.bed / _downregulated.bed (the per-line
#   names).
# * Not ported (not an output of this pipeline): the design_matrix file, the replicate tracker across comparisons
#   (cell_line_replicates_summary.csv), the per-cell-line <line>_promoters_high_vs_low.csv table, the cross-cell-line consensus
#   tables and the top-10 printouts. The printed legacy counts go to stdout and NAME.summary.tsv. The legacy filter histogram
#   (hist of the promoter background abundance with the cut) is drawn as NAME.diag_promoter_filter.svg.
#   Plots: svg (default) replaces the legacy single pdf page set; one svg per diagnostic (svglite writes real <text> elements
#   with the requested font family; close_svg rewrites the family name that fontconfig resolves, e.g. Nimbus Sans, to --font).
#
options(stringsAsFactors = FALSE)

die <- function(msg, status = 3L) {
  cat("ERROR: ", msg, "\n", sep = "", file = stderr())
  quit(save = "no", status = status)
}
parse_args <- function(argv, allowed) {
  out <- list()
  for (a in argv) {
    m <- regmatches(a, regexec("^--([^=]+)=(.*)$", a))[[1]]
    if (length(m) != 3L) die(paste0("cannot parse argument '", a, "' (expected --key=value)"), 2L)
    if (!(m[2] %in% allowed)) die(paste0("unknown argument --", m[2]), 2L)
    out[[m[2]]] <- m[3]
  }
  out
}
get_req <- function(args, key) {
  v <- args[[key]]
  if (is.null(v) || !nzchar(v)) die(paste0("missing required argument --", key), 2L)
  v
}
get_opt <- function(args, key, default) {
  v <- args[[key]]
  if (is.null(v) || !nzchar(v)) default else v
}
split_list <- function(x) { v <- trimws(strsplit(x, ",", fixed = TRUE)[[1]]); v[nzchar(v)] }
as_num <- function(x, key) {
  v <- suppressWarnings(as.numeric(x))
  if (length(v) != 1L || is.na(v)) die(paste0("--", key, " must be a number, got '", x, "'"), 2L)
  v
}

# ---- plots: one svg per diagnostic (svglite, editable <text>, font family = --font); --plot-format=pdf = one legacy pdf ----
# open_svg / close_svg are copied from the RNASEQ module (bin/deg_common.R). close_svg rewrites the family names that
# svglite resolves through fontconfig (e.g. "Nimbus Sans" for Helvetica) to the requested --font; widths are kept.
open_svg <- function(path, width, height) {
  svglite::svglite(filename = path, width = width, height = height,
    system_fonts = list(sans = plot_font, serif = plot_font, mono = plot_font, symbol = plot_font))
  invisible(path)
}
close_svg <- function(path) {
  grDevices::dev.off()
  txt <- readLines(path, warn = FALSE, encoding = "UTF-8")
  txt <- gsub("font-family: ('[^']*'|\"[^\"]*\"|[^;]*);", paste0("font-family: ", plot_font, ";"), txt)
  writeLines(txt, path, useBytes = TRUE)
  invisible(NULL)
}
# One diagnostic. Returns its svg path, or NULL in pdf mode (the plot then goes to the shared pdf). Pair with diag_close().
diag_device <- function(stub, width, height) {
  if (plot_format == "pdf") return(NULL)
  path <- paste0(name, ".", stub, ".svg")
  open_svg(path, width, height)
  path
}
diag_close <- function(path) {
  if (!is.null(path)) close_svg(path)
  invisible(NULL)
}

detect_peak_format <- function(fname) {
  if (grepl("\\.narrowPeak$", fname, ignore.case = TRUE)) return("narrowPeak")
  if (grepl("\\.broadPeak$",  fname, ignore.case = TRUE)) return("broadPeak")
  if (grepl("\\.bed$",        fname, ignore.case = TRUE)) return("BED")
  die(paste0("unknown peak format (expected .narrowPeak, .broadPeak or .bed): ", fname), 3L)
}

# Legacy utils.R import_and_merge_peaks, with full file paths instead of (peaks, peakdir)
import_and_merge_peaks <- function(peaks) {

  gr_list <- vector("list", length(peaks))

  for (i in seq_along(peaks)) {

    format <- detect_peak_format(peaks[i])

    gr <- rtracklayer::import(
      con    = peaks[i],
      format = format
    )

    gr <- GenomeInfoDb::keepStandardChromosomes(
      gr, pruning.mode = "coarse"
    )

    gr_list[[i]] <- gr
  }

  # Merge all GRanges into a single union set
  Reduce(union, gr_list)
}

# Summary rows of one region set: legacy printed counts (table() calls) and logFC summaries, written as rows
summary_rows <- function(tabcom, set_label) {
  is.sig <- tabcom$FDR <= qvalue
  counts <- c(
    n_regions               = nrow(tabcom),
    n_increased             = sum(tabcom$rep.logFC > lfc, na.rm = TRUE),
    n_decreased             = sum(tabcom$rep.logFC < -lfc, na.rm = TRUE),
    n_fdr_significant       = sum(is.sig, na.rm = TRUE),
    n_significant_increased = sum((tabcom$rep.logFC > lfc)[is.sig], na.rm = TRUE),
    n_significant_decreased = sum((tabcom$rep.logFC < -lfc)[is.sig], na.rm = TRUE)
  )
  stat_names <- c("min", "q1", "median", "mean", "q3", "max")
  q_all <- suppressWarnings(as.numeric(summary(tabcom$rep.logFC)))[1:6]
  q_sig <- suppressWarnings(as.numeric(summary(tabcom$rep.logFC[is.sig])))[1:6]
  data.frame(
    region_set = set_label,
    metric = c(names(counts), paste0("logFC_", stat_names), paste0("significant_logFC_", stat_names)),
    value = c(counts, q_all, q_sig),
    stringsAsFactors = FALSE
  )
}

# Shared statistics: model fit and test for one count matrix (legacy "STATISTICAL MODELING" chunk of each loop)
test_model <- function(y, tag = "") {
  # Estimate dispersions
  y <- estimateDisp(y, design)

  cat("\nSummary of trended dispersions:\n")
  y$trended.dispersion %>% summary %>% print()

  if (ncol(y) > 2) {
    # No residual df (only reachable with a raised --qlf-min-libs): the dispersions are NA and plotBCV() fails
    if (ncol(y) - design_rank >= 1L) {
      bcv_path <- diag_device(paste0("diag_", tag, "bcv"), plot_width, plot_height)
      plotBCV(y)
      diag_close(bcv_path)
    }

    # MDS Plot
    cat("PCA plot:\n")
    mds_path <- diag_device(paste0("diag_", tag, "mds"), plot_width, plot_height)
    plotMDS(cpm(y, log = TRUE),
            top = 10000,
            gene.selection = "common", # creates PCA plot
            labels = ids,
            col = c("black", "purple")[as.integer(grouping)])
    diag_close(mds_path)
  }

  # Estimate quasi-likelihood dispersions
  if (ncol(y) >= qlf_min) {
    cat("Estimating quasi-likelihood dispersions...\n")
    fit <- glmQLFit(y, design, robust = TRUE)

    ql_path <- diag_device(paste0("diag_", tag, "qldisp"), plot_width, plot_height)
    plotQLDisp(fit)
    diag_close(ql_path)

    if (use_blocks) {
      res <- glmQLFTest(fit, coef = coef_name)
    } else {
      res <- glmQLFTest(fit, contrast = contrast)
    }
  } else {
    cat("Not enough libraries to estimate dispersion (", ncol(y), " < ", qlf_min, ")\n", sep = "")
    cat("Setting dispersion manually to ", fixed_disp, "\n", sep = "")
    fit <- glmFit(y, design, dispersion = fixed_disp)

    if (use_blocks) {
      res <- glmLRT(fit, coef = coef_name)
    } else {
      res <- glmLRT(fit, contrast = contrast)
    }
  }

  cat("Summary of results:\n")
  res$table$logFC %>% range %>% print()
  res$table$logFC %>% summary %>% print()

  res
}

allowed <- c("name", "bams", "sample-ids", "groups", "blocks", "control-peaks", "case-peaks", "regions", "txdb",
             "genome-build", "style", "minq", "pe", "max-frag", "dedup", "bin-size", "lfc", "fdr", "table",
             "qlf-min-libs", "fixed-dispersion", "promoter-up", "promoter-down", "promoter-min-fc", "cpus",
             "case-label", "control-label", "font", "plot-format", "plot-width", "plot-height")
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) == 0L || any(argv %in% c("-h", "--help"))) {
  cat("usage: diff_accessibility.R --name=NAME --bams=a.bam,b.bam --sample-ids=x,y --groups=Control,Condition",
      "[--blocks=..] --control-peaks=.. --case-peaks=.. --regions=peaks[,promoters] [--txdb=F] --genome-build=hs1",
      "[--style=UCSC --minq=10 --pe=both --max-frag=500 --dedup=auto --bin-size=10000 --lfc=0.585 --fdr=0.1",
      "--table=best --qlf-min-libs=3 --fixed-dispersion=0.05 --promoter-up=2000 --promoter-down=200",
      "--promoter-min-fc=4 --cpus=1 --case-label= --control-label=]\n", file = stderr())
  quit(save = "no", status = 2L)
}
args <- parse_args(argv, allowed)
name          <- get_req(args, "name")
bams          <- split_list(get_req(args, "bams"))
ids           <- split_list(get_req(args, "sample-ids"))
groups        <- split_list(get_req(args, "groups"))
blocks_arg    <- get_opt(args, "blocks", "")
blocks        <- if (nzchar(trimws(blocks_arg))) split_list(blocks_arg) else NULL
ctrl.peaks    <- split_list(get_req(args, "control-peaks"))
cond.peaks    <- split_list(get_req(args, "case-peaks"))
region_sets   <- unique(split_list(get_opt(args, "regions", "peaks")))
txdb_file     <- get_opt(args, "txdb", "")
build         <- get_req(args, "genome-build")
style         <- get_opt(args, "style", "UCSC")
minq          <- as_num(get_opt(args, "minq", "10"), "minq")
pe            <- get_opt(args, "pe", "both")
max_frag      <- as_num(get_opt(args, "max-frag", "500"), "max-frag")
dedup_arg     <- tolower(get_opt(args, "dedup", "auto"))
bin.size      <- as_num(get_opt(args, "bin-size", "10000"), "bin-size")
lfc           <- as_num(get_opt(args, "lfc", "0.585"), "lfc")
qvalue        <- as_num(get_opt(args, "fdr", "0.1"), "fdr")
table_sel     <- tolower(get_opt(args, "table", "best"))
qlf_min       <- as_num(get_opt(args, "qlf-min-libs", "3"), "qlf-min-libs")
fixed_disp    <- as_num(get_opt(args, "fixed-dispersion", "0.05"), "fixed-dispersion")
prom_up       <- as_num(get_opt(args, "promoter-up", "2000"), "promoter-up")
prom_down     <- as_num(get_opt(args, "promoter-down", "200"), "promoter-down")
min.fc        <- as_num(get_opt(args, "promoter-min-fc", "4"), "promoter-min-fc")
cpus          <- as_num(get_opt(args, "cpus", "1"), "cpus")
case_label    <- get_opt(args, "case-label", "Condition")
control_label <- get_opt(args, "control-label", "Control")
plot_font     <- get_opt(args, "font", "Helvetica")
plot_format   <- tolower(get_opt(args, "plot-format", "svg"))
plot_width    <- as_num(get_opt(args, "plot-width", "8"), "plot-width")
plot_height   <- as_num(get_opt(args, "plot-height", "4"), "plot-height")

# ---- argument checks (exit 2 usage, exit 3 data) --------------------------------------------------------------------------
if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", name)) die("--name must match ^[A-Za-z0-9][A-Za-z0-9_.-]*$", 2L)
if (length(ids) != length(bams) || length(groups) != length(bams)) {
  die(sprintf("--bams (%d), --sample-ids (%d) and --groups (%d) must have the same length",
              length(bams), length(ids), length(groups)), 2L)
}
if (!is.null(blocks) && length(blocks) != length(bams)) {
  die(sprintf("--blocks (%d) must have one label per BAM (%d)", length(blocks), length(bams)), 2L)
}
if (!all(groups %in% c("Control", "Condition"))) die("--groups must be Control or Condition for every BAM", 2L)
if (anyDuplicated(ids)) die(paste0("duplicated sample ids: ", paste(ids[duplicated(ids)], collapse = ", ")), 2L)
if (anyDuplicated(bams)) die("duplicated BAM paths in --bams", 2L)
if (length(region_sets) == 0L || !all(region_sets %in% c("peaks", "promoters"))) {
  die("--regions must be peaks, promoters or peaks,promoters", 2L)
}
if (!(pe %in% c("both", "none"))) die("--pe must be both or none", 2L)
if (!(dedup_arg %in% c("true", "false", "auto"))) die("--dedup must be true, false or auto", 2L)
if (!(table_sel %in% c("best", "combined"))) die("--table must be best or combined", 2L)
if (bin.size < 1000 || bin.size != round(bin.size)) die("--bin-size must be an integer >= 1000 bp", 2L)
if (cpus < 1 || cpus != round(cpus)) die("--cpus must be a positive integer", 2L)
if (qlf_min < 1) die("--qlf-min-libs must be >= 1", 2L)
if (!(plot_format %in% c("svg", "pdf"))) die("--plot-format must be svg or pdf", 2L)
if (plot_width <= 0 || plot_height <= 0) die("--plot-width and --plot-height must be > 0", 2L)

n_ctrl <- sum(groups == "Control")
n_case <- sum(groups == "Condition")
if (n_ctrl < 1L || n_case < 1L) {
  die(sprintf("each side needs at least one BAM (Control: %d, Condition: %d)", n_ctrl, n_case))
}
miss <- bams[!file.exists(bams)]
if (length(miss)) die(paste0("BAM file(s) not found: ", paste(miss, collapse = ", ")))
has_idx <- file.exists(paste0(bams, ".bai")) | file.exists(sub("\\.bam$", ".bai", bams))
if (!all(has_idx)) die(paste0("BAM index (.bai) missing for: ", paste(bams[!has_idx], collapse = ", "),
                              " (index the BAMs with samtools index)"))
if (length(ctrl.peaks) == 0L || length(cond.peaks) == 0L) {
  die("both --control-peaks and --case-peaks need at least one peak file")
}
peak_files <- c(ctrl.peaks, cond.peaks)
miss <- peak_files[!file.exists(peak_files)]
if (length(miss)) die(paste0("peak file(s) not found: ", paste(miss, collapse = ", ")))
for (p in peak_files) detect_peak_format(p)
if ("promoters" %in% region_sets &&
    (!nzchar(txdb_file) || !file.exists(txdb_file) || file.size(txdb_file) == 0)) {
  die("promoters requested but no usable TxDb: give --txdb (txdb_build.R output)", 3L)
}

# Dedup rule (legacy): TRUE only for comparisons with <= 2 BAMs
dedup <- if (dedup_arg == "auto") length(bams) <= 2L else dedup_arg == "true"

# Attach order as in the legacy (dplyr first, csaw last): dplyr masks union()/filter() otherwise and the GRanges union breaks
suppressPackageStartupMessages({
  library(dplyr)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(GenomicFeatures)
  library(rtracklayer)
  library(BiocParallel)
  library(edgeR)
  library(csaw)
})

BPPARAM <- if (cpus > 1) MulticoreParam(cpus) else SerialParam()

# ---- design (legacy "Generate design matrix") -----------------------------------------------------------------------------
grouping <- factor(groups, levels = c("Control", "Condition"))
use_blocks <- !is.null(blocks) && length(unique(blocks)) >= 2L
if (use_blocks) {
  block_factor <- factor(blocks)
  design <- model.matrix(~ block_factor + grouping)
  coef_name <- "groupingCondition"
} else {
  design <- model.matrix(~0 + grouping)
  colnames(design) <- levels(grouping)
  contrast <- makeContrasts(Condition - Control, levels = design)
}

# ---- design checks (exit 3, before any read is counted) --------------------------------------------------------------------
# The group indicator is in the span of the block indicators exactly when every block holds libraries of ONE group only; then
# groupingCondition is not estimable (design rank < ncol). One informative block (both groups present) is enough for full rank.
# A block with one group only is allowed (full rank) but carries no information on the Condition effect: warn.
n_lib <- length(bams)
design_rank <- qr(design)$rank
if (use_blocks) {
  blk_tab <- table(block_factor, grouping)
  blk_desc <- paste0(rownames(blk_tab), " (Control ", blk_tab[, "Control"], ", Condition ", blk_tab[, "Condition"], ")",
                     collapse = "; ")
  informative <- rownames(blk_tab)[blk_tab[, "Control"] > 0L & blk_tab[, "Condition"] > 0L]
  if (design_rank < ncol(design) || length(informative) == 0L) {
    die(paste0("--blocks are confounded with --groups: no block holds both Control and Condition libraries, so the Condition ",
               "effect cannot be separated from the block effect (design matrix rank ", design_rank, " < ", ncol(design),
               " columns). Libraries per block: ", blk_desc, ". Give every analysed block at least one Control and one ",
               "Condition library, or leave the block column empty."))
  }
  single_group <- setdiff(rownames(blk_tab), informative)
  if (length(single_group)) {
    cat("WARNING: block(s) with libraries of one group only (no information on the Condition effect): ",
        paste(single_group, collapse = ", "), ". Informative blocks: ", length(informative), " (", paste(informative, collapse = ", "),
        ")\n", sep = "", file = stderr())
  }
  if (length(informative) == 1L) {
    cat("WARNING: only one block (", informative, ") holds both groups: the Condition effect is estimated from that block alone.\n",
        sep = "", file = stderr())
  }
  cat("Blocks: ", blk_desc, "; informative blocks (effective n for the group effect): ", length(informative), "\n", sep = "")
}
if (n_lib >= qlf_min && n_lib - design_rank < 1L) {
  die(paste0("no residual degrees of freedom for the quasi-likelihood test: ", n_lib, " libraries and a design of rank ",
             design_rank, ". Add libraries (one more per block and group is enough), or drop --blocks. (--qlf-min-libs=", qlf_min,
             "; raising it would switch to the fixed-dispersion LRT, which is not a valid test: see the header.)"))
}

# ---- read counting setup (legacy csaw_parameters.R, dedup per comparison) --------------------------------------------------
param <- readParam(minq = minq,
                   pe = pe,
                   max.frag = max_frag,
                   dedup = dedup)
cat("Using the following parameters for read counting: \n")
print(param)
cat("Dedup for this comparison (", length(bams), " BAMs): ", dedup, "\n", sep = "")

if (plot_format == "pdf") {
  if (tolower(plot_font) %in% c("helvetica", "arial", "sans")) {
    grDevices::pdf(file = paste0(name, ".diagnostics.pdf"), width = plot_width, height = plot_height, family = "Helvetica")
  } else {
    grDevices::cairo_pdf(filename = paste0(name, ".diagnostics.pdf"), width = plot_width, height = plot_height, family = plot_font)
  }
}

# ---- regions ------------------------------------------------------------------------------------------------------------
# Contigs of the BAM headers (csaw counts regions on a contig that is not in the BAM as zero without a message, and an
# all-zero table then fails inside edgeR with an unreadable error).
bam_contigs <- lapply(bams, function(b) Rsamtools::scanBamHeader(b)[[1]]$targets)
if (any(vapply(bam_contigs, function(x) !identical(sort(names(x)), sort(names(bam_contigs[[1]]))), logical(1)))) {
  die("the BAM headers do not list the same contigs; all BAMs of a comparison must come from the same reference")
}
check_regions <- function(gr, what) {
  if (length(gr) == 0L) {
    die(paste0("the ", what, " set is empty (no ", what, " on a standard chromosome in the peak files / TxDb)"))
  }
  inbam <- as.character(seqnames(gr)) %in% names(bam_contigs[[1]])
  if (!any(inbam)) {
    die(paste0("none of the ", length(gr), " ", what, " regions is on a contig of the BAM header (peak/TxDb seqnames e.g. '",
               as.character(seqnames(gr))[1], "', BAM contigs e.g. '", names(bam_contigs[[1]])[1],
               "'): chromosome naming style differs"))
  }
  if (!all(inbam)) {
    cat("WARNING: ", sum(!inbam), " of ", length(gr), " ", what, " regions are on contigs absent from the BAM header ",
        "(their counts are zero)\n", sep = "", file = stderr())
  }
  invisible(NULL)
}

if ("peaks" %in% region_sets) {
  cat("Using the following files for control peaks: \n")
  print(ctrl.peaks)
  cat("Using the following files for treatment peaks: \n")
  print(cond.peaks)

  ctrl.peaks.gr <- import_and_merge_peaks(ctrl.peaks)
  cond.peaks.gr <- import_and_merge_peaks(cond.peaks)

  # Get the union set between condition and control peaksets
  my.regions <- union(ctrl.peaks.gr, cond.peaks.gr)
  check_regions(my.regions, "peak")
}

if ("promoters" %in% region_sets) {
  txdb <- AnnotationDbi::loadDb(file = txdb_file)

  # Generate Promoter Regions (T2T/hs1 compatible)
  all_promoters <- promoters(genes(txdb), upstream = prom_up, downstream = prom_down)
  all_promoters <- keepStandardChromosomes(all_promoters, pruning.mode = "coarse")
  check_regions(all_promoters, "promoter")
}

# ---- bins: counted once for the whole comparison (library sizes and TMM normalisation factors) -------------------------
cat("\nCounting reads into bins\n")
bins <- windowCounts(bams,
                     bin = TRUE,
                     width = bin.size,
                     param = param,
                     BPPARAM = BPPARAM)

# ---- one test per region set ------------------------------------------------------------------------------------------------
summary_list <- list()
norm_factors <- NULL
set_index <- 0L

for (set_label in region_sets) {
  set_index <- set_index + 1L

  if (set_label == "peaks") {
    cat("\nCounting reads over union set of peaks\n")
    data <- regionCounts(bams,
                         my.regions,
                         ext = NA,
                         param = param,
                         BPPARAM = BPPARAM)

    ## Compute normalization factors
    data <- normFactors(bins,
                        se.out = data,
                        weighted = F
                        )
  } else {
    cat("\nCounting reads over promoters\n")
    data <- regionCounts(bams, all_promoters, ext = NA, param = param, BPPARAM = BPPARAM)

    # Filter low-abundance regions
    filter.global <- filterWindowsGlobal(data, bins)
    cat("\nFiltering statistics for promoters:\n")
    summary(filter.global$filter) %>% print()

    fh_path <- diag_device("diag_promoter_filter", plot_width, plot_height)
    hist(filter.global$filter,
         main = name,
         breaks = 50,
         xlab = expression("Background"~"abundance" ~ (log[2]~"CPM")))
    abline(v = log2(min.fc), col = "red", lwd = 2)
    diag_close(fh_path)

    keep <- filter.global$filter > log2(min.fc)
    ### To investigate the effectiveness of our filtering strategy:
    cat("\nNumber of regions kept after filtering:\n")
    summary(keep) %>% print

    data <- data[keep, ]

    cat("\nComputing normalization factors for promoters\n\n")
    data <- normFactors(bins, se.out = data, weighted = F)
  }

  if (set_index == 1L) {
    normfacs <- data$norm.factors
    norm_factors <- normfacs

    ## Visualize the effects of normalization
    bin.ab <- scaledAverage(bins)
    adjc <- calculateCPM(bins, use.norm.factors = F)

    # Grid as legacy par(mfrow = c(1, 2)) per pair in pdf mode; in svg mode one svg with all pairs in a 2-column grid
    npairs <- dim(bins)[2] - 1L
    ncol_ma <- if (plot_format == "pdf") 2L else min(2L, npairs)
    nrow_ma <- if (plot_format == "pdf") 1L else ceiling(npairs / ncol_ma)
    ma_path <- diag_device("diag_ma", plot_width, plot_height * nrow_ma)
    par(mfrow = c(nrow_ma, ncol_ma), mar = c(5, 4, 5, 4))
    for (i in 2:dim(bins)[2]) {
      smoothScatter(bin.ab, adjc[, 1] - adjc[, i],
                    ylim = c(-6, 6),
                    main = paste(c(ids[1], "vs", ids[i]), sep = "\n"),
                    xlab = "Mean signal (logCPM) across samples (A)",
                    ylab = "Log2 Fold-Change (M)")

      abline(h = log2(normfacs[1] / normfacs[i]), col = "red")
    }
    diag_close(ma_path)
  }

  # STATISTICAL MODELING

  ## Convert RangedSummarizedExperiment object into a DGEList for modelling with edgeR.
  y <- asDGEList(data)

  cat("\nModel for ", set_label, " (", ncol(y), " libraries): ", if (use_blocks) "~ block + grouping" else "~ 0 + grouping",
      "\n", sep = "")
  tag <- if (set_label == "peaks") "" else "promoters_"
  res <- test_model(y, tag)

  ## Consolidate results
  if (set_label == "peaks") {
    merged <- overlapResults(data,
                             res$table,
                             regions = my.regions
                             )
  } else {
    merged <- overlapResults(data,
                             res$table,
                             regions = all_promoters
                             )
  }

  ## Inspect results
  if (table_sel == "best") {
    tabcom <- merged$best
  } else {
    tabcom <- merged$combined
  }

  summary_list[[length(summary_list) + 1L]] <- summary_rows(tabcom, set_label)

  ## Label regions

  ### Get regions
  ranges <- merged$regions

  ### Annotate
  ranges$PValue <- tabcom$PValue
  ranges$FDR <- tabcom$FDR
  ranges$logFC <- tabcom$rep.logFC
  ranges$score <- tabcom$rep.logFC

  ### Annotate genome build
  genome(ranges) <- build

  if (set_label == "peaks") {
    ### Make sure seqnames are in correct style
    seqnames(ranges) <- mapSeqlevels(as.character(seqnames(ranges)),
                                     style = style)

    ## Export peaks

    ### Convert GRanges to data frame
    bed <- ranges %>% as.data.frame()

    # Name peaks
    bed$name <- paste0("region", 1:nrow(bed))

    # Export as CSV
    write.csv(bed, paste0(name, ".csv"))

    # Get up and downregulated regions
    up <- filter(bed, logFC > lfc)
    down <- filter(bed, logFC < -lfc)

    # Convert the filtered data frame back to a GRanges object
    bed <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)

    # Export as BED
    rtracklayer::export.bed(bed, paste0(name, ".bed"))

    # Convert the filtered data frame back to a GRanges object
    up <- makeGRangesFromDataFrame(up, keep.extra.columns = TRUE)

    # Export as BED
    rtracklayer::export.bed(up, paste0(name, "_upregulated.bed"))

    # Convert the filtered data frame back to a GRanges object
    down <- makeGRangesFromDataFrame(down, keep.extra.columns = TRUE)

    # Export as BED
    rtracklayer::export.bed(down, paste0(name, "_downregulated.bed"))
  } else {
    # We extract it specifically so it survives the data.frame conversion
    gene_ids <- mcols(ranges)$gene_id

    ## Convert to data frame and fix columns
    bed <- as.data.frame(ranges)

    # Ensure gene_id is a named column in the data frame
    if (!"gene_id" %in% colnames(bed)) {
      bed$gene_id <- gene_ids
    }

    # Add a specific name column for BED format compatibility
    bed$name <- bed$gene_id

    # Remove promoters with NA logFCs
    bed <- filter(bed, !is.na(score))

    # Export as CSV (row names = gene ids, as legacy)
    write.csv(bed, paste0(name, "_promoters.csv"))

    # Get up and downregulated regions
    up <- filter(bed, logFC > lfc)
    down <- filter(bed, logFC < -lfc)

    # Convert the filtered data frame back to a GRanges object
    bed <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)

    # Export as BED
    rtracklayer::export.bed(bed, paste0(name, "_promoters.bed"))

    up <- makeGRangesFromDataFrame(up, keep.extra.columns = TRUE)
    rtracklayer::export.bed(up, paste0(name, "_promoters_up.bed"))

    down <- makeGRangesFromDataFrame(down, keep.extra.columns = TRUE)
    rtracklayer::export.bed(down, paste0(name, "_promoters_down.bed"))
  }
}

if (plot_format == "pdf") grDevices::dev.off()

# ---- outputs --------------------------------------------------------------------------------------------------------------
cat("\nNormalisation factors (TMM, weighted = FALSE):\n")
nf_df <- data.frame(
  sample_id = ids,
  norm_factor = as.numeric(norm_factors),
  lib_size = as.numeric(bins$totals),
  group = groups,
  block = if (use_blocks) blocks else "NA"
)
print(nf_df)
write.table(nf_df, paste0(name, ".norm_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summary_df <- do.call(rbind, summary_list)
write.table(summary_df, paste0(name, ".summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nSummary:\n")
print(summary_df)

rep_df <- data.frame(
  Comparison = name,
  Case = case_label,
  Control = control_label,
  Case_Replicates = as.integer(n_case),
  Control_Replicates = as.integer(n_ctrl),
  Total_Samples = length(bams),
  stringsAsFactors = FALSE
)
write.csv(rep_df, paste0(name, ".replicates.csv"), row.names = FALSE)
cat("\nReplicates written to: ", paste0(name, ".replicates.csv"), "\n", sep = "")
