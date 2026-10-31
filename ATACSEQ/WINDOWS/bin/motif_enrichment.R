#!/usr/bin/env Rscript
# motif_enrichment.R -- binned motif enrichment (monaLisa, binomial test) for ONE analysis: the differential regions of a
# comparison (diff), the regions annotated to a target gene vs all other regions (gene), or the regions overlapping a target
# set of regions (BED) vs the other regions (regions).
# Port of the legacy motif_analysis.Rmd (Scripts/, AV_ATACSEQ_2026, March 2026): chunk "loop and do motif analysis" (diff; the
# pooled motif_analysis_pooled.Rmd is the same code with comparison gd2_high_vs_low), chunk "motif analysis ST8SIA1 vs all
# other peaks" (gene), chunks "overlap atac peaks with ccres" and "motif analysis ccre vs non ccre" (regions), plus utils.R
# (compositionDiagonisticPlot and nucFreqDiagnosticPlot, copied verbatim below). The computation (sequences, composition
# plots, calcBinnedMotifEnrR, the motif selection, the heatmaps, the merge into the full table) is the legacy code; the I/O is
# the CLI, and the three modes share one copy of the selection and plotting code.
#
# CLI contract (called by the MOTIF_ENRICHMENT module; --key=value):
#   --mode=diff|gene|regions      REQUIRED. Which analysis (see the three sections below).
#   --name=NAME                   REQUIRED. Output prefix, ^[A-Za-z0-9][A-Za-z0-9_.-]*$. diff: the COMPARISON name (e.g.
#                                 BT54_high_vs_low; the region goes into the summary row). gene/regions: the motif_targets name.
#   --table=F.csv                 REQUIRED. diff: the differential table (NAME.csv or NAME_promoters.csv; needs logFC).
#                                 gene/regions: the all_sites annotation table of the universe comparison (DIFF_ANNOTATE).
#   --species="Homo sapiens"      OPTIONAL, default "Homo sapiens". species= of the JASPAR2024 query; default of --bsgenome.
#   --genome-build=hs1            OPTIONAL, default hs1. Used in the default --bsgenome.
#   --bsgenome=PKG                OPTIONAL. BSgenome package for getSeq. Default BSgenome.Hsapiens.UCSC.<build> (Mus musculus:
#                                 BSgenome.Mmusculus.UCSC.<build>). Missing package: exit 3.
#   --pwm-rds=F                   OPTIONAL. Cached PWMList (legacy Importable_Data/rds/pwms.rds). Blank = JASPAR2024 from the
#                                 package: getMatrixSet(JASPAR2024::db(JASPAR2024()), opts = list(matrixtype = "PWM", species)).
#   --min-score=10                OPTIONAL, default 10. min.score of calcBinnedMotifEnrR (PWM hit threshold).
#   --test=binomial               OPTIONAL, default binomial. test of calcBinnedMotifEnrR.
#   --max-kmer=3                  OPTIONAL, default 3. maxKmerSize of calcBinnedMotifEnrR (composition correction).
#   --background=otherBins        OPTIONAL, default otherBins. background of calcBinnedMotifEnrR (the other bin is the background).
#   --seed=123                    OPTIONAL, default 123. RNGseed of the BiocParallel parameter.
#   --top=40                      OPTIONAL, default 40. Motifs shown in the heatmaps when <= 2 or > top pass the primary threshold.
#   --padj-primary=2              OPTIONAL, default 2. -log10 adjusted P for the primary selection (legacy padj < 0.01).
#   --padj-relaxed=1              OPTIONAL, default 1. -log10 adjusted P for the relaxed selection (legacy padj < 0.1), used when <= 2
#                                 motifs pass the primary threshold.
#   --cpus=1                      OPTIONAL, default 1. 1: BiocParallel SerialParam(RNGseed = seed), exactly the legacy. > 1:
#                                 MulticoreParam(cpus, RNGseed = seed). NOT bit-reproducible against the legacy serial run.
#   --font=Helvetica              OPTIONAL, default Helvetica. Font family of every plot (svg: font-family of the text elements).
#   --plot-format=svg|pdf         OPTIONAL, default svg. svg: one svg per plot (svglite, editable <text>). pdf: the legacy pdf files
#                                 (same names with .pdf, legacy device arguments).
#   --heatmap-width=11            OPTIONAL, default 11 (inches). Width of both motif heatmaps.
#   --heatmap-height=9            OPTIONAL, default 9 (inches). Height of both motif heatmaps.
#   --gc-width=6                  OPTIONAL, default 6 (inches). Width of the GC fraction plot.
#   --gc-height=5                 OPTIONAL, default 5 (inches). Height of the GC fraction plot.
#   --dinuc-width=6               OPTIONAL, default 6 (inches). Width of the dinucleotide frequency plot.
#   --dinuc-height=6              OPTIONAL, default 6 (inches). Height of the dinucleotide frequency plot.
# diff mode only:
#   --region=promoters            REQUIRED. Region key (e.g. promoters, peaks); used in messages and in the skipped report.
#   --region-label=Promoters      OPTIONAL, default = --region. Text of the region column of the summary row (legacy "Promoters").
#   --lfc=0.585                   OPTIONAL, default 0.585. Bins: logFC > lfc (increased) and logFC < -lfc (decreased).
#   --min-bin-peaks=10            OPTIONAL, default 10. Skip rule: fewer regions than this in either bin.
# gene mode only:
#   --target=ST8SIA1              REQUIRED. Gene symbol matched in --gene-column (legacy: grepl, fixed, on geneId).
#   --gene-column=geneId          OPTIONAL, default geneId (annotation column).
#   --match=substring|exact       OPTIONAL, default substring (legacy). exact: the whole cell equals --target.
#   --target-label=               OPTIONAL, default "<target> peaks" (legacy bin label).
#   --other-label="Other peaks"   OPTIONAL, default "Other peaks".
#   --min-peaks=15                OPTIONAL, default 15. Skip rule: fewer regions than this in the target or the other bin.
# regions mode only:
#   --target=F.bed               REQUIRED. BED of the target regions (read with rtracklayer::import). The legacy set is the
#                                 ST8SIA1 cCREs lifted over to hs1 (st8sia1_ccres_hs1.rds); the liftover is a project pre-step.
#   --target-label="Target regions" OPTIONAL. Bin label (legacy "ST8SIA1 cCRE").
#   --other-label="Other peaks"   OPTIONAL.
#   --min-peaks=15                OPTIONAL, default 15. Skip rule: fewer overlapping regions than this in the target or other bin.
#
# Outputs (cwd, prefix = --name):
#   NAME_motif_analysis.csv          full table, one row per motif (motif.id, motif.name, then the enrichment and adjusted P
#                                    columns). diff: log2enr_down, log2enr_up, padj_down, padj_up (legacy merge). gene/regions:
#                                    log2enr_<level>, negLog10Padj_<level> (levels = alphabetical factor levels of the bins).
#   NAME_motif_heatmap_padj.svg      heatmap, rows alphabetical (no clustering) (--plot-format=pdf: .pdf)
#   NAME_motif_similarity_padj.svg   heatmap, rows clustered by PWM similarity (1 - Pearson correlation, average linkage)
#   NAME.gc_fraction.svg             GC fraction per bin (violin + box)
#   NAME.dinucleotide_frequency.svg  dinucleotide frequency heatmap (centred per bin)
#   NAME.motif_summary.csv           one row. diff: comparison, region, total_peaks, n_peaks_up, n_peaks_down, n_motifs_tested,
#                                    n_sig_padj_0.01, n_sig_padj_0.1, n_motifs_plotted (diff: the motifs passing the threshold,
#                                    BEFORE the topX cap; gene/regions: after the cap, as legacy). gene: comparison, region (all_sites),
#                                    total_peaks, n_peaks_target, n_peaks_other, n_motifs_tested, n_sig_padj_0.01, n_sig_padj_0.1,
#                                    n_motifs_plotted. regions: comparison, n_ccre_peaks, n_other_peaks, peak_width_bp, n_ccres_hs1,
#                                    n_motifs_tested, n_sig_padj_0.01, n_sig_padj_0.1, n_motifs_plotted (legacy column names).
#   NAME.skipped.txt                 ONLY when a bin has fewer regions than the minimum: the reason and the counts. No other output
#                                    is written in that case (exit 0).
#
# Exit codes: 0 ok (also when skipped); 2 usage error (also --min-bin-peaks / --min-peaks < 1); 3 known data problem: table or BED
# missing, a column missing, a BSgenome package or the JASPAR2024 package not installed, a --pwm-rds file missing, no region on
# a contig of the genome (seqname style), or a bin left empty after the out-of-bounds regions are dropped (see below).
# Regions that leave their chromosome after the resize to the common width (or are on a contig the genome lacks) cannot be read by
# getSeq (the legacy stopped with "trying to load regions beyond the boundaries"); they are dropped with a WARNING that gives the
# count, and total_peaks / n_peaks_* of the summary row count the regions actually analysed.
#
# ---- Question, model, assumptions -------------------------------------------------------------------------------------------
# Question: which TF motifs occur more often in one bin of regions than in the other bin? diff: increased accessibility
# (logFC > lfc) vs decreased (logFC < -lfc); gene: regions of the target gene vs all other regions; regions: regions overlapping
# the target set vs the rest.
# Model: for each motif, PWM hits (min.score) are counted over the regions of each bin and compared with a BINOMIAL test
# (test = "binomial": total hit counts, so several hits per region count), using the other bin as the background
# (background = "otherBins"; monaLisa calcBinnedMotifEnrR). log2enr = log2 enrichment of the bin; negLog10Padj = -log10 of the
# Benjamini-Hochberg adjusted P value (monaLisa). Significant: negLog10Padj > 2 (adjusted P < 0.01), relaxed > 1 (< 0.1).
# Composition control: maxKmerSize = 3, so the background sequences are sampled to match GC, dinucleotide and trinucleotide
# composition of the foreground (monaLisa). All regions are resized to one width (the mean width, centred), so the width does
# not confound the test.
# Assumptions: the other bin is an adequate background for the bin (no other confounder besides composition); the sampling of
# the background is fixed by --seed with a serial BiocParallel parameter (legacy RNGseed = 123), so the serial run is repeatable.
# The gene and regions modes test a set of regions, not a set of genes: a gene's regions are not independent of each other.
#
# ---- Statistical limits of the legacy method, ported unchanged (state them when reporting results) ---------------------------
# * Bins (diff mode) are cut on the estimated logFC only (> lfc, < -lfc); there is NO FDR filter. The bins therefore hold regions
#   with a large estimated fold change, significant or not. For comparisons without replication (fixed-dispersion LRT in
#   diff_accessibility.R) the extremes of logFC are partly noise. The result describes "regions that gained versus regions that
#   lost accessibility by more than lfc", not "significantly differential regions". Small bins (the skip rule is only 10
#   regions) give little power.
# * The reference is the other bin ("otherBins"): every statement is relative (Increased versus Decreased; target versus all other
#   regions of the universe table), never "enriched in accessible chromatin" in general.
# * monaLisa's binomial test is ONE-SIDED (P(hits >= observed), pbinom lower.tail = FALSE): negLog10P and negLog10Padj are
#   >= 0 and carry no sign. A large padj_up means enriched in the Increased bin; a motif depleted there shows up in the other
#   column. Padj is Benjamini-Hochberg over all motifs and both bins jointly (2 x number of motifs tests). The two columns of a
#   two-bin analysis mirror each other, so they are not independent tests.
# * Every PWM hit counts as one independent Bernoulli trial (trials = number of regions, success probability = background hits per
#   region; capped below 1 for motifs with more than one hit per region on average). Hits are clustered inside a region
#   (overlapping hits, homotypic clusters), and JASPAR holds many correlated PWMs of one TF family, so the p-values are
#   anti-conservative and the number of "significant motifs" counts redundant matrices. Read the ranking and log2enr, not the
#   count. A motif with no hit in either bin has NA (kept as NA in the table).
# * The heatmaps show the topX motifs by -log10 padj when <= 2 motifs pass the relaxed threshold or more than topX pass the primary
#   one, whether or not they are significant; n_sig_padj_* of the summary row give the number that actually passed.
# * Gene mode: "regions of the gene" are the regions whose ChIPseeker nearest-TSS gene (no distance limit) matches --target. The
#   default --match=substring also matches longer symbols (FOS: FOSB, FOSL1; JUN: JUNB, JUND; MYC: MYCN), use --match=exact for
#   such targets. The target regions of one gene are few and not independent (nearby regions share sequence and context).
# * Width: all regions are cut or extended to the mean width, so long peaks lose their flanks and short peaks gain sequence that is
#   not accessible. The composition control (kmers up to 3) removes GC / di- / tri-nucleotide differences, not other confounders.
# * motif.name is not unique in JASPAR2024 (84 of the 720 PWMs in the cached legacy set share one of 39 names): motif.id is the
#   key. motif_common.R still matches on motif.name (legacy) and warns.
#
# ---- Where the port deviates from the legacy, and why -----------------------------------------------------------------------
# * BPPARAM: --cpus = 1 gives exactly the legacy SerialParam(RNGseed). --cpus > 1 uses MulticoreParam(cpus, RNGseed); the
#   random streams are then not the legacy serial ones, so the output is not bit-reproducible against the legacy run.
# * The legacy diff chunk ended the calcBinnedMotifEnrR call with an empty trailing argument ("BPPARAM = ..., )"); it is dropped
#   (no effect on the call).
# * The legacy diff chunk drew the dinucleotide heatmap with draw() at construction time and printed the stored value later; here
#   nucFreqDiagnosticPlot() is called inside its own pdf() device, the same drawing.
# * The legacy diff chunk printed "range(SimMatSel)" (no effect) and the commented-out relaxed-LFC and skip code. Not ported.
# * The legacy diff chunk's min-bin check was commented out; the skip rule (--min-bin-peaks) is the one this port defines.
# * Gene and regions: the GC subtitle uses an ASCII hyphen ("All regions - N peaks") instead of the legacy em dash.
# * n_motifs_plotted (diff): the legacy summary row is built before the topX cap, so the value is the number of motifs passing
#   the threshold (e.g. 145 for BT54), although the legacy comment says "after topX cap". Ported as the code computes it; the
#   heatmaps themselves show at most topX motifs. gene and regions build their row after the cap (value = topX), as legacy.
# * PWMs: --pwm-rds gives the legacy cached PWMList (720 motifs). Without it, the JASPAR2024 package is queried (no network in
#   the test container: see the report for the outcome).
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
as_num <- function(x, key) {
  v <- suppressWarnings(as.numeric(x))
  if (length(v) != 1L || is.na(v)) die(paste0("--", key, " must be a number, got '", x, "'"), 2L)
  v
}

# ---- plots: one svg per plot (svglite, editable <text>, font family = --font); --plot-format=pdf = the legacy pdf ----
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
# One plot written to BASE.svg (--plot-format=svg) or BASE.pdf (pdf: the legacy device). Pair with plot_close(BASE).
plot_open <- function(base, width, height) {
  if (plot_format == "svg") {
    open_svg(paste0(base, ".svg"), width, height)
  } else if (tolower(plot_font) %in% c("helvetica", "arial", "sans")) {
    grDevices::pdf(file = paste0(base, ".pdf"), width = width, height = height, family = "Helvetica")
  } else {
    grDevices::cairo_pdf(filename = paste0(base, ".pdf"), width = width, height = height, family = plot_font)
  }
  invisible(base)
}
plot_close <- function(base) {
  if (plot_format == "svg") close_svg(paste0(base, ".svg")) else grDevices::dev.off()
  invisible(NULL)
}

# Legacy utils.R (verbatim): GC composition diagnostic plot
compositionDiagonisticPlot <- function(bins, seqs, title = "", subtitle = ""){
  # Compute GC fraction for each sequence
  gc_frac <- letterFrequency(seqs, letters = c("G","C"), as.prob = TRUE)
  gc_frac <- rowSums(gc_frac)  # fraction of G or C per sequence

  df_gc <- data.frame(
    Bin = bins,
    GC = gc_frac * 100 # convert to percentage
  )

  binCols <- getColsByBin(bins)


  p <- ggplot(data = df_gc,
              mapping = aes(x = GC, y = Bin)) +
    geom_violin(trim = TRUE,
                mapping = aes(fill = Bin),
                show.legend = FALSE) +
    geom_boxplot(width = 0.05,
                 outlier.size = 0.8,
                 show.legend = FALSE) +
    ggtitle(label = title, subtitle = subtitle) +
    scale_fill_manual(values = attr(binCols, "cols")) +
    scale_x_continuous(name = "GC fraction",
                       labels = scales::percent_format(scale = 1)) +
    labs(y = element_blank()) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))

  return(p)
}

# Legacy utils.R (verbatim): dinucleotide frequency plot
nucFreqDiagnosticPlot <- function(bins, seqs, title = "", threshold = 0.005){
  # compute dinucleotide frequencies
  dnf <- Biostrings::oligonucleotideFrequency(seqs, width = 2, as.prob = TRUE)
  dnfsplit <- lapply(split.data.frame(dnf, bins), colMeans)
  dnfmat <- t(do.call(rbind, dnfsplit))

  # scale row-wise
  dnfmat_scaled <- t(scale(t(dnfmat), center = TRUE, scale = FALSE))

  # determine symmetric color scale
  max_obs <- max(abs(dnfmat_scaled))
  lim <- max(max_obs, threshold)

  cols <- circlize::colorRamp2(
    breaks = c(-lim, 0, lim),
    colors = c("blue", "white", "red")
  )

  hm <- ComplexHeatmap::Heatmap(
    dnfmat_scaled,
    column_names_centered = TRUE,
    column_names_rot = 0,
    column_title = title,
    right_annotation = ComplexHeatmap::rowAnnotation(
      OverallFreq = rowMeans(dnfmat),
      bar2 = anno_barplot(rowMeans(dnfmat)),
      col = list(OverallFreq = circlize::colorRamp2(
        range(rowMeans(dnfmat)), c("white", "grey50"))),
      show_annotation_name = FALSE,
      annotation_label = list(
        OverallFreq = "Average\nrelative\nfrequency\nacross bins")),
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    col = cols,
    name = "Δ mean\nrelative\nfrequency"
  )

  draw(hm)
}

allowed <- c("mode", "name", "table", "species", "genome-build", "bsgenome", "pwm-rds", "min-score", "test",
             "max-kmer", "background", "seed", "top", "padj-primary", "padj-relaxed", "cpus",
             "font", "plot-format", "heatmap-width", "heatmap-height", "gc-width", "gc-height",
             "dinuc-width", "dinuc-height",
             "region", "region-label", "lfc", "min-bin-peaks",
             "target", "gene-column", "match", "target-label", "other-label", "min-peaks")
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) == 0L || any(argv %in% c("-h", "--help"))) {
  cat("usage: motif_enrichment.R --mode=diff|gene|regions --name=NAME --table=F.csv [mode options] [--pwm-rds=F]\n",
      file = stderr())
  quit(save = "no", status = 2L)
}
args <- parse_args(argv, allowed)

mode <- get_req(args, "mode")
if (!(mode %in% c("diff", "gene", "regions"))) die("--mode must be diff, gene or regions", 2L)
name <- get_req(args, "name")
if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", name)) die(paste0("--name must match ^[A-Za-z0-9][A-Za-z0-9_.-]*$, got '", name, "'"), 2L)
table_file <- get_req(args, "table")

species   <- get_opt(args, "species", "Homo sapiens")
build     <- get_opt(args, "genome-build", "hs1")
if (species == "Mus musculus") {
  default_bsgenome <- paste0("BSgenome.Mmusculus.UCSC.", build)
} else {
  default_bsgenome <- paste0("BSgenome.Hsapiens.UCSC.", build)
}
bsgenome_pkg <- get_opt(args, "bsgenome", default_bsgenome)
pwm_rds      <- get_opt(args, "pwm-rds", "")

min_score   <- as_num(get_opt(args, "min-score", "10"), "min-score")
test        <- get_opt(args, "test", "binomial")
max_kmer    <- as.integer(as_num(get_opt(args, "max-kmer", "3"), "max-kmer"))
background  <- get_opt(args, "background", "otherBins")
seed        <- as.integer(as_num(get_opt(args, "seed", "123"), "seed"))
topX        <- as.integer(as_num(get_opt(args, "top", "40"), "top"))
padj_primary <- as_num(get_opt(args, "padj-primary", "2"), "padj-primary")
padj_relaxed <- as_num(get_opt(args, "padj-relaxed", "1"), "padj-relaxed")
cpus        <- as.integer(as_num(get_opt(args, "cpus", "1"), "cpus"))
if (is.na(cpus) || cpus < 1L) die("--cpus must be a positive integer", 2L)
plot_font    <- get_opt(args, "font", "Helvetica")
plot_format  <- tolower(get_opt(args, "plot-format", "svg"))
heat_width   <- as_num(get_opt(args, "heatmap-width", "11"), "heatmap-width")
heat_height  <- as_num(get_opt(args, "heatmap-height", "9"), "heatmap-height")
gc_width     <- as_num(get_opt(args, "gc-width", "6"), "gc-width")
gc_height    <- as_num(get_opt(args, "gc-height", "5"), "gc-height")
dinuc_width  <- as_num(get_opt(args, "dinuc-width", "6"), "dinuc-width")
dinuc_height <- as_num(get_opt(args, "dinuc-height", "6"), "dinuc-height")
if (!(plot_format %in% c("svg", "pdf"))) die("--plot-format must be svg or pdf", 2L)
if (min(heat_width, heat_height, gc_width, gc_height, dinuc_width, dinuc_height) <= 0) {
  die("--heatmap-*, --gc-* and --dinuc-* sizes must be > 0", 2L)
}

if (mode == "diff") {
  region       <- get_req(args, "region")
  region_label <- get_opt(args, "region-label", region)
  lfc_thres    <- as_num(get_opt(args, "lfc", "0.585"), "lfc")
  min_bin      <- as.integer(as_num(get_opt(args, "min-bin-peaks", "10"), "min-bin-peaks"))
} else {
  target       <- get_req(args, "target")
  min_peaks    <- as.integer(as_num(get_opt(args, "min-peaks", "15"), "min-peaks"))
  if (mode == "gene") {
    gene_column  <- get_opt(args, "gene-column", "geneId")
    match_mode   <- get_opt(args, "match", "substring")
    if (!(match_mode %in% c("substring", "exact"))) die("--match must be substring or exact", 2L)
    target_label <- get_opt(args, "target-label", paste0(target, " peaks"))
  } else {
    target_label <- get_opt(args, "target-label", "Target regions")
  }
  other_label  <- get_opt(args, "other-label", "Other peaks")
}

if (mode == "diff") {
  if (is.na(min_bin) || min_bin < 1L) die("--min-bin-peaks must be an integer >= 1 (a bin cannot be empty)", 2L)
} else if (is.na(min_peaks) || min_peaks < 1L) {
  die("--min-peaks must be an integer >= 1 (a bin cannot be empty)", 2L)
}
if (!file.exists(table_file)) die(paste0("table not found: ", table_file), 3L)
if (mode == "regions" && !file.exists(target)) die(paste0("target BED not found: ", target), 3L)
if (nzchar(pwm_rds) && !file.exists(pwm_rds)) die(paste0("PWM rds not found: ", pwm_rds), 3L)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(SummarizedExperiment)
  library(Biostrings)
  library(BiocParallel)
  library(TFBSTools)
  library(monaLisa)
  library(ComplexHeatmap)
  library(circlize)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(rtracklayer)
  library(grid)
})
for (p in c(bsgenome_pkg, "JASPAR2024")) {
  if (!requireNamespace(p, quietly = TRUE)) die(paste0("package not installed in this image: ", p), 3L)
}
genome <- BSgenome::getBSgenome(bsgenome_pkg)

# PWMs: cached PWMList (legacy rds/pwms.rds), else the JASPAR2024 package (legacy query)
if (nzchar(pwm_rds)) {
  pwm <- readRDS(file = pwm_rds)
} else {
  jaspar_obj <- JASPAR2024::JASPAR2024()
  pwm <- getMatrixSet(JASPAR2024::db(jaspar_obj),
                      opts = list(matrixtype = "PWM",
                                  species = species)
  )
}

# BiocParallel: serial (legacy) unless --cpus > 1
bpparam <- if (cpus > 1L) MulticoreParam(cpus, RNGseed = seed) else SerialParam(RNGseed = seed)

# Drop the (resized) regions that getSeq cannot read: contig not in the genome, start < 1 or end beyond the chromosome. Returns
# the kept regions and their bins (levels not used any more are dropped; nothing changes when no region is dropped).
drop_out_of_bounds <- function(gr, bins) {
  sl <- seqlengths(genome)[as.character(seqnames(gr))]
  ok <- !is.na(sl) & start(gr) >= 1L & end(gr) <= sl
  if (!any(ok)) {
    die(paste0("none of the ", length(gr), " regions lies on a contig of ", bsgenome_pkg, " (seqnames e.g. '",
               as.character(seqnames(gr))[1], "', genome contigs e.g. '", seqnames(genome)[1], "'): chromosome naming style differs"), 3L)
  }
  if (!all(ok)) {
    cat("WARNING: ", sum(!ok), " of ", length(gr), " regions are dropped: after the resize to ", unique(width(gr)),
        " bp they leave their chromosome or lie on a contig absent from ", bsgenome_pkg, "\n", sep = "", file = stderr())
  }
  bins_kept <- droplevels(bins[ok])
  if (nlevels(bins_kept) < 2L) {
    die(paste0("after dropping the out-of-bounds regions one bin is empty (bins left: ", paste(levels(bins_kept), collapse = ", "), ")"), 3L)
  }
  list(gr = gr[ok], bins = bins_kept)
}

skip_run <- function(reason) {
  writeLines(c(paste0("SKIPPED: ", reason),
               "Motif enrichment was not run for this analysis (fewer regions than the minimum in one bin)."),
             paste0(name, ".skipped.txt"))
  cat(paste0("\nSkipped ", name, ": ", reason, "\n"))
  quit(save = "no", status = 0L)
}

# Shared: significance counts, motif selection (legacy rules) and the heatmaps
sig_counts <- function(se) {
  c(primary = sum(apply(assay(se, "negLog10Padj"), 1,
                        function(x) max(abs(x), 0, na.rm = TRUE)) > padj_primary),
    relaxed = sum(apply(assay(se, "negLog10Padj"), 1,
                        function(x) max(abs(x), 0, na.rm = TRUE)) > padj_relaxed))
}

select_motifs <- function(se) {
  # Primary threshold: padj < 0.01 (-log10 > 2)
  sel <- apply(assay(se, "negLog10Padj"), 1,
               function(x) max(abs(x), 0, na.rm = TRUE)) > padj_primary
  cat(paste0("\n# of motifs with adj. P-value below the primary threshold: ", sum(sel), "\n"))
  seSel <- se[sel, ]

  # Relax to padj < 0.1 if too few pass the primary threshold
  if (sum(sel) <= 2) {
    sel <- apply(assay(se, "negLog10Padj"), 1,
                 function(x) max(abs(x), 0, na.rm = TRUE)) > padj_relaxed
    cat(paste0("# of motifs with adj. P-value below the relaxed threshold: ", sum(sel), "\n"))
    seSel <- se[sel, ]
  }
  # Legacy diff summary (n_motifs_plotted) counts these, before the topX cap below
  n_selected <- nrow(seSel)

  # Cap at topX if too many pass or too few pass either threshold
  if (sum(sel) <= 2 | sum(sel) > topX) {
    cat(paste0("\nPlotting Top ", topX, " significant motifs\n"))
    maxVals <- apply(assay(se, "negLog10Padj"), 1,
                     function(x) max(abs(x), na.rm = TRUE))
    ord <- order(maxVals, decreasing = TRUE)
    seSel <- se[ord[seq_len(min(topX, length(ord)))], ]
  }
  list(seSel = seSel, n_selected = n_selected)
}

plot_motif_heatmaps <- function(seSel, name) {
  # Heatmap 1: alphabetical order, no clustering
  plot_open(paste0(name, "_motif_heatmap_padj"), heat_width, heat_height)
  plotMotifHeatmaps(x = seSel,
                    which.plots = c("log2enr", "negLog10Padj"),
                    width = 1.5,
                    maxEnr = max(max(abs(seSel@assays@data@listData$log2enr)),
                                 2),
                    maxSig = max(quantile(seSel@assays@data@listData$negLog10Padj,
                                          0.9),
                                 10),
                    row_names_gp = gpar(fontsize = 9),
                    # Sort rows alphabetically by motif name within each cluster
                    row_order = order(rowData(seSel)$motif.name),
                    cluster = FALSE,
                    show_bin_legend = TRUE,
                    show_seqlogo = TRUE)
  plot_close(paste0(name, "_motif_heatmap_padj"))

  # Heatmap 2: clustered by motif similarity
  SimMatSel <- motifSimilarity(rowData(seSel)$motif.pfm)

  # Create hclust object, similarity defined by 1 - Pearson correlation
  hcl <- hclust(as.dist(1 - SimMatSel), method = "average")

  plot_open(paste0(name, "_motif_similarity_padj"), heat_width, heat_height)
  plotMotifHeatmaps(x = seSel,
                    which.plots = c("log2enr", "negLog10Padj"),
                    width = 1.5,
                    cluster = hcl,
                    maxEnr = max(max(abs(seSel@assays@data@listData$log2enr)),
                                 2),
                    maxSig = max(quantile(seSel@assays@data@listData$negLog10Padj,
                                          0.9),
                                 10),
                    row_names_gp = gpar(fontsize = 9),
                    # Sort rows alphabetically by motif name within each cluster
                    row_order = order(rowData(seSel)$motif.name),
                    show_dendrogram = TRUE,
                    show_seqlogo = TRUE,
                    show_bin_legend = TRUE,
                    show_motif_GC = TRUE)
  plot_close(paste0(name, "_motif_similarity_padj"))
}

write_composition_plots <- function(bins, seqs, title, subtitle, name) {
  gc_plot <- compositionDiagonisticPlot(bins,
                                        seqs,
                                        title = title,
                                        subtitle = subtitle)
  plot_open(paste0(name, ".gc_fraction"), gc_width, gc_height)
  print(gc_plot)
  plot_close(paste0(name, ".gc_fraction"))

  plot_open(paste0(name, ".dinucleotide_frequency"), dinuc_width, dinuc_height)
  invisible(nucFreqDiagnosticPlot(bins, seqs, title = title))
  plot_close(paste0(name, ".dinucleotide_frequency"))
}

run_enrichment <- function(seqs, bins) {
  calcBinnedMotifEnrR(
    seqs  = seqs,
    bins  = bins,
    pwmL  = pwm,
    min.score = min_score,
    test  = test,           # "binomial": total hits per motif, allows >1 hit per region
    maxKmerSize = max_kmer, # corrects up to tri-nucleotides (GC/dinuc/tri)
    background = background,
    BPPARAM = bpparam
  )
}

# ---- diff mode: legacy "loop and do motif analysis" --------------------------------------------------------------------
if (mode == "diff") {
  peak_df <- read.csv(table_file)
  if (!("logFC" %in% colnames(peak_df))) die(paste0("no logFC column in ", table_file), 3L)

  # Subset peaks
  peaks_up <- filter(peak_df, logFC > lfc_thres)
  peaks_down <- filter(peak_df, logFC < -lfc_thres)

  cat(paste0("\nBinning peaks based on LFC threshold: ", lfc_thres, "\n"))
  down_tit <- paste0("LFC < ", -lfc_thres)
  up_tit <- paste0("LFC > ", lfc_thres)
  peak_counts <- table(rep(c(down_tit, up_tit),
            c(nrow(peaks_down), nrow(peaks_up))))
  print(peak_counts)

  if (nrow(peaks_up) < min_bin || nrow(peaks_down) < min_bin) {
    skip_run(paste0(name, " ", region, ": increased ", nrow(peaks_up),
                    ", decreased ", nrow(peaks_down), " regions (minimum ", min_bin, " per bin)"))
  }

  # Convert to GRanges
  peaks_up <- GenomicRanges::makeGRangesFromDataFrame(peaks_up, keep.extra.columns = TRUE)
  peaks_down <- GenomicRanges::makeGRangesFromDataFrame(peaks_down, keep.extra.columns = TRUE)

  peaks_fixed <- c(peaks_down, peaks_up)
  peaks_fixed <- resize(peaks_fixed,
                        width = round(mean(width(peaks_fixed)),0),
                        fix = "center")

  # Define binning vector
  bins <- rep(c("Decreased Accessibility", "Increased Accessibility"),
               c(length(peaks_down), length(peaks_up)))
  bins <- factor(bins)

  bounded <- drop_out_of_bounds(peaks_fixed, bins)
  peaks_fixed <- bounded$gr
  bins <- bounded$bins

  cat("\n")
  print(table(bins))

  # Check for differences in GC composition between bins
  seqs <- getSeq(genome, peaks_fixed)
  seqs_clean <- Biostrings::replaceAmbiguities(seqs, new="N")
  seqs <- seqs_clean

  write_composition_plots(bins, seqs,
                          title = gsub("_", " ", name, fixed = TRUE),
                          subtitle = paste0(region_label, " - ", length(peaks_fixed), "peaks"),
                          name = name)

  ## Run binned motif enrichment (composition-controlled)
  se <- run_enrichment(seqs, bins)

  # Get full table of motifs
  df <- rowData(se)
  df <- df[,c("motif.id", "motif.name")]

  # Get log2 enrichment for all motifs
  logdf <- se@assays@data@listData$log2enr
  colnames(logdf) <- c("log2enr_down", "log2enr_up")

  # Get negLog10 padj for all motifs
  padjdf <- se@assays@data@listData$negLog10Padj
  colnames(padjdf) <- c("padj_down", "padj_up")

  # Merge enrichment data and padj data
  df_full <- merge(df, logdf, by = "row.names")
  df_full <- df_full %>% column_to_rownames(var = "Row.names")
  df_full <- merge(df_full, padjdf, by = "row.names")[,-1]

  # Export full table
  write.csv(x = df_full, file = paste0(name, "_motif_analysis.csv"), row.names = FALSE)

  sel_res <- select_motifs(se)
  seSel <- sel_res$seSel
  sig <- sig_counts(se)

  # --- Accumulate a summary row for this comparison x region ---
  summary_row <- data.frame(
    comparison       = name,
    region           = region_label,
    # Total peaks analysed (both bins combined)
    total_peaks      = length(peaks_fixed),
    n_peaks_up       = sum(bins == "Increased Accessibility"),
    n_peaks_down     = sum(bins == "Decreased Accessibility"),
    # Total motifs tested
    n_motifs_tested  = nrow(se),
    # How many motifs cleared each significance bar
    n_sig_padj_0.01  = unname(sig["primary"]),
    n_sig_padj_0.1   = unname(sig["relaxed"]),
    # Legacy: the count of motifs passing the threshold, before the topX cap (see select_motifs)
    n_motifs_plotted = sel_res$n_selected,
    stringsAsFactors = FALSE
  )

  plot_motif_heatmaps(seSel, name)
}

# ---- gene mode: legacy "motif analysis ST8SIA1 vs all other peaks" ---------------------------------------------------
if (mode == "gene") {
  anno_df <- read.csv(table_file, stringsAsFactors = FALSE)
  if (!(gene_column %in% colnames(anno_df))) die(paste0("no column '", gene_column, "' in ", table_file), 3L)
  cat("Loaded annotation table:", nrow(anno_df), "regions\n")

  gene_vals <- as.character(anno_df[[gene_column]])
  is_target <- if (match_mode == "exact") {
    !is.na(gene_vals) & gene_vals == target
  } else {
    grepl(target, gene_vals, fixed = TRUE)
  }

  peaks_target <- anno_df[is_target, ]
  peaks_other  <- anno_df[!is_target, ]

  cat(paste0("\nPeaks mapped to ", target, ": ", nrow(peaks_target), "\n"))
  cat(paste0("All other peaks: ",               nrow(peaks_other),  "\n"))

  if (nrow(peaks_target) < min_peaks | nrow(peaks_other) < min_peaks) {
    skip_run(paste0("target ", target, " has ", nrow(peaks_target), " regions and the other bin ",
                    nrow(peaks_other), " (minimum ", min_peaks, " per bin)"))
  }

  # Convert both groups to GRanges and resize to uniform width
  peaks_target_gr <- GenomicRanges::makeGRangesFromDataFrame(
    peaks_target, keep.extra.columns = TRUE
  )
  peaks_other_gr  <- GenomicRanges::makeGRangesFromDataFrame(
    peaks_other,  keep.extra.columns = TRUE
  )

  # Combine and resize all peaks to the same width (centred on peak midpoint)
  peaks_all <- c(peaks_target_gr, peaks_other_gr)
  peaks_all <- resize(peaks_all,
                      width = round(mean(width(peaks_all)), 0),
                      fix   = "center")

  # Bin vector ordered to match the concatenation above (target first)
  bins <- factor(
    rep(c(target_label, other_label),
        c(length(peaks_target_gr), length(peaks_other_gr)))
  )

  bounded <- drop_out_of_bounds(peaks_all, bins)
  peaks_all <- bounded$gr
  bins <- bounded$bins

  cat("\nBin counts:\n")
  print(table(bins))

  seqs       <- getSeq(genome, peaks_all)
  seqs_clean <- Biostrings::replaceAmbiguities(seqs, new = "N")
  seqs       <- seqs_clean

  write_composition_plots(bins, seqs,
                          title    = paste0(target, " vs Other Peaks"),
                          subtitle = paste0("All regions - ", length(peaks_all), " peaks"),
                          name     = name)

  se <- run_enrichment(seqs, bins)

  # Export the full motif results table
  df      <- as.data.frame(rowData(se)[, c("motif.id", "motif.name")])
  logdf   <- as.data.frame(se@assays@data@listData$log2enr)
  padjdf  <- as.data.frame(se@assays@data@listData$negLog10Padj)

  # Column names reflect the two bins (monaLisa orders columns to match the alphabetical factor levels)
  colnames(logdf)  <- paste0("log2enr_",  levels(bins))
  colnames(padjdf) <- paste0("negLog10Padj_", levels(bins))

  df_full <- cbind(df, logdf, padjdf)

  write.csv(x = df_full, file = paste0(name, "_motif_analysis.csv"), row.names = FALSE)

  seSel <- select_motifs(se)$seSel
  sig <- sig_counts(se)

  summary_row <- data.frame(
    comparison       = name,
    region           = "all_sites",
    total_peaks      = length(peaks_all),
    n_peaks_target   = sum(bins == target_label),
    n_peaks_other    = sum(bins == other_label),
    n_motifs_tested  = nrow(se),
    n_sig_padj_0.01  = unname(sig["primary"]),
    n_sig_padj_0.1   = unname(sig["relaxed"]),
    n_motifs_plotted = nrow(seSel),
    stringsAsFactors = FALSE
  )

  plot_motif_heatmaps(seSel, name)
}

# ---- regions mode: legacy "overlap atac peaks with ccres" + "motif analysis ccre vs non ccre" -------------------------
if (mode == "regions") {
  target_gr <- rtracklayer::import(target)
  cat("Target regions:", length(target_gr), "\n")

  anno_df <- read.csv(table_file, stringsAsFactors = FALSE)

  # Remove columns whose names clash with reserved GRanges slot names before
  # passing to makeGRangesFromDataFrame ("width", "strand")
  cols_to_drop <- c("width", "strand")

  peaks_gr <- makeGRangesFromDataFrame(
    anno_df[, !colnames(anno_df) %in% cols_to_drop],
    keep.extra.columns = TRUE,
    seqnames.field     = "seqnames",
    start.field        = "start",
    end.field          = "end",
    ignore.strand      = TRUE
  )

  # findOverlaps returns a Hits object: queryHits = peak index, subjectHits = target index
  hits <- findOverlaps(peaks_gr, target_gr, ignore.strand = TRUE)

  # Boolean vector: TRUE for peaks that overlap at least one target region
  overlaps_ccre <- seq_along(peaks_gr) %in% queryHits(hits)

  peaks_ccre       <- peaks_gr[ overlaps_ccre]
  peaks_non_ccre   <- peaks_gr[!overlaps_ccre]

  cat(paste0("Peaks overlapping the target regions:     ", length(peaks_ccre),     "\n",
             "Peaks NOT overlapping the target regions: ", length(peaks_non_ccre), "\n"))

  if (length(peaks_ccre) < min_peaks | length(peaks_non_ccre) < min_peaks) {
    skip_run(paste0(length(peaks_ccre), " regions overlap the target and ", length(peaks_non_ccre),
                    " do not (minimum ", min_peaks, " per bin)"))
  }

  # Resize peaks to uniform width
  peaks_all_gr <- c(peaks_ccre, peaks_non_ccre)
  mean_width   <- round(mean(width(peaks_all_gr)), 0)

  peaks_all_fixed <- resize(peaks_all_gr, width = mean_width, fix = "center")

  cat("Peaks resized to uniform width:", mean_width, "bp\n")

  # Bin labels
  bins <- factor(
    rep(c(target_label, other_label),
        c(length(peaks_ccre), length(peaks_non_ccre)))
  )

  bounded <- drop_out_of_bounds(peaks_all_fixed, bins)
  peaks_all_fixed <- bounded$gr
  bins <- bounded$bins

  cat("\nBin counts:\n")
  print(table(bins))

  # Extract sequences
  seqs <- getSeq(genome, peaks_all_fixed)
  seqs <- Biostrings::replaceAmbiguities(seqs, new = "N")

  write_composition_plots(bins, seqs,
                          title    = paste0(target_label, " vs ", other_label),
                          subtitle = paste0(length(peaks_all_fixed), " total peaks"),
                          name     = name)

  se <- run_enrichment(seqs, bins)

  # Export full results table
  df     <- as.data.frame(rowData(se)[, c("motif.id", "motif.name")])
  logdf  <- as.data.frame(se@assays@data@listData$log2enr)
  padjdf <- as.data.frame(se@assays@data@listData$negLog10Padj)

  # monaLisa orders columns by factor level; label them explicitly
  colnames(logdf)  <- paste0("log2enr_",      levels(bins))
  colnames(padjdf) <- paste0("negLog10Padj_", levels(bins))

  df_full <- cbind(df, logdf, padjdf)

  write.csv(df_full,
            file      = paste0(name, "_motif_analysis.csv"),
            row.names = FALSE)

  seSel <- select_motifs(se)$seSel
  sig <- sig_counts(se)

  summary_row <- data.frame(
    comparison          = name,
    n_ccre_peaks        = sum(bins == target_label),
    n_other_peaks       = sum(bins == other_label),
    peak_width_bp       = mean_width,
    n_ccres_hs1         = length(target_gr),
    n_motifs_tested     = nrow(se),
    n_sig_padj_0.01     = unname(sig["primary"]),
    n_sig_padj_0.1      = unname(sig["relaxed"]),
    n_motifs_plotted    = nrow(seSel),
    stringsAsFactors    = FALSE
  )

  plot_motif_heatmaps(seSel, name)
}

write.csv(summary_row, file = paste0(name, ".motif_summary.csv"), row.names = FALSE)
cat(paste0("\nAll outputs written for: ", name, " (mode ", mode, ")\n"))
