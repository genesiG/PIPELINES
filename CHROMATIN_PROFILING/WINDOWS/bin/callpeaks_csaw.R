#!/usr/bin/env Rscript
# callpeaks_csaw.R -- target-vs-control (IgG) peak calling with csaw/edgeR for ONE peak set.
#
# CLI contract (called by the CSAW module; thresholds come from peaksets.csv cells, else config defaults, resolved by the pipeline):
#   --peakset=<name>              peak set label (BED names, logs)
#   --target-bams=<comma list>    target (antibody) BAMs, final filtered, coordinate-sorted + indexed
#   --control-bams=<comma list>   control (IgG) BAMs
#   --target-ids=<comma list>     sample ids of the target BAMs (same order)
#   --control-ids=<comma list>    sample ids of the control BAMs
#   --window=<int> --spacing=<int>   small window width / spacing (config csaw_window / csaw_spacing; 150 / 50)
#   --bin=<int>                   background bin width (config csaw_bin; 10000)
#   --small-filt=<num> --large-filt=<num>   log2 enrichment over global background; windows are KEPT when score > cutoff
#   --merge-bp=<int>              merge tolerance (bp) of windows into regions
#   --max-width=<int|NA>          maximum merged region width (mergeWindows max.width); NA/0/empty = no limit
#   --fdr=<num>                   region-level FDR threshold (peak when FDR <= fdr)
#   --lfc=<num>                   minimum log2 fold change of the region's representative window (see "lfc decision")
#   --mapq=<int> --max-frag=<int> csaw readParam(minq=, max.frag=); pass bowtie2_maxins as --max-frag (500 if unset)
#   --norm-bins=background|all    bins used to compute the TMM factors (see "Normalisation bins")
#   --dedup=TRUE|FALSE            csaw readParam(dedup=); TRUE drops reads flagged 0x400. With duplicate-removed BAMs (pipeline
#                                 default) both values give identical counts. NEVER use FALSE on markdup-flagged BAMs in which
#                                 only the IgG was deduplicated (targets would keep duplicates, IgG not).
#   --cpus=<int>                  workers for read counting (BiocParallel MulticoreParam over BAM files)
#   --out-prefix=<name>
#   OPTIONAL: --large-window=2000 --large-spacing=500 (the second, wide window scale of the legacy two-scale design),
#             --pe=both|none (default both; none = single-end with cross-correlation fragment length),
#             --bg-max-enrich=1 (log2 units, background-bin definition below),
#             --count-filter=10 (csaw windowCounts filter = minimum summed count over all libraries for a small or large WINDOW;
#                 the background bins keep the csaw default 10 so the TMM reference does not move with this knob),
#             --lfc-from=combined|best (default combined = legacy rule; see "lfc decision"),
#             --disjoint-regions=true|false (default false = csaw's own merge, regions may overlap; see "Disjoint regions")
# Outputs:
#   <name>.peaks.bed                 chr, start, end, peak_<i>, score, ".", logFC   (score = min(1000, round(-10*log10(FDR))),
#                                    FDR floored at 1e-300; logFC = the log2 FC (Target/Control) used by the peak rule:
#                                    best-window with --lfc-from=best, representative window of the combined test with =combined)
#   <name>.csaw_results.tsv.gz       ALL merged regions: chr,start,end,width,num_tests,num_up,num_down,PValue,FDR,direction,
#                                    rep_logFC,best_logFC,best_PValue,best_FDR,is_fdr,is_peak
#   <name>.mds.svg                   MDS of log-CPM (replicate QC): top 100/500/1000/5000 windows, small and large scale
#   <name>.rds                       atomic RDS: parameters, norm factors, library sizes, all merged regions
#   <name>.norm_sensitivity.tsv      TMM factors, #regions, #FDR-significant, #peaks under norm-bins=background AND all
#                                    (always computed so the effect of the setting can be judged); one row per
#                                    (setting, sample); run-level columns repeat on each row. Last columns: lfc_from, disjoint_regions,
#                                    count_filter and the overlap diagnostic of csaw's OWN merge (n_regions_csaw_merge,
#                                    n_overlapping_regions_csaw_merge, n_peaks_csaw_merge, n_overlapping_peaks_csaw_merge: number of
#                                    regions / peaks that overlap at least one other region / peak); with --disjoint-regions=true
#                                    n_regions / n_peaks describe the disjoint regions and the *_csaw_merge columns the csaw ones
# Exit codes: 0 ok; 2 usage error; 3 data problem (missing BAM/index, no controls/targets, too few background bins).
#
# ---- Question, model, assumptions -------------------------------------------------------------------------------
# Question: which genomic regions carry more fragments in the target libraries than in the matched IgG libraries, beyond
# replicate variability? Counts in sliding windows (two scales: window/spacing and large-window/large-spacing) are modelled
# jointly for all libraries with an edgeR quasi-likelihood NB GLM, ~0 + group (Target, Control), test Target - Control
# (H0: log2 FC = 0). Windows are pre-filtered by their average abundance relative to the global background
# (csaw::filterWindowsGlobal against the `--bin` bins; independent of the test statistic), tested separately per scale, and
# combined into regions by csaw::mergeResultsList (equiweight Simes combination, region FDR by BH over regions).
# Needs >= 1 target and >= 1 control; with only one library per group (no residual df) the NB dispersion cannot be
# estimated: a fixed common dispersion 0.05 with a likelihood-ratio test is used (legacy behaviour; the FDR is then NOT
# valid and a WARNING is printed). Replicates of both groups are strongly recommended.
#
# ---- Normalisation bins (what "background" means) -----------------------------------------------------------------
# TMM factors are computed from `--bin`-wide bins counted in ALL libraries (target + control) and applied as edgeR
# norm.factors to the window counts; lib.size is ALWAYS the full csaw library total (counts$totals), never the column sum
# of retained bins or windows. TMM assumes most bins are not differentially abundant between any two libraries. For
# target-vs-IgG this is false for the ENRICHED bins by construction, and the target libraries differ from IgG in the share of
# their reads that fall in enriched domains (broad marks such as H3K27me3 can cover a large fraction of bins).
#   all         = every bin that passed csaw's count filter (legacy round-9 behaviour: csaw::normFactors(bins)). The
#                 trimmed mean of M-values then sits between background and enriched bins, so the target factor is
#                 pulled toward the enrichment: effective target library size is overestimated, logFC and peak numbers
#                 are biased DOWN (conservative).
#   background  = bins whose average abundance over ALL libraries is not elevated:  aveLogCPM(bin) <= median(aveLogCPM(all
#                 bins)) + bg-max-enrich (default 1 log2 unit, i.e. < 2-fold above the typical bin). The median is the same
#                 global-background estimate csaw::filterWindowsGlobal uses, and the cut-off is applied to the mean of
#                 target AND control libraries, so selection is symmetric in the groups and independent of the contrast
#                 (under the null A and M are uncorrelated, selecting on A does not bias M). This is the composition-bias
#                 reference used for the peak call: it aligns the NON-enriched genome of the target and IgG libraries and
#                 lets enriched regions show their enrichment. Weakly enriched bins (< 2-fold) remain in the reference.
# Weighted TMM is OFF (doWeighting = FALSE) in both settings, as in csaw::normFactors(weighted = FALSE).
# --norm-bins picks which setting drives <name>.peaks.bed etc.; the other one is evaluated and reported in
# <name>.norm_sensitivity.tsv. Expected direction: background >= all in peak number and logFC; large differences mean the
# peak set depends on how composition is handled and should be reported as such.
#
# ---- lfc decision (inferential meaning) ---------------------------------------------------------------------------------
# Decision: the post-hoc filter is KEPT (round-9 behaviour): a region is a peak when FDR <= --fdr (test of H0: log2 FC = 0,
# BH over regions) AND the region's estimated log2 FC > --lfc. Which window's estimate: --lfc-from=combined (default, the legacy
# rule: GG_EED_CUTANDRUN / GG_EPICYPHER_CBX2 callpeaks_csaw.R filter the exported BED on the rep.logFC column of combineTests;
# their getBestTest logFC only feeds the "Enriched" log line; verified 2026-10-23, identical peak sets) uses the rep.logFC of
# csaw::combineTests (the window that represents the Simes combination); --lfc-from=best uses csaw::getBestTest rep.logFC, the
# logFC of the window with the smallest (weighted-Holm) p-value in the region. Both columns are always written (rep_logFC, best_logFC). Region FDR is always the COMBINED (Simes) FDR. The best
# window is selected on the same p-values it is then filtered on, so its logFC is more winner's-curse-biased than the combined
# one; the combined estimate is the legacy definition and the default, best is one flag away. Meaning: the FDR controls the
# proportion of false "Target differs from IgG" calls among the tested regions; the lfc cut is a DESCRIPTIVE effect-size
# filter applied after selection. It does NOT control the FDR for the claim "enrichment exceeds 2^lfc", and the logFC of
# the representative window is a selected (winner's-curse-biased) estimate. A true test against the threshold would use
# edgeR::glmTreat(lfc=) inside the merge; it was not adopted because it changes every round-9 peak set, is much more
# conservative with 2-3 replicates, and could not be validated here. <name>.norm_sensitivity.tsv reports n_fdr (FDR only)
# and n_peaks (FDR + lfc) so the effect of the filter is visible. Do not write "peaks are significantly >2-fold enriched".
#
# ---- Disjoint regions (--disjoint-regions) ------------------------------------------------------------------------------
# csaw::mergeWindows sorts windows by (chr, start, end) and starts a new region when start[i] - end[i-1] - 1 > tol, i.e. it chains
# on the end of the PREVIOUS window in that order, not on the running maximum end. With two window widths (150 and 2000 bp) a wide
# window is followed by narrow windows whose ends are far smaller, so a later region can start inside the span of an earlier one:
# the merged regions overlap (a region is then counted, ranked and tested more than once; reads/peaks in the overlap are shared).
# The log line "csaw merge: ... overlap" gives the number of overlapping regions/peaks for the run; it is also written to
# <name>.norm_sensitivity.tsv. Default (false) keeps csaw's regions (legacy-compatible). With true:
#   1. all windows kept at both scales are reduced with GenomicRanges::reduce(min.gapwidth = merge_bp + 1) (gap <= tol merges, the
#      csaw criterion, but chained on the running maximum end): connected components of the union, no overlaps by construction;
#   2. max.width finite: a region wider than --max-width is cut exactly as csaw's split_cluster does (ceiling(width/max.width)
#      equal tiles, windows by midpoint) and each tile is trimmed to its own tile interval (disjoint, <= max.width);
#   3. every window belongs to the region that contains its midpoint (= the region containing the whole window, except split tiles);
#   4. csaw::combineTests (Simes, region FDR by BH over the NEW regions) and csaw::getBestTest are recomputed over the concatenated
#      small + large tables with the weights mergeWindowsList would assign: 1 / (number of windows of the same scale in the region)
#      (equiweight), so each scale carries the same total weight per region.
# Expected effect: fewer, wider regions; region statistics and BH ranks change; peak numbers are NOT expected to be the same as with
# csaw's merge. The window-level fits, TMM factors and the abundance filters are unaffected. No FDR claim is added or lost: both
# definitions combine p-values of windows that are independent of the p-value (abundance filter), but only the disjoint definition has
# regions that are non-overlapping tests.
#
# ---- Other fixes relative to the legacy script -------------------------------------------------------------------------------
# readParam max.frag = --max-frag (legacy: csaw default 500 regardless of bowtie2 -X); mergeResultsList gets max.width
# from --max-width; the peak rule is selectable (--lfc-from; default combined = legacy); BED score is -10 log10(FDR) (CBX2 version), logFC in column 7 as a real number (legacy: integer-truncated);
# no reticulate/config.py; no install.packages; thresholds only from the command line; RDS of the (large) window counts is no
# longer written (results RDS only); chromosome names are kept as in the BAM header (no style conversion).

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
as_bool <- function(x, key) {
  u <- toupper(trimws(x))
  if (u %in% c("TRUE", "T", "1", "YES")) TRUE
  else if (u %in% c("FALSE", "F", "0", "NO")) FALSE
  else die(paste0("--", key, " must be TRUE or FALSE, got '", x, "'"), 2L)
}
check_bams <- function(bams) {
  miss <- bams[!file.exists(bams)]
  if (length(miss)) die(paste0("BAM file(s) not found: ", paste(miss, collapse = ", ")))
  has_idx <- file.exists(paste0(bams, ".bai")) | file.exists(sub("\\.bam$", ".bai", bams))
  if (!all(has_idx)) die(paste0("BAM index (.bai) missing for: ", paste(bams[!has_idx], collapse = ", ")))
}
save_rds_atomic <- function(obj, path) {
  tmp <- paste0(path, ".", Sys.getpid(), ".tmp")
  saveRDS(obj, tmp)
  if (!file.rename(tmp, path)) die(paste0("could not write ", path))
}

suppressPackageStartupMessages({
  library(csaw)
  library(edgeR)
  library(GenomicRanges)
  library(BiocParallel)
})

# >>> CSAW_CORE (identical in callpeaks_csaw.R and csaw_cutoff_scan.R)
make_bp <- function(cpus) {
  if (cpus > 1) BiocParallel::MulticoreParam(workers = as.integer(cpus), progressbar = FALSE) else BiocParallel::SerialParam()
}
# count small windows, large windows and background bins in all libraries. count_filter = minimum summed count over all
# libraries for a WINDOW to be kept (csaw default 10); the background bins always use the csaw default 10, so the TMM
# reference does not move when count_filter is changed.
count_all <- function(bams, ids, pe, mapq, max_frag, dedup, window, spacing, large_window, large_spacing, bin, bp, cc_bams,
                      count_filter = 10) {
  param <- csaw::readParam(minq = mapq, pe = pe, max.frag = max_frag, dedup = dedup)
  ext <- NA
  if (pe == "none") {
    cat("Single-end: estimating fragment length by cross-correlation (first target + first control)\n")
    cc <- csaw::correlateReads(cc_bams, max.dist = 500,
                               param = csaw::readParam(minq = mapq, pe = "none", dedup = TRUE))
    ext <- csaw::maximizeCcf(cc)
    cat("  estimated fragment length:", ext, "bp\n")
  }
  cat("Counting small windows (", window, " bp, spacing ", spacing, ", count filter ", count_filter, ")\n", sep = "")
  small <- csaw::windowCounts(bams, spacing = spacing, width = window, ext = ext, filter = count_filter, param = param,
                              BPPARAM = bp)
  cat("Counting large windows (", large_window, " bp, spacing ", large_spacing, ", count filter ", count_filter, ")\n", sep = "")
  large <- csaw::windowCounts(bams, spacing = large_spacing, width = large_window, ext = ext, filter = count_filter,
                              param = param, BPPARAM = bp)
  cat("Counting background bins (", bin, " bp)\n", sep = "")
  bins <- if (pe == "none") {
    csaw::windowCounts(bams, bin = TRUE, width = bin, ext = ext, filter = 10, param = param, BPPARAM = bp)
  } else {
    csaw::windowCounts(bams, bin = TRUE, width = bin, filter = 10, param = param, BPPARAM = bp)
  }
  colnames(small) <- ids; colnames(large) <- ids; colnames(bins) <- ids
  list(small = small, large = large, bins = bins)
}
# TMM factors from bins. setting = "all" (every counted bin) or "background" (bins with aveLogCPM over ALL libraries
# <= median + bg_max). lib.size is the FULL library total in both cases (never the column sums of the subset).
tmm_factors <- function(bins, setting, bg_max) {
  m <- SummarizedExperiment::assay(bins)
  totals <- bins$totals
  ab <- edgeR::aveLogCPM(m, lib.size = totals)
  ref <- stats::median(ab)
  keep <- if (setting == "background") ab <= ref + bg_max else rep(TRUE, nrow(m))
  if (sum(keep) < 1000L) {
    die(sprintf("only %d bins qualify as '%s' normalisation bins (need >= 1000); check bin size / depth", sum(keep), setting))
  }
  nf <- edgeR::calcNormFactors(m[keep, , drop = FALSE], lib.size = totals, method = "TMM", doWeighting = FALSE)
  nf <- as.numeric(nf)
  if (any(!is.finite(nf)) || any(nf <= 0)) die(paste0("TMM (", setting, ") gave non-finite factors"))
  list(nf = nf, n_bins = sum(keep), n_total = nrow(m), ref_logCPM = ref)
}
# edgeR QL fit + test of Target - Control on a window-count matrix; effective library size = totals x nf
fit_test <- function(counts, totals, nf, design, has_reps, return_y = FALSE) {
  y <- edgeR::DGEList(counts = counts, lib.size = totals, norm.factors = nf)
  if (has_reps) {
    y <- edgeR::estimateDisp(y, design)
    fit <- edgeR::glmQLFit(y, design)
    res <- edgeR::glmQLFTest(fit, contrast = c(1, -1))
  } else {
    y$common.dispersion <- 0.05
    fit <- edgeR::glmFit(y, design)
    res <- edgeR::glmLRT(fit, contrast = c(1, -1))
  }
  list(table = res$table, y = if (return_y) y else NULL)
}
# csaw's own merge (mergeResultsList): windows of both scales are sorted by (chr, start, end) and chained on the END OF THE
# PREVIOUS WINDOW IN THAT ORDER (merge_windows.cpp: start[i] - end[i-1] - 1 > tol opens a new region), so a 2000 bp window
# followed by 150 bp windows can open a region that lies inside the span of an earlier one: the merged regions can OVERLAP.
merge_call <- function(gr_list, tab_list, merge_bp, max_width) {
  merge_args <- if (is.finite(max_width)) list(max.width = max_width) else list()
  csaw::mergeResultsList(gr_list, tab.list = tab_list, equiweight = TRUE, tol = merge_bp, merge.args = merge_args)
}
# DISJOINT regions (--disjoint-regions=true). Same statistics as merge_call, other region definition:
#  1. all windows of all scales are reduced with GenomicRanges::reduce(min.gapwidth = tol + 1): two windows are in one region
#     when the gap between them is <= tol bp (csaw's rule gap <= tol) but chaining uses the running MAXIMUM end, so regions
#     can no longer overlap (connected components of the union of windows);
#  2. if max_width is finite, a region wider than max_width is cut as csaw does (mult = ceiling(width / max_width) equal
#     tiles, windows assigned to a tile by their midpoint) and each tile is trimmed to its own tile interval, so tiles stay
#     disjoint and <= max_width;
#  3. every window is assigned to exactly one region by its midpoint (floor((start + end) / 2));
#  4. combineTests (Simes) and getBestTest (weighted Holm) are recomputed over the concatenated per-scale tables with the same
#     equiweight weights as mergeWindowsList: weight = 1 / (number of windows of THAT scale in the region), recomputed for the
#     new regions; region FDR = BH over the disjoint regions. Returns the same fields as mergeResultsList (regions, combined, best).
merge_disjoint <- function(gr_list, tab_list, merge_bp, max_width) {
  if (length(gr_list) != length(tab_list)) die("internal: ranges and tables differ in length")
  n_each <- vapply(gr_list, length, FUN.VALUE = 0L)
  if (!identical(n_each, vapply(tab_list, nrow, FUN.VALUE = 0L))) die("internal: ranges and result tables differ in rows")
  if (merge_bp < 0) die("--disjoint-regions=true needs --merge-bp >= 0", 2L)
  all_gr <- do.call(c, lapply(gr_list, function(g) { S4Vectors::mcols(g) <- NULL; g }))
  origin <- rep(seq_along(gr_list), n_each)
  tol <- as.integer(merge_bp)
  red <- reduce(all_gr, min.gapwidth = tol + 1L, ignore.strand = TRUE)
  mid <- as.integer(floor((start(all_gr) + end(all_gr)) / 2))
  mid_gr <- GRanges(seqnames(all_gr), IRanges(mid, mid), seqinfo = seqinfo(all_gr))
  rid <- findOverlaps(mid_gr, red, select = "first")
  if (anyNA(rid)) die("internal: window midpoint outside every disjoint region")
  r_start <- start(red); r_width <- width(red)
  mult <- if (is.finite(max_width)) pmax(1, ceiling(r_width / max_width)) else rep(1, length(red))
  if (all(mult == 1)) {
    ids <- as.integer(rid)
    reg <- red
  } else {
    # same partition as csaw split_cluster: tile k of a region of width W cut into m tiles = [s + ceil(k W / m), s + ceil((k + 1) W / m) - 1]
    m_w <- mult[rid]
    k <- floor(((start(all_gr) + end(all_gr)) - 2 * r_start[rid]) * m_w / (2 * r_width[rid]))   # csaw: int(((s+e)/2 - region start) / subwidth)
    k <- pmin(pmax(k, 0), m_w - 1)
    key <- as.numeric(rid) * (max(mult) + 1) + k
    uk <- sort(unique(key))
    ids <- match(key, uk)
    first <- match(uk, key)                       # one window per tile to read region and tile index
    t_rid <- rid[first]; t_k <- k[first]; t_m <- mult[t_rid]; t_W <- r_width[t_rid]; t_s <- r_start[t_rid]
    lo <- t_s + ceiling(t_k * t_W / t_m)
    hi <- t_s + ceiling((t_k + 1) * t_W / t_m) - 1
    w_start <- as.numeric(tapply(start(all_gr), ids, min))
    w_end <- as.numeric(tapply(end(all_gr), ids, max))
    new_start <- as.integer(pmax(w_start, lo)); new_end <- as.integer(pmin(w_end, hi))
    if (any(new_start > new_end)) die("internal: empty tile after max.width split")
    reg <- GRanges(seqnames(red)[t_rid], IRanges(new_start, new_end), seqinfo = seqinfo(red))
  }
  weights <- numeric(length(ids))
  for (i in seq_along(gr_list)) {
    sel <- which(origin == i)
    cur <- ids[sel]
    weights[sel] <- (1 / tabulate(cur, nbins = length(reg)))[cur]
  }
  tab <- do.call(rbind, tab_list)
  combined <- csaw::combineTests(ids, tab, weights = weights)
  best <- csaw::getBestTest(ids, tab, weights = weights)
  want <- as.character(seq_along(reg))
  io <- match(want, rownames(combined)); ib <- match(want, rownames(best))
  if (anyNA(io) || anyNA(ib)) die("internal: combined statistics do not cover every disjoint region")
  list(regions = reg, combined = combined[io, , drop = FALSE], best = best[ib, , drop = FALSE])
}
merge_any <- function(gr_list, tab_list, merge_bp, max_width, disjoint) {
  if (disjoint) merge_disjoint(gr_list, tab_list, merge_bp, max_width) else merge_call(gr_list, tab_list, merge_bp, max_width)
}
# number of ranges that overlap at least one OTHER range of the same set (1 bp overlap counts; adjacency does not)
count_overlapping <- function(gr) {
  if (length(gr) < 2L) return(0L)
  ov <- findOverlaps(gr, gr, ignore.strand = TRUE)
  length(unique(S4Vectors::queryHits(ov)[S4Vectors::queryHits(ov) != S4Vectors::subjectHits(ov)]))
}
df_to_gr <- function(d) GRanges(seqnames = d$chr, ranges = IRanges(start = d$start, end = d$end))
# one row per merged region. is_fdr = combined (Simes) region FDR <= fdr. is_peak = is_fdr AND the effect-size filter:
#   lfc_from = "best":     best-window rep.logFC (getBestTest, lowest weighted p in the region) > lfc
#   lfc_from = "combined": rep.logFC of combineTests (representative window of the Simes combination) > lfc   (legacy rule, default)
summarise_merged <- function(merged, fdr, lfc, lfc_from) {
  if (!(lfc_from %in% c("best", "combined"))) die("lfc_from must be 'best' or 'combined'", 2L)
  tabcom <- merged$combined; tabbest <- merged$best; reg <- merged$regions
  getc <- function(tab, nm) if (nm %in% colnames(tab)) tab[[nm]] else rep(NA, nrow(tab))
  rep_lfc <- if ("rep.logFC" %in% colnames(tabcom)) tabcom[["rep.logFC"]] else tabbest[["rep.logFC"]]
  best_lfc <- getc(tabbest, "rep.logFC")
  fdr_v <- tabcom[["FDR"]]
  is_fdr <- !is.na(fdr_v) & fdr_v <= fdr
  used_lfc <- if (lfc_from == "best") best_lfc else rep_lfc
  data.frame(chr = as.character(seqnames(reg)), start = start(reg), end = end(reg), width = width(reg),
             num_tests = getc(tabcom, "num.tests"), num_up = getc(tabcom, "num.up.logFC"),
             num_down = getc(tabcom, "num.down.logFC"), PValue = tabcom[["PValue"]], FDR = fdr_v,
             direction = as.character(getc(tabcom, "direction")), rep_logFC = rep_lfc,
             best_logFC = best_lfc, best_PValue = getc(tabbest, "PValue"),
             best_FDR = getc(tabbest, "FDR"), is_fdr = is_fdr,
             is_peak = is_fdr & !is.na(used_lfc) & used_lfc > lfc)
}
# <<< CSAW_CORE

# ---- arguments ----------------------------------------------------------------------------------------------------------
allowed <- c("peakset", "target-bams", "control-bams", "target-ids", "control-ids", "window", "spacing", "bin",
             "small-filt", "large-filt", "merge-bp", "max-width", "fdr", "lfc", "mapq", "max-frag", "norm-bins", "dedup",
             "cpus", "out-prefix", "large-window", "large-spacing", "pe", "bg-max-enrich", "count-filter", "lfc-from",
             "disjoint-regions")
args <- parse_args(commandArgs(trailingOnly = TRUE), allowed)
peakset    <- get_req(args, "peakset")
t_bams     <- split_list(get_req(args, "target-bams"))
c_bams     <- split_list(get_req(args, "control-bams"))
t_ids      <- split_list(get_req(args, "target-ids"))
c_ids      <- split_list(get_req(args, "control-ids"))
window     <- as_num(get_req(args, "window"), "window")
spacing    <- as_num(get_req(args, "spacing"), "spacing")
bin        <- as_num(get_req(args, "bin"), "bin")
small_filt <- as_num(get_req(args, "small-filt"), "small-filt")
large_filt <- as_num(get_req(args, "large-filt"), "large-filt")
merge_bp   <- as_num(get_req(args, "merge-bp"), "merge-bp")
mw_raw     <- toupper(get_opt(args, "max-width", "NA"))
max_width  <- if (mw_raw %in% c("NA", "NONE", "INF", "0")) Inf else as_num(mw_raw, "max-width")
fdr        <- as_num(get_req(args, "fdr"), "fdr")
lfc        <- as_num(get_req(args, "lfc"), "lfc")
mapq       <- as_num(get_req(args, "mapq"), "mapq")
max_frag   <- as_num(get_req(args, "max-frag"), "max-frag")
norm_bins  <- get_req(args, "norm-bins")
dedup      <- as_bool(get_req(args, "dedup"), "dedup")
cpus       <- as_num(get_req(args, "cpus"), "cpus")
out_prefix <- get_req(args, "out-prefix")
large_window  <- as_num(get_opt(args, "large-window", "2000"), "large-window")
large_spacing <- as_num(get_opt(args, "large-spacing", "500"), "large-spacing")
pe         <- get_opt(args, "pe", "both")
bg_max     <- as_num(get_opt(args, "bg-max-enrich", "1"), "bg-max-enrich")
count_filter <- as_num(get_opt(args, "count-filter", "10"), "count-filter")
lfc_from   <- tolower(get_opt(args, "lfc-from", "combined"))
disjoint   <- as_bool(get_opt(args, "disjoint-regions", "false"), "disjoint-regions")

if (!(norm_bins %in% c("background", "all"))) die("--norm-bins must be 'background' or 'all'", 2L)
if (!(pe %in% c("both", "none"))) die("--pe must be 'both' or 'none'", 2L)
if (!(lfc_from %in% c("best", "combined"))) die("--lfc-from must be 'best' or 'combined'", 2L)
if (count_filter < 0) die("--count-filter must be >= 0", 2L)
if (disjoint && merge_bp < 0) die("--disjoint-regions=true needs --merge-bp >= 0", 2L)
if (fdr <= 0 || fdr > 1) die("--fdr must be in (0, 1]", 2L)
if (lfc < 0) die("--lfc must be >= 0", 2L)
if (window < 1 || spacing < 1 || bin < 1000 || large_window <= window || large_spacing < 1) {
  die("need window >= 1, spacing >= 1, bin >= 1000, large-window > window", 2L)
}
if (cpus < 1) die("--cpus must be >= 1", 2L)
if (length(t_bams) < 1L || length(c_bams) < 1L) die("need at least one target BAM and one control BAM")
if (length(t_bams) != length(t_ids)) die("--target-bams and --target-ids differ in length")
if (length(c_bams) != length(c_ids)) die("--control-bams and --control-ids differ in length")
bams <- c(t_bams, c_bams); ids <- c(t_ids, c_ids)
if (anyDuplicated(ids)) die(paste0("duplicated sample ids: ", paste(ids[duplicated(ids)], collapse = ", ")))
check_bams(bams)

n_t <- length(t_bams); n_c <- length(c_bams)
group <- factor(c(rep("Target", n_t), rep("Control", n_c)), levels = c("Target", "Control"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
has_reps <- (nrow(design) - ncol(design)) > 0
if (!has_reps) {
  cat("WARNING: no replicates (", n_t, " target, ", n_c, " control): dispersion fixed at 0.05, LRT; FDR is NOT valid\n",
      sep = "", file = stderr())
}
cat("=== callpeaks_csaw.R: ", peakset, " ===\n", sep = "")
cat("targets: ", paste(t_ids, collapse = ", "), "\ncontrols: ", paste(c_ids, collapse = ", "), "\n", sep = "")
cat(sprintf("window=%g/%g large=%g/%g bin=%g small.filt=%g large.filt=%g merge.bp=%g max.width=%s fdr=%g lfc=%g mapq=%g max.frag=%g dedup=%s pe=%s norm.bins=%s bg.max.enrich=%g count.filter=%g lfc.from=%s disjoint.regions=%s\n",
            window, spacing, large_window, large_spacing, bin, small_filt, large_filt, merge_bp, max_width, fdr, lfc, mapq,
            max_frag, dedup, pe, norm_bins, bg_max, count_filter, lfc_from, disjoint))

# ---- counting, filtering --------------------------------------------------------------------------------------------------
bp <- make_bp(cpus)
cnt <- count_all(bams, ids, pe, mapq, max_frag, dedup, window, spacing, large_window, large_spacing, bin, bp,
                 cc_bams = c(t_bams[1], c_bams[1]), count_filter = count_filter)
totals <- cnt$bins$totals
names(totals) <- ids
cat("Library sizes (csaw totals): ", paste(ids, totals, sep = "=", collapse = ", "), "\n", sep = "")

filt_small <- csaw::filterWindowsGlobal(cnt$small, cnt$bins)$filter
filt_large <- csaw::filterWindowsGlobal(cnt$large, cnt$bins)$filter
keep_small <- filt_small > small_filt
keep_large <- filt_large > large_filt
cat("Small windows retained: ", sum(keep_small), " / ", length(keep_small), "\n", sep = "")
cat("Large windows retained: ", sum(keep_large), " / ", length(keep_large), "\n", sep = "")
if (sum(keep_small) == 0L || sum(keep_large) == 0L) die("no windows pass the abundance filter at one scale; lower --small-filt/--large-filt")
small_f <- cnt$small[keep_small, ]
large_f <- cnt$large[keep_large, ]
bins <- cnt$bins
rm(cnt); invisible(gc())

# ---- normalisation factors under both settings -----------------------------------------------------------------------------
tmm <- list(background = tmm_factors(bins, "background", bg_max), all = tmm_factors(bins, "all", bg_max))
for (k in names(tmm)) {
  cat(sprintf("TMM factors [%s, %d/%d bins]: %s\n", k, tmm[[k]]$n_bins, tmm[[k]]$n_total,
              paste(ids, signif(tmm[[k]]$nf, 4), sep = "=", collapse = ", ")))
}
other <- setdiff(c("background", "all"), norm_bins)

# ---- model + merge for one normalisation setting ---------------------------------------------------------------------------
analyse <- function(nf, keep_y) {
  fs <- fit_test(SummarizedExperiment::assay(small_f), small_f$totals, nf, design, has_reps, return_y = keep_y)
  fl <- fit_test(SummarizedExperiment::assay(large_f), large_f$totals, nf, design, has_reps, return_y = keep_y)
  grs <- list(SummarizedExperiment::rowRanges(small_f), SummarizedExperiment::rowRanges(large_f))
  tabs <- list(fs$table, fl$table)
  # csaw's own merge is always computed: it is the default result and the overlap diagnostic (see "Disjoint regions")
  std <- merge_call(grs, tabs, merge_bp, max_width)
  res_std <- summarise_merged(std, fdr, lfc, lfc_from)
  dgn <- list(n_regions = nrow(res_std), n_ov_regions = count_overlapping(std$regions), n_peaks = sum(res_std$is_peak),
               n_ov_peaks = count_overlapping(df_to_gr(res_std[res_std$is_peak, , drop = FALSE])))
  if (disjoint) {
    res <- summarise_merged(merge_disjoint(grs, tabs, merge_bp, max_width), fdr, lfc, lfc_from)
    if (count_overlapping(df_to_gr(res)) != 0L) die("internal: disjoint regions overlap")
  } else {
    res <- res_std
  }
  list(res = res, diag = dgn, y_small = fs$y, y_large = fl$y)
}
cat("Fitting edgeR models with norm-bins = ", norm_bins, " (primary)\n", sep = "")
prim <- analyse(tmm[[norm_bins]]$nf, keep_y = TRUE)
cat("Fitting edgeR models with norm-bins = ", other, " (sensitivity)\n", sep = "")
alt <- analyse(tmm[[other]]$nf, keep_y = FALSE)
runs <- setNames(list(prim$res, alt$res), c(norm_bins, other))
diags <- setNames(list(prim$diag, alt$diag), c(norm_bins, other))

res <- prim$res
cat("Merged regions: ", nrow(res), "; FDR <= ", fdr, ": ", sum(res$is_fdr), "; peaks (FDR and ", lfc_from, " logFC > ", lfc, "): ",
    sum(res$is_peak), "\n", sep = "")
dg <- prim$diag
cat(sprintf("csaw merge: %d regions, %d (%.1f%%) overlap at least one other region; %d peaks, %d (%.1f%%) overlap another peak%s\n",
            dg$n_regions, dg$n_ov_regions, 100 * dg$n_ov_regions / max(1, dg$n_regions), dg$n_peaks, dg$n_ov_peaks,
            100 * dg$n_ov_peaks / max(1, dg$n_peaks),
            if (disjoint) sprintf("; reported regions are DISJOINT: %d regions, %d peaks", nrow(res), sum(res$is_peak)) else
              " (--disjoint-regions=false: reported regions are csaw's, overlaps kept)"))

# ---- BED ------------------------------------------------------------------------------------------------------------------------
pk <- res[res$is_peak, , drop = FALSE]
bed <- data.frame(chr = pk$chr, start = as.integer(pk$start - 1L), end = as.integer(pk$end),
                  name = if (nrow(pk)) paste0("peak_", seq_len(nrow(pk))) else character(0),
                  score = as.integer(pmin(1000, round(-10 * log10(pmax(pk$FDR, 1e-300))))),
                  strand = rep(".", nrow(pk)), logFC = round(if (lfc_from == "best") pk$best_logFC else pk$rep_logFC, 4))
write.table(bed, paste0(out_prefix, ".peaks.bed"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ---- all regions, gzip ------------------------------------------------------------------------------------------------------------
gz <- gzfile(paste0(out_prefix, ".csaw_results.tsv.gz"), "w")
write.table(res, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)

# ---- MDS (replicate QC) -------------------------------------------------------------------------------------------------------------
cols <- ifelse(group == "Target", "#1f77b4", "#d62728")
svg(paste0(out_prefix, ".mds.svg"), width = 14, height = 7)
par(mfrow = c(2, 4), mar = c(5, 4, 3, 2), cex.main = 0.8)
for (sc in c("small", "large")) {
  y <- if (sc == "small") prim$y_small else prim$y_large
  lc <- edgeR::cpm(y, log = TRUE)
  for (top in c(100, 500, 1000, 5000)) {
    ok <- tryCatch({
      limma::plotMDS(lc, top = top, labels = ids, col = cols, main = sprintf("Top %d %s windows (blue Target, red Control)", top, sc))
      TRUE
    }, error = function(e) {
      plot.new(); title(main = sprintf("MDS top %d %s: not drawn", top, sc))
      mtext(conditionMessage(e), cex = 0.5)
      FALSE
    })
  }
}
invisible(dev.off())

# ---- sensitivity table ------------------------------------------------------------------------------------------------------------------
as_gr <- df_to_gr
n_overlap <- function(a, b) {
  pa <- a[a$is_peak, , drop = FALSE]; pb <- b[b$is_peak, , drop = FALSE]
  if (nrow(pa) == 0L || nrow(pb) == 0L) return(0L)
  sum(suppressWarnings(IRanges::overlapsAny(as_gr(pa), as_gr(pb))))
}
sens <- do.call(rbind, lapply(names(runs), function(k) {
  r <- runs[[k]]; o <- runs[[setdiff(names(runs), k)]]
  data.frame(peakset = peakset, norm_bins = k, is_primary = (k == norm_bins), sample_id = ids, group = as.character(group),
             lib_size = as.numeric(totals), n_bins_used = tmm[[k]]$n_bins, n_bins_total = tmm[[k]]$n_total,
             bg_ref_logCPM = tmm[[k]]$ref_logCPM, bg_max_enrich = bg_max, tmm_factor = tmm[[k]]$nf,
             log2_tmm_factor = log2(tmm[[k]]$nf), n_small_kept = sum(keep_small), n_large_kept = sum(keep_large),
             n_regions = nrow(r), n_fdr = sum(r$is_fdr), n_peaks = sum(r$is_peak),
             n_peaks_overlapping_other_setting = n_overlap(r, o), fdr = fdr, lfc = lfc,
             lfc_from = lfc_from, disjoint_regions = disjoint, count_filter = count_filter,
             n_regions_csaw_merge = diags[[k]]$n_regions, n_overlapping_regions_csaw_merge = diags[[k]]$n_ov_regions,
             n_peaks_csaw_merge = diags[[k]]$n_peaks, n_overlapping_peaks_csaw_merge = diags[[k]]$n_ov_peaks)
}))
write.table(sens, paste0(out_prefix, ".norm_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- RDS (results only; window count matrices are not stored) -------------------------------------------------------------------------------
save_rds_atomic(list(
  peakset = peakset, primary_norm_bins = norm_bins, ids = ids, group = as.character(group), lib_size = totals,
  params = list(window = window, spacing = spacing, large_window = large_window, large_spacing = large_spacing, bin = bin,
                small_filt = small_filt, large_filt = large_filt, merge_bp = merge_bp, max_width = max_width, fdr = fdr,
                lfc = lfc, mapq = mapq, max_frag = max_frag, dedup = dedup, pe = pe, bg_max_enrich = bg_max,
                count_filter = count_filter, lfc_from = lfc_from, disjoint_regions = disjoint),
  tmm = lapply(tmm, function(x) x[c("nf", "n_bins", "n_total", "ref_logCPM")]),
  results = runs, session = utils::sessionInfo()
), paste0(out_prefix, ".rds"))
cat("Done: ", out_prefix, ".peaks.bed (", nrow(bed), " peaks)\n", sep = "")
