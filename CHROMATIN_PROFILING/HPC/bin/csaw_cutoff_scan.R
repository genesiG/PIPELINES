#!/usr/bin/env Rscript
# csaw_cutoff_scan.R -- POST-HOC DIAGNOSTIC: number of csaw peaks over a grid of window-filter cut-offs. Never feeds thresholds.
#
# CLI contract (called by the optional CSAW_CUTOFF_SCAN module): the same arguments as callpeaks_csaw.R, plus
#   --small-grid=<comma list>   candidate --small-filt values (log2 enrichment over background), e.g. 2,2.5,3,...
#   --large-grid=<comma list>   candidate --large-filt values
# i.e. --peakset --target-bams --control-bams --target-ids --control-ids --window --spacing --bin --small-filt --large-filt
#      --merge-bp --max-width --fdr --lfc --mapq --max-frag --norm-bins --dedup --cpus --out-prefix
#      [--large-window=2000 --large-spacing=500 --pe=both --bg-max-enrich=1 --count-filter=10 --lfc-from=combined
#       --disjoint-regions=false]
# --count-filter, --lfc-from and --disjoint-regions mean exactly what they mean in callpeaks_csaw.R (they define the peak set that
# is being scanned, so they must match the peak-calling run); --disjoint-regions=true replaces csaw's overlapping merge by the
# disjoint regions documented there. The csaw overlap diagnostic is NOT computed here (see callpeaks_csaw.R).
# The current --small-filt / --large-filt are ALWAYS added to the grids (row flagged is_current = TRUE) so the current peak set
# can be located in the scan. See callpeaks_csaw.R for the meaning of every shared argument.
# Output:
#   <name>.cutoff_scan.tsv   one row per (c_small, c_large): peakset, label ("POST-HOC DIAGNOSTIC ..."), norm_bins, c_small,
#                            c_large, is_current, n_small_kept, n_large_kept, n_regions, n_fdr, n_peaks, width_min/q25/median/
#                            q75/max (of peak regions), fdr, lfc, merge_bp, max_width, lfc_from, disjoint_regions, count_filter
# Exit codes: 0 ok; 2 usage error; 3 data problem.
#
# ---- Interpretation (read before using the table) ---------------------------------------------------------------------
# This is a descriptive sensitivity scan, NOT an optimiser. Picking the grid cell with the most peaks chooses the thresholds
# with the same data that are then tested: the number of peaks is not an estimate of anything and no error rate holds for the
# selected cell. Report it as "peak number as a function of the abundance filter", and keep the thresholds in peaksets.csv fixed
# beforehand. The `label` column says so on every row; nothing downstream reads this table.
#
# ---- What differs from the legacy optimize_csaw_cutoffs.R (fixes) ---------------------------------------------------------
# 1. DGEList library sizes are the full csaw totals (bins$totals) -- legacy used the column sums of the RETAINED windows, which
#    makes every library's size depend on the cut-off and silently re-normalises composition.
# 2. Normalisation is exactly callpeaks_csaw.R's TMM-on-bins (selected by --norm-bins, same background-bin definition,
#    lib.size = totals). Legacy used 1/SF of ChIPseqSpikeInFree/spike-in tables with the IgG factors mixed in (wrong direction
#    for edgeR, which wants divisors of the library size, and IgG must never be scaled by another run's factors). No external
#    factors are accepted here: factors from norm_chipseqspikeinfree.R / norm_spikein.py apply to bigWigs, not to this model.
# 3. Window filter is `score > cut-off` as in callpeaks_csaw.R (legacy scan used >=); peak rule is FDR <= --fdr AND the logFC chosen
#    by --lfc-from > --lfc (default combined = representative window of the combined test, as in the legacy scripts; `best` =
#    best-window logFC); legacy scan used FDR <; max.width is honoured.
# 4. Efficiency: windows of one scale do not depend on the other scale's cut-off, so each cut-off is fitted once per scale
#    (|small-grid| + |large-grid| edgeR fits, not 2 x |small-grid| x |large-grid|); only the merge is done per grid cell.
#    Dispersion/QL fits on a retained set are identical to what callpeaks_csaw.R does for that cut-off.

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
             "disjoint-regions", "small-grid", "large-grid")
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
parse_grid <- function(x, key) {
  v <- suppressWarnings(as.numeric(split_list(x)))
  if (length(v) == 0L || anyNA(v)) die(paste0("--", key, " must be a comma list of numbers"), 2L)
  v
}
small_grid <- sort(unique(c(parse_grid(get_req(args, "small-grid"), "small-grid"), small_filt)))
large_grid <- sort(unique(c(parse_grid(get_req(args, "large-grid"), "large-grid"), large_filt)))

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
  cat("WARNING: no replicates: dispersion fixed at 0.05, LRT; FDR is NOT valid\n", file = stderr())
}
cat("=== csaw_cutoff_scan.R (POST-HOC DIAGNOSTIC): ", peakset, " ===\n", sep = "")
cat(sprintf("peak definition: lfc.from=%s disjoint.regions=%s count.filter=%g fdr=%g lfc=%g merge.bp=%g\n", lfc_from, disjoint,
            count_filter, fdr, lfc, merge_bp))
cat("small grid: ", paste(small_grid, collapse = ","), "\nlarge grid: ", paste(large_grid, collapse = ","), "\n", sep = "")

# ---- counting, scores, normalisation (once) ---------------------------------------------------------------------------------
bp <- make_bp(cpus)
cnt <- count_all(bams, ids, pe, mapq, max_frag, dedup, window, spacing, large_window, large_spacing, bin, bp,
                 cc_bams = c(t_bams[1], c_bams[1]), count_filter = count_filter)
totals <- cnt$bins$totals
cat("Library sizes (csaw totals): ", paste(ids, totals, sep = "=", collapse = ", "), "\n", sep = "")
filt_small <- csaw::filterWindowsGlobal(cnt$small, cnt$bins)$filter
filt_large <- csaw::filterWindowsGlobal(cnt$large, cnt$bins)$filter
tm <- tmm_factors(cnt$bins, norm_bins, bg_max)
nf <- tm$nf
cat(sprintf("TMM factors [%s, %d/%d bins]: %s\n", norm_bins, tm$n_bins, tm$n_total,
            paste(ids, signif(nf, 4), sep = "=", collapse = ", ")))

m_small <- SummarizedExperiment::assay(cnt$small); gr_small <- SummarizedExperiment::rowRanges(cnt$small)
m_large <- SummarizedExperiment::assay(cnt$large); gr_large <- SummarizedExperiment::rowRanges(cnt$large)
tot_small <- totals; tot_large <- totals   # bins$totals: full library sizes for every cut-off
rm(cnt); invisible(gc())

# ---- one edgeR fit per cut-off and scale (library size = totals, factors = nf) ----------------------------------------------
fit_scale <- function(m, tot, filt, cuts) {
  res <- BiocParallel::bplapply(cuts, function(ct) {
    keep <- filt > ct
    n <- sum(keep)
    if (n == 0L) return(list(n = 0L, table = NULL))
    ft <- fit_test(m[keep, , drop = FALSE], tot, nf, design, has_reps)
    list(n = n, table = ft$table[, c("logFC", "PValue"), drop = FALSE])
  }, BPPARAM = bp)
  names(res) <- as.character(cuts)
  res
}
cat("Fitting ", length(small_grid), " small-window and ", length(large_grid), " large-window cut-offs\n", sep = "")
fits_small <- fit_scale(m_small, tot_small, filt_small, small_grid)
fits_large <- fit_scale(m_large, tot_large, filt_large, large_grid)

# ---- merge per grid cell ---------------------------------------------------------------------------------------------------------
cells <- expand.grid(i = seq_along(small_grid), j = seq_along(large_grid))
cat("Merging ", nrow(cells), " grid cells\n", sep = "")
rows <- BiocParallel::bplapply(seq_len(nrow(cells)), function(r) {
  cs <- small_grid[cells$i[r]]; cl <- large_grid[cells$j[r]]
  fs <- fits_small[[cells$i[r]]]; fl <- fits_large[[cells$j[r]]]
  grs <- list(); tabs <- list()
  if (fs$n > 0L) { grs <- c(grs, list(gr_small[filt_small > cs])); tabs <- c(tabs, list(fs$table)) }
  if (fl$n > 0L) { grs <- c(grs, list(gr_large[filt_large > cl])); tabs <- c(tabs, list(fl$table)) }
  out <- data.frame(c_small = cs, c_large = cl, n_small_kept = fs$n, n_large_kept = fl$n)
  if (length(grs) == 0L) {
    return(cbind(out, n_regions = 0L, n_fdr = 0L, n_peaks = 0L, width_min = NA_real_, width_q25 = NA_real_,
                 width_median = NA_real_, width_q75 = NA_real_, width_max = NA_real_))
  }
  sm <- summarise_merged(merge_any(grs, tabs, merge_bp, max_width, disjoint), fdr, lfc, lfc_from)
  w <- sm$width[sm$is_peak]
  q <- if (length(w)) as.numeric(stats::quantile(w, c(0, 0.25, 0.5, 0.75, 1))) else rep(NA_real_, 5)
  cbind(out, n_regions = nrow(sm), n_fdr = sum(sm$is_fdr), n_peaks = sum(sm$is_peak), width_min = q[1], width_q25 = q[2],
        width_median = q[3], width_q75 = q[4], width_max = q[5])
}, BPPARAM = bp)
scan <- do.call(rbind, rows)
scan <- data.frame(peakset = peakset, label = "POST-HOC DIAGNOSTIC: data-dependent threshold scan, not an estimate; do not select thresholds from it",
                   norm_bins = norm_bins, scan, is_current = (scan$c_small == small_filt & scan$c_large == large_filt),
                   fdr = fdr, lfc = lfc, merge_bp = merge_bp, max_width = max_width, lfc_from = lfc_from,
                   disjoint_regions = disjoint, count_filter = count_filter)
write.table(scan, paste0(out_prefix, ".cutoff_scan.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cur <- scan[scan$is_current, ]
cat("Current setting (", small_filt, ", ", large_filt, "): ", cur$n_peaks, " peaks; scan range of n_peaks: ", min(scan$n_peaks), " - ",
    max(scan$n_peaks), " (post-hoc diagnostic)\n", sep = "")
cat("Done: ", out_prefix, ".cutoff_scan.tsv\n", sep = "")
