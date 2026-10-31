#!/usr/bin/env Rscript
# norm_qc_plot.R -- plot step of the normalisation QC: one MA panel per pair, one reference line per normalisation method.
# Reads only tables (pairs table of norm_qc.R, the wide scale-factor table, bin counts, libsizes). The line positions are read
# from the pairs table (`line`); only the cloud (A, M per pair) is recomputed from the counts, with the same formulas as norm_qc.R.
#
# CLI contract (called by the NORM_QC_PLOT module; R style --key=value):
#   --pairs=<run>.norm_qc_pairs.tsv        pairs table from norm_qc.R (pair_id, pair_type, ref, other, method, is_primary, status, line,
#                                          n_bins, min_count, prior, M_peak, core_lo, core_hi, ...)
#   --scale-factors-all=<run>.scale_factors_all.tsv   wide table (used only for the reason text of unavailable methods)
#   --counts=<bins_<size>.counts.tsv.gz> --libsizes=<bins_<size>.libsizes.tsv>   BIN_COUNTS output (same files as NORM_QC)
#   --out-prefix=<run>
#   OPTIONAL: --panels-per-page=12 (3 columns x up to 4 rows, 3.9 x 4.4 in per panel);
#             --lines-out=<tsv>  writes page, panel, ref, other, method, y_drawn, ylim_lo, ylim_hi;
#             --theme=<file>  accepted and ignored (base graphics; kept so a generic module can pass it)
# Output: <run>.ma.svg (page 1), <run>.ma.p2.svg, <run>.ma.p3.svg, ...
#   Panel: x = A = mean log2 CPM (0.1-99.9 % quantiles), y = M = log2 CPM(other) - log2 CPM(ref) (0.5-99.5 % quantiles, extended to
#   include every drawn line + 10 %; negative = less signal in `other`). smoothScatter density as ONE raster image
#   (options(preferRaster = TRUE)), nrpoints = 100. Horizontal lines: dashed = status ok, dotted = qc_fail or imputed (factor of a
#   QC-failed or imputed library); primary lwd 2.5, others 1.5; Okabe-Ito colours (chipseqspikeinfree #D55E00, tmm #0072B2,
#   cpm #000000, spikein #009E73). Corner text = the line values (3 decimals). Black bar at the right margin = background-cloud
#   core (density peak .. density-weighted centroid), white dot = density peak. One legend per page in the bottom outer margin;
#   unavailable methods in grey; not_run methods are omitted. A pair with < 100 bins gets a text panel "too few bins (n)".
# Exit codes: 0 ok; 2 usage; 3 data problem.
# Reading the plot: a line on the dense cloud means that method keeps the bulk of the bins unchanged. This is by construction for TMM
# and CPM and is not evidence that they are right when the treatment changes the mark globally; between replicates every reliable
# method should lie on the cloud.

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

allowed <- c("pairs", "scale-factors-all", "counts", "libsizes", "out-prefix", "panels-per-page", "lines-out", "theme")
args <- parse_args(commandArgs(trailingOnly = TRUE), allowed)
pairs_f    <- get_req(args, "pairs")
wide_f     <- get_req(args, "scale-factors-all")
counts_f   <- get_req(args, "counts")
libsize_f  <- get_req(args, "libsizes")
out_prefix <- get_req(args, "out-prefix")
per_page   <- suppressWarnings(as.integer(get_opt(args, "panels-per-page", "12")))
lines_out  <- get_opt(args, "lines-out", "")
if (is.na(per_page) || per_page < 1L) die("--panels-per-page must be a positive integer", 2L)
for (f in c(pairs_f, wide_f, counts_f, libsize_f)) if (!file.exists(f)) die(paste0("file not found: ", f))

METHODS <- c("chipseqspikeinfree", "tmm", "cpm", "spikein")
MCOL <- c(chipseqspikeinfree = "#D55E00", tmm = "#0072B2", cpm = "#000000", spikein = "#009E73")
MABB <- c(chipseqspikeinfree = "CSIF", tmm = "TMM", cpm = "CPM", spikein = "SPIKE")
MLAB <- c(chipseqspikeinfree = "ChIPseqSpikeInFree", tmm = "TMM", cpm = "CPM", spikein = "Spike-in")

pairs <- read.delim(pairs_f, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
need <- c("pair_id", "pair_type", "ref", "other", "method", "is_primary", "status", "line", "n_bins", "min_count", "prior", "M_peak", "core_lo", "core_hi")
if (!all(need %in% colnames(pairs))) die(paste0("pairs table lacks columns: ", paste(setdiff(need, colnames(pairs)), collapse = ", ")))
if (!nrow(pairs)) die("pairs table is empty")
pairs$is_primary <- as.logical(pairs$is_primary)
pair_ids <- unique(pairs$pair_id)
ids <- unique(c(pairs$ref, pairs$other))

# counts: read only the blacklist flag and the run's columns
hdr <- strsplit(readLines(gzfile(counts_f), n = 1L), "\t", fixed = TRUE)[[1]]
if (!all(c("blacklisted", ids) %in% hdr)) die(paste0("counts table lacks columns: ", paste(setdiff(c("blacklisted", ids), hdr), collapse = ", ")))
cc <- ifelse(hdr %in% c("start", "end", "blacklisted", ids), NA, "NULL"); cc[hdr == "blacklisted"] <- "logical"
cn <- read.delim(counts_f, check.names = FALSE, stringsAsFactors = FALSE, colClasses = cc)
blk <- cn$blacklisted; blk[is.na(blk)] <- FALSE
cm <- as.matrix(cn[, ids, drop = FALSE])
bin_bp <- median(cn$end - cn$start + 1)
ls0 <- read.delim(libsize_f, check.names = FALSE, stringsAsFactors = FALSE)
Tt <- setNames(as.numeric(ls0$csaw_total[match(ids, ls0$sample_id)]), ids)
if (anyNA(Tt)) die("libsizes table lacks some of the pair libraries")

# unavailable reasons from the wide table
w <- read.delim(wide_f, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character", na.strings = c("NA", ""))
unavail_reason <- function(m) {
  nt <- if (paste0(m, "_note") %in% colnames(w)) w[[paste0(m, "_note")]] else character(0)
  nt <- nt[!is.na(nt) & nzchar(nt)]
  if (!length(nt)) return("no reason given")
  sub("^UNAVAILABLE:\\s*", "", nt[grepl("^UNAVAILABLE", nt)][1] %||% nt[1])
}
`%||%` <- function(a, b) if (length(a) && !is.na(a)) a else b

ma_coords <- function(c_ref, c_oth, T_ref, T_oth, prior) {   # same formulas as norm_qc.R
  mT <- mean(c(T_ref, T_oth))
  p_r <- prior * T_ref / mT; p_o <- prior * T_oth / mT
  xr <- log2((c_ref + p_r) / (T_ref + 2 * p_r) * 1e6)
  xo <- log2((c_oth + p_o) / (T_oth + 2 * p_o) * 1e6)
  list(A = (xo + xr) / 2, M = xo - xr)
}

present <- METHODS[METHODS %in% pairs$method]
primary_m <- unique(pairs$method[pairs$is_primary])
unavail_m <- present[vapply(present, function(m) all(pairs$status[pairs$method == m] == "unavailable"), logical(1))]

# ---- legend (one per page, bottom outer margin) ----------------------------------------------------------------------------------
draw_legend <- function() {
  ok_m <- setdiff(present, unavail_m)
  lab <- character(0); col <- character(0); lty <- integer(0); lwd <- numeric(0)
  for (m in ok_m) {
    lab <- c(lab, paste0(MLAB[[m]], if (m %in% primary_m) " (primary)" else ""))
    col <- c(col, MCOL[[m]]); lty <- c(lty, 2L); lwd <- c(lwd, if (m %in% primary_m) 2.5 else 1.5)
  }
  lab <- c(lab, "dotted = factor from a QC-failed or imputed library"); col <- c(col, "black"); lty <- c(lty, 3L); lwd <- c(lwd, 1.5)
  lab <- c(lab, "black bar = background-cloud core, white dot = density peak"); col <- c(col, "black"); lty <- c(lty, 1L); lwd <- c(lwd, 5)
  for (m in unavail_m) {
    lab <- c(lab, paste0(m, ": unavailable (", unavail_reason(m), ")")); col <- c(col, "grey60"); lty <- c(lty, 1L); lwd <- c(lwd, 1.5)
  }
  legend(x = grconvertX(0.5, "ndc", "user"), y = grconvertY(0.005, "ndc", "user"), xjust = 0.5, yjust = 0, xpd = NA, bty = "n",
         cex = 0.75, ncol = 2, legend = lab, col = col, lty = lty, lwd = lwd)
}

# ---- pages ------------------------------------------------------------------------------------------------------------------------
options(preferRaster = TRUE)  # smoothScatter -> image(): raster layer instead of ~16k vector rects per panel
n_pan <- length(pair_ids)
n_page <- ceiling(n_pan / per_page)
lines_log <- list()
for (pg in seq_len(n_page)) {
  idx <- ((pg - 1L) * per_page + 1L):min(n_pan, pg * per_page)
  ncol_p <- min(3L, length(idx)); nrow_p <- ceiling(length(idx) / ncol_p)
  fn <- if (pg == 1L) paste0(out_prefix, ".ma.svg") else paste0(out_prefix, ".ma.p", pg, ".svg")
  svg(fn, width = 3.9 * ncol_p + 0.4, height = 4.4 * nrow_p + 0.4 + 1.0)
  op <- par(mfrow = c(nrow_p, ncol_p), mar = c(4.5, 4, 5, 1.5), oma = c(5, 0, 0, 0))
  for (j in seq_along(idx)) {
    pid <- pair_ids[idx[j]]
    pr <- pairs[pairs$pair_id == pid, , drop = FALSE]
    r <- pr$ref[1]; o <- pr$other[1]
    ttl <- paste0(o, " vs ", r)
    typ <- if (pr$pair_type[1] == "replicate") "replicates" else "treated vs reference"
    ln <- pr[is.finite(pr$line), , drop = FALSE]
    ln <- ln[order(ln$is_primary), , drop = FALSE]  # primary drawn last
    sel <- !blk & (cm[, r] + cm[, o] >= pr$min_count[1])
    n <- sum(sel)
    sub_t <- sprintf("%g kb bins, pair sum >= %g, prior %g, n = %d", bin_bp / 1000, pr$min_count[1], pr$prior[1], n)
    corner <- function(cex = 0.6) {
      u <- par("usr"); xs <- u[1] + 0.02 * (u[2] - u[1]); ytop <- u[4] - 0.02 * (u[4] - u[3])
      pm <- ln[!duplicated(ln$method), , drop = FALSE]
      pm <- pm[match(intersect(METHODS, pm$method), pm$method), , drop = FALSE]
      for (q in seq_len(nrow(pm))) {
        tk <- sprintf("%s %+.3f", MABB[[pm$method[q]]], pm$line[q])
        text(xs, ytop, tk, adj = c(0, 1), cex = cex, col = MCOL[[pm$method[q]]], font = 2, xpd = NA)
        xs <- xs + strwidth(tk, cex = cex, font = 2) + strwidth("  ", cex = cex)
      }
    }
    if (n < 100L) {
      plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "", main = "")
      title(main = ttl, cex.main = 0.8, line = 3); mtext(typ, line = 2, cex = 0.7); mtext(sub_t, line = 1, cex = 0.6)
      text(0.5, 0.55, sprintf("too few bins (%d)", n), cex = 1.1)
      if (nrow(ln)) {
        tk <- paste(sprintf("%s %+.3f", MABB[ln$method[!duplicated(ln$method)]], ln$line[!duplicated(ln$method)]), collapse = "  ")
        text(0.5, 0.4, tk, cex = 0.6)
      }
      for (q in seq_len(nrow(pr))) lines_log[[length(lines_log) + 1L]] <- data.frame(page = pg, panel = j, ref = r, other = o,
        method = pr$method[q], y_drawn = NA_real_, ylim_lo = NA_real_, ylim_hi = NA_real_)
      next
    }
    ma <- ma_coords(cm[sel, r], cm[sel, o], Tt[[r]], Tt[[o]], pr$prior[1])
    xl <- as.numeric(quantile(ma$A, c(0.001, 0.999)))
    yq <- as.numeric(quantile(ma$M, c(0.005, 0.995)))
    yr <- range(c(yq, ln$line)); pad <- if (diff(yr) > 0) 0.1 * diff(yr) else 0.1
    yl <- c(yr[1] - pad, yr[2] + pad)
    inb <- ma$A >= xl[1] & ma$A <= xl[2] & ma$M >= yl[1] & ma$M <= yl[2]
    smoothScatter(x = ma$A[inb], y = ma$M[inb], xlab = "Mean log2 CPM (A)", ylab = "log2 CPM difference, other - ref (M)",
                  main = "", xlim = xl, ylim = yl, nbin = 128, nrpoints = 100)
    title(main = ttl, cex.main = 0.8, line = 3); mtext(typ, line = 2, cex = 0.7); mtext(sub_t, line = 1, cex = 0.6)
    for (q in seq_len(nrow(ln))) {
      m <- ln$method[q]
      abline(h = ln$line[q], col = MCOL[[m]], lwd = if (isTRUE(ln$is_primary[q])) 2.5 else 1.5,
             lty = if (ln$status[q] %in% c("qc_fail", "imputed")) 3 else 2)
    }
    corner()
    if (is.finite(pr$core_lo[1])) {
      u <- par("usr"); xr <- u[2] - 0.03 * (u[2] - u[1])
      segments(xr, pr$core_lo[1], xr, pr$core_hi[1], lwd = 5, lend = 1)
      points(xr, pr$M_peak[1], pch = 16, col = "white", cex = 0.5)
    }
    for (q in seq_len(nrow(pr))) {
      yd <- if (is.finite(pr$line[q])) pr$line[q] else NA_real_
      lines_log[[length(lines_log) + 1L]] <- data.frame(page = pg, panel = j, ref = r, other = o, method = pr$method[q],
                                                        y_drawn = yd, ylim_lo = yl[1], ylim_hi = yl[2])
    }
  }
  draw_legend()
  par(op)
  invisible(dev.off())
  cat("Done: ", fn, "\n", sep = "")
}
if (nzchar(lines_out)) {
  write.table(do.call(rbind, lines_log), lines_out, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
