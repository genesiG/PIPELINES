#!/usr/bin/env Rscript
# classify_regions_plot.R -- plot step of the region classification: MA plot of the classified regions.
# Reads ONLY the tables written by classify_regions.R (<id>.per_region.tsv.gz and <id>.summary.tsv); recomputes no statistics.
#
# CLI contract (called by the CLASSIFY_REGIONS_PLOT module; R style --key=value):
#   --per-region=<id>.per_region.tsv.gz   columns chr, start0, end, name, logFC, logCPM, PValue, FDR, class
#   --summary=<id>.summary.tsv            class, n, Mb, lfc_retained, lfc_lost, fdr, n_ref, n_treat, dispersion, ... (first row's
#                                         thresholds are used)
#   --id=<name>                           classification id (plot title)
#   --out=<file.svg>                      output SVG path (e.g. <id>.ma.svg)
# Output: one SVG. x = average log2-CPM (A), y = log2 fold change treated / reference (M), one colour per class
#   (retained blue, lost red, unclassified grey), dashed horizontal lines at lfc_retained and lfc_lost, dotted line at 0.
#   Classes with > 2000 regions are drawn as a density layer = ONE embedded raster image (options(preferRaster = TRUE) makes image()
#   use rasterImage(); otherwise svg() writes one vector <rect> per density cell) plus the 100 sparsest points as vector dots;
#   smaller classes are drawn as vector points. Axes, text, lines and legend stay vector.
# Exit codes: 0 ok; 2 usage; 3 data problem (missing/empty table).
# Reading the plot: the y position of the cloud depends on the external scale factor applied upstream (no re-centring is done here).
#   "retained" includes every region without significant loss (FDR >= threshold), whatever its logFC.

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

args <- parse_args(commandArgs(trailingOnly = TRUE), c("per-region", "summary", "id", "out", "theme"))
per_f  <- get_req(args, "per-region")
sum_f  <- get_req(args, "summary")
id     <- get_req(args, "id")
out_f  <- get_req(args, "out")
if (!file.exists(per_f)) die(paste0("per-region table not found: ", per_f))
if (!file.exists(sum_f)) die(paste0("summary table not found: ", sum_f))

per <- read.delim(per_f, check.names = FALSE, stringsAsFactors = FALSE)
summ <- read.delim(sum_f, check.names = FALSE, stringsAsFactors = FALSE)
need <- c("logFC", "logCPM", "class")
if (!all(need %in% colnames(per))) die(paste0(per_f, " lacks columns: ", paste(setdiff(need, colnames(per)), collapse = ", ")))
need_s <- c("class", "n", "Mb", "lfc_retained", "lfc_lost", "fdr")
if (!all(need_s %in% colnames(summ))) die(paste0(sum_f, " lacks columns: ", paste(setdiff(need_s, colnames(summ)), collapse = ", ")))
if (nrow(per) == 0L) die("per-region table has no rows")
ok <- is.finite(per$logFC) & is.finite(per$logCPM)
if (!any(ok)) die("no region has a finite logFC and logCPM")
per <- per[ok, , drop = FALSE]
lfc_ret <- summ$lfc_retained[1]; lfc_lost <- summ$lfc_lost[1]; fdr_thr <- summ$fdr[1]
fixed_disp <- "dispersion" %in% colnames(summ) && isTRUE(summ$dispersion[1] == "fixed_0.05")

classes <- c("unclassified", "retained", "lost")          # draw order: lost on top
cols <- c(retained = "#4575b4", lost = "#d73027", unclassified = "#9aa1a9")
xlim <- range(per$logCPM); ylim <- range(c(per$logFC, lfc_ret, lfc_lost, 0))
xpad <- 0.03 * diff(xlim); ypad <- 0.03 * diff(ylim)
xlim <- xlim + c(-xpad, xpad); ylim <- ylim + c(-ypad, ypad)

# density of one class on its own grid; returns expected number of regions per cell and the density at each point
dens_class <- function(x, y) {
  bw <- c(diff(stats::quantile(x, c(0.05, 0.95))), diff(stats::quantile(y, c(0.05, 0.95)))) / 25
  bw[!is.finite(bw) | bw <= 0] <- c(diff(xlim), diff(ylim))[!is.finite(bw) | bw <= 0] / 50
  est <- KernSmooth::bkde2D(cbind(x, y), bandwidth = bw, gridsize = c(256L, 256L))
  cell <- diff(est$x1[1:2]) * diff(est$x2[1:2])
  cnt <- est$fhat * length(x) * cell
  ix <- pmin(pmax(findInterval(x, est$x1), 1L), length(est$x1)); iy <- pmin(pmax(findInterval(y, est$x2), 1L), length(est$x2))
  list(x1 = est$x1, x2 = est$x2, cnt = cnt, at_pt = cnt[cbind(ix, iy)])
}
parts <- list()
for (k in classes) {
  d <- per[per$class == k, , drop = FALSE]
  if (nrow(d) == 0L) next
  parts[[k]] <- list(d = d, dens = if (nrow(d) > 2000L) tryCatch(dens_class(d$logCPM, d$logFC), error = function(e) NULL) else NULL)
}
cmax <- max(c(1, vapply(parts, function(p) if (is.null(p$dens)) 0 else max(p$dens$cnt), 0)))

options(preferRaster = TRUE)
svg(out_f, width = 6.2, height = 5.8)
op <- par(mar = c(4.6, 4.6, 4.2, 1.2))
plot(NA, xlim = xlim, ylim = ylim, xlab = "Average log2-CPM (A)", ylab = "log2 fold change, treated / reference (M)",
     main = id, cex.main = 1)
for (k in names(parts)) {
  p <- parts[[k]]; d <- p$d; col <- cols[[k]]
  if (is.null(p$dens)) {
    points(d$logCPM, d$logFC, pch = 16, cex = 0.4, col = adjustcolor(col, 0.7))
  } else {
    rgbv <- grDevices::col2rgb(col) / 255
    ramp <- grDevices::rgb(rgbv[1], rgbv[2], rgbv[3], alpha = seq(0, 1, length.out = 100))
    # pmin with the matrix FIRST: pmin keeps the attributes (dim) of its first argument
    z <- sqrt(pmin(p$dens$cnt / (0.3 * cmax), 1)); z[p$dens$cnt < 0.01 * cmax] <- 0
    image(p$dens$x1, p$dens$x2, z, col = ramp, zlim = c(0, 1), add = TRUE, useRaster = TRUE)
    sparse <- order(p$dens$at_pt)[seq_len(min(100L, nrow(d)))]
    points(d$logCPM[sparse], d$logFC[sparse], pch = 16, cex = 0.4, col = adjustcolor(col, 0.7))
  }
}
abline(h = 0, col = "grey60", lty = 3)
abline(h = lfc_ret, col = cols[["retained"]], lty = 2, lwd = 1.2)
abline(h = lfc_lost, col = cols[["lost"]], lty = 2, lwd = 1.2)
srow <- function(k) summ[summ$class == k, , drop = FALSE]
lab <- vapply(c("retained", "lost", "unclassified"), function(k) {
  r <- srow(k)
  if (nrow(r) == 0L) return(k)
  sprintf("%s: n = %s (%.2f Mb)", k, format(r$n[1], big.mark = ","), r$Mb[1])
}, "")
legend("topright", bty = "n", cex = 0.75, pch = 16, col = cols[c("retained", "lost", "unclassified")], legend = lab)
mtext(sprintf("lost: logFC <= %g and FDR < %g;  retained: round(logFC, 1) >= %g or FDR >= %g", lfc_lost, fdr_thr, lfc_ret, fdr_thr),
      side = 3, line = 0.4, cex = 0.65)
if (fixed_disp) mtext("no replicates: fixed dispersion, FDR not valid", side = 3, line = 1.4, cex = 0.65, col = "#d73027")
par(op)
invisible(dev.off())
cat("Done: ", out_f, "\n", sep = "")
