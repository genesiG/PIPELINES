#!/usr/bin/env Rscript
# qc_replicate_correlation_plot.R -- plots of the replicate concordance step.
#
# Called by REPLICATE_CORRELATION_PLOT (modules/local/replicate_correlation_plot.nf). Draws from the two tables of
# qc_replicate_correlation.R; the bin counts are read only for the scatter plots (x = log2 CPM of the two replicates).
#
# CLI (--key=value):
#   --table=<replicate_correlation.tsv>         replicate_correlation.tsv from qc_replicate_correlation.R
#   --pairs-all=<replicate_correlation_pairs_all.tsv>
#   --counts=<bins_<size>.counts.tsv.gz> --libsizes=<bins_<size>.libsizes.tsv>   the counts the table was computed from
#   --theme=my_theme.R                          (optional) lab theme; minimal built-in theme when absent
#   --pairs-per-page=12                         scatter panels per SVG page (3 columns)
#   --outdir=.
# Outputs:
#   replicate_correlation.svg    one panel column per antibody, rows pearson_all and pearson_signal. x = group; filled points =
#                                replicate pairs; grey jittered points = non-replicate same-antibody pairs involving that group.
#                                No reference lines.
#   replicate_scatter.svg        smoothScatter of log2 CPM of sample_1 (x) vs sample_2 (y) per replicate pair, identity line, text
#                                with pearson_all / spearman_all / pearson_signal; further pages replicate_scatter.p2.svg, ...
#                                The density is one embedded raster (options(preferRaster = TRUE)).
# Exit codes: 0 ok; 2 usage; 3 data problem.

options(stringsAsFactors = FALSE, preferRaster = TRUE)
suppressPackageStartupMessages(library(ggplot2))

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

allowed <- c("table", "pairs-all", "counts", "libsizes", "theme", "pairs-per-page", "outdir")
args <- parse_args(commandArgs(trailingOnly = TRUE), allowed)
tab_file   <- get_req(args, "table")
pairs_file <- get_req(args, "pairs-all")
counts_file <- get_req(args, "counts")
lib_file   <- get_req(args, "libsizes")
per_page   <- as.integer(get_opt(args, "pairs-per-page", "12"))
outdir     <- get_opt(args, "outdir", ".")
if (is.na(per_page) || per_page < 1L) die("--pairs-per-page must be a positive integer", 2L)
for (f in c(tab_file, pairs_file, counts_file, lib_file)) if (!file.exists(f)) die(paste0("file not found: ", f))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
if (!is.null(args$theme) && file.exists(args$theme)) source(args$theme)
base_theme <- if (exists("theme_paper")) theme_paper() else theme_classic(base_size = 11)
pal_get <- function(k, d) if (exists("pal") && is.list(pal) && !is.null(pal[[k]])) pal[[k]] else d

tab <- read.delim(tab_file, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
pa  <- read.delim(pairs_file, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))
if (nrow(tab) == 0L) die("empty table ", tab_file)

# ---- replicate_correlation.svg ------------------------------------------------
grp_of <- character(0)
for (i in seq_len(nrow(tab))) {
  grp_of[tab$sample_1[i]] <- tab$group[i]
  if (!is.na(tab$sample_2[i])) grp_of[tab$sample_2[i]] <- tab$group[i]
}
metrics <- c(pearson_all = "Pearson, all bins", pearson_signal = "Pearson, top-signal bins")
rep_rows <- tab[tab$status == "ok", , drop = FALSE]
mk_rep <- function(mt) if (nrow(rep_rows)) data.frame(antibody = rep_rows$antibody, group = rep_rows$group,
                                                      value = rep_rows[[mt]], metric = metrics[[mt]]) else NULL
rep_df <- do.call(rbind, lapply(names(metrics), mk_rep))
non <- pa[!pa$same_group, , drop = FALSE]
mk_non <- function(mt) {
  if (!nrow(non)) return(NULL)
  rbind(data.frame(antibody = non$antibody, group = grp_of[non$sample_1], value = non[[mt]], metric = metrics[[mt]]),
        data.frame(antibody = non$antibody, group = grp_of[non$sample_2], value = non[[mt]], metric = metrics[[mt]]))
}
non_df <- do.call(rbind, lapply(names(metrics), mk_non))
ab_levels <- unique(tab$antibody)
grp_levels <- unique(tab$group)
fix_df <- function(d) {
  if (is.null(d)) return(NULL)
  d <- d[!is.na(d$value) & !is.na(d$group), , drop = FALSE]
  d$antibody <- factor(d$antibody, levels = ab_levels)
  d$group <- factor(d$group, levels = grp_levels)
  d$metric <- factor(d$metric, levels = unname(metrics))
  d
}
rep_df <- fix_df(rep_df); non_df <- fix_df(non_df)
p <- ggplot() + facet_grid(metric ~ antibody, scales = "free_x", space = "free_x") + base_theme +
  labs(x = NULL, y = "Pearson correlation of log2 CPM") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
if (!is.null(non_df) && nrow(non_df))
  p <- p + geom_jitter(data = non_df, aes(group, value), width = 0.18, height = 0, colour = "#a7b4c2", size = 1.4, alpha = 0.8)
if (!is.null(rep_df) && nrow(rep_df))
  p <- p + geom_point(data = rep_df, aes(group, value), shape = 21, size = 2.6, stroke = 0.6,
                      colour = pal_get("ink", "#222222"), fill = pal_get("purple", "#86187f"))
if (is.null(rep_df) || !nrow(rep_df))
  p <- p + annotate("text", x = 1, y = 0.5, label = "no replicate pairs", size = 4)
n_grp <- length(grp_levels)
grDevices::svg(filename = file.path(outdir, "replicate_correlation.svg"), width = max(5, 0.45 * n_grp + 1.5 * length(ab_levels) + 1.5), height = 6.5)
print(p)
invisible(dev.off())

# ---- replicate_scatter*.svg ---------------------------------------------------
pw <- 3.9; ph <- 3.9; ncol_p <- 3L
lib <- read.delim(lib_file, check.names = FALSE, stringsAsFactors = FALSE)
Tl <- setNames(as.numeric(lib$csaw_total), lib$sample_id)
prior <- tab$prior[1]
if (nrow(rep_rows) == 0L) {
  grDevices::svg(file.path(outdir, "replicate_scatter.svg"), width = 6, height = 2.5)
  plot.new(); text(0.5, 0.5, "no replicate pairs (every group has one library)")
  invisible(dev.off())
  quit(save = "no", status = 0L)
}
cnt <- read.delim(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
cnt <- cnt[!as.logical(cnt$blacklisted), , drop = FALSE]
need_ids <- unique(c(rep_rows$sample_1, rep_rows$sample_2))
if (!all(need_ids %in% colnames(cnt))) die("counts lack columns for some replicate-pair libraries")
if (!all(need_ids %in% names(Tl))) die("libsizes lack some replicate-pair libraries")

n_pairs <- nrow(rep_rows)
n_page <- ceiling(n_pairs / per_page)
fmt <- function(v) if (is.na(v)) "NA" else sprintf("%.3f", v)
for (pg in seq_len(n_page)) {
  idx <- ((pg - 1L) * per_page + 1L):min(n_pairs, pg * per_page)
  nc <- min(ncol_p, length(idx)); nr <- ceiling(length(idx) / ncol_p)
  fname <- if (pg == 1L) "replicate_scatter.svg" else sprintf("replicate_scatter.p%d.svg", pg)
  grDevices::svg(file.path(outdir, fname), width = pw * nc, height = ph * nr)
  par(mfrow = c(nr, nc), mar = c(4, 4, 3, 1), mgp = c(2.3, 0.8, 0))
  for (i in idx) {
    r <- rep_rows[i, ]
    c1 <- cnt[[r$sample_1]]; c2 <- cnt[[r$sample_2]]
    keep <- (c1 + c2) > 0
    c1 <- c1[keep]; c2 <- c2[keep]
    T1 <- Tl[[r$sample_1]]; T2 <- Tl[[r$sample_2]]
    pm <- mean(c(T1, T2)); p1 <- prior * T1 / pm; p2 <- prior * T2 / pm
    x1 <- log2((c1 + p1) / (T1 + 2 * p1) * 1e6); x2 <- log2((c2 + p2) / (T2 + 2 * p2) * 1e6)
    rng <- range(c(x1, x2))
    smoothScatter(x1, x2, nrpoints = 100, xlim = rng, ylim = rng, xlab = paste0(r$sample_1, " (log2 CPM)"),
                  ylab = paste0(r$sample_2, " (log2 CPM)"), main = r$group, cex.main = 0.9, cex.lab = 0.8)
    abline(0, 1, lty = 2, col = "black")
    legend("topleft", bty = "n", cex = 0.8, text.col = "black",
           legend = c(paste0("pearson_all ", fmt(r$pearson_all)), paste0("spearman_all ", fmt(r$spearman_all)),
                      paste0("pearson_signal ", fmt(r$pearson_signal))))
  }
  invisible(dev.off())
}
