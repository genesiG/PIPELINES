#!/usr/bin/env Rscript
# qc_alignment_plot.R -- alignment QC plots, drawn only from alignment_stats.tsv.
#
# Called by QC_ALIGNMENT_PLOT (modules/local/qc_alignment_plot.nf):
#   qc_alignment_plot.R --table=alignment_stats.tsv --theme=my_theme.R
#
# Arguments (--key=value): --table (alignment_stats.tsv from qc_alignment_table.R), --theme (my_theme.R;
#   minimal built-in theme when absent), --outdir (default ".")
# Outputs (always all seven; a plot without data becomes a one-line "not available" figure):
#   alignment_read_composition.svg      raw reads (millions): target / chrM / spike-in / unmapped, stacked
#   alignment_read_composition_pct.svg  same, as % of all raw reads
#   alignment_pct_spikein.svg           % spike-in of mapped reads, one panel per antibody
#   alignment_n_spikein.svg             spike-in read count per sample (reference line at 1,000)
#   alignment_pct_chrM.svg              % chrM of mapped reads
#   alignment_duplicates.svg            Picard duplicate fractions (unique / PCR / optical read pairs)
#   alignment_final_reads.svg           reads in the final BAM (millions)
# Colours: categories use the theme palette; antibodies cycle through the theme palette (from the data).

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))
parse_args <- function(args) {
  out <- list()
  for (m in regmatches(args, regexec("^--([^=]+)=(.*)$", args))) if (length(m) == 3) out[[m[2]]] <- m[3]
  out
}
die <- function(...) { message("ERROR: ", ...); quit(save = "no", status = 3) }
opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt$table)) die("usage: qc_alignment_plot.R --table=alignment_stats.tsv [--theme=my_theme.R] [--outdir=.]")
outdir <- if (is.null(opt$outdir)) "." else opt$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
if (!is.null(opt$theme) && file.exists(opt$theme)) source(opt$theme)
base_theme <- if (exists("theme_paper")) theme_paper() else theme_classic(base_size = 11)
pal_get <- function(k, d) if (exists("pal") && is.list(pal) && !is.null(pal[[k]])) pal[[k]] else d
cycle_cols <- function(keys) {
  base <- c(pal_get("purple", "#86187f"), pal_get("steel", "#6b93ad"), pal_get("slate", "#425b74"),
            pal_get("grey", "#a7b4c2"), pal_get("navy", "#1d3557"))
  setNames(rep_len(base, length(keys)), keys)
}

d <- read.delim(opt$table, check.names = FALSE, na.strings = c("NA", ""))
if (nrow(d) == 0) die("empty table ", opt$table)
is_ctrl <- grepl("^igg$|^input$", d$antibody, ignore.case = TRUE)
d <- d[order(is_ctrl, d$antibody, d$group, d$sample_id), ]
ab_levels <- unique(d$antibody)
ab_cols <- cycle_cols(ab_levels)
d$sample_f <- factor(d$sample_id, levels = rev(d$sample_id))   # first sample at the top of horizontal plots
d$antibody <- factor(d$antibody, levels = ab_levels)
n <- nrow(d)
hgt <- max(4, 0.22 * n + 1.8)

save_svg <- function(p, name, w, h) {
  grDevices::svg(filename = file.path(outdir, name), width = w, height = h)
  print(p)
  invisible(dev.off())
}
unavailable <- function(name, msg, w = 7, h = 3) {
  p <- ggplot() + annotate("text", x = 0, y = 0, label = msg, size = 4) + theme_void()
  save_svg(p, name, w, h)
}

# ---- composition (raw) ----------------------------------------------------------
dr <- d[!is.na(d$raw_total_reads), , drop = FALSE]
cat_cols <- c("Target (mapped)" = pal_get("purple", "#86187f"), "chrM" = pal_get("steel", "#6b93ad"),
              "Spike-in" = pal_get("grey", "#a7b4c2"), "Unmapped" = pal_get("ink", "#222222"))
if (nrow(dr) > 0) {
  long <- do.call(rbind, lapply(seq_len(nrow(dr)), function(i) data.frame(
    sample_id = dr$sample_id[i],
    category = factor(names(cat_cols), levels = rev(names(cat_cols))),
    reads = c(dr$n_target_mapped[i], dr$n_chrM[i], dr$n_spike_in[i], dr$raw_total_reads[i] - dr$raw_mapped_reads[i]))))
  long$reads[is.na(long$reads)] <- 0
  long$sample_id <- factor(long$sample_id, levels = dr$sample_id)
  long$pct <- 100 * long$reads / ave(long$reads, long$sample_id, FUN = sum)
  w <- max(6, 0.28 * nrow(dr) + 3)
  p1 <- ggplot(long, aes(sample_id, reads / 1e6, fill = category)) +
    geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = cat_cols, name = NULL, breaks = names(cat_cols)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Raw reads (millions)") + base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  save_svg(p1, "alignment_read_composition.svg", w, 5)
  p2 <- ggplot(long, aes(sample_id, pct, fill = category)) +
    geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = cat_cols, name = NULL, breaks = names(cat_cols)) +
    scale_y_continuous(breaks = seq(0, 100, 20), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.02))) +
    labs(x = NULL, y = "Share of raw reads") + base_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  save_svg(p2, "alignment_read_composition_pct.svg", w, 5)
} else {
  for (nm in c("alignment_read_composition.svg", "alignment_read_composition_pct.svg"))
    unavailable(nm, "Not available: no raw alignment statistics (all samples are pre-aligned BAMs).")
}

# ---- per-sample single-metric horizontal bars -----------------------------------
hbar <- function(df, col, xlab, name, vline = NULL, vlabel = NULL, facet = FALSE, w = 7) {
  df <- df[!is.na(df[[col]]), , drop = FALSE]
  if (nrow(df) == 0) return(unavailable(name, paste0("Not available: no values for '", col, "'.")))
  df$val <- df[[col]]
  p <- ggplot(df, aes(val, sample_f, fill = antibody)) + geom_col(width = 0.7) +
    scale_fill_manual(values = ab_cols, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = xlab, y = NULL) + base_theme + theme(axis.text.y = element_text(size = 7))
  if (facet) p <- p + facet_wrap(~antibody, scales = "free_y", ncol = 1)
  if (!is.null(vline)) {
    p <- p + geom_vline(xintercept = vline, linetype = "dashed", colour = pal_get("ink", "#222222"), linewidth = 0.5) +
      annotate("text", x = vline, y = Inf, label = vlabel, hjust = -0.05, vjust = 1.5, size = 3)
  }
  h <- if (facet) max(4, 0.22 * nrow(df) + 1.2 * length(unique(df$antibody)) + 1.5) else hgt
  save_svg(p, name, w, h)
}
hbar(d, "pct_spike_in", "Spike-in reads (% of mapped reads)", "alignment_pct_spikein.svg", facet = TRUE)
hbar(d, "n_spike_in", "Spike-in reads (n)", "alignment_n_spikein.svg", vline = 1000, vlabel = "1,000 reads")
hbar(d, "pct_chrM", "chrM reads (% of mapped reads)", "alignment_pct_chrM.svg")
hbar(d, "final_mapped_reads", "Reads in final BAM", "alignment_final_reads.svg")

# ---- duplicates (Picard) -----------------------------------------------------------
dd <- d[!is.na(d$dup_pairs_examined) & d$dup_pairs_examined > 0, , drop = FALSE]
if (nrow(dd) > 0) {
  dup_cols <- c("Unique pairs" = pal_get("purple", "#86187f"), "PCR duplicates" = pal_get("steel", "#6b93ad"),
                "Optical duplicates" = pal_get("ink", "#222222"))
  long <- do.call(rbind, lapply(seq_len(nrow(dd)), function(i) {
    ex <- dd$dup_pairs_examined[i]; du <- dd$dup_pairs_duplicate[i]; op <- dd$dup_pairs_optical[i]
    if (is.na(op)) op <- 0
    data.frame(sample_f = dd$sample_f[i], category = factor(names(dup_cols), levels = rev(names(dup_cols))),
               pct = 100 * c(ex - du, du - op, op) / ex)
  }))
  p <- ggplot(long, aes(pct, sample_f, fill = category)) + geom_col(width = 0.75, colour = "white", linewidth = 0.2) +
    scale_fill_manual(values = dup_cols, name = NULL, breaks = names(dup_cols)) +
    scale_x_continuous(breaks = seq(0, 100, 20), labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Share of read pairs examined by Picard MarkDuplicates", y = NULL) + base_theme +
    theme(axis.text.y = element_text(size = 7))
  save_svg(p, "alignment_duplicates.svg", 7.5, max(4, 0.22 * nrow(dd) + 1.8))
} else {
  unavailable("alignment_duplicates.svg", "Not available: no Picard MarkDuplicates metrics.")
}
cat("alignment plots written to", normalizePath(outdir), "\n")
