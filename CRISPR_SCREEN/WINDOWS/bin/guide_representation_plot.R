#!/usr/bin/env Rscript
#' ##########################################################################
#' guide_representation_plot.R
#'
#' Plotting half of the guide-representation QC: renders every plot from the
#' tables written by guide_representation_wrangle.R. No statistics are
#' computed here beyond what ggplot2 itself does (ECDF steps, densities,
#' loess smooth) — change styling, add or drop plots, and re-run without
#' recomputing any tables.
#'
#' USAGE:
#'   Rscript guide_representation_plot.R --tables-dir=<dir> --prefix=<prefix> \
#'       [--theme=<theme.R>] [--outdir=<dir>] [--dpi=150]
#'
#'   --tables-dir  directory holding <prefix>_*.csv from the wrangle step
#'   --prefix      count-table prefix (e.g. MYSCREEN_library)
#'   --theme       R file defining a ggplot2 theme object named `my_theme`
#'                 (default: my_theme.R next to this script)
#'   --outdir      where to write the plots (default: current directory)
#'   --dpi         PNG resolution (default: 150)
#'
#' Plots: guide representation (ranked), ECDF + mean-read-count bars, Lorenz,
#' density, control vs targeting, 90/10 ratio, undetected guides, mapping rate,
#' read classification, GC vs depth, PCA of samples (>= 3 samples), Pearson and
#' Spearman correlation heatmaps (ordered by the wrangle step's hclust order),
#' read depth (total vs library-matched reads per sample).
#'
#' Every plot is written as an SVG and a PNG twin (same base name, default 150 dpi,
#' same width/height in inches), so viewers that cannot render SVG can open it.
#'
#' Any plot whose input table is absent (e.g. no mapping-rate table because
#' no MAGeCK log was given) is skipped with a note, not an error.
#' ##########################################################################

options(stringsAsFactors = FALSE)

suppressMessages({
  library(dplyr, quietly = TRUE)
  library(tidyr, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
})

# ---- Arguments (--key=value) ----------------------------------------------
parse_args <- function(args) {
  out <- list()
  for (m in regmatches(args, regexec("^--([^=]+)=(.*)$", args))) {
    if (length(m) == 3) out[[m[2]]] <- m[3]
  }
  out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt[["tables-dir"]]) || is.null(opt$prefix)) {
  stop("Usage: guide_representation_plot.R --tables-dir=<dir> --prefix=<prefix> [--theme=<theme.R>] [--outdir=<dir>] [--dpi=150]")
}
tables_dir <- opt[["tables-dir"]]
prefix <- opt$prefix
outdir <- if (is.null(opt$outdir)) "." else opt$outdir
dpi <- as.numeric(if (is.null(opt$dpi)) "150" else opt$dpi)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

theme_file <- opt$theme
if (is.null(theme_file)) {
  self <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  theme_file <- file.path(dirname(normalizePath(self)), "my_theme.R")
}
source(theme_file) # defines `my_theme`

# ---- Helpers --------------------------------------------------------------
short_label <- function(x) {
  out <- character(length(x))
  for (i in seq_along(x)) {
    v <- x[i]
    out[i] <- if (is.na(v)) NA_character_
      else if (abs(v) >= 1e9) paste0(format(v / 1e9, trim = TRUE, scientific = FALSE, digits = 3), "B")
      else if (abs(v) >= 1e6) paste0(format(v / 1e6, trim = TRUE, scientific = FALSE, digits = 3), "M")
      else if (abs(v) >= 1e3) paste0(format(v / 1e3, trim = TRUE, scientific = FALSE, digits = 3), "K")
      else format(v, trim = TRUE, scientific = FALSE, digits = 3)
  }
  out
}

read_table <- function(name) {
  path <- file.path(tables_dir, paste0(prefix, "_", name, ".csv"))
  if (!file.exists(path)) {
    return(NULL)
  }
  read.csv(path, check.names = FALSE)
}
# A ggplot is print()ed; anything else (gtable/grid object) is drawn on a new page.
draw_plot <- function(p) {
  if (inherits(p, "ggplot")) print(p) else { grid::grid.newpage(); grid::grid.draw(p) }
}
# PNG twin of every SVG, for viewers that cannot render SVG (e.g. AI agents).
# ragg (ships with the tidyverse image) if present, else base png() on cairo.
save_png <- function(p, path, width, height, dpi = 150) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = path, width = width, height = height, units = "in", res = dpi)
  } else {
    grDevices::png(filename = path, width = width, height = height, units = "in", res = dpi, type = "cairo")
  }
  draw_plot(p)
  invisible(grDevices::dev.off())
}
save_plot <- function(p, name, width, height) {
  path <- file.path(outdir, paste0(prefix, "_", name, ".svg"))
  grDevices::svg(filename = path, width = width, height = height)
  draw_plot(p)
  invisible(dev.off())
  cat("  Saved:", path, "\n")
  png_path <- sub("\\.svg$", ".png", path)
  save_png(p, png_path, width, height, dpi = dpi)
  cat("  Saved:", png_path, "\n")
}
skip <- function(what, table) cat("  [skip]", what, "— no", paste0(prefix, "_", table, ".csv"), "\n")
centered_title <- theme(plot.title = element_text(hjust = 0.5))
angled_x <- theme(
  plot.title = element_text(hjust = 0.5),
  axis.text.x = element_text(angle = 35, hjust = 1)
)

counts_long <- read_table("counts_long")
if (is.null(counts_long)) stop("Missing ", prefix, "_counts_long.csv in ", tables_dir)
n_samples <- length(unique(counts_long$sample))
bar_width <- max(6, n_samples * 1.2 + 1.5)

# ---- 1. Ranked guide representation ---------------------------------------
ranked <- read_table("ranked")
if (!is.null(ranked)) {
  p <- ggplot(ranked, aes(x = rank, y = log10(avg_cpm + 1))) +
    geom_area(fill = "purple4", alpha = 0.3) +
    geom_line(color = "purple4", linewidth = 1) +
    labs(
      x = "gRNA Rank (lowest to highest representation)",
      y = "Mean CPM + 1 (log10)",
      title = paste("Guide Representation —\n", prefix)
    ) +
    my_theme +
    centered_title
  save_plot(p, "guide_representation", 6, 6)
} else {
  skip("guide representation", "ranked")
}

# ---- 2. ECDF + mean read count per guide ------------------------------------
# Two stacked panels sharing one log10 x axis: the ECDF curves (top) and, below,
# one bar per sample = its mean read count per guide (zeros included), in the
# same colour as its curve. The mean is the count every guide would have if
# reads were split perfectly evenly across the library.
sample_levels <- sort(unique(counts_long$sample))
sample_cols <- setNames(scales::hue_pal()(length(sample_levels)), sample_levels)
sample_means <- counts_long |> group_by(sample) |> summarize(mean_count = mean(count), .groups = "drop")
sample_means$pos <- match(sample_means$sample, sample_levels)
sample_means$pos <- length(sample_levels) + 1 - sample_means$pos # first sample at the top
xl <- c(1, max(max(counts_long$count) + 1, max(sample_means$mean_count + 1) * 4))

p1 <- ggplot(counts_long, aes(x = count + 1, color = sample)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  scale_x_log10(limits = xl, expand = expansion(mult = c(0.02, 0.05))) +
  scale_color_manual(values = sample_cols) +
  labs(
    x = NULL, y = "Cumulative Fraction of gRNAs",
    title = paste("ECDF of gRNA Representation —\n", prefix)
  ) +
  my_theme +
  centered_title +
  theme(legend.position = "right")

p2 <- ggplot(sample_means) +
  geom_rect(aes(xmin = 1, xmax = mean_count + 1, ymin = pos - 0.35, ymax = pos + 0.35, fill = sample)) +
  geom_text(aes(x = mean_count + 1, y = pos,
                label = trimws(format(round(mean_count, 1), big.mark = ",", nsmall = 1))),
            hjust = -0.1, size = 3.2) +
  scale_x_log10(limits = xl, expand = expansion(mult = c(0.02, 0.05))) +
  scale_y_continuous(breaks = sample_means$pos, labels = sample_means$sample,
                     expand = c(0.03, 0.03)) +
  scale_fill_manual(values = sample_cols, guide = "none") +
  labs(x = "gRNA Read Count + 1 (log10)   (bars: mean read count per guide)", y = NULL) +
  my_theme +
  theme(legend.position = "none")

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g <- rbind(g1, g2, size = "max") # gtable method (gtable_rbind is not exported in this gtable version)
panel_rows <- sort(g$layout$t[grepl("^panel", g$layout$name)])
bar_in <- max(1.3, 0.22 * n_samples + 0.4)
g$heights[panel_rows[1]] <- unit(4.5, "in")
g$heights[panel_rows[2]] <- unit(bar_in, "in")
save_plot(g, "ecdf", 8, 4.5 + bar_in + 2.4)

# ---- 3. Lorenz curves -------------------------------------------------------
lorenz <- read_table("lorenz")
if (!is.null(lorenz)) {
  p <- ggplot(lorenz, aes(x = cum_gRNAs, y = cum_reads, color = sample)) +
    geom_line(linewidth = 1.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      x = "Cumulative Fraction of gRNAs", y = "Cumulative Fraction of Total Reads",
      title = paste("Lorenz Curves —\n", prefix)
    ) +
    my_theme +
    centered_title
  save_plot(p, "lorenz", 8, 6)
} else {
  skip("Lorenz curves", "lorenz")
}

# ---- 4. Read-count density overlay -----------------------------------------
p <- ggplot(counts_long, aes(x = count + 1, fill = sample)) +
  geom_density(alpha = 0.3) +
  scale_x_log10() +
  labs(
    x = "Read count + 1 (log10)", y = "Density",
    title = paste("Read Count Density Overlay —\n", prefix)
  ) +
  my_theme +
  centered_title
save_plot(p, "density", 7.5, 6)

# ---- 5. Control vs targeting (faceted by sample) ---------------------------
if (any(as.logical(counts_long$is_control))) {
  p <- ggplot(counts_long, aes(x = count + 1, color = as.logical(is_control))) +
    stat_ecdf(geom = "step", linewidth = 1) +
    scale_x_log10() +
    facet_wrap(~sample) +
    scale_color_manual(
      values = c("FALSE" = "gray70", "TRUE" = "red"),
      labels = c("FALSE" = "Targeting", "TRUE" = "Control")
    ) +
    labs(
      x = "Read count + 1 (log10)", y = "Cumulative Fraction",
      title = paste("Control vs Targeting gRNA Representation —\n", prefix),
      color = "gRNA Type"
    ) +
    my_theme +
    theme(legend.position = "bottom")
  save_plot(p, "controls_ecdf", 8, 8)
} else {
  cat("  [skip] control vs targeting — no control guides (Control/NTC/Non-targeting)\n")
}

# ---- 6. 90/10 percentile ratio ---------------------------------------------
ratio <- read_table("percentile_ratio")
if (!is.null(ratio)) {
  p <- ggplot(ratio, aes(x = sample, y = ratio_90_10, fill = ratio_90_10 > 10)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "forestgreen", linewidth = 0.8) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
    scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"), guide = "none") +
    annotate("text",
      x = -Inf, y = 5, label = "Ideal ≤ 5", hjust = -0.1, vjust = -0.4,
      color = "forestgreen", size = 3.5
    ) +
    annotate("text",
      x = -Inf, y = 10, label = "Flag  ≥ 10", hjust = -0.1, vjust = -0.4,
      color = "firebrick", size = 3.5
    ) +
    labs(
      x = "Sample", y = "90th / 10th Percentile Read Count Ratio",
      title = paste("Library Skew: 90/10 Percentile Ratio —\n", prefix)
    ) +
    my_theme +
    angled_x
  save_plot(p, "90_10_ratio", bar_width, 5)
} else {
  skip("90/10 ratio", "percentile_ratio")
}

# ---- 7. % undetected guides -------------------------------------------------
undetected <- read_table("undetected")
if (!is.null(undetected)) {
  p <- ggplot(undetected, aes(x = sample, y = pct_undetected, fill = pct_undetected > 5)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "forestgreen", linewidth = 0.8) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
    geom_text(aes(label = paste0(pct_undetected, "%\n(", zero_guides, " guides)")),
      vjust = -0.4, size = 3.5
    ) +
    scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"), guide = "none") +
    annotate("text",
      x = -Inf, y = 1, label = "Good ≤ 1%", hjust = -0.1, vjust = -0.4,
      color = "forestgreen", size = 3.5
    ) +
    annotate("text",
      x = -Inf, y = 5, label = "Flag  ≥ 5%", hjust = -0.1, vjust = -0.4,
      color = "firebrick", size = 3.5
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(
      x = "Sample", y = "% Guides with Zero Reads",
      title = paste("Undetected Guides (Zero Counts) —\n", prefix)
    ) +
    my_theme +
    angled_x
  save_plot(p, "undetected_guides", bar_width, 5)
} else {
  skip("undetected guides", "undetected")
}

# ---- 8. Mapping rate --------------------------------------------------------
mapping <- read_table("mapping_rate")
if (!is.null(mapping)) {
  p <- ggplot(mapping, aes(x = label)) +
    geom_col(aes(y = 100), fill = "gray85", width = 0.6) +
    geom_col(aes(y = mapping_pct, fill = mapping_pct < 30), width = 0.6) +
    geom_text(
      aes(
        y = mapping_pct,
        label = paste0(mapping_pct, "%\n(", format(mappedreads, big.mark = ","), " reads)")
      ),
      vjust = -0.4, size = 3.5
    ) +
    geom_hline(yintercept = 30, linetype = "dashed", color = "firebrick", linewidth = 0.8) +
    scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"), guide = "none") +
    annotate("text",
      x = -Inf, y = 30, label = "Flag < 30%", hjust = -0.1, vjust = -0.4,
      color = "firebrick", size = 3.5
    ) +
    scale_y_continuous(limits = c(0, 110), expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    labs(
      x = "Sample", y = "% Reads Mapped to Library (= % Exact Guide Match)",
      title = paste("Mapping Rate / % Exact Guide Match —\n", prefix)
    ) +
    my_theme +
    angled_x
  save_plot(p, "mapping_rate", max(6, nrow(mapping) * 1.2 + 2), 5)
} else {
  skip("mapping rate", "mapping_rate")
}

# ---- 8b. Read classification: exact / near-miss / far-miss -----------------
# mageck count's own matching is always exact, so this breaks the *unmapped*
# reads (mapping_rate's complement) down further: near-miss (1 substitution
# from a real guide -- likely sequencing/PCR error) vs far-miss (not close to
# any guide -- contamination, wrong library, adapter reads, ...).
classification <- read_table("read_classification")
if (!is.null(classification)) {
  class_long <- classification |>
    select(sample, pct_exact_match, pct_near_miss, pct_far_miss) |>
    pivot_longer(-sample, names_to = "category", values_to = "pct") |>
    mutate(category = factor(
      category,
      levels = c("pct_far_miss", "pct_near_miss", "pct_exact_match"),
      labels = c("Far-miss (unrelated)", "Near-miss (1 mismatch)", "Exact match")
    ))
  p <- ggplot(class_long, aes(x = sample, y = pct, fill = category)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = ifelse(pct >= 3, paste0(round(pct), "%"), "")),
      position = position_stack(vjust = 0.5), size = 3, color = "white"
    ) +
    scale_fill_manual(values = c(
      "Exact match"            = "steelblue",
      "Near-miss (1 mismatch)" = "goldenrod2",
      "Far-miss (unrelated)"   = "gray60"
    )) +
    scale_y_continuous(limits = c(0, 100.5), expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    labs(
      x = "Sample", y = "% of Reads", fill = NULL,
      title = paste("Read Classification —\n", prefix)
    ) +
    my_theme +
    angled_x +
    theme(legend.position = "bottom")
  save_plot(p, "read_classification", max(6, nrow(classification) * 1.2 + 2), 5.5)
} else {
  skip("read classification", "read_classification")
}

# ---- 9. GC content vs read depth -------------------------------------------
gc <- read_table("gc_vs_depth")
if (!is.null(gc)) {
  p <- ggplot(gc, aes(x = gc_pct, y = log10(norm_count + 1), color = sample)) +
    geom_point(alpha = 0.15, size = 0.6) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 1.2) +
    labs(
      x = "Spacer GC Content (%)", y = "Median-Normalised Read Count + 1 (log10)",
      title = paste("GC Content vs Read Depth —\n", prefix), color = "Sample"
    ) +
    my_theme +
    centered_title
  save_plot(p, "gc_vs_depth", 7, 6)
} else {
  skip("GC content vs depth", "gc_vs_depth")
}

# ---- 10. PCA of samples -----------------------------------------------------

pca <- read_table("pca")
pca_var <- read_table("pca_variance")
if (!is.null(pca) && !is.null(pca_var)) {
  pc_lab <- function(pc) {
    v <- pca_var$pct_variance[pca_var$PC == pc]
    sprintf("%s (%.1f%% of variance)", pc, v)
  }
  pad <- function(r) {
    d <- diff(r)
    if (d == 0) d <- max(abs(r[1]), 1)
    r + c(-1, 1) * 0.14 * d
  }
  pca_xlim <- pad(range(pca$PC1))
  pca_ylim <- pad(range(pca$PC2))
  p <- ggplot(pca, aes(x = PC1, y = PC2)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    geom_point(aes(fill = sample), shape = 21, size = 4.5, alpha = 0.7, colour = "grey20", stroke = 0.3) +
    scale_fill_manual(values = sample_cols, name = "Sample") +
    guides(fill = guide_legend(override.aes = list(alpha = 0.7, size = 4))) +
    scale_x_continuous(limits = pca_xlim, expand = c(0, 0)) +
    scale_y_continuous(limits = pca_ylim, expand = c(0, 0)) +
    labs(
      x = pc_lab("PC1"), y = pc_lab("PC2"),
      title = paste("PCA of Sample Counts —\n", prefix)
    ) +
    my_theme +
    centered_title
  save_plot(p, "pca", 7, 6)
} else {
  skip("PCA", "pca")
}

# ---- 11. Replicate correlation heatmaps ------------------------------------
cor_long <- read_table("replicate_correlations")
if (!is.null(cor_long)) {
  cor_order <- read_table("correlation_order")
  sample_order <- if (!is.null(cor_order)) {
    cor_order$sample[order(cor_order$order)]
  } else {
    sort(unique(c(cor_long$Sample_A, cor_long$Sample_B)))
  }
  for (method in c("Pearson", "Spearman")) {
    d <- cor_long |>
      rename(Correlation = all_of(method)) |>
      mutate(
        Sample_A = factor(Sample_A, levels = sample_order),
        Sample_B = factor(Sample_B, levels = rev(sample_order))
      )
    cor_lo <- min(0, min(d$Correlation, na.rm = TRUE))
    cor_mid <- (cor_lo + 1) / 2
    p <- ggplot(d, aes(x = Sample_A, y = Sample_B, fill = Correlation)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.2f", Correlation),
                    colour = Correlation > cor_mid), size = 3.5, show.legend = FALSE) +
      scale_colour_manual(values = c("FALSE" = "gray20", "TRUE" = "white")) +
      scale_fill_gradient(
        low = "#f7fbff", high = "#08519c", limits = c(cor_lo, 1),
        name = paste(method, "\ncorrelation")
      ) +
      labs(title = paste(method, "Correlation —\n", prefix), x = NULL, y = NULL) +
      my_theme +
      angled_x
    save_plot(
      p, paste0(tolower(method), "_correlation"),
      max(5, n_samples + 1.5), max(4.5, n_samples + 1)
    )
  }
} else {
  skip("replicate correlation", "replicate_correlations")
}

# ---- 12. Read depth ----------------------------------------------------------
# geom_col has no valid baseline on a log axis, so bars are geom_rect from a
# shared constant floor.
read_depth <- read_table("read_depth")
if (!is.null(read_depth)) {
  cat_levels <- c("Total reads (sequenced)", "Reads matched to this library")
  rd_pos <- data.frame(sample = read_depth$sample, pos = rev(seq_len(nrow(read_depth))))
  rd_long <- read_depth |>
    pivot_longer(c(total_reads, reads_matched_to_library), names_to = "category", values_to = "value") |>
    mutate(category = factor(category, levels = c("total_reads", "reads_matched_to_library"),
                             labels = cat_levels)) |>
    filter(!is.na(value), value > 0) |>
    left_join(rd_pos, by = "sample") |>
    mutate(ymid = pos + ifelse(category == cat_levels[1], 0.2, -0.2))
  if (nrow(rd_long) > 0) {
    rd_min <- min(rd_long$value)
    rd_floor <- if (rd_min < 10) 1 else 10^floor(log10(rd_min))
    p <- ggplot(rd_long) +
      geom_rect(aes(xmin = rd_floor, xmax = value, ymin = ymid - 0.18, ymax = ymid + 0.18, fill = category)) +
      geom_text(aes(x = value, y = ymid, label = trimws(format(round(value), big.mark = ","))),
                hjust = -0.1, size = 3) +
      scale_x_log10(
        limits = c(rd_floor, max(rd_long$value) * 6),
        labels = short_label,
        expand = expansion(mult = c(0, 0.02))
      ) +
      scale_y_continuous(breaks = rd_pos$pos, labels = rd_pos$sample, expand = c(0.03, 0.03)) +
      scale_fill_manual(values = setNames(c("gray65", "steelblue"), cat_levels), drop = FALSE) +
      labs(
        x = "Reads (log10 scale)", y = NULL, fill = NULL,
        title = paste("Read Depth per Sample —\n", prefix)
      ) +
      my_theme +
      centered_title +
      theme(legend.position = "bottom")
    save_plot(p, "read_depth", 8, max(3.5, 0.7 * n_samples + 2))
  } else {
    cat("  [skip] read depth — no positive read counts\n")
  }
} else {
  skip("read depth", "read_depth")
}

cat("\nGuide representation plots complete.\n")
