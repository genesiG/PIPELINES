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
#'       [--theme=<theme.R>] [--outdir=<dir>] [--dpi=300]
#'
#'   --tables-dir  directory holding <prefix>_*.csv from the wrangle step
#'   --prefix      count-table prefix (e.g. MYSCREEN_library)
#'   --theme       R file defining a ggplot2 theme object named `my_theme`
#'                 (default: my_theme.R next to this script)
#'   --outdir      where to write PNGs (default: current directory)
#'   --dpi         PNG resolution (default: 300)
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
  stop("Usage: guide_representation_plot.R --tables-dir=<dir> --prefix=<prefix> [--theme=<theme.R>] [--outdir=<dir>] [--dpi=300]")
}
tables_dir <- opt[["tables-dir"]]
prefix <- opt$prefix
outdir <- if (is.null(opt$outdir)) "." else opt$outdir
dpi <- as.numeric(if (is.null(opt$dpi)) "300" else opt$dpi)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

theme_file <- opt$theme
if (is.null(theme_file)) {
  self <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  theme_file <- file.path(dirname(normalizePath(self)), "my_theme.R")
}
source(theme_file) # defines `my_theme`

# ---- Helpers --------------------------------------------------------------
read_table <- function(name) {
  path <- file.path(tables_dir, paste0(prefix, "_", name, ".csv"))
  if (!file.exists(path)) {
    return(NULL)
  }
  read.csv(path, check.names = FALSE)
}
save_plot <- function(p, name, width, height) {
  path <- file.path(outdir, paste0(prefix, "_", name, ".svg"))
  grDevices::svg(filename = path, width = width, height = height)
  print(p)
  invisible(dev.off())
  cat("  Saved:", path, "\n")
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
  p <- ggplot(ranked, aes(x = rank, y = log10(avg_count + 1))) +
    geom_area(fill = "purple4", alpha = 0.3) +
    geom_line(color = "purple4", linewidth = 1) +
    labs(
      x = "gRNA Rank (lowest to highest representation)",
      y = "Average Read Count + 1 (log10)",
      title = paste("Guide Representation —\n", prefix)
    ) +
    my_theme +
    centered_title
  save_plot(p, "guide_representation", 6, 6)
} else {
  skip("guide representation", "ranked")
}

# ---- 2. ECDF ----------------------------------------------------------------
# The dashed vertical line is each sample's own mean read count — the count
# every guide would have if reads were split perfectly evenly across the
# library. It's the sample's own value (not a fixed cutoff) because total
# depth varies between samples/pools, so only a sample's own mean is a fair
# reference for it. This also makes the plot interpretable with a single
# sample: without it, one lone curve has no reference point at all.
sample_means <- counts_long |> group_by(sample) |> summarize(mean_count = mean(count), .groups = "drop")
p <- ggplot(counts_long, aes(x = count + 1, color = sample)) +
  stat_ecdf(geom = "step", linewidth = 1.2) +
  geom_vline(data = sample_means, aes(xintercept = mean_count + 1, color = sample),
             linetype = "dashed", linewidth = 0.7, alpha = 0.8, show.legend = FALSE) +
  scale_x_log10() +
  labs(
    x = "gRNA Read Count + 1 (log10)", y = "Cumulative Fraction of gRNAs",
    title = paste("ECDF of gRNA Representation —\n", prefix)
  ) +
  my_theme +
  centered_title
save_plot(p, "ecdf", 8, 6)

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

# ---- 10. Replicate correlation heatmaps ------------------------------------
cor_long <- read_table("replicate_correlations")
if (!is.null(cor_long)) {
  for (method in c("Pearson", "Spearman")) {
    d <- cor_long |> rename(Correlation = all_of(method))
    p <- ggplot(d, aes(x = Sample_A, y = Sample_B, fill = Correlation)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = round(Correlation, 3)), size = 3.5) +
      scale_fill_gradient2(
        low = "firebrick", mid = "white", high = "steelblue", midpoint = 0.9,
        limits = c(max(0, min(d$Correlation, na.rm = TRUE) - 0.02), 1),
        name = paste(method, "\nCorrelation")
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

cat("\nGuide representation plots complete.\n")
