#!/usr/bin/env Rscript
#' ##########################################################################
#' plot_barcode_position.R
#'
#' Density plot of where in the read each sample's barcode was matched,
#' from the per-(sample, position) count table written by
#' barcode_position_qc.py. A shared, sharp peak across samples confirms a
#' consistent read structure (stagger + barcode length); a shifted or
#' multi-modal peak flags a misconfigured barcode length or a read-structure
#' problem.
#'
#' The input is already binned (one row per sample x position, with a read
#' count), so each position is weighted by its share of that sample's reads
#' rather than expanding millions of reads back into rows.
#'
#' USAGE:
#'   Rscript plot_barcode_position.R --counts=<counts.csv> --project=<name> \
#'       [--theme=<theme.R>] [--bw=<bandwidth in bp>] [--outdir=<dir>]
#'
#'   --theme   R file defining a ggplot2 theme object named `my_theme`
#'             (default: my_theme.R next to this script)
#'   --bw      kernel bandwidth in bp (default: 0.5 — narrow, since positions
#'             are integers and the expected peak is a single position)
#'   --outdir  where to write the plot (default: current directory); an SVG
#'             and a PNG twin (same base name, 150 dpi, same size in inches)
#' ##########################################################################

options(stringsAsFactors = FALSE)

suppressMessages({
  library(dplyr,   quietly = TRUE)
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
if (is.null(opt$counts) || is.null(opt$project)) {
  stop("Usage: plot_barcode_position.R --counts=<counts.csv> --project=<name> [--theme=<theme.R>] [--bw=<bp>] [--outdir=<dir>]")
}
bw     <- as.numeric(if (is.null(opt$bw)) "0.5" else opt$bw)
outdir <- if (is.null(opt$outdir)) "." else opt$outdir

# Default theme: my_theme.R next to this script (found via Rscript's own
# --file= argument; --theme overrides it).
theme_file <- opt$theme
if (is.null(theme_file)) {
  self <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  theme_file <- file.path(dirname(normalizePath(self)), "my_theme.R")
}
source(theme_file)   # defines `my_theme`

# ---- Data -----------------------------------------------------------------
df <- read.csv(opt$counts) |>
  group_by(sample) |>
  mutate(weight = count / sum(count)) |>
  ungroup()

peaks <- df |>
  group_by(sample) |>
  slice_max(count, n = 1, with_ties = FALSE) |>
  ungroup()

cat("Peak barcode position per sample:\n")
for (i in seq_len(nrow(peaks))) {
  cat(sprintf("  %s: %dbp (%.1f%% of reads)\n",
              peaks$sample[i], peaks$position[i], peaks$weight[i] * 100))
}
if (length(unique(peaks$position)) > 1) {
  cat("  [NOTE] Peak position varies across samples — check barcode length/stagger assumptions if unexpected.\n")
}

# ---- Plot -----------------------------------------------------------------
p <- ggplot(df, aes(x = position, weight = weight, color = sample, fill = sample)) +
  geom_density(bw = bw, alpha = 0.15, linewidth = 0.8) +
  labs(
    title    = "Barcode Match Position Within Read",
    subtitle = "A shared, sharp peak across samples confirms a consistent read structure",
    x        = "Position in read (bp, 0-based)",
    y        = "Density",
    color    = "Sample",
    fill     = "Sample"
  ) +
  my_theme

# PNG twin of every SVG, for viewers that cannot render SVG (e.g. AI agents).
# ragg (ships with the tidyverse image) if present, else base png() on cairo.
save_png <- function(p, path, width, height, dpi = 150) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename = path, width = width, height = height, units = "in", res = dpi)
  } else {
    grDevices::png(filename = path, width = width, height = height, units = "in", res = dpi, type = "cairo")
  }
  print(p)
  invisible(grDevices::dev.off())
}

out_path <- file.path(outdir, paste0(opt$project, "_barcode_position_density.svg"))
# ggsave(device = "svg") needs svglite, which rocker/tidyverse lacks — use
# base R's own SVG device instead.
svg(out_path, width = 9, height = 5)
print(p)
invisible(dev.off())
cat("Saved:", out_path, "\n")
png_path <- sub("\\.svg$", ".png", out_path)
save_png(p, png_path, 9, 5)
cat("Saved:", png_path, "\n")
