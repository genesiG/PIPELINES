# =============================================================================
# my_theme.R - Shared figure style for the paper (ggplot2), matching my_theme_py.py
#
# Canonical copy: /path/to/utils/my_theme.R
# Use from any project (source the canonical file; do not copy it):
#
#   source("/path/to/utils/my_theme.R")
#   p <- ggplot(d, aes(group, fold, colour = line)) +
#     geom_no_preference() + geom_line() + geom_point() +
#     scale_x_retention_group() + scale_y_fold_vs_even() +
#     theme_paper()                                # panel.border on by default
#   save_paper(p, "plots/my_figure", width = 6, height = 4)
#
# Example: /path/to/utils/my_theme_demo.R (panel A of erasure_axis_lines_fold_vs_even in ggplot2).
# =============================================================================

suppressPackageStartupMessages(library(ggplot2))

# --- Palette (identical hex codes to my_theme_py.py) ----------------------------
pal <- list(
  purple   = "#86187f",  # H3K27M line, EZH2i, most retained end of the diverging scale
  steel    = "#6b93ad",  # untreated / DMSO reference, IDH-mutant lines
  slate    = "#425b74",  # K27M-KO DMSO genome track
  grey     = "#a7b4c2",  # K27M-KO EZH2i genome track
  navy     = "#1d3557",  # least retained end of the diverging scale
  ink      = "#222222",  # text, panel borders
  muted    = "#5f6670",  # secondary text
  rule     = "#9aa1a9",  # guides, baselines
  unscored = "#dfe3e7"   # bins excluded from scoring
)

# Treatment comparisons (validated all-pairs: worst CVD dE 14.3, contrast >= 3:1).
# Inhibitor colours match the inhibitor viability figure (2026-10-09): EPZ6438 black, EEDi-5285 dark purple,
# MAK683 blue-lilac; DMSO neutral grey
treatment_colors    <- c("DMSO" = "#808080", "EZH2i" = "#000000", "EEDi 1" = "#76069a", "EEDi 2" = "#8d8dff")
treatment_linetypes <- c("DMSO" = "dashed", "EZH2i" = "solid", "EEDi 1" = "solid", "EEDi 2" = "solid")

# --- Labels and sizes ---------------------------------------------------------
n_groups        <- 20
retention_xlab  <- "Retention Group"
end_labels      <- c("Least retained", "Most retained")
fold_breaks     <- c(0.25, 0.5, 1, 2, 4)
fold_labels     <- c("0.25x", "0.5x", "1x", "2x", "4x")
size_pt         <- list(base = 11, axis_title = 13, legend = 9, annotation = 8, panel_letter = 16)
lw_pt           <- list(axis = 0.6, border = 0.8, data_main = 2.5, data_other = 1.8, reference = 0.8)

pt2lw <- function(pt) pt / .pt          # line width in points -> ggplot2 linewidth units
pt2sz <- function(pt) pt / .pt          # font size in points  -> geom_text/annotate size

# --- Theme --------------------------------------------------------------------
theme_paper <- function(base_size = size_pt$base, border = TRUE) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line         = element_line(colour = pal$ink, linewidth = pt2lw(lw_pt$axis)),
      axis.ticks        = element_line(colour = pal$ink, linewidth = pt2lw(lw_pt$axis)),
      axis.text         = element_text(colour = pal$ink, size = base_size),
      axis.title        = element_text(colour = pal$ink, size = size_pt$axis_title),
      legend.title      = element_blank(),
      legend.text       = element_text(colour = pal$ink, size = size_pt$legend),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      panel.border      = if (border) element_rect(colour = pal$ink, fill = NA,
                                                   linewidth = pt2lw(lw_pt$border)) else element_blank(),
      plot.tag          = element_text(size = size_pt$panel_letter, face = "bold"),
      strip.background  = element_blank(),
      strip.text        = element_text(colour = pal$ink, size = base_size)
    )
}

# --- Scales and layers ---------------------------------------------------------
scale_x_retention_group <- function(name = retention_xlab, ...) {
  scale_x_continuous(name = name, breaks = c(1, 5, 10, 15, 20), limits = c(0.5, n_groups + 0.5),
                     expand = c(0, 0), ...)
}

scale_y_fold_vs_even <- function(name = "Share of H3K27me3\n(fold over expected 5%)", limits = c(0.4, 4.5), ...) {
  scale_y_continuous(name = name, trans = "log2", breaks = fold_breaks, labels = fold_labels,
                     limits = limits, ...)
}

# fill aesthetic must be log2(share fold over expected 5%)
scale_fill_share_vs_even <- function(name = "Share of signal\n(fold over expected 5%)", ...) {
  scale_fill_gradientn(name = name, colours = c(pal$navy, "#f7f7f7", pal$purple),
                       limits = c(-2.3, 2.3), breaks = log2(fold_breaks), labels = fold_labels,
                       oob = scales::squish, ...)
}

geom_no_preference <- function(label = TRUE) {
  layers <- list(geom_hline(yintercept = 1, linetype = "dashed", colour = "#333333",
                            linewidth = pt2lw(lw_pt$reference)))
  if (label) layers <- c(layers, list(annotate("text", x = 1, y = 1.05, label = "Expected share (5%)",
                                               hjust = 0, vjust = 0, size = pt2sz(size_pt$annotation))))
  layers
}

# --- Saving ---------------------------------------------------------------------
save_paper <- function(plot, file_base, width, height, units = "in", dpi = 300) {
  dir.create(dirname(file_base), recursive = TRUE, showWarnings = FALSE)
  svg_dev <- if (requireNamespace("svglite", quietly = TRUE)) svglite::svglite else grDevices::svg
  ggsave(paste0(file_base, ".svg"), plot, width = width, height = height, units = units, device = svg_dev)
  ggsave(paste0(file_base, ".png"), plot, width = width, height = height, units = units, dpi = dpi)
  message("Saved: ", file_base, ".svg/.png")
}
