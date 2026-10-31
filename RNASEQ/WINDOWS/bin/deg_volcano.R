#!/usr/bin/env Rscript
# deg_volcano.R -- volcano plots of all comparisons of a model (individual svg, multi-page pdfs, combined figure).
#
# Legacy mapping (step_6_deg_analysis.Rmd): the volcano parts of the loop in "Statistical test (Walds Test)" (individual
#   svgs, volcano_plots_{yscaled,xscaled,freescale}.pdf), "Publication-Ready Volcano Plots" (combined figure). The plot function is makeVolcanoPlot() of deg_common.R.
#
# CLI (all --key=value; booleans TRUE/FALSE):
#   --common=deg_common.R
#   --model=<name>
#   --comparisons=comparisons.tsv     normalised table; rows of --model, in table order
#   --alpha=0.05 --lfc=1              dashed guide lines: horizontal at -log10(alpha) in the raw plot, vertical at +-lfc
#   --sval-threshold=0.005            horizontal guide line of the SHRUNK plot at -log10(sval-threshold) (the classification
#                                     itself is read from the csv). The Rmd drew the shrunk plot's line at -log10(alpha) =
#                                     1.30 although shrunk_effect uses svalue < 0.005 (2.30), so s-values between 0.005 and
#                                     0.05 appeared above the line but were coloured NS.
#   --legacy-shrunk-guide=FALSE       TRUE reproduces the Rmd line (-log10(alpha)) in the shrunk plot
#   --shrink=TRUE                     also draw the shrunk plots (needs shrunkLFC / svalue / shrunk_effect columns)
#   --color-down=#46627A --color-up=#8c6bb1 --color-ns=gray
#                                     (legacy colours were #4d004b / goldenrod / gray)
#   --label-genes="G1,G2"             genes that get a label and a black ring (comma list; empty = none)
#   --xlim="" --ylim=""               default per-plot axis limits "min,max"; a comparison row's xlim / ylim column wins;
#                                     empty = automatic
#   --common-xlim="" --common-ylim="" limits of the x-/y-scaled pdfs and the combined figure; empty = computed for the
#                                     pdfs as in the commented-out Rmd code (x: range of all log2FoldChange, y: 0 to
#                                     -log10 of the smallest padj of the model) and automatic for the combined figure
#   --combined-layout=stacked|side
#   --svg-device=svglite              svglite (real <text>, default) or cairo (legacy, glyph outlines)
#   --font=Helvetica                  font family of every svg / pdf (see deg_common.R set_font)
#   --legend=TRUE --legend-position=bottom   legend on/off and position (bottom|top|left|right); the legacy plot had it at the bottom
#   --label-size=4 --label-nudge-x=1.5 --label-nudge-y=1.5 --label-force=1 --label-force-pull=1 --label-fontface=italic
#                                     ggrepel label settings (Rmd values). The label start is nudged by (nudge-x, nudge-y)
#   --label-keep-in-bounds=TRUE       a label whose nudged start would fall outside the axis limits is nudged to the other
#                                     side / shortened (FALSE = legacy: ggplot2 dropped that label, e.g. CDKN2A at -log10(padj)
#                                     25.6 under ylim 0,26, leaving only the ring)
#   --point-size=2 --ring-stroke=0.5  point size and the width of the black ring around labelled genes
#   --width=4.5 --height=5            size (in) of the individual svgs; the combined figure uses width x height per panel
#   --pdf-width=4 --pdf-height=5      size of the y-scaled / free-scale pdfs (x-scaled pdf: pdf-width + 0.5 as the Rmd 4.5)
#   --individual=TRUE --pdfs=TRUE --combined=TRUE   switch each output group
#   --results-dir=.                   folder with the staged results_<comparison>.csv files (default: current directory)
#   --out-dir=.
# Rmd constants kept: force = 1, force_pull = 1, label_size = 4, nudgex = nudgey = 1.5, legend at the bottom, point
#   colour order (down, NS, up) with legend labels Downregulated / NS / Upregulated, set.seed(1) before each draw.
# Outputs (in <out-dir>/<model>/):
#   <comparison>/volcano_all_genes_<comparison>.svg          4.5 x 5 in; x = log2FoldChange, y = -log10(padj), colour = effect
#   <comparison>/volcano_shrunk_all_genes_<comparison>.svg   4.5 x 5 in; x = shrunkLFC, y = -log10(svalue), colour = shrunk_effect
#   volcano_plots_yscaled.pdf (4 x 5 in, + scale_y_continuous(limits = common ylim)), volcano_plots_xscaled.pdf
#   (4.5 x 5 in, + scale_x_continuous(limits = common xlim)), volcano_plots_freescale.pdf (4 x 5 in): every plot
#   (comparison, then its shrunk plot, ...) one per page
#   volcano_combined_<n>_comparisons_<stacked|side>.svg   comparisons with a combined_order, ordered by it, titles =
#   combined_title (fallback label); stacked 4.5 x 4.5n in, side 4.5n x 5 in; needs >= 2 such comparisons, else skipped
# Exit codes: 0 ok, 3 data problem (missing results file / columns), 1 unexpected error.

.argv <- commandArgs(trailingOnly = FALSE)
.cm <- sub("^--common=", "", grep("^--common=", .argv, value = TRUE)[1])
if (is.na(.cm) || !file.exists(.cm)) .cm <- file.path(dirname(sub("^--file=", "", grep("^--file=", .argv, value = TRUE)[1])), "deg_common.R")
if (!file.exists(.cm)) {
  cat("ERROR: --common=<deg_common.R> not found\n", file = stderr())
  quit(save = "no", status = 3)
}
source(.cm)
suppressPackageStartupMessages(library(patchwork))

opt <- parse_args()
model <- opt_need(opt, "model")
comp_path <- opt_need(opt, "comparisons")
alpha <- as_num(opt_get(opt, "alpha", "0.05"), "alpha")
lfc <- as_num(opt_get(opt, "lfc", "1"), "lfc")
do_shrink <- as_bool(opt_get(opt, "shrink", "TRUE"), "shrink")
sval_threshold <- as_num(opt_get(opt, "sval-threshold", "0.005"), "sval-threshold")
legacy_shrunk_guide <- as_bool(opt_get(opt, "legacy-shrunk-guide", "FALSE"), "legacy-shrunk-guide")
color_down <- opt_get(opt, "color-down", "#46627A")
color_up <- opt_get(opt, "color-up", "#8c6bb1")
color_ns <- opt_get(opt, "color-ns", "gray")
label_genes <- trimws(strsplit(opt_get(opt, "label-genes", ""), ",", fixed = TRUE)[[1]])
label_genes <- label_genes[nzchar(label_genes)]
xlim_global <- parse_limits(opt_get(opt, "xlim", ""), "--xlim")
ylim_global <- parse_limits(opt_get(opt, "ylim", ""), "--ylim")
common_xlim <- parse_limits(opt_get(opt, "common-xlim", ""), "--common-xlim")
common_ylim <- parse_limits(opt_get(opt, "common-ylim", ""), "--common-ylim")
layout <- opt_get(opt, "combined-layout", "stacked")
set_font(opt)
show_legend <- as_bool(opt_get(opt, "legend", "TRUE"), "legend")
legend_pos <- if (show_legend) opt_get(opt, "legend-position", "bottom") else "none"
label_size <- as_num(opt_get(opt, "label-size", "4"), "label-size")
label_nudge_x <- as_num(opt_get(opt, "label-nudge-x", "1.5"), "label-nudge-x")
label_nudge_y <- as_num(opt_get(opt, "label-nudge-y", "1.5"), "label-nudge-y")
label_force <- as_num(opt_get(opt, "label-force", "1"), "label-force")
label_force_pull <- as_num(opt_get(opt, "label-force-pull", "1"), "label-force-pull")
label_fontface <- opt_get(opt, "label-fontface", "italic")
label_in_bounds <- as_bool(opt_get(opt, "label-keep-in-bounds", "TRUE"), "label-keep-in-bounds")
point_size <- as_num(opt_get(opt, "point-size", "2"), "point-size")
ring_stroke <- as_num(opt_get(opt, "ring-stroke", "0.5"), "ring-stroke")
svg_w <- as_num(opt_get(opt, "width", "4.5"), "width")
svg_h <- as_num(opt_get(opt, "height", "5"), "height")
pdf_w <- as_num(opt_get(opt, "pdf-width", "4"), "pdf-width")
pdf_h <- as_num(opt_get(opt, "pdf-height", "5"), "pdf-height")
do_individual <- as_bool(opt_get(opt, "individual", "TRUE"), "individual")
do_pdfs <- as_bool(opt_get(opt, "pdfs", "TRUE"), "pdfs")
do_combined <- as_bool(opt_get(opt, "combined", "TRUE"), "combined")
res_dir <- opt_get(opt, "results-dir", ".")
out_dir <- opt_get(opt, "out-dir", ".")
if (!layout %in% c("stacked", "side")) die("--combined-layout must be stacked or side")

comp <- read_comparisons_tsv(comp_path)
comp <- comp[comp$model == model, , drop = FALSE]
if (nrow(comp) == 0) die(paste0("no comparison rows for model '", model, "' in ", comp_path))
comp$label_resolved <- comparison_label(comp)
mdir <- file.path(out_dir, model)
dir.create(mdir, showWarnings = FALSE, recursive = TRUE)

# ---- read the results tables (staged flat) and add the Label column ----
res <- list()
for (cmp in comp$comparison) {
  f <- file.path(res_dir, paste0("results_", cmp, ".csv"))
  if (!file.exists(f)) die(paste0("results file not found: ", f))
  d <- read.csv(f, colClasses = c(Gene = "character"))
  need <- c("Gene", "log2FoldChange", "padj", "effect", if (do_shrink) c("shrunkLFC", "svalue", "shrunk_effect"))
  if (!all(need %in% names(d))) die(paste0(f, " lacks column(s): ", paste(setdiff(need, names(d)), collapse = ", "), if (do_shrink) " (use --shrink=FALSE for a model without shrinkage)"))
  d$Label <- NA
  d$Label[d$Gene %in% label_genes] <- d$Gene[d$Gene %in% label_genes]
  res[[cmp]] <- d
}
color_values <- c(color_down, color_ns, color_up)
color_labels <- c("Downregulated", "NS", "Upregulated")

# one volcano (raw or shrunk) with the Rmd settings; xlim / ylim / title are the only varying parts
vp_make <- function(df, title, shrunk = FALSE, xlim = NULL, ylim = NULL) {
  makeVolcanoPlot(df,
    lfc.threshold = lfc,
    alpha = if (shrunk && !legacy_shrunk_guide) sval_threshold else alpha,
    color_by = if (shrunk) "shrunk_effect" else "effect",
    label_by = "Label",
    highlight_by = "Label",
    font_color = "black",
    title = if (shrunk) paste0(title, " (Shrunk)") else title,
    xlim = xlim,
    ylim = ylim,
    color_values = color_values,
    color_labels = color_labels,
    x_axis = if (shrunk) "shrunkLFC" else "log2FoldChange",
    y_axis = if (shrunk) "svalue" else "padj",
    x_lab = if (shrunk) "Shrunk Log2 Fold-Change" else "Log2 Fold-Change",
    force = label_force,
    force_pull = label_force_pull,
    label_size = label_size,
    nudgex = label_nudge_x,
    nudgey = label_nudge_y,
    point_size = point_size,
    highlight_stroke = ring_stroke,
    label_fontface = label_fontface,
    keep_labels_in_bounds = label_in_bounds,
    y_lab = if (shrunk) "-Log10(s-value)" else "-Log10(adjusted P-value)",
    leg.pos = legend_pos
  )
}

# ---- individual plots + the list for the pdfs ----
all_volcano_plots <- list()
for (i in seq_len(nrow(comp))) {
  cmp <- comp$comparison[i]
  lab <- comp$label_resolved[i]
  xl <- parse_limits(comp$xlim[i], "xlim"); if (is.null(xl)) xl <- xlim_global
  yl <- parse_limits(comp$ylim[i], "ylim"); if (is.null(yl)) yl <- ylim_global
  cdir <- file.path(mdir, cmp)
  dir.create(cdir, showWarnings = FALSE, recursive = TRUE)

  vp2 <- vp_make(res[[cmp]], lab, FALSE, xl, yl)
  all_volcano_plots[[cmp]] <- vp2
  if (do_individual) draw_svg(vp2, file.path(cdir, paste0("volcano_all_genes_", cmp, ".svg")), svg_w, svg_h)
  if (do_shrink) {
    vp_shrunk <- vp_make(res[[cmp]], lab, TRUE, xl, yl)
    all_volcano_plots[[paste0(cmp, "_shrunk")]] <- vp_shrunk
    if (do_individual) draw_svg(vp_shrunk, file.path(cdir, paste0("volcano_shrunk_all_genes_", cmp, ".svg")), svg_w, svg_h)
  }
}

# ---- multi-page pdfs ----
if (do_pdfs && length(all_volcano_plots) > 0) {
  # default common limits = the commented-out Rmd code
  all_lfc <- unlist(lapply(res, function(d) d$log2FoldChange))
  all_padj <- unlist(lapply(res, function(d) d$padj))
  x_lim <- common_xlim
  if (is.null(x_lim)) x_lim <- range(all_lfc, na.rm = TRUE)
  y_lim <- common_ylim
  if (is.null(y_lim)) {
    minp <- min(all_padj[!is.na(all_padj) & all_padj > 0])
    y_lim <- c(0, -log10(minp))
  }
  write_pdf <- function(name, width, height, add = NULL) {
    open_pdf(file.path(mdir, name), width, height)
    for (comp_name in names(all_volcano_plots)) {
      plot <- all_volcano_plots[[comp_name]]
      if (!is.null(add)) plot <- plot + add
      set.seed(1)
      print(plot)
    }
    invisible(dev.off())
  }
  write_pdf("volcano_plots_yscaled.pdf", pdf_w, pdf_h, scale_y_continuous(limits = y_lim))
  write_pdf("volcano_plots_xscaled.pdf", pdf_w + 0.5, pdf_h, scale_x_continuous(limits = x_lim))
  write_pdf("volcano_plots_freescale.pdf", pdf_w, pdf_h)
}

# ---- combined figure ----
if (do_combined) {
  ord <- suppressWarnings(as.numeric(comp$combined_order))
  sel <- which(!is.na(ord))
  sel <- sel[order(ord[sel])]
  if (length(sel) < 2) {
    message("Combined volcano figure skipped: fewer than 2 comparisons with a combined_order")
  } else {
    stack_plots <- layout == "stacked"
    combined_plots <- list()
    for (k in seq_along(sel)) {
      i <- sel[k]
      cmp <- comp$comparison[i]
      title <- if (is_blank(comp$combined_title[i])) comp$label_resolved[i] else comp$combined_title[i]
      vp <- vp_make(res[[cmp]], title, FALSE, common_xlim, common_ylim)
      # If stacked, remove x-axis title for all but the last plot
      if (stack_plots && k < length(sel)) vp <- vp + theme(axis.title.x = element_blank())
      combined_plots[[cmp]] <- vp
    }
    n <- length(combined_plots)
    if (stack_plots) {
      final_plot <- Reduce(function(a, b) a / b, combined_plots)
      w <- svg_w
      h <- svg_w * n
    } else {
      final_plot <- Reduce(function(a, b) a | b, combined_plots)
      w <- svg_w * n
      h <- svg_h
    }
    final_plot <- final_plot + plot_layout(guides = "collect") & theme(legend.position = legend_pos)
    draw_svg(final_plot, file.path(mdir, paste0("volcano_combined_", n, "_comparisons_", layout, ".svg")), w, h)
  }
}
cat("deg_volcano: model ", model, " done\n", sep = "")
