#!/usr/bin/env Rscript
# deg_ma.R -- MA plots (mean expression vs log2 fold change) of all comparisons of a model.
#
# NEW (not in the legacy Rmd). Same theme and colours as the volcano plots (deg_common.R).
#
# CLI (all --key=value):
#   --common=deg_common.R
#   --model=<name>
#   --comparisons=comparisons.tsv     normalised table; rows of --model
#   --alpha=0.05                      accepted for symmetry with the other scripts (classification is read from the csv)
#   --lfc=1                           dashed horizontal guide lines at +-lfc
#   --color-down=#46627A --color-up=#8c6bb1 --color-ns=gray
#   --svg-device=svglite              svglite (real <text>, default) or cairo (legacy, glyph outlines)
#   --font=Helvetica --legend=TRUE --legend-position=bottom --point-size=2 --width=4.5 --height=5
#   --results-dir=.                   folder with the staged results_<comparison>.csv (default: current directory)
#   --out-dir=.
# Output: <out-dir>/<model>/<comparison>/ma_<comparison>.svg (4.5 x 5 in): x = baseMean (log10 scale), y = log2FoldChange,
#   colour = effect (Downregulated / NS / Upregulated). Genes with baseMean = 0 or NA values cannot be drawn and are dropped.
# Exit codes: 0 ok, 3 data problem, 1 unexpected error.

.argv <- commandArgs(trailingOnly = FALSE)
.cm <- sub("^--common=", "", grep("^--common=", .argv, value = TRUE)[1])
if (is.na(.cm) || !file.exists(.cm)) .cm <- file.path(dirname(sub("^--file=", "", grep("^--file=", .argv, value = TRUE)[1])), "deg_common.R")
if (!file.exists(.cm)) {
  cat("ERROR: --common=<deg_common.R> not found\n", file = stderr())
  quit(save = "no", status = 3)
}
source(.cm)

opt <- parse_args()
model <- opt_need(opt, "model")
comp_path <- opt_need(opt, "comparisons")
lfc <- as_num(opt_get(opt, "lfc", "1"), "lfc")
color_down <- opt_get(opt, "color-down", "#46627A")
color_up <- opt_get(opt, "color-up", "#8c6bb1")
color_ns <- opt_get(opt, "color-ns", "gray")
res_dir <- opt_get(opt, "results-dir", ".")
out_dir <- opt_get(opt, "out-dir", ".")
set_font(opt)
show_legend <- as_bool(opt_get(opt, "legend", "TRUE"), "legend")
legend_pos <- if (show_legend) opt_get(opt, "legend-position", "bottom") else "none"
point_size <- as_num(opt_get(opt, "point-size", "2"), "point-size")
svg_w <- as_num(opt_get(opt, "width", "4.5"), "width")
svg_h <- as_num(opt_get(opt, "height", "5"), "height")

comp <- read_comparisons_tsv(comp_path)
comp <- comp[comp$model == model, , drop = FALSE]
if (nrow(comp) == 0) die(paste0("no comparison rows for model '", model, "' in ", comp_path))
color_map <- c(Downregulated = color_down, NS = color_ns, Upregulated = color_up)

for (i in seq_len(nrow(comp))) {
  cmp <- comp$comparison[i]
  f <- file.path(res_dir, paste0("results_", cmp, ".csv"))
  if (!file.exists(f)) die(paste0("results file not found: ", f))
  d <- read.csv(f, colClasses = c(Gene = "character"))
  if (!all(c("baseMean", "log2FoldChange", "effect") %in% names(d))) die(paste0(f, " lacks baseMean / log2FoldChange / effect"))
  d <- d[!is.na(d$baseMean) & d$baseMean > 0 & !is.na(d$log2FoldChange), , drop = FALSE]
  d$effect <- factor(d$effect, levels = names(color_map))
  # draw NS first so that the significant genes stay on top
  d <- d[order(d$effect != "NS", decreasing = FALSE), , drop = FALSE]

  p <- ggplot(d, aes(x = baseMean, y = log2FoldChange, color = effect)) +
    geom_point(shape = 16, size = point_size) +
    my_theme +
    xlab("Mean of normalized counts") +
    ylab("Log2 Fold-Change") +
    ggtitle(comparison_label(comp[i, , drop = FALSE])) +
    geom_hline(yintercept = c(-lfc, lfc), lty = "dashed", alpha = 0.7) +
    scale_x_log10(labels = scales::label_comma()) +
    theme(legend.title = element_blank(), legend.position = legend_pos) +
    scale_color_manual(values = color_map, drop = FALSE) +
    guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))
  draw_svg(p, file.path(out_dir, model, cmp, paste0("ma_", cmp, ".svg")), svg_w, svg_h)
}
cat("deg_ma: model ", model, " done\n", sep = "")
