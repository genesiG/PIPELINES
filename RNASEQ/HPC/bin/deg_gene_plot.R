#!/usr/bin/env Rscript
# deg_gene_plot.R -- transformed expression of selected genes, one jitter plot per gene.
#
# Generalises the Rmd chunk "Compare rlog counts for specific genes" (the Rmd matched sample names against hard-coded
# patterns to derive the groups; here the groups come from metadata columns).
#
# CLI (all --key=value):
#   --common=deg_common.R
#   --counts=<rlogCounts.rds>         data.frame with a Gene column + one column per sample (deg_transform.R, one model)
#   --metadata=sample_metadata_edited.rds
#   --genes="CDKN2A,EYA4"             comma-separated gene names; genes missing from the table are skipped (message)
#   --x-by=<metadata column>          groups on the x axis (factor level order of the metadata)
#   --color-by=<metadata column>      point colour (empty = same as --x-by)
#   --colors=<csv group,color>        colours keyed by the levels of --color-by; default palette when empty / missing
#   --svg-device=svglite              svglite (real <text>, default) or cairo (legacy, glyph outlines)
#   --font=Helvetica --width=4 --height=4 --point-size=2.5 --alpha=0.75 --jitter-width=0.2 --legend=FALSE
#   --out-dir=.
# Output: <out-dir>/gene_plots/<gene>.svg (4 x 4 in): geom_jitter(width = jitter_width, size = point_size, alpha = point_alpha), theme_classic,
#   y axis "Counts\n(rlog scale)", title = gene name, x labels at 45 degrees, no legend. set.seed(1) before each draw.
# Exit codes: 0 ok (also when no requested gene is found), 3 data problem (unknown column), 1 unexpected error.

.argv <- commandArgs(trailingOnly = FALSE)
.cm <- sub("^--common=", "", grep("^--common=", .argv, value = TRUE)[1])
if (is.na(.cm) || !file.exists(.cm)) .cm <- file.path(dirname(sub("^--file=", "", grep("^--file=", .argv, value = TRUE)[1])), "deg_common.R")
if (!file.exists(.cm)) {
  cat("ERROR: --common=<deg_common.R> not found\n", file = stderr())
  quit(save = "no", status = 3)
}
source(.cm)

opt <- parse_args()
counts_path <- opt_need(opt, "counts")
meta_path <- opt_need(opt, "metadata")
genes <- trimws(strsplit(opt_need(opt, "genes"), ",", fixed = TRUE)[[1]])
genes <- unique(genes[nzchar(genes)])
x_by <- opt_need(opt, "x-by")
color_by <- opt_get(opt, "color-by", "")
if (!nzchar(color_by)) color_by <- x_by
colors_csv <- opt_get(opt, "colors", "")
out_dir <- opt_get(opt, "out-dir", ".")
set_font(opt)
svg_w <- as_num(opt_get(opt, "width", "4"), "width")
svg_h <- as_num(opt_get(opt, "height", "4"), "height")
point_size <- as_num(opt_get(opt, "point-size", "2.5"), "point-size")
point_alpha <- as_num(opt_get(opt, "alpha", "0.75"), "alpha")
jitter_width <- as_num(opt_get(opt, "jitter-width", "0.2"), "jitter-width")
legend_pos <- if (as_bool(opt_get(opt, "legend", "FALSE"), "legend")) "right" else "none"

for (p in c(counts_path, meta_path)) if (!file.exists(p)) die(paste0("file not found: ", p))
rlogCounts <- readRDS(counts_path)
sample_metadata <- readRDS(meta_path)
if (!"Gene" %in% names(rlogCounts)) die("--counts table needs a Gene column")
for (cn in c(x_by, color_by)) if (!cn %in% names(sample_metadata)) die(paste0("metadata has no column '", cn, "'"))
samples <- setdiff(names(rlogCounts), "Gene")
if (!all(samples %in% rownames(sample_metadata))) die("sample columns of --counts are not all in the metadata")

x_fac <- droplevels(as.factor(sample_metadata[samples, x_by]))
col_fac <- droplevels(as.factor(sample_metadata[samples, color_by]))
color_mapping <- group_colors(levels(col_fac), colors_csv)

found <- 0
for (g in genes) {
  row <- rlogCounts[rlogCounts$Gene == g, samples, drop = FALSE]
  if (nrow(row) == 0) {
    message("Gene ", g, " not found in the table; skipped")
    next
  }
  found <- found + 1
  df <- data.frame(
    Sample = samples,
    RlogCounts = as.numeric(row[1, ]),
    X = x_fac,
    Color = col_fac
  )
  b <- ggplot(df, aes(x = X, y = RlogCounts, color = Color)) +
    geom_jitter(width = 0.2, size = 2.5, alpha = 0.75) +
    scale_color_manual(values = color_mapping) +
    theme_classic() +
    labs(
      title = g,
      y = "Counts\n(rlog scale)"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      legend.position = legend_pos,
      title = element_text(hjust = 0.5)
    )
  draw_svg(b, file.path(out_dir, "gene_plots", paste0(gsub("[^A-Za-z0-9._+-]", "_", g), ".svg")), svg_w, svg_h)
}
if (found == 0) message("WARNING: none of the requested genes is in the table; no plot written")
cat("deg_gene_plot: ", found, " of ", length(genes), " gene(s) plotted\n", sep = "")
