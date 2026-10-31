#!/usr/bin/env Rscript
# deg_qc.R -- sample-level QC of the prepared counts: reads per sample, MDS plot, sample correlation heatmap.
#
# Legacy mapping (step_6_deg_analysis.Rmd): chunks "QC - Check counts per sample" (counts_per_index.pdf),
#   "MDS plots (by genotypes)" (pca.pdf) and "Sample Correlation Heatmap" (spearman_correlation_heatmap.svg).
#
# CLI (all --key=value; booleans TRUE/FALSE):
#   --common=deg_common.R
#   --counts=counts_edited.rds        data.frame genes x samples (from deg_prepare.R)
#   --metadata=sample_metadata_edited.rds
#   --comparisons=comparisons.tsv     only used to find the default --color-by
#   --color-by=<metadata column>      group / colour column; empty = first design variable of the first model
#   --colors=<csv group,color>        colours keyed by the levels of --color-by; empty / missing groups = default
#                                     palette (Set2, or hue palette above 8 groups)
#   --barplot=TRUE --mds=TRUE --correlation=TRUE   switch each plot
#   --mds-input=logcpm|counts         logcpm (default) = edgeR::cpm(log = TRUE, prior.count = 2); counts = raw count matrix
#                                     as the legacy (dominated by library size and the most highly expressed genes)
#   --mds-label-by=<column>           text drawn at each point; empty = --color-by
#   --mds-legend-title=<text>         legend title; empty = the --color-by column name
#   --legacy-mds-legend-col=<empty|column>   reproduce the LEGACY legend (see below)
#   --legacy-color-index=FALSE        TRUE reproduces the LEGACY per-sample colours (see below)
#   --cor-method=spearman             method of cor() for the correlation heatmap
#   --svg-device=svglite              svglite (real <text>, default) or cairo (legacy, glyph outlines)
#   --font=Helvetica --barplot-title="Reads per sample" --mds-cex=0.75 --mds-legend-pos=topleft
#   --cor-numbers=TRUE --cor-width=8.5 --cor-height=6
#   --out-dir=.
# Outputs (in --out-dir):
#   counts_per_index.pdf  ggplot bar plot (theme_classic, pdf 4 x 5 in, title "Reads per sample"), one bar per sample in
#                         metadata order, filled with the group colour, black outline
#   counts_per_sample.tsv sample, total_counts, group
#   pca.pdf               limma::plotMDS(gene.selection = "common", cex = mds_cex), pdf 5 x 5 in, legend top left. The
#                         file is called pca.pdf as in the legacy project but it is an MDS plot, not a PCA.
#   <cor-method>_correlation_heatmap.svg   (spearman_correlation_heatmap.svg by default) pheatmap of
#                         cor(counts, method = cor-method): palette #313695-#ffffbf-#a50026 (200), group annotation on
#                         rows and columns, euclidean distance on the correlations, complete linkage, numbers
#                         shown; svg 8.5 x 6 in (width/height kept as written in the Rmd)
# LEGACY BUGS (defaults are the corrected behaviour, each can be switched back):
#   1. MDS legend: the legacy legend LABELS were unique(gRNA) while the colours came from the treatment column, so
#      labels and colours did not match. Default: legend = levels of --color-by with their own colours. With
#      --legacy-mds-legend-col=<col> the legacy legend is reproduced (labels = gsub("_|-", " ", unique(col)),
#      fill = unique(sample colours), title = --mds-legend-title or the column name).
#   2. Colour indexing: the Rmd did treatment_colors[sample_metadata$treatment], i.e. it indexed the colour vector
#      with the INTEGER CODES of the (re-levelled) factor instead of its names, so bars / MDS points of groups whose
#      factor level position differs from their position in the colour vector got another group's colour (in the GG
#      project CBX2_1uM and LUC_1uM were swapped). Default: colours by group name. --legacy-color-index=TRUE
#      reproduces the legacy indexing (colour vector in csv order indexed by factor codes). The correlation heatmap
#      annotation always used names (correct).
# Exit codes: 0 ok, 3 data problem, 1 unexpected error.

.argv <- commandArgs(trailingOnly = FALSE)
.cm <- sub("^--common=", "", grep("^--common=", .argv, value = TRUE)[1])
if (is.na(.cm) || !file.exists(.cm)) .cm <- file.path(dirname(sub("^--file=", "", grep("^--file=", .argv, value = TRUE)[1])), "deg_common.R")
if (!file.exists(.cm)) {
  cat("ERROR: --common=<deg_common.R> not found\n", file = stderr())
  quit(save = "no", status = 3)
}
source(.cm)
suppressPackageStartupMessages({
  library(limma)
  library(pheatmap)
})

opt <- parse_args()
counts_path <- opt_need(opt, "counts")
meta_path <- opt_need(opt, "metadata")
comp_path <- opt_need(opt, "comparisons")
color_by <- opt_get(opt, "color-by", "")
colors_csv <- opt_get(opt, "colors", "")
do_bar <- as_bool(opt_get(opt, "barplot", "TRUE"), "barplot")
do_mds <- as_bool(opt_get(opt, "mds", "TRUE"), "mds")
do_cor <- as_bool(opt_get(opt, "correlation", "TRUE"), "correlation")
mds_input <- opt_get(opt, "mds-input", "logcpm")
mds_label_by <- opt_get(opt, "mds-label-by", "")
mds_legend_title <- opt_get(opt, "mds-legend-title", "")
legacy_legend_col <- opt_get(opt, "legacy-mds-legend-col", "")
legacy_color_index <- as_bool(opt_get(opt, "legacy-color-index", "FALSE"), "legacy-color-index")
cor_method <- opt_get(opt, "cor-method", "spearman")
out_dir <- opt_get(opt, "out-dir", ".")
set_font(opt)
bar_title <- opt_get(opt, "barplot-title", "Reads per sample")
mds_cex <- as_num(opt_get(opt, "mds-cex", "0.75"), "mds-cex")
mds_legend_pos <- opt_get(opt, "mds-legend-pos", "topleft")
cor_numbers <- as_bool(opt_get(opt, "cor-numbers", "TRUE"), "cor-numbers")
cor_w <- as_num(opt_get(opt, "cor-width", "8.5"), "cor-width")
cor_h <- as_num(opt_get(opt, "cor-height", "6"), "cor-height")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
if (!mds_input %in% c("counts", "logcpm")) die("--mds-input must be counts or logcpm")
if (!cor_method %in% c("spearman", "pearson", "kendall")) die("--cor-method must be spearman, pearson or kendall")

for (p in c(counts_path, meta_path)) if (!file.exists(p)) die(paste0("file not found: ", p))
counts <- as.data.frame(readRDS(counts_path))
sample_metadata <- readRDS(meta_path)
if (!identical(colnames(counts), rownames(sample_metadata))) die("counts columns and metadata rows differ; run deg_prepare.R first")

comp <- read_comparisons_tsv(comp_path)
if (!nzchar(color_by)) color_by <- all.vars(as.formula(comp$design[1]))[1]
if (!color_by %in% names(sample_metadata)) die(paste0("--color-by column '", color_by, "' not in the metadata"))
if (!nzchar(mds_label_by)) mds_label_by <- color_by
if (!mds_label_by %in% names(sample_metadata)) die(paste0("--mds-label-by column '", mds_label_by, "' not in the metadata"))
if (nzchar(legacy_legend_col) && !legacy_legend_col %in% names(sample_metadata)) die(paste0("--legacy-mds-legend-col column '", legacy_legend_col, "' not in the metadata"))
if (!nzchar(mds_legend_title)) mds_legend_title <- if (nzchar(legacy_legend_col)) legacy_legend_col else color_by

grp <- as.factor(sample_metadata[[color_by]])
grp <- droplevels(grp)
grp_cols <- group_colors(levels(grp), colors_csv) # named by level
if (legacy_color_index) {
  # legacy: colour vector (csv order) indexed with the factor's integer codes
  legacy_vec <- read_colors_csv(colors_csv)
  if (is.null(legacy_vec)) legacy_vec <- grp_cols
  cols <- unname(legacy_vec[as.integer(grp)])
  message("NOTE: --legacy-color-index=TRUE: sample colours taken by factor code, as in the Rmd")
} else {
  cols <- unname(grp_cols[as.character(grp)])
}

# ------------------------------------------------------------------ counts per sample ----
df_counts <- data.frame(
  Sample = colnames(counts),
  Counts = colSums(counts)
)
write.table(data.frame(sample = df_counts$Sample, total_counts = df_counts$Counts, group = as.character(grp)),
  file.path(out_dir, "counts_per_sample.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (do_bar) {
  # Ensure Sample is a factor with levels in the desired order
  df_counts$Sample <- factor(df_counts$Sample, levels = df_counts$Sample)

  bp <- ggplot(df_counts, aes(x = Sample, y = Counts, fill = Sample)) +
    geom_bar(stat = "identity", show.legend = FALSE, color = "black") +
    scale_fill_manual(values = cols) +
    labs(title = bar_title) +
    xlab("sAMPLES") +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, size = 11),
      axis.text.y = element_text(size = 11),
      axis.title.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  open_pdf(file.path(out_dir, "counts_per_index.pdf"), 4, 5)
  print(bp)
  invisible(dev.off())
}

# -------------------------------------------------------------------------------- MDS ----
if (do_mds) {
  counts_df <- counts
  if (mds_input == "counts") {
    message(
      "NOTE: MDS on raw counts: distances are dominated by library size and a few highly expressed genes. Use --mds-input=logcpm to compare samples."
    )
  }
  if (mds_input == "logcpm") {
    counts_df <- as.data.frame(edgeR::cpm(as.matrix(counts), log = TRUE, prior.count = 2))
  }
  labels1 <- as.character(sample_metadata[[mds_label_by]])

  open_pdf(file.path(out_dir, "pca.pdf"), 5, 5)
  plotMDS(counts_df,
    labels = labels1,
    main = "",
    col = cols, # colors
    gene.selection = "common", cex = 0.75
  )
  if (nzchar(legacy_legend_col)) {
    # legacy legend: labels from another column than the colours (kept on request)
    leg1 <- unique(sample_metadata[[legacy_legend_col]])
    leg1_sub <- gsub("_|-", " ", leg1)
    legend(mds_legend_pos,
      fill = unique(cols),
      legend = leg1_sub,
      title = mds_legend_title
    )
  } else {
    legend(mds_legend_pos,
      fill = unname(grp_cols[levels(grp)]),
      legend = levels(grp),
      title = mds_legend_title
    )
  }
  invisible(dev.off())
}

# --------------------------------------------------------------- correlation heatmap ----
if (do_cor) {
  cor_matrix <- cor(counts, method = cor_method)

  # sample annotations (column name = capitalised --color-by, as the legacy "Treatment")
  ann_name <- paste0(toupper(substr(color_by, 1, 1)), substring(color_by, 2))
  annotation_df <- data.frame(grp, row.names = rownames(sample_metadata))
  names(annotation_df) <- ann_name
  # legend order = order of the colour vector: colours csv order first (as the Rmd), then groups missing from the csv
  ann_cols <- grp_cols
  given <- read_colors_csv(colors_csv)
  if (!is.null(given)) ann_cols <- c(given[names(given) %in% names(grp_cols)], grp_cols[setdiff(names(grp_cols), names(given))])
  annotation_colors <- setNames(list(ann_cols), ann_name)

  hm_path <- file.path(out_dir, paste0(cor_method, "_correlation_heatmap.svg"))
  open_svg(hm_path, cor_w, cor_h)
  cor_hm <- pheatmap(
    cor_matrix,
    main = paste0("Sample Correlation (", paste0(toupper(substr(cor_method, 1, 1)), substring(cor_method, 2)), ")"),
    color = colorRampPalette(c("#313695", "#ffffbf", "#a50026"))(200),
    annotation_col = annotation_df,
    annotation_row = annotation_df,
    annotation_colors = annotation_colors,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    display_numbers = cor_numbers,
    number_color = "black",
    fontsize_number = 8,
    silent = TRUE
  )
  grid::grid.newpage()
  grid::grid.draw(cor_hm$gtable)
  close_svg()
}
cat("deg_qc: done (color-by = ", color_by, ")\n", sep = "")
