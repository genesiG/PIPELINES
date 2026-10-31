#!/usr/bin/env Rscript
# deg_heatmap.R -- heatmap of the differentially expressed genes of one comparison (transformed counts, row z-score).
#
# Legacy mapping (step_6_deg_analysis.Rmd): chunk "Gene Expression Heatmap" (treatment model, cutree_rows = 3, no branch
#   swap); the reverse_clusters block of that chunk is --col-swap=scaled.
#
# CLI (all --key=value):
#   --common=deg_common.R
#   --name=<heatmap name>             output file name part
#   --model=<name>                    output sub-folder
#   --comparison=<name>               comparison whose DEGs are shown (only used in messages)
#   --counts=<rlogCounts.rds>         data.frame with a Gene column + one column per sample (deg_transform.R)
#   --results=results_<comparison>.csv   results table of the comparison (deg_results.R)
#   --gene-set=<file>                 optional text file, one gene per line: only DEGs that are in the set are drawn
#   --label-rules=<csv>               optional table pattern,replacement applied in order to the column names
#                                     (typed \n in the replacement = newline); empty = column names as they are
#   --clust-method=ward.D2            hclust method of rows and columns
#   --cutree-rows=<int>               number of row clusters to cut (empty = none)
#   --col-swap=none|scaled|unscaled   none: pheatmap clusters the columns itself; scaled / unscaled: hclust(dist(t(M))) on
#                                     the row-scaled / the unscaled matrix, the two top branches are swapped, and the result
#                                     is passed as cluster_cols (legacy "reverse_clusters")
#   --palette-n=1000 --fontsize=9 --fontsize-row=6.5 --width=6.5 --height=4.5
#   --title="DEGs (N = {n})"          {n} is replaced by the number of genes drawn
#   --svg-device=svglite              svglite (real <text>, default) or cairo (legacy, glyph outlines)
#   --font=Helvetica --show-rownames=FALSE --angle-col=90 --cluster-rows=TRUE --treeheight-row= --treeheight-col=
#                                     (empty treeheight = pheatmap default)
#   --color-down=#46627A --color-mid=white --color-up=#88419d    (legacy palette: #4d004b - white - goldenrod)
#   --out-dir=.
# Genes: effect != "NS" in the results table (and in --gene-set when given), in the row order of the transformed counts
#   table (NOT the results order, as the Rmd). pheatmap(scale = "row", clustering_distance_rows = "euclidean",
#   show_rownames = show_rownames, angle_col = 90).
# Output: <out-dir>/<model>/heatmap_<name>.svg. With fewer than 2 genes a message is printed, nothing is written, exit 0.
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
  library(pheatmap)
  library(dplyr)
  library(tibble)
})

opt <- parse_args()
name <- opt_need(opt, "name")
model <- opt_need(opt, "model")
comparison <- opt_get(opt, "comparison", "")
counts_path <- opt_need(opt, "counts")
results_path <- opt_need(opt, "results")
gene_set_path <- opt_get(opt, "gene-set", "")
rules_path <- opt_get(opt, "label-rules", "")
clust_method <- opt_get(opt, "clust-method", "ward.D2")
cutree_arg <- opt_get(opt, "cutree-rows", "")
col_swap <- opt_get(opt, "col-swap", "none")
palette_n <- as_num(opt_get(opt, "palette-n", "1000"), "palette-n")
fontsize <- as_num(opt_get(opt, "fontsize", "9"), "fontsize")
fontsize_row <- as_num(opt_get(opt, "fontsize-row", "6.5"), "fontsize-row")
width <- as_num(opt_get(opt, "width", "6.5"), "width")
height <- as_num(opt_get(opt, "height", "4.5"), "height")
title <- opt_get(opt, "title", "DEGs (N = {n})")
color_down <- opt_get(opt, "color-down", "#46627A")
color_mid <- opt_get(opt, "color-mid", "white")
color_up <- opt_get(opt, "color-up", "#88419d")
out_dir <- opt_get(opt, "out-dir", ".")
set_font(opt)
show_rownames <- as_bool(opt_get(opt, "show-rownames", "FALSE"), "show-rownames")
angle_col <- as_num(opt_get(opt, "angle-col", "90"), "angle-col")
cluster_rows <- as_bool(opt_get(opt, "cluster-rows", "TRUE"), "cluster-rows")
th_row <- opt_get(opt, "treeheight-row", "")
th_col <- opt_get(opt, "treeheight-col", "")
if (!col_swap %in% c("none", "scaled", "unscaled")) die("--col-swap must be none, scaled or unscaled")
cutree_rows <- if (nzchar(trimws(cutree_arg))) as_num(cutree_arg, "cutree-rows") else NA

for (p in c(counts_path, results_path)) if (!file.exists(p)) die(paste0("file not found: ", p))

# Load table of DEG
results_df <- read.csv(results_path, header = TRUE, colClasses = c(Gene = "character"))
if (!all(c("Gene", "effect") %in% names(results_df))) die("results table needs the columns Gene and effect")
results_df <- filter(results_df, effect != "NS")
deg_genes <- results_df$Gene

# Subset the transformed count table by DEG (row order of the counts table)
rlogCounts <- readRDS(counts_path)
if (!"Gene" %in% names(rlogCounts)) die("--counts table needs a Gene column (rlogCounts.rds from deg_transform.R)")
keep <- rlogCounts$Gene %in% deg_genes
if (nzchar(gene_set_path)) {
  if (!file.exists(gene_set_path)) die(paste0("gene set file not found: ", gene_set_path))
  gs <- trimws(readLines(gene_set_path, warn = FALSE))
  gs <- gs[nzchar(gs)]
  keep <- keep & rlogCounts$Gene %in% gs
}
rlogCounts <- rlogCounts[keep, , drop = FALSE]
rownames(rlogCounts) <- NULL

# Convert rlogCounts to a numeric matrix
rlogCounts <- rlogCounts %>% column_to_rownames(var = "Gene")
rlogCounts_mat <- as.matrix(rlogCounts)

# Get the number of rows (genes)
n_all <- nrow(rlogCounts_mat)
if (n_all < 2) {
  message("Heatmap '", name, "' (", comparison, "): ", n_all, " gene(s) to draw (need >= 2); nothing written")
  quit(save = "no", status = 0)
}
if (!is.na(cutree_rows) && cutree_rows > n_all) {
  message("NOTE: --cutree-rows=", cutree_rows, " > ", n_all, " genes; using ", n_all)
  cutree_rows <- n_all
}

labels_col <- colnames(rlogCounts)
if (nzchar(rules_path)) labels_col <- apply_label_rules(labels_col, read_csv_clean(rules_path, trim = FALSE))

if (col_swap != "none") {
  # scaled: cluster on the row-scaled data, which is what pheatmap does by default with scale = "row"
  mat_for_clust <- if (col_swap == "scaled") t(scale(t(rlogCounts_mat))) else rlogCounts_mat
  # Swap top-level dendrogram branches
  col_dend_raw <- as.dendrogram(hclust(dist(t(mat_for_clust)), method = clust_method))
  col_dend_swap <- col_dend_raw
  col_dend_swap[[1]] <- col_dend_raw[[2]]
  col_dend_swap[[2]] <- col_dend_raw[[1]]
  col_dend <- as.hclust(col_dend_swap)
} else {
  # Don't swap the dendrogram branches
  col_dend <- TRUE
}

hm_dir <- file.path(out_dir, model)
dir.create(hm_dir, showWarnings = FALSE, recursive = TRUE)
open_svg(file.path(hm_dir, paste0("heatmap_", name, ".svg")), width, height)
hm <- pheatmap(rlogCounts_mat,
  color = colorRampPalette(c(color_down, color_mid, color_up))(palette_n),
  scale = "row",
  clustering_distance_rows = "euclidean",
  cluster_rows = cluster_rows,
  cluster_cols = col_dend,
  clustering_method = clust_method,
  show_rownames = FALSE,
  main = gsub("{n}", n_all, title, fixed = TRUE),
  fontsize_row = fontsize_row,
  fontsize = fontsize,
  labels_col = labels_col,
  cutree_rows = cutree_rows,
  angle_col = angle_col,
  treeheight_row = if (nzchar(th_row)) as_num(th_row, "treeheight-row") else if (cluster_rows) 50 else 0,
  treeheight_col = if (nzchar(th_col)) as_num(th_col, "treeheight-col") else 50,
  silent = TRUE
)
grid::grid.newpage()
grid::grid.draw(hm$gtable)
close_svg()
cat("deg_heatmap: ", name, " (", n_all, " genes) written\n", sep = "")
