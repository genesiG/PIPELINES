#!/usr/bin/env Rscript
#' ##########################################################################
#' guide_representation_compare_report.R
#'
#' One HTML report comparing every library (pool) in this run side by side:
#' the distributional plots faceted one panel per pool (3 columns), plus the
#' headline scalar metrics (Gini, 90/10 ratio, % undetected, mapping rate) as
#' combined bar charts across pools (one full-width chart per metric, one facet
#' per pool). Also: an ECDF figure with one block per pool (curves on top, mean
#' read count per guide bars below), a PCA of samples per pool, and the total
#' reads per sample. Every section is titled with the question it answers.
#' Built only from the tables
#' guide_representation_wrangle.R already wrote for each library — it never
#' reads a count table directly or recomputes anything.
#'
#' Only useful with more than one library; main.nf only runs this step when
#' the run has more than one.
#'
#' If --tables-dir also holds <project>_library_breakdown.csv
#' (bin/library_breakdown.py), the report opens with a "Reads per Library"
#' section: how each sample's reads split across the libraries it was counted
#' against, plus reads matching no library.
#'
#' USAGE:
#'   Rscript guide_representation_compare_report.R --tables-dir=<dir> \
#'       --project=<name> [--theme=<theme.R>] [--outdir=<dir>]
#'
#'   --tables-dir  directory holding every library's <prefix>_*.csv from the
#'                 wrangle step (one prefix per library)
#'   --project     project_name (used to strip the common prefix off each
#'                 library's name for compact facet labels, and shown in the
#'                 report header)
#'   --theme       R file defining a ggplot2 theme object named `my_theme`
#'                 (default: my_theme.R next to this script)
#'   --outdir      where to write the report (default: current directory)
#'
#' OUTPUT: <project>_pools_comparison_qc_report.html (self-contained — every
#' plot is embedded inline as SVG). Besides the HTML, each plot is also saved
#' as SVG + PNG under <outdir>/pools_comparison/ (<project>_pools_<plot>.svg/.png,
#' e.g. _pools_bar_Gini, _pools_ecdf, _pools_pca, _pools_read_depth).
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
if (is.null(opt[["tables-dir"]]) || is.null(opt$project)) {
  stop("Usage: guide_representation_compare_report.R --tables-dir=<dir> --project=<name> [--theme=<theme.R>] [--outdir=<dir>]")
}
tables_dir <- opt[["tables-dir"]]
project    <- opt$project
outdir     <- if (is.null(opt$outdir)) "." else opt$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

theme_file <- opt$theme
if (is.null(theme_file)) {
  self <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  theme_file <- file.path(dirname(normalizePath(self)), "my_theme.R")
}
source(theme_file) # defines `my_theme`

# ---- Helpers ----------------------------------------------------------------
html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x
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

# Same fix as guide_representation_report.R: R's svg() device restarts its
# internal <symbol>/<clipPath> id counters (glyph0-0, clip1, ...) in every
# file, so embedding several in one page makes them collide and render with
# each other's glyph shapes. Namespacing every id fixes it.
namespace_svg_ids <- function(svg_text, ns) {
  ids <- unique(regmatches(svg_text, gregexpr('(?<=id=")[^"]+', svg_text, perl = TRUE))[[1]])
  for (id in ids) {
    new_id <- paste0(ns, "-", id)
    svg_text <- gsub(paste0('id="', id, '"'),    paste0('id="', new_id, '"'),    svg_text, fixed = TRUE)
    svg_text <- gsub(paste0('href="#', id, '"'), paste0('href="#', new_id, '"'), svg_text, fixed = TRUE)
    svg_text <- gsub(paste0('url(#', id, ')'),   paste0('url(#', new_id, ')'),   svg_text, fixed = TRUE)
  }
  svg_text
}

# Renders a ggplot (or grob) to SVG and returns it as a namespaced, embeddable string.
plot_to_svg <- function(p, ns, width, height) {
  tmp <- tempfile(fileext = ".svg")
  grDevices::svg(filename = tmp, width = width, height = height)
  draw_plot(p)
  invisible(grDevices::dev.off())
  raw <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
  plots_dir <- file.path(outdir, "pools_comparison")
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  base <- paste0(project, "_pools_", gsub("[^A-Za-z0-9]+", "_", sub("^cmp-", "", ns)))
  file.copy(tmp, file.path(plots_dir, paste0(base, ".svg")), overwrite = TRUE)
  save_png(p, file.path(plots_dir, paste0(base, ".png")), width, height)
  unlink(tmp)
  start <- regexpr("<svg[ >]", raw)
  end   <- regexpr("</svg>", raw)
  svg   <- substr(raw, start, end + attr(end, "match.length") - 1)
  namespace_svg_ids(svg, ns)
}

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

fmt_int <- function(x) ifelse(is.na(x), "—", format(round(x), big.mark = ",", scientific = FALSE))
fmt_pct <- function(x, digits = 1) ifelse(is.na(x), "—", paste0(round(x, digits), "%"))
fmt_num <- function(x, digits = 3) ifelse(is.na(x), "—", format(round(x, digits), nsmall = digits))
badge   <- function(text, kind) sprintf('<span class="badge badge-%s">%s</span>', kind, html_escape(text))

# ---- Discover libraries and load their tables ------------------------------
prefixes <- sub("_qc_summary\\.csv$", "", list.files(tables_dir, pattern = "_qc_summary\\.csv$"))
if (length(prefixes) < 2) {
  stop("Found ", length(prefixes), " library qc_summary.csv file(s) in ", tables_dir,
       " — need at least two to build a cross-pool comparison report.")
}
# Compact facet/axis label: strip the shared "<project>_" prefix, e.g.
# "AM_POOL_2026_pool_1" -> "pool_1".
pool_of <- function(px) {
  sub(paste0("^", project, "_"), "", px)
}
pools <- vapply(prefixes, pool_of, character(1))
n_pools     <- length(pools)
facet_ncol  <- min(3, n_pools)
facet_nrow  <- ceiling(n_pools / facet_ncol)
facet_w     <- 3.3 * facet_ncol + 1
facet_h     <- 2.9 * facet_nrow + 1

read_all <- function(name, required = TRUE) {
  parts <- lapply(seq_along(prefixes), function(i) {
    path <- file.path(tables_dir, paste0(prefixes[i], "_", name, ".csv"))
    if (!file.exists(path)) return(NULL)
    df <- read.csv(path, check.names = FALSE)
    df$pool <- pools[i]
    df
  })
  parts <- Filter(Negate(is.null), parts)
  if (length(parts) == 0) {
    if (required) stop("No ", name, ".csv found for any library in ", tables_dir)
    return(NULL)
  }
  dplyr::bind_rows(parts)
}

qc_all     <- read_all("qc_summary")
counts_all <- read_all("counts_long")
ranked_all <- read_all("ranked")
lorenz_all <- read_all("lorenz")
gc_all     <- read_all("gc_vs_depth", required = FALSE)
class_all  <- read_all("read_classification", required = FALSE)
pca_all    <- read_all("pca", required = FALSE)
pcavar_all <- read_all("pca_variance", required = FALSE)
depth_all  <- read_all("read_depth", required = FALSE)

# One colour per sample, shared by the ECDF and PCA figures (same sample, same
# colour in every pool).
all_samples <- sort(unique(counts_all$sample))
sample_cols <- setNames(scales::hue_pal()(length(all_samples)), all_samples)

breakdown_path <- file.path(tables_dir, paste0(project, "_library_breakdown.csv"))
breakdown_all  <- if (file.exists(breakdown_path)) read.csv(breakdown_path, check.names = FALSE) else NULL
if (!is.null(breakdown_all) && nrow(breakdown_all) == 0) breakdown_all <- NULL

# All-vs-all mode (no samplesheet `library` column): every sample is counted
# against every library, so samples have several real library rows. Detected
# from the data: the breakdown file exists and at least one sample has more
# than one real library row (assigned mode: exactly one per sample).
breakdown_real <- if (!is.null(breakdown_all)) breakdown_all[!(breakdown_all$library %in% c("Any library", "No library guide")), ] else NULL
all_vs_all <- !is.null(breakdown_real) && any(table(breakdown_real$sample) > 1)

has_mapping <- "mapping_pct" %in% names(qc_all)
has_rpg     <- "reads_per_guide" %in% names(qc_all)

# One x-axis label per (pool, sample): just the pool name when it has a
# single sample (the common case), "pool: sample" when it has several.
# Ordered by pool then sample, so the bar charts below group each pool's
# bar(s) together instead of sorting alphabetically by label.
qc_all <- qc_all |>
  group_by(pool) |>
  mutate(n_in_pool = n_distinct(sample)) |>
  ungroup() |>
  mutate(x_label = ifelse(n_in_pool > 1, paste0(pool, ": ", sample), pool)) |>
  arrange(pool, sample) |>
  mutate(x_label = factor(x_label, levels = unique(x_label)))

# Low mapping and the "Mainly <pool>" note, one value per qc_all row. Computed
# after the sort above so each value lines up with its row.
low_mapping <- if (has_mapping) qc_all$mapping_pct < 30 else rep(FALSE, nrow(qc_all))
mainly_pool <- rep(NA_character_, nrow(qc_all))
if (all_vs_all) {
  # The sample's "Any library" % of reads: below 30 means it matches no pool
  # well. Falls back to the mapping_pct rule for a sample missing from the breakdown.
  any_rows    <- breakdown_all[breakdown_all$library == "Any library", ]
  any_pct     <- setNames(suppressWarnings(as.numeric(any_rows$pct_of_sample_reads)), any_rows$sample)
  any_row_pct <- any_pct[qc_all$sample]
  low_mapping <- ifelse(is.na(any_row_pct), low_mapping, any_row_pct < 30)
  # The sample's main library: the real library with the most mapped reads
  # (same rule as the breakdown table below), named like qc_all$pool.
  main_pool_of <- function(smp) {
    pl <- breakdown_real[breakdown_real$sample == smp, ]
    mr <- suppressWarnings(as.numeric(pl$mapped_reads))
    if (nrow(pl) == 0 || all(is.na(mr)) || max(mr, na.rm = TRUE) <= 0) return(NA_character_)
    pool_of(paste0(project, "_", pl$library[which.max(mr)]))
  }
  mainly_pool <- vapply(qc_all$sample, main_pool_of, character(1))
}

# Plain-language notes shown only in all-vs-all mode.
aa_note_glance    <- "All-vs-all mode: each Mapping % is the share of that sample's reads that match THIS pool, so low values for samples of other pools are expected. Low mapping is only flagged when a sample matches no pool well (Any library below 30%); a sample whose reads mostly match another pool is marked &ldquo;Mainly &lt;pool&gt;&rdquo;."
aa_note_map       <- "All-vs-all mode: a pool's mapping rate bar is the share of each sample's reads that match THAT pool, so low bars for samples that belong to other pools are expected. The Low mapping flag is only raised when a sample matches no pool well (Any library below 30%); a sample whose reads mostly match another pool is labelled &ldquo;Mainly &lt;pool&gt;&rdquo; instead."
aa_note_breakdown <- "All-vs-all mode: each sample is counted against every pool, so in the per-pool panels below a sample's low mapping or dropout in a pool it does not belong to is expected. Its main library here shows the pool it belongs to. Low mapping is only flagged when a sample matches no pool well (Any library below 30%)."

# ---- "At a glance" table: every pool, every sample, one row each ----------
flags_col <- vapply(seq_len(nrow(qc_all)), function(i) {
  b <- character(0)
  if (isTRUE(qc_all$library_skew_flag[i]))    b <- c(b, badge("Library skew", "flag"))
  if (isTRUE(qc_all$undetected_flag[i]))      b <- c(b, badge("High dropout", "flag"))
  if (has_mapping && isTRUE(low_mapping[i]))  b <- c(b, badge("Low mapping",  "flag"))
  if (all_vs_all && !isTRUE(low_mapping[i]) && !is.na(mainly_pool[i]) && mainly_pool[i] != qc_all$pool[i])
    b <- c(b, badge(paste0("Mainly ", mainly_pool[i]), "info"))
  if (length(b) == 0) b <- badge("OK", "ok")
  paste(b, collapse = " ")
}, character(1))

summary_rows <- paste0(
  "<tr>",
  "<td>", html_escape(qc_all$pool), "</td>",
  "<td>", html_escape(qc_all$sample), "</td>",
  if (has_mapping) paste0("<td class='num'>", fmt_pct(qc_all$mapping_pct), "</td>") else "",
  if (has_rpg) paste0("<td class='num'>", fmt_num(qc_all$reads_per_guide, 1), "</td>") else "",
  "<td class='num'>", fmt_num(qc_all$Gini, 4), "</td>",
  "<td class='num'>", fmt_num(qc_all$ratio_90_10, 2), "</td>",
  "<td class='num'>", fmt_pct(qc_all$pct_undetected, 2), "</td>",
  "<td>", flags_col, "</td>",
  "</tr>",
  collapse = "\n"
)
summary_table <- paste0(
  "<table class='qc-table'><thead><tr><th>Pool</th><th>Sample</th>",
  if (has_mapping) "<th>Mapping % (Exact)</th>" else "",
  if (has_rpg) "<th>Reads per Guide</th>" else "",
  "<th>Gini</th><th>90/10 Ratio</th><th>% Undetected</th><th>Flags</th>",
  "</tr></thead><tbody>", summary_rows, "</tbody></table>"
)

# ---- Faceted distributional plots (one panel per pool) ---------------------
angled <- theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 35, hjust = 1))
centered <- theme(plot.title = element_text(hjust = 0.5))
many_samples <- length(unique(counts_all$sample)) > 12
legend_guard <- if (many_samples) guides(color = "none", fill = "none") else NULL

svg_guide_rep <- {
  p <- ggplot(ranked_all, aes(x = rank, y = log10(avg_cpm + 1))) +
    geom_area(fill = "purple4", alpha = 0.3) +
    geom_line(color = "purple4", linewidth = 0.7) +
    facet_wrap(~pool, ncol = facet_ncol, scales = "free_x") +
    labs(x = "gRNA Rank (lowest to highest)", y = "Mean CPM + 1 (log10)",
         title = paste("Guide Representation by Pool —", project)) +
    my_theme + centered
  plot_to_svg(p, "cmp-guide-rep", facet_w, facet_h)
}

# One figure, one block per pool: ECDF curves on top and, below on the same
# log10 x axis, one bar per sample = its mean read count per guide (zeros
# included), coloured like its curve. The mean is the count every guide would
# have if reads were split perfectly evenly. Blocks are gtables placed in a
# facet_ncol-column grid of absolute sizes (no patchwork/cowplot needed).
ecdf_block <- function(pl) {
  d <- counts_all[counts_all$pool == pl, ]
  smp <- sort(unique(d$sample))
  m <- d |> group_by(sample) |> summarize(mean_count = mean(count), .groups = "drop")
  m$pos <- length(smp) + 1 - match(m$sample, smp) # first sample at the top
  xl <- c(1, max(max(d$count) + 1, max(m$mean_count + 1) * 4))
  p1 <- ggplot(d, aes(x = count + 1, color = sample)) +
    stat_ecdf(geom = "step", linewidth = 0.9) +
    scale_x_log10(limits = xl, expand = expansion(mult = c(0.02, 0.05)),
                  labels = short_label) +
    scale_color_manual(values = sample_cols) +
    labs(x = NULL, y = "Cumulative Fraction", title = pl) +
    my_theme + centered + theme(legend.position = "none")
  p2 <- ggplot(m) +
    geom_rect(aes(xmin = 1, xmax = mean_count + 1, ymin = pos - 0.35, ymax = pos + 0.35, fill = sample)) +
    geom_text(aes(x = mean_count + 1, y = pos,
                  label = trimws(format(round(mean_count, 1), big.mark = ",", nsmall = 1))),
              hjust = -0.1, size = 2.8) +
    scale_x_log10(limits = xl, expand = expansion(mult = c(0.02, 0.05)),
                  labels = short_label) +
    scale_y_continuous(breaks = m$pos, labels = m$sample, expand = c(0.03, 0.03)) +
    scale_fill_manual(values = sample_cols) +
    labs(x = "Read Count + 1 (log10)", y = NULL) +
    my_theme + theme(legend.position = "none")
  g <- rbind(ggplotGrob(p1), ggplotGrob(p2), size = "max") # gtable method (gtable_rbind is not exported here)
  panel_rows <- sort(g$layout$t[grepl("^panel", g$layout$name)])
  g$heights[panel_rows[1]] <- unit(2.6, "in")
  g$heights[panel_rows[2]] <- unit(max(1.1, 0.2 * length(smp) + 0.4), "in")
  g
}
svg_ecdf <- {
  blocks <- lapply(pools, ecdf_block)
  # Total height of each block in inches (absolute panel heights + text rows).
  grDevices::pdf(NULL)
  blk_h <- vapply(blocks, function(g) grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE), numeric(1))
  invisible(grDevices::dev.off())
  ecdf_rows <- ceiling(n_pools / facet_ncol)
  row_h <- vapply(seq_len(ecdf_rows), function(r) max(blk_h[((r - 1) * facet_ncol + 1):min(r * facet_ncol, n_pools)]), numeric(1))
  col_w <- rep(3.6, facet_ncol)
  title_h <- 0.45
  outer <- gtable::gtable(widths = unit(col_w, "in"), heights = unit(c(title_h, row_h), "in"))
  outer <- gtable::gtable_add_grob(outer,
    grid::textGrob(paste("ECDF by Pool —", project), gp = grid::gpar(fontface = "bold", fontsize = 12)),
    t = 1, l = 1, r = facet_ncol, name = "title")
  for (i in seq_len(n_pools)) {
    outer <- gtable::gtable_add_grob(outer, blocks[[i]],
      t = 1 + ceiling(i / facet_ncol), l = ((i - 1) %% facet_ncol) + 1, clip = "off", name = paste0("block", i))
  }
  plot_to_svg(outer, "cmp-ecdf", sum(col_w), title_h + sum(row_h))
}

svg_lorenz <- {
  p <- ggplot(lorenz_all, aes(x = cum_gRNAs, y = cum_reads, color = sample)) +
    geom_line(linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~pool, ncol = facet_ncol) +
    labs(x = "Cumulative % of Guides", y = "Cumulative % of Reads",
         title = paste("Lorenz Curves by Pool —", project)) +
    my_theme + centered + theme(legend.position = "bottom") + legend_guard
  plot_to_svg(p, "cmp-lorenz", facet_w, facet_h + 0.6)
}

svg_density <- {
  p <- ggplot(counts_all, aes(x = count + 1, fill = sample)) +
    geom_density(alpha = 0.3) +
    scale_x_log10() +
    facet_wrap(~pool, ncol = facet_ncol) +
    labs(x = "Read Count + 1 (log10)", y = "Density",
         title = paste("Read Count Density by Pool —", project)) +
    my_theme + centered + theme(legend.position = "bottom") + legend_guard
  plot_to_svg(p, "cmp-density", facet_w, facet_h + 0.6)
}

ctrl_pools <- counts_all |> group_by(pool) |> summarize(has_ctrl = any(as.logical(is_control)), .groups = "drop") |>
  filter(has_ctrl) |> pull(pool)
svg_controls <- if (length(ctrl_pools) > 0) {
  d <- counts_all |> filter(pool %in% ctrl_pools)
  ctrl_ncol <- min(3, length(ctrl_pools))
  p <- ggplot(d, aes(x = count + 1, color = as.logical(is_control))) +
    stat_ecdf(geom = "step", linewidth = 0.9) +
    scale_x_log10() +
    facet_wrap(~pool, ncol = ctrl_ncol) +
    scale_color_manual(values = c("FALSE" = "gray70", "TRUE" = "red"),
                        labels = c("FALSE" = "Targeting", "TRUE" = "Control")) +
    labs(x = "Read Count + 1 (log10)", y = "Cumulative Fraction", color = "gRNA Type",
         title = paste("Control vs Targeting by Pool —", project)) +
    my_theme + centered + theme(legend.position = "bottom")
  plot_to_svg(p, "cmp-controls", 3.3 * ctrl_ncol + 1, 2.9 * ceiling(length(ctrl_pools) / ctrl_ncol) + 1.6)
} else NULL

svg_gc <- if (!is.null(gc_all)) {
  p <- ggplot(gc_all, aes(x = gc_pct, y = log10(norm_count + 1), color = sample)) +
    geom_point(alpha = 0.12, size = 0.5) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 1) +
    facet_wrap(~pool, ncol = facet_ncol) +
    labs(x = "Spacer GC Content (%)", y = "Normalized Read Count + 1 (log10)",
         title = paste("GC Content vs Read Depth by Pool —", project)) +
    my_theme + centered + theme(legend.position = "bottom") + legend_guard
  plot_to_svg(p, "cmp-gc", facet_w, facet_h + 0.6)
} else NULL

# mageck count's own matching is always exact, so this breaks each pool's
# UNMAPPED reads down further: near-miss (1 substitution from a real guide)
# vs far-miss (not close to any guide) — see bin/guide_representation_wrangle.R.
svg_class <- if (!is.null(class_all)) {
  class_long <- class_all |>
    select(pool, sample, pct_exact_match, pct_near_miss, pct_far_miss) |>
    pivot_longer(c(pct_exact_match, pct_near_miss, pct_far_miss), names_to = "category", values_to = "pct") |>
    mutate(category = factor(
      category,
      levels = c("pct_far_miss", "pct_near_miss", "pct_exact_match"),
      labels = c("Far-miss (unrelated)", "Near-miss (1 mismatch)", "Exact match")
    ))
  p <- ggplot(class_long, aes(x = sample, y = pct, fill = category)) +
    geom_col(width = 0.65) +
    scale_fill_manual(values = c(
      "Exact match"            = "steelblue",
      "Near-miss (1 mismatch)" = "goldenrod2",
      "Far-miss (unrelated)"   = "gray60"
    )) +
    scale_y_continuous(limits = c(0, 100.5), expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    facet_wrap(~pool, ncol = facet_ncol) +
    labs(x = NULL, y = "% of Reads", fill = NULL,
         title = paste("Read Classification by Pool —", project)) +
    my_theme + centered + theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  plot_to_svg(p, "cmp-class", facet_w, facet_h + 1.4)
} else NULL

# ---- Reads per library (bin/library_breakdown.py) --------------------------
# One stacked bar per sample: % of its reads matching each library, plus reads
# matching no library. Libraries are shown under their compact pool name.
breakdown_html <- NULL
if (!is.null(breakdown_all)) {
  any_lab <- "Any library"; none_lab <- "No library guide"
  bd <- breakdown_all
  bd$pct_of_sample_reads <- suppressWarnings(as.numeric(bd$pct_of_sample_reads))
  bd$mapped_reads <- suppressWarnings(as.numeric(bd$mapped_reads))
  per_lib <- bd[!(bd$library %in% c(any_lab, none_lab)), ]
  lib_levels <- unique(per_lib$library)
  sample_levels <- unique(bd$sample)
  stack <- rbind(
    per_lib[, c("sample", "library", "pct_of_sample_reads")],
    bd[bd$library == none_lab, c("sample", "library", "pct_of_sample_reads")]
  )
  stack$library <- factor(stack$library, levels = c(lib_levels, none_lab))
  stack$sample  <- factor(stack$sample, levels = rev(sample_levels))
  lib_cols <- c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948")
  fills <- if (length(lib_levels) <= length(lib_cols)) lib_cols[seq_along(lib_levels)] else scales::hue_pal()(length(lib_levels))
  fills <- c(setNames(fills, lib_levels), setNames("gray75", none_lab))
  p <- ggplot(stack, aes(x = sample, y = pct_of_sample_reads, fill = library)) +
    geom_col(width = 0.7, color = "white", linewidth = 0.3) +
    coord_flip() +
    scale_fill_manual(values = fills, breaks = c(lib_levels, none_lab)) +
    scale_y_continuous(limits = c(0, 100.5), expand = c(0, 0), labels = function(x) paste0(x, "%")) +
    labs(x = NULL, y = "% of the sample's reads", fill = NULL,
         title = paste("Reads per Library —", project)) +
    my_theme + centered + theme(legend.position = "bottom")
  svg_breakdown <- plot_to_svg(p, "cmp-breakdown", 9, max(3.2, 0.42 * length(sample_levels) + 2.2))

  # Table: one row per sample, % of reads per library (main library in bold).
  cell <- function(smp, lib) {
    r <- bd[bd$sample == smp & bd$library == lib, ]
    if (nrow(r) == 0) return("<td class='num'>&mdash;</td>")
    v <- fmt_pct(r$pct_of_sample_reads[1], 1)
    if (identical(r$status[1], "dropped_zero_mapped")) v <- "0% <span class='hint'>(not in table)</span>"
    paste0("<td class='num'>", v, "</td>")
  }
  rows <- vapply(sample_levels, function(smp) {
    reads <- bd$reads[bd$sample == smp][1]
    lib_cells <- vapply(lib_levels, function(l) cell(smp, l), character(1))
    pl <- per_lib[per_lib$sample == smp, ]
    main <- if (nrow(pl) > 0 && max(pl$mapped_reads, na.rm = TRUE) > 0) pl$library[which.max(pl$mapped_reads)] else NA
    lib_cells[!is.na(main) & lib_levels == main] <- sub("<td class='num'>", "<td class='num main-lib'>", lib_cells[!is.na(main) & lib_levels == main])
    any_pct  <- bd$pct_of_sample_reads[bd$sample == smp & bd$library == any_lab][1]
    none_pct <- bd$pct_of_sample_reads[bd$sample == smp & bd$library == none_lab][1]
    share <- if (!is.na(main)) round(100 * max(pl$mapped_reads) / max(1, bd$mapped_reads[bd$sample == smp & bd$library == any_lab][1])) else NA
    b <- if (is.na(any_pct) || any_pct < 5) badge("No library", "flag")
         else if (!is.na(share) && share < 80) badge(paste0("Mixed: ", share, "% ", pool_of(paste0(project, "_", main))), "flag")
         else badge(pool_of(paste0(project, "_", main)), "ok")
    paste0("<tr><td>", html_escape(smp), "</td><td class='num'>", fmt_int(reads), "</td>",
           paste(lib_cells, collapse = ""), "<td class='num'>", fmt_pct(any_pct, 1), "</td><td class='num'>",
           fmt_pct(none_pct, 1), "</td><td>", b, "</td></tr>")
  }, character(1))
  breakdown_table <- paste0(
    "<div class='table-scroll'><table class='qc-table'><thead><tr><th>Sample</th><th>Reads</th>",
    paste0("<th>", html_escape(vapply(lib_levels, function(l) pool_of(paste0(project, "_", l)), character(1))), "</th>", collapse = ""),
    "<th>Any library</th><th>No library guide</th><th>Main library</th></tr></thead><tbody>",
    paste(rows, collapse = "\n"), "</tbody></table></div>")
  breakdown_html <- paste0("<div class='plot-row'><div class='plot plot-wide'>", svg_breakdown, "</div></div>", breakdown_table)
}

# ---- Headline metrics: combined bar charts (one bar per pool/sample) ------
bar_compare <- function(y, y_lab, title, hline_good = NULL, hline_flag = NULL) {
  p <- ggplot(qc_all, aes(x = sample, y = .data[[y]], fill = pool)) +
    geom_col(width = 0.7)
  if (!is.null(hline_good)) p <- p + geom_hline(yintercept = hline_good, linetype = "dashed", color = "forestgreen", linewidth = 0.7)
  if (!is.null(hline_flag)) p <- p + geom_hline(yintercept = hline_flag, linetype = "dashed", color = "firebrick", linewidth = 0.7)
  p <- p + facet_grid(. ~ pool, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = y_lab, title = title) + my_theme + centered +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text.x = element_text(size = 10, face = "bold"),
          panel.spacing = unit(0.3, "lines"))
  plot_to_svg(p, paste0("cmp-bar-", y), max(8, 0.32 * nrow(qc_all) + 2.5 + 0.35 * n_pools), 5.2)
}
svg_gini  <- bar_compare("Gini", "Gini Coefficient", "Gini by Pool")
svg_ratio <- bar_compare("ratio_90_10", "90/10 Ratio", "Library Skew by Pool", hline_good = 5, hline_flag = 10)
svg_undet <- bar_compare("pct_undetected", "% Undetected", "Undetected Guides by Pool", hline_good = 1, hline_flag = 5)
svg_map   <- if (has_mapping) bar_compare("mapping_pct", "% Mapped (Exact)", "Mapping Rate / % Exact Match by Pool", hline_flag = 30) else NULL

# ---- PCA of samples, one panel per pool -------------------------------------
# Greedy deterministic label placement (no ggrepel), same helper as
# guide_representation_plot.R: points sorted by x each take the first of
# above / below / right / left (then diagonals, then farther rings) whose label
# box overlaps no box already placed (other labels, all points) and stays
# inside the axis limits; otherwise the candidate with the least overlap.
place_labels <- function(x, y, labels, xlim, ylim) {
  xs <- diff(xlim); ys <- diff(ylim)
  w <- nchar(labels) * 0.015 * xs
  h <- 0.05 * ys
  ox <- 0.025 * xs; oy <- 0.03 * ys
  pr_x <- 0.012 * xs; pr_y <- 0.012 * ys
  # boxes: xmin, xmax, ymin, ymax; all points are obstacles from the start
  boxes <- cbind(x - pr_x, x + pr_x, y - pr_y, y + pr_y)
  overlap <- function(b, others) {
    if (nrow(others) == 0) return(0)
    ow <- pmax(0, pmin(b[2], others[, 2]) - pmax(b[1], others[, 1]))
    oh <- pmax(0, pmin(b[4], others[, 4]) - pmax(b[3], others[, 3]))
    sum(ow * oh)
  }
  out <- data.frame(x = x, y = y, hjust = 0.5, vjust = 0, label = labels)
  placed <- boxes[0, , drop = FALSE]
  for (i in order(x, y)) {
    # candidate anchors: above, below, right, left, then the four diagonals;
    # tried at increasing distances (rings) from the point
    cands <- list()
    for (ring in c(1, 2.2, 3.6)) {
      dx <- ox * ring; dy <- oy * ring
      cands <- c(cands, list(
        list(px = x[i],      py = y[i] + dy, hj = 0.5, vj = 0,
             b = c(x[i] - w[i] / 2, x[i] + w[i] / 2, y[i] + dy, y[i] + dy + h)),
        list(px = x[i],      py = y[i] - dy, hj = 0.5, vj = 1,
             b = c(x[i] - w[i] / 2, x[i] + w[i] / 2, y[i] - dy - h, y[i] - dy)),
        list(px = x[i] + dx, py = y[i],      hj = 0,   vj = 0.5,
             b = c(x[i] + dx, x[i] + dx + w[i], y[i] - h / 2, y[i] + h / 2)),
        list(px = x[i] - dx, py = y[i],      hj = 1,   vj = 0.5,
             b = c(x[i] - dx - w[i], x[i] - dx, y[i] - h / 2, y[i] + h / 2)),
        list(px = x[i] + dx, py = y[i] + dy, hj = 0,   vj = 0,
             b = c(x[i] + dx, x[i] + dx + w[i], y[i] + dy, y[i] + dy + h)),
        list(px = x[i] - dx, py = y[i] + dy, hj = 1,   vj = 0,
             b = c(x[i] - dx - w[i], x[i] - dx, y[i] + dy, y[i] + dy + h)),
        list(px = x[i] + dx, py = y[i] - dy, hj = 0,   vj = 1,
             b = c(x[i] + dx, x[i] + dx + w[i], y[i] - dy - h, y[i] - dy)),
        list(px = x[i] - dx, py = y[i] - dy, hj = 1,   vj = 1,
             b = c(x[i] - dx - w[i], x[i] - dx, y[i] - dy - h, y[i] - dy))
      ))
    }
    others <- rbind(boxes[-i, , drop = FALSE], placed)
    score <- vapply(cands, function(cd) {
      inside <- cd$b[1] >= xlim[1] && cd$b[2] <= xlim[2] && cd$b[3] >= ylim[1] && cd$b[4] <= ylim[2]
      overlap(cd$b, others) + if (inside) 0 else 1e6 * xs * ys
    }, numeric(1))
    k <- which(score == 0)[1]
    if (is.na(k)) k <- which.min(score)
    cd <- cands[[k]]
    out[i, c("x", "y", "hjust", "vjust")] <- c(cd$px, cd$py, cd$hj, cd$vj)
    placed <- rbind(placed, cd$b)
  }
  out
}

svg_pca <- if (!is.null(pca_all) && !is.null(pcavar_all) && nrow(pca_all) > 0) {
  pad <- function(r) {
    d <- diff(r)
    if (d == 0) d <- max(abs(r[1]), 1)
    r + c(-1, 1) * 0.22 * d
  }
  pca_pools <- pools[pools %in% pca_all$pool]
  lab_list <- list(); blank_list <- list(); pt_list <- list()
  for (pl in pca_pools) {
    d <- pca_all[pca_all$pool == pl, ]
    v <- pcavar_all[pcavar_all$pool == pl, ]
    v1 <- v$pct_variance[v$PC == "PC1"][1]; v2 <- v$pct_variance[v$PC == "PC2"][1]
    lbl <- paste0(pl, "\nPC1 ", sprintf("%.1f", v1), "%, PC2 ", sprintf("%.1f", v2), "%")
    xlim <- pad(range(d$PC1)); ylim <- pad(range(d$PC2))
    ld <- place_labels(d$PC1, d$PC2, d$sample, xlim, ylim)
    ld$sample <- d$sample; ld$pool_label <- lbl
    lab_list[[pl]] <- ld
    blank_list[[pl]] <- data.frame(PC1 = xlim, PC2 = ylim, pool_label = lbl)
    d$pool_label <- lbl
    pt_list[[pl]] <- d
  }
  lab_df <- dplyr::bind_rows(lab_list); blank_df <- dplyr::bind_rows(blank_list); pt_df <- dplyr::bind_rows(pt_list)
  lvls <- unique(pt_df$pool_label)
  lab_df$pool_label <- factor(lab_df$pool_label, levels = lvls)
  blank_df$pool_label <- factor(blank_df$pool_label, levels = lvls)
  pt_df$pool_label <- factor(pt_df$pool_label, levels = lvls)
  p <- ggplot(pt_df, aes(x = PC1, y = PC2, color = sample)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    geom_blank(data = blank_df, aes(x = PC1, y = PC2), inherit.aes = FALSE) +
    geom_point(size = 2.5) +
    geom_text(data = lab_df, aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
              size = 2.8, show.legend = FALSE) +
    scale_color_manual(values = sample_cols, guide = "none") +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~pool_label, scales = "free", ncol = facet_ncol) +
    labs(x = "PC1", y = "PC2", title = paste("PCA of Samples by Pool —", project)) +
    my_theme + centered
  plot_to_svg(p, "cmp-pca", facet_w, facet_h + 0.4)
} else NULL

# ---- Read depth per sample ---------------------------------------------------
# Total reads per sample: from the library breakdown if present, else the
# qc_summary reads column, else the read_depth tables.
depth_df <- NULL
one_per_sample <- function(smp, tot) {
  tot <- suppressWarnings(as.numeric(tot))
  ok <- !is.na(tot)
  if (!any(ok)) return(NULL)
  d <- data.frame(sample = smp[ok], reads = tot[ok])
  d <- d |> group_by(sample) |> summarize(reads = max(reads), .groups = "drop")
  d[match(unique(smp[ok]), d$sample), ]
}
if (!is.null(breakdown_all) && "reads" %in% names(breakdown_all)) depth_df <- one_per_sample(breakdown_all$sample, breakdown_all$reads)
if (is.null(depth_df) && "reads" %in% names(qc_all)) depth_df <- one_per_sample(qc_all$sample, qc_all$reads)
if (is.null(depth_df) && !is.null(depth_all) && "total_reads" %in% names(depth_all)) depth_df <- one_per_sample(depth_all$sample, depth_all$total_reads)
if (!is.null(depth_df)) depth_df <- depth_df[depth_df$reads > 0, ]
svg_depth <- if (!is.null(depth_df) && nrow(depth_df) > 0) {
  n <- nrow(depth_df)
  depth_df$pos <- n + 1 - seq_len(n) # first sample at the top
  xmin_d <- 10^floor(log10(min(depth_df$reads)))
  p <- ggplot(depth_df) +
    geom_rect(aes(xmin = xmin_d, xmax = reads, ymin = pos - 0.35, ymax = pos + 0.35), fill = "steelblue") +
    geom_text(aes(x = reads, y = pos, label = format(round(reads), big.mark = ",")), hjust = -0.1, size = 3.2) +
    scale_x_log10(labels = short_label,
                  expand = expansion(mult = c(0, 0.15))) +
    scale_y_continuous(breaks = depth_df$pos, labels = depth_df$sample, expand = c(0.02, 0.02)) +
    labs(x = "Total Reads (log10)", y = NULL, title = paste("Read Depth per Sample —", project)) +
    my_theme + centered
  plot_to_svg(p, "cmp-read-depth", 9, max(4, 0.28 * n + 1.8))
} else NULL

# ---- Sections ---------------------------------------------------------------
plot_div <- function(svg) paste0("<div class='plot-row'><div class='plot plot-wide'>", svg, "</div></div>")
missing_div <- function(reason) paste0("<div class='missing-box'><strong>Not generated.</strong> ", reason, "</div>")

# Each section is titled with the question it answers (num. question) and keeps
# its old name as a subtitle. c(question, subtitle) per section id.
sec_head <- list(
  breakdown = c("0. Which library does each sample's DNA belong to, and what share of its reads matches no library?", "Reads per Library"),
  headline = c("1. Do any pools stand out on the key quality numbers?", "Headline Metrics"),
  guide_representation = c("2. Are each pool's guides represented evenly in each sample?", "Guide Representation (Ranked)"),
  ecdf = c("3. What fraction of each pool's guides received few reads, and how deep is each sample on average?", "ECDF of gRNA Representation"),
  lorenz = c("4. How unequally are reads shared among each pool's guides?", "Lorenz Curves"),
  density = c("5. Do each pool's guides form one population or several (e.g. a dropout group near zero)?", "Read-Count Density"),
  controls = c("6. Do control guides behave differently from targeting guides in each pool?", "Control vs Targeting gRNA Representation"),
  class = c("7. Why didn't the other reads match a guide: sequencing errors or unrelated sequences?", "Read Classification: Exact / Near-Miss / Far-Miss"),
  gc = c("8. Does a guide's GC content bias how many reads it gets?", "GC Content vs Read Depth"),
  pca = c("9. Which samples resemble each other overall in each pool, and do replicates group together?", "Principal Component Analysis (PCA) by Pool"),
  depth = c("10. Was each sample sequenced deeply enough?", "Read Depth per Sample")
)

section <- function(id, description, body, interpretation, troubleshooting) {
  troubleshooting_html <- if (grepl("<li>", troubleshooting, fixed = TRUE)) paste0("<ul>", troubleshooting, "</ul>") else paste0("<p>", troubleshooting, "</p>")
  paste0(
    "<section id='", id, "'>", "<h2>", html_escape(sec_head[[id]][1]), "</h2>",
    "<p class='subtitle'>", html_escape(sec_head[[id]][2]), "</p>", "<p class='desc'>", description, "</p>", body,
    "<div class='callout callout-interp'><h3>How to interpret it</h3><p>", interpretation, "</p></div>",
    "<div class='callout callout-trouble'><h3>If one pool looks off</h3>", troubleshooting_html, "</div>",
    "</section>"
  )
}

sections_html <- paste(c(
  if (!is.null(breakdown_html)) {
    section("breakdown", paste0("For every sample: the share of its reads that match a guide of each library it was counted against, and the share matching no library at all. A read matching a guide shared by several libraries counts once in &ldquo;Any library&rdquo;. With no <code>library</code> column in <code>samplesheet.csv</code>, every sample is counted against every library, so this shows which pool(s) each sample really holds. Numbers: <code>results/counts/qc/tables/", project, "_library_breakdown.csv</code> (one row per sample &times; library) and <code>_library_breakdown_wide.csv</code> (one row per sample)."),
      breakdown_html,
      paste0(r"(A sample screened with one pool shows one dominant colour (its main library, in bold in the
      table) and little else. Several sizeable colours mean the sample holds a mix of pools. A sample
      that is almost entirely grey matches no library: it is empty, failed, or from a different screen.
      In the per-pool sections below, every sample appears in every pool's panel; a sample that holds
      few reads of a pool will look like heavy dropout there. Read those panels together with this one.)",
      if (all_vs_all) paste0(" ", aa_note_breakdown) else ""),
      r"(<li><b>No library</b>: check the sample's FastQC report (short or poly-G reads usually mean an
      empty or failed library) and whether it belongs to this screen at all.</li>
      <li><b>Mixed</b>: expected if the sample was meant to hold several pools; otherwise check for a
      pooling or sample-sheet mix-up.</li>
      <li>&ldquo;0% (not in table)&rdquo;: the sample had no read of that library, so it was left out of that
      library's count table (only done when every sample is counted against every library).</li>)")
  },
  section("headline", "The four numbers behind each pool's flags, one bar chart per metric (one panel per pool, one bar per sample). Dashed lines mark the same good/flagged thresholds used throughout this report.",
    paste0(plot_div(svg_gini), plot_div(svg_ratio), plot_div(svg_undet),
           if (!is.null(svg_map)) plot_div(svg_map) else ""),
    paste0("A pool above the flagged line on any of these has that specific problem — see that metric's own section below for what it means and what to check.",
           if (all_vs_all) paste0(" ", aa_note_map) else ""),
    r"(<li>A pool flagged on several metrics at once usually indicates a labeling or sample mix-up
    in <code>samplesheet.csv</code>, rather than several unrelated problems.</li>
    <li>% Undetected near 100% for a pool is a different, more severe problem than ordinary
    dropout: it usually means that sample wasn't actually screened with this library, rather than
    its guides genuinely failing one by one.</li>)"),

  section("guide_representation", "Each panel shows one pool's guides, sorted from lowest to highest mean CPM. Each guide's value is its mean counts per million reads (CPM) across the pool's samples, which puts samples of different depth on one scale.",
    plot_div(svg_guide_rep),
    r"(A steep drop at the low-rank end of a panel indicates a subset of that pool's guides are
    barely present or missing. A steep rise concentrated at the high-rank end indicates a small
    number of guides amplified far beyond the rest. A distinct step or plateau in the middle of a
    panel's curve indicates two sub-populations of guides at different typical abundance in that
    pool — the same signal the Density panel below shows as two peaks.)",
    r"(<li>Cross-check that pool's % Undetected and 90/10 ratio in the headline metrics above.</li>
    <li>A drop already present at T0/plasmid indicates a library or transduction coverage problem
    specific to that pool.</li>)"),

  section("ecdf", "One block per pool. Top: the cumulative fraction of that pool's guides at or below a given read count, one curve per sample. Bottom: one bar per sample showing its mean read count per guide (zeros included) on the same axis &mdash; what every guide would have if reads were split perfectly evenly. Bars and curves share the sample's colour and the bar value is printed at its end.",
    plot_div(svg_ecdf),
    "Compare a curve with its bar below: the more of the curve sitting to the left of the end of its bar, the larger the fraction of that pool's guides falling below their fair share of reads. A short bar means the sample is shallow for that pool (few reads per guide); bars of very different length between replicates mean uneven depth.",
    r"(<li>A curve that rises far to the left of its bar end indicates that pool needs closer
    attention — cross-check its mapping rate and % undetected in the headline metrics above.</li>
    <li>A bar far shorter than the others in the same pool: that sample has few reads per guide in
    this pool (shallow, or mostly from another pool in all-vs-all mode) — see the Reads per Library and
    Read Depth sections.</li>)"),

  section("lorenz", "Each panel shows cumulative % of reads against cumulative % of guides for one pool. The diagonal represents an ideal pool with perfectly even representation.",
    plot_div(svg_lorenz),
    "The further a panel's curve bows away from the diagonal, the less evenly that pool's guides are represented — that pool's Gini coefficient in the headline metrics above will be higher too.",
    "<li>Confirm that pool's 90/10 ratio above agrees with what the curve shows.</li>"),

  section("density", "Each panel shows the read-count distribution for one pool.",
    plot_div(svg_density),
    r"(A single, smooth peak in a panel is the ideal case. A second peak near zero read counts
    indicates gRNA dropout in that pool; the height of that peak reflects severity — a small bump
    affects a handful of guides, a tall peak affects a large fraction of the library. A second peak
    away from zero indicates two distinct sub-populations of guides at different typical
    abundance, rather than dropout — check the Control vs Targeting panel for that pool below to
    see if it splits along the same lines (often intentional, e.g. controls spiked in at a
    different ratio).)",
    r"(<li>The exact dropped-out guides for a pool are listed in that pool's
    <code>zero_count_guides.csv</code> (in <code>results/counts/qc/tables/</code>), useful for
    checking whether they cluster by gene.</li>
    <li>A peak away from zero that doesn't line up with controls: check that pool's own
    samplesheet.csv assignment for a possible sample mix-up.</li>)"),

  if (!is.null(svg_controls)) {
    section("controls", "For pools with control/NTC/non-targeting guides: control vs targeting read-count distributions, one panel per pool.",
      plot_div(svg_controls),
      "Divergence between a panel's two curves already present at a reference/T0 timepoint indicates a technical difference between how that pool's controls and targeting guides were made, not biology.",
      "<li>Compare a pool's own T0/plasmid panel against its treated timepoint, if this run has both.</li>")
  } else {
    section("controls", "For pools with control/NTC/non-targeting guides: control vs targeting read-count distributions, one panel per pool.",
      missing_div("None of this run's libraries have control, NTC, or non-targeting guides in their Gene column."),
      "", "")
  },

  if (!is.null(svg_class)) {
    section("class", "Every pool's reads, split into exact matches (= that pool's Mapping Rate above), near-misses (1 substitution from a guide), and far-misses (not close to any guide), one panel per pool.",
      plot_div(svg_class),
      r"(A large exact-match share with little near- or far-miss is the ideal case. A pool with a
      substantial near-miss share is losing real signal to sequencing/PCR error &mdash; those reads
      almost certainly came from a real guide. A pool with a substantial far-miss share has
      something else going on: contamination, adapter-only reads, or a sample that doesn't belong
      to that library at all.)",
      r"(<li>A pool that's mostly near-miss: check its FastQC/MultiQC base-quality metrics, and
      confirm <code>trim5</code>/<code>sgrna_len</code> in <code>config.yaml</code> exactly match
      the amplicon design for that pool.</li>
      <li>A pool that's mostly far-miss: check the Undetected Guides and Mapping Rate bars for that
      pool in the headline metrics above — the same signal, another angle.</li>)")
  } else {
    section("class", "Every pool's reads, split into exact matches, near-misses, and far-misses, one panel per pool.",
      missing_div("None of this run's libraries had a usable combination of library sequences, mageck log, and unmapped-reads file to classify reads with."),
      "", "")
  },

  if (!is.null(svg_gc)) {
    section("gc", "Guide GC content vs. normalized read count, one panel per pool.",
      plot_div(svg_gc),
      r"(A flat trend line is the ideal case. The steeper a monotonic slope, the stronger the
      GC-dependent amplification bias for that pool. A U-shaped or dome-shaped trend &mdash; dipping
      or peaking in the middle of the GC range instead of running in one direction &mdash; is the
      classic PCR-bias signature, with both GC extremes disadvantaged (or advantaged) relative to
      the middle. Some slope is common even in good pools; it's the pools with a distinctly steeper
      slope or a pronounced U/dome shape, compared to the rest, that are worth a closer look.)",
      r"(<li>A pool with a distinctly steeper slope than its neighbors is worth checking for a
      library-PCR difference specific to it (cycle number, polymerase, annealing temperature).</li>
      <li>mageck's normalization doesn't correct for GC bias, so it's worth noting even when it
      doesn't invalidate the screen.</li>)")
  } else {
    section("gc", "Guide GC content vs. normalized read count, one panel per pool.",
      missing_div("None of this run's libraries had the guide id/sequence columns needed to compute GC content."),
      "", "")
  },

  if (!is.null(svg_pca)) {
    section("pca", "One panel per pool, one point per sample. Each pool's PCA uses log2(counts per million + 1) of that pool's guides, centred but not scaled; PC1 is plotted against PC2 and the percentage of variance each explains is in the panel title. Pools with fewer than 3 samples with reads in them are left out (no PCA).",
      plot_div(svg_pca),
      paste0("Points close together have similar guide profiles; replicates should cluster. Separation along PC1 by condition or timepoint is expected after selection. A sample far from its replicates is suspect.",
             if (all_vs_all) " In all-vs-all mode every sample is counted against every pool, so samples that belong to other pools have very few reads in a pool and dominate or separate in that pool's PCA &mdash; read it together with section 0 (Reads per Library)." else ""),
      "<li>A sample far from its replicates: check its mapping rate and depth in the sections above, and its row in <code>samplesheet.csv</code>.</li>")
  } else {
    section("pca", "One panel per pool, one point per sample (PCA of log2 counts per million of that pool's guides).",
      missing_div("No pool had at least 3 samples with reads in it, or the PCA tables were not found."),
      "", "")
  },

  if (!is.null(svg_depth)) {
    section("depth", "One bar per sample: its total sequenced reads, on a log10 axis, with the number at the end of the bar. The reads you need are roughly the number of guides times the reads per guide you want. A common rule of thumb is at least 200-500 reads per guide at the end point of a screen (plasmid/T0 samples often need more). The Reads per Guide column in At a Glance gives each pool's matched reads per guide (matched reads / guides). This is a rule of thumb, not a hard cutoff.",
      plot_div(svg_depth),
      "Samples of similar depth are directly comparable. A sample far below the others, or below guides x coverage, has noisy counts for many guides.",
      r"(<li>A sample far below the others (or below guides &times; the coverage you need): too shallow, so noisy counts &mdash; consider re-sequencing it.</li>
      <li>Very uneven depth between replicates lowers their correlation; compare with the PCA above.</li>)")
  } else {
    section("depth", "One bar per sample: its total sequenced reads.",
      missing_div("No total read counts were found (no library breakdown, no reads column in the qc summaries, no read_depth tables)."),
      "", "")
  }
), collapse = "\n")

toc_ids    <- c(if (!is.null(breakdown_html)) "breakdown", "headline", "guide_representation", "ecdf", "lorenz", "density", "controls", "class", "gc", "pca", "depth")
toc_titles <- vapply(toc_ids, function(i) sec_head[[i]][1], character(1))
toc_items  <- paste0("<li><a href='#", toc_ids, "'>", html_escape(toc_titles), "</a></li>", collapse = "\n")

# ---- Assemble the page -------------------------------------------------------
css <- r"(
:root {
  --bg: #ffffff; --fg: #1a1a1a; --muted: #6b6b6b; --border: #e2e2e2;
  --accent: #2f5d8a; --ok-bg: #e6f4ea; --ok-fg: #1e7a34;
  --flag-bg: #fdecea; --flag-fg: #b3261e; --card-bg: #f7f8fa;
}
* { box-sizing: border-box; }
body {
  font-family: -apple-system, "Helvetica Neue", Helvetica, Arial, sans-serif;
  color: var(--fg); background: var(--bg); margin: 0; line-height: 1.55;
}
.wrap { max-width: 1180px; margin: 0 auto; padding: 2rem 1.5rem 4rem; }
header.page-header { border-bottom: 3px solid var(--accent); padding-bottom: 1rem; margin-bottom: 1.5rem; }
header.page-header h1 { margin: 0 0 0.25rem; font-size: 1.6rem; }
header.page-header .meta { color: var(--muted); font-size: 0.92rem; }
.about { background: var(--card-bg); border: 1px solid var(--border); border-radius: 8px; padding: 1rem 1.25rem; margin-bottom: 1.5rem; font-size: 0.95rem; }
nav.toc { background: var(--card-bg); border: 1px solid var(--border); border-radius: 8px; padding: 1rem 1.25rem; margin-bottom: 2rem; }
nav.toc h2 { margin-top: 0; font-size: 1rem; text-transform: uppercase; letter-spacing: 0.04em; color: var(--muted); }
nav.toc ul { columns: 2; margin: 0; padding-left: 1.2rem; }
nav.toc a { color: var(--accent); text-decoration: none; }
nav.toc a:hover { text-decoration: underline; }
h2 { font-size: 1.3rem; margin-top: 3rem; border-top: 1px solid var(--border); padding-top: 1.5rem; }
.subtitle { margin: -0.3rem 0 1rem; font-size: 0.95rem; color: var(--muted); font-weight: 400; }
.desc { color: #333; }
.plot-row { display: flex; flex-wrap: wrap; gap: 1.5rem; justify-content: center; margin: 1rem 0; }
.plot { flex: 1 1 300px; max-width: 100%; text-align: center; }
.plot-wide { flex-basis: 100%; }
.plot svg { max-width: 100%; height: auto; }
.missing-box { background: var(--card-bg); border: 1px dashed var(--border); border-radius: 8px; padding: 1rem 1.25rem; color: var(--muted); font-size: 0.92rem; }
.callout { border-radius: 8px; padding: 0.9rem 1.25rem; margin-top: 0.9rem; font-size: 0.93rem; }
.callout h3 { margin: 0 0 0.35rem; font-size: 0.85rem; text-transform: uppercase; letter-spacing: 0.03em; }
.callout p { margin: 0; }
.callout ul { margin: 0; padding-left: 1.2rem; }
.callout li { margin: 0.15rem 0; }
.callout-interp { background: #eef3f8; }
.callout-interp h3 { color: var(--accent); }
.callout-trouble { background: #fff8ec; }
.callout-trouble h3 { color: #8a5a00; }
table.qc-table { border-collapse: collapse; width: 100%; font-size: 0.88rem; margin: 1rem 0; }
table.qc-table th, table.qc-table td { border: 1px solid var(--border); padding: 0.45rem 0.6rem; text-align: left; }
table.qc-table th { background: var(--card-bg); }
table.qc-table td.num { text-align: right; font-variant-numeric: tabular-nums; }
.badge { display: inline-block; border-radius: 999px; padding: 0.1rem 0.6rem; font-size: 0.78rem; font-weight: 600; white-space: nowrap; }
.badge-ok { background: var(--ok-bg); color: var(--ok-fg); }
.badge-flag { background: var(--flag-bg); color: var(--flag-fg); margin-right: 0.2rem; }
.badge-info { background: #e8eef6; color: #2b4a6f; margin-right: 0.2rem; }
.table-scroll { overflow-x: auto; }
td.main-lib { font-weight: 700; }
.hint { color: var(--muted); font-size: 0.8em; }
footer { color: var(--muted); font-size: 0.82rem; margin-top: 3rem; border-top: 1px solid var(--border); padding-top: 1rem; }
@media print { nav.toc { display: none; } h2 { break-before: page; } }
)"

html <- paste0(
  "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width, initial-scale=1'>",
  "<title>", html_escape(project), " &mdash; pool comparison QC</title>",
  "<style>", css, "</style></head><body><div class='wrap'>",
  "<header class='page-header'><h1>Guide Representation QC &mdash; Pool Comparison</h1>",
  "<div class='meta'>Project: <strong>", html_escape(project), "</strong> &nbsp;&middot;&nbsp; ",
  "Pools: <strong>", n_pools, "</strong> (", html_escape(paste(pools, collapse = ", ")), ") &nbsp;&middot;&nbsp; ",
  "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), "</div></header>",
  "<div class='about'>This report puts every pool in this run side by side, using the same plots ",
  "and metrics as each pool's own <code>qc_report.html</code> in this same folder. Replicate ",
  "correlation is the one thing left out here, since it can't be compared meaningfully across ",
  "pools that used different guides — see a pool's own report for that.</div>",
  "<nav class='toc'><h2>Sections</h2><ul>", toc_items, "</ul></nav>",
  "<section id='summary'><h2>How does each sample look overall, and which ones need attention?</h2><p class='subtitle'>At a Glance</p>", summary_table,
  if (all_vs_all) paste0("<p class='desc'>", aa_note_glance, "</p>") else "",
  "</section>",
  sections_html,
  "<footer>Generated by the CRISPR_SCREEN pipeline's guide_representation_compare_report.R from tables in ",
  html_escape(tables_dir), ".</footer>",
  "</div></body></html>"
)

out_path <- file.path(outdir, paste0(project, "_pools_comparison_qc_report.html"))
writeLines(html, out_path)
cat("  Saved:", out_path, "\n")
