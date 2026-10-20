#!/usr/bin/env Rscript
#' ##########################################################################
#' guide_representation_compare_report.R
#'
#' One HTML report comparing every library (pool) in this run side by side:
#' the distributional plots faceted one panel per pool (3 columns), plus the
#' headline scalar metrics (Gini, 90/10 ratio, % undetected, mapping rate) as
#' combined bar charts across pools. Built only from the tables
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
#' plot is embedded inline as SVG).
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

# Renders a ggplot to SVG and returns it as a namespaced, embeddable string.
plot_to_svg <- function(p, ns, width, height) {
  tmp <- tempfile(fileext = ".svg")
  grDevices::svg(filename = tmp, width = width, height = height)
  print(p)
  invisible(grDevices::dev.off())
  raw <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
  unlink(tmp)
  start <- regexpr("<svg[ >]", raw)
  end   <- regexpr("</svg>", raw)
  svg   <- substr(raw, start, end + attr(end, "match.length") - 1)
  namespace_svg_ids(svg, ns)
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

svg_ecdf <- {
  # Dashed line: each sample's own mean read count (see guide_representation_plot.R
  # for why a per-sample mean, not a fixed cutoff, is the right reference).
  means <- counts_all |> group_by(pool, sample) |> summarize(mean_count = mean(count), .groups = "drop")
  p <- ggplot(counts_all, aes(x = count + 1, color = sample)) +
    stat_ecdf(geom = "step", linewidth = 0.9) +
    geom_vline(data = means, aes(xintercept = mean_count + 1, color = sample),
               linetype = "dashed", linewidth = 0.6, alpha = 0.8, show.legend = FALSE) +
    scale_x_log10() +
    facet_wrap(~pool, ncol = facet_ncol) +
    labs(x = "gRNA Read Count + 1 (log10)", y = "Cumulative Fraction",
         title = paste("ECDF by Pool —", project)) +
    my_theme + centered + theme(legend.position = "bottom") + legend_guard
  plot_to_svg(p, "cmp-ecdf", facet_w, facet_h + 0.6)
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
    my_theme + centered + theme(legend.position = "bottom", axis.text.x = element_text(angle = 35, hjust = 1))
  plot_to_svg(p, "cmp-class", facet_w, facet_h + 0.8)
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
  p <- ggplot(qc_all, aes(x = x_label, y = .data[[y]], fill = pool)) +
    geom_col(width = 0.7)
  if (!is.null(hline_good)) p <- p + geom_hline(yintercept = hline_good, linetype = "dashed", color = "forestgreen", linewidth = 0.7)
  if (!is.null(hline_flag)) p <- p + geom_hline(yintercept = hline_flag, linetype = "dashed", color = "firebrick", linewidth = 0.7)
  p <- p + labs(x = NULL, y = y_lab, title = title) + my_theme + angled + theme(legend.position = "none")
  plot_to_svg(p, paste0("cmp-bar-", y), max(6, n_pools * 0.9 + 1.5), 4.6)
}
svg_gini  <- bar_compare("Gini", "Gini Coefficient", "Gini by Pool")
svg_ratio <- bar_compare("ratio_90_10", "90/10 Ratio", "Library Skew by Pool", hline_good = 5, hline_flag = 10)
svg_undet <- bar_compare("pct_undetected", "% Undetected", "Undetected Guides by Pool", hline_good = 1, hline_flag = 5)
svg_map   <- if (has_mapping) bar_compare("mapping_pct", "% Mapped (Exact)", "Mapping Rate / % Exact Match by Pool", hline_flag = 30) else NULL

# ---- Sections ---------------------------------------------------------------
plot_div <- function(svg) paste0("<div class='plot-row'><div class='plot plot-wide'>", svg, "</div></div>")
missing_div <- function(reason) paste0("<div class='missing-box'><strong>Not generated.</strong> ", reason, "</div>")

section <- function(id, title, description, body, interpretation, troubleshooting) {
  troubleshooting_html <- if (grepl("<li>", troubleshooting, fixed = TRUE)) paste0("<ul>", troubleshooting, "</ul>") else paste0("<p>", troubleshooting, "</p>")
  paste0(
    "<section id='", id, "'>", "<h2>", title, "</h2>", "<p class='desc'>", description, "</p>", body,
    "<div class='callout callout-interp'><h3>How to interpret it</h3><p>", interpretation, "</p></div>",
    "<div class='callout callout-trouble'><h3>If one pool looks off</h3>", troubleshooting_html, "</div>",
    "</section>"
  )
}

sections_html <- paste(c(
  if (!is.null(breakdown_html)) {
    section("breakdown", "0. Reads per Library", paste0("For every sample: the share of its reads that match a guide of each library it was counted against, and the share matching no library at all. A read matching a guide shared by several libraries counts once in &ldquo;Any library&rdquo;. With no <code>library</code> column in <code>samplesheet.csv</code>, every sample is counted against every library, so this shows which pool(s) each sample really holds. Numbers: <code>results/counts/qc/tables/", project, "_library_breakdown.csv</code> (one row per sample &times; library) and <code>_library_breakdown_wide.csv</code> (one row per sample)."),
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
  section("headline", "1. Headline Metrics", "The four numbers behind each pool's flags, one bar chart per metric across every pool. Dashed lines mark the same good/flagged thresholds used throughout this report.",
    paste0("<div class='plot-row'>",
           "<div class='plot'>", svg_gini, "</div>",
           "<div class='plot'>", svg_ratio, "</div>",
           "<div class='plot'>", svg_undet, "</div>",
           if (!is.null(svg_map)) paste0("<div class='plot'>", svg_map, "</div>") else "",
           "</div>"),
    paste0("A pool above the flagged line on any of these has that specific problem — see that metric's own section below for what it means and what to check.",
           if (all_vs_all) paste0(" ", aa_note_map) else ""),
    r"(<li>A pool flagged on several metrics at once usually indicates a labeling or sample mix-up
    in <code>samplesheet.csv</code>, rather than several unrelated problems.</li>
    <li>% Undetected near 100% for a pool is a different, more severe problem than ordinary
    dropout: it usually means that sample wasn't actually screened with this library, rather than
    its guides genuinely failing one by one.</li>)"),

  section("guide_representation", "2. Guide Representation (Ranked)", "Each panel shows one pool's guides, sorted from lowest to highest mean CPM. Each guide's value is its mean counts per million reads (CPM) across the pool's samples, which puts samples of different depth on one scale.",
    plot_div(svg_guide_rep),
    r"(A steep drop at the low-rank end of a panel indicates a subset of that pool's guides are
    barely present or missing. A steep rise concentrated at the high-rank end indicates a small
    number of guides amplified far beyond the rest. A distinct step or plateau in the middle of a
    panel's curve indicates two sub-populations of guides at different typical abundance in that
    pool — the same signal the Density panel below shows as two peaks.)",
    r"(<li>Cross-check that pool's % Undetected and 90/10 ratio in the headline metrics above.</li>
    <li>A drop already present at T0/plasmid indicates a library or transduction coverage problem
    specific to that pool.</li>)"),

  section("ecdf", "3. ECDF of gRNA Representation", "Each panel shows one pool's cumulative fraction of guides at or below a given read count. The dashed line marks that pool's mean read count — what every guide would have if reads were split perfectly evenly.",
    plot_div(svg_ecdf),
    "The more of a panel's curve sitting to the left of its mean line, the larger the fraction of that pool's guides falling below their fair share of reads.",
    r"(<li>A wide gap between a pool's curve and its mean line indicates that pool needs closer
    attention — cross-check its mapping rate and % undetected in the headline metrics above.</li>)"),

  section("lorenz", "4. Lorenz Curves", "Each panel shows cumulative % of reads against cumulative % of guides for one pool. The diagonal represents an ideal pool with perfectly even representation.",
    plot_div(svg_lorenz),
    "The further a panel's curve bows away from the diagonal, the less evenly that pool's guides are represented — that pool's Gini coefficient in the headline metrics above will be higher too.",
    "<li>Confirm that pool's 90/10 ratio above agrees with what the curve shows.</li>"),

  section("density", "5. Read-Count Density", "Each panel shows the read-count distribution for one pool.",
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
    section("controls", "6. Control vs Targeting gRNA Representation", "For pools with control/NTC/non-targeting guides: control vs targeting read-count distributions, one panel per pool.",
      plot_div(svg_controls),
      "Divergence between a panel's two curves already present at a reference/T0 timepoint indicates a technical difference between how that pool's controls and targeting guides were made, not biology.",
      "<li>Compare a pool's own T0/plasmid panel against its treated timepoint, if this run has both.</li>")
  } else {
    section("controls", "6. Control vs Targeting gRNA Representation", "For pools with control/NTC/non-targeting guides: control vs targeting read-count distributions, one panel per pool.",
      missing_div("None of this run's libraries have control, NTC, or non-targeting guides in their Gene column."),
      "", "")
  },

  if (!is.null(svg_class)) {
    section("class", "7. Read Classification: Exact / Near-Miss / Far-Miss", "Every pool's reads, split into exact matches (= that pool's Mapping Rate above), near-misses (1 substitution from a guide), and far-misses (not close to any guide), one panel per pool.",
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
    section("class", "7. Read Classification: Exact / Near-Miss / Far-Miss", "Every pool's reads, split into exact matches, near-misses, and far-misses, one panel per pool.",
      missing_div("None of this run's libraries had a usable combination of library sequences, mageck log, and unmapped-reads file to classify reads with."),
      "", "")
  },

  if (!is.null(svg_gc)) {
    section("gc", "8. GC Content vs Read Depth", "Guide GC content vs. normalized read count, one panel per pool.",
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
    section("gc", "8. GC Content vs Read Depth", "Guide GC content vs. normalized read count, one panel per pool.",
      missing_div("None of this run's libraries had the guide id/sequence columns needed to compute GC content."),
      "", "")
  }
), collapse = "\n")

toc_ids    <- c(if (!is.null(breakdown_html)) "breakdown", "headline", "guide_representation", "ecdf", "lorenz", "density", "controls", "class", "gc")
toc_titles <- c(if (!is.null(breakdown_html)) "0. Reads per Library", "1. Headline Metrics", "2. Guide Representation", "3. ECDF", "4. Lorenz Curves", "5. Read-Count Density",
                "6. Control vs Targeting", "7. Read Classification", "8. GC Content vs Read Depth")
toc_items  <- paste0("<li><a href='#", toc_ids, "'>", toc_titles, "</a></li>", collapse = "\n")

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
  "<section id='summary'><h2>At a Glance</h2>", summary_table,
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
