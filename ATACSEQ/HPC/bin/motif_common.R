#!/usr/bin/env Rscript
# motif_common.R -- motifs that are significantly enriched in several comparisons (one motif table per comparison, one region
# set). Port of the legacy find_common_motifs.Rmd (Scripts/, AV_ATACSEQ_2026, March 2026): chunks "load cell lines",
# "merge data frames and summarize motifs", "Find common motifs" and "venn diagram of significant motifs". The computation
# (the reshape() tables, the significance flag, the >= min-lines count, the Venn sets) is the legacy code, verbatim, except that the
# key is motif.id (deviation, see the header below); only the input discovery (a file list instead of list.files on motif/) and
# the output paths are replaced by the CLI.
#
# CLI contract (called by the MOTIF_COMMON module; --key=value):
#   --inputs=a_motif_analysis.csv,b_motif_analysis.csv   REQUIRED. Diff-mode motif tables of ONE region set (MOTIF_ENRICHMENT,
#                                 NAME_motif_analysis.csv). Columns used: motif.id, motif.name, padj_up, log2enr_up, padj_down,
#                                 log2enr_down. Missing file: exit 3.
#   --labels=lab1,lab2                                  REQUIRED. One label per input, same order, unique. The legacy label is the
#                                 comparison directory name (e.g. BT54_high_vs_low) with the --label-strip regex removed.
#   --label-strip=_high_vs_low.*  OPTIONAL. Regex removed from every label with gsub (legacy: "_high_vs_low.*" -> "BT54").
#                                 Blank = labels used as given.
#   --padj=0.10                   OPTIONAL, default 0.10. A motif is significant in a comparison when padj_up >= -log10(padj),
#                                 i.e. adjusted P <= padj (legacy padj_thres).
#   --min-lines=2                 OPTIONAL, default 2. Common motifs: significant in at least this many comparisons.
#   --summaries=a.csv,b.csv       OPTIONAL. MOTIF_ENRICHMENT NAME.motif_summary.csv files of the same region set. Their rows are
#                                 stacked (dplyr::bind_rows, so differing columns are allowed) into motif_analysis_summary.csv.
#   --out-prefix=.                OPTIONAL, default "." (current directory). Directory the three outputs are written to.
#   --font=Helvetica              OPTIONAL, default Helvetica. Font family of the plot (svg: font-family of the text elements).
#   --plot-format=svg|pdf         OPTIONAL, default svg. svg: venn_common_motifs.svg (svglite, editable <text>). pdf: the legacy
#                                 venn_common_motifs.pdf (ggsave, as before).
#   --venn-width=8                OPTIONAL, default 8 (inches). Width of the Venn diagram.
#   --venn-height=7               OPTIONAL, default 7 (inches). Height of the Venn diagram.
#
# Outputs (in --out-prefix):
#   common_significant_motifs.csv  write.csv(row.names = FALSE). One row per motif significant in >= min-lines comparisons,
#                                  sorted by the count (descending). Columns: motif.id, motif.name, label (the motif.name, or
#                                  "motif.name (motif.id)" when several matrices share the name), one TRUE/FALSE column per label
#                                  (significant), n_cell_lines_significant (legacy name), wide columns log2enr_up.<label> and
#                                  padj_up.<label> per label, avg_log2enr, max_log2enr, max_padj.
#   venn_common_motifs.svg         Venn diagram of the significant motifs per label (--plot-format=pdf: venn_common_motifs.pdf).
#                                  Written ONLY when at least one motif is significant in all inputs (legacy); otherwise not written.
#   motif_analysis_summary.csv     rbind of the --summaries files (only when --summaries is given)
#
# Exit codes: 0 ok; 2 usage error; 3 known data problem (an input file is missing or has none of the required columns).
# With fewer than 2 inputs nothing is written: a message is printed and the exit code is 0 (legacy: needs several cell lines).
#
# Question, model, assumptions: in how many comparisons is a motif significantly enriched in the increased-accessibility bin?
# A motif is called significant in a comparison when its adjusted P value (negLog10, from motif_enrichment.R: binomial test
# with composition control) passes --padj. The comparisons are then counted; no new test is run, so the count is a set-overlap
# statistic, not a test of common enrichment. The motifs are matched by motif.id (the JASPAR matrix ID).
#
# Not ported: the legacy printouts of top-10 tables and the Venn print() to the console device (no file effect).
#
# ---- Statistical limits of the legacy method, ported unchanged (state them when reporting results) ---------------------------
# * Significance uses padj_up only. monaLisa's binomial test is one-sided (enrichment), so padj_up is a sign-free enrichment
#   score of the Increased bin versus the Decreased bin and the call "significant" means "enriched in the Increased bin of that
#   comparison"; motifs enriched in the Decreased bin (padj_down) are not counted. See the limits listed in motif_enrichment.R
#   (bins cut on logFC without an FDR filter, hits treated as independent, anti-conservative p-values).
# * The count of comparisons in which a motif is significant is a set-overlap descriptor, not a test of common enrichment, and
#   the thresholds (padj <= --padj) are applied to p-values that are anti-conservative. Comparisons that rest on one library per
#   side (fixed-dispersion LRT) should not be mixed with replicated ones without saying so.
# * motif.name is not unique in JASPAR2024 (84 of 720 PWMs share one of 39 names, e.g. CTCF, FOS, ASCL1), so the script keys
#   on motif.id; the label column adds the id to names shared by several matrices.
#
options(stringsAsFactors = FALSE)

die <- function(msg, status = 3L) {
  cat("ERROR: ", msg, "\n", sep = "", file = stderr())
  quit(save = "no", status = status)
}
parse_args <- function(argv, allowed) {
  out <- list()
  for (a in argv) {
    m <- regmatches(a, regexec("^--([^=]+)=(.*)$", a))[[1]]
    if (length(m) != 3L) die(paste0("cannot parse argument '", a, "' (expected --key=value)"), 2L)
    if (!(m[2] %in% allowed)) die(paste0("unknown argument --", m[2]), 2L)
    out[[m[2]]] <- m[3]
  }
  out
}
get_req <- function(args, key) {
  v <- args[[key]]
  if (is.null(v) || !nzchar(v)) die(paste0("missing required argument --", key), 2L)
  v
}
get_opt <- function(args, key, default) {
  v <- args[[key]]
  if (is.null(v) || !nzchar(v)) default else v
}
split_list <- function(x) { v <- trimws(strsplit(x, ",", fixed = TRUE)[[1]]); v[nzchar(v)] }
as_num <- function(x, key) {
  v <- suppressWarnings(as.numeric(x))
  if (length(v) != 1L || is.na(v)) die(paste0("--", key, " must be a number, got '", x, "'"), 2L)
  v
}

# ---- plots: one svg per plot (svglite, editable <text>, font family = --font); --plot-format=pdf = the legacy pdf ----
# open_svg / close_svg are copied from the RNASEQ module (bin/deg_common.R). close_svg rewrites the family names that
# svglite resolves through fontconfig (e.g. "Nimbus Sans" for Helvetica) to the requested --font; widths are kept.
open_svg <- function(path, width, height) {
  svglite::svglite(filename = path, width = width, height = height,
    system_fonts = list(sans = plot_font, serif = plot_font, mono = plot_font, symbol = plot_font))
  invisible(path)
}
close_svg <- function(path) {
  grDevices::dev.off()
  txt <- readLines(path, warn = FALSE, encoding = "UTF-8")
  txt <- gsub("font-family: ('[^']*'|\"[^\"]*\"|[^;]*);", paste0("font-family: ", plot_font, ";"), txt)
  writeLines(txt, path, useBytes = TRUE)
  invisible(NULL)
}
# One plot written to BASE.svg (--plot-format=svg). Pair with plot_close(BASE). The pdf branch uses ggsave (legacy).
plot_open <- function(base, width, height) {
  open_svg(paste0(base, ".svg"), width, height)
  invisible(base)
}
plot_close <- function(base) {
  close_svg(paste0(base, ".svg"))
  invisible(NULL)
}

allowed <- c("inputs", "labels", "label-strip", "padj", "min-lines", "summaries", "out-prefix",
             "font", "plot-format", "venn-width", "venn-height")
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) == 0L || any(argv %in% c("-h", "--help"))) {
  cat("usage: motif_common.R --inputs=a.csv,b.csv --labels=a,b [--label-strip=REGEX] [--padj=0.10] [--min-lines=2]\n",
      "       [--summaries=a.csv,b.csv] [--out-prefix=.] [--font=Helvetica] [--plot-format=svg|pdf]\n",
      "       [--venn-width=8] [--venn-height=7]\n", file = stderr())
  quit(save = "no", status = 2L)
}
args <- parse_args(argv, allowed)

results    <- split_list(get_req(args, "inputs"))
labels     <- split_list(get_req(args, "labels"))
label_strip <- get_opt(args, "label-strip", "")
padj_thres <- as_num(get_opt(args, "padj", "0.10"), "padj")
min_lines  <- as.integer(as_num(get_opt(args, "min-lines", "2"), "min-lines"))
summaries  <- split_list(get_opt(args, "summaries", ""))
out_prefix <- get_opt(args, "out-prefix", ".")
plot_font   <- get_opt(args, "font", "Helvetica")
plot_format <- tolower(get_opt(args, "plot-format", "svg"))
venn_width  <- as_num(get_opt(args, "venn-width", "8"), "venn-width")
venn_height <- as_num(get_opt(args, "venn-height", "7"), "venn-height")

if (length(results) != length(labels)) die(paste0("--inputs (", length(results), ") and --labels (", length(labels),
                                                  ") must have the same length"), 2L)
if (anyDuplicated(labels)) die("--labels must be unique", 2L)
if (!(plot_format %in% c("svg", "pdf"))) die("--plot-format must be svg or pdf", 2L)
if (venn_width <= 0 || venn_height <= 0) die("--venn-width and --venn-height must be > 0", 2L)
if (length(results) < 2L) {
  cat("Fewer than 2 motif tables: nothing to compare, no output written.\n")
  quit(save = "no", status = 0L)
}
if (!all(file.exists(results))) die(paste0("motif table not found: ", results[!file.exists(results)][1]), 3L)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggVennDiagram)
})

if (nzchar(label_strip)) labels <- gsub(label_strip, "", labels)
if (anyDuplicated(labels)) {
  die(paste0("the labels are not unique after --label-strip (", paste(unique(labels[duplicated(labels)]), collapse = ", "),
             "): two comparisons would be merged into one column"), 2L)
}
cell_lines <- labels

# ── 1. Load each results table and tag it with its cell line name ──────────────

# Keep only the columns we need
motif_list <- setNames(
  lapply(seq_along(results), function(i) {
    df <- read.csv(results[i], stringsAsFactors = FALSE)
    needed <- c("motif.id", "motif.name", "padj_up", "log2enr_up", "padj_down", "log2enr_down")
    if (!all(needed %in% colnames(df))) die(paste0("missing columns in ", results[i], ": ",
                                                   paste(setdiff(needed, colnames(df)), collapse = ", ")), 3L)

    # Retain motif identifiers and the significance score for the "up" bin
    df <- df[, c("motif.id", "motif.name",
                 "padj_up", "log2enr_up",
                 "padj_down", "log2enr_down")]

    # Add a column recording which cell line this row came from
    df$cell_line <- cell_lines[i]

    return(df)
  }),
  cell_lines
)

cat("Loaded results for cell lines:", paste(cell_lines, collapse = ", "), "\n")

# ── 2. Flag significant motifs in each cell line ───────────────────────────────

# A motif is considered significant when its adjusted P value passes the threshold
motif_list_sig <- lapply(motif_list, function(df) {
  df$significant <- NA
  df$significant <- df$padj_up >= -log10(padj_thres)
  return(df)
})

# ── 3. Build a wide summary: one row per motif, one column per cell line ────────

# Combine all data frames into one long table — must be done before reshaping
long_df <- do.call(rbind, motif_list_sig)

# The key is motif.id. id_map gives the motif.name of every id in the inputs; label is the name, with the id added when
# the name belongs to more than one matrix (motif.name is not unique in JASPAR2024)
id_map <- unique(do.call(rbind, lapply(motif_list, function(df) df[, c("motif.id", "motif.name")])))
id_map <- id_map[!duplicated(id_map$motif.id), ]
shared_names <- id_map$motif.name[duplicated(id_map$motif.name)]
id_map$label <- ifelse(id_map$motif.name %in% shared_names,
                       paste0(id_map$motif.name, " (", id_map$motif.id, ")"), id_map$motif.name)

# --- Significance wide table (TRUE/FALSE per cell line) ---
# Using base R reshape to avoid tidyr version incompatibilities with values_fill

wide_sig <- reshape(
  data      = long_df[, c("motif.id", "cell_line", "significant")],
  idvar     = "motif.id",
  timevar   = "cell_line",
  direction = "wide"
)
# reshape() prefixes column names with "significant." — strip that prefix
colnames(wide_sig) <- gsub("^significant\\.", "", colnames(wide_sig))

# Fill any NA (motif not tested in a given cell line) with FALSE
wide_sig[is.na(wide_sig)] <- FALSE

# --- log2enr wide table (numeric per cell line) ---
log2_wide <- reshape(
  data      = long_df[, c("motif.id", "cell_line",
                          "log2enr_up", "padj_up")],
  idvar     = "motif.id",
  timevar   = "cell_line",
  direction = "wide"
)

print(head(log2_wide))

# ── 4. Find common significant motifs ──

# Fill NA (motif not tested) with NA_real_ — kept as NA rather than 0 to avoid
# implying the motif was tested and had zero enrichment
log2_wide[is.na(log2_wide)] <- NA_real_

# Get column indexes for log enrichment data
log2_columns <- which(grepl("log2enr", colnames(log2_wide)))

# Compute average enrichment for each motif
log2_wide$avg_log2enr <- rowMeans(log2_wide[,log2_columns], na.rm = TRUE)

# Use apply to find the max of each row (MARGIN = 1 means rows)
log2_wide$max_log2enr <- apply(log2_wide[,log2_columns],
                               MARGIN = 1,
                               FUN = max, na.rm = TRUE)

# Get column indexes for significance data
padj_columns <- which(grepl("padj", colnames(log2_wide)))

# Use apply to find the max significance of each row
log2_wide$max_padj <- apply(log2_wide[,padj_columns],
                               MARGIN = 1,
                               FUN = max, na.rm = TRUE)

print(head(log2_wide))

# Count in how many cell lines each motif was significant
wide_sig$n_cell_lines_significant <- rowSums(wide_sig[, cell_lines], na.rm = TRUE)

# Keep only motifs significant in at least min_lines cell lines, sorted by breadth
common_motifs <- wide_sig[wide_sig$n_cell_lines_significant >= min_lines, ]
common_motifs <- common_motifs[order(-common_motifs$n_cell_lines_significant), ]

common_motifs <- merge(common_motifs, log2_wide, by = "motif.id")
common_motifs <- common_motifs[order(-common_motifs$n_cell_lines_significant), ]

# Add the names (motif.name and label) after the id, then the per-label columns in their order
common_motifs$motif.name <- id_map$motif.name[match(common_motifs$motif.id, id_map$motif.id)]
common_motifs$label      <- id_map$label[match(common_motifs$motif.id, id_map$motif.id)]
common_motifs <- common_motifs[, c("motif.id", "motif.name", "label",
                                   setdiff(colnames(common_motifs), c("motif.id", "motif.name", "label")))]

# Print a compact preview to the console
cat("\nMotifs significant in >= ", min_lines, " cell lines: ", nrow(common_motifs), "\n", sep = "")

# Export the full table
write.csv(
  x         = common_motifs,
  file      = file.path(out_prefix, "common_significant_motifs.csv"),
  row.names = FALSE
)
cat("\nFull table written to: common_significant_motifs.csv\n")

# ── 5. Venn diagram — only when at least one motif is significant in ALL cell lines ──

# Build a named list of significant motif names per cell line (for the Venn input)
sig_sets <- lapply(motif_list_sig, function(df) {
  df$motif.id[df$significant]
})

# Check whether the intersection across all cell lines is non-empty
universal_motifs <- Reduce(intersect, sig_sets)
universal_motifs <- na.exclude(universal_motifs)
cat("\nMotifs significant in ALL cell lines:", length(universal_motifs), "\n")

if (length(universal_motifs) > 0) {
  cat("\n They are:\n")
  print(universal_motifs)

  venn_plot <- ggVennDiagram(
    x            = sig_sets,
    label_alpha  = 0,          # transparent label background
    label        = "count"     # show counts in each region
  ) +
    scale_fill_gradient(low = "white", high = "#2166AC") +
    scale_color_manual(values = rep("grey30", length(sig_sets))) +
    ggtitle(
      label    = "Shared significant Motifs Across Cell Lines",
      subtitle = paste0("Threshold: FDR < ", padj_thres,
                        "  |  Motifs in all lines: ", length(universal_motifs))
    ) +
    theme(plot.title    = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5))

  # Save to disk: svg (svglite, Helvetica text) or, with --plot-format=pdf, the legacy ggsave pdf
  if (plot_format == "pdf") {
    ggsave(filename = file.path(out_prefix, "venn_common_motifs.pdf"),
           plot     = venn_plot,
           width    = venn_width,
           height   = venn_height)
  } else {
    venn_base <- file.path(out_prefix, "venn_common_motifs")
    plot_open(venn_base, venn_width, venn_height)
    print(venn_plot)
    plot_close(venn_base)
  }

  cat(paste0("Venn diagram written to: venn_common_motifs.", plot_format, "\n"))

} else {
  cat("No motif was significant in every cell line — Venn diagram skipped.\n")
}

# ── 6. Stack the per-comparison summaries of this region set ──

if (length(summaries) > 0) {
  summary_df <- dplyr::bind_rows(lapply(summaries, read.csv, stringsAsFactors = FALSE))
  write.csv(summary_df, file = file.path(out_prefix, "motif_analysis_summary.csv"), row.names = FALSE)
  cat("Summary file written to: motif_analysis_summary.csv\n")
}
