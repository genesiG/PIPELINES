#!/usr/bin/env Rscript
# heatmap_plot_enriched.R -- plot step of the heatmap module (PLOT_ENRICHED process, container cutandrun-r_1.1.sif).
# Merge of the legacy helpers Scripts/utils/plot_enriched_heatmaps.R (non-clustered) and
# plot_enriched_heatmaps_clustered.R (cluster row-split). Reads the objects written by heatmap_build_matrix.R and
# recomputes only colour ranges, profile limits and row order. No JSON.
#
# CLI contract (--key=value; comma lists):
#   --norm-mats=<NAME.norm_mats.rds>     named list of normalizedMatrix objects (names = sample ids)
#   --targets=<NAME.targets.rds>         targets GRanges (used for the row count check)
#   --window=W                           half-window in bp (axis labels: -Wkb / center / +Wkb; bp when W < 1000)
#   --out-prefix=NAME                    output prefix; file written to the working directory
#   --names=id1,id2                      sample ids to plot, in column order (keys of the norm_mats list)
#   --labels="a,b"                       display label per sample (surrounding quotes stripped, literal \n -> newline)
#   --colors="#hex,#hex"                 one hue per sample; colour scale white -> hue
#   --legend-titles="t1,t2"              optional legend title per sample (1 entry = all samples; n entries = one each;
#                                        empty = the display label). Passed as text, never hard-coded.
#   --z-min=, --z-max=                   optional colour range: empty | one value | one per sample.
#                                        z-min default 0; z-max default = max column mean of the sample (data-driven)
#   --y-max=                             optional profile annotation ylim: empty | one value | one per sample.
#                                        Independent of z-max. Default = 1.08 x max over row groups of the column means
#   --assignments=<NAME.cluster_assignments.csv>  one integer per matrix row (1..k), no header. Row count must equal the
#                                        matrix rows (exit 3). Empty/absent = non-clustered plot.
#   --cluster-labels="a,b"               k labels; used only if the count == k, else "cluster 1".."cluster k" (warning)
#   --sort-by=<sample id>                sample whose row means order the rows within a cluster (default: first of --names);
#                                        without --assignments it is ignored (warning)
#   --format=svg                         optional; only svg is written (default svg)
#   --heatmap-width-cm=3.5               width per heatmap column (cm)
#   --column-gap-cm=1                    gap between columns (cm)
#   --anno-height-cm=2.5                 profile annotation height (cm)
#
# Row order: without --assignments, decreasing mean signal over all plotted samples.
# With --assignments: row_split by cluster (levels in numeric order), decreasing mean of the --sort-by sample within a cluster.
# Cluster profile colours: RColorBrewer Set1 (k <= 9; recycled beyond 9). Clustered legends add a "Cluster profile" key.
# Output: NAME.enriched_heatmap.svg (svg device) in the working directory.
# Exit codes: 0 ok; 2 usage (bad or missing argument); 3 data problem (missing file, length mismatch, row count mismatch,
#   draw() failure).

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(EnrichedHeatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
})

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
# Display text: strip surrounding quotes, literal backslash-n -> newline
clean_label <- function(x) {
  x <- gsub("^[\"']+|[\"']+$", "", x)
  gsub("\\n", "\n", x, fixed = TRUE)
}
# Numeric per-sample list: empty -> all NA (data-driven); one value -> recycled; n values -> one each
num_list <- function(args, key, n) {
  raw <- get_opt(args, key, "")
  if (!nzchar(raw)) return(rep(NA_real_, n))
  toks <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
  x    <- suppressWarnings(as.numeric(toks))
  if (any(is.na(x) & nzchar(toks))) die(paste0("--", key, " must hold numbers (comma list)"), 2L)
  if (length(x) == 1L) x <- rep(x, n)
  if (length(x) != n) die(sprintf("--%s must have 1 or %d entries (got %d)", key, n, length(x)), 3L)
  x
}
num_opt <- function(args, key, default) {
  v <- get_opt(args, key, NULL)
  if (is.null(v)) return(default)
  x <- suppressWarnings(as.numeric(v))
  if (is.na(x)) die(paste0("--", key, " must be numeric"), 2L)
  x
}

allowed <- c("norm-mats", "targets", "window", "out-prefix", "names", "labels", "colors", "legend-titles",
             "z-min", "z-max", "y-max", "assignments", "cluster-labels", "sort-by", "format",
             "heatmap-width-cm", "column-gap-cm", "anno-height-cm")
args <- parse_args(commandArgs(trailingOnly = TRUE), allowed)

norm_file   <- get_req(args, "norm-mats")
tgt_file    <- get_req(args, "targets")
out_prefix  <- get_req(args, "out-prefix")
sample_ids  <- split_list(get_req(args, "names"))
labels_raw  <- split_list(get_req(args, "labels"))
colors      <- split_list(get_req(args, "colors"))
window_bp   <- suppressWarnings(as.integer(get_req(args, "window")))
if (is.na(window_bp) || window_bp < 1L) die("--window must be a positive integer (bp)", 2L)
out_format  <- tolower(get_opt(args, "format", "svg"))
if (out_format != "svg") die("--format must be svg", 2L)
hw_cm       <- num_opt(args, "heatmap-width-cm", 3.5)
gap_cm      <- num_opt(args, "column-gap-cm", 1)
anno_cm     <- num_opt(args, "anno-height-cm", 2.5)

n <- length(sample_ids)
if (n < 1L) die("--names is empty", 2L)
if (anyDuplicated(sample_ids))
  die(paste0("duplicate --names: ", paste(unique(sample_ids[duplicated(sample_ids)]), collapse = ", ")), 3L)
if (length(labels_raw) != n)
  die(sprintf("--labels (%d) and --names (%d) must have the same count", length(labels_raw), n), 3L)
if (length(colors) != n)
  die(sprintf("--colors (%d) and --names (%d) must have the same count", length(colors), n), 3L)

for (f in c(norm_file, tgt_file))
  if (!file.exists(f)) die(paste0("file not found: ", f), 3L)

norm_all <- tryCatch(readRDS(norm_file),
                     error = function(e) die(paste0("cannot read norm_mats: ", conditionMessage(e)), 3L))
targets  <- tryCatch(readRDS(tgt_file),
                     error = function(e) die(paste0("cannot read targets: ", conditionMessage(e)), 3L))

missing_s <- setdiff(sample_ids, names(norm_all))
if (length(missing_s) > 0L)
  die(paste0("sample ids not found in norm_mats: ", paste(missing_s, collapse = ", ")), 3L)
norm_mats <- norm_all[sample_ids]
n_rows    <- length(targets)
for (s in sample_ids) {
  if (nrow(norm_mats[[s]]) != n_rows)
    die(sprintf("norm_mats '%s' has %d rows but targets has %d", s, nrow(norm_mats[[s]]), n_rows), 3L)
}

cat("=== heatmap_plot_enriched.R ===\n")
cat(sprintf("  Prefix  : %s\n", out_prefix))
cat(sprintf("  Window  : +/- %d bp\n", window_bp))
cat(sprintf("  Samples : %s\n", paste(sample_ids, collapse = ", ")))

# -----------------------------------------------------------------------------
# Labels, legend titles, per-sample scales
# -----------------------------------------------------------------------------
disp <- vapply(labels_raw, clean_label, character(1), USE.NAMES = FALSE)

lt_raw <- split_list(get_opt(args, "legend-titles", ""))
if (length(lt_raw) == 0L) {
  legend_titles <- disp
} else if (length(lt_raw) == 1L) {
  legend_titles <- rep(clean_label(lt_raw), n)
} else if (length(lt_raw) == n) {
  legend_titles <- vapply(lt_raw, clean_label, character(1), USE.NAMES = FALSE)
} else {
  die(sprintf("--legend-titles must have 1 or %d entries (got %d)", n, length(lt_raw)), 3L)
}

z_min_v <- num_list(args, "z-min", n)
z_max_v <- num_list(args, "z-max", n)
y_max_v <- num_list(args, "y-max", n)

window_kb   <- window_bp / 1000
axis_labels <- if (window_bp >= 1000L) {
  c(paste0("-", window_kb, "kb"), "center", paste0("+", window_kb, "kb"))
} else {
  c(paste0("-", window_bp, "bp"), "center", paste0("+", window_bp, "bp"))
}

mat_list <- lapply(sample_ids, function(s) as.matrix(norm_mats[[s]]))

sample_cols <- vector("list", n)
for (i in seq_len(n)) {
  zmin <- if (is.na(z_min_v[i])) 0 else z_min_v[i]
  zmax <- if (is.na(z_max_v[i])) max(colMeans(mat_list[[i]], na.rm = TRUE), 1e-9) else z_max_v[i]
  sample_cols[[i]] <- colorRamp2(c(zmin, zmax), c("white", colors[i]))
  cat(sprintf("  Color scale [%s]: %.3f - %.3f  hue=%s\n", disp[i], zmin, zmax, colors[i]))
}

# -----------------------------------------------------------------------------
# Clusters (optional), row order and row split
# -----------------------------------------------------------------------------
assign_file <- get_opt(args, "assignments", "")
clustered   <- nzchar(assign_file)
cl          <- NULL
cl_names    <- NULL
sort_by     <- get_opt(args, "sort-by", sample_ids[1])
if (!(sort_by %in% sample_ids))
  die(paste0("--sort-by '", sort_by, "' is not one of --names"), 3L)

if (clustered) {
  if (!file.exists(assign_file)) die(paste0("file not found: ", assign_file), 3L)
  txt <- trimws(readLines(assign_file, warn = FALSE))
  txt <- txt[nzchar(txt)]
  cl  <- suppressWarnings(as.integer(txt))
  if (length(cl) == 0L || anyNA(cl) || any(cl < 1L))
    die("--assignments must hold one integer >= 1 per row", 3L)
  if (length(cl) != n_rows)
    die(sprintf("--assignments has %d rows but the matrices have %d rows", length(cl), n_rows), 3L)
  k <- max(cl)

  cl_raw <- split_list(get_opt(args, "cluster-labels", ""))
  if (length(cl_raw) == k) {
    cl_names <- cl_raw
  } else {
    if (length(cl_raw) > 0L)
      message(sprintf("WARNING: --cluster-labels has %d entries but k = %d; using 'cluster 1'..'cluster %d'",
                      length(cl_raw), k, k))
    cl_names <- paste0("cluster ", seq_len(k))
  }
  cat(sprintf("  k clusters : %d  (%s)\n", k, paste(cl_names, collapse = " | ")))
  cat(sprintf("  Sort by    : %s\n", sort_by))

  sort_mat  <- mat_list[[match(sort_by, sample_ids)]]
  row_order_vec <- order(cl, -rowMeans(sort_mat, na.rm = TRUE))
  row_split     <- factor(cl_names[cl], levels = cl_names)

  pal        <- brewer.pal(max(3L, min(k, 9L)), "Set1")
  prof_cols  <- pal[(seq_len(k) - 1L) %% length(pal) + 1L]
  group_ids  <- cl
} else {
  if (nzchar(get_opt(args, "sort-by", ""))) message("WARNING: --sort-by ignored without --assignments")
  combined_mat  <- do.call(cbind, mat_list)
  row_order_vec <- order(-rowMeans(combined_mat, na.rm = TRUE))
  row_split     <- NULL
  group_ids     <- rep(1L, n_rows)
}

# Profile ylim: per-group worst-case max column mean, 8% headroom (legacy rule)
profile_ymax <- function(m) {
  mx <- vapply(sort(unique(group_ids)),
               function(g) max(colMeans(m[group_ids == g, , drop = FALSE], na.rm = TRUE)),
               numeric(1))
  max(mx) * 1.08
}
profile_ylims <- lapply(seq_len(n), function(i) {
  ymax <- if (is.na(y_max_v[i])) profile_ymax(mat_list[[i]]) else y_max_v[i]
  c(0, max(ymax, 1e-9))
})

# -----------------------------------------------------------------------------
# Legends (packed in one horizontal row at the bottom)
# -----------------------------------------------------------------------------
font_legend_title <- 9
font_legend_text  <- 8
font_col_title    <- 11
font_axis_name    <- 8
font_anno_axis    <- 7
font_row_title    <- 9
profile_lwd       <- 1.5

lgd_list <- lapply(seq_len(n), function(i) {
  Legend(
    col_fun    = sample_cols[[i]],
    title      = legend_titles[i],
    direction  = "horizontal",
    title_gp   = gpar(fontsize = font_legend_title, fontface = "bold"),
    labels_gp  = gpar(fontsize = font_legend_text)
  )
})
if (clustered) {
  lgd_cluster <- Legend(
    title     = "Cluster profile",
    labels    = cl_names,
    legend_gp = gpar(col = prof_cols, lwd = profile_lwd),
    type      = "lines",
    direction = "horizontal",
    title_gp  = gpar(fontsize = font_legend_title, fontface = "bold"),
    labels_gp = gpar(fontsize = font_legend_text)
  )
  lgd_list <- c(lgd_list, list(lgd_cluster))
}
combined_legend <- do.call(packLegend,
                           c(lgd_list, list(direction = "horizontal", gap = unit(4, "mm"))))

# -----------------------------------------------------------------------------
# Heatmap list
# -----------------------------------------------------------------------------
ht_list <- NULL
for (i in seq_len(n)) {
  s <- sample_ids[i]
  anno_gp <- if (clustered) gpar(col = prof_cols, lwd = profile_lwd) else gpar(col = colors[i], lwd = profile_lwd)

  top_anno <- HeatmapAnnotation(
    enriched = anno_enriched(
      gp         = anno_gp,
      ylim       = profile_ylims[[i]],
      axis_param = list(gp = gpar(fontsize = font_anno_axis))
    ),
    height = unit(anno_cm, "cm")
  )

  ht <- EnrichedHeatmap(
    norm_mats[[s]],
    col                 = sample_cols[[i]],
    name                = s,
    width               = unit(hw_cm, "cm"),
    column_title        = disp[i],
    column_title_gp     = gpar(fontsize = font_col_title, fontface = "bold"),
    axis_name           = axis_labels,
    axis_name_gp        = gpar(fontsize = font_axis_name),
    row_order           = row_order_vec,
    row_split           = row_split,
    row_title_side      = "left",
    row_title_rot       = 90,
    row_title_gp        = gpar(fontsize = font_row_title, fontface = "bold"),
    row_gap             = if (clustered) unit(2, "mm") else unit(0, "mm"),
    border              = clustered,
    top_annotation      = top_anno,
    pos_line            = FALSE,
    use_raster          = TRUE,
    show_heatmap_legend = FALSE
  )
  ht_list <- if (is.null(ht_list)) ht else ht_list + ht
}

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
out_file       <- paste0(out_prefix, ".enriched_heatmap.svg")
total_width_cm <- (n * hw_cm) + ((max(1L, n) - 1L) * gap_cm) + 5.0
svg_w          <- total_width_cm / 2.54
svg_h          <- 10

cat(sprintf("  Saving heatmap: %s (width=%.2f in, height=%.2f in)\n",
            basename(out_file), svg_w, svg_h))
svg(out_file, width = svg_w, height = svg_h)
invisible(tryCatch(
  draw(ht_list,
       gap                    = unit(gap_cm, "cm"),
       annotation_legend_list = list(combined_legend),
       annotation_legend_side = "bottom",
       merge_legend           = FALSE,
       padding                = unit(c(5, 5, 5, 5), "mm")),
  error = function(e) {
    invisible(try(dev.off(), silent = TRUE))
    die(paste0("draw() failed: ", conditionMessage(e)), 3L)
  }
))
invisible(dev.off())

cat(sprintf("  Saved: %s\n", basename(out_file)))
cat(sprintf("=== heatmap_plot_enriched.R complete (%s) ===\n\n", out_prefix))
