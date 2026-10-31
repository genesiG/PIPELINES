#!/usr/bin/env Rscript
# qc_preseq_plot.R -- preseq c_curve complexity curves, one line per sample.
#
# Called by QC_PRESEQ_PLOT (modules/local/qc_preseq_plot.nf):
#   qc_preseq_plot.R --curves-dir=. --theme=my_theme.R
#
# Arguments (--key=value):
#   --curves-dir  directory with <sample_id>.preseq_ccurve.txt files (preseq c_curve output: total_reads, distinct_reads)
#   --theme       my_theme.R (palette); a minimal built-in palette is used when absent
#   --outdir      default "."
# Output: preseq_ccurve.svg
#   x = fragments sequenced (millions), y = distinct fragments (millions), one line per sample, the dashed diagonal
#   marks a library without any duplicate (distinct = total). Colours cycle through the theme palette; when there
#   are more samples than palette colours the line type changes too, so each sample stays distinguishable.
#   The sample name is the file name without ".preseq_ccurve.txt" (the plot has no samplesheet input).
# Exit 3: no curve file, or none parsable.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))
parse_args <- function(args) {
  out <- list()
  for (m in regmatches(args, regexec("^--([^=]+)=(.*)$", args))) if (length(m) == 3) out[[m[2]]] <- m[3]
  out
}
die <- function(...) { message("ERROR: ", ...); quit(save = "no", status = 3) }
opt <- parse_args(commandArgs(trailingOnly = TRUE))
dir_in <- if (is.null(opt[["curves-dir"]])) "." else opt[["curves-dir"]]
outdir <- if (is.null(opt$outdir)) "." else opt$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
if (!is.null(opt$theme) && file.exists(opt$theme)) source(opt$theme)
base_theme <- if (exists("theme_paper")) theme_paper() else theme_classic(base_size = 11)
pal_get <- function(k, d) if (exists("pal") && is.list(pal) && !is.null(pal[[k]])) pal[[k]] else d
base_cols <- c(pal_get("purple", "#86187f"), pal_get("steel", "#6b93ad"), pal_get("slate", "#425b74"),
               pal_get("grey", "#a7b4c2"), pal_get("navy", "#1d3557"))

files <- sort(list.files(dir_in, pattern = "\\.preseq_ccurve\\.txt$", full.names = TRUE))
if (length(files) == 0) die("no *.preseq_ccurve.txt in ", dir_in)
read_curve <- function(f) {
  l <- readLines(f, warn = FALSE); l <- l[grepl("^[0-9]", trimws(l))]
  if (length(l) < 2) return(NULL)
  v <- do.call(rbind, strsplit(trimws(l), "[[:space:]]+"))
  data.frame(sample_id = sub("\\.preseq_ccurve\\.txt$", "", basename(f)),
             total = as.numeric(v[, 1]), distinct = as.numeric(v[, 2]))
}
cur <- Filter(Negate(is.null), lapply(files, function(f) tryCatch(read_curve(f), error = function(e) NULL)))
if (length(cur) == 0) die("no parsable preseq c_curve file")
d <- do.call(rbind, cur)
ids <- unique(d$sample_id)
d$sample_id <- factor(d$sample_id, levels = ids)
cols <- setNames(rep_len(base_cols, length(ids)), ids)
lts <- setNames(rep_len(c("solid", "dashed", "dotdash", "longdash", "twodash"), ceiling(length(ids) / length(base_cols)))[ceiling(seq_along(ids) / length(base_cols))], ids)
lim <- max(d$total, d$distinct, na.rm = TRUE)
p <- ggplot(d, aes(total / 1e6, distinct / 1e6, colour = sample_id, linetype = sample_id)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = pal_get("rule", "#9aa1a9"), linewidth = 0.6) +
  geom_line(linewidth = 0.9) +
  scale_colour_manual(values = cols, name = NULL) + scale_linetype_manual(values = lts, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.03))) + scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
  coord_cartesian(xlim = c(0, lim / 1e6), ylim = c(0, lim / 1e6)) +
  labs(x = "Fragments sequenced (millions)", y = "Distinct fragments (millions)") + base_theme
grDevices::svg(filename = file.path(outdir, "preseq_ccurve.svg"), width = 8, height = 6)
print(p)
invisible(dev.off())
cat("preseq plot written for", length(ids), "samples\n")
