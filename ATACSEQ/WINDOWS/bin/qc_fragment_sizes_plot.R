#!/usr/bin/env Rscript
# qc_fragment_sizes_plot.R -- fragment size plots, drawn only from the two fragment-size tables.
#
# Called by QC_FRAGMENT_SIZES_PLOT (modules/local/qc_fragment_sizes_plot.nf):
#   qc_fragment_sizes_plot.R --summary=fragment_sizes.tsv --distribution=fragment_size_distribution.tsv --theme=my_theme.R
#
# Arguments (--key=value): --summary, --distribution (tables from qc_fragment_sizes_table.R), --theme (my_theme.R;
#   minimal built-in theme when absent), --outdir (default "."), --xmax (x-axis limit of the density plot, default 1000)
# Colour group = the condition column; a blank condition is replaced by the sample_id, so each such sample gets its own colour.
# Outputs:
#   fragment_sizes_density.svg  one smoothed size-distribution line per sample, coloured by condition (colours cycle
#                               through the theme palette in order of appearance in the data); dashed guides at
#                               120/200/400/600 bp. Smoothing: Gaussian kernel, bandwidth 4 bp, applied to the exact counts.
#   fragment_sizes_boxplot.svg  median fragment size per sample, grouped by condition (blank: sample_id)

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))
parse_args <- function(args) {
  out <- list()
  for (m in regmatches(args, regexec("^--([^=]+)=(.*)$", args))) if (length(m) == 3) out[[m[2]]] <- m[3]
  out
}
die <- function(...) { message("ERROR: ", ...); quit(save = "no", status = 3) }
opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt$summary) || is.null(opt$distribution)) die("usage: qc_fragment_sizes_plot.R --summary=<tsv> --distribution=<tsv> [--theme=my_theme.R]")
outdir <- if (is.null(opt$outdir)) "." else opt$outdir
xmax <- if (is.null(opt$xmax)) 1000 else as.numeric(opt$xmax)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
if (!is.null(opt$theme) && file.exists(opt$theme)) source(opt$theme)
base_theme <- if (exists("theme_paper")) theme_paper() else theme_classic(base_size = 11)
pal_get <- function(k, d) if (exists("pal") && is.list(pal) && !is.null(pal[[k]])) pal[[k]] else d
cycle_cols <- function(keys) {
  base <- c(pal_get("purple", "#86187f"), pal_get("steel", "#6b93ad"), pal_get("slate", "#425b74"),
            pal_get("grey", "#a7b4c2"), pal_get("navy", "#1d3557"))
  setNames(rep_len(base, length(keys)), keys)
}
save_svg <- function(p, name, w, h) {
  grDevices::svg(filename = file.path(outdir, name), width = w, height = h)
  print(p)
  invisible(dev.off())
}

s <- read.delim(opt$summary, check.names = FALSE)
dist <- read.delim(opt$distribution, check.names = FALSE)
if (nrow(s) == 0 || nrow(dist) == 0) die("empty fragment size table(s)")
# read.delim turns an all-blank column into NA: blank or NA condition -> sample_id
cond_raw <- if ("condition" %in% names(s)) as.character(s$condition) else rep("", nrow(s))
cond_raw[is.na(cond_raw)] <- ""
s$cond <- ifelse(cond_raw == "", as.character(s$sample_id), cond_raw)
s <- s[order(s$cond, s$sample_id), ]
cond_levels <- unique(s$cond)
cond_cols <- cycle_cols(cond_levels)
s$cond <- factor(s$cond, levels = cond_levels)

# smoothed density per sample from exact counts (dist has only sizes with count > 0)
grid <- 1:xmax
dens <- do.call(rbind, lapply(s$sample_id, function(id) {
  d <- dist[dist$sample_id == id, ]
  if (nrow(d) == 0) return(NULL)
  tot <- sum(d$count)
  y <- vapply(grid, function(g) sum(d$count * dnorm(g, mean = d$size, sd = 4)) / tot, numeric(1))
  data.frame(sample_id = id, cond = s$cond[s$sample_id == id], size = grid, density = y)
}))
vl <- c(120, 200, 400, 600)
p1 <- ggplot(dens, aes(size, density, colour = cond, group = sample_id)) +
  geom_line(linewidth = 0.6, alpha = 0.8) +
  geom_vline(xintercept = vl, linetype = "dashed", colour = pal_get("rule", "#9aa1a9"), linewidth = 0.4) +
  annotate("text", x = vl + 8, y = Inf, label = paste0(vl, " bp"), angle = 90, vjust = 1.2, hjust = 1.1, size = 3,
           colour = pal_get("muted", "#5f6670")) +
  scale_colour_manual(values = cond_cols, name = NULL) +
  scale_x_continuous(limits = c(0, xmax), expand = expansion(mult = c(0, 0.01))) +
  labs(x = "Fragment size (bp)", y = "Density") + base_theme
save_svg(p1, "fragment_sizes_density.svg", 8, 4.5)

p2 <- ggplot(s, aes(cond, median_size, colour = cond, fill = cond)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, linewidth = 0.5) +
  geom_jitter(width = 0.15, height = 0, size = 2, alpha = 0.9) +
  scale_colour_manual(values = cond_cols, guide = "none") + scale_fill_manual(values = cond_cols, guide = "none") +
  labs(x = NULL, y = "Median fragment size (bp)") + base_theme +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_svg(p2, "fragment_sizes_boxplot.svg", max(4, 0.7 * length(cond_levels) + 2), 4.5)
cat("fragment size plots written\n")
