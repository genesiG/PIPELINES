# deg_common.R -- helpers shared by the deg_*.R scripts of the DEG module (sourced with --common=deg_common.R).
#
# Not a command line script: every deg_*.R script does `source(<--common value>)` right after start-up.
# Legacy mapping: parse/IO helpers are new; `my_theme` and `makeVolcanoPlot()` are ported from the legacy
#   Scripts/utils/makeVolcanoPlot.R of the GG_EZH2i_RNAseq project. The only change in makeVolcanoPlot() is that
#   `alpha` (FDR cut-off of the horizontal dashed line) is an ARGUMENT; the legacy function read a global `alpha`.
#
# Exit codes used by all deg_*.R scripts: 0 ok, 3 known data/usage problem (message on stderr, no retry),
#   1 = any other (unexpected) R error.
#
# Functions
#   parse_args(args)               --key=value -> named list (keys keep dashes, e.g. opt[["sample-col"]])
#   opt_get(opt, key, default)     value of key, default when absent
#   die(msg, status = 3)           "ERROR: msg" on stderr, quit(status)
#   as_bool(x, name)               TRUE/FALSE (also T/F/yes/no/1/0, any case)
#   parse_limits(x)                "a,b" -> numeric(2); "" / "auto" -> NULL
#   read_csv_clean(path, trim)     csv or tsv as all-character data.frame (BOM removed, blank rows dropped)
#   apply_label_rules(labels, rules_df)   ordered regex gsub, columns pattern,replacement ("\n" typed -> newline)
#   read_colors_csv(path)          csv group,color -> named character vector in file order (NULL when empty)
#   group_colors(levels, csv)      named colours for `levels` (csv first, Set2 / hue palette for the missing ones)
#   read_comparisons_tsv(path)     normalised comparisons table with every contract column (blank = "")
#   comparison_label(df)           label column, or comparison with "_" replaced by " " when blank
#   set_font(opt) / open_svg() / open_pdf()  module font (--font=Helvetica) and svg device (--svg-device=svglite|cairo)
#   draw_svg(plot, path, w, h)     set.seed(1) + cairo svg device (module font) + print
#   my_theme, makeVolcanoPlot()    volcano theme / plot (see below)

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))

# ---------------------------------------------------------------- command line / messages ----
die <- function(msg, status = 3) {
  cat("ERROR: ", msg, "\n", sep = "", file = stderr())
  quit(save = "no", status = status)
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  for (m in regmatches(args, regexec("^--([^=]+)=(.*)$", args))) {
    if (length(m) == 3) out[[m[2]]] <- m[3]
  }
  bad <- args[!grepl("^--[^=]+=", args)]
  if (length(bad)) message("WARNING: ignoring argument(s) not of the form --key=value: ", paste(bad, collapse = " "))
  out
}

opt_get <- function(opt, key, default = "") {
  v <- opt[[key]]
  if (is.null(v)) default else v
}

# a required argument; dies with a usage hint when absent or empty
opt_need <- function(opt, key) {
  v <- opt[[key]]
  if (is.null(v) || !nzchar(trimws(v))) die(paste0("missing required argument --", key))
  v
}

is_blank <- function(x) is.na(x) | !nzchar(trimws(as.character(x)))

as_bool <- function(x, name = "value") {
  if (is.logical(x)) return(x)
  v <- tolower(trimws(as.character(x)))
  if (v %in% c("true", "t", "yes", "y", "1")) return(TRUE)
  if (v %in% c("false", "f", "no", "n", "0")) return(FALSE)
  die(paste0("--", name, " must be TRUE or FALSE (got '", x, "')"))
}

as_num <- function(x, name = "value") {
  v <- suppressWarnings(as.numeric(x))
  if (length(v) != 1 || is.na(v)) die(paste0("--", name, " must be a number (got '", x, "')"))
  v
}

# "a,b" -> c(a, b); "" or "auto" -> NULL (= automatic limits)
parse_limits <- function(x, name = "limits") {
  if (is.null(x) || is.na(x)) return(NULL)
  x <- trimws(x)
  if (!nzchar(x) || tolower(x) == "auto") return(NULL)
  v <- suppressWarnings(as.numeric(trimws(strsplit(x, ",")[[1]])))
  if (length(v) != 2 || anyNA(v) || v[1] >= v[2]) die(paste0(name, " must be 'min,max' with min < max (got '", x, "')"))
  v
}

# --------------------------------------------------------------------------------- tables ----
# Reads a csv (or tsv: detected from a tab in the header line) as an ALL-CHARACTER data.frame.
# UTF-8 BOM removed, cells trimmed when trim = TRUE, fully blank rows and unnamed empty columns dropped.
# "NA" is NOT turned into NA (a gene/label called NA stays text); blank stays "".
read_csv_clean <- function(path, trim = TRUE) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) die(paste0("file not found: '", path, "'"))
  first <- readLines(path, n = 1, warn = FALSE, encoding = "UTF-8")
  first <- sub("^﻿", "", first)
  tsv <- length(first) == 1 && grepl("\t", first, fixed = TRUE)
  df <- tryCatch(
    read.table(path,
      header = TRUE, sep = if (tsv) "\t" else ",", quote = if (tsv) "" else "\"",
      colClasses = "character", na.strings = character(0), check.names = FALSE,
      strip.white = FALSE, comment.char = "", blank.lines.skip = TRUE, fill = TRUE,
      fileEncoding = "UTF-8-BOM", stringsAsFactors = FALSE
    ),
    error = function(e) {
      if (grepl("no lines available", conditionMessage(e))) {
        data.frame()
      } else {
        die(paste0("cannot read '", path, "': ", conditionMessage(e)))
      }
    }
  )
  names(df) <- trimws(names(df))
  if (trim) df[] <- lapply(df, trimws)
  keep_col <- nzchar(names(df)) | vapply(df, function(v) any(nzchar(trimws(v))), logical(1))
  df <- df[, keep_col, drop = FALSE]
  if (nrow(df) > 0) {
    blank_row <- rowSums(as.data.frame(lapply(df, function(v) nzchar(trimws(v))))) == 0
    df <- df[!blank_row, , drop = FALSE]
  }
  rownames(df) <- NULL
  df
}

# Ordered regex substitution of labels. rules_df: columns `pattern`, `replacement` (file order = application
# order, like the legacy chain of gsub() calls; R extended regex, not perl, not fixed). A replacement typed as
# the two characters \n becomes a real newline.
apply_label_rules <- function(labels, rules_df) {
  if (is.null(rules_df) || nrow(rules_df) == 0) return(labels)
  if (!all(c("pattern", "replacement") %in% names(rules_df))) {
    die("label rules table needs the columns 'pattern' and 'replacement'")
  }
  out <- as.character(labels)
  for (i in seq_len(nrow(rules_df))) {
    rep_txt <- gsub("\\n", "\n", rules_df$replacement[i], fixed = TRUE)
    out <- tryCatch(gsub(rules_df$pattern[i], rep_txt, out),
      error = function(e) die(paste0("bad label rule regex '", rules_df$pattern[i], "': ", conditionMessage(e)))
    )
  }
  out
}

# csv with columns group,color -> named character vector in file order; NULL when path empty
read_colors_csv <- function(path) {
  if (is.null(path) || !nzchar(trimws(path))) return(NULL)
  df <- read_csv_clean(path)
  if (!all(c("group", "color") %in% names(df))) die("colours table needs the columns 'group' and 'color'")
  df <- df[nzchar(df$group) & nzchar(df$color), , drop = FALSE]
  if (nrow(df) == 0) return(NULL)
  setNames(df$color, df$group)
}

default_palette <- function(n) {
  if (n <= 0) return(character(0))
  if (n <= 8) RColorBrewer::brewer.pal(max(3, n), "Set2")[seq_len(n)] else scales::hue_pal()(n)
}

# colours for `levels` (named by level). The csv wins; groups missing from it get the default palette colour
# of their position (Set2 up to 8 groups, scales::hue_pal beyond).
group_colors <- function(levels, colors_csv = "") {
  levels <- as.character(levels)
  cols <- setNames(default_palette(length(levels)), levels)
  given <- read_colors_csv(colors_csv)
  if (!is.null(given)) {
    hit <- intersect(levels, names(given))
    cols[hit] <- given[hit]
    miss <- setdiff(levels, names(given))
    if (length(miss)) message("NOTE: no colour given for group(s) ", paste(miss, collapse = ", "), "; default palette used")
  }
  cols
}

# Comparisons table (normalised comparisons.tsv written by deg_prepare.R, or the raw user csv): all contract
# columns present, blank = "".
comparison_columns <- c(
  "model", "design", "comparison", "factor", "numerator", "denominator", "coef", "add", "subtract",
  "label", "combined_order", "combined_title", "xlim", "ylim"
)
read_comparisons_tsv <- function(path) {
  df <- read_csv_clean(path)
  need <- c("model", "design", "comparison")
  if (!all(need %in% names(df))) die(paste0("comparisons table needs the columns ", paste(need, collapse = ", ")))
  for (cn in setdiff(comparison_columns, names(df))) df[[cn]] <- ""
  df
}
comparison_label <- function(df) {
  ifelse(is_blank(df$label), gsub("_", " ", df$comparison), df$label)
}

# ---- font: every plot of the module is drawn with one family (--font=Helvetica, the default) ----
# The family is set on the DEVICE (svg / pdf), so ggplot2, ggrepel, pheatmap and base graphics all use it; the theme
# keeps family "" (= the device family). svglite (default) writes real <text> elements (see open_svg); the cairo svg
# device (--svg-device=cairo) converts text to glyph outlines. Text is measured through fontconfig: where Helvetica is
# not installed or aliased (Nimbus Sans / Liberation Sans) another sans font is used for the layout.
set_font <- function(opt) {
  options(deg.font = opt_get(opt, "font", "Helvetica"))
  dev <- tolower(opt_get(opt, "svg-device", "svglite"))
  if (!dev %in% c("svglite", "cairo")) die("--svg-device must be svglite or cairo")
  options(deg.svg.device = dev)
  invisible(getOption("deg.font"))
}
deg_font <- function() getOption("deg.font", "Helvetica")

# svg device. svglite (default) writes real <text> elements with font-family = the module font (editable in
# Illustrator / Inkscape); --svg-device=cairo gives the legacy grDevices::svg() device (text converted to glyph outlines).
# svglite measures text through systemfonts / fontconfig, so layout follows the font that fontconfig resolves.
open_svg <- function(path, width, height) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (getOption("deg.svg.device", "svglite") == "cairo") {
    grDevices::svg(filename = path, width = width, height = height, family = deg_font())
  } else {
    f <- deg_font()
    svglite::svglite(filename = path, width = width, height = height,
      system_fonts = list(sans = f, serif = f, mono = f, symbol = f))
  }
  options(deg.svg.path = path)
  invisible(path)
}

# Closes the device opened by open_svg(). svglite writes the font family that fontconfig resolved (e.g. "Nimbus Sans"
# when Helvetica is an alias of it); the family names of the svg text are rewritten to the requested font (default
# exactly `Helvetica`) so that Illustrator / Inkscape / a browser use the user's own Helvetica. Widths are kept: every
# text element carries textLength.
close_svg <- function() {
  grDevices::dev.off()
  path <- getOption("deg.svg.path")
  if (getOption("deg.svg.device", "svglite") == "svglite" && !is.null(path) && file.exists(path)) {
    txt <- readLines(path, warn = FALSE, encoding = "UTF-8")
    txt <- gsub("font-family: ('[^']*'|\"[^\"]*\"|[^;]*);", paste0("font-family: ", deg_font(), ";"), txt)
    writeLines(txt, path, useBytes = TRUE)
  }
  options(deg.svg.path = NULL)
  invisible(NULL)
}

# pdf device with the module font. Helvetica / Arial / sans use the pdf() built-in Helvetica (legacy: pdf(file, width,
# height), whose default family is Helvetica); any other family name goes through cairo_pdf (fontconfig).
open_pdf <- function(path, width, height) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  if (tolower(deg_font()) %in% c("helvetica", "arial", "sans")) {
    grDevices::pdf(file = path, width = width, height = height, family = "Helvetica")
  } else {
    grDevices::cairo_pdf(filename = path, width = width, height = height, family = deg_font())
  }
}

# set.seed(1) before drawing (ggrepel label placement is random) + cairo svg device
draw_svg <- function(plot, path, width, height) {
  open_svg(path, width, height)
  set.seed(1)
  on.exit(close_svg())
  print(plot)
  invisible(path)
}

# ------------------------------------------------------------------------ volcano plot ----
# Ported verbatim from the legacy utils/makeVolcanoPlot.R (theme and function), except `alpha` is an argument.
my_theme <- theme_bw(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2),
    strip.background = element_blank(),
    strip.text = element_text(size = 17),
    plot.title = element_text(size = 17, hjust = 0.5),
    legend.position = "none",
    axis.text = element_text(color = "black", size = 14),
    axis.title = element_text(color = "black", size = 15),
    axis.line = element_blank()
  )

# Takes the results data frame (Gene, log2FoldChange, padj, effect, shrunkLFC, ..., Label) and returns a ggplot.
makeVolcanoPlot <- function(results_df,
                            lfc.threshold = 1,
                            alpha = 0.05,
                            ylim = NULL,
                            xlim = NULL,
                            color_by = "effect",
                            label_by = "Label",
                            highlight_by = NULL,
                            title = "",
                            y_lab = (expression(bold(bolditalic(-Log[10]) ~ "FDR"))),
                            x_lab = (expression(bold(bolditalic(Log[2]) ~ "Fold-Change"))),
                            color_labels = c("Downregulated", "NS", "Upregulated"),
                            color_values = c("navyblue", "gray", "red"),
                            x_axis = "log2FoldChange",
                            y_axis = "padj",
                            nbreaks = 5,
                            nudgex = 0,
                            nudgey = 0,
                            force = 2,
                            force_pull = 4,
                            label_size = 3,
                            font_color = NULL,
                            leg.pos = "bottom",
                            maxoverlaps = 10,
                            point_size = 2,
                            highlight_stroke = 0.5,
                            label_fontface = "italic",
                            keep_labels_in_bounds = TRUE) {
  # Ensure that the number of labels matches the number of colors provided
  if (length(color_labels) != length(color_values)) {
    stop("The number of color labels must match the number of color values.")
  }

  # Convert `color_labels` and `color_values` into a named vector for scale_color_manual
  color_map <- setNames(color_values, color_labels)

  # Check if all labels are NA - if so, don't use geom_text_repel
  all_labels_na <- all(is.na(results_df[[label_by]]))

  volcanoPlot <- ggplot(results_df, aes(
    x = .data[[x_axis]], # use .data pronoun with strings
    y = -log10(.data[[y_axis]]),
    color = .data[[color_by]]
  ))

  volcanoPlot <- volcanoPlot + geom_point(shape = 16, size = point_size)

  # Add border to highlighted points
  if (!is.null(highlight_by)) {
    highlight_df <- results_df[!is.na(results_df[[highlight_by]]), ]
    if (nrow(highlight_df) > 0) {
      volcanoPlot <- volcanoPlot +
        geom_point(
          data = highlight_df,
          shape = 21,
          size = point_size,
          stroke = highlight_stroke,
          color = "black", fill = NA
        )
    }
  }

  volcanoPlot <- volcanoPlot +
    my_theme +
    # Add axes labels
    ylab(y_lab) +
    xlab(x_lab) +
    # Add plot tile based on which comparison was made
    ggtitle(title) +
    # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold
    geom_vline(xintercept = -lfc.threshold, lty = "dashed", alpha = 0.7) +
    geom_vline(xintercept = lfc.threshold, lty = "dashed", alpha = 0.7) +
    geom_hline(yintercept = -log10(alpha), lty = "dashed", alpha = 0.7) +
    scale_x_continuous(
      n.breaks = nbreaks,
      limits = xlim
    ) +
    scale_y_continuous(limits = ylim, n.breaks = (nbreaks - 2)) +
    # Modify plotting theme
    theme(
      legend.title = element_blank(),
      legend.position = leg.pos,
    ) +
    # Apply custom color scale based on user-defined labels and colors
    scale_color_manual(values = color_map) +
    guides(color = guide_legend(override.aes = list(shape = 16, size = 3)))

  # Only add the label layer if there are non-NA labels
  if (!all_labels_na) {
    # Only the labelled rows go to the label layer (ggrepel dropped the NA labels anyway).
    label_df <- results_df[!is.na(results_df[[label_by]]), , drop = FALSE]
    # nudge_x / nudge_y shift the label start from its point. ggplot2 applies the scale limits (oob = censor) AFTER
    # the nudge, so a point near the top / right edge whose nudged start is outside the limits lost its label (legacy
    # bug: CDKN2A at -log10(padj) = 25.6, ylim 0-26, nudge 1.5 -> start 27.1 -> label dropped, only the ring drawn).
    # With keep_labels_in_bounds the nudge is flipped to the other side, or shortened, when the nudged start would
    # leave the axis limits; labels that fit keep exactly the legacy nudge.
    keep_in <- function(v, nudge, lim) {
      out <- rep(nudge, length(v))
      if (is.null(lim) || !keep_labels_in_bounds) return(out)
      bad <- is.finite(v) & (v + nudge < lim[1] | v + nudge > lim[2])
      flip_ok <- v - nudge >= lim[1] & v - nudge <= lim[2]
      out[bad & flip_ok] <- -nudge
      out[bad & !flip_ok] <- pmin(pmax(v[bad & !flip_ok] + nudge, lim[1]), lim[2]) - v[bad & !flip_ok]
      out
    }
    nudge_x_v <- keep_in(label_df[[x_axis]], nudgex, xlim)
    nudge_y_v <- keep_in(-log10(label_df[[y_axis]]), nudgey, ylim)
    volcanoPlot <- volcanoPlot +
      ggrepel::geom_label_repel(
        data = label_df,
        aes(label = .data[[label_by]]),
        color = font_color,
        nudge_x = nudge_x_v,
        nudge_y = nudge_y_v,
        max.overlaps = maxoverlaps,
        max.time = 2,
        force = force,
        force_pull = force_pull,
        size = label_size,
        show.legend = FALSE,
        fontface = label_fontface,
        box.padding = 0.4,
        point.padding = 0.3,
        segment.color = "black", # font_color,
        segment.size = 0.35,
        min.segment.length = 0.1
      )
  } else {
    cat("\nSkipping labels\n")
  }

  return(volcanoPlot)
}
