#!/usr/bin/env Rscript
#' ##########################################################################
#' guide_representation_report.R
#'
#' Assembles one self-contained HTML QC report per library from the tables
#' guide_representation_wrangle.R wrote and the SVGs guide_representation_plot.R
#' rendered from them. For every plot: what it measures, how to read it, and
#' what to check next if it looks poor. Every section opens with the plain-language
#' question it answers (the plot name is the subtitle under it). Sections: guide
#' representation, ECDF, Lorenz, density, controls, 90/10, undetected, mapping, read
#' classification, GC, sample correlation heatmaps, PCA and read depth.
#' No statistics are computed here — this only reads existing tables/plots and writes HTML.
#'
#' USAGE:
#'   Rscript guide_representation_report.R --tables-dir=<dir> --plots-dir=<dir> \
#'       --prefix=<prefix> [--project=<name>] [--outdir=<dir>] [--breakdown=<file>]
#'
#'   --tables-dir  directory holding <prefix>_*.csv from the wrangle step
#'   --plots-dir   directory holding <prefix>_*.svg from the plot step
#'   --prefix      count-table prefix (e.g. MYSCREEN_library)
#'   --project     project_name, shown in the report header (optional)
#'   --outdir      where to write the report (default: current directory)
#'   --breakdown   <project>_library_breakdown_wide.csv from bin/library_breakdown.py
#'                 (optional; used in all-vs-all runs to show the reads mapped to
#'                 each library. Without it the report is unchanged)
#'
#' OUTPUT: <prefix>_qc_report.html (a single file — every plot is embedded
#' inline as SVG, so it can be emailed or opened from any machine with no
#' other files alongside it).
#' ##########################################################################

options(stringsAsFactors = FALSE)

# ---- Arguments (--key=value) ----------------------------------------------
parse_args <- function(args) {
  out <- list()
  for (m in regmatches(args, regexec("^--([^=]+)=(.*)$", args))) {
    if (length(m) == 3) out[[m[2]]] <- m[3]
  }
  out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt[["tables-dir"]]) || is.null(opt[["plots-dir"]]) || is.null(opt$prefix)) {
  stop("Usage: guide_representation_report.R --tables-dir=<dir> --plots-dir=<dir> --prefix=<prefix> [--project=<name>] [--outdir=<dir>] [--breakdown=<file>]")
}
tables_dir <- opt[["tables-dir"]]
plots_dir  <- opt[["plots-dir"]]
prefix     <- opt$prefix
project    <- if (is.null(opt$project)) prefix else opt$project
outdir     <- if (is.null(opt$outdir)) "." else opt$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- Helpers ----------------------------------------------------------------
html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x
}

read_table_opt <- function(name) {
  path <- file.path(tables_dir, paste0(prefix, "_", name, ".csv"))
  if (!file.exists(path)) return(NULL)
  read.csv(path, check.names = FALSE)
}

# R's svg() device (Cairo) numbers its internal <symbol>/<clipPath> ids from
# scratch in every file it writes (glyph0-0, clip1, ...). Embedding several
# such SVGs in one HTML page makes those ids collide: the browser resolves
# each <use xlink:href="#glyph0-0"> to whichever plot's <symbol id="glyph0-0">
# happens to win, so letters render with another plot's glyph shapes ("H"
# from one title, "s" from another, ...). Namespacing every id (and every
# reference to it) with something unique to this plot fixes it. `ns` must be
# unique among everything embedded in one page — the plot's own name suffix
# (e.g. "ecdf", "pearson_correlation") is enough within one report.
namespace_svg_ids <- function(svg_text, ns) {
  ids <- unique(regmatches(svg_text, gregexpr('(?<=id=")[^"]+', svg_text, perl = TRUE))[[1]])
  for (id in ids) {
    new_id <- paste0(ns, "-", id)
    svg_text <- gsub(paste0('id="', id, '"'),  paste0('id="', new_id, '"'),  svg_text, fixed = TRUE)
    # Matches both href="#id" and xlink:href="#id" (substring match).
    svg_text <- gsub(paste0('href="#', id, '"'), paste0('href="#', new_id, '"'), svg_text, fixed = TRUE)
    # fill="url(#id)", clip-path="url(#id)", mask="url(#id)", ...
    svg_text <- gsub(paste0('url(#', id, ')'), paste0('url(#', new_id, ')'), svg_text, fixed = TRUE)
  }
  svg_text
}

# Reads a plot script SVG file and returns just its <svg ...>...</svg> element
# (stripping the XML declaration / DOCTYPE) with its internal ids namespaced,
# so several can be embedded inline in one HTML document. NULL if missing.
read_svg <- function(name) {
  path <- file.path(plots_dir, paste0(prefix, "_", name, ".svg"))
  if (!file.exists(path)) return(NULL)
  raw <- paste(readLines(path, warn = FALSE), collapse = "\n")
  start <- regexpr("<svg[ >]", raw)
  if (start < 0) return(NULL)
  end <- regexpr("</svg>", raw)
  if (end < 0) return(NULL)
  svg <- substr(raw, start, end + attr(end, "match.length") - 1)
  namespace_svg_ids(svg, name)
}

fmt_int <- function(x) ifelse(is.na(x), "—", format(round(x), big.mark = ",", scientific = FALSE))
fmt_pct <- function(x, digits = 1) ifelse(is.na(x), "—", paste0(round(x, digits), "%"))
fmt_num <- function(x, digits = 3) ifelse(is.na(x), "—", format(round(x, digits), nsmall = digits))

badge <- function(text, kind) sprintf('<span class="badge badge-%s">%s</span>', kind, html_escape(text))

# ---- Load the tables every report needs ------------------------------------
counts_long <- read_table_opt("counts_long")
if (is.null(counts_long)) {
  stop("Missing ", prefix, "_counts_long.csv in ", tables_dir,
       " — guide_representation_wrangle.R must run (successfully) before this report.")
}
qc <- read_table_opt("qc_summary")
if (is.null(qc)) {
  stop("Missing ", prefix, "_qc_summary.csv in ", tables_dir,
       " — guide_representation_wrangle.R must run (successfully) before this report.")
}

has_controls <- any(as.logical(counts_long$is_control))
n_samples    <- length(unique(counts_long$sample))
one_sample   <- n_samples == 1
has_mapping  <- "mapping_pct" %in% names(qc)
low_mapping  <- if (has_mapping) qc$mapping_pct < 30 else rep(FALSE, nrow(qc))
has_class    <- "pct_near_miss" %in% names(qc)
has_rpg      <- "reads_per_guide" %in% names(qc)

# ---- Optional library breakdown (bin/library_breakdown.py, wide file) ------
# Used in all-vs-all runs (several libraries, no samplesheet library column),
# where every sample was counted against every library. Without a usable file
# none of the breakdown logic below applies.
bd_wide <- NULL
if (!is.null(opt[["breakdown"]]) && file.exists(opt[["breakdown"]]) &&
    file.size(opt[["breakdown"]]) > 0) {
  bd_wide <- tryCatch(read.csv(opt[["breakdown"]], check.names = FALSE),
                      error = function(e) NULL)
}
bd_pct_cols <- character(0)
if (!is.null(bd_wide) && "sample" %in% names(bd_wide)) {
  bd_pct_cols <- setdiff(grep("^pct_", names(bd_wide), value = TRUE),
                         c("pct_any_library", "pct_no_library_guide"))
}
has_breakdown <- length(bd_pct_cols) > 0
# This report's library name: the prefix without its leading "<project>_".
project_stub <- paste0(project, "_")
lib_name <- if (startsWith(prefix, project_stub)) substring(prefix, nchar(project_stub) + 1) else prefix
this_col <- paste0("pct_", lib_name)
# Row of the breakdown for each sample of this report (NA if absent).
bd_idx <- if (has_breakdown) match(qc$sample, bd_wide$sample) else rep(NA_integer_, nrow(qc))
bd_in  <- !is.na(bd_idx)
bd_pick <- function(col, numeric = TRUE) {
  x <- if (has_breakdown && col %in% names(bd_wide)) bd_wide[[col]][bd_idx] else rep(NA, nrow(qc))
  if (numeric) suppressWarnings(as.numeric(x)) else as.character(x)
}
bd_any_pct <- bd_pick("pct_any_library")
bd_main    <- bd_pick("main_library", numeric = FALSE)

# ---- "At a glance" summary table (every headline metric, one row/sample) ---
flags_col <- vapply(seq_len(nrow(qc)), function(i) {
  b <- character(0)
  if (isTRUE(qc$library_skew_flag[i]))          b <- c(b, badge("Library skew",  "flag"))
  if (isTRUE(qc$undetected_flag[i]))            b <- c(b, badge("High dropout",  "flag"))
  if (bd_in[i]) {
    # All-vs-all run: "Low mapping" means the sample matches no library well.
    if (isTRUE(bd_any_pct[i] < 30)) b <- c(b, badge("Low mapping", "flag"))
    # A sample mainly matching another library gets a neutral note only,
    # and only when it is not already flagged Low mapping.
    if (!isTRUE(bd_any_pct[i] < 30) && !is.na(bd_main[i]) && nzchar(bd_main[i]) &&
        bd_main[i] != "none" && bd_main[i] != lib_name) {
      b <- c(b, badge(paste0("Mainly ", bd_main[i]), "info"))
    }
  } else if (has_mapping && isTRUE(low_mapping[i])) {
    b <- c(b, badge("Low mapping",   "flag"))
  }
  if (length(b) == 0) b <- badge("OK", "ok")
  paste(b, collapse = " ")
}, character(1))

summary_rows <- paste0(
  "<tr>",
  "<td>", html_escape(qc$sample), "</td>",
  if (has_mapping) paste0("<td class='num'>", fmt_int(qc$reads), "</td><td class='num'>", fmt_int(qc$mappedreads), "</td><td class='num'>", fmt_pct(qc$mapping_pct), "</td>") else "",
  if (has_rpg) paste0("<td class='num'>", fmt_num(qc$reads_per_guide, 1), "</td>") else "",
  if (has_class) paste0("<td class='num'>", fmt_pct(qc$pct_near_miss), "</td><td class='num'>", fmt_pct(qc$pct_far_miss), "</td>") else "",
  "<td class='num'>", fmt_num(qc$Gini, 4), "</td>",
  "<td class='num'>", fmt_num(qc$ratio_90_10, 2), "</td>",
  "<td class='num'>", fmt_pct(qc$pct_undetected, 2), " (", fmt_int(qc$zero_guides), ")</td>",
  "<td class='num'>", fmt_pct(qc$pct_within_2_log2_mean, 1), "</td>",
  "<td>", flags_col, "</td>",
  "</tr>",
  collapse = "\n"
)

summary_table <- paste0(
  "<table class='qc-table'><thead><tr>",
  "<th>Sample</th>",
  if (has_mapping) "<th>Total Reads</th><th>Mapped Reads</th><th>Mapping % (Exact)</th>" else "",
  if (has_rpg) "<th>Reads per Guide</th>" else "",
  if (has_class) "<th>% Near-Miss</th><th>% Far-Miss</th>" else "",
  "<th>Gini</th><th>90/10 Ratio</th><th>% Undetected (n)</th><th>% Within 4-fold of Mean</th><th>Flags</th>",
  "</tr></thead><tbody>",
  summary_rows,
  "</tbody></table>"
)

# ---- "Reads mapped to each library" table (only with a breakdown) ----------
# Shown under the Mapping Rate plot. One row per sample of this report.
mapping_table_html <- ""
if (has_breakdown && any(bd_in)) {
  lib_heads  <- sub("^pct_", "", bd_pct_cols)
  pct_vals   <- lapply(bd_pct_cols, bd_pick)
  bd_reads   <- bd_pick("reads")
  bd_none    <- bd_pick("pct_no_library_guide")
  head_cells <- paste0("<th", ifelse(bd_pct_cols == this_col, " class='this-lib'", ""),
                       ">", html_escape(lib_heads), "</th>", collapse = "")
  rows_bd <- vapply(which(bd_in), function(i) {
    cells <- vapply(seq_along(bd_pct_cols), function(k) {
      cls <- c("num",
               if (bd_pct_cols[k] == this_col) "this-lib",
               if (!is.na(bd_main[i]) && bd_main[i] == lib_heads[k]) "main-lib")
      paste0("<td class='", paste(cls, collapse = " "), "'>", fmt_pct(pct_vals[[k]][i], 1), "</td>")
    }, character(1))
    main_txt <- if (is.na(bd_main[i]) || !nzchar(bd_main[i])) "&mdash;" else html_escape(bd_main[i])
    paste0("<tr><td>", html_escape(qc$sample[i]), "</td><td class='num'>", fmt_int(bd_reads[i]), "</td>",
           paste(cells, collapse = ""),
           "<td class='num'>", fmt_pct(bd_any_pct[i], 1), "</td>",
           "<td class='num'>", fmt_pct(bd_none[i], 1), "</td>",
           "<td>", main_txt, "</td></tr>")
  }, character(1))
  mapping_table_html <- paste0(
    "<h3>Reads mapped to each library</h3>",
    "<div class='table-scroll'><table class='qc-table'><thead><tr><th>Sample</th><th>Total Reads</th>",
    head_cells,
    "<th>Any library</th><th>No library guide</th><th>Main library</th></tr></thead><tbody>",
    paste(rows_bd, collapse = "\n"), "</tbody></table></div>"
  )
}

# ---- Plot sections: description / interpretation / troubleshooting --------
# Plain, short sentences — this report is read by wet-lab scientists, not
# statisticians. Each entry: id (report anchor), svgs (one or more
# <prefix>_<svg>.svg names), title, description (what it is), interpretation
# (what good/bad looks like, 2-4 short sentences), troubleshooting (what to
# check next, as a short bullet list), and an optional missing_reason shown
# in plain English instead of the plot when it wasn't generated.
# `avg_word`/`ecdf_*`/`density_*` below adapt wording for a single-sample
# library, where "average across samples" and "compare curves" don't apply.
avg_word <- if (one_sample) "counts per million reads (CPM)" else "average counts per million reads (CPM) across this library's samples"
density_desc <- if (one_sample) {
  "The shape of this sample's read-count distribution: where most guides sit, and whether there's a separate cluster elsewhere in the range."
} else {
  "The shape of each sample's read-count distribution, overlaid: where most guides sit, and whether there's a separate cluster elsewhere in the range."
}
# Mapping Rate wording: the breakdown paragraphs are added only when present.
mapping_desc_base <- r"(The percent of sequenced reads that matched a library guide, from the mageck
    count log. mageck count's own matching is always exact &mdash; it has no mismatch-tolerance
    option &mdash; so this is also, exactly, the percentage of reads that exactly matched a guide.
    The Read Classification plot below breaks the rest (the non-matching reads) down further.)"
mapping_desc <- if (has_breakdown) {
  paste0(mapping_desc_base, r"( In this run every sample was counted against every library in the
    pool, so the Mapping Rate shown here is the share of each sample's reads that match THIS library.
    A low value just means the sample belongs to another library. The table below shows, for each
    sample, the percent of its reads that matched each library. "Low mapping" is flagged only when a
    sample matches no library well (Any library under 30%).)")
} else mapping_desc_base
mapping_interp_base <- r"(A healthy pooled screen typically maps 60&ndash;70% or more of its reads;
    this report flags anything under 30%. Mapping rate is independent of the representation metrics
    above: a sample can map poorly yet still show even representation among the reads that did map
    — that just means most of its sequencing depth was wasted, not misrepresented.)"
mapping_interp <- if (has_breakdown) {
  paste0(mapping_interp_base, r"( Red bars in the plot above (under 30%) are expected for samples of
    other libraries; read the plot together with the table below.)")
} else mapping_interp_base

# PCA needs at least 3 samples; the reason shown when it is missing depends on that.
pca_missing_reason <- if (n_samples < 3) {
  paste0("PCA needs at least 3 samples with reads in this library; this library has ", n_samples, ".")
} else "the PCA table was not found (the tables come from the wrangle step)."
# Guides in this library (for the Read Depth text), if the summary table has it.
n_guides <- if ("total_guides" %in% names(qc)) qc$total_guides[1] else NA
rd_guides_txt <- if (is.na(n_guides)) "" else paste0(" This library has ", fmt_int(n_guides), " guides.")

pca_desc <- r"(Each sample is one point. The plot is computed from log2(counts per million + 1) of
    this library's guides, using the samples that have reads in this library. Guides with no reads in
    any sample are removed. The data are centred but not scaled (PCA via <code>prcomp</code>). The axis
    labels give the percent of the total variance each principal component (PC) explains. Only PC1
    vs PC2 is plotted; PC3 is in the <code>pca.csv</code> table.)"
pca_interp_base <- r"(Points close together have similar guide profiles. Replicates should cluster
    together. Separation along PC1 between conditions or timepoints is expected after selection: it
    is the screen working. A sample far from its replicates is suspect (low depth, mix-up,
    contamination). The percentages say how much of the variation the plot shows.)"
pca_interp <- if (has_breakdown) {
  paste0(pca_interp_base, r"( Caveat for this all-vs-all run: every sample is counted against every
    library, so samples that belong to OTHER libraries have very few reads in this library and
    dominate the PCA. They separate from the rest simply because they hold other libraries. Read the
    PCA together with the Mapping Rate table above.)")
} else pca_interp_base

corr_desc <- r"(Samples are ordered by hierarchical clustering (average linkage on 1 minus the
    Pearson correlation of log10(count + 1)), so similar samples sit next to each other. Each cell
    shows the correlation of one pair of samples, to 2 decimals. The colour scale (0 to 1) is the
    same in every library, so libraries can be compared by eye. <strong>Pearson</strong> uses log
    counts and is sensitive to overall magnitude. <strong>Spearman</strong> uses ranks and is robust
    to outliers.)"
corr_interp_base <- r"(Replicates should form blocks of high values, typically above 0.9 for a
    technically sound screen. Lower correlation between conditions or timepoints is expected after
    selection. Low correlation between nominal replicates suggests a sample mix-up, a failed or
    low-depth replicate, or heavy dropout in one of them.)"
corr_interp <- if (has_breakdown) {
  paste0(corr_interp_base, r"( In all-vs-all runs, samples that belong to other libraries have very
    few reads here and will correlate poorly.)")
} else corr_interp_base

sections <- list(
  list(
    id = "guide_representation", svgs = "guide_representation",
    question = "1. Are the library's guides represented evenly in each sample?",
    title = "Guide Representation (Ranked)",
    description = paste0("Every guide's ", avg_word, ", sorted from lowest to highest, on a log scale (log10 of CPM + 1). CPM scales each sample to its own total reads, so samples of different sequencing depth sit on the same scale. This is the standard check for how evenly a library is represented."),
    interpretation = r"(The curve's slope reflects representation evenness. A shallow, gradual rise
    across the whole curve is the ideal case: most guides sit within a similar count range. A steep
    drop at the low-rank end indicates a subset of guides are barely present or missing entirely. A
    steep rise concentrated at the high-rank end — a small number of guides far above the rest of
    the curve — indicates PCR over-amplification of specific sequences. A distinct step or plateau
    in the middle of the curve, rather than a smooth rise, indicates two sub-populations of guides
    at different typical abundance — the same signal the Density plot below shows as two peaks.)",
    troubleshooting = r"(<li>Cross-check "% Undetected Guides" below to see whether the low-rank
    guides are truly at zero, or just underrepresented.</li>
    <li>A drop already present in a plasmid/T0 sample indicates a library or transduction coverage
    problem — aim for 300&ndash;500&times; coverage per guide.</li>
    <li>A drop appearing only after treatment or selection indicates real biology: guides targeting
    essential or fitness genes dropping out.</li>
    <li>A mid-curve step: see the Density and Control vs Targeting plots below for whether it lines
    up with control vs. targeting guides (often by design) or something else.</li>)"
  ),
  list(
    id = "ecdf", svgs = "ecdf",
    question = "2. What fraction of guides received few reads, and how deep is each sample on average?",
    title = "ECDF of gRNA Representation",
    description = r"(The top panel shows the cumulative fraction of guides at or below a given read
    count, one curve per sample. The bottom panel uses the same read-count axis and has one bar per
    sample, in the same colour as its curve. Each bar shows that sample's mean read count per guide
    &mdash; what every guide would have if reads were split perfectly evenly. The value is printed
    next to the bar.)",
    interpretation = r"(A curve that rises steeply close to the end of its bar (the mean) means even
    coverage. The more of the curve sitting far to the left of the mean, the more guides fall below
    their fair share of reads. The mean also tells you how deep each sample is per guide. Compare
    the bars across samples to spot a sample with far fewer reads per guide than its replicates.)",
    troubleshooting = r"(<li>Most of the curve sitting well to the left of its bar's end indicates a
    large fraction of guides are underrepresented. Cross-check "% Undetected Guides" below for how
    many sit at exactly zero.</li>
    <li>A bar much shorter than the bars of its replicates means that sample has far fewer reads per
    guide. Check the Read Depth section below.</li>
    <li>An unusual step or kink in the curve's shape indicates the sample may be a mixture of two
    populations, e.g. cross-contamination between pools.</li>)"
  ),
  list(
    id = "lorenz", svgs = "lorenz",
    question = "3. How unequally are reads shared among guides?",
    title = "Lorenz Curves (Gini)",
    description = r"(Guides sorted from least to most reads, plotted as cumulative % of reads
    against cumulative % of guides.)",
    interpretation = r"(The diagonal line represents an ideal sample with perfectly even gRNA
    representation: every guide contributing an identical share of reads. The further the curve
    bows away from the diagonal, the less evenly the library is represented. This is the same
    thing the Gini coefficient in the table above measures as one number &mdash; Gini is defined as
    twice the area between the curve and the diagonal, computed over all of this library's guides,
    including guides with zero reads.)",
    troubleshooting = r"(<li>A curve already bowed at T0/plasmid indicates a library or coverage
    problem — see Guide Representation above.</li>
    <li>Bowing that increases from T0 through treated/selected timepoints reflects selection
    pressure concentrating reads on surviving guides — expected, not a technical artifact.</li>
    <li>Confirm the Gini and 90/10 ratio (below) tell the same story.</li>)"
  ),
  list(
    id = "density", svgs = "density",
    question = "4. Do guides form one population or several (e.g. a dropout group near zero)?",
    title = "Read-Count Density",
    description = density_desc,
    interpretation = r"(A single, smooth peak is the ideal case: one uniform population of
    well-represented guides. A second peak near zero read counts indicates gRNA dropout: guides
    with almost no reads. The height of that second peak reflects how severe the dropout is — a
    small bump means a handful of guides dropped out; a tall peak means a large fraction of the
    library did. A second peak away from zero, at a substantial read count, is a different signal:
    two distinct sub-populations of guides with different typical abundance, rather than one
    population plus dropout. More than two peaks indicates more than two such sub-populations.)",
    troubleshooting = r"(<li>A second peak near zero, already present in a plasmid/T0 sample,
    indicates a subset of guides were cloned or synthesized less efficiently — check GC Content vs
    Read Depth below for a possible cause.</li>
    <li>A second peak near zero appearing only after selection indicates real dropout. The wrangle
    step's <code>zero_count_guides.csv</code> lists exactly which guides, useful for checking
    whether they cluster by gene.</li>
    <li>A second peak away from zero: check the Control vs Targeting plot below. If it splits
    cleanly into controls forming one peak and targeting guides the other, this is often
    intentional — some libraries include control guides at a different molar ratio than targeting
    guides for statistical power — not a problem. If the split doesn't line up with controls, it
    more likely indicates uneven pooling of sub-libraries or synthesis batches, or two samples'
    reads merged by mistake — check <code>samplesheet.csv</code> for a possible mix-up.</li>)"
  ),
  list(
    id = "controls_ecdf", svgs = "controls_ecdf",
    question = "5. Do control guides behave differently from targeting guides?",
    title = "Control vs Targeting gRNA Representation",
    description = r"(The same per-sample count distribution as the ECDF above, split into control
    guides (Control/NTC/Non-targeting) and targeting guides.)",
    interpretation = r"(Control guides aren't expected to affect cell fitness, so their curve is a
    built-in reference for what unaffected representation looks like in this exact sample. The
    ideal case is the control and targeting curves overlapping closely at early/reference
    timepoints, since neither has been selected on yet. The targeting curve pulling away from the
    control curve after treatment or selection indicates real selection: a subset of targeting
    guides genuinely dropping out or growing out.)",
    troubleshooting = r"(<li>Divergence already present at T0/plasmid, before any selection,
    indicates a technical difference in how the controls and targeting guides were made or
    sequenced, not biology.</li>
    <li>Divergence appearing only after selection supports a real, working screen — worth citing
    alongside <code>mageck test</code> results.</li>
    <li>A library with only a handful of control guides will show a jagged, noisy control curve
    from having so few data points — treat small deviations in that curve as less reliable than in
    a library with many control guides.</li>)",
    missing_reason = "this library has no control, NTC, or non-targeting guides in its Gene column."
  ),
  list(
    id = "90_10_ratio", svgs = "90_10_ratio",
    question = "6. How many times more reads do the best-covered guides get than the worst-covered?",
    title = "90/10 Percentile Ratio (Library Skew)",
    description = r"(The ratio of the top 10% of guides' read counts to the bottom 10%, across all
    of this library's guides, including guides with zero reads. A sample with more than 10% zero-count
    guides has a bottom-10% value of zero, which is counted as 1 read, so its ratio is the 90th
    percentile divided by 1: dropout raises the ratio.)",
    interpretation = r"(A ratio near 1 indicates very even representation. This report flags a
    ratio of 10 or higher as meaningfully skewed, and treats 5 or below as good — rules of thumb,
    not hard cutoffs.)",
    troubleshooting = r"(<li>A ratio already high at T0 indicates a library or coverage problem
    (see Guide Representation and Lorenz above); a ratio that's only high after selection is
    expected.</li>
    <li>A high ratio driven by many near-zero guides indicates a coverage problem; a high ratio
    driven by a few extremely high-count guides indicates PCR over-amplification.</li>)"
  ),
  list(
    id = "undetected_guides", svgs = "undetected_guides",
    question = "7. How many guides received no reads at all?",
    title = "Undetected Guides (Zero Counts)",
    description = "The percent (and count) of guides with exactly zero reads in each sample.",
    interpretation = r"(Some zero counts occur by chance even in a well-represented library, so the
    ideal is close to 0%, not necessarily exactly 0%. This report flags 5% or more as high dropout;
    under about 1% is a good target for a plasmid/reference sample. A value near 100% is a
    different, more severe problem than ordinary dropout: it usually means this sample doesn't
    belong to this library at all, rather than its guides genuinely failing one by one.)",
    troubleshooting = r"(<li>Confirm sequencing depth was sufficient for the library size first —
    aim for 300&ndash;500&times; coverage per guide.</li>
    <li>High dropout despite adequate depth at T0 indicates a problem upstream of sequencing:
    reconsider transduction MOI/coverage.</li>
    <li>A value near 100%: check the Mapping Rate plot below and confirm this sample was actually
    screened with this library, not a different one.</li>
    <li>The wrangle step's <code>zero_count_guides.csv</code> lists exactly which guides dropped
    out. Dropout clustered in one gene or pathway after treatment usually indicates real
    essential-gene dropout rather than a technical problem.</li>)"
  ),
  list(
    id = "mapping_rate", svgs = "mapping_rate",
    question = "8. What share of each sample's reads matches a guide of this library?",
    title = "Mapping Rate (% Exact Guide Match)",
    description = mapping_desc,
    interpretation = mapping_interp,
    troubleshooting = r"(<li>A low rate in one sample among otherwise-good ones indicates a
    sample-specific issue: check its barcode/index assignment, or whether <code>trim5</code>
    matches its actual read structure.</li>
    <li>A uniformly low rate across every sample indicates a pipeline-wide parameter mismatch:
    check <code>trim5</code> and <code>sgrna_len</code> in <code>config.yaml</code> against the
    real amplicon design, and check the barcode-position plot from demultiplexing QC.</li>
    <li>A nonzero-but-low rate that isn't a config problem usually indicates the sample was
    screened with a different library than the one it's assigned to — the Read Classification plot
    below (a high far-miss share) is a good way to confirm that.</li>)",
    missing_reason = "no usable mapping numbers could be read from the mageck count log for this library.",
    extra_html = mapping_table_html
  ),
  list(
    id = "read_classification", svgs = "read_classification",
    question = "9. Why didn't the other reads match a guide: sequencing errors or unrelated sequences?",
    title = "Read Classification: Exact / Near-Miss / Far-Miss",
    description = r"(Every sequenced read, split into three categories: an <strong>exact match</strong>
    to a library guide (what mageck count itself counts &mdash; the same number as Mapping Rate
    above), a <strong>near-miss</strong> (exactly one base different from some guide &mdash; not
    counted by mageck, but recovered here), or a <strong>far-miss</strong> (not within one base of
    any guide in the library). Each percentage is of that sample's own total reads. The three can add
    up to a little less than 100%: reads that mageck counted in the sample's total but listed neither
    as matched nor as unmapped (for example, reads too short to hold a whole guide) are not shown.)",
    interpretation = r"(A large exact-match share with very little near- or far-miss is the ideal
    case. A substantial near-miss share indicates real signal being lost to sequencing or PCR
    error: those reads almost certainly came from a real guide, just with one wrong base call. A
    substantial far-miss share indicates something else entirely &mdash; contamination, adapter-only
    reads, or a sample that doesn't belong to this library at all.)",
    troubleshooting = r"(<li>Mostly near-miss: check base-quality metrics in the FastQC/MultiQC
    report for this sample, and confirm <code>trim5</code>/<code>sgrna_len</code> exactly match the
    amplicon design (an off-by-one trim shifts every read by one base, turning what should be exact
    matches into near-misses across the whole library).</li>
    <li>Mostly far-miss: check whether this sample was actually screened with this library — see
    the Undetected Guides and Mapping Rate plots above for the same signal from another angle.</li>
    <li>A near-miss share that's concentrated in a few guides rather than spread evenly across the
    library can indicate a specific guide's sequence in <code>library.csv</code> doesn't match what
    was actually synthesized or cloned.</li>)",
    missing_reason = "this library's guides aren't all the same length, or one of --library/--mageck-log/--unmapped wasn't available, so reads couldn't be classified."
  ),
  list(
    id = "gc_vs_depth", svgs = "gc_vs_depth",
    question = "10. Does a guide's GC content bias how many reads it gets?",
    title = "GC Content vs Read Depth",
    description = "Each guide's GC content plotted against its normalized read count, with a smoothed trend line.",
    interpretation = r"(A flat trend line is the ideal case: read count doesn't depend on GC
    content. The steeper a monotonic slope, the stronger the GC-dependent amplification bias in
    that direction (high-GC or low-GC guides systematically over- or under-represented). A
    U-shaped or dome-shaped trend &mdash; dipping or peaking in the middle of the GC range instead
    of running in one direction &mdash; is the classic PCR-bias signature: guides at both GC
    extremes are disadvantaged (or both advantaged) relative to guides in the middle, rather than
    one extreme simply amplifying better than the other. Some slope is common even in good
    libraries; a very steep slope or a pronounced U/dome shape is the more specific signal to act
    on.)",
    troubleshooting = r"(<li>A trend already present at T0/plasmid, and consistent across samples,
    usually originates in the library PCR step (cycle number, polymerase, annealing temperature)
    rather than the biology — true whether the shape is a simple slope or a U/dome.</li>
    <li>mageck's normalization doesn't correct for GC bias, so it's worth noting even when it
    doesn't invalidate the screen.</li>
    <li>A trend that differs between control and treatment samples indicates GC content happens to
    correlate with a gene family in this library's design — worth a closer look.</li>)",
    missing_reason = "this library's MAGeCK-format file didn't include the guide id/sequence columns needed to compute GC content."
  ),
  list(
    id = "replicate_correlation", svgs = c("pearson_correlation", "spearman_correlation"),
    question = "11. Do samples, and especially replicates, agree with each other?",
    title = "Sample Correlation Heatmap (Pearson and Spearman)",
    description = corr_desc,
    interpretation = corr_interp,
    troubleshooting = r"(<li>Low correlation between nominal replicates indicates a possible sample
    mix-up in <code>samplesheet.csv</code> — check that before suspecting anything biological.</li>
    <li>Low correlation is very often explained by one replicate having much lower depth than its
    partner — check its mapping rate and % undetected guides above.</li>
    <li>High Spearman with lower Pearson for a pair indicates the guide ranking is preserved but
    the magnitude of a subset of guides differs — check GC Content and Guide Representation for
    that sample.</li>
    <li>A correlation at or extremely close to 1.0 between two nominally independent replicates is
    itself worth checking — it can indicate the same data was counted twice (a duplicated FASTQ, or
    two sample names pointing at the same file), not genuinely perfect reproducibility.</li>)",
    missing_reason = "this library has only one sample, so there's nothing to correlate it against."
  ),
  list(
    id = "pca", svgs = "pca",
    question = "12. Which samples resemble each other overall, and do replicates group together?",
    title = "Principal Component Analysis (PCA)",
    description = pca_desc,
    interpretation = pca_interp,
    troubleshooting = r"(<li>A sample far from its replicates: check its mapping rate, % undetected
    guides and read depth in the sections of this report, and check <code>samplesheet.csv</code>
    for a mix-up.</li>
    <li>All samples of one condition split into two groups: check for a batch effect or an index
    mix-up.</li>)",
    missing_reason = pca_missing_reason
  ),
  list(
    id = "read_depth", svgs = "read_depth",
    question = "13. Was each sample sequenced deeply enough?",
    title = "Read Depth",
    description = r"(For each sample, the grey bar is the total number of sequenced reads (from the
    mageck count log). The blue bar is the number of reads that matched a guide of this library. The
    scale is log10, with abbreviated labels (K = thousand, M = million).)",
    interpretation = paste0(r"(Depth should give enough reads per guide. The reads needed are roughly
    the number of guides times the reads wanted per guide. A common rule of thumb is at least
    200&ndash;500 reads per guide at the end point of a screen. Plasmid and T0 samples are often
    sequenced deeper. This report gives reads per guide as matched reads divided by the number of
    guides (column "Reads per Guide" in the At a Glance table). It is a rule of thumb, not a hard
    cutoff.)", rd_guides_txt,
    if (has_breakdown) r"( In this all-vs-all run the grey bar is the whole sample and the blue bar
    is the part that belongs to this library.)" else ""),
    troubleshooting = r"(<li>Low reads per guide at the end point: expect noisy counts with more
    dropout. Re-sequence deeper, or reduce the size of the pooled library.</li>
    <li>A large gap between the grey and blue bars: see Mapping Rate and Read Classification.</li>
    <li>Very different depth between replicates: expect a lower correlation between them.</li>)",
    missing_reason = "no read-depth table was found (the tables come from the wrangle step)."
  )
)

render_section <- function(sec) {
  found <- Filter(Negate(is.null), lapply(sec$svgs, read_svg))
  body <- if (length(found) > 0) {
    paste0("<div class='plot-row'>", paste(vapply(found, function(s) paste0("<div class='plot'>", s, "</div>"), character(1)), collapse = "\n"), "</div>")
  } else if (!is.null(sec$missing_reason)) {
    paste0("<div class='missing-box'><strong>Not generated for this library</strong> &mdash; ", sec$missing_reason, "</div>")
  } else {
    "<div class='missing-box'><strong>Not generated for this library.</strong></div>"
  }
  # Optional extra block (e.g. the per-library table) placed after the plot.
  extra <- if (is.null(sec$extra_html)) "" else sec$extra_html
  troubleshooting_html <- if (grepl("<li>", sec$troubleshooting, fixed = TRUE)) {
    paste0("<ul>", sec$troubleshooting, "</ul>")
  } else {
    paste0("<p>", sec$troubleshooting, "</p>")
  }
  paste0(
    "<section id='", sec$id, "'>",
    "<h2>", sec$question, "</h2><p class='subtitle'>", sec$title, "</p>",
    "<p class='desc'>", sec$description, "</p>",
    body,
    extra,
    "<div class='callout callout-interp'><h3>How to interpret it</h3><p>", sec$interpretation, "</p></div>",
    "<div class='callout callout-trouble'><h3>If this looks poor</h3>", troubleshooting_html, "</div>",
    "</section>"
  )
}

toc_items <- paste0("<li><a href='#", vapply(sections, `[[`, character(1), "id"), "'>",
                     vapply(sections, `[[`, character(1), "question"), "</a></li>", collapse = "\n")
sections_html <- paste(vapply(sections, render_section, character(1)), collapse = "\n")

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
.wrap { max-width: 980px; margin: 0 auto; padding: 2rem 1.5rem 4rem; }
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
.plot { flex: 1 1 420px; max-width: 100%; text-align: center; }
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
table.qc-table td.this-lib, table.qc-table th.this-lib { background: #eef3f8; }
footer { color: var(--muted); font-size: 0.82rem; margin-top: 3rem; border-top: 1px solid var(--border); padding-top: 1rem; }
@media print { nav.toc { display: none; } h2 { break-before: page; } }
)"

html <- paste0(
  "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width, initial-scale=1'>",
  "<title>", html_escape(project), " &mdash; ", html_escape(prefix), " guide-representation QC</title>",
  "<style>", css, "</style></head><body><div class='wrap'>",
  "<header class='page-header'><h1>Guide Representation QC Report</h1>",
  "<div class='meta'>Project: <strong>", html_escape(project), "</strong> &nbsp;&middot;&nbsp; ",
  "Library: <strong>", html_escape(prefix), "</strong> &nbsp;&middot;&nbsp; ",
  "Samples: <strong>", n_samples, "</strong> &nbsp;&middot;&nbsp; ",
  "Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), "</div></header>",
  "<div class='about'>This report covers sgRNA count-level QC for this library only: how evenly ",
  "its guides are represented, and (when it has more than one sample) how well they agree with ",
  "each other. It doesn't cover raw read QC (see <code>results/demultiplexing/qc/multiqc/</code>) ",
  "or statistical hit calling (see <code>results/mageck_test/</code>, if <code>contrasts.csv</code> ",
  "was set). The <strong>% Within 4-fold of Mean</strong> column in the table below is a quick ",
  "uniformity check: what percent of guides fall within a 4-fold range of the sample's average ",
  "count. Higher is more uniform &mdash; another way of looking at the same evenness the Gini ",
  "coefficient and 90/10 ratio measure.</div>",
  "<nav class='toc'><h2>Sections</h2><ul>", toc_items, "</ul></nav>",
  "<section id='summary'><h2>How does each sample look overall, and which ones need attention?</h2>",
  "<p class='subtitle'>At a Glance</p>", summary_table, "</section>",
  sections_html,
  "<footer>Generated by the CRISPR_SCREEN pipeline's guide_representation_report.R from tables in ",
  html_escape(tables_dir), " and plots in ", html_escape(plots_dir), ".</footer>",
  "</div></body></html>"
)

out_path <- file.path(outdir, paste0(prefix, "_qc_report.html"))
writeLines(html, out_path)
cat("  Saved:", out_path, "\n")
