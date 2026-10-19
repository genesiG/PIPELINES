#!/usr/bin/env Rscript
#' ##########################################################################
#' guide_representation_report.R
#'
#' Assembles one self-contained HTML QC report per library from the tables
#' guide_representation_wrangle.R wrote and the SVGs guide_representation_plot.R
#' rendered from them. For every plot: what it measures, how to read it, and
#' what to check next if it looks poor. No statistics are computed here —
#' this only reads existing tables/plots and writes HTML.
#'
#' USAGE:
#'   Rscript guide_representation_report.R --tables-dir=<dir> --plots-dir=<dir> \
#'       --prefix=<prefix> [--project=<name>] [--outdir=<dir>]
#'
#'   --tables-dir  directory holding <prefix>_*.csv from the wrangle step
#'   --plots-dir   directory holding <prefix>_*.svg from the plot step
#'   --prefix      count-table prefix (e.g. MYSCREEN_library)
#'   --project     project_name, shown in the report header (optional)
#'   --outdir      where to write the report (default: current directory)
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
  stop("Usage: guide_representation_report.R --tables-dir=<dir> --plots-dir=<dir> --prefix=<prefix> [--project=<name>] [--outdir=<dir>]")
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

# ---- "At a glance" summary table (every headline metric, one row/sample) ---
flags_col <- vapply(seq_len(nrow(qc)), function(i) {
  b <- character(0)
  if (isTRUE(qc$library_skew_flag[i]))          b <- c(b, badge("Library skew",  "flag"))
  if (isTRUE(qc$undetected_flag[i]))            b <- c(b, badge("High dropout",  "flag"))
  if (has_mapping && isTRUE(low_mapping[i]))    b <- c(b, badge("Low mapping",   "flag"))
  if (length(b) == 0) b <- badge("OK", "ok")
  paste(b, collapse = " ")
}, character(1))

summary_rows <- paste0(
  "<tr>",
  "<td>", html_escape(qc$sample), "</td>",
  if (has_mapping) paste0("<td class='num'>", fmt_int(qc$reads), "</td><td class='num'>", fmt_int(qc$mappedreads), "</td><td class='num'>", fmt_pct(qc$mapping_pct), "</td>") else "",
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
  if (has_class) "<th>% Near-Miss</th><th>% Far-Miss</th>" else "",
  "<th>Gini</th><th>90/10 Ratio</th><th>% Undetected (n)</th><th>% Within 4-fold of Mean</th><th>Flags</th>",
  "</tr></thead><tbody>",
  summary_rows,
  "</tbody></table>"
)

# ---- Plot sections: description / interpretation / troubleshooting --------
# Plain, short sentences — this report is read by wet-lab scientists, not
# statisticians. Each entry: id (report anchor), svgs (one or more
# <prefix>_<svg>.svg names), title, description (what it is), interpretation
# (what good/bad looks like, 2-4 short sentences), troubleshooting (what to
# check next, as a short bullet list), and an optional missing_reason shown
# in plain English instead of the plot when it wasn't generated.
# `avg_word`/`ecdf_*`/`density_*` below adapt wording for a single-sample
# library, where "average across samples" and "compare curves" don't apply.
avg_word <- if (one_sample) "read count" else "average read count across this library's samples"
density_desc <- if (one_sample) {
  "The shape of this sample's read-count distribution: where most guides sit, and whether there's a separate cluster elsewhere in the range."
} else {
  "The shape of each sample's read-count distribution, overlaid: where most guides sit, and whether there's a separate cluster elsewhere in the range."
}

sections <- list(
  list(
    id = "guide_representation", svgs = "guide_representation",
    title = "1. Guide Representation (Ranked)",
    description = paste0("Every guide's ", avg_word, ", sorted from lowest to highest, on a log scale. This is the standard check for how evenly a library is represented."),
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
    title = "2. ECDF of gRNA Representation",
    description = r"(The cumulative fraction of guides with a read count at or below a given value.
    The dashed vertical line marks the sample's mean read count &mdash; what every guide would have
    if reads were split perfectly evenly.)",
    interpretation = r"(In a perfectly even library, every guide has the same count, so the curve
    would rise in one sharp step right at the mean line. The more of the curve sitting to the left
    of the mean line, the larger the fraction of guides falling below their fair share of reads.)",
    troubleshooting = r"(<li>A wide gap between the curve and the mean line &mdash; most of the
    curve well to the left of it &mdash; indicates a large fraction of guides are underrepresented.
    Cross-check "% Undetected Guides" below for how many sit at exactly zero.</li>
    <li>An unusual step or kink in the curve's shape indicates the sample may be a mixture of two
    populations, e.g. cross-contamination between pools.</li>)"
  ),
  list(
    id = "lorenz", svgs = "lorenz",
    title = "3. Lorenz Curves",
    description = r"(Guides sorted from least to most reads, plotted as cumulative % of reads
    against cumulative % of guides.)",
    interpretation = r"(The diagonal line represents an ideal sample with perfectly even gRNA
    representation: every guide contributing an identical share of reads. The further the curve
    bows away from the diagonal, the less evenly the library is represented. This is the same
    thing the Gini coefficient in the table above measures as one number &mdash; Gini is defined as
    twice the area between the curve and the diagonal.)",
    troubleshooting = r"(<li>A curve already bowed at T0/plasmid indicates a library or coverage
    problem — see Guide Representation above.</li>
    <li>Bowing that increases from T0 through treated/selected timepoints reflects selection
    pressure concentrating reads on surviving guides — expected, not a technical artifact.</li>
    <li>Confirm the Gini and 90/10 ratio (below) tell the same story.</li>)"
  ),
  list(
    id = "density", svgs = "density",
    title = "4. Read-Count Density",
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
    title = "5. Control vs Targeting gRNA Representation",
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
    title = "6. 90/10 Percentile Ratio (Library Skew)",
    description = r"(The ratio of the top 10% of guides' read counts to the bottom 10%, for
    guides with at least one read.)",
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
    title = "7. Undetected Guides (Zero Counts)",
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
    title = "8. Mapping Rate (% Exact Guide Match)",
    description = r"(The percent of sequenced reads that matched a library guide, from the mageck
    count log. mageck count's own matching is always exact &mdash; it has no mismatch-tolerance
    option &mdash; so this is also, exactly, the percentage of reads that exactly matched a guide.
    The Read Classification plot below breaks the rest (the non-matching reads) down further.)",
    interpretation = r"(A healthy pooled screen typically maps 60&ndash;70% or more of its reads;
    this report flags anything under 30%. Mapping rate is independent of the representation metrics
    above: a sample can map poorly yet still show even representation among the reads that did map
    — that just means most of its sequencing depth was wasted, not misrepresented.)",
    troubleshooting = r"(<li>A low rate in one sample among otherwise-good ones indicates a
    sample-specific issue: check its barcode/index assignment, or whether <code>trim5</code>
    matches its actual read structure.</li>
    <li>A uniformly low rate across every sample indicates a pipeline-wide parameter mismatch:
    check <code>trim5</code> and <code>sgrna_len</code> in <code>config.yaml</code> against the
    real amplicon design, and check the barcode-position plot from demultiplexing QC.</li>
    <li>A nonzero-but-low rate that isn't a config problem usually indicates the sample was
    screened with a different library than the one it's assigned to — the Read Classification plot
    below (a high far-miss share) is a good way to confirm that.</li>)",
    missing_reason = "no usable mapping numbers could be read from the mageck count log for this library."
  ),
  list(
    id = "read_classification", svgs = "read_classification",
    title = "9. Read Classification: Exact / Near-Miss / Far-Miss",
    description = r"(Every sequenced read, split into three categories: an <strong>exact match</strong>
    to a library guide (what mageck count itself counts &mdash; the same number as Mapping Rate
    above), a <strong>near-miss</strong> (exactly one base different from some guide &mdash; not
    counted by mageck, but recovered here), or a <strong>far-miss</strong> (not within one base of
    any guide in the library).)",
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
    title = "10. GC Content vs Read Depth",
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
    title = "11. Replicate Correlation",
    description = r"(How well every pair of samples in this pool agree with each other:
    <strong>Pearson</strong> correlation (on log-transformed counts, sensitive to overall count
    magnitude) and <strong>Spearman</strong> correlation (rank-based, robust to outliers).)",
    interpretation = r"(True replicates correlate highest of any pair in the pool — above 0.9 is
    typical for a technically sound screen. Lower correlation between different conditions or
    timepoints is expected, not a problem, provided replicates of the same condition still agree
    closely with each other.)",
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
  troubleshooting_html <- if (grepl("<li>", sec$troubleshooting, fixed = TRUE)) {
    paste0("<ul>", sec$troubleshooting, "</ul>")
  } else {
    paste0("<p>", sec$troubleshooting, "</p>")
  }
  paste0(
    "<section id='", sec$id, "'>",
    "<h2>", sec$title, "</h2>",
    "<p class='desc'>", sec$description, "</p>",
    body,
    "<div class='callout callout-interp'><h3>How to interpret it</h3><p>", sec$interpretation, "</p></div>",
    "<div class='callout callout-trouble'><h3>If this looks poor</h3>", troubleshooting_html, "</div>",
    "</section>"
  )
}

toc_items <- paste0("<li><a href='#", vapply(sections, `[[`, character(1), "id"), "'>",
                     vapply(sections, `[[`, character(1), "title"), "</a></li>", collapse = "\n")
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
  "<section id='summary'><h2>At a Glance</h2>", summary_table, "</section>",
  sections_html,
  "<footer>Generated by the CRISPR_SCREEN pipeline's guide_representation_report.R from tables in ",
  html_escape(tables_dir), " and plots in ", html_escape(plots_dir), ".</footer>",
  "</div></body></html>"
)

out_path <- file.path(outdir, paste0(prefix, "_qc_report.html"))
writeLines(html, out_path)
cat("  Saved:", out_path, "\n")
