#!/usr/bin/env Rscript
# diff_annotate.R -- annotate the regions of ONE differential table with ChIPseeker (TxDb + org.db) and write the all /
# promoter / non-promoter / distal tables and BEDs, the logFC histogram and a per-set count summary.
# Port of the legacy peak_annotation.Rmd chunk "loop and annotate" (Scripts/peak_annotation.Rmd, AV_ATACSEQ_2026, March 2026).
# The computation is the legacy code (read.csv, makeGRangesFromDataFrame, keepStandardChromosomes, seqlengths from the
# BSgenome, genome<-, seqlevelsStyle<-, annotatePeak with tssRegion and annoDb, as.data.frame, the grepl subsets, the BED
# export); only the I/O is replaced by the CLI. The TxDb is the one built by txdb_build.R from the same GTF (the legacy read
# rds/txdb.hs1.sqlite, built from the same file).
#
# CLI contract (called by the DIFF_ANNOTATE module; --key=value):
#   --name=NAME                   REQUIRED. Comparison name = output prefix (e.g. BT54_high_vs_low). ^[A-Za-z0-9][A-Za-z0-9_.-]*$
#   --table=F.csv                 REQUIRED. Differential table written by diff_accessibility.R (NAME.csv, or NAME_promoters.csv).
#                                 Needs the columns seqnames, start, end, strand, logFC and FDR. Zero rows: exit 3.
#   --txdb=F                      REQUIRED. TxDb SQLite (txdb_build.R). Missing file: exit 3.
#   --species="Homo sapiens"      OPTIONAL, default "Homo sapiens". Only sets the defaults of --bsgenome and --org-db.
#   --genome-build=hs1            OPTIONAL, default hs1. Written to the GRanges with genome<- (legacy build label).
#   --bsgenome=PKG                OPTIONAL. BSgenome package that gives the seqlengths. Default BSgenome.Hsapiens.UCSC.<build>
#                                 (species Mus musculus: BSgenome.Mmusculus.UCSC.<build>). Missing package: exit 3.
#   --org-db=PKG                  OPTIONAL. annoDb of annotatePeak (gene symbols). Default org.Hs.eg.db (Mus musculus: org.Mm.eg.db).
#   --tss-region=-2000,200        OPTIONAL, default -2000,200. ChIPseeker tssRegion: two integers, upstream and downstream of the TSS.
#   --fdr=                        OPTIONAL. Blank (default) = annotate ALL regions of the table (legacy AV_ATACSEQ_2026 behaviour).
#                                 A number x: annotate the regions with FDR < x; when none pass, all regions are annotated
#                                 with a warning (older legacy behaviour, which filtered FDR < 0.1 before annotating).
#   --cpus=1                      OPTIONAL, default 1. Accepted for the module contract; annotatePeak is serial (no parallel loop).
#   --font=Helvetica              OPTIONAL, default Helvetica. Font family of the plot (svg: font-family of the text elements).
#   --plot-format=svg|pdf         OPTIONAL, default svg. svg: NAME_peak_distribution.svg (svglite, editable <text>). pdf: the
#                                 legacy NAME_peak_distribution.pdf (legacy device arguments).
#   --dist-width=5                OPTIONAL, default 5 (inches). Width of the logFC histogram.
#   --dist-height=5               OPTIONAL, default 5 (inches). Height of the logFC histogram.
#
# Outputs (cwd, prefix = --name):
#   NAME_all_sites.csv            every annotated region (write.csv, row names as legacy; columns of annotatePeak as.data.frame)
#   NAME_all_sites.bed            the same regions as BED (rtracklayer::export.bed)
#   NAME_promoter_sites.csv/.bed  regions whose annotation contains "Promoter"
#   NAME_non-promoter_sites.csv/.bed  regions whose annotation does not contain "Promoter"
#   NAME_distal_sites.csv/.bed    regions whose annotation contains "Distal Intergenic"
#   NAME_peak_distribution.svg    histogram of logFC (100 bins; --dist-width x --dist-height in, default 5 x 5) of the annotated
#                                 regions (--plot-format=pdf: NAME_peak_distribution.pdf). The legacy called hist() outside
#                                 pdf(), so its PDF was empty; here hist() is drawn inside the device. Same histogram.
#   NAME.annotation_summary.tsv   set, n_regions: all_sites, promoter_sites, non-promoter_sites, distal_sites
#
# Exit codes: 0 ok; 2 usage error (incl. a bad --plot-format or plot size); 3 known data problem: table or TxDb missing,
# table with no rows, a BSgenome or org package not installed, or the TxDb cannot be loaded.
#
# Interpretation limits (no statistic is computed here; the table's PValue / FDR are those of diff_accessibility.R):
# * By default ALL regions of the table are annotated, significant or not, so the promoter / non-promoter / distal sets and
#   their counts describe the tested region set, not differentially accessible regions. Use the FDR and logFC columns of the
#   tables for that (or --fdr, which drops non-significant regions before annotating).
# * geneId / SYMBOL are the NEAREST TSS gene as ChIPseeker assigns it, with no distance limit (read distanceToTSS): a distal
#   intergenic region hundreds of kb from a gene still carries that gene. With the RefSeq GTF TxDb the gene_id is the gene symbol.
#   A region that overlaps a promoter and another feature is labelled by ChIPseeker's priority order (Promoter first). An
#   annotation is not evidence that the region regulates that gene.
# * Regions on contigs that are not standard chromosomes of the TxDb style are dropped by keepStandardChromosomes (coarse).
#
# Not ported: the legacy final chunk "annotate genome" (promoters/non-promoters of the whole genome BED). It never ran
# (it referenced non_promoters_gr, while the object was named non_promoter_gr) and wrote no file that any later step reads.
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
# One plot written to BASE.svg (--plot-format=svg) or BASE.pdf (pdf: the legacy device). Pair with plot_close(BASE).
plot_open <- function(base, width, height) {
  if (plot_format == "svg") {
    open_svg(paste0(base, ".svg"), width, height)
  } else if (tolower(plot_font) %in% c("helvetica", "arial", "sans")) {
    grDevices::pdf(file = paste0(base, ".pdf"), width = width, height = height, family = "Helvetica")
  } else {
    grDevices::cairo_pdf(filename = paste0(base, ".pdf"), width = width, height = height, family = plot_font)
  }
  invisible(base)
}
plot_close <- function(base) {
  if (plot_format == "svg") close_svg(paste0(base, ".svg")) else grDevices::dev.off()
  invisible(NULL)
}

allowed <- c("name", "table", "txdb", "species", "genome-build", "bsgenome", "org-db", "tss-region", "fdr", "cpus",
             "font", "plot-format", "dist-width", "dist-height")
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) == 0L || any(argv %in% c("-h", "--help"))) {
  cat("usage: diff_annotate.R --name=NAME --table=NAME.csv --txdb=txdb.sqlite [--species=...] [--genome-build=hs1]\n",
      "       [--bsgenome=PKG] [--org-db=PKG] [--tss-region=-2000,200] [--fdr=X] [--cpus=1]\n",
      "       [--font=Helvetica] [--plot-format=svg|pdf] [--dist-width=5] [--dist-height=5]\n", file = stderr())
  quit(save = "no", status = 2L)
}
args <- parse_args(argv, allowed)
name       <- get_req(args, "name")
table_file <- get_req(args, "table")
txdb_file  <- get_req(args, "txdb")
species    <- get_opt(args, "species", "Homo sapiens")
build      <- get_opt(args, "genome-build", "hs1")
if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", name)) die(paste0("--name must match ^[A-Za-z0-9][A-Za-z0-9_.-]*$, got '", name, "'"), 2L)

if (species == "Mus musculus") {
  default_bsgenome <- paste0("BSgenome.Mmusculus.UCSC.", build)
  default_org_db   <- "org.Mm.eg.db"
} else {
  default_bsgenome <- paste0("BSgenome.Hsapiens.UCSC.", build)
  default_org_db   <- "org.Hs.eg.db"
}
bsgenome_pkg <- get_opt(args, "bsgenome", default_bsgenome)
annoDb       <- get_opt(args, "org-db", default_org_db)

tss_region <- suppressWarnings(as.integer(split_list(get_opt(args, "tss-region", "-2000,200"))))
if (length(tss_region) != 2L || anyNA(tss_region)) die("--tss-region must be two integers, e.g. -2000,200", 2L)

fdr_arg <- get_opt(args, "fdr", "")
fdr_cut <- NA_real_
if (nzchar(fdr_arg)) {
  fdr_cut <- suppressWarnings(as.numeric(fdr_arg))
  if (is.na(fdr_cut)) die(paste0("--fdr must be a number or blank, got '", fdr_arg, "'"), 2L)
}
cpus <- suppressWarnings(as.integer(get_opt(args, "cpus", "1")))
if (is.na(cpus) || cpus < 1L) die("--cpus must be a positive integer", 2L)
plot_font   <- get_opt(args, "font", "Helvetica")
plot_format <- tolower(get_opt(args, "plot-format", "svg"))
dist_width  <- suppressWarnings(as.numeric(get_opt(args, "dist-width", "5")))
dist_height <- suppressWarnings(as.numeric(get_opt(args, "dist-height", "5")))
if (!(plot_format %in% c("svg", "pdf"))) die("--plot-format must be svg or pdf", 2L)
if (is.na(dist_width) || is.na(dist_height) || dist_width <= 0 || dist_height <= 0) {
  die("--dist-width and --dist-height must be numbers > 0", 2L)
}

if (!file.exists(table_file)) die(paste0("differential table not found: ", table_file), 3L)
if (!file.exists(txdb_file)) die(paste0("TxDb not found: ", txdb_file, " (run TXDB_BUILD)"), 3L)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(GenomicFeatures)
  library(ChIPseeker)
  library(rtracklayer)
})
for (p in c(bsgenome_pkg, annoDb)) {
  if (!requireNamespace(p, quietly = TRUE)) die(paste0("package not installed in this image: ", p), 3L)
}
txdb <- tryCatch(
  AnnotationDbi::loadDb(file = txdb_file),
  error = function(e) die(paste0("could not load TxDb ", txdb_file, ": ", conditionMessage(e)), 3L)
)
genome <- BSgenome::getBSgenome(bsgenome_pkg)

comp <- gsub("_", " ", name, fixed = TRUE)
cat(paste0("\nPerforming annotation: ", comp, "\n"))

# LOAD PEAKS
peak_df <- read.csv(table_file, header = TRUE)
if (nrow(peak_df) == 0L) die(paste0("differential table has no rows: ", table_file), 3L)

# FILTER PEAKS by FDR (optional; blank --fdr annotates all regions, as the legacy AV_ATACSEQ_2026 code did)
peakdf <- peak_df
if (!is.na(fdr_cut)) {
  peakdf <- subset(peak_df, FDR < fdr_cut)
  if (nrow(peakdf) == 0L) {
    message(paste0("WARNING: no region with FDR < ", fdr_cut, " for ", comp,
                   "; annotating all regions instead"))
    peakdf <- peak_df
  } else {
    cat(paste0("Regions with FDR < ", fdr_cut, ": ", nrow(peakdf), " of ", nrow(peak_df), "\n"))
  }
}

# Plot histogram of logFC (drawn inside the plot device: NAME_peak_distribution.svg, or .pdf with --plot-format=pdf)
xlab = paste0("Differential Accessibility\n(log2 fold-change ", comp, ")")
dist_base <- paste0(name, "_peak_distribution")
plot_open(dist_base, dist_width, dist_height)
hist(peakdf$logFC,
     100,
     col = "gray",
     main = "",
     xlab = xlab,
     ylab = "Number of Peaks")
plot_close(dist_base)

# Convert to GRanges
peaks <- makeGRangesFromDataFrame(peakdf, keep.extra.columns = TRUE)

# ANNOTATE PEAKS
peaks <- GenomeInfoDb::keepStandardChromosomes(peaks, pruning.mode = "coarse")
common <- intersect(seqlevels(peaks), seqlevels(genome))
seqlengths(peaks)[common] <- seqlengths(genome)[common]
genome(peaks) <- build
seqlevelsStyle(peaks) <- seqlevelsStyle(txdb)[1]

peakAnno <- annotatePeak(peaks,
                         TxDb      = txdb,
                         tssRegion = tss_region,
                         annoDb    = annoDb)

# Convert annotations to data frame
anno_df <- as.data.frame(peakAnno)

write.csv(anno_df, file = paste0(name, "_all_sites.csv"))
rtracklayer::export.bed(makeGRangesFromDataFrame(anno_df, keep.extra.columns = TRUE),
                        con = paste0(name, "_all_sites.bed"))

## Promoter Peaks
promoter_df <- subset(anno_df, grepl("Promoter", annotation))
promoter_peaks <- GenomicRanges::makeGRangesFromDataFrame(promoter_df,
                                                          keep.extra.columns = TRUE)
write.csv(promoter_df, file = paste0(name, "_promoter_sites.csv"))
rtracklayer::export.bed(promoter_peaks,
                        con = paste0(name, "_promoter_sites.bed"))

## Non-Promoter Peaks
nonpromoter_df <- subset(anno_df, !grepl("Promoter", annotation))
nonpromoter_peaks <- GenomicRanges::makeGRangesFromDataFrame(nonpromoter_df,
                                                             keep.extra.columns = TRUE)
write.csv(nonpromoter_df, file = paste0(name, "_non-promoter_sites.csv"))
rtracklayer::export.bed(nonpromoter_peaks,
                        con = paste0(name, "_non-promoter_sites.bed"))

## Distal Intergenic Peaks
distal_df <- subset(anno_df, grepl("Distal Intergenic", annotation))
distal_peaks <- GenomicRanges::makeGRangesFromDataFrame(distal_df,
                                                        keep.extra.columns = TRUE)
write.csv(distal_df, file = paste0(name, "_distal_sites.csv"))
rtracklayer::export.bed(distal_peaks, con = paste0(name, "_distal_sites.bed"))

# Report peak counts
cat("\nPromoter peaks:", length(promoter_peaks), "\n")
cat("Non-promoter peaks:", length(nonpromoter_peaks), "\n")
cat("Distal peaks:", length(distal_peaks), "\n")

summary_df <- data.frame(
  set       = c("all_sites", "promoter_sites", "non-promoter_sites", "distal_sites"),
  n_regions = c(nrow(anno_df), length(promoter_peaks), length(nonpromoter_peaks), length(distal_peaks))
)
write.table(summary_df, file = paste0(name, ".annotation_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(paste0("\nAnnotation complete: ", name, "\n"))
