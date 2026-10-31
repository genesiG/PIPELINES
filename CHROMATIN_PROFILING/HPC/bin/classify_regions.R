#!/usr/bin/env Rscript
# classify_regions.R -- classify pre-specified regions as RETAINED / LOST / UNCLASSIFIED between a reference and a treated group
# from normalised read counts (csaw::regionCounts + edgeR quasi-likelihood), using scale factors from the normalisation step.
#
# CLI contract (called by the CLASSIFY_REGIONS module; R style --key=value):
#   --regions=<bed>             regions to classify: BED3+ (4th column = name, used when present; else region_<i>), 0-based starts
#   --ref-bams=<comma list>     reference-group BAMs (final filtered, coordinate-sorted + indexed)
#   --treat-bams=<comma list>   treated-group BAMs
#   --ref-ids=<comma list>      sample ids of the reference BAMs (same order); match sample_id of the scale-factor table
#   --treat-ids=<comma list>    sample ids of the treated BAMs (same order)
#   --scale-factors=<tsv>       <norm_run>.scale_factors.tsv; ONLY the columns sample_id and edger_norm_factor are read
#   --id=<name>                 classification id = prefix of every output file
#   --lfc-retained=<num>        retained when round(logFC, 1) >= this (e.g. -1)
#   --lfc-lost=<num>            lost when logFC <= this AND FDR < --fdr (e.g. -2)
#   --fdr=<num>                 FDR threshold (BH over all input regions)
#   --paired-end=TRUE|FALSE     csaw readParam(pe = "both" | "none")
#   --minq=<int>                csaw readParam(minq=) (pipeline `mapq`)
#   --dedup=TRUE|FALSE          csaw readParam(dedup=) (pipeline `remove_duplicates`; TRUE drops reads flagged 0x400)
#   --max-frag=<int>            csaw readParam(max.frag=) (same value as the peak calling: csaw_max_frag, else bowtie2_maxins (500 if unset))
#   OPTIONAL: --outdir=.        output directory (default current directory)
#             --ext=<int>       single-end only: fixed read extension in bp (default 100 = csaw/legacy default); ignored for paired-end
#             --norm-run=<name> recorded in the summary only
#             --cpus=<int>      BiocParallel workers for counting (default 1)
# Outputs (all in --outdir, prefix = --id):
#   <id>.retained.bed | <id>.lost.bed | <id>.unclassified.bed
#       BED6 with 0-BASED starts (chr, start0, end, name, score = logFC rounded to 4 decimals, "."), input order. Always written; an
#       empty class gives an EMPTY file (0 bytes). A logFC of NA (no estimate) is written as score 0 in the BED (per_region keeps NA).
#   <id>.per_region.tsv.gz   one row per input region, input order: chr, start0, end, name, logFC, logCPM, PValue, FDR, class
#   <id>.summary.tsv         one row per class (retained, lost, unclassified) + "total": id, class, n, Mb, lfc_retained, lfc_lost, fdr,
#                            n_ref, n_treat, ref_ids, treat_ids, norm_run, n_overlapping_regions, dispersion, paired_end, minq, dedup,
#                            max_frag
# Exit codes: 0 ok; 2 usage error; 3 data problem (missing/NaN scale factor, empty BAM list, missing BAM/index, regions on chromosomes
#   absent from a BAM, malformed BED, empty library).
#
# ---- Question, model, assumptions --------------------------------------------------------------------------------------
# Question: by how much does the signal in each pre-specified region change between treated and reference libraries, and is the
# change distinguishable from replicate variability? Reads (fragments) are counted per region and library; an edgeR NB GLM
# ~ group (reference first, intercept + treated effect) with quasi-likelihood F-test of the treated coefficient (coef = 2); logFC =
# log2(treated / reference), FDR = Benjamini-Hochberg over ALL input regions.
# Normalisation: edgeR norm.factors = edger_norm_factor of the normalisation run (spike-in, ChIPseqSpikeInFree or TMM), library size =
# csaw's full BAM total (counts$totals, never the column sums over the regions). The factor is an external, composition-aware
# scaling: it must not be re-estimated from these regions (which are enriched by construction, so TMM on them would assume what is
# being tested). The statistics are only as valid as that factor: under a global loss of signal a data-driven factor would remove the
# effect; this script never computes one.
# Invariance to the centring constant: the pipeline's edger_norm_factor = SF / geomean(SF) (geometric mean over the run = 1), the
# legacy classifier used the raw SF. Multiplying every norm.factor by one constant c multiplies every effective library size by c,
# i.e. adds the constant -log(c) to every offset of the log-linear model; the intercept absorbs it, so fitted counts, logFC, p-values,
# dispersions and QL statistics are unchanged (up to numerical error ~1e-6); only logCPM shifts by -log2(c). Classes depend only on
# logFC and FDR, so the classification is invariant. (Per-sample ratios between factors are what matters; direction: divisor of the
# library size, as edgeR expects.)
# Rules (legacy classify_regions_lfc.R, GG_EED_CUTANDRUN): lost = logFC <= lfc_lost & FDR < fdr; retained = round(logFC, 1) >=
# lfc_retained | FDR >= fdr (non-significant regions count as retained: "not shown to be lost", NOT "shown to be retained"); unclassified
# = the rest; a region that is both lost and retained is LOST. NA logFC/FDR never satisfy a condition (an NA FDR with logFC >=
# lfc_retained is retained, as in the legacy script). Inference: the lfc cut-offs are descriptive; FDR controls the false calls for
# "differs from reference", not for "|change| exceeds 2^lfc". `retained` is an absence of evidence of loss and depends on power
# (replicates, depth); report it as "no significant loss" unless equivalence is shown separately.
# Dispersion: needs >= 1 residual degree of freedom (n_ref + n_treat > 2). With 1 + 1 libraries the NB dispersion cannot be estimated:
# fixed common dispersion 0.05 with a likelihood-ratio test (the fallback of callpeaks_csaw.R); the FDR is then NOT valid and a WARNING
# is printed. With 2 + 1 (or 1 + 2) it is estimated from the replicated group only (WARNING).
# Overlapping input regions: a read in an overlap is counted in each region; the tests are then not independent (BH FDR is
# approximate). A warning gives the number of regions that overlap another one (merge to disjoint regions upstream if it is large).
#
# ---- Differences from the legacy classify_regions_lfc.R ------------------------------------------------------------------
# * BED starts are 0-based (legacy wrote as.data.frame(gr)$start, 1-based: every BED was shifted by 1 bp).
# * regionCounts uses readParam(minq, dedup, max.frag, pe) given on the command line (legacy: csaw defaults minq 0, dedup FALSE,
#   max.frag 500). With duplicate-removed BAMs dedup is a no-op; minq/max.frag matter only if reads/fragments are filtered by them.
# * norm.factors from edger_norm_factor (centred, see above) instead of raw SF; inputs by id, not by BAM filename.
# * lost wins when both rules are true (legacy would list the region in both BEDs); empty class -> empty file (legacy: no file);
#   per-region table and summary are written; no reticulate/config.py, no dplyr/rtracklayer.

options(stringsAsFactors = FALSE)

die <- function(msg, status = 3L) {
  cat("ERROR: ", msg, "\n", sep = "", file = stderr())
  quit(save = "no", status = status)
}
warn <- function(msg) cat("WARNING: ", msg, "\n", sep = "", file = stderr())
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
as_bool <- function(x, key) {
  u <- toupper(trimws(x))
  if (u %in% c("TRUE", "T", "1", "YES")) TRUE
  else if (u %in% c("FALSE", "F", "0", "NO")) FALSE
  else die(paste0("--", key, " must be TRUE or FALSE, got '", x, "'"), 2L)
}
check_bams <- function(bams) {
  miss <- bams[!file.exists(bams)]
  if (length(miss)) die(paste0("BAM file(s) not found: ", paste(miss, collapse = ", ")))
  has_idx <- file.exists(paste0(bams, ".bai")) | file.exists(sub("\\.bam$", ".bai", bams))
  if (!all(has_idx)) die(paste0("BAM index (.bai) missing for: ", paste(bams[!has_idx], collapse = ", ")))
}

suppressPackageStartupMessages({
  library(csaw)
  library(edgeR)
  library(GenomicRanges)
  library(BiocParallel)
})

# ---- arguments ----------------------------------------------------------------------------------------------------------
allowed <- c("regions", "ref-bams", "treat-bams", "ref-ids", "treat-ids", "scale-factors", "id", "lfc-retained", "lfc-lost",
             "fdr", "paired-end", "minq", "dedup", "max-frag", "outdir", "ext", "norm-run", "cpus")
args <- parse_args(commandArgs(trailingOnly = TRUE), allowed)
regions_f   <- get_req(args, "regions")
# an absent list argument is a usage error (2); a present but empty list is a data problem (3: the group has no libraries)
get_list <- function(args, key) {
  if (is.null(args[[key]])) die(paste0("missing required argument --", key), 2L)
  v <- split_list(args[[key]])
  if (!length(v)) die(paste0("--", key, " is an empty list: the group has no libraries"))
  v
}
ref_bams    <- get_list(args, "ref-bams")
treat_bams  <- get_list(args, "treat-bams")
ref_ids     <- get_list(args, "ref-ids")
treat_ids   <- get_list(args, "treat-ids")
sf_file     <- get_req(args, "scale-factors")
id          <- get_req(args, "id")
lfc_ret     <- as_num(get_req(args, "lfc-retained"), "lfc-retained")
lfc_lost    <- as_num(get_req(args, "lfc-lost"), "lfc-lost")
fdr_thr     <- as_num(get_req(args, "fdr"), "fdr")
paired      <- as_bool(get_req(args, "paired-end"), "paired-end")
minq        <- as_num(get_req(args, "minq"), "minq")
dedup       <- as_bool(get_req(args, "dedup"), "dedup")
max_frag    <- as_num(get_req(args, "max-frag"), "max-frag")
outdir      <- get_opt(args, "outdir", ".")
ext_se      <- as_num(get_opt(args, "ext", "100"), "ext")
norm_run    <- get_opt(args, "norm-run", "")
cpus        <- as_num(get_opt(args, "cpus", "1"), "cpus")
if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", id)) die("--id must match ^[A-Za-z0-9][A-Za-z0-9_.-]*$", 2L)
if (fdr_thr <= 0 || fdr_thr > 1) die("--fdr must be in (0, 1]", 2L)
if (cpus < 1) die("--cpus must be >= 1", 2L)
if (ext_se < 1) die("--ext must be >= 1", 2L)
if (lfc_lost > lfc_ret) {
  warn(sprintf("--lfc-lost (%g) is above --lfc-retained (%g): regions between them are both lost and retained and are called lost", lfc_lost, lfc_ret))
}
if (length(ref_bams) < 1L || length(treat_bams) < 1L) die("need at least one reference BAM and one treated BAM (empty BAM list)")
if (length(ref_bams) != length(ref_ids)) die("--ref-bams and --ref-ids differ in length")
if (length(treat_bams) != length(treat_ids)) die("--treat-bams and --treat-ids differ in length")
bams <- c(ref_bams, treat_bams); ids <- c(ref_ids, treat_ids)
if (anyDuplicated(ids)) die(paste0("sample ids used more than once (a sample cannot be in both groups): ", paste(unique(ids[duplicated(ids)]), collapse = ", ")))
if (anyDuplicated(bams)) die("the same BAM is listed more than once")
check_bams(bams)
n_ref <- length(ref_bams); n_treat <- length(treat_bams)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
prefix <- file.path(outdir, id)

# ---- scale factors (only sample_id and edger_norm_factor are read) ----------------------------------------------------------
if (!file.exists(sf_file)) die(paste0("scale-factor table not found: ", sf_file))
sf <- read.delim(sf_file, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character")
if (!all(c("sample_id", "edger_norm_factor") %in% colnames(sf))) {
  die(paste0(sf_file, " lacks the column(s): ", paste(setdiff(c("sample_id", "edger_norm_factor"), colnames(sf)), collapse = ", ")))
}
pos <- match(ids, sf$sample_id)
if (anyNA(pos)) {
  die(paste0("no scale factor for sample(s) ", paste(ids[is.na(pos)], collapse = ", "), " in ", sf_file,
             ": the normalisation run must contain every reference and treated sample"))
}
nf <- suppressWarnings(as.numeric(sf$edger_norm_factor[pos]))
if (any(!is.finite(nf)) || any(nf <= 0)) {
  bad <- ids[!is.finite(nf) | nf <= 0]
  die(paste0("edger_norm_factor is missing, NaN or not positive for: ", paste(bad, collapse = ", "),
             " (normalisation QC failed for these libraries?); fix the normalisation run, no default of 1 is used"))
}

# ---- regions --------------------------------------------------------------------------------------------------------------------
if (!file.exists(regions_f)) die(paste0("regions BED not found: ", regions_f))
ln <- readLines(regions_f, warn = FALSE)
ln <- ln[nzchar(trimws(ln)) & !grepl("^(#|track|browser)", ln)]
if (!length(ln)) die(paste0("regions BED is empty: ", regions_f))
bd <- tryCatch(utils::read.table(text = ln, sep = "\t", quote = "", comment.char = "", colClasses = "character", fill = TRUE,
                                 header = FALSE, stringsAsFactors = FALSE),
               error = function(e) NULL)
if (is.null(bd) || ncol(bd) < 3L) {
  bd <- tryCatch(utils::read.table(text = ln, sep = "", quote = "", comment.char = "", colClasses = "character", fill = TRUE,
                                   header = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
}
if (is.null(bd) || ncol(bd) < 3L) die(paste0("cannot parse ", regions_f, " as BED (need chr, start, end)"))
b_chr <- bd[[1]]
b_start <- suppressWarnings(as.numeric(bd[[2]])); b_end <- suppressWarnings(as.numeric(bd[[3]]))
if (anyNA(b_start) || anyNA(b_end) || any(b_start < 0) || any(b_end <= b_start)) {
  die(paste0("malformed BED coordinates in ", regions_f, " (need integers with 0 <= start < end)"))
}
n_reg <- length(b_chr)
b_name <- if (ncol(bd) >= 4L) bd[[4]] else rep(NA_character_, n_reg)
b_name[is.na(b_name) | !nzchar(b_name) | b_name == "."] <- paste0("region_", seq_len(n_reg))[is.na(b_name) | !nzchar(b_name) | b_name == "."]
if (n_reg < 10L) die(paste0("only ", n_reg, " regions in ", regions_f, ": too few to fit dispersion trends (need >= 10)"))
gr <- GRanges(seqnames = b_chr, ranges = IRanges(start = as.integer(b_start) + 1L, end = as.integer(b_end)))  # BED 0-based -> 1-based
n_ov <- {
  ov <- findOverlaps(gr, gr, ignore.strand = TRUE)
  length(unique(S4Vectors::queryHits(ov)[S4Vectors::queryHits(ov) != S4Vectors::subjectHits(ov)]))
}
if (n_ov > 0L) {
  warn(sprintf("%d of %d input regions (%.1f%%) overlap at least one other input region: reads in an overlap are counted in each region and the tests are not independent",
               n_ov, n_reg, 100 * n_ov / n_reg))
}
for (b in bams) {
  tg <- names(Rsamtools::scanBamHeader(b)[[1]]$targets)
  miss <- setdiff(unique(b_chr), tg)
  if (length(miss)) {
    die(paste0("regions on chromosome(s) not in the header of ", b, ": ", paste(head(miss, 5), collapse = ", "),
               if (length(miss) > 5) ", ..." else "", " (chromosome naming of the regions vs the BAMs differs?)"))
  }
}

# ---- counting ------------------------------------------------------------------------------------------------------------------
cat("=== classify_regions.R: ", id, " ===\n", sep = "")
cat("reference: ", paste(ref_ids, collapse = ", "), "\ntreated:   ", paste(treat_ids, collapse = ", "), "\n", sep = "")
cat("edger_norm_factor: ", paste(ids, signif(nf, 5), sep = "=", collapse = ", "), "\n", sep = "")
cat(sprintf("regions=%d lfc.retained=%g lfc.lost=%g fdr=%g pe=%s minq=%g dedup=%s max.frag=%g\n", n_reg, lfc_ret, lfc_lost, fdr_thr,
            if (paired) "both" else "none", minq, dedup, max_frag))
pe <- if (paired) "both" else "none"
param <- csaw::readParam(minq = minq, pe = pe, max.frag = max_frag, dedup = dedup)
bp <- if (cpus > 1) BiocParallel::MulticoreParam(workers = as.integer(cpus), progressbar = FALSE) else BiocParallel::SerialParam()
counts <- csaw::regionCounts(bams, gr, ext = if (paired) NA else ext_se, param = param, BPPARAM = bp)
totals <- as.numeric(counts$totals)
names(totals) <- ids
cat("Library sizes (csaw totals): ", paste(ids, totals, sep = "=", collapse = ", "), "\n", sep = "")
if (any(!is.finite(totals)) || any(totals <= 0)) die(paste0("empty library after read filtering: ", paste(ids[!(totals > 0)], collapse = ", ")))
mat <- SummarizedExperiment::assay(counts)
if (nrow(mat) != n_reg) die("internal: regionCounts returned a different number of regions")

# ---- edgeR ----------------------------------------------------------------------------------------------------------------------
group <- factor(c(rep("ref", n_ref), rep("treat", n_treat)), levels = c("ref", "treat"))
design <- model.matrix(~ group)
has_reps <- (nrow(design) - ncol(design)) > 0
if (n_ref < 2L || n_treat < 2L) {
  warn(sprintf("%d reference and %d treated library(ies): at least one group has no replicate", n_ref, n_treat))
}
y <- edgeR::DGEList(counts = mat, lib.size = totals, norm.factors = nf, group = group)
if (has_reps) {
  if (n_ref < 2L || n_treat < 2L) warn("dispersion is estimated from the replicated group only (weakly determined)")
  y <- edgeR::estimateDisp(y, design)
  fit <- edgeR::glmQLFit(y, design)
  tt <- edgeR::glmQLFTest(fit, coef = 2)$table
  dispersion <- "estimated"
} else {
  warn("no replicates (1 reference, 1 treated): dispersion fixed at 0.05, LRT; the FDR is NOT valid")
  y$common.dispersion <- 0.05
  fit <- edgeR::glmFit(y, design)
  tt <- edgeR::glmLRT(fit, coef = 2)$table
  dispersion <- "fixed_0.05"
}
logcpm <- as.numeric(edgeR::aveLogCPM(y))
logfc <- tt$logFC
pval <- tt$PValue
fdr <- stats::p.adjust(pval, method = "BH")

# ---- classes -------------------------------------------------------------------------------------------------------------------
is_lost <- (logfc <= lfc_lost) & (fdr < fdr_thr)
is_lost[is.na(is_lost)] <- FALSE
is_ret <- (round(logfc, 1) >= lfc_ret) | (fdr >= fdr_thr)
is_ret[is.na(is_ret)] <- FALSE
n_both <- sum(is_lost & is_ret)
cls <- ifelse(is_lost, "lost", ifelse(is_ret, "retained", "unclassified"))
if (n_both > 0L) cat(n_both, " region(s) satisfy both the lost and the retained rule: classified as lost\n", sep = "")

per <- data.frame(chr = b_chr, start0 = as.integer(b_start), end = as.integer(b_end), name = b_name, logFC = logfc, logCPM = logcpm,
                  PValue = pval, FDR = fdr, class = cls)
gz <- gzfile(paste0(prefix, ".per_region.tsv.gz"), "w")
write.table(per, gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)

write_bed <- function(d, path) {
  if (nrow(d) == 0L) { file.create(path); return(invisible(0L)) }
  sc <- round(d$logFC, 4); sc[is.na(sc)] <- 0
  write.table(data.frame(d$chr, d$start0, d$end, d$name, sc, rep(".", nrow(d))), path, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  invisible(nrow(d))
}
classes <- c("retained", "lost", "unclassified")
for (k in classes) write_bed(per[per$class == k, , drop = FALSE], paste0(prefix, ".", k, ".bed"))

mb <- function(d) sum(as.numeric(d$end) - as.numeric(d$start0)) / 1e6
rows <- lapply(classes, function(k) { d <- per[per$class == k, , drop = FALSE]; data.frame(class = k, n = nrow(d), Mb = mb(d)) })
rows <- do.call(rbind, c(rows, list(data.frame(class = "total", n = nrow(per), Mb = mb(per)))))
summ <- data.frame(id = id, rows, lfc_retained = lfc_ret, lfc_lost = lfc_lost, fdr = fdr_thr, n_ref = n_ref, n_treat = n_treat,
                   ref_ids = paste(ref_ids, collapse = ","), treat_ids = paste(treat_ids, collapse = ","), norm_run = norm_run,
                   n_overlapping_regions = n_ov, dispersion = dispersion, paired_end = paired, minq = minq, dedup = dedup,
                   max_frag = max_frag)
write.table(summ, paste0(prefix, ".summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Total input regions: %d\n  retained     (round(logFC,1) >= %g or FDR >= %g): %d\n  lost         (logFC <= %g and FDR < %g): %d\n  unclassified: %d\n",
            n_reg, lfc_ret, fdr_thr, sum(cls == "retained"), lfc_lost, fdr_thr, sum(cls == "lost"), sum(cls == "unclassified")))
cat("Done: ", prefix, ".{retained,lost,unclassified}.bed, .per_region.tsv.gz, .summary.tsv\n", sep = "")
