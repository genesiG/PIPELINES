#!/usr/bin/env Rscript
# norm_tmm.R -- TMM (composition-bias) scale factors on large genomic bins for ONE normalisation run.
#
# CLI contract (called by the NORM_TMM module):
#   Two input modes (exactly one of --bams / --counts):
#   BAM mode (original; kept for CLI compatibility):
#   --bams=<comma list>        final target BAMs (coordinate-sorted + indexed)
#   --sample-ids=<comma list>  same order/length as --bams
#   Counts mode (the pipeline: counts come from bin_counts.R, filter 0, all active libraries of the project):
#   --counts=<tsv.gz>          bin_counts.R <prefix>.counts.tsv.gz (chr, start, end, blacklisted, one column per sample id)
#   --libsizes=<tsv>           bin_counts.R <prefix>.libsizes.tsv; T_i = csaw_total of each id of the run. --sample-ids selects
#                              the run's columns (verbatim ids; the file may hold more). The counts' bin size must equal --bin-size.
#                              Bins with rowSums over the run's libraries < 10 are dropped (= windowCounts(filter = 10) on the
#                              run alone); the note's max.frag/dedup come from the libsizes file. Factors are numerically
#                              identical to BAM mode on the same BAMs and parameters (bin_counts.R uses the same readParam).
#                              --mapq, --dedup, --max-frag and --pe are not needed (ignored with a WARNING if given); --bams is not
#                              allowed. The factors use T_i = csaw_total as lib.size (NOT the column sums of the retained bins), the
#                              same csaw::normFactors call as BAM mode on a SummarizedExperiment with colData totals = T.
#   --sample-ids=<comma list>  (counts mode) the run's ids
#   --norm-run=<name>          run label (note column)
#   --bin-size=<int>           bin width in bp (config tmm_bin_size, default 10000)
#   --mapq=<int>               csaw readParam(minq=)
#   --dedup=TRUE|FALSE         csaw readParam(dedup=): TRUE drops reads carrying the duplicate flag (0x400). Use TRUE only when
#                              the BAMs are markdup-flagged and NOT already duplicate-removed; with duplicate-removed BAMs
#                              (pipeline default) TRUE and FALSE give identical counts. Caution: counting flagged duplicates
#                              in targets but not in IgG (dedup=FALSE on flagged BAMs) biases target-vs-IgG comparisons.
#   --out=<name>.scale_factors.tsv
#   --fail-soft=TRUE|FALSE     OPTIONAL, default FALSE. TRUE (diagnostic use): non-finite or non-positive TMM factors do not
#                              stop the script; the table is written with NA factors, qc_pass FALSE and note
#                              "UNAVAILABLE: <reason>", exit 0. FALSE (primary use): exit 3 as before.
#   --max-frag=<int>           OPTIONAL, default 500 (csaw default): fragments longer than this are dropped (PE). Pass the
#                              same value as the peak caller (bowtie2_maxins) (500 if unset). Legacy used the maximum observed size (no cap).
#   --pe=both|none             OPTIONAL, default both
# Output: <name>.scale_factors.tsv  (sample_id, method, raw_factor, bigwig_multiplier, edger_norm_factor, qc_pass, note)
# Exit codes: 0 ok; 2 usage error; 3 data problem (missing BAM/index, <2 libraries, non-finite factors unless --fail-soft=TRUE).
#
# ---- Question, assumptions, and where they fail -------------------------------------------------------------------
# Question: by what factor must each library's size be adjusted so that the BULK of the genome (large bins) has the same
# abundance in every library? This removes composition/depth bias (e.g. a library with a larger share of its reads in a few
# enriched regions). Assumption of TMM (Robinson & Oshlack 2010): MOST bins are not differentially abundant between any two
# libraries (a trimmed mean of the M-values, 30 percent per side, estimates the factor). VIOLATED by a global loss/gain of
# the mark: after an EZH2/EED inhibitor most H3K27me3 domain bins drop, and TMM then aligns the unaffected remainder and
# thereby ERASES the global change. Consequence: TMM factors are suitable for comparisons of libraries of the same
# condition (replicates; target vs IgG with background-bin normalisation, see callpeaks_csaw.R) and must NOT be used to claim
# or deny a global change between treated and untreated cells -- use norm_spikein.py or norm_chipseqspikeinfree.R for that.
#
# ---- Implementation (csaw/edgeR of the run's image; cluster: 1.36.0/4.0.16, Docker omics image: 1.40.0/4.4.0) --------
# counts   = csaw::windowCounts(bin = TRUE, width = bin_size, filter = 10, param = readParam(minq, pe, max.frag, dedup))
#            (filter = 10 = the csaw default: bins with < 10 fragments summed over all libraries are dropped BEFORE
#             normalisation, as in the legacy script; irrelevant for 10 kb bins at normal depth)
# nf       = csaw::normFactors(counts, weighted = FALSE, se.out = FALSE)
#            `weighted` IS a documented argument of csaw normFactors (checked in 1.36.0 and 1.40.0) (default FALSE; forwarded to
#            edgeR::calcNormFactors(doWeighting = weighted)); the legacy call `weighted = EFFICIENCY_BIAS` (= FALSE) was
#            therefore valid. weighted = FALSE keeps the plain TMM trimmed mean, as the legacy comment intended for composition
#            bias. csaw uses the library totals (counts$totals), NOT the column sums of the retained bins, as lib.size.
# nf has geometric mean 1 (edgeR convention).
#
# ---- Algebra ----------------------------------------------------------------------------------------------------------
#   edgeR effective library size = lib.size x nf;  normalised signal = count / (lib.size x nf) = RPM-type / nf
#   raw_factor         = nf                                   (TMM norm.factor; the legacy table stored nf x lib.size / 1e6)
#   edger_norm_factor  = nf                                   (already geometric-mean-centred by edgeR)
#   bigwig_multiplier  = 1 / nf                               (bamCoverage --scaleFactor with --normalizeUsing CPM|RPKM:
#                        deepTools multiplies --scaleFactor into the CPM/RPKM value, i.e. count/lib x 1e6 x (1/nf) =
#                        count / (lib x nf) x 1e6.) Legacy: --scaleFactor 1/(nf x lib/1e6) without --normalizeUsing; the two
#                        agree up to the constant deepTools-total / csaw-total of each BAM.
#   Legacy comparison value is stored in the note: legacy_SF = nf x lib.size / 1e6.
#   Invariant: bigwig_multiplier x edger_norm_factor = 1.
# qc_pass = FALSE when a library has < 100000 fragments (factor imprecise) or fewer than 1000 bins are retained.

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
write_sf_table <- function(df, path) {
  cols <- c("sample_id", "method", "raw_factor", "bigwig_multiplier", "edger_norm_factor", "qc_pass", "note")
  stopifnot(identical(colnames(df), cols))
  df$note <- gsub("[\t\r\n]+", " ", df$note)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

allowed <- c("bams", "sample-ids", "norm-run", "bin-size", "mapq", "dedup", "out", "max-frag", "pe",
             "counts", "libsizes", "fail-soft")
args <- parse_args(commandArgs(trailingOnly = TRUE), allowed)
counts_mode <- !is.null(args[["counts"]]) && nzchar(args[["counts"]])
ids      <- split_list(get_req(args, "sample-ids"))
norm_run <- get_req(args, "norm-run")
bin_size <- as_num(get_req(args, "bin-size"), "bin-size")
out_file <- get_req(args, "out")
fail_soft <- as_bool(get_opt(args, "fail-soft", "FALSE"), "fail-soft")
pe       <- get_opt(args, "pe", "both")
if (counts_mode) {
  if (!is.null(args[["bams"]]) && nzchar(args[["bams"]])) die("--counts and --bams are mutually exclusive", 2L)
  ignored <- intersect(c("mapq", "dedup", "max-frag", "pe"), names(args))
  if (length(ignored)) cat("WARNING: counts mode ignores ", paste0("--", ignored, collapse = ", "),
                           " (mapq, max.frag and dedup come from the libsizes file)
", sep = "", file = stderr())
  counts_file  <- get_req(args, "counts")
  libsize_file <- get_req(args, "libsizes")
  bams <- character(0)
} else {
  bams     <- split_list(get_req(args, "bams"))
  mapq     <- as_num(get_req(args, "mapq"), "mapq")
  dedup    <- as_bool(get_req(args, "dedup"), "dedup")
  max_frag <- as_num(get_opt(args, "max-frag", "500"), "max-frag")
}
if (!(pe %in% c("both", "none"))) die("--pe must be 'both' or 'none'", 2L)
if (bin_size < 1000) die("--bin-size < 1000 bp is not composition-bias normalisation (use >= 10 kb)", 2L)

if (!counts_mode && length(bams) != length(ids)) die(sprintf("--bams has %d entries but --sample-ids has %d", length(bams), length(ids)))
if (length(ids) < 2L) die("TMM needs at least 2 libraries per run")
if (anyDuplicated(ids)) die(paste0("duplicated sample ids: ", paste(ids[duplicated(ids)], collapse = ", ")))
if (counts_mode) {
  if (!file.exists(counts_file)) die(paste0("counts file not found: ", counts_file))
  if (!file.exists(libsize_file)) die(paste0("libsizes file not found: ", libsize_file))
} else check_bams(bams)

suppressPackageStartupMessages({
  library(csaw)
  library(edgeR)
})

if (counts_mode) {
  ct <- read.delim(counts_file, check.names = FALSE, colClasses = c(chr = "character"))
  ls <- read.delim(libsize_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("sample_id", "csaw_total", "bin_size", "max_frag", "dedup") %in% colnames(ls))) die("--libsizes lacks sample_id/csaw_total/bin_size/max_frag/dedup columns")
  miss <- setdiff(ids, colnames(ct))
  if (length(miss)) die(paste0("sample id(s) not in the counts table: ", paste(miss, collapse = ", ")))
  miss <- setdiff(ids, ls$sample_id)
  if (length(miss)) die(paste0("sample id(s) not in the libsizes table: ", paste(miss, collapse = ", ")))
  li <- match(ids, ls$sample_id)
  if (any(as.numeric(ls$bin_size[li]) != bin_size)) die(sprintf("the counts were made with bin size %s, --bin-size is %s", paste(unique(ls$bin_size[li]), collapse = ","), format(bin_size)))
  max_frag <- as.numeric(ls$max_frag[li][1])
  dedup <- as_bool(as.character(ls$dedup[li][1]), "libsizes dedup")
  totals <- as.numeric(ls$csaw_total[li])
  cm <- as.matrix(ct[, ids, drop = FALSE])
  keep <- rowSums(cm) >= 10
  cm <- cm[keep, , drop = FALSE]
  n_inform <- nrow(cm)
  cat("Counts mode: ", ncol(cm), " libraries from ", counts_file, "\n", sep = "")
  cat("Library sizes (csaw totals): ", paste(ids, totals, sep = "=", collapse = ", "), "\n", sep = "")
  cat("Bins retained: ", n_inform, "\n", sep = "")
} else {
  param <- readParam(minq = mapq, pe = pe, max.frag = max_frag, dedup = dedup)
  cat("Counting ", bin_size, " bp bins in ", length(bams), " libraries (minq=", mapq, ", pe=", pe, ", max.frag=", max_frag,
      ", dedup=", dedup, ")\n", sep = "")
  counts <- windowCounts(bams, bin = TRUE, width = bin_size, filter = 10, param = param)
  colnames(counts) <- ids
  totals <- counts$totals
  n_inform <- nrow(counts)
  cat("Library sizes (csaw totals): ", paste(ids, totals, sep = "=", collapse = ", "), "\n", sep = "")
  cat("Bins retained: ", nrow(counts), "\n", sep = "")
}

unavailable <- function(reason) {
  cat("WARNING: TMM unavailable (--fail-soft): ", reason, "\n", sep = "", file = stderr())
  write_sf_table(data.frame(sample_id = ids, method = "tmm", raw_factor = NA_real_, bigwig_multiplier = NA_real_,
                            edger_norm_factor = NA_real_, qc_pass = FALSE,
                            note = paste0("UNAVAILABLE: ", reason, "; run=", norm_run)), out_file)
  cat("Done: ", out_file, " (UNAVAILABLE)\n", sep = "")
  quit(save = "no", status = 0L)
}

nf <- if (counts_mode) {
  tryCatch({
    # same code path as BAM mode: csaw reads the library sizes from colData(se)$totals (a plain matrix would ignore lib.sizes)
    se <- SummarizedExperiment::SummarizedExperiment(list(counts = cm), colData = S4Vectors::DataFrame(totals = totals))
    csaw::normFactors(se, weighted = FALSE, se.out = FALSE)
  }, error = function(e) e)
} else csaw::normFactors(counts, weighted = FALSE, se.out = FALSE)
if (inherits(nf, "error")) {
  if (fail_soft) unavailable(paste0("TMM failed: ", conditionMessage(nf)))
  die(paste0("TMM failed: ", conditionMessage(nf)))
}
nf <- as.numeric(nf)
if (length(nf) != length(ids) || any(!is.finite(nf)) || any(nf <= 0)) {
  msg <- paste0("TMM produced non-finite or non-positive factors: ", paste(signif(nf, 4), collapse = ", "))
  if (fail_soft) unavailable(msg)
  die(msg)
}
if (abs(exp(mean(log(nf))) - 1) > 1e-6) {
  cat("NOTE: TMM factors are not geometric-mean 1 (", exp(mean(log(nf))), "); re-centring edger_norm_factor\n", sep = "", file = stderr())
}
nf_edger <- nf / exp(mean(log(nf)))

qc_pass <- totals >= 1e5 & n_inform >= 1e3
if (!all(qc_pass)) {
  cat("WARNING: qc_pass=FALSE for: ", paste(ids[!qc_pass], collapse = ", "), "\n", sep = "", file = stderr())
}
out <- data.frame(
  sample_id = ids, method = "tmm", raw_factor = nf, bigwig_multiplier = 1 / nf, edger_norm_factor = nf_edger,
  qc_pass = qc_pass,
  note = sprintf("run=%s; TMM weighted=FALSE on %d bins of %d bp; lib_size=%d; legacy_SF=%.6g; max.frag=%d; dedup=%s; assumes most bins unchanged (not valid for global loss)",
                 norm_run, n_inform, as.integer(bin_size), as.integer(totals), nf * totals / 1e6,
                 as.integer(max_frag), dedup)
)
write_sf_table(out, out_file)
cat("Done: ", out_file, "\n", sep = "")
