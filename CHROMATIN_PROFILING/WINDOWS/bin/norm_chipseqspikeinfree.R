#!/usr/bin/env Rscript
# norm_chipseqspikeinfree.R -- ChIPseqSpikeInFree scale factors for ONE normalisation run.
#
# CLI contract (--key=value; called by the NORM_CHIPSEQSPIKEINFREE module):
#   --bams=<comma list>          final (filtered, deduplicated) target BAMs, coordinate-sorted + indexed
#   --sample-ids=<comma list>    sample ids, same order/length as --bams
#   --norm-run=<name>            run label (written to the package metadata as ANTIBODY = one scaling family)
#   --bin-size=<int>             ChIPseqSpikeInFree binSize (bp)
#   --cutoff=<num>               cutoff_QC     (package default 1.2)
#   --max-turns=<num>            maxLastTurn   (package default 0.99)
#   --min-first-turn=<num>       OPTIONAL minFirstTurn (package default 0.5; allowed range 0-1). When the flag is ABSENT the
#                                argument is not passed to the package at all, so the run (SF table, note, every output) is
#                                byte-identical to the previous version of this script. When given (even as 0.5) it is passed
#                                as minFirstTurn = and recorded in the `note` column ("minFirstTurn=<x>"). Legacy use:
#                                CBX2 H3K27me2 (0.7, GG_EPICYPHER_CBX2 config.py; wrapper ChIPseqSpikeInFree.R:17,274).
#   --chrom-sizes=<path>         chromosome sizes (2 columns, no header; a .fai is accepted, first 2 columns used)
#   --cpus=<int>                 workers for the read counting (the registered BiocParallel default, see below) and the
#                                package's ncores (its ParseReadCounts PSOCK cluster)
#   --out-prefix=<name>          output prefix
#   --spikein-table=<tsv>        OPTIONAL scale-factor table from norm_spikein.py (same samples): overlay of the
#                                physical spike-in factors (dashed line) on the MA plot
#   --impute-failed=TRUE|FALSE   OPTIONAL, default TRUE (see "QC failures" below)
#   --sf-out=<file>              OPTIONAL: name of the scale-factor table. Default <out-prefix>.scale_factors.tsv (the previous
#                                name and behaviour). The pipeline passes <run>.chipseqspikeinfree.scale_factors.tsv.
#   --fail-soft=TRUE|FALSE       OPTIONAL, default FALSE. TRUE (diagnostic use, CSIF is not the run's primary method): a package
#                                error or ALL libraries failing the package QC do not stop the script; the table is written
#                                with NA factors, qc_pass FALSE and a note starting "UNAVAILABLE: <reason>", exit 0 (no SVG,
#                                no package tables). FALSE (primary use): exit 3 as described below (imputation unchanged).
# Outputs:
#   <sf-out>                            (default <name>.scale_factors.tsv) sample_id, method, raw_factor, bigwig_multiplier, edger_norm_factor, qc_pass, note
#   <name>.chipseqspikeinfree.svg       MA plots (first library vs each other), 10 kb-style bins = the package's own
#                                       count table; solid red = ChIPseqSpikeInFree line, dashed red = physical spike-in
#                                       (options(preferRaster = TRUE): the density is a raster inside the SVG, so the file is
#                                       small; it was ~45 MB per run with vector density cells)
#   <name>.rawCounts.tsv                the package's bin x library count table (the table the factors were derived from)
#   <name>.parsedMatrix.tsv             the package's parsed (cumulative-fraction vs CPM cutoff) matrix
#   <name>.package_SF.tsv               the package's own SF table, unmodified (QC string, TURNS, SF)
#   <name>.chipseqspikeinfree_distribution.pdf / _boxplot.pdf   the package's own diagnostic PDFs (if produced)
# Package pin: ChIPseqSpikeInFree commit f9cc88f (do not upgrade). That commit has the minFirstTurn argument and ACCEPTS bin sizes
# outside 100-10000 bp; later package commits stop() on such bin sizes, which would break the legacy 5000/30000 bp runs.
# Exit codes: 0 success (INCLUDING runs where some samples fail the package QC: see below), 2 usage error,
#   3 data/package error (missing BAM/index, all samples failed QC, package error, unparsable output); with
#   --fail-soft=TRUE a package error or all-samples-failed QC exits 0 with the UNAVAILABLE table (other errors still exit 3).
#
# ---- What the statistic answers and what it assumes -------------------------------------------------------
# ChIPseqSpikeInFree (Jin et al. 2020) estimates, WITHOUT an exogenous spike-in, a per-library scaling factor SF from the
# shape of the genome-wide bin-count distribution (the position of the "enriched" turning point of the cumulative
# read-fraction curve). It assumes (i) the libraries of one ANTIBODY/run have comparable background and (ii) the curve has a
# detectable enrichment turn. It FAILS QC for a library with complete loss of enrichment or poor enrichment (SF = NA).
# It is therefore a *relative global-change* estimator that can be fragile (K27M-KO K27me3 round 9: only max_turns
# 0.988-0.9895 gave plausible factors). Always inspect the MA plot and <name>.chipseqspikeinfree_distribution.pdf.
#
# ---- Algebra (direction of every factor) -------------------------------------------------------------------
# Package semantics (verified on round-9 NUC data, K27M K27me3: EZH2i SF 2.19/1.85, DMSO 1.00/1.08, and step_5 used
# --scaleFactor 1/SF): SF is a DIVISOR of the depth-normalised signal:   normalised signal = RPM / SF.
#   raw_factor          = SF as written by the package (2 decimals, package precision; NA when the package QC failed)
#   SF_used             = raw_factor, or for QC failures the run's maximum finite SF when --impute-failed=TRUE
#                         (legacy round-9 rule, kept so round-9 numbers stay comparable; the sample is flagged
#                          qc_pass = FALSE and the note says "imputed"). NA when --impute-failed=FALSE.
#   bigwig_multiplier   = 1 / SF_used     -> value for bamCoverage --scaleFactor *together with* --normalizeUsing CPM|RPKM
#                         (deepTools multiplies --scaleFactor into the CPM/RPKM value: signal = RPM_or_RPKM x scaleFactor
#                          = RPM / SF).
#   edger_norm_factor   = SF_used / geometric_mean(SF_used)   (edgeR convention: norm.factors have geometric mean 1.
#                         edgeR effective library size = lib.size x norm.factors, CPM = count / (lib.size x nf) x 1e6 =
#                         RPM / nf, i.e. the same divisor semantics as SF. Centring only rescales all libraries of the
#                         run by the same constant, so no test or ratio is affected.)
#   Invariant: bigwig_multiplier x edger_norm_factor = 1 / geometric_mean(SF_used) is the same for every sample of the run.
#
# ---- QC failures -------------------------------------------------------------------------------------------
# qc_pass = FALSE for every library the package marks "fail" (SF = NA) or whose SF is not finite. The table is still written
# and the script exits 0 (a QC failure is a scientific result, not a crash); a WARNING listing the failed samples is printed to
# stderr. The script exits 3 when ALL samples fail (no information to scale on) or on any package/parsing error.
# Counting (package-internal, round-9 behaviour): every alignment (each mate) is counted once per bin by its own position
# (singleEnd = TRUE is the package default and is passed explicitly); no MAPQ / duplicate / fragment-size filtering is done
# by the package, so the BAMs must already be the final filtered set. No cache: the package is run in a fresh, empty work
# directory so an old *_rawCounts.txt / *_parsedMatrix.txt can never be reused (legacy bug).

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
write_sf_table <- function(df, path) {
  cols <- c("sample_id", "method", "raw_factor", "bigwig_multiplier", "edger_norm_factor", "qc_pass", "note")
  stopifnot(identical(colnames(df), cols))
  df$note <- gsub("[\t\r\n]+", " ", df$note)
  write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
gmean <- function(x) exp(mean(log(x)))
# --fail-soft: write the table with NA factors and exit 0. Reads ids/norm_run/sf_out/pkg_dir from the script environment.
unavailable <- function(reason) {
  cat("WARNING: ChIPseqSpikeInFree unavailable (--fail-soft): ", reason, "\n", sep = "", file = stderr())
  write_sf_table(data.frame(sample_id = ids, method = "chipseqspikeinfree", raw_factor = NA_real_,
                            bigwig_multiplier = NA_real_, edger_norm_factor = NA_real_, qc_pass = FALSE,
                            note = paste0("UNAVAILABLE: ", reason, "; run=", norm_run)), sf_out)
  unlink(pkg_dir, recursive = TRUE, force = TRUE)
  cat("Done: ", sf_out, " (UNAVAILABLE)\n", sep = "")
  quit(save = "no", status = 0L)
}

allowed <- c("bams", "sample-ids", "norm-run", "bin-size", "cutoff", "max-turns", "chrom-sizes", "cpus",
             "out-prefix", "spikein-table", "impute-failed", "min-first-turn", "sf-out", "fail-soft")
args <- parse_args(commandArgs(trailingOnly = TRUE), allowed)
bams       <- split_list(get_req(args, "bams"))
ids        <- split_list(get_req(args, "sample-ids"))
norm_run   <- get_req(args, "norm-run")
bin_size   <- as_num(get_req(args, "bin-size"), "bin-size")
cutoff     <- as_num(get_req(args, "cutoff"), "cutoff")
max_turns  <- as_num(get_req(args, "max-turns"), "max-turns")
chrom_file <- get_req(args, "chrom-sizes")
cpus       <- as_num(get_req(args, "cpus"), "cpus")
out_prefix <- get_req(args, "out-prefix")
spike_file <- get_opt(args, "spikein-table", "")
impute     <- as_bool(get_opt(args, "impute-failed", "TRUE"), "impute-failed")
sf_out     <- get_opt(args, "sf-out", paste0(out_prefix, ".scale_factors.tsv"))
fail_soft  <- as_bool(get_opt(args, "fail-soft", "FALSE"), "fail-soft")
min_first  <- if (is.null(args[["min-first-turn"]]) || !nzchar(args[["min-first-turn"]])) NULL else
  as_num(args[["min-first-turn"]], "min-first-turn")   # NULL = flag absent: package default, nothing recorded

# ---- validation -------------------------------------------------------------------------------------------------
if (length(bams) != length(ids)) die(sprintf("--bams has %d entries but --sample-ids has %d", length(bams), length(ids)))
if (length(bams) < 2L) die("ChIPseqSpikeInFree needs at least 2 libraries per run")
if (anyDuplicated(ids)) die(paste0("duplicated sample ids: ", paste(ids[duplicated(ids)], collapse = ", ")))
if (anyDuplicated(basename(bams))) die("BAM basenames must be unique (the package identifies libraries by file name)")
miss <- bams[!file.exists(bams)]
if (length(miss)) die(paste0("BAM file(s) not found: ", paste(miss, collapse = ", ")))
if (!file.exists(chrom_file)) die(paste0("chromosome sizes file not found: ", chrom_file))
if (bin_size < 1 || bin_size != round(bin_size)) die("--bin-size must be a positive integer", 2L)
if (cpus < 1 || cpus != round(cpus)) die("--cpus must be a positive integer", 2L)
if (!is.null(min_first) && (min_first < 0 || min_first > 1)) die("--min-first-turn must be between 0 and 1", 2L)
if (nzchar(spike_file) && !file.exists(spike_file)) die(paste0("--spikein-table not found: ", spike_file))
no_idx <- bams[!(file.exists(paste0(bams, ".bai")) | file.exists(sub("\\.bam$", ".bai", bams)))]
if (length(no_idx)) cat("WARNING: no .bai next to: ", paste(no_idx, collapse = ", "), "\n", sep = "", file = stderr())

suppressPackageStartupMessages({
  library(edgeR)
  library(ChIPseqSpikeInFree)
})

# The package's summarizeOverlaps uses the REGISTERED BiocParallel default (detectCores() - 2 forked workers): make it
# follow --cpus. Counting is integer arithmetic, so the factors do not depend on the worker count.
options(mc.cores = as.integer(cpus))
BiocParallel::register(if (cpus > 1) BiocParallel::MulticoreParam(workers = as.integer(cpus), progressbar = FALSE)
                       else BiocParallel::SerialParam(), default = TRUE)

# ---- fresh work directory (no cache reuse) ------------------------------------------------------------------------
pkg_dir <- paste0(out_prefix, ".pkgwork")
unlink(pkg_dir, recursive = TRUE, force = TRUE)
dir.create(pkg_dir, recursive = TRUE)
pkg_prefix <- file.path(pkg_dir, "run")

# chromosome file: keep the first two columns only (a .fai has five)
cs <- tryCatch(read.table(chrom_file, sep = "\t", header = FALSE, quote = "", comment.char = "", stringsAsFactors = FALSE),
               error = function(e) NULL)
if (is.null(cs) || ncol(cs) < 2L || !is.numeric(cs[[2]])) die("--chrom-sizes must be a tab-separated file: contig<TAB>length")
chrom_use <- file.path(pkg_dir, "chrom.sizes")
write.table(cs[, 1:2], chrom_use, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

meta <- data.frame(ID = basename(bams), ANTIBODY = gsub("[^A-Za-z0-9_.-]", "_", norm_run),
                   GROUP = gsub("[^A-Za-z0-9_.-]", "_", norm_run))
meta_file <- file.path(pkg_dir, "meta.txt")
write.table(meta, meta_file, sep = "\t", quote = FALSE, row.names = FALSE)

# ---- run the package (abort on any error) ---------------------------------------------------------------------------
cat("Running ChIPseqSpikeInFree: ", length(bams), " libraries, binSize=", bin_size, ", cutoff_QC=", cutoff,
    ", maxLastTurn=", max_turns, if (!is.null(min_first)) paste0(", minFirstTurn=", min_first), ", ncores=", cpus,
    "\n", sep = "")
if (!is.null(min_first) && min_first >= max_turns)
  cat("WARNING: --min-first-turn (", min_first, ") >= --max-turns (", max_turns, "): no admissible turning point range\n",
      sep = "", file = stderr())
# minFirstTurn is added to the call only when requested (absent flag = the previous call exactly)
pkg_extra <- if (is.null(min_first)) list() else list(minFirstTurn = min_first)
run <- tryCatch({
  do.call(ChIPseqSpikeInFree, c(list(bamFiles = bams, chromFile = chrom_use, metaFile = meta_file, prefix = pkg_prefix,
                     binSize = bin_size, cutoff_QC = cutoff, maxLastTurn = max_turns,
                     singleEnd = TRUE, disableOverwritten = FALSE, ncores = cpus), pkg_extra))
  NULL
}, error = function(e) e)
if (inherits(run, "error")) {
  if (fail_soft) unavailable(paste0("ChIPseqSpikeInFree failed: ", gsub("[\t\r\n]+", " ", conditionMessage(run))))
  die(paste0("ChIPseqSpikeInFree failed: ", conditionMessage(run)))
}

sf_path <- paste0(pkg_prefix, "_SF.txt")
raw_path <- paste0(pkg_prefix, "_rawCounts.txt")
pm_path <- paste0(pkg_prefix, "_parsedMatrix.txt")
for (p in c(sf_path, raw_path, pm_path)) if (!file.exists(p)) {
  if (fail_soft) unavailable(paste0("ChIPseqSpikeInFree did not produce ", basename(p)))
  die(paste0("ChIPseqSpikeInFree did not produce ", basename(p)))
}

sf_tab <- read.table(sf_path, sep = "\t", header = TRUE, fill = TRUE, quote = "", check.names = FALSE,
                     stringsAsFactors = FALSE)
if (!all(c("ID", "QC", "SF") %in% colnames(sf_tab))) die("package SF table lacks ID/QC/SF columns")
idx <- match(basename(bams), sf_tab$ID)
if (anyNA(idx)) die(paste0("package SF table lacks libraries: ", paste(basename(bams)[is.na(idx)], collapse = ", ")))
sf_pkg <- suppressWarnings(as.numeric(sf_tab$SF[idx]))
sf_pkg[!is.finite(sf_pkg) | sf_pkg <= 0] <- NA_real_
qc_str <- as.character(sf_tab$QC[idx])
turns  <- if ("TURNS" %in% colnames(sf_tab)) as.character(sf_tab$TURNS[idx]) else rep("", length(idx))
pkg_pass <- grepl("^pass", qc_str) & !is.na(sf_pkg)

# ---- imputation of QC failures (legacy rule, explicit) ---------------------------------------------------------------
sf_used <- sf_pkg
sf_used[!pkg_pass] <- NA_real_
fail_ids <- ids[!pkg_pass]
if (length(fail_ids) == length(ids)) {
  if (fail_soft) unavailable("all libraries failed the ChIPseqSpikeInFree QC; no scale factors can be estimated")
  die("all libraries failed the ChIPseqSpikeInFree QC; no scale factors can be estimated")
}
imputed_value <- NA_real_
if (length(fail_ids)) {
  cat("WARNING: libraries failing ChIPseqSpikeInFree QC (qc_pass=FALSE): ", paste(fail_ids, collapse = ", "), "\n",
      sep = "", file = stderr())
  if (impute) {
    imputed_value <- max(sf_used[pkg_pass])
    sf_used[!pkg_pass] <- imputed_value
  }
}
nf <- rep(NA_real_, length(sf_used))
fin <- is.finite(sf_used)
nf[fin] <- sf_used[fin] / gmean(sf_used[fin])
mft_note <- if (is.null(min_first)) "" else paste0("; minFirstTurn=", min_first)   # empty when the flag is absent
note <- ifelse(pkg_pass,
               paste0("package_QC=pass; TURNS=", turns, "; SF precision = package (2 dp)", mft_note),
               paste0("package_QC=", qc_str, "; TURNS=", turns, "; ",
                      if (impute) paste0("SF IMPUTED as run maximum finite SF (", imputed_value,
                                         "), legacy rule: no global scaling information for this library")
                      else "not imputed (--impute-failed=FALSE): factors NA, downstream must stop", mft_note))
out <- data.frame(sample_id = ids, method = "chipseqspikeinfree", raw_factor = sf_pkg,
                  bigwig_multiplier = 1 / sf_used, edger_norm_factor = nf, qc_pass = pkg_pass, note = note)
write_sf_table(out, sf_out)

# ---- copy the package's tables / PDFs ----------------------------------------------------------------------------------
invisible(file.copy(raw_path, paste0(out_prefix, ".rawCounts.tsv"), overwrite = TRUE))
invisible(file.copy(pm_path, paste0(out_prefix, ".parsedMatrix.tsv"), overwrite = TRUE))
invisible(file.copy(sf_path, paste0(out_prefix, ".package_SF.tsv"), overwrite = TRUE))
for (suf in c("distribution", "boxplot")) {
  f <- paste0(pkg_prefix, "_", suf, ".pdf")
  if (file.exists(f)) invisible(file.copy(f, paste0(out_prefix, ".chipseqspikeinfree_", suf, ".pdf"), overwrite = TRUE))
}

# ---- MA plot on the package's own count table ---------------------------------------------------------------------------
raw <- read.table(raw_path, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE, quote = "", comment.char = "")
if (!all(basename(bams) %in% colnames(raw))) die("rawCounts table lacks some libraries (column names differ from BAM basenames)")
counts <- as.matrix(raw[, basename(bams), drop = FALSE])
colnames(counts) <- ids
keep <- rowSums(counts) > 0
lcpm <- edgeR::cpm(counts[keep, , drop = FALSE], log = TRUE, prior.count = 2)

nf_phys <- NULL
if (nzchar(spike_file)) {
  sp <- read.delim(spike_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("sample_id", "edger_norm_factor") %in% colnames(sp))) die("--spikein-table lacks sample_id / edger_norm_factor columns")
  nf_phys <- sp$edger_norm_factor[match(ids, sp$sample_id)]
  if (anyNA(nf_phys) || any(nf_phys <= 0)) {
    cat("WARNING: --spikein-table does not cover every library with a positive factor; physical spike-in line omitted\n",
        file = stderr())
    nf_phys <- NULL
  }
}

options(preferRaster = TRUE)   # density cells as a raster in the SVG (file size); no effect on any factor
n <- ncol(lcpm); npan <- n - 1L
ncol_p <- min(3L, npan); nrow_p <- ceiling(npan / ncol_p)
svg(paste0(out_prefix, ".chipseqspikeinfree.svg"), width = 3.8 * ncol_p + 0.4, height = 4.3 * nrow_p + 0.4)
op <- par(mfrow = c(nrow_p, ncol_p), mar = c(5, 4, 5, 2))
for (i in seq_len(npan)) {
  x <- lcpm[, 1]; y <- lcpm[, i + 1]
  smoothScatter(x = (x + y) / 2, y = x - y,
                xlab = "Mean signal (logCPM) (A)", ylab = "log2 ratio of CPM (M)",
                main = paste0(ids[1], "\nvs\n", ids[i + 1]), cex.main = 0.8)
  abline(h = 0, col = "grey60", lty = 3)
  if (is.finite(sf_used[1]) && is.finite(sf_used[i + 1])) abline(h = log2(sf_used[1] / sf_used[i + 1]), col = "red")
  if (!is.null(nf_phys)) abline(h = log2(nf_phys[1] / nf_phys[i + 1]), col = "red", lty = 2)
  if (i == 1L) {
    legend("topleft", bty = "n", cex = 0.7, col = "red", lty = c(1, if (!is.null(nf_phys)) 2),
           legend = c("ChIPseqSpikeInFree", if (!is.null(nf_phys)) "Physical spike-in"))
  }
}
par(op)
invisible(dev.off())

unlink(pkg_dir, recursive = TRUE, force = TRUE)
cat("Done: ", sf_out, " (", sum(pkg_pass), "/", length(ids), " libraries pass package QC)\n", sep = "")
