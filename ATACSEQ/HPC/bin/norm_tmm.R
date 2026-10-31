#!/usr/bin/env Rscript
# norm_tmm.R -- TMM (composition-bias) scale factors on large genomic bins for ONE ATAC-seq normalisation group.
#
# CLI contract (called by the NORM_TMM module, one task per TMM group):
#   --bams=<comma list>          final BAMs of the group (coordinate-sorted + indexed, paired-end)
#   --sample-ids=<comma list>    same order/length as --bams
#   --norm-run=<name>            group label (note column; MA plot title)
#   --bin-size=<int>             bin width in bp (config tmm_bin_size, default 10000; must be >= 1000)
#   --mapq=<int>                 csaw readParam(minq=)
#   --dedup=TRUE|FALSE           csaw readParam(dedup=): TRUE drops reads carrying the duplicate flag (0x400). Pipeline default
#                                FALSE (BAMs are duplicate-REMOVED). Use TRUE only for markdup-flagged, not-removed BAMs.
#   --max-frag=auto|<int>        fragments longer than this are dropped by csaw (PE). auto = the maximum observed proper-pair
#                                fragment size over ALL BAMs of the group (see "max.frag" below); an integer is used as given.
#   --pe=both|none               OPTIONAL, default both (both datasets are paired-end; none = single-end counting, then
#                                --max-frag must be an integer and is ignored by csaw)
#   --bigwig-units=none|cpm|rpkm how bigwig_multiplier will be used in bamCoverage (see Algebra)
#   --out=<GROUP>.scale_factors.tsv
# Outputs (prefix P = --out without ".scale_factors.tsv"):
#   P.scale_factors.tsv  (sample_id, method, raw_factor, bigwig_multiplier, edger_norm_factor, qc_pass, note;
#                         schema /project/rephilllab/CUTANDRUN/docs/scale_factor_schema.md)
#   P.ma_plot.png, P.ma_plot.pdf   MA plots of library 1 vs every other library, unnormalised logCPM, red line =
#                         log2(nf_1/nf_j) (legacy scaleFactors.R plot; PDF = all comparisons, 6 per page; PNG = first 12)
#   P.max_frag.txt       the max.frag value used (one integer; for --max-frag=auto the computed value)
# Exit codes: 0 ok; 2 usage error; 3 data problem (missing BAM/index, <2 libraries, no proper pairs, non-finite factors).
#
# ---- Question, assumptions, and where they fail -------------------------------------------------------------------
# Question: by what factor must each library's size be adjusted so that the BULK of the genome (10 kb bins) has the same
# accessibility signal in every library of the group? This removes composition/depth bias (a library with a larger share of
# its fragments in a few strongly open regions, e.g. different signal-to-noise or FRiP).
# Assumption of TMM (Robinson & Oshlack 2010): MOST bins are NOT differentially accessible between any two libraries of the
# group (a 30 percent per side trimmed mean of the bin log-ratios estimates the factor). Where it FAILS: libraries of one
# group that differ globally in accessibility (a drug/knock-out/differentiation that opens or closes most of the genome, or
# sorted fractions of different global state). TMM then aligns the unchanged remainder and ERASES the global difference.
# Consequence: the factors are legitimate for display-scaling of replicates and for groups whose members are expected to
# share the bulk of the genome; they must NOT be used to claim or deny a GLOBAL change in accessibility between members of a
# group. Group composition (config tmm_group_column) therefore matters: the exact groups keep K27M and K27MKO separate, the
# legacy-membership column (tmm_group_legacy) pools them and so also forces their bulk to be equal.
# Violations are visible in the MA plot: the red line must pass through the centre of the dense low-A background cloud; if
# the cloud is clearly displaced from the line, or is not horizontal, the assumption is violated for that pair.
#
# ---- Implementation (csaw 1.36.0 / edgeR 4.0.16 / Rsamtools 2.18.0; R 4.3.3) ------------------------------------------
# counts   = csaw::windowCounts(bin = TRUE, width = bin_size, filter = 10, param = readParam(minq, pe, max.frag, dedup))
#            filter = 10 (the csaw default, as legacy): bins with < 10 fragments summed over ALL libraries are removed BEFORE
#            normalisation; irrelevant for 10 kb bins at ATAC depth. All contigs present in the BAM header are counted, so
#            chrM must be removed upstream (pipeline remove_chroms) or it forms one extreme bin (TMM trims it, but it still
#            enters lib.size).
# nf       = csaw::normFactors(counts, weighted = FALSE, se.out = FALSE)
#            weighted = FALSE: plain trimmed mean of M-values (composition bias; legacy EFFICIENCY_BIAS = FALSE).
#            csaw's library size is counts$totals (total fragments = proper pairs counted in the BAM that pass minq/max.frag/
#            dedup, NOT the column sums of the retained bins) -- this is the lib.size legacy used and edgeR uses.
#            nf has geometric mean 1 (edgeR convention); a deviation > 1e-6 is treated as a data problem (exit 3).
#
# ---- max.frag --------------------------------------------------------------------------------------------------------
# Legacy getPEsizes.R took max(csaw::getPESizes(bam)$sizes) per sample but stored `max.frag[i] = max.frag %>% append(...)`
# (line 78), which keeps only the first sample's maximum at every index; maxFrag was therefore the first sample's maximum
# (legacy value 1999 with bowtie2 -X 2000). Here auto = the true maximum over all BAMs of the group.
# Efficient computation: getPESizes loads every pair into memory (tens of millions for a deep ATAC BAM); instead each BAM is
# streamed with Rsamtools::scanBam(what = flag, mapq, isize) in chunks of 2e6 records, keeping only a running maximum, so
# memory is O(chunk). A record counts when: paired, proper-pair flag (0x2), mapped, primary (no 0x100/0x800), MAPQ >= minq
# (MAPQ 255 = unavailable passes, as in samtools -q), not duplicate-flagged when --dedup=TRUE; fragment size = |isize| (TLEN),
# both mates carry the same |TLEN| so the maximum is unaffected by counting both. csaw's own pair test is orientation based
# (opposite strands, forward mate upstream) rather than the 0x2 flag; for the pipeline's final BAMs (samtools -f 3, bowtie2
# --no-discordant --no-mixed) the two coincide. Because auto is the largest observed size, max.frag then drops nothing
# (legacy behaviour); pass an integer (e.g. the bowtie2 -X value or ~1000) to actually trim long fragments.
#
# ---- Algebra ----------------------------------------------------------------------------------------------------------
#   edgeR effective library size = lib.size x nf;  normalised signal = count / (lib.size x nf) = RPM-type / nf
#   raw_factor         = nf                                   (TMM norm.factor; the legacy table stored nf x lib.size / 1e6)
#   edger_norm_factor  = nf                                   (geometric mean 1; use with DGEList(lib.size = counts$totals))
#   legacy_SF          = nf x lib.size / 1e6                  (legacy tmm_scaling_factors_*.txt column SF; stored in `note`)
#   bigwig_multiplier  = depends on --bigwig-units:
#     none     1e6 / (nf x lib.size) = 1 / legacy_SF  -> bamCoverage --scaleFactor <m> WITHOUT --normalizeUsing (exact legacy:
#              step_5 passed 1/SF; the track is then reads-per-bin x 1e6 / (fragments x nf)).
#     cpm|rpkm 1 / nf                                  -> bamCoverage --scaleFactor <m> WITH --normalizeUsing CPM|RPKM (deepTools
#              multiplies --scaleFactor into the CPM/RPKM value = count / total_deepTools x 1e6 x (1/nf); equals the "none"
#              track up to the per-BAM constant csaw-total / deepTools-total, ~1/2 for PE because deepTools counts mates).
#   Direction: all methods express the factor as a DIVISOR of the depth-normalised signal, signal = CPM / nf.
#   Invariants: cpm|rpkm: bigwig_multiplier x edger_norm_factor = 1;  none: bigwig_multiplier x nf x lib.size / 1e6 = 1.
#   MA reference line between libraries i and j: log2(nf_i / nf_j) (M = logCPM_i - logCPM_j).
# qc_pass = FALSE when a library has < 100000 fragments (factor imprecise) or fewer than 1000 bins are retained.
#
# ---- Not done here -----------------------------------------------------------------------------------------------------
# The MA plot uses edgeR::cpm(log = TRUE, prior.count = 2) with lib.size = counts$totals and no normalisation factors, as
# legacy; the legacy "+ 6 log2(10)" shift of the A axis (cosmetic) is dropped. No blacklist/greylist filtering of bins.

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

# Maximum proper-pair fragment size |TLEN| of one BAM, streamed in chunks (memory O(chunk)); see header "max.frag".
max_proper_pair_size <- function(bam, minq, dedup, chunk = 2e6) {
  bf <- Rsamtools::BamFile(bam, yieldSize = chunk)
  open(bf)
  on.exit(close(bf))
  sp <- Rsamtools::ScanBamParam(what = c("flag", "mapq", "isize"))
  best <- 0L
  n_used <- 0
  n_seen <- 0
  repeat {
    x <- Rsamtools::scanBam(bf, param = sp)[[1]]
    n <- length(x$flag)
    if (n == 0L) break
    n_seen <- n_seen + n
    fl <- x$flag
    keep <- bitwAnd(fl, 0x1L) != 0L & bitwAnd(fl, 0x2L) != 0L & bitwAnd(fl, 0x4L) == 0L & bitwAnd(fl, 0x900L) == 0L
    if (dedup) keep <- keep & bitwAnd(fl, 0x400L) == 0L
    if (minq > 0) keep <- keep & (is.na(x$mapq) | x$mapq >= minq)
    sz <- abs(x$isize[keep])
    sz <- sz[!is.na(sz) & sz > 0L]
    if (length(sz)) {
      best <- max(best, max(sz))
      n_used <- n_used + length(sz)
    }
  }
  list(max = as.integer(best), n_records = n_used, n_seen = n_seen)
}

draw_ma <- function(lcpm, ids, nf, title, which_j) {
  for (j in which_j) {
    smoothScatter(x = (lcpm[, 1] + lcpm[, j]) / 2, y = lcpm[, 1] - lcpm[, j],
                  xlab = "Average signal (logCPM)", ylab = "Difference in signal (log2 ratio, unnormalised)",
                  main = paste0(title, "\n", ids[1], " vs ", ids[j]), cex.main = 0.8)
    abline(h = 0, col = "grey50", lty = 2)
    abline(h = log2(nf[1] / nf[j]), col = "red", lwd = 2)  # log2 ratio of norm factors: should cross the background cloud
  }
}

allowed <- c("bams", "sample-ids", "norm-run", "bin-size", "mapq", "dedup", "out", "max-frag", "pe", "bigwig-units")
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) == 0L || any(argv %in% c("-h", "--help"))) {
  cat("usage: norm_tmm.R --bams=a.bam,b.bam --sample-ids=x,y --norm-run=GROUP --bin-size=10000 --mapq=30 --dedup=FALSE",
      "--max-frag=auto|INT --pe=both --bigwig-units=none|cpm|rpkm --out=GROUP.scale_factors.tsv\n", file = stderr())
  quit(save = "no", status = 2L)
}
args <- parse_args(argv, allowed)
bams      <- split_list(get_req(args, "bams"))
ids       <- split_list(get_req(args, "sample-ids"))
norm_run  <- get_req(args, "norm-run")
bin_size  <- as_num(get_req(args, "bin-size"), "bin-size")
mapq      <- as_num(get_req(args, "mapq"), "mapq")
dedup     <- as_bool(get_req(args, "dedup"), "dedup")
out_file  <- get_req(args, "out")
mf_arg    <- get_opt(args, "max-frag", "auto")
pe        <- get_opt(args, "pe", "both")
units     <- get_req(args, "bigwig-units")
if (!(pe %in% c("both", "none"))) die("--pe must be 'both' or 'none'", 2L)
if (!(units %in% c("none", "cpm", "rpkm"))) die("--bigwig-units must be none, cpm or rpkm", 2L)
if (bin_size < 1000 || bin_size != round(bin_size)) die("--bin-size must be an integer >= 1000 bp (composition-bias normalisation needs large bins; >= 10 kb advised)", 2L)
if (mapq < 0 || mapq != round(mapq)) die("--mapq must be a non-negative integer", 2L)
mf_auto <- tolower(trimws(mf_arg)) == "auto"
if (mf_auto) {
  if (pe != "both") die("--max-frag=auto needs --pe=both; give an integer", 2L)
  max_frag <- NA_real_
} else {
  max_frag <- as_num(mf_arg, "max-frag")
  if (max_frag < 1 || max_frag != round(max_frag)) die("--max-frag must be 'auto' or a positive integer", 2L)
}
prefix <- sub("\\.scale_factors\\.tsv$", "", out_file)
if (identical(prefix, out_file)) die("--out must end in .scale_factors.tsv", 2L)

if (length(bams) != length(ids)) die(sprintf("--bams has %d entries but --sample-ids has %d", length(bams), length(ids)))
if (length(bams) < 2L) die("TMM needs at least 2 libraries per group")
if (anyDuplicated(ids)) die(paste0("duplicated sample ids: ", paste(ids[duplicated(ids)], collapse = ", ")))
if (anyDuplicated(bams)) die("duplicated BAM paths in --bams")
check_bams(bams)

suppressPackageStartupMessages({
  library(csaw)
  library(edgeR)
})

# ---- max.frag -------------------------------------------------------------------------------------------------------
sample_max <- rep(NA_integer_, length(bams))
if (mf_auto) {
  for (i in seq_along(bams)) {
    r <- tryCatch(max_proper_pair_size(bams[i], mapq, dedup), error = function(e) die(paste0("cannot read ", bams[i], ": ", conditionMessage(e))))
    if (r$max < 1L) die(paste0("no proper-paired fragments found in ", bams[i], " (is it a paired-end, filtered final BAM?)"))
    sample_max[i] <- r$max
    cat("max proper-pair fragment size ", ids[i], ": ", r$max, " (", r$n_records, " of ", r$n_seen, " records used)\n", sep = "")
  }
  max_frag <- max(sample_max)   # true maximum over ALL libraries (legacy getPEsizes.R:78 kept only the first sample's)
}
max_frag <- as.integer(max_frag)
mf_src <- if (mf_auto) "auto" else "given"
writeLines(as.character(max_frag), paste0(prefix, ".max_frag.txt"))

# ---- counting and TMM ------------------------------------------------------------------------------------------------
param <- readParam(minq = mapq, pe = pe, max.frag = max_frag, dedup = dedup)
cat("Counting ", bin_size, " bp bins in ", length(bams), " libraries (minq=", mapq, ", pe=", pe, ", max.frag=", max_frag,
    " [", mf_src, "], dedup=", dedup, ")\n", sep = "")
counts <- tryCatch(windowCounts(bams, bin = TRUE, width = as.integer(bin_size), filter = 10, param = param),
                   error = function(e) die(paste0("windowCounts failed: ", conditionMessage(e))))
colnames(counts) <- ids
totals <- as.numeric(counts$totals)
n_inform <- nrow(counts)
cat("Library sizes (csaw totals): ", paste(ids, totals, sep = "=", collapse = ", "), "\n", sep = "")
cat("Bins retained: ", n_inform, "\n", sep = "")
if (any(!is.finite(totals)) || any(totals < 1)) die("a library has no counted fragments (csaw totals < 1)")

nf <- as.numeric(csaw::normFactors(counts, weighted = FALSE, se.out = FALSE))
if (length(nf) != length(ids) || any(!is.finite(nf)) || any(nf <= 0)) {
  die(paste0("TMM produced non-finite or non-positive factors: ", paste(signif(nf, 4), collapse = ", ")))
}
gm <- exp(mean(log(nf)))
if (abs(gm - 1) > 1e-6) die(sprintf("TMM factors are not geometric-mean 1 (%.8g); refusing to guess a centring", gm))

if (units == "none") {
  bigwig_multiplier <- 1e6 / (nf * totals)   # = 1 / legacy_SF; bamCoverage WITHOUT --normalizeUsing
} else {
  bigwig_multiplier <- 1 / nf                # bamCoverage WITH --normalizeUsing CPM|RPKM
}
qc_pass <- totals >= 1e5 & n_inform >= 1e3
if (!all(qc_pass)) cat("WARNING: qc_pass=FALSE for: ", paste(ids[!qc_pass], collapse = ", "), "\n", sep = "", file = stderr())

if (units == "none") {
  unit_note <- "bigwig_multiplier=1e6/(nf*lib_size) for bamCoverage --scaleFactor WITHOUT --normalizeUsing (legacy)"
} else {
  unit_note <- sprintf("bigwig_multiplier=1/nf for bamCoverage --scaleFactor WITH --normalizeUsing %s", toupper(units))
}
out <- data.frame(
  sample_id = ids, method = "tmm", raw_factor = nf, bigwig_multiplier = bigwig_multiplier, edger_norm_factor = nf,
  qc_pass = qc_pass,
  note = sprintf("run=%s; TMM weighted=FALSE on %d bins of %d bp; lib_size=%d; legacy_SF=%.8g; max.frag=%d (%s%s); mapq=%d; dedup=%s; bigwig_units=%s (%s); assumes most bins not differentially accessible within the group",
                 norm_run, n_inform, as.integer(bin_size), as.integer(totals), nf * totals / 1e6, max_frag, mf_src,
                 ifelse(is.na(sample_max), "", paste0("; sample max ", sample_max)), as.integer(mapq), dedup, units, unit_note)
)
write_sf_table(out, out_file)
cat("Wrote ", out_file, "\n", sep = "")

# ---- MA plots ---------------------------------------------------------------------------------------------------------
y <- DGEList(counts = SummarizedExperiment::assay(counts), lib.size = totals)
lcpm <- cpm(y, log = TRUE, prior.count = 2)       # unnormalised (norm.factors = 1), as legacy
ttl <- sprintf("%s (%d bins, %d bp)", norm_run, n_inform, as.integer(bin_size))
js <- seq_along(ids)[-1]
tryCatch({
  pdf(paste0(prefix, ".ma_plot.pdf"), width = 11, height = 7.5)
  par(mfrow = c(2, 3), mar = c(5, 4, 5, 2))
  draw_ma(lcpm, ids, nf, ttl, js)
  dev.off()
  js_png <- head(js, 12L)
  ncol_p <- min(3L, length(js_png)); nrow_p <- ceiling(length(js_png) / ncol_p)
  ptype <- if (isTRUE(capabilities("cairo"))) "cairo" else getOption("bitmapType")
  png(paste0(prefix, ".ma_plot.png"), width = 4.2 * ncol_p, height = 4.2 * nrow_p, units = "in", res = 150, type = ptype)
  par(mfrow = c(nrow_p, ncol_p), mar = c(5, 4, 5, 2))
  draw_ma(lcpm, ids, nf, ttl, js_png)
  dev.off()
}, error = function(e) die(paste0("MA plot failed: ", conditionMessage(e))))
cat("Done\n")
