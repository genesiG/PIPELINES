#!/usr/bin/env Rscript
# qc_fragment_sizes_table.R -- fragment (insert) size distribution of each final BAM, computed with Rsamtools.
#
# Called by QC_FRAGMENT_SIZES_TABLE (modules/local/qc_fragment_sizes_table.nf):
#   qc_fragment_sizes_table.R --bam-dir=. --samplesheet=samplesheet.csv
#
# Arguments (--key=value):
#   --bam-dir      directory with <sample_id>.final.bam (+ .bai), or <sample_id>.tn5.bam when tn5_shift is on;
#                  sample_id = file name minus ".final.bam" / ".tn5.bam"; one BAM per sample_id
#   --samplesheet  samplesheet.csv; column sample_id is required. condition, replicate, tmm_group and peak_group
#                  are copied from the sheet; a column missing from the sheet is filled with blank (metadata only)
#   --max-size     largest fragment size counted, bp (default 2000); larger fragments are counted in n_over_max
#   --nfr-max      nucleosome-free boundary, bp (default 100): frac_nfr = fragments with size < nfr-max
#   --mono-min     mononucleosome lower bound, bp (default 180)
#   --mono-max     mononucleosome upper bound, bp (default 247): frac_mono = mono-min <= size <= mono-max
#   --di-min       dinucleosome lower bound, bp (default 315)
#   --di-max       dinucleosome upper bound, bp (default 473): frac_di = di-min <= size <= di-max
#                  (boundaries as in ATACseqQC: NFR < 100, mono 180-247, di 315-473)
#   --out-prefix   default "" (outputs fragment_sizes.tsv / fragment_size_distribution.tsv)
# Method: one count per proper pair = read1 of every primary, mapped, properly paired alignment; fragment size =
#   |TLEN| (distance between the outer ends of the two mates). Read in chunks, so memory use is independent of BAM
#   size. There is NO cache between runs: the BAMs listed on disk are the only input.
#   The frac_* columns are fractions of the counted fragments (sizes 1..max-size); fragments above max-size are
#   not in the denominator (they are in n_over_max).
# Outputs:
#   fragment_sizes.tsv              sample_id, condition, replicate, tmm_group, peak_group, n_fragments, mean_size,
#                                   sd_size, median_size, mode_size, max_size, pct_lt_120 (<120 bp),
#                                   pct_120_300 (120-300 bp), pct_gt_300 (>300 bp), n_over_max, frac_nfr, frac_mono,
#                                   frac_di (fractions 0-1, 4 decimals)
#   fragment_size_distribution.tsv  sample_id, size, count, fraction (fraction of fragments <= --max-size), rows with count > 0
# Exit 3: no BAM, BAM without a samplesheet row, duplicate BAM for one sample, a BAM without any proper pair,
#         or a non-integer boundary.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(Rsamtools))
parse_args <- function(args) {
  out <- list()
  for (m in regmatches(args, regexec("^--([^=]+)=(.*)$", args))) if (length(m) == 3) out[[m[2]]] <- m[3]
  out
}
die <- function(...) { message("ERROR: ", ...); quit(save = "no", status = 3) }
opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt$samplesheet)) die("usage: qc_fragment_sizes_table.R --samplesheet=<csv> [--bam-dir=.] [--max-size=2000]")
bam_dir <- if (is.null(opt[["bam-dir"]])) "." else opt[["bam-dir"]]
int_opt <- function(key, default) if (is.null(opt[[key]])) as.integer(default) else suppressWarnings(as.integer(opt[[key]]))
cap      <- int_opt("max-size", 2000)
nfr_max  <- int_opt("nfr-max", 100)
mono_min <- int_opt("mono-min", 180)
mono_max <- int_opt("mono-max", 247)
di_min   <- int_opt("di-min", 315)
di_max   <- int_opt("di-max", 473)
if (anyNA(c(cap, nfr_max, mono_min, mono_max, di_min, di_max))) die("--max-size and the fraction boundaries must be integers")
prefix  <- if (is.null(opt[["out-prefix"]])) "" else opt[["out-prefix"]]

ss <- read.csv(opt$samplesheet, colClasses = "character", check.names = FALSE, fileEncoding = "UTF-8-BOM")
if (!("sample_id" %in% names(ss))) die("samplesheet lacks column sample_id")
meta_cols <- c("condition", "replicate", "tmm_group", "peak_group")
for (cc in meta_cols) if (!(cc %in% names(ss))) ss[[cc]] <- ""
meta <- ss[, c("sample_id", meta_cols), drop = FALSE]
meta[is.na(meta)] <- ""

bams <- sort(list.files(bam_dir, pattern = "\\.(final|tn5)\\.bam$", full.names = TRUE))
if (length(bams) == 0) die("no *.final.bam or *.tn5.bam in ", bam_dir)
ids <- sub("\\.(final|tn5)\\.bam$", "", basename(bams))
if (anyDuplicated(ids)) die("more than one BAM for sample(s): ", paste(unique(ids[duplicated(ids)]), collapse = ", "))
unknown <- setdiff(ids, ss$sample_id)
if (length(unknown)) die("BAM(s) without a samplesheet row: ", paste(unknown, collapse = ", "))

count_sizes <- function(bam) {
  param <- ScanBamParam(what = "isize", flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE, isFirstMateRead = TRUE,
                        isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE, isUnmappedQuery = FALSE))
  bf <- BamFile(bam, yieldSize = 1e6)
  open(bf); on.exit(close(bf))
  counts <- numeric(cap + 1)           # bins 1..cap, bin cap+1 = overflow
  maxsz <- 0
  repeat {
    x <- scanBam(bf, param = param)[[1]]$isize
    if (length(x) == 0) break
    x <- abs(x[!is.na(x)]); x <- x[x > 0]
    if (length(x) == 0) next
    maxsz <- max(maxsz, x)
    counts <- counts + tabulate(pmin(x, cap + 1), nbins = cap + 1)
  }
  list(counts = counts, max = maxsz)
}

summ <- list(); dist <- list()
for (i in seq_along(bams)) {
  message("fragment sizes: ", ids[i])
  r <- count_sizes(bams[i])
  cnt <- r$counts[1:cap]; over <- r$counts[cap + 1]
  n <- sum(cnt)
  if (n == 0) die("no properly paired read1 alignments in ", bams[i])
  sz <- seq_len(cap)
  mu <- sum(sz * cnt) / n
  sdv <- sqrt(sum(cnt * (sz - mu)^2) / max(n - 1, 1))
  med <- sz[which(cumsum(cnt) >= n / 2)[1]]
  summ[[i]] <- data.frame(sample_id = ids[i], n_fragments = n, mean_size = round(mu, 2), sd_size = round(sdv, 2),
    median_size = med, mode_size = sz[which.max(cnt)], max_size = r$max,
    pct_lt_120 = round(100 * sum(cnt[sz < 120]) / n, 2), pct_120_300 = round(100 * sum(cnt[sz >= 120 & sz <= 300]) / n, 2),
    pct_gt_300 = round(100 * sum(cnt[sz > 300]) / n, 2), n_over_max = over,
    frac_nfr = round(sum(cnt[sz < nfr_max]) / n, 4),
    frac_mono = round(sum(cnt[sz >= mono_min & sz <= mono_max]) / n, 4),
    frac_di = round(sum(cnt[sz >= di_min & sz <= di_max]) / n, 4))
  keep <- cnt > 0
  dist[[i]] <- data.frame(sample_id = ids[i], size = sz[keep], count = cnt[keep], fraction = signif(cnt[keep] / n, 6))
}
summ <- merge(meta, do.call(rbind, summ), by = "sample_id")
summ <- summ[match(ids, summ$sample_id), ]
write.table(summ, paste0(prefix, "fragment_sizes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, dist), paste0(prefix, "fragment_size_distribution.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote fragment_sizes.tsv and fragment_size_distribution.tsv for", length(ids), "samples\n")
