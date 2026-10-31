#!/usr/bin/env Rscript
# qc_fragment_sizes_table.R -- fragment (insert) size distribution of each final BAM, computed with Rsamtools.
#
# Called by QC_FRAGMENT_SIZES_TABLE (modules/local/qc_fragment_sizes_table.nf):
#   qc_fragment_sizes_table.R --bam-dir=. --samplesheet=samplesheet.csv --threads=<n>
#
# Arguments (--key=value):
#   --bam-dir      directory with <sample_id>.final.bam (+ .bai); sample_id = file name minus ".final.bam"
#   --samplesheet  samplesheet.csv; columns sample_id, antibody, condition, replicate, group (metadata only from here)
#   --max-size     largest fragment size counted, bp (default 2000); larger fragments are counted in n_over_max
#   --threads      BAMs counted in parallel, one forked worker per BAM (parallel::mclapply; default 1)
#   --out-prefix   default "" (outputs fragment_sizes.tsv / fragment_size_distribution.tsv)
# Method: one count per proper pair = read1 of every primary, mapped, properly paired alignment; fragment size =
#   |TLEN| (distance between the outer ends of the two mates). Read in chunks, so memory use is independent of BAM
#   size. There is NO cache between runs: the BAMs listed on disk are the only input.
# Outputs:
#   fragment_sizes.tsv              sample_id, antibody, condition, replicate, group, n_fragments, mean_size, sd_size,
#                                   median_size, mode_size, max_size, pct_lt_120 (<120 bp), pct_120_300 (120-300 bp),
#                                   pct_gt_300 (>300 bp), n_over_max
#   fragment_size_distribution.tsv  sample_id, size, count, fraction (fraction of fragments <= --max-size), rows with count > 0
# Exit 3: no BAM, BAM without a samplesheet row, or a BAM without any proper pair.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(Rsamtools))
library(parallel)
parse_args <- function(args) {
  out <- list()
  for (m in regmatches(args, regexec("^--([^=]+)=(.*)$", args))) if (length(m) == 3) out[[m[2]]] <- m[3]
  out
}
die <- function(...) { message("ERROR: ", ...); quit(save = "no", status = 3) }
opt <- parse_args(commandArgs(trailingOnly = TRUE))
bam_dir <- if (is.null(opt[["bam-dir"]])) "." else opt[["bam-dir"]]
cap     <- if (is.null(opt[["max-size"]])) 2000L else as.integer(opt[["max-size"]])
prefix  <- if (is.null(opt[["out-prefix"]])) "" else opt[["out-prefix"]]
threads <- if (is.null(opt[["threads"]])) 1L else max(1L, as.integer(opt[["threads"]]))
if (is.null(opt$samplesheet)) die("usage: qc_fragment_sizes_table.R --samplesheet=<csv> [--bam-dir=.] [--max-size=2000] [--threads=1]")

ss <- read.csv(opt$samplesheet, colClasses = "character", check.names = FALSE, fileEncoding = "UTF-8-BOM")
need <- c("sample_id", "antibody", "condition", "replicate", "group")
if (!all(need %in% names(ss))) die("samplesheet lacks column(s): ", paste(setdiff(need, names(ss)), collapse = ", "))

bams <- sort(list.files(bam_dir, pattern = "\\.final\\.bam$", full.names = TRUE))
if (length(bams) == 0) die("no *.final.bam in ", bam_dir)
ids <- sub("\\.final\\.bam$", "", basename(bams))
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

# one forked worker per BAM (mc.preschedule = FALSE), at most `threads` at a time; a failed worker returns a "try-error"
counted <- mclapply(seq_along(bams), function(i) {
  message("fragment sizes: ", ids[i])
  count_sizes(bams[i])
}, mc.cores = threads, mc.preschedule = FALSE)
for (i in seq_along(bams)) if (inherits(counted[[i]], "try-error")) die("fragment size counting failed for ", ids[i], ": ", as.character(counted[[i]]))

summ <- list(); dist <- list()
for (i in seq_along(bams)) {
  r <- counted[[i]]
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
    pct_gt_300 = round(100 * sum(cnt[sz > 300]) / n, 2), n_over_max = over)
  keep <- cnt > 0
  dist[[i]] <- data.frame(sample_id = ids[i], size = sz[keep], count = cnt[keep], fraction = signif(cnt[keep] / n, 6))
}
summ <- merge(ss[, need], do.call(rbind, summ), by = "sample_id")
summ <- summ[match(ids, summ$sample_id), ]
write.table(summ, paste0(prefix, "fragment_sizes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, dist), paste0(prefix, "fragment_size_distribution.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote fragment_sizes.tsv and fragment_size_distribution.tsv for", length(ids), "samples\n")
