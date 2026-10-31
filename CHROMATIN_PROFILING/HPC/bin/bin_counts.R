#!/usr/bin/env Rscript
# bin_counts.R -- one csaw bin-count pass over all active libraries, shared by TMM, the normalisation QC and the replicate
# concordance steps (counting is done once; every consumer subsets the columns it needs).
#
# CLI contract (--key=value; called by the BIN_COUNTS module, one task per distinct bin size):
#   --bams=<comma list>          final target BAMs (coordinate-sorted + indexed)
#   --sample-ids=<comma list>    sample ids, same order/length as --bams (used verbatim as count column names)
#   --idxstats=<comma list>      samtools idxstats files of the same final BAMs, same order/length as --bams
#                                (bam_mapped_reads = sum of column 3 over all contigs)
#   --bin-size=<int>             bin width in bp
#   --mapq=<int>                 csaw readParam(minq=)
#   --dedup=TRUE|FALSE           csaw readParam(dedup=): TRUE drops reads carrying the duplicate flag (0x400)
#   --max-frag=<int>             csaw readParam(max.frag=): fragments longer than this are dropped (paired-end, pe = "both")
#   --blacklist=<bed>            OPTIONAL 0-based BED (chr, start, end). A bin overlapping it by >= 1 bp gets blacklisted = TRUE.
#                                Absent/empty: all bins FALSE. Converted with IRanges(start + 1, end) (csaw bins are 1-based)
#   --cpus=<int>                 OPTIONAL, default 1: BiocParallel MulticoreParam workers for windowCounts (counts are
#                                integers, so results do not depend on the worker count)
#   --out-prefix=<prefix>        output prefix, by convention bins_<bin size>
# Outputs:
#   <prefix>.counts.tsv.gz   tab-separated, header chr, start, end, blacklisted, then one integer column per sample id.
#                            ALL bins (csaw filter = 0; no bin is dropped, no blacklist exclusion). start is the 1-based
#                            csaw start, end the csaw end (inclusive); rows follow the csaw genome order.
#   <prefix>.libsizes.tsv    sample_id, csaw_total, bam_mapped_reads, bam_fragments, bin_size, mapq, max_frag, dedup
#                            csaw_total       = csaw totals (fragments counted under the readParam above) = edgeR lib.size
#                            bam_mapped_reads = sum of idxstats column 3
#                            bam_fragments    = bam_mapped_reads / 2
# Exit codes: 0 ok; 2 usage error; 3 data problem (missing BAM/index/idxstats, length mismatch, duplicated ids).
#
# Consumers subset the columns of their run and apply their own bin filter (e.g. TMM: rowSums >= 10), which gives the same
# counts and totals as a per-run windowCounts(filter = 10).

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
idxstats_mapped <- function(path) {
  if (!file.exists(path)) die(paste0("idxstats file not found: ", path))
  d <- tryCatch(read.table(path, sep = "\t", header = FALSE, quote = "", comment.char = "", fill = TRUE,
                           colClasses = "character"), error = function(e) NULL)
  if (is.null(d) || ncol(d) < 3L) die(paste0(path, " is empty or not a samtools idxstats file"))
  v <- suppressWarnings(as.numeric(d[[3]]))
  if (anyNA(v)) die(paste0(path, ": non-numeric mapped count in column 3"))
  sum(v)
}

allowed <- c("bams", "sample-ids", "idxstats", "bin-size", "mapq", "dedup", "max-frag", "blacklist", "cpus", "out-prefix")
args <- parse_args(commandArgs(trailingOnly = TRUE), allowed)
bams       <- split_list(get_req(args, "bams"))
ids        <- split_list(get_req(args, "sample-ids"))
idx_files  <- split_list(get_req(args, "idxstats"))
bin_size   <- as_num(get_req(args, "bin-size"), "bin-size")
mapq       <- as_num(get_req(args, "mapq"), "mapq")
dedup      <- as_bool(get_req(args, "dedup"), "dedup")
max_frag   <- as_num(get_req(args, "max-frag"), "max-frag")
bl_file    <- get_opt(args, "blacklist", "")
cpus       <- as_num(get_opt(args, "cpus", "1"), "cpus")
out_prefix <- get_req(args, "out-prefix")

if (bin_size < 1 || bin_size != round(bin_size)) die("--bin-size must be a positive integer", 2L)
if (cpus < 1 || cpus != round(cpus)) die("--cpus must be a positive integer", 2L)
if (length(bams) != length(ids)) die(sprintf("--bams has %d entries but --sample-ids has %d", length(bams), length(ids)))
if (length(idx_files) != length(ids)) die(sprintf("--idxstats has %d entries but --sample-ids has %d", length(idx_files), length(ids)))
if (length(bams) < 1L) die("no libraries given")
if (anyDuplicated(ids)) die(paste0("duplicated sample ids: ", paste(ids[duplicated(ids)], collapse = ", ")))
check_bams(bams)
if (nzchar(bl_file) && !file.exists(bl_file)) die(paste0("blacklist file not found: ", bl_file))
mapped <- vapply(idx_files, idxstats_mapped, numeric(1))

suppressPackageStartupMessages({
  library(csaw)
  library(GenomicRanges)
  library(IRanges)
  library(BiocParallel)
})

param <- readParam(minq = mapq, pe = "both", max.frag = max_frag, dedup = dedup)
bp <- if (cpus > 1) MulticoreParam(workers = as.integer(cpus), progressbar = FALSE) else SerialParam()
cat("Counting ", bin_size, " bp bins in ", length(bams), " libraries (minq=", mapq, ", pe=both, max.frag=", max_frag,
    ", dedup=", dedup, ", cpus=", cpus, ")\n", sep = "")
counts <- windowCounts(bams, bin = TRUE, width = bin_size, filter = 0, param = param, BPPARAM = bp)
totals <- counts$totals
gr <- rowRanges(counts)

blk <- rep(FALSE, length(gr))
if (nzchar(bl_file) && file.size(bl_file) > 0) {
  bl <- read.table(bl_file, sep = "\t", header = FALSE, quote = "", comment.char = "", colClasses = c("character", "numeric", "numeric"),
                   fill = TRUE)[, 1:3]
  if (nrow(bl)) {
    bl_gr <- GRanges(bl[[1]], IRanges(start = bl[[2]] + 1, end = bl[[3]]))   # BED is 0-based half-open
    blk <- overlapsAny(gr, bl_gr, ignore.strand = TRUE)
  }
}
cat("Bins: ", length(gr), "; blacklisted: ", sum(blk), "\n", sep = "")

cm <- assay(counts, "counts")
storage.mode(cm) <- "integer"
colnames(cm) <- ids
tab <- data.frame(chr = as.character(seqnames(gr)), start = start(gr), end = end(gr), blacklisted = blk,
                  check.names = FALSE)
tab <- cbind(tab, as.data.frame(cm, check.names = FALSE))
con <- gzfile(paste0(out_prefix, ".counts.tsv.gz"), "w")
write.table(tab, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
close(con)

fmt0 <- function(x) sprintf("%.0f", x)   # never scientific notation (as.character(1e5) is "1e+05")
lib <- data.frame(sample_id = ids, csaw_total = fmt0(totals), bam_mapped_reads = fmt0(mapped),
                  bam_fragments = vapply(mapped / 2, format, character(1), scientific = FALSE), bin_size = as.integer(bin_size), mapq = as.integer(mapq),
                  max_frag = as.integer(max_frag), dedup = dedup)
write.table(lib, paste0(out_prefix, ".libsizes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Done: ", out_prefix, ".counts.tsv.gz, ", out_prefix, ".libsizes.tsv\n", sep = "")
