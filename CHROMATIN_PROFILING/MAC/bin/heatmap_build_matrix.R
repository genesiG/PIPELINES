#!/usr/bin/env Rscript
# heatmap_build_matrix.R -- build step of the heatmap module (BUILD_MATRIX process, container cutandrun-r_1.1.sif).
# Port of the legacy helper Scripts/utils/build_matrix.R to a --key=value CLI.
#
# CLI contract (--key=value; comma lists):
#   --regions-bed=<F.bed>        BED with >= 3 columns (chr start end [name score strand ...]). Lines starting with '#',
#                                'track' or 'browser', and blank lines, are skipped. Row order is kept (no sorting).
#                                name = column 4 if present, else region_<row index> (1-based over the data rows).
#                                Names must be unique (exit 3 listing the duplicates).
#                                strand = column 6 if present (+, -, . ; '.' -> '*'), else '*'.
#   --bigwigs=a.bw,b.bw          bigWig file per sample (same count as --names)
#   --names=id1,id2              sample ids (keys of the norm_mats list), same order as --bigwigs
#   --window=W                   half-window in bp; matrices span -W..+W around the centre point
#   --bin-size=B                 bin size in bp; W must be an integer multiple of B (exit 3 otherwise), n_bins = W/B per side
#   --out-prefix=NAME            output prefix; files are written to the working directory
#   --cluster-samples=id1,id2    optional subset of --names used for the cluster matrix; empty = none
#   --cpus=1                     parallel workers (parallel::mclapply, one sample per worker)
#
# Method: centre point = resize(targets, width = 1, fix = "center") keeping strand (EnrichedHeatmap's normalizeToMatrix
# flips minus-strand regions). Reference-point mode only: include_target = FALSE, mean_mode = "w0", smooth = FALSE,
# value_column = "score", w = B, extend = W.
#
# Outputs (NAME = --out-prefix):
#   NAME.norm_mats.rds           named list (names = --names) of normalizedMatrix objects; one row per BED row, BED order
#   NAME.targets.rds             GRanges of the BED rows with original names; mcols name, score (NA if absent),
#                                bed_strand (column 6 as read, NA if absent); strand slot = '+', '-' or '*'
#   NAME.cluster_matrix.tsv.gz   only when --cluster-samples is non-empty: cbind of those samples' matrices, NA -> 0,
#                                tab-separated, no header, no row names, row order = BED order
#
# Exit codes: 0 ok; 2 usage (bad or missing argument); 3 data problem (missing file, --names/--bigwigs count mismatch,
#   duplicate names, bad BED record, W not a multiple of B, sample failed, matrix row count != BED row count).

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(rtracklayer)
  library(EnrichedHeatmap)
  library(GenomicRanges)
  library(parallel)     # mclapply: fork-based parallelism on Linux
})

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

allowed <- c("regions-bed", "bigwigs", "names", "window", "bin-size", "out-prefix", "cluster-samples", "cpus")
args <- parse_args(commandArgs(trailingOnly = TRUE), allowed)
bed_file   <- get_req(args, "regions-bed")
bw_files   <- split_list(get_req(args, "bigwigs"))
sample_ids <- split_list(get_req(args, "names"))
out_prefix <- get_req(args, "out-prefix")
cl_samples <- split_list(get_opt(args, "cluster-samples", ""))

window_bp <- suppressWarnings(as.integer(get_req(args, "window")))
bin_bp    <- suppressWarnings(as.integer(get_req(args, "bin-size")))
cpus      <- suppressWarnings(as.integer(get_opt(args, "cpus", "1")))
if (is.na(window_bp) || window_bp < 1L) die("--window must be a positive integer (bp)", 2L)
if (is.na(bin_bp) || bin_bp < 1L) die("--bin-size must be a positive integer (bp)", 2L)
if (is.na(cpus) || cpus < 1L) die("--cpus must be a positive integer", 2L)

if (length(bw_files) != length(sample_ids))
  die(sprintf("--names (%d) and --bigwigs (%d) must have the same count",
              length(sample_ids), length(bw_files)), 3L)
if (anyDuplicated(sample_ids))
  die(paste0("duplicate --names: ", paste(unique(sample_ids[duplicated(sample_ids)]), collapse = ", ")), 3L)
if (window_bp %% bin_bp != 0L)
  die(sprintf("--window (%d) is not a multiple of --bin-size (%d): n_bins would not be an integer",
              window_bp, bin_bp), 3L)
n_bins <- window_bp %/% bin_bp

if (!file.exists(bed_file)) die(paste0("regions BED not found: ", bed_file), 3L)
for (bw in bw_files) if (!file.exists(bw)) die(paste0("bigWig not found: ", bw), 3L)
missing_cs <- setdiff(cl_samples, sample_ids)
if (length(missing_cs) > 0L)
  die(paste0("--cluster-samples not found in --names: ", paste(missing_cs, collapse = ", ")), 3L)

cat("=== heatmap_build_matrix.R ===\n")
cat(sprintf("  Prefix     : %s\n", out_prefix))
cat(sprintf("  Regions    : %s\n", bed_file))
cat(sprintf("  Window     : +/- %d bp   Bin: %d bp   Bins per side: %d\n", window_bp, bin_bp, n_bins))
cat(sprintf("  Samples    : %s\n", paste(sample_ids, collapse = ", ")))

# -----------------------------------------------------------------------------
# Load regions (BED order kept)
# -----------------------------------------------------------------------------
read_bed <- function(path) {
  lines <- readLines(path, warn = FALSE)
  skip  <- !nzchar(trimws(lines)) | grepl("^(#|track|browser)", lines)
  lines <- lines[!skip]
  if (length(lines) == 0L) die(paste0("no BED records in ", path), 3L)

  fld <- strsplit(lines, "\t", fixed = TRUE)
  n_fld <- vapply(fld, length, integer(1))
  if (any(n_fld < 3L))
    die(sprintf("BED data record %d has fewer than 3 columns", which(n_fld < 3L)[1]), 3L)
  col_of <- function(j) vapply(fld, function(f) if (length(f) >= j) f[j] else NA_character_, character(1))

  chr   <- col_of(1L)
  start <- suppressWarnings(as.integer(col_of(2L)))
  end   <- suppressWarnings(as.integer(col_of(3L)))
  if (anyNA(start) || anyNA(end)) die("BED start/end must be integers", 3L)

  reg_names <- col_of(4L)
  reg_names <- ifelse(is.na(reg_names) | !nzchar(reg_names),
                      paste0("region_", seq_along(reg_names)), reg_names)
  dup <- unique(reg_names[duplicated(reg_names)])
  if (length(dup) > 0L)
    die(paste0("duplicate region names in BED (names must be unique): ",
               paste(dup, collapse = ", ")), 3L)

  score <- suppressWarnings(as.numeric(col_of(5L)))

  strand_raw <- col_of(6L)
  if (all(is.na(strand_raw))) {
    strand_raw <- rep(NA_character_, length(lines))
    gr_strand  <- rep("*", length(lines))
  } else {
    if (anyNA(strand_raw)) die("BED column 6 (strand) present for some records only", 3L)
    bad_strand <- setdiff(unique(strand_raw), c("+", "-", "."))
    if (length(bad_strand) > 0L)
      die(paste0("BED strand must be +, - or . ; found: ", paste(bad_strand, collapse = ", ")), 3L)
    gr_strand <- ifelse(strand_raw == ".", "*", strand_raw)
  }

  bed <- data.frame(chr = chr, start = start, end = end, stringsAsFactors = FALSE)
  targets <- makeGRangesFromDataFrame(bed,
               seqnames.field          = "chr",
               start.field             = "start",
               end.field               = "end",
               keep.extra.columns      = FALSE,
               starts.in.df.are.0based = TRUE)
  strand(targets) <- gr_strand
  names(targets)  <- reg_names
  targets$name    <- reg_names
  targets$score   <- score
  targets$bed_strand <- strand_raw
  targets
}

targets <- read_bed(bed_file)
cat(sprintf("  Regions loaded: %d\n", length(targets)))

# -----------------------------------------------------------------------------
# Reference points (centre, strand kept) and import windows
# -----------------------------------------------------------------------------
targets_ref <- resize(targets, width = 1L, fix = "center")
names(targets_ref) <- names(targets)
import_windows <- resize(targets_ref, width = 2L * window_bp, fix = "center")

# -----------------------------------------------------------------------------
# Build normalizedMatrix per sample
# -----------------------------------------------------------------------------
build_nm <- function(label, bw_path) {
  cat(sprintf("  [%s] Importing bigWig ...\n", label))
  bw_signal <- import(bw_path,
                      format    = "BigWig",
                      selection = BigWigSelection(import_windows))
  cat(sprintf("  [%s] Building normalizedMatrix ...\n", label))
  mat <- normalizeToMatrix(
    signal         = bw_signal,
    target         = targets_ref,
    extend         = window_bp,
    w              = bin_bp,
    value_column   = "score",
    mean_mode      = "w0",
    include_target = FALSE,
    smooth         = FALSE
  )
  cat(sprintf("  [%s] Dimensions: %d x %d\n", label, nrow(mat), ncol(mat)))
  mat
}

cat(sprintf("  Parallelism : cpus=%d (mclapply fork)\n", cpus))
raw_results <- mclapply(
  seq_along(sample_ids),
  function(i) try(build_nm(sample_ids[i], bw_files[i]), silent = TRUE),
  mc.cores  = cpus,
  mc.silent = FALSE
)

failed <- vapply(raw_results, inherits, logical(1), what = "try-error")
if (any(failed)) {
  for (i in which(failed)) cat(as.character(raw_results[[i]]), file = stderr())
  die(paste0("the following samples failed: ", paste(sample_ids[failed], collapse = ", ")), 3L)
}
for (i in seq_along(sample_ids)) {
  if (nrow(raw_results[[i]]) != length(targets))
    die(sprintf("sample '%s': matrix has %d rows but the BED has %d rows",
                sample_ids[i], nrow(raw_results[[i]]), length(targets)), 3L)
}

norm_mats <- setNames(raw_results, sample_ids)

# -----------------------------------------------------------------------------
# Save RDS outputs
# -----------------------------------------------------------------------------
mats_file <- paste0(out_prefix, ".norm_mats.rds")
tgt_file  <- paste0(out_prefix, ".targets.rds")
saveRDS(norm_mats, mats_file)
saveRDS(targets,   tgt_file)
cat(sprintf("  Saved: %s / %s\n", basename(mats_file), basename(tgt_file)))

# -----------------------------------------------------------------------------
# Optional: cluster matrix (subset of samples) as gzip TSV, NA -> 0
# -----------------------------------------------------------------------------
if (length(cl_samples) > 0L) {
  cat(sprintf("  Building cluster matrix from: %s\n", paste(cl_samples, collapse = " + ")))
  cluster_mat <- do.call(
    cbind,
    lapply(cl_samples, function(s) as.matrix(norm_mats[[s]]))
  )
  cluster_mat[is.na(cluster_mat)] <- 0
  cat(sprintf("  Cluster matrix dimensions: %d x %d\n",
              nrow(cluster_mat), ncol(cluster_mat)))

  gz_path <- paste0(out_prefix, ".cluster_matrix.tsv.gz")
  gz_con  <- gzfile(gz_path, "w")
  write.table(cluster_mat, gz_con, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  close(gz_con)
  cat(sprintf("  Cluster matrix saved: %s\n", basename(gz_path)))
}

cat(sprintf("=== heatmap_build_matrix.R complete (%s) ===\n\n", out_prefix))
