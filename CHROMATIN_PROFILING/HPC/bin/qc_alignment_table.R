#!/usr/bin/env Rscript
# qc_alignment_table.R -- per-sample alignment statistics table (compute half of the alignment QC).
#
# Called by QC_ALIGNMENT_TABLE (modules/local/qc_alignment_table.nf):
#   qc_alignment_table.R --stats-dir=. --samplesheet=samplesheet.csv --spike-prefix=ecoli_chr --mito-chrom=chrM
#
# Arguments (--key=value):
#   --stats-dir     directory holding the per-sample files (default "."), named by sample_id:
#                     <S>.raw.flagstat.txt  <S>.raw.idxstats.txt     RAW stats: the sorted, unfiltered alignment
#                                                                     (absent for pre-aligned BAM samples -> NA)
#                     <S>.final.flagstat.txt <S>.final.idxstats.txt  FINAL stats: the BAM used downstream
#                     <S>.markdup_metrics.txt                         Picard MarkDuplicates metrics (optional)
#   --samplesheet   samplesheet.csv; columns sample_id, antibody, condition, replicate, group. Sample metadata
#                   comes only from here, never from file names. Samples without any stats file are not reported.
#   --spike-prefix  contigs whose name starts with this prefix count as spike-in (empty = no spike-in)
#   --mito-chrom    mitochondrial contig name (default chrM)
#   --out           output table (default alignment_stats.tsv)
# Output: alignment_stats.tsv (tab separated, one row per sample; NA where a source file is absent)
#   raw_*        counted on the RAW alignment (before MAPQ/flag filtering, duplicate removal, chrM/spike removal)
#   n_target_mapped = raw mapped reads on contigs that are neither chrM nor spike-in
#   final_*      counted on the FINAL BAM; pct_final_of_raw_mapped = 100 * final_mapped_reads / raw_mapped_reads
#   dup_*        Picard MarkDuplicates read pairs (examined, duplicate incl. optical, optical)
# Exit 3 on data problems (missing columns, no samples with stats).

options(stringsAsFactors = FALSE)
parse_args <- function(args) {
  out <- list()
  for (m in regmatches(args, regexec("^--([^=]+)=(.*)$", args))) if (length(m) == 3) out[[m[2]]] <- m[3]
  out
}
die <- function(...) { message("ERROR: ", ...); quit(save = "no", status = 3) }
opt <- parse_args(commandArgs(trailingOnly = TRUE))
stats_dir    <- if (is.null(opt[["stats-dir"]])) "." else opt[["stats-dir"]]
ss_path      <- opt[["samplesheet"]]
spike_prefix <- if (is.null(opt[["spike-prefix"]])) "" else opt[["spike-prefix"]]
mito_chrom   <- if (is.null(opt[["mito-chrom"]])) "chrM" else opt[["mito-chrom"]]
out_path     <- if (is.null(opt[["out"]])) "alignment_stats.tsv" else opt[["out"]]
if (is.null(ss_path)) die("usage: qc_alignment_table.R --samplesheet=<csv> [--stats-dir=.] [--spike-prefix=P] [--mito-chrom=chrM]")

ss <- read.csv(ss_path, colClasses = "character", check.names = FALSE, fileEncoding = "UTF-8-BOM")
need <- c("sample_id", "antibody", "condition", "replicate", "group")
if (!all(need %in% names(ss))) die("samplesheet lacks column(s): ", paste(setdiff(need, names(ss)), collapse = ", "))
ss <- ss[, need]

# ---- parsers (all return NA when a file is absent or unparsable) -------------
first_num <- function(line) if (length(line) == 0 || is.na(line)) NA_real_ else suppressWarnings(as.numeric(strsplit(trimws(line), " ")[[1]][1]))
parse_flagstat <- function(path) {
  if (!file.exists(path)) return(list(total = NA_real_, mapped = NA_real_))
  l <- readLines(path, warn = FALSE)
  list(total  = first_num(l[grepl("in total", l)][1]),
       mapped = first_num(l[grepl("^[0-9]+ \\+ [0-9]+ mapped \\(", l)][1]))
}
parse_idxstats <- function(path, spike_prefix, mito) {
  na <- list(chrM = NA_real_, spike = NA_real_, mapped = NA_real_)
  if (!file.exists(path) || file.info(path)$size == 0) return(na)
  df <- tryCatch(read.delim(path, header = FALSE, col.names = c("chr", "len", "mapped", "unmapped"), quote = ""),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(na)
  df <- df[df$chr != "*", ]
  is_spike <- if (nzchar(spike_prefix)) startsWith(df$chr, spike_prefix) else rep(FALSE, nrow(df))
  list(chrM = sum(df$mapped[df$chr == mito], na.rm = TRUE),
       spike = sum(df$mapped[is_spike], na.rm = TRUE),
       mapped = sum(df$mapped, na.rm = TRUE))
}
parse_picard <- function(path) {
  na <- list(examined = NA_real_, dup = NA_real_, optical = NA_real_)
  if (!file.exists(path)) return(na)
  l <- readLines(path, warn = FALSE)
  h <- which(grepl("READ_PAIRS_EXAMINED", l))
  if (length(h) == 0 || h[1] + 1 > length(l) || !nzchar(trimws(l[h[1] + 1]))) return(na)
  hdr <- strsplit(l[h[1]], "\t")[[1]]; v <- strsplit(l[h[1] + 1], "\t")[[1]]
  if (length(v) < length(hdr)) return(na)
  names(v) <- hdr
  num <- function(k) suppressWarnings(as.numeric(v[k]))
  list(examined = num("READ_PAIRS_EXAMINED"), dup = num("READ_PAIR_DUPLICATES"), optical = num("READ_PAIR_OPTICAL_DUPLICATES"))
}
pct <- function(a, b) ifelse(is.na(a) | is.na(b) | b == 0, NA_real_, round(100 * a / b, 2))

has_stats <- function(s) any(file.exists(file.path(stats_dir, paste0(s, c(".raw.flagstat.txt", ".final.flagstat.txt")))))
ss <- ss[vapply(ss$sample_id, has_stats, logical(1)), , drop = FALSE]
if (nrow(ss) == 0) die("no sample in ", ss_path, " has a <sample_id>.raw.flagstat.txt or <sample_id>.final.flagstat.txt in ", stats_dir)

rows <- lapply(ss$sample_id, function(s) {
  f <- function(x) file.path(stats_dir, paste0(s, x))
  rf <- parse_flagstat(f(".raw.flagstat.txt"));   ri <- parse_idxstats(f(".raw.idxstats.txt"), spike_prefix, mito_chrom)
  ff <- parse_flagstat(f(".final.flagstat.txt")); fi <- parse_idxstats(f(".final.idxstats.txt"), spike_prefix, mito_chrom)
  pd <- parse_picard(f(".markdup_metrics.txt"))
  raw_mapped <- if (!is.na(rf$mapped)) rf$mapped else ri$mapped
  if (is.na(rf$total) && is.na(ff$total)) message("WARNING: no parsable flagstat for ", s)
  data.frame(
    raw_total_reads = rf$total, raw_mapped_reads = raw_mapped,
    pct_mapped = pct(raw_mapped, rf$total),
    n_chrM = ri$chrM, pct_chrM = pct(ri$chrM, raw_mapped),
    n_spike_in = ri$spike, pct_spike_in = pct(ri$spike, raw_mapped),
    n_target_mapped = pmax(raw_mapped - ri$chrM - ri$spike, 0),
    pct_unmapped = pct(rf$total - raw_mapped, rf$total),
    final_total_reads = ff$total, final_mapped_reads = ff$mapped,
    final_n_chrM = fi$chrM, final_n_spike_in = fi$spike,
    pct_final_of_raw_mapped = pct(ff$mapped, raw_mapped),
    dup_pairs_examined = pd$examined, dup_pairs_duplicate = pd$dup, dup_pairs_optical = pd$optical,
    pct_duplication = pct(pd$dup, pd$examined))
})
res <- cbind(ss, do.call(rbind, rows))
write.table(res, out_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat("Wrote", out_path, "for", nrow(res), "samples;", sum(is.na(res$raw_total_reads)), "without raw stats\n")
