#!/usr/bin/env Rscript
#' ##########################################################################
#' guide_representation_wrangle.R
#'
#' Computation half of the guide-representation QC: reads one MAGeCK count
#' table and writes every table the plots (and labmates) need — no plotting.
#' Its sibling, guide_representation_plot.R, renders the plots from these
#' tables, so plot styling can be changed and re-run without recomputing
#' anything here.
#'
#' USAGE:
#'   Rscript guide_representation_wrangle.R --count-table=<X.count.txt> \
#'       [--library=<mageck_library.txt>] [--mageck-log=<X_mageck_count.log>] \
#'       [--unmapped=<X.unmapped.txt>] [--outdir=<dir>]
#'
#'   --library     MAGeCK-format library (sgRNA_id, sequence, gene); enables
#'                 the GC-content table and (with --unmapped) read
#'                 classification. Optional.
#'   --mageck-log  mageck count log; enables the mapping-rate table and (with
#'                 --library and --unmapped) read classification. Optional.
#'   --unmapped    mageck count's --unmapped-to-file output (sequence,
#'                 sequence, one count column per sample); enables read
#'                 classification. Optional.
#'   --outdir      where to write tables (default: current directory)
#'
#' OUTPUTS (<prefix> = count-table name without .count.txt):
#'   <prefix>_counts_long.csv            sgRNA, Gene, sample, count, is_control
#'                                       (feeds ECDF, density, control-vs-targeting)
#'   <prefix>_ranked.csv                 sgRNA, Gene, avg_count, rank
#'   <prefix>_lorenz.csv                 sample, count, cum_gRNAs, cum_reads
#'   <prefix>_gini.csv                   sample, Gini
#'   <prefix>_percentile_ratio.csv       sample, p10, p90, ratio_90_10
#'   <prefix>_undetected.csv             sample, total_guides, zero_guides, pct_undetected
#'   <prefix>_mapping_rate.csv           label, reads, mappedreads, mapping_pct   (if --mageck-log)
#'   <prefix>_gc_vs_depth.csv            sgRNA, sample, count, norm_count, gc_pct (if --library)
#'   <prefix>_replicate_correlations.csv Sample_A, Sample_B, Pearson, Spearman   (if >= 2 samples)
#'   <prefix>_dropout_summary.csv        per-sample dropout thresholds
#'   <prefix>_zero_count_guides.csv      sample, sgRNA, Gene                     (if any)
#'   <prefix>_read_classification.csv    sample, reads, mappedreads, pct_exact_match,
#'                                       near_miss_reads, pct_near_miss, far_miss_reads,
#'                                       pct_far_miss (if --library, --mageck-log, --unmapped)
#'   <prefix>_qc_summary.csv             one row per sample, all headline metrics + flags
#' ##########################################################################

options(stringsAsFactors = FALSE)

suppressMessages({
  library(dplyr, quietly = TRUE)
  library(tidyr, quietly = TRUE)
})

# ---- Arguments (--key=value) ----------------------------------------------
parse_args <- function(args) {
  out <- list()
  for (m in regmatches(args, regexec("^--([^=]+)=(.*)$", args))) {
    if (length(m) == 3) out[[m[2]]] <- m[3]
  }
  out
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(opt[["count-table"]])) {
  stop("Usage: guide_representation_wrangle.R --count-table=<X.count.txt> [--library=<lib.txt>] [--mageck-log=<log>] [--outdir=<dir>]")
}
count_path <- opt[["count-table"]]
outdir     <- if (is.null(opt$outdir)) "." else opt$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

pool_prefix <- sub("\\.count\\.txt$", "", basename(count_path))
out_file    <- function(name) file.path(outdir, paste0(pool_prefix, "_", name, ".csv"))
save_table  <- function(df, name) {
  write.csv(df, out_file(name), row.names = FALSE)
  cat("  Saved:", out_file(name), "\n")
}

# ---- Helpers --------------------------------------------------------------
calc_gini <- function(x) {
  x <- sort(x[x > 0])
  n <- length(x)
  if (n == 0 || sum(x) == 0) return(NA_real_)
  sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
}

gc_content <- function(seq) {
  seq <- toupper(seq)
  gc  <- nchar(gsub("[^GC]", "", seq))
  gc / nchar(seq) * 100
}

# Parse a MAGeCK count log for per-file mapping statistics. Anchors to the
# structured "Summary of file" block so it never picks up fields from the
# Parameters: line (which contains --sample-label, etc.).
parse_mageck_log <- function(log_path) {
  if (is.null(log_path) || !file.exists(log_path)) return(NULL)
  lines <- readLines(log_path, warn = FALSE)
  summary_starts <- grep("Summary of file", lines, fixed = TRUE)
  if (length(summary_starts) == 0) return(NULL)

  get_field <- function(block, key) {
    pattern <- paste0(":\\s+", key, "\\s+(\\S+)")
    m <- regmatches(block, regexec(pattern, block))
    hits <- Filter(function(x) length(x) == 2, m)
    if (length(hits) == 0) return(NA_character_)
    trimws(hits[[1]][2])
  }

  do.call(rbind, lapply(summary_starts, function(start_idx) {
    block       <- lines[start_idx:min(start_idx + 10L, length(lines))]
    reads       <- suppressWarnings(as.numeric(get_field(block, "reads")))
    mappedreads <- suppressWarnings(as.numeric(get_field(block, "mappedreads")))
    data.frame(
      label       = get_field(block, "label"),
      reads       = reads,
      mappedreads = mappedreads,
      mapping_pct = round(mappedreads / max(reads, 1) * 100, 2)
    )
  }))
}

# ---- Load -----------------------------------------------------------------
cat("Processing:", count_path, "\n")
df <- read.table(count_path, header = TRUE, sep = "\t", check.names = FALSE)
if (!all(c("sgRNA", "Gene") %in% colnames(df))) {
  stop("Count table ", count_path, " is missing the expected sgRNA/Gene columns.")
}

sample_cols    <- setdiff(colnames(df), c("sgRNA", "Gene"))
active_samples <- sample_cols[colSums(df[, sample_cols, drop = FALSE], na.rm = TRUE) > 0]
if (length(active_samples) == 0) stop("No active samples (>0 reads) in ", count_path)
cat("  Active samples:", paste(active_samples, collapse = ", "), "\n")

df_long <- df |>
  pivot_longer(cols = all_of(active_samples), names_to = "sample", values_to = "count") |>
  mutate(is_control = grepl("Control|NTC|Non-targeting", Gene, ignore.case = TRUE)) |>
  select(sgRNA, Gene, sample, count, is_control)
save_table(df_long, "counts_long")

# ---- Ranked guide representation ------------------------------------------
df_ranked <- df |>
  mutate(avg_count = rowMeans(across(all_of(active_samples)), na.rm = TRUE)) |>
  arrange(avg_count) |>
  mutate(rank = row_number()) |>
  select(sgRNA, Gene, avg_count, rank)
save_table(df_ranked, "ranked")

# ---- Lorenz curves + Gini -------------------------------------------------
lorenz_df <- df_long |>
  group_by(sample) |>
  arrange(count, .by_group = TRUE) |>
  mutate(
    cum_gRNAs = row_number() / n(),
    cum_reads = cumsum(count) / max(sum(count), 1)
  ) |>
  ungroup() |>
  select(sample, count, cum_gRNAs, cum_reads)
save_table(lorenz_df, "lorenz")

gini_df <- df_long |>
  group_by(sample) |>
  summarize(Gini = round(calc_gini(count), 4), .groups = "drop")
save_table(gini_df, "gini")

# ---- 90/10 percentile ratio -----------------------------------------------
pct_ratio_df <- df_long |>
  group_by(sample) |>
  summarize(
    p10         = quantile(count[count > 0], 0.10, na.rm = TRUE),
    p90         = quantile(count[count > 0], 0.90, na.rm = TRUE),
    ratio_90_10 = round(
      quantile(count[count > 0], 0.90, na.rm = TRUE) /
        max(quantile(count[count > 0], 0.10, na.rm = TRUE), 1),
      2
    ),
    .groups = "drop"
  )
save_table(pct_ratio_df, "percentile_ratio")

# ---- Undetected guides ----------------------------------------------------
undetected_df <- df_long |>
  group_by(sample) |>
  summarize(
    total_guides   = n(),
    zero_guides    = sum(count == 0),
    pct_undetected = round(sum(count == 0) / n() * 100, 2),
    .groups = "drop"
  )
save_table(undetected_df, "undetected")

# ---- Mapping rate (from MAGeCK log) ---------------------------------------
mapping_df <- parse_mageck_log(opt[["mageck-log"]])
if (!is.null(mapping_df) && nrow(mapping_df) > 0) {
  save_table(mapping_df, "mapping_rate")
} else {
  cat("  [skip] mapping rate: no usable --mageck-log\n")
  mapping_df <- NULL
}

# ---- Library sequences (used by GC content below and read classification) -
lib_seqs <- NULL
lib_path <- opt$library
if (!is.null(lib_path) && file.exists(lib_path)) {
  lib_raw <- read.table(lib_path, header = TRUE, sep = "\t", check.names = FALSE)
  id_col  <- intersect(c("sgRNA_id", "sgRNA", "id"), colnames(lib_raw))[1]
  if ("sequence" %in% colnames(lib_raw) && !is.na(id_col)) {
    lib_seqs <- lib_raw |> rename(sgRNA = all_of(id_col)) |> select(sgRNA, sequence)
  } else {
    cat("  [skip] library sequences: --library lacks sgRNA id / sequence columns\n")
  }
} else {
  cat("  [skip] library sequences: no --library\n")
}

# ---- GC content vs median-normalised depth --------------------------------
if (!is.null(lib_seqs)) {
  seq_df <- lib_seqs |>
    mutate(gc_pct = sapply(sequence, gc_content)) |>
    select(sgRNA, gc_pct)
  gc_count_df <- df_long |>
    select(sgRNA, sample, count) |>
    group_by(sample) |>
    mutate(norm_count = count / max(median(count[count > 0], na.rm = TRUE), 1)) |>
    ungroup() |>
    inner_join(seq_df, by = "sgRNA")
  save_table(gc_count_df, "gc_vs_depth")
} else {
  cat("  [skip] GC content: no usable library sequences\n")
}

# ---- Read classification: exact / near-miss / far-miss --------------------
# mageck count's own guide matching is always exact -- there is no
# mismatch-tolerance option in mageck count itself (verified against
# `mageck count --help`), so every mapped read is by construction an exact
# guide match (mapping_rate above already *is* "% exact matches"). This
# classifies the reads mageck COULDN'T map (its --unmapped-to-file output)
# into "near-miss" (exactly 1 substitution away from some guide -- most
# likely lost to sequencing/PCR error, and would map if mismatch-tolerant
# counting existed) vs "far-miss" (not close to any guide in the library --
# contamination, wrong library, adapter-only reads, ...).
unmapped_path      <- opt$unmapped
classification_df  <- NULL
if (is.null(lib_seqs) || is.null(unmapped_path) || !file.exists(unmapped_path) || is.null(mapping_df)) {
  cat("  [skip] read classification: needs --library, --mageck-log, and --unmapped\n")
} else {
  seq_lens <- nchar(lib_seqs$sequence)
  L <- seq_lens[1]
  if (!all(seq_lens == L)) {
    cat("  [skip] read classification: library guides aren't all the same length\n")
  } else if (file.info(unmapped_path)$size == 0) {
    # Empty file: mageck wrote no unmapped reads at all (every sample mapped 100%).
    classification_df <- data.frame(sample = active_samples, unmapped_reads = 0L,
                                     near_miss_reads = 0L, far_miss_reads = 0L)
  } else {
    unmapped_raw  <- read.table(unmapped_path, header = FALSE, sep = "\t", check.names = FALSE)
    n_sample_cols <- ncol(unmapped_raw) - 2
    if (n_sample_cols != length(sample_cols)) {
      cat("  [skip] read classification: --unmapped has", n_sample_cols,
          "sample column(s), expected", length(sample_cols), "(from the count table)\n")
    } else {
      colnames(unmapped_raw) <- c("sequence", "sequence_dup", sample_cols)
      # Hamming<=1 lookup via the standard wildcard-neighborhood trick: every
      # guide contributes, per position, itself with that position blanked;
      # a query sequence hits if any of ITS OWN blanked variants matches one
      # of these. O(n_guides * L) to build, O(n_query * L) to check -- both
      # vectorized and linear, unlike an all-pairs distance comparison.
      index_keys <- unique(unlist(lapply(seq_len(L), function(i) {
        paste0(substr(lib_seqs$sequence, 1, i - 1), ".", substr(lib_seqs$sequence, i + 1, L))
      })))
      query <- unmapped_raw$sequence
      hit <- rep(FALSE, length(query))
      for (i in seq_len(L)) {
        k   <- paste0(substr(query, 1, i - 1), ".", substr(query, i + 1, L))
        hit <- hit | (k %in% index_keys)
      }
      classification_df <- do.call(rbind, lapply(active_samples, function(s) {
        counts <- unmapped_raw[[s]]
        data.frame(sample = s, unmapped_reads = sum(counts),
                   near_miss_reads = sum(counts[hit]), far_miss_reads = sum(counts[!hit]))
      }))
    }
  }
  if (!is.null(classification_df)) {
    classification_df <- classification_df |>
      left_join(mapping_df |> rename(sample = label) |> select(sample, reads, mappedreads), by = "sample") |>
      mutate(
        pct_exact_match = round(mappedreads      / pmax(reads, 1) * 100, 2),
        pct_near_miss   = round(near_miss_reads   / pmax(reads, 1) * 100, 2),
        pct_far_miss    = round(far_miss_reads    / pmax(reads, 1) * 100, 2)
      ) |>
      select(sample, reads, mappedreads, pct_exact_match, near_miss_reads, pct_near_miss,
             far_miss_reads, pct_far_miss)
    save_table(classification_df, "read_classification")
  }
}

# ---- Replicate correlation (needs >= 2 samples) ---------------------------
if (length(active_samples) >= 2) {
  count_mat  <- as.matrix(df[, active_samples, drop = FALSE])
  pearson_m  <- cor(log10(count_mat + 1), method = "pearson",  use = "complete.obs")
  spearman_m <- cor(count_mat,            method = "spearman", use = "complete.obs")
  cor_long <- merge(
    as.data.frame(as.table(pearson_m))  |> rename(Sample_A = Var1, Sample_B = Var2, Pearson = Freq),
    as.data.frame(as.table(spearman_m)) |> rename(Sample_A = Var1, Sample_B = Var2, Spearman = Freq),
    by = c("Sample_A", "Sample_B")
  ) |> mutate(across(c(Pearson, Spearman), \(x) round(x, 4)))
  save_table(cor_long, "replicate_correlations")
} else {
  cat("  [skip] replicate correlation: only 1 active sample\n")
}

# ---- Dropout summary ------------------------------------------------------
dropout_table <- lorenz_df |>
  group_by(sample) |>
  summarize(
    dropout_1pct      = max(count[cum_reads <= 0.01], default = 0L),
    dropout_5pct      = max(count[cum_reads <= 0.05], default = 0L),
    dropout_10pct     = max(count[cum_reads <= 0.10], default = 0L),
    zero_count_guides = sum(count == 0),
    .groups = "drop"
  )
save_table(dropout_table, "dropout_summary")

# ---- Zero-count guides ----------------------------------------------------
zero_guides_df <- df_long |>
  filter(count == 0) |>
  select(sample, sgRNA, Gene) |>
  arrange(sample, Gene, sgRNA)
if (nrow(zero_guides_df) > 0) save_table(zero_guides_df, "zero_count_guides")

# ---- Log2 uniformity (% of guides within 2 of log2(mean)) -----------------
log2_metrics <- df_long |>
  group_by(sample) |>
  summarize(
    pct_within_2_log2_mean = {
      m <- mean(count)
      l <- ifelse(m > 0, log2(m), 0)
      round(sum(abs(log2(count + 1) - l) <= 2) / n() * 100, 2)
    },
    .groups = "drop"
  )

# ---- Per-pool QC summary --------------------------------------------------
summary_df <- gini_df |>
  left_join(pct_ratio_df  |> select(sample, ratio_90_10), by = "sample") |>
  left_join(undetected_df |> select(sample, total_guides, zero_guides, pct_undetected), by = "sample") |>
  left_join(log2_metrics, by = "sample")
if (!is.null(mapping_df)) {
  summary_df <- summary_df |>
    left_join(mapping_df |> rename(sample = label) |> select(sample, reads, mappedreads, mapping_pct),
              by = "sample")
}
if (!is.null(classification_df)) {
  summary_df <- summary_df |>
    left_join(classification_df |> select(sample, near_miss_reads, pct_near_miss, far_miss_reads, pct_far_miss),
              by = "sample")
}
summary_df <- summary_df |>
  mutate(
    library_skew_flag = ratio_90_10 >= 10,
    undetected_flag   = pct_undetected >= 5
  )
save_table(summary_df, "qc_summary")

cat("\n=== QC summary:", pool_prefix, "===\n")
print(summary_df, row.names = FALSE)
cat("\nGuide representation tables complete.\n")
