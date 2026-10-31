#!/usr/bin/env Rscript
# qc_replicate_correlation.R -- replicate concordance of the final BAMs, from the shared csaw bin counts (no BAM is read).
#
# Called by REPLICATE_CORRELATION (modules/local/replicate_correlation.nf). Descriptive values only: no pass/fail threshold.
#
# CLI (--key=value):
#   --counts=<bins_<size>.counts.tsv.gz>  bin_counts.R output: tab separated, header chr, start, end, blacklisted, then one
#                                         integer count column per sample_id; ALL bins (filter 0)
#   --libsizes=<bins_<size>.libsizes.tsv> bin_counts.R output; columns sample_id, csaw_total (= edgeR lib.size), ...
#   --samplesheet=<samplesheet.csv>       columns sample_id, antibody, condition, group are read
#   --sample-ids=<id1,id2,...>            active libraries (default: every count column)
#   --prior=2                             prior count of the log2 CPM used by the Pearson statistics
#   --signal-frac=0.10                    fraction of "all" bins (highest mean log2 CPM of the pair) forming the signal universe
#   --bin-size=<bp>                       reported in the table (default: bin_size column of --libsizes, else median bin width)
#   --out-prefix=replicate_correlation    output file prefix
#   --outdir=.                            output directory
# Outputs:
#   <prefix>.tsv             one row per pair of libraries of the same samplesheet group (groups with one library: one row,
#                            status single_library, NA metrics). Columns: group, antibody, condition, sample_1, sample_2, status,
#                            n_bins_all, pearson_all, spearman_all, n_bins_signal, median_pair_count_signal, pearson_signal,
#                            robust_sd_log2ratio_signal, nonrep_pearson_all_median, nonrep_n, bin_size, prior, signal_frac
#   <prefix>_pairs_all.tsv   every same-antibody pair of active libraries: sample_1, sample_2, antibody, same_group,
#                            pearson_all, spearman_all, pearson_signal (the last column feeds the signal row of the plot)
# Definitions (per pair, library sizes T = csaw_total; blacklisted bins dropped):
#   universe "all"   = bins with c1 + c2 > 0
#   x_i = log2((c_i + p_i) / (T_i + 2 p_i) * 1e6), p_i = prior * T_i / mean(T_1, T_2)   (edgeR cpm, log = TRUE, on the two columns)
#   pearson_all = cor(x_1, x_2); spearman_all = Spearman of the raw counts (average ranks for ties)
#   universe "signal" = the top signal_frac of the "all" bins by A = (x_1 + x_2) / 2, ties at the cut-off included
#   robust_sd_log2ratio_signal = 1.4826 * median(|M - median(M)|), M = x_2 - x_1 on the signal bins (median-centred)
#   nonrep_pearson_all_median = median pearson_all between either library of the pair and every library of the same antibody
#                               in OTHER groups; nonrep_n = number of those pairs
# status: ok | single_library | too_few_bins (fewer than 3 bins in the universe; metrics NA)
# A library without a group value in the samplesheet is treated as its own group.
# Exit codes: 0 ok; 2 usage; 3 data problem.

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

# log2 CPM of two libraries, prior scaled by library size (edgeR cpm(log = TRUE) on the two columns)
pair_x <- function(c1, c2, T1, T2, prior) {
  pm <- mean(c(T1, T2))
  p1 <- prior * T1 / pm
  p2 <- prior * T2 / pm
  list(x1 = log2((c1 + p1) / (T1 + 2 * p1) * 1e6), x2 = log2((c2 + p2) / (T2 + 2 * p2) * 1e6))
}
safe_cor <- function(a, b, method = "pearson") {
  if (length(a) < 3L) return(NA_real_)
  suppressWarnings(cor(a, b, method = method))
}
# all metrics of one pair; counts are the NON-blacklisted bins
pair_metrics <- function(c1, c2, T1, T2, prior, frac) {
  keep <- (c1 + c2) > 0
  c1 <- c1[keep]; c2 <- c2[keep]
  n <- length(c1)
  out <- list(n_bins_all = n, pearson_all = NA_real_, spearman_all = NA_real_, n_bins_signal = 0L,
              median_pair_count_signal = NA_real_, pearson_signal = NA_real_, robust_sd = NA_real_)
  if (n < 3L) return(out)
  xx <- pair_x(c1, c2, T1, T2, prior)
  out$pearson_all <- safe_cor(xx$x1, xx$x2)
  out$spearman_all <- safe_cor(c1, c2, "spearman")
  A <- (xx$x1 + xx$x2) / 2
  k <- min(n, max(1L, ceiling(frac * n)))
  thr <- sort(A, decreasing = TRUE)[k]
  sel <- A >= thr
  out$n_bins_signal <- sum(sel)
  out$median_pair_count_signal <- median(c1[sel] + c2[sel])
  out$pearson_signal <- safe_cor(xx$x1[sel], xx$x2[sel])
  M <- xx$x2[sel] - xx$x1[sel]
  out$robust_sd <- 1.4826 * median(abs(M - median(M)))
  out
}

allowed <- c("counts", "libsizes", "samplesheet", "sample-ids", "prior", "signal-frac", "bin-size", "out-prefix", "outdir")
args <- parse_args(commandArgs(trailingOnly = TRUE), allowed)
counts_file <- get_req(args, "counts")
lib_file    <- get_req(args, "libsizes")
ss_file     <- get_req(args, "samplesheet")
prior       <- as.numeric(get_opt(args, "prior", "2"))
frac        <- as.numeric(get_opt(args, "signal-frac", "0.10"))
out_prefix  <- get_opt(args, "out-prefix", "replicate_correlation")
outdir      <- get_opt(args, "outdir", ".")
if (!is.finite(prior) || prior < 0) die("--prior must be a number >= 0", 2L)
if (!is.finite(frac) || frac <= 0 || frac > 1) die("--signal-frac must be in (0, 1]", 2L)
for (f in c(counts_file, lib_file, ss_file)) if (!file.exists(f)) die(paste0("file not found: ", f))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cnt <- read.delim(counts_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("chr", "start", "end", "blacklisted") %in% colnames(cnt))) die("counts lack chr, start, end, blacklisted columns")
all_ids <- colnames(cnt)[-(1:4)]
ids <- if (is.null(args[["sample-ids"]]) || !nzchar(args[["sample-ids"]])) all_ids else split_list(args[["sample-ids"]])
if (length(ids) < 1L) die("no sample ids")
miss <- setdiff(ids, all_ids)
if (length(miss)) die(paste0("sample ids not in the counts file: ", paste(miss, collapse = ", ")))
lib <- read.delim(lib_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("sample_id", "csaw_total") %in% colnames(lib))) die("libsizes lack sample_id, csaw_total")
Tl <- setNames(as.numeric(lib$csaw_total[match(ids, lib$sample_id)]), ids)
if (anyNA(Tl) || any(Tl <= 0)) die("csaw_total missing or not positive for some sample ids")
bin_size <- get_opt(args, "bin-size", NA)
if (is.na(bin_size)) {
  bin_size <- if ("bin_size" %in% colnames(lib)) lib$bin_size[1] else median(cnt$end - cnt$start)
}
bin_size <- as.numeric(bin_size)

keep_bins <- !(as.logical(cnt$blacklisted))
if (anyNA(keep_bins)) die("blacklisted column has missing values")
C <- as.matrix(cnt[keep_bins, ids, drop = FALSE])
storage.mode(C) <- "double"
rm(cnt)
if (nrow(C) == 0L) die("no non-blacklisted bins")

ss <- read.csv(ss_file, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM", na.strings = character(0))
need <- c("sample_id", "antibody", "condition", "group")
if (!all(need %in% colnames(ss))) die(paste0("samplesheet lacks columns: ", paste(setdiff(need, colnames(ss)), collapse = ", ")))
ss <- ss[match(ids, ss$sample_id), need, drop = FALSE]
if (anyNA(ss$sample_id)) die(paste0("sample ids not in the samplesheet: ", paste(ids[is.na(ss$sample_id)], collapse = ", ")))
ss$group <- ifelse(nzchar(trimws(ss$group)), ss$group, ss$sample_id)

# ---- every same-antibody pair -------------------------------------------------
pr <- list(); met <- list(); k <- 0L
for (i in seq_len(length(ids) - 1L)) for (j in (i + 1L):length(ids)) {
  if (ss$antibody[i] != ss$antibody[j]) next
  m <- pair_metrics(C[, i], C[, j], Tl[i], Tl[j], prior, frac)
  k <- k + 1L
  met[[k]] <- m
  pr[[k]] <- data.frame(sample_1 = ids[i], sample_2 = ids[j], antibody = ss$antibody[i],
                        same_group = ss$group[i] == ss$group[j], pearson_all = m$pearson_all,
                        spearman_all = m$spearman_all, pearson_signal = m$pearson_signal)
}
pairs_all <- if (k > 0L) do.call(rbind, pr) else
  data.frame(sample_1 = character(0), sample_2 = character(0), antibody = character(0), same_group = logical(0),
             pearson_all = numeric(0), spearman_all = numeric(0), pearson_signal = numeric(0))

# ---- table: one row per replicate pair, or single_library --------------------
rows <- list()
for (g in unique(ss$group)) {
  mem <- which(ss$group == g)
  ab <- ss$antibody[mem[1]]; cond <- ss$condition[mem[1]]
  if (length(mem) == 1L) {
    rows[[length(rows) + 1L]] <- data.frame(group = g, antibody = ab, condition = cond, sample_1 = ids[mem], sample_2 = NA_character_,
                                            status = "single_library", n_bins_all = NA_integer_, pearson_all = NA_real_,
                                            spearman_all = NA_real_, n_bins_signal = NA_integer_, median_pair_count_signal = NA_real_,
                                            pearson_signal = NA_real_, robust_sd_log2ratio_signal = NA_real_,
                                            nonrep_pearson_all_median = NA_real_, nonrep_n = NA_integer_)
    next
  }
  for (a in seq_len(length(mem) - 1L)) for (b in (a + 1L):length(mem)) {
    s1 <- ids[mem[a]]; s2 <- ids[mem[b]]
    idx <- which(pairs_all$sample_1 == s1 & pairs_all$sample_2 == s2)
    if (length(idx) != 1L) die(paste0("group ", g, " mixes antibodies; pair ", s1, "/", s2, " is not a same-antibody pair"))
    m <- met[[idx]]
    other <- !pairs_all$same_group & pairs_all$antibody == ab &
      (pairs_all$sample_1 %in% c(s1, s2) | pairs_all$sample_2 %in% c(s1, s2))
    nr <- pairs_all$pearson_all[other]
    rows[[length(rows) + 1L]] <- data.frame(group = g, antibody = ab, condition = cond, sample_1 = s1, sample_2 = s2,
                                            status = if (m$n_bins_all >= 3L) "ok" else "too_few_bins",
                                            n_bins_all = m$n_bins_all, pearson_all = m$pearson_all, spearman_all = m$spearman_all,
                                            n_bins_signal = m$n_bins_signal, median_pair_count_signal = m$median_pair_count_signal,
                                            pearson_signal = m$pearson_signal, robust_sd_log2ratio_signal = m$robust_sd,
                                            nonrep_pearson_all_median = if (any(!is.na(nr))) median(nr, na.rm = TRUE) else NA_real_,
                                            nonrep_n = sum(!is.na(nr)))
  }
}
tab <- do.call(rbind, rows)
tab$bin_size <- bin_size
tab$prior <- prior
tab$signal_frac <- frac

write.table(tab, file.path(outdir, paste0(out_prefix, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(pairs_all, file.path(outdir, paste0(out_prefix, "_pairs_all.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
cat(sprintf("replicate concordance: %d libraries, %d non-blacklisted bins, %d table rows (%d single_library), %d same-antibody pairs\n",
            length(ids), nrow(C), nrow(tab), sum(tab$status == "single_library"), nrow(pairs_all)))
