#!/usr/bin/env Rscript
# norm_qc.R -- compute step of the normalisation QC: MA pairs, one reference line per normalisation method, numeric cloud summary.
# Reads bin counts from BIN_COUNTS (no BAMs). The plot is made by norm_qc_plot.R from the tables written here.
#
# CLI contract (called by the NORM_QC module; R style --key=value):
#   --counts=<bins_<size>.counts.tsv.gz>   bin_counts.R counts: chr, start, end, blacklisted, one integer column per sample id (ALL bins)
#   --libsizes=<bins_<size>.libsizes.tsv>  bin_counts.R libsizes: sample_id, csaw_total, bam_mapped_reads, bam_fragments, ...
#   --sample-ids=<comma list>              libraries of the norm_run, in run order (subset of the counts columns)
#   --samplesheet=<csv>                    project samplesheet (columns sample_id, antibody, condition, replicate, group are used)
#   --scale-factors-all=<tsv>              wide table <run>.scale_factors_all.tsv from norm_combine.py
#                                          (sample_id, primary_method, then <m>_raw_factor, _bigwig_multiplier, _edger_norm_factor,
#                                          _qc_pass, _status, _note for m in chipseqspikeinfree, tmm, cpm, spikein)
#   --primary=<method>                     primary method of the run (is_primary column)
#   --ref-condition=<condition or ''>      MA reference condition; '' = condition of the run's first library (warning)
#   --min-count=20                         bins with c_ref + c_other < min-count are dropped from a pair
#   --prior=0.5                            prior count of the MA log2 CPM, computed on the two columns of the pair
#   --out-prefix=<run>
# Outputs:
#   <run>.norm_qc_pairs.tsv     one row per pair x method (methods with status not_run in every library are omitted):
#                               pair_id, pair_type (replicate|treatment), ref, other, ref_group, other_group, method, is_primary,
#                               status, line, bigwig_offset, line_se, log2_factor_ratio_check, n_bins, min_count, prior, A_peak,
#                               M_peak, M_centroid, core_lo, core_hi, line_minus_centroid
#   <run>.norm_qc_libsizes.tsv  sample_id, csaw_total, bam_mapped_reads, bam_fragments, then <m>_edger_norm_factor, <m>_status per method
#   <run>.norm_qc.tsv           DEPRECATED, nothing new reads it: chr, start, end, one column per library with edgeR
#                               cpm(log = TRUE, prior.count = 2), lib.size = csaw_total, bins with run sum >= 10
# Exit codes: 0 ok; 2 usage error; 3 data problem (missing file/column, ref condition not in the run, line vs multiplier check fails).
#
# ---- What the lines mean ------------------------------------------------------------------------------------------------
# Pair (ref, other): x_i = log2((c_i + p_i)/(T_i + 2 p_i) * 1e6), p_i = prior * T_i / mean(T_ref, T_other);
#   M = x_other - x_ref (negative = less signal in `other`, the sign of the classification logFC), A = (x_other + x_ref)/2.
# Method m declares "no change" at line = log2(f_other / f_ref), f = edger_norm_factor (divisor of CPM); cpm is 0 exactly.
# It is a function of the factors only, not estimated from the cloud. log2_factor_ratio_check recomputes it from
# bigwig_multiplier, log2(mult_ref / mult_other), and the script exits 3 if the two differ by >= 1e-9.
# bigwig_offset = log2(B_other/T_other) - log2(B_ref/T_ref) (B = final-BAM fragments, T = csaw total): the bigWigs of method m
# agree at M = line + bigwig_offset. line_se (spikein only) = sqrt(2/S_ref + 2/S_other)/ln 2 with S = spike_reads parsed from the
# spikein_note of the wide table (NA when not parseable).
# Cloud centre: 2-D kernel density as drawn by smoothScatter (KernSmooth::bkde2D, bandwidth = 5-95 percent range / 25):
# M_peak = density peak, M_centroid = density-weighted centroid above half the peak, core = [min, max] of the two.
# CAVEAT: a line on the cloud is by construction for TMM and CPM and is not evidence that they are right when the mark changes globally.

options(stringsAsFactors = FALSE)

METHODS <- c("chipseqspikeinfree", "tmm", "cpm", "spikein")

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

# ---- pairing (pure function; spec 2.1) ----------------------------------------------------------------------------------------
# meta: data.frame in RUN ORDER with columns sample_id, group, condition, antibody, replicate.
# Returns data.frame(pair_id, pair_type, ref, other, ref_group, other_group): replicate-null pairs first (groups in order of first
# appearance, unordered pairs of one group, ref = earlier library), then treatment pairs (run order): each library whose condition is
# not ref_condition against the ref_condition library of the same antibody and replicate label, else the first ref_condition library
# of that antibody, else skipped with a warning.
make_pairs <- function(meta, ref_condition) {
  rows <- list()
  add <- function(type, r, o) {
    rows[[length(rows) + 1L]] <<- data.frame(
      pair_id = paste0(meta$sample_id[o], "_vs_", meta$sample_id[r]), pair_type = type,
      ref = meta$sample_id[r], other = meta$sample_id[o], ref_group = meta$group[r], other_group = meta$group[o])
  }
  for (g in unique(meta$group)) {
    ix <- which(meta$group == g)
    if (length(ix) < 2L) next
    for (a in seq_len(length(ix) - 1L)) for (b in (a + 1L):length(ix)) add("replicate", ix[a], ix[b])
  }
  for (i in seq_len(nrow(meta))) {
    if (meta$condition[i] == ref_condition) next
    cand <- which(meta$condition == ref_condition & meta$antibody == meta$antibody[i])
    if (!length(cand)) {
      cat("WARNING: no '", ref_condition, "' library of antibody '", meta$antibody[i], "' in the run; treatment pair for ",
          meta$sample_id[i], " skipped\n", sep = "", file = stderr())
      next
    }
    same <- cand[meta$replicate[cand] == meta$replicate[i]]
    add("treatment", if (length(same)) same[1] else cand[1], i)
  }
  if (length(rows)) do.call(rbind, rows) else
    data.frame(pair_id = character(0), pair_type = character(0), ref = character(0), other = character(0),
               ref_group = character(0), other_group = character(0))
}

# ---- MA coordinates of one pair (spec 2.2); c_ref, c_oth counts, T_ref, T_oth csaw totals -------------------------------------
ma_coords <- function(c_ref, c_oth, T_ref, T_oth, prior) {
  mT <- mean(c(T_ref, T_oth))
  p_r <- prior * T_ref / mT; p_o <- prior * T_oth / mT
  xr <- log2((c_ref + p_r) / (T_ref + 2 * p_r) * 1e6)
  xo <- log2((c_oth + p_o) / (T_oth + 2 * p_o) * 1e6)
  list(A = (xo + xr) / 2, M = xo - xr)
}

# ---- cloud centre (same density as smoothScatter) ------------------------------------------------------------------------------
cloud_center <- function(A, M, nbin = 256L) {
  na <- c(A_peak = NA_real_, M_peak = NA_real_, M_centroid = NA_real_)
  if (length(A) < 100L) return(na)
  X <- cbind(A, M)
  bw <- as.numeric(diff(apply(X, 2, quantile, probs = c(0.05, 0.95)))) / 25
  if (any(!is.finite(bw)) || any(bw <= 0)) return(na)
  d <- KernSmooth::bkde2D(X, bandwidth = bw, gridsize = c(nbin, nbin))
  f <- d$fhat
  ij <- which(f == max(f), arr.ind = TRUE)[1, ]
  Mg <- matrix(d$x2, nrow = nrow(f), ncol = ncol(f), byrow = TRUE)
  hi <- f >= 0.5 * max(f)
  c(A_peak = d$x1[ij[1]], M_peak = d$x2[ij[2]], M_centroid = sum(f[hi] * Mg[hi]) / sum(f[hi]))
}

# ---- status of a pair from the two libraries' status for a method --------------------------------------------------------------
pair_status <- function(s_ref, s_oth) {
  s <- c(s_ref, s_oth)
  if (any(s == "unavailable")) "unavailable"
  else if (any(s == "not_run")) "not_run"
  else if (any(s == "imputed")) "imputed"
  else if (any(s == "qc_fail")) "qc_fail"
  else "ok"
}

spike_reads <- function(note) {
  if (is.na(note)) return(NA_real_)
  m <- regmatches(note, regexec("spike_reads=([0-9]+)", note))[[1]]
  if (length(m) == 2L) as.numeric(m[2]) else NA_real_
}

main <- function() {
  allowed <- c("counts", "libsizes", "sample-ids", "samplesheet", "scale-factors-all", "primary", "ref-condition", "min-count",
               "prior", "out-prefix")
  args <- parse_args(commandArgs(trailingOnly = TRUE), allowed)
  counts_f   <- get_req(args, "counts")
  libsize_f  <- get_req(args, "libsizes")
  ids        <- split_list(get_req(args, "sample-ids"))
  ss_f       <- get_req(args, "samplesheet")
  wide_f     <- get_req(args, "scale-factors-all")
  primary    <- get_req(args, "primary")
  ref_cond   <- trimws(get_opt(args, "ref-condition", ""))
  min_count  <- as_num(get_opt(args, "min-count", "20"), "min-count")
  prior      <- as_num(get_opt(args, "prior", "0.5"), "prior")
  out_prefix <- get_req(args, "out-prefix")
  if (length(ids) < 2L) die("need at least 2 libraries")
  if (anyDuplicated(ids)) die("duplicated sample ids")
  for (f in c(counts_f, libsize_f, ss_f, wide_f)) if (!file.exists(f)) die(paste0("file not found: ", f))

  suppressPackageStartupMessages({ library(edgeR); library(KernSmooth) })

  # samplesheet -> per-library metadata in run order
  ss <- read.csv(ss_f, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character", na.strings = character(0))
  colnames(ss)[1] <- sub("^﻿", "", colnames(ss)[1])
  need_ss <- c("sample_id", "antibody", "condition", "replicate", "group")
  if (!all(need_ss %in% colnames(ss))) die(paste0("samplesheet lacks columns: ", paste(setdiff(need_ss, colnames(ss)), collapse = ", ")))
  ix <- match(ids, ss$sample_id)
  if (anyNA(ix)) die(paste0("sample ids not in the samplesheet: ", paste(ids[is.na(ix)], collapse = ", ")))
  meta <- ss[ix, need_ss, drop = FALSE]
  if (!nzchar(ref_cond)) {
    ref_cond <- meta$condition[1]
    cat("WARNING: norm_run ", out_prefix, ": MA reference condition = ", ref_cond, " (first library of the run; no --ref-condition)\n",
        sep = "", file = stderr())
  } else if (!(ref_cond %in% meta$condition)) {
    die(paste0("reference condition '", ref_cond, "' is not a condition of the run's libraries"))
  }

  # counts and library sizes
  cn <- read.delim(counts_f, check.names = FALSE, stringsAsFactors = FALSE, colClasses = c(chr = "character"))
  need_cn <- c("chr", "start", "end", "blacklisted", ids)
  if (!all(need_cn %in% colnames(cn))) die(paste0("counts table lacks columns: ", paste(setdiff(need_cn, colnames(cn)), collapse = ", ")))
  ls0 <- read.delim(libsize_f, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("sample_id", "csaw_total", "bam_mapped_reads", "bam_fragments") %in% colnames(ls0))) die("libsizes table lacks columns")
  li <- match(ids, ls0$sample_id)
  if (anyNA(li)) die(paste0("sample ids not in the libsizes table: ", paste(ids[is.na(li)], collapse = ", ")))
  Tt <- setNames(as.numeric(ls0$csaw_total[li]), ids)
  Bf <- setNames(as.numeric(ls0$bam_fragments[li]), ids)
  cm <- as.matrix(cn[, ids, drop = FALSE])
  blk <- as.logical(cn$blacklisted)
  if (anyNA(blk)) blk[is.na(blk)] <- FALSE

  # wide scale-factor table
  w <- read.delim(wide_f, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character", na.strings = c("NA", ""))
  if (!("sample_id" %in% colnames(w))) die("scale_factors_all table lacks sample_id")
  wi <- match(ids, w$sample_id)
  if (anyNA(wi)) die(paste0("sample ids not in the scale_factors_all table: ", paste(ids[is.na(wi)], collapse = ", ")))
  getcol <- function(m, field) if (paste0(m, "_", field) %in% colnames(w)) w[[paste0(m, "_", field)]][wi] else rep(NA_character_, length(ids))
  fac <- list()
  for (m in METHODS) {
    st <- getcol(m, "status"); st[is.na(st)] <- "not_run"
    fac[[m]] <- list(status = setNames(st, ids),
                     nf = setNames(suppressWarnings(as.numeric(getcol(m, "edger_norm_factor"))), ids),
                     mult = setNames(suppressWarnings(as.numeric(getcol(m, "bigwig_multiplier"))), ids),
                     note = setNames(getcol(m, "note"), ids))
  }
  run_methods <- METHODS[vapply(fac, function(z) any(z$status != "not_run"), logical(1))]
  if (!length(run_methods)) die("no method with a status other than not_run in the scale_factors_all table")

  # deprecated binned logCPM table (prior 2, run sum >= 10)
  keep10 <- rowSums(cm) >= 10
  if (sum(keep10) >= 1L) {
    lcpm <- edgeR::cpm(cm[keep10, , drop = FALSE], lib.size = Tt, log = TRUE, prior.count = 2)
    tab <- data.frame(chr = cn$chr[keep10], start = cn$start[keep10], end = cn$end[keep10], round(lcpm, 4), check.names = FALSE)
  } else {
    tab <- data.frame(chr = character(0), start = numeric(0), end = numeric(0))
  }
  write.table(tab, paste0(out_prefix, ".norm_qc.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  # pairs and lines
  prs <- make_pairs(meta, ref_cond)
  if (!nrow(prs)) die("no MA pairs (no replicate group with 2 libraries and no treated library with a reference library)")
  rows <- list(); max_off <- 0
  for (k in seq_len(nrow(prs))) {
    r <- prs$ref[k]; o <- prs$other[k]
    sel <- !blk & (cm[, r] + cm[, o] >= min_count)
    ma <- ma_coords(cm[sel, r], cm[sel, o], Tt[[r]], Tt[[o]], prior)
    cc <- cloud_center(ma$A, ma$M)
    lo <- suppressWarnings(min(cc[c("M_peak", "M_centroid")])); hi <- suppressWarnings(max(cc[c("M_peak", "M_centroid")]))
    if (!is.finite(lo)) { lo <- NA_real_; hi <- NA_real_ }
    off <- log2(Bf[[o]] / Tt[[o]]) - log2(Bf[[r]] / Tt[[r]])
    if (is.finite(off)) max_off <- max(max_off, abs(off))
    for (m in run_methods) {
      z <- fac[[m]]
      st <- pair_status(z$status[[r]], z$status[[o]])
      fr <- z$nf[[r]]; fo <- z$nf[[o]]
      line <- if (st %in% c("unavailable", "not_run") || !all(is.finite(c(fr, fo))) || fr <= 0 || fo <= 0) NA_real_ else log2(fo / fr)
      chk <- if (is.finite(line) && all(is.finite(c(z$mult[[r]], z$mult[[o]]))) && z$mult[[r]] > 0 && z$mult[[o]] > 0)
        log2(z$mult[[r]] / z$mult[[o]]) else NA_real_
      if (is.finite(chk) && abs(chk - line) >= 1e-9)
        die(sprintf("method %s, pair %s vs %s: line %.12g from edger_norm_factor differs from %.12g from bigwig_multiplier (|diff| >= 1e-9); multiplier x factor is not constant within the run",
                    m, o, r, line, chk))
      se <- NA_real_
      if (m == "spikein" && is.finite(line)) {
        S <- c(spike_reads(z$note[[r]]), spike_reads(z$note[[o]]))
        if (all(is.finite(S)) && all(S > 0)) se <- sqrt(2 / S[1] + 2 / S[2]) / log(2)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        pair_id = prs$pair_id[k], pair_type = prs$pair_type[k], ref = r, other = o, ref_group = prs$ref_group[k],
        other_group = prs$other_group[k], method = m, is_primary = (m == primary), status = st, line = line,
        bigwig_offset = off, line_se = se, log2_factor_ratio_check = chk, n_bins = sum(sel), min_count = min_count, prior = prior,
        A_peak = unname(cc["A_peak"]), M_peak = unname(cc["M_peak"]), M_centroid = unname(cc["M_centroid"]),
        core_lo = lo, core_hi = hi, line_minus_centroid = line - unname(cc["M_centroid"]))
    }
  }
  pairs <- do.call(rbind, rows)
  if (!(primary %in% run_methods)) cat("WARNING: primary method '", primary, "' has no data in the scale_factors_all table\n", sep = "", file = stderr())
  if (max_off > 0.01)
    cat(sprintf("WARNING: bigWigs and edgeR steps differ by up to %.4f log2 in run %s (bigwig_offset): fragments longer than max.frag or orphan reads\n",
                max_off, out_prefix), file = stderr())
  write.table(pairs, paste0(out_prefix, ".norm_qc_pairs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

  lib <- data.frame(sample_id = ids, csaw_total = Tt, bam_mapped_reads = as.numeric(ls0$bam_mapped_reads[li]), bam_fragments = Bf)
  for (m in METHODS) {
    lib[[paste0(m, "_edger_norm_factor")]] <- unname(fac[[m]]$nf)
    lib[[paste0(m, "_status")]] <- unname(fac[[m]]$status)
  }
  write.table(lib, paste0(out_prefix, ".norm_qc_libsizes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  cat("Done: ", out_prefix, ".norm_qc_pairs.tsv (", nrow(prs), " pairs, ", length(run_methods), " methods, ", nrow(pairs), " rows)\n", sep = "")
}

if (sys.nframe() == 0L) main()   # not run when the file is source()d
