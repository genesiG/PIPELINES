#!/usr/bin/env Rscript
# peak_annotation.R -- annotate one or more peak sets against a GTF (ChIPseeker) and write per-region metrics,
# per-label summaries, a pairwise overlap QC table and annotation plots.
# Port of the legacy utils/peak_annotation.R (ChIPseeker/clusterProfiler annotation engine). The legacy driver
# (bsub submission script) is replaced by the Nextflow PEAK_ANNOTATION module; all settings come from the CLI.
#
# CLI contract (called by the PEAK_ANNOTATION module; --key=value):
#   --bed-files=a.bed,b.bed        REQUIRED. Comma list of BED files, one per label. Tab-separated; >= 3 columns
#                                  (chr, start, end; 0-based half-open). Column 4, if present, is read as 'name'.
#                                  Lines starting with '#', 'track' or 'browser' are skipped. A BED with 0 rows
#                                  is dropped from all analyses with a warning; exit 3 if all are empty.
#   --labels=l1,l2                 REQUIRED. Comma list, same count as --bed-files, unique, no commas.
#   --gtf=F                        REQUIRED. Gene annotation GTF (gene_id = gene symbol). The TxDb is built on
#                                  every run with makeTxDbFromGFF(format = "gtf"); no TxDb cache file is written
#                                  (the Nextflow task cache covers reruns of the whole task).
#   --out-prefix=NAME              REQUIRED. Outputs go to the CURRENT directory as NAME.<name> (see below).
#   --colors="l1=#hex,l2=#hex"     OPTIONAL. Colours for labels. Labels without a colour get colours from the
#                                  RColorBrewer Set2 then Dark2 palettes, assigned in label order.
#   --tss-region=3000              OPTIONAL. bp on each side of the TSS: ChIPseeker tssRegion = c(-r, r) and
#                                  getPromoters(upstream = r, downstream = r). Default 3000.
#   --kegg=true|false              OPTIONAL. Default true. false (or a failed enrichment) writes a placeholder page.
#   --cpus=1                       OPTIONAL. Sets options(mc.cores); this script has no parallel loop.
#
# Outputs (cwd, prefix NAME):
#   NAME.region_metrics.tsv              one row per region: label, seqnames, start, end, width_bp, gc_content,
#                                        cpg_oe, annotation, dist_to_tss_bp, nearest_gene_id, gene_biotype
#   NAME.summary_stats.tsv               one row per label: counts, width statistics, GC, CpG o/e, annotation %, TSS distance
#   NAME.overlap_matrix.tsv              all ordered label pairs: set_A, set_B, n_overlapping_pairs, jaccard_bp
#   NAME.tss_windows.bed                 TSS +/- tss-region windows (0-based BED) from getPromoters
#   NAME.annotation_pie_<slug>.svg       one per label; slug = label with [^A-Za-z0-9_] replaced by _
#   NAME.annotation_bar.svg              ChIPseeker annotation bars, all labels
#   NAME.dist_to_tss.svg                 ChIPseeker distance-to-TSS bars, all labels
#   NAME.width_distribution.svg          peak width violin + box (log10 scale)
#   NAME.width_ecdf.svg                  peak width ECDF with density panel
#   NAME.gc_cpgoe.svg                    GC fraction and CpG o/e distributions (message page if sequences unavailable)
#   NAME.chromosome_distribution.svg     peak count per chromosome (chr1-22, X, Y)
#   NAME.tss_dist_ecdf.svg               distance to nearest TSS, cumulative distribution
#   NAME.chromosome_coverage.svg         total bp per chromosome (chr1-22, X, Y)
#   NAME.chromosome_coverage_<slug>.svg  one per label, chromosomes ordered high to low
#   NAME.kegg.svg                        compareCluster KEGG dotplot, or a placeholder page giving the reason
#
# Legacy plot numbering (Scripts/utils/peak_annotation.R) -> new files:
#   03 -> annotation_pie_<slug>; 04 -> annotation_bar; 05 -> dist_to_tss; 06 -> width_distribution;
#   07 -> width_ecdf (title no longer names a mark); 08 -> gc_cpgoe; 09 -> chromosome_distribution;
#   10 -> tss_dist_ecdf; 11 -> chromosome_coverage; 12 -> chromosome_coverage_<slug>; 13 -> kegg.
#   01a/01b (EnrichedHeatmap TSS / gene-body) and 02 (average profile) are NOT produced here: heatmaps are
#   made by the heatmap module, and the legacy heatmaps needed bigWigs.
#
# Genome: BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0 only. GC/CpG metrics are NA with a warning when the package
# or the sequences are not found. Chromosome names are matched to the GTF/BSgenome, converting seqlevelsStyle
# when none match directly; a label with no matching chromosome names is an error (exit 3).
#
# Exit codes: 0 ok; 2 usage error; 3 data or environment problem (missing BED/GTF, no usable peaks, chromosome
# names not matching the GTF, tables not written). A failing plot does not stop the run: it is written as a
# message page and the run continues.

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
slug <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

allowed <- c("bed-files", "labels", "gtf", "out-prefix", "colors", "tss-region", "kegg", "cpus")
args <- parse_args(commandArgs(trailingOnly = TRUE), allowed)
bed_files  <- split_list(get_req(args, "bed-files"))
labels     <- split_list(get_req(args, "labels"))
gtf_file   <- get_req(args, "gtf")
out_prefix <- get_req(args, "out-prefix")
colors_arg <- get_opt(args, "colors", "")
tss_r      <- suppressWarnings(as.integer(get_opt(args, "tss-region", "3000")))
kegg_arg   <- tolower(get_opt(args, "kegg", "true"))
cpus       <- suppressWarnings(as.integer(get_opt(args, "cpus", "1")))

if (length(bed_files) == 0L || length(labels) == 0L) die("--bed-files and --labels must not be empty", 2L)
if (length(bed_files) != length(labels)) {
  die(sprintf("--bed-files (%d) and --labels (%d) must have the same length",
              length(bed_files), length(labels)), 2L)
}
if (anyDuplicated(labels)) die("--labels must be unique", 2L)
if (anyDuplicated(slug(labels))) die("--labels must stay unique after replacing characters outside [A-Za-z0-9_] by _", 2L)
if (is.na(tss_r) || tss_r < 1L) die("--tss-region must be a positive integer (bp)", 2L)
if (!(kegg_arg %in% c("true", "false"))) die("--kegg must be true or false", 2L)
kegg_on <- kegg_arg == "true"
if (is.na(cpus) || cpus < 1L) die("--cpus must be a positive integer", 2L)

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicFeatures)
  library(Biostrings)
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(gridExtra)
  library(RColorBrewer)
})
ok_pkg <- function(p) requireNamespace(p, quietly = TRUE)
if (!ok_pkg("ChIPseeker")) die("ChIPseeker is not installed in this image", 3L)
suppressPackageStartupMessages(library(ChIPseeker))
CLUSTERP_AVAILABLE <- ok_pkg("clusterProfiler") && ok_pkg("org.Hs.eg.db")
if (CLUSTERP_AVAILABLE) {
  suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
  })
}
options(mc.cores = cpus)

cat("=== peak_annotation.R ===\n")
cat(sprintf("  Labels    : %s\n", paste(labels, collapse = ", ")))
cat(sprintf("  Output    : %s.*  (current directory)\n\n", out_prefix))

# =============================================================================
# 1. Colours: user colours first, then Set2 / Dark2 in label order
# =============================================================================
user_cols <- character(0)
for (p in split_list(colors_arg)) {
  kv <- regmatches(p, regexec("^([^=]+)=(.+)$", p))[[1]]
  if (length(kv) != 3L) die(paste0("cannot parse --colors entry '", p, "' (expected label=#hex)"), 2L)
  user_cols[trimws(kv[2])] <- trimws(kv[3])
}
if (!all(grepl("^#[0-9A-Fa-f]{6}$", user_cols))) die("--colors values must be #RRGGBB", 2L)
unknown_cols <- setdiff(names(user_cols), labels)
if (length(unknown_cols) > 0L) {
  cat("WARNING: --colors given for labels not in --labels (ignored): ",
      paste(unknown_cols, collapse = ", "), "\n", sep = "", file = stderr())
}

colour_pool <- c(brewer.pal(8, "Set2"), brewer.pal(8, "Dark2"))
COLS <- character(0)
n_auto <- 0L
for (lbl in labels) {
  if (lbl %in% names(user_cols)) {
    COLS[lbl] <- user_cols[[lbl]]
  } else {
    COLS[lbl] <- colour_pool[(n_auto %% length(colour_pool)) + 1L]
    n_auto <- n_auto + 1L
  }
}

# Purple ramp for categorical slices (pie charts); fixed, no palette package needed
purple_ramp <- function(n) colorRampPalette(c("#efedf5", "#bcbddc", "#807dba", "#54278f"))(n)

# =============================================================================
# 2. Read BED files as GRanges (>= 3 columns; empty BEDs dropped)
# =============================================================================
cat("=== Reading BED files ===\n")

read_bed3 <- function(path, label) {
  if (!file.exists(path)) die(paste0("BED file not found: ", path), 3L)
  lines <- readLines(path, warn = FALSE)
  lines <- lines[nzchar(trimws(lines)) & !grepl("^(#|track|browser)", lines)]
  if (length(lines) == 0L) return(GRanges())
  f <- strsplit(lines, "\t", fixed = TRUE)
  if (!all(lengths(f) >= 3L)) die(paste0("BED file has fewer than 3 tab-separated columns: ", path), 3L)
  chr <- vapply(f, `[`, "", 1L)
  st  <- suppressWarnings(as.numeric(vapply(f, `[`, "", 2L)))
  en  <- suppressWarnings(as.numeric(vapply(f, `[`, "", 3L)))
  if (anyNA(st) || anyNA(en)) die(paste0("non-numeric start/end in: ", path), 3L)
  if (any(en <= st)) die(paste0("BED interval with end <= start in: ", path), 3L)
  gr <- GRanges(seqnames = chr, ranges = IRanges(start = st + 1, end = en))  # BED 0-based -> 1-based
  if (all(lengths(f) >= 4L)) gr$name <- vapply(f, `[`, "", 4L)
  gr$label <- label
  cat(sprintf("  [%s] %d regions (%.2f Mb)\n", label, length(gr), sum(width(gr)) / 1e6))
  gr
}

gr_all <- lapply(seq_along(bed_files), function(i) read_bed3(bed_files[i], labels[i]))
names(gr_all) <- labels
n_regions <- vapply(gr_all, length, 0L)
for (lbl in labels[n_regions == 0L]) {
  cat(sprintf("WARNING: no regions in BED for label '%s'; dropped from all analyses\n", lbl), file = stderr())
}
labels <- labels[n_regions > 0L]
if (length(labels) == 0L) die("no BED file contains any region", 3L)
gr_list <- gr_all[labels]

# =============================================================================
# 3. TxDb from the GTF (built per run; no cache file)
# =============================================================================
cat("\n=== Building TxDb from GTF ===\n")
cat(sprintf("  GTF: %s\n", gtf_file))
if (!file.exists(gtf_file)) die(paste0("GTF file not found: ", gtf_file), 3L)
txdb <- tryCatch(
  suppressWarnings(makeTxDbFromGFF(gtf_file, format = "gtf", taxonomyId = 9606L)),
  error = function(e) die(paste0("could not build TxDb from GTF: ", conditionMessage(e)), 3L)
)
cat("  TxDb built.\n")

# Chromosome names: use as is when they match the GTF, else convert seqlevelsStyle, else stop
txdb_seqs <- tryCatch(seqlevels(txdb), error = function(e) character(0))
for (lbl in labels) {
  g <- gr_list[[lbl]]
  if (length(intersect(seqlevels(g), txdb_seqs)) == 0L) {
    g2 <- g
    tryCatch(seqlevelsStyle(g2) <- seqlevelsStyle(txdb)[1], error = function(e) NULL)
    if (length(intersect(seqlevels(g2), txdb_seqs)) == 0L) {
      die(sprintf("no chromosome name of label '%s' (e.g. %s) matches the GTF",
                  lbl, paste(head(seqlevels(g), 3), collapse = ", ")), 3L)
    }
    cat(sprintf("  [%s] chromosome names converted to GTF style\n", lbl))
    gr_list[[lbl]] <- g2
  }
}

# =============================================================================
# 4. BSgenome sequence metrics (GC, CpG o/e)
# =============================================================================
cat("\n=== BSgenome: GC content and CpG O/E ===\n")
BSG_PKG <- "BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0"
bsg <- NULL
if (ok_pkg(BSG_PKG)) {
  suppressPackageStartupMessages(library(BSG_PKG, character.only = TRUE))
  bsg <- get(BSG_PKG)
  cat(sprintf("  Using BSgenome package: %s\n", BSG_PKG))
} else {
  cat(sprintf("WARNING: %s not installed; GC/CpG metrics will be NA\n", BSG_PKG), file = stderr())
}

std_chrs <- paste0("chr", c(1:22, "X", "Y", "M"))

compute_seq_metrics <- function(gr) {
  res <- data.frame(gc = rep(NA_real_, length(gr)), cpgoe = rep(NA_real_, length(gr)))
  if (is.null(bsg)) return(res)
  g <- gr
  bsg_seqs <- seqnames(bsg)
  if (length(intersect(seqlevels(g), bsg_seqs)) == 0L) {
    tryCatch(seqlevelsStyle(g) <- seqlevelsStyle(bsg)[1], error = function(e) NULL)
  }
  idx <- which(as.character(seqnames(g)) %in% intersect(bsg_seqs, std_chrs))
  if (length(idx) == 0L) {
    cat("WARNING: no peak chromosome found in BSgenome; GC/CpG metrics will be NA\n", file = stderr())
    return(res)
  }
  sub_gr <- g[idx]
  sub_gr <- keepSeqlevels(sub_gr, unique(as.character(seqnames(sub_gr))), pruning.mode = "coarse")
  seqs <- tryCatch(getSeq(bsg, sub_gr), error = function(e) {
    cat("WARNING: getSeq failed; GC/CpG metrics will be NA: ", conditionMessage(e), "\n", sep = "", file = stderr())
    NULL
  })
  if (is.null(seqs)) return(res)
  bases <- letterFrequency(seqs, letters = c("C", "G"), as.prob = FALSE)
  cpg_f <- as.numeric(dinucleotideFrequency(seqs, as.prob = FALSE)[, "CG"])
  len   <- as.numeric(width(sub_gr))
  nC    <- as.numeric(bases[, "C"])
  nG    <- as.numeric(bases[, "G"])
  res$gc[idx]    <- as.numeric(rowSums(bases)) / len
  res$cpgoe[idx] <- ifelse(nC > 0 & nG > 0, (cpg_f * len) / (nC * nG), NA_real_)
  res
}

seq_metrics_list <- lapply(gr_list, compute_seq_metrics)

# =============================================================================
# 5. ChIPseeker annotation
# =============================================================================
cat("\n=== ChIPseeker annotation ===\n")
anno_list <- lapply(labels, function(lbl) {
  cat(sprintf("  Annotating [%s] ...\n", lbl))
  tryCatch(
    annotatePeak(gr_list[[lbl]], tssRegion = c(-tss_r, tss_r), TxDb = txdb, verbose = FALSE),
    error = function(e) die(paste0("annotatePeak failed for '", lbl, "': ", conditionMessage(e)), 3L)
  )
})
names(anno_list) <- labels
anno_df_list <- lapply(anno_list, as.data.frame)

# =============================================================================
# 6. Tables: region_metrics.tsv, summary_stats.tsv, overlap_matrix.tsv
# =============================================================================
cat("\n=== Assembling tables ===\n")

anno_fractions <- function(ann) {
  list(
    pct_promoter   = mean(grepl("Promoter",          ann)) * 100,
    pct_exon       = mean(grepl("Exon",              ann)) * 100,
    pct_intron     = mean(grepl("Intron",            ann)) * 100,
    pct_downstream = mean(grepl("Downstream",        ann)) * 100,
    pct_intergenic = mean(grepl("Intergenic|Distal", ann)) * 100
  )
}

jaccard_bp <- function(A, B) {
  inter_bp <- sum(width(GenomicRanges::intersect(A, B, ignore.strand = TRUE)))
  union_bp <- sum(width(GenomicRanges::union(A, B, ignore.strand = TRUE)))
  if (union_bp == 0L) return(0)
  inter_bp / union_bp
}

build_tables <- function() {
  metrics_df_list <- mapply(function(lbl, gr, anno_df, seq_df) {
    data.frame(
      label           = lbl,
      seqnames        = as.character(seqnames(gr)),
      start           = start(gr),
      end             = end(gr),
      width_bp        = width(gr),
      gc_content      = seq_df$gc,
      cpg_oe          = seq_df$cpgoe,
      annotation      = anno_df$annotation,
      dist_to_tss_bp  = anno_df$distanceToTSS,
      nearest_gene_id = anno_df$geneId,
      gene_biotype    = anno_df$geneBiotype %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }, labels, gr_list, anno_df_list, seq_metrics_list, SIMPLIFY = FALSE)
  metrics_all <- bind_rows(metrics_df_list)
  metrics_all$label <- factor(metrics_all$label, levels = labels)

  summary_list <- mapply(function(lbl, gr, anno_df, seq_df) {
    w  <- width(gr)
    af <- anno_fractions(anno_df$annotation)
    data.frame(
      label           = lbl,
      n_regions       = length(gr),
      total_bp        = sum(w),
      total_Mb        = round(sum(w) / 1e6, 3),
      median_width_bp = as.integer(median(w)),
      mean_width_bp   = round(mean(w), 1),
      sd_width_bp     = round(sd(w), 1),
      width_q25_bp    = as.integer(quantile(w, 0.25)),
      width_q75_bp    = as.integer(quantile(w, 0.75)),
      width_max_bp    = max(w),
      mean_gc         = round(mean(seq_df$gc,    na.rm = TRUE), 4),
      mean_cpg_oe     = round(mean(seq_df$cpgoe, na.rm = TRUE), 4),
      pct_promoter    = round(af$pct_promoter,   2),
      pct_exon        = round(af$pct_exon,       2),
      pct_intron      = round(af$pct_intron,     2),
      pct_downstream  = round(af$pct_downstream, 2),
      pct_intergenic  = round(af$pct_intergenic, 2),
      n_unique_chroms = length(unique(as.character(seqnames(gr)))),
      n_unique_genes  = length(unique(na.omit(anno_df$geneId))),
      median_tss_dist = as.integer(median(abs(anno_df$distanceToTSS), na.rm = TRUE)),
      mean_tss_dist   = round(mean(abs(anno_df$distanceToTSS), na.rm = TRUE), 1),
      stringsAsFactors = FALSE
    )
  }, labels, gr_list, anno_df_list, seq_metrics_list, SIMPLIFY = FALSE)
  summary_all <- bind_rows(summary_list)

  pair_grid <- expand.grid(set_A = labels, set_B = labels, stringsAsFactors = FALSE)
  overlap_df <- do.call(rbind, lapply(seq_len(nrow(pair_grid)), function(i) {
    A <- gr_list[[pair_grid$set_A[i]]]
    B <- gr_list[[pair_grid$set_B[i]]]
    data.frame(
      set_A               = pair_grid$set_A[i],
      set_B               = pair_grid$set_B[i],
      n_overlapping_pairs = length(findOverlaps(A, B, ignore.strand = TRUE)),
      jaccard_bp          = round(jaccard_bp(A, B), 5),
      stringsAsFactors    = FALSE
    )
  }))

  list(metrics = metrics_all, summary = summary_all, overlap = overlap_df)
}

tables <- tryCatch(build_tables(), error = function(e) {
  die(paste0("could not assemble the tables: ", conditionMessage(e)), 3L)
})
metrics_all <- tables$metrics

write_tsv <- function(df, name) {
  path <- paste0(out_prefix, ".", name)
  res <- tryCatch({
    write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
    TRUE
  }, error = function(e) {
    cat("ERROR: cannot write ", path, ": ", conditionMessage(e), "\n", sep = "", file = stderr())
    FALSE
  })
  if (res) cat(sprintf("  Written: %s (%d rows)\n", path, nrow(df)))
  res
}
tbl_ok <- c(
  write_tsv(tables$metrics,  "region_metrics.tsv"),
  write_tsv(tables$summary,  "summary_stats.tsv"),
  write_tsv(tables$overlap,  "overlap_matrix.tsv")
)
if (!all(tbl_ok)) die("could not write region_metrics.tsv, summary_stats.tsv and overlap_matrix.tsv", 3L)

# TSS windows (0-based BED) from getPromoters, same window as tssRegion
tryCatch({
  promoter_windows <- getPromoters(TxDb = txdb, upstream = tss_r, downstream = tss_r)
  tss_df <- as.data.frame(promoter_windows)[, c("seqnames", "start", "end")]
  tss_df$start <- tss_df$start - 1L
  write.table(tss_df, paste0(out_prefix, ".tss_windows.bed"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  cat(sprintf("  Written: %s.tss_windows.bed (%d regions)\n", out_prefix, nrow(tss_df)))
}, error = function(e) {
  cat("WARNING: tss_windows.bed not written: ", conditionMessage(e), "\n", sep = "", file = stderr())
})

# =============================================================================
# 7. Plots (every plot in svg(); a failing plot becomes a message page)
# =============================================================================
cat("\n=== Generating plots ===\n")

draw_message <- function(msg) {
  plot.new()
  text(0.5, 0.5, msg, cex = 0.9, col = "grey40")
}

save_plot <- function(name, width, height, draw) {
  path <- paste0(out_prefix, ".", name, ".svg")
  svg(path, width = width, height = height)
  err <- tryCatch({ draw(); NULL }, error = function(e) conditionMessage(e))
  if (!is.null(err)) {
    dev.off()
    cat("WARNING: plot failed [", name, "]: ", err, "\n", sep = "", file = stderr())
    svg(path, width = width, height = height)
    draw_message(paste0("Plot failed: ", name, "\n", err))
  }
  dev.off()
  cat(sprintf("  Saved: %s\n", basename(path)))
}

# --- annotation pie, one per label (ChIPseeker base graphics) ---
for (lbl in labels) {
  n_slices <- tryCatch(nrow(ChIPseeker:::getAnnoStat(anno_list[[lbl]])), error = function(e) 8L)
  n_slices <- max(n_slices, 2L)
  pie_cols <- rev(purple_ramp(n_slices))
  save_plot(paste0("annotation_pie_", slug(lbl)), 6, 6, function() {
    plotAnnoPie(anno_list[[lbl]], col = pie_cols)
    title(main = paste0(lbl, ": genomic annotation"), cex.main = 1.1)
  })
}

# --- annotation bar, all labels ---
save_plot("annotation_bar", 9, 5, function() {
  p <- plotAnnoBar(anno_list) +
    scale_fill_manual(values = COLS) +
    labs(title = "Genomic feature distribution by label",
         x = "Percentage (%)", fill = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "right",
          plot.title   = element_text(face = "bold"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  print(p)
})

# --- distance to TSS bars (ChIPseeker), all labels ---
save_plot("dist_to_tss", 9, 5, function() {
  p <- plotDistToTSS(anno_list, title = "Distance to nearest TSS") +
    scale_fill_manual(values = COLS) +
    labs(fill = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "right",
          plot.title   = element_text(face = "bold", hjust = 0.5),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  print(p)
})

# --- width distribution (violin + box, log scale) ---
save_plot("width_distribution", max(5, 2 * length(labels) + 2), 6, function() {
  df_w <- bind_rows(lapply(labels, function(lbl) {
    data.frame(label = lbl, width_bp = width(gr_list[[lbl]]))
  }))
  df_w$label <- factor(df_w$label, levels = labels)
  p <- ggplot(df_w, aes(x = label, y = width_bp, fill = label)) +
    geom_violin(trim = TRUE, alpha = 0.75, color = NA) +
    geom_boxplot(width = 0.10, outlier.shape = NA,
                 fill = "white", color = "grey30", linewidth = 0.5) +
    scale_y_log10(labels = label_comma(), name = "Peak width (bp, log10 scale)") +
    scale_fill_manual(values = COLS) +
    labs(title = "Peak width distribution", x = "Label", fill = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
  print(p)
})

# --- peak width ECDF with density panel ---
save_plot("width_ecdf", 10, 6, function() {
  df_w <- bind_rows(lapply(labels, function(lbl) {
    data.frame(label = lbl, width_kb = width(gr_list[[lbl]]) / 1e3)
  }))
  df_w$label <- factor(df_w$label, levels = labels)

  p_ecdf <- ggplot(df_w, aes(x = width_kb, colour = label)) +
    stat_ecdf(geom = "step", linewidth = 1.1, pad = FALSE) +
    scale_x_log10(name = "Peak width (kb, log10)", labels = label_comma(suffix = " kb")) +
    scale_colour_manual(values = COLS) +
    scale_y_continuous(name = "Cumulative fraction", labels = label_percent()) +
    labs(title = "Peak width ECDF",
         subtitle = "Cumulative distribution of peak widths",
         colour = NULL) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))

  p_dens <- ggplot(df_w, aes(x = width_kb, fill = label)) +
    geom_density(alpha = 0.55, linewidth = 0.4) +
    scale_x_log10(name = "Peak width (kb)", labels = label_comma(suffix = " kb")) +
    scale_fill_manual(values = COLS) +
    labs(y = "Density", fill = NULL) +
    theme_bw(base_size = 11) +
    theme(legend.position = "none")

  grid.arrange(p_ecdf, p_dens, ncol = 2, widths = c(1.6, 1))
})

# --- GC content and CpG o/e ---
save_plot("gc_cpgoe", 11, 5, function() {
  df_gc <- metrics_all %>% filter(!is.na(gc_content))
  if (nrow(df_gc) == 0L) {
    draw_message(paste0("GC and CpG o/e not computed.\n",
                        "BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0 unavailable or sequences not found."))
  } else {
    p_gc <- ggplot(df_gc, aes(x = gc_content, fill = label)) +
      geom_density(alpha = 0.5) +
      scale_fill_manual(values = COLS) +
      labs(title = "GC content per peak", x = "GC fraction", y = "Density", fill = NULL) +
      theme_classic(base_size = 14) +
      theme(legend.position = "top",
            plot.title = element_text(face = "bold", size = 13),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

    p_cpg <- ggplot(df_gc %>% filter(!is.na(cpg_oe)), aes(x = cpg_oe, fill = label)) +
      geom_density(alpha = 0.5) +
      coord_cartesian(xlim = c(0, 1.5)) +
      scale_fill_manual(values = COLS) +
      labs(title = "CpG O/E ratio per peak",
           subtitle = "(nCpG x length) / (nC x nG)",
           x = "CpG O/E", y = "Density", fill = NULL) +
      theme_classic(base_size = 14) +
      theme(legend.position = "top",
            plot.title = element_text(face = "bold", size = 13),
            panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))

    grid.arrange(p_gc, p_cpg, ncol = 2)
  }
})

# --- peak count per chromosome ---
save_plot("chromosome_distribution", 12, 6, function() {
  df_chr <- bind_rows(lapply(labels, function(lbl) {
    data.frame(label = lbl, chr = as.character(seqnames(gr_list[[lbl]])))
  })) %>%
    filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%
    mutate(chr   = factor(chr, levels = paste0("chr", c(as.character(1:22), "X", "Y"))),
           label = factor(label, levels = labels)) %>%
    group_by(label, chr) %>%
    summarise(n_peaks = n(), .groups = "drop")

  p <- ggplot(df_chr, aes(x = chr, y = n_peaks, fill = label)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
    scale_fill_manual(values = COLS) +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Peak count by chromosome",
         x = "Chromosome", y = "Number of peaks", fill = NULL) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9),
          legend.position = "top",
          plot.title = element_text(face = "bold"))
  print(p)
})

# --- distance to nearest TSS, cumulative distribution ---
save_plot("tss_dist_ecdf", 5.25, 5, function() {
  df_tss <- metrics_all %>%
    mutate(abs_dist_kb = abs(dist_to_tss_bp) / 1e3) %>%
    filter(!is.na(abs_dist_kb))
  df_tss$label <- factor(df_tss$label, levels = labels)

  p <- ggplot(df_tss, aes(x = abs_dist_kb, colour = label)) +
    stat_ecdf(geom = "step", linewidth = 1.1, pad = FALSE) +
    scale_x_log10(name = "Distance to nearest TSS (kb, log10)",
                  labels = label_comma(suffix = " kb"),
                  limits = c(0.01, NA)) +
    scale_y_continuous(name = "Cumulative fraction", labels = label_percent()) +
    scale_colour_manual(values = COLS) +
    labs(title = "Distance to nearest TSS: cumulative distribution", colour = NULL) +
    theme_classic(base_size = 14) +
    theme(legend.position = "top",
          plot.title = element_text(face = "bold", hjust = 0.5),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1))
  print(p)
})

# --- total coverage per chromosome, all labels ---
save_plot("chromosome_coverage", 12, 6, function() {
  df_cov <- bind_rows(lapply(labels, function(lbl) {
    gr <- gr_list[[lbl]]
    data.frame(label = lbl, chr = as.character(seqnames(gr)), width_bp = width(gr))
  })) %>%
    filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%
    mutate(chr   = factor(chr, levels = paste0("chr", c(as.character(1:22), "X", "Y"))),
           label = factor(label, levels = labels)) %>%
    group_by(label, chr) %>%
    summarise(total_Mb = sum(width_bp) / 1e6, .groups = "drop")

  p <- ggplot(df_cov, aes(x = chr, y = total_Mb, fill = label)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.85) +
    scale_fill_manual(values = COLS) +
    labs(title = "Genomic coverage by chromosome",
         x = "Chromosome", y = "Total coverage (Mb)", fill = NULL) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 9),
          legend.position = "top",
          plot.title = element_text(face = "bold"))
  print(p)
})

# --- coverage per chromosome, one per label (highest to lowest) ---
for (lbl in labels) {
  save_plot(paste0("chromosome_coverage_", slug(lbl)), 10, 5, function() {
    df_lbl <- data.frame(chr = as.character(seqnames(gr_list[[lbl]])),
                         width_bp = width(gr_list[[lbl]])) %>%
      filter(grepl("^chr([0-9]+|X|Y)$", chr)) %>%
      group_by(chr) %>%
      summarise(total_Mb = sum(width_bp) / 1e6, .groups = "drop") %>%
      arrange(desc(total_Mb)) %>%
      mutate(chr = factor(chr, levels = chr))

    p <- ggplot(df_lbl, aes(x = chr, y = total_Mb)) +
      geom_bar(stat = "identity", fill = COLS[[lbl]], alpha = 0.87, width = 0.7) +
      scale_y_continuous(labels = label_comma(suffix = " Mb"), name = "Total coverage (Mb)") +
      labs(title = paste0(lbl, ": coverage by chromosome"),
           subtitle = "Chromosomes ordered by total coverage (high to low)",
           x = NULL) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 9),
            plot.title = element_text(face = "bold"),
            panel.grid.major.x = element_blank())
    print(p)
  })
}

# --- comparative KEGG enrichment (compareCluster + dotplot) ---
# GTF gene_id is a gene symbol, so symbols are converted to ENTREZ IDs (bitr, SYMBOL -> ENTREZID).
run_kegg <- function() {
  if (!kegg_on) return(list(obj = NULL, reason = "KEGG enrichment disabled by --kegg=false."))
  if (!CLUSTERP_AVAILABLE) {
    return(list(obj = NULL, reason = "clusterProfiler and/or org.Hs.eg.db not installed."))
  }
  genes_by_label <- lapply(anno_list, function(x) unique(na.omit(as.character(as.data.frame(x)$geneId))))
  genes_entrez <- lapply(genes_by_label, function(sym) {
    if (length(sym) == 0L) return(character(0))
    tryCatch(
      bitr(sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID,
      error = function(e) {
        cat("WARNING: SYMBOL to ENTREZID conversion failed: ", conditionMessage(e), "\n", sep = "", file = stderr())
        character(0)
      }
    )
  })
  if (all(lengths(genes_entrez) == 0L)) {
    return(list(obj = NULL, reason = "no gene symbols could be converted to ENTREZ IDs."))
  }
  ck <- compareCluster(geneClusters  = genes_entrez,
                       fun           = "enrichKEGG",
                       organism      = "hsa",
                       pvalueCutoff  = 0.05,
                       pAdjustMethod = "BH")
  list(obj = ck, reason = NULL)
}

kegg <- tryCatch(run_kegg(), error = function(e) {
  cat("WARNING: KEGG enrichment failed: ", conditionMessage(e), "\n", sep = "", file = stderr())
  list(obj = NULL, reason = paste0("KEGG enrichment failed (needs internet access for KEGG): ",
                                   conditionMessage(e)))
})

save_plot("kegg", 10, 8, function() {
  if (is.null(kegg$obj)) {
    draw_message(paste0("KEGG pathway enrichment not drawn.\n", kegg$reason))
  } else if (nrow(as.data.frame(kegg$obj)) == 0L) {
    draw_message("KEGG pathway enrichment: no terms enriched at p < 0.05.")
  } else {
    p <- dotplot(kegg$obj, showCategory = 15, title = "Comparative KEGG pathway enrichment") +
      theme_classic(base_size = 13) +
      theme(plot.title  = element_text(face = "bold", hjust = 0.5),
            axis.text.x = element_text(angle = 30, hjust = 1))
    print(p)
  }
})

# =============================================================================
# 8. Done
# =============================================================================
cat(sprintf("\n=== peak_annotation.R complete (prefix: %s) ===\n", out_prefix))
