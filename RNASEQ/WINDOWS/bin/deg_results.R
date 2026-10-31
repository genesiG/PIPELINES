#!/usr/bin/env Rscript
# deg_results.R -- Wald-test results, LFC shrinkage, up/down classification and summary for all comparisons of a model.
#
# Legacy mapping (step_6_deg_analysis.Rmd): chunk "Statistical test (Walds Test)" (results table, shrunk LFC, effect
#   classification, analysis_summary.csv / .txt).
#
# CLI (all --key=value; booleans TRUE/FALSE):
#   --common=deg_common.R
#   --deseq=DESeq.rds                 fitted DESeqDataSet of the model (deg_deseq2.R)
#   --comparisons=comparisons.tsv     normalised table (deg_prepare.R); rows of --model are processed in table order
#   --model=<name>
#   --alpha=0.05                      FDR cut-off (results(alpha =)) and padj threshold of the classification
#   --lfc=1                           |log2FC| threshold of the classification
#   --lfc-threshold-in-results=FALSE  TRUE passes lfcThreshold = lfc to results() AND to lfcShrink() (FSOS s-value, T = lfc;
#                                     the Rmd never ran this branch and did not pass it to lfcShrink) and classifies by sign only
#   --shrink-type=ashr|apeglm|none    lfcShrink(type =, svalue = TRUE); apeglm only works with `coef` rows
#                                     (alpha is never passed to lfcShrink(), so ashr always uses its default model)
#   --sval-threshold=0.005            s-value threshold of shrunk_effect (0.005 is equivalent to alpha 0.05)
#   --legacy-summary-counts=FALSE     TRUE reproduces the legacy summary bug (see below)
#   --out-dir=.
# Per comparison row the DESeq2 call is
#   factor/numerator/denominator -> results(dds, contrast = c(factor, numerator, denominator), alpha =)
#   coef                         -> results(dds, name = coef, alpha =)               (lfcShrink: coef =)
#   add[;..] / subtract[;..]     -> results(dds, contrast = list(add, subtract), alpha =)  (subtract empty: list(add))
# and lfcShrink(dds, contrast = / coef = (same form), res = res, type =, svalue = TRUE) unless --shrink-type=none.
# Outputs (in <out-dir>/<model>/):
#   <comparison>/results_<comparison>.csv   write.csv(quote = FALSE, row.names = FALSE); columns exactly
#       Gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,shrunkLFC,svalue,effect,shrunk_effect
#       (sorted by padj, NA last; shrunkLFC,svalue,shrunk_effect omitted with --shrink-type=none).
#       effect / shrunk_effect = Upregulated | Downregulated | NS; NA padj / svalue = NS.
#   analysis_summary.csv, analysis_summary.txt   as the Rmd (comparison = label column, or comparison with "_" -> " ")
#   normalized_counts.csv   Gene + counts(dds, normalized = TRUE) (new, not in the Rmd)
# LEGACY BUG: the Rmd counted up / down genes as effect == "Enriched" / "Depleted" (labels that never occur), so its
#   analysis_summary showed 0 up / 0 down and non_sig = the "NS" count only. Default: counts "Upregulated" /
#   "Downregulated" (correct). --legacy-summary-counts=TRUE reproduces the zeros. The "downregualted" typo of the
#   legacy txt report is corrected in both modes.
# Exit codes: 0 ok, 3 data problem (unknown coefficient / level, apeglm with a non-coef row), 1 unexpected error.

.argv <- commandArgs(trailingOnly = FALSE)
.cm <- sub("^--common=", "", grep("^--common=", .argv, value = TRUE)[1])
if (is.na(.cm) || !file.exists(.cm)) .cm <- file.path(dirname(sub("^--file=", "", grep("^--file=", .argv, value = TRUE)[1])), "deg_common.R")
if (!file.exists(.cm)) {
  cat("ERROR: --common=<deg_common.R> not found\n", file = stderr())
  quit(save = "no", status = 3)
}
source(.cm)
suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tibble)
})

opt <- parse_args()
deseq_path <- opt_need(opt, "deseq")
comp_path <- opt_need(opt, "comparisons")
model <- opt_need(opt, "model")
alpha <- as_num(opt_get(opt, "alpha", "0.05"), "alpha")
lfc <- as_num(opt_get(opt, "lfc", "1"), "lfc")
use_lfc_threshold_in_results <- as_bool(opt_get(opt, "lfc-threshold-in-results", "FALSE"), "lfc-threshold-in-results")
shrink_type <- opt_get(opt, "shrink-type", "ashr")
sval_threshold <- as_num(opt_get(opt, "sval-threshold", "0.005"), "sval-threshold")
legacy_summary <- as_bool(opt_get(opt, "legacy-summary-counts", "FALSE"), "legacy-summary-counts")
out_dir <- opt_get(opt, "out-dir", ".")
if (!shrink_type %in% c("ashr", "apeglm", "none")) die("--shrink-type must be ashr, apeglm or none")
do_shrink <- shrink_type != "none"

if (!file.exists(deseq_path)) die(paste0("file not found: ", deseq_path))
DESeq <- readRDS(deseq_path)
comp <- read_comparisons_tsv(comp_path)
comp <- comp[comp$model == model, , drop = FALSE]
if (nrow(comp) == 0) die(paste0("no comparison rows for model '", model, "' in ", comp_path))
mdir <- file.path(out_dir, model)
dir.create(mdir, showWarnings = FALSE, recursive = TRUE)

split_names <- function(x) {
  v <- trimws(strsplit(x, ";", fixed = TRUE)[[1]])
  v[nzchar(v)]
}
coef_names <- resultsNames(DESeq)
check_coefs <- function(v, comparison) {
  bad <- setdiff(v, coef_names)
  if (length(bad)) {
    die(paste0(
      "comparison '", comparison, "': coefficient(s) not in resultsNames(dds): ", paste(bad, collapse = ", "),
      "\n  available: ", paste(coef_names, collapse = ", ")
    ))
  }
}

de_summary <- list()
for (i in seq_len(nrow(comp))) {
  comparisons <- comp$comparison[i]
  comp_label <- comparison_label(comp[i, , drop = FALSE])
  cat(paste0("\nPerforming DEG analysis for comparison: ", comp_label, "\n"))

  # --- the contrast / coefficient of this row (one kind per row) ---
  res_args <- list(object = DESeq, alpha = alpha)
  shrink_args <- list(dds = DESeq, type = shrink_type, svalue = TRUE)
  if (!is_blank(comp$coef[i])) {
    check_coefs(comp$coef[i], comparisons)
    res_args$name <- comp$coef[i]
    shrink_args$coef <- comp$coef[i]
  } else if (!is_blank(comp$add[i])) {
    add <- split_names(comp$add[i])
    subt <- split_names(comp$subtract[i])
    check_coefs(c(add, subt), comparisons)
    contrast <- if (length(subt)) list(add, subt) else list(add)
    res_args$contrast <- contrast
    shrink_args$contrast <- contrast
  } else {
    f <- comp$factor[i]
    if (!f %in% colnames(colData(DESeq))) die(paste0("comparison '", comparisons, "': factor '", f, "' not in the DESeq2 colData"))
    lv <- levels(colData(DESeq)[[f]])
    bad <- setdiff(c(comp$numerator[i], comp$denominator[i]), lv)
    if (length(bad)) die(paste0("comparison '", comparisons, "': level(s) ", paste(bad, collapse = ", "), " not in factor '", f, "'"))
    contrast <- c(f, comp$numerator[i], comp$denominator[i])
    res_args$contrast <- contrast
    shrink_args$contrast <- contrast
  }
  if (shrink_type == "apeglm" && is.null(shrink_args$coef)) {
    die(paste0("comparison '", comparisons, "': --shrink-type=apeglm needs a `coef` row (character / add-subtract contrasts are not supported); use ashr"))
  }

  # --- Building the results table ---
  if (use_lfc_threshold_in_results) res_args$lfcThreshold <- lfc
  results <- do.call(DESeq2::results, res_args)

  # --- Calculate shrunk LFCs ---
  if (do_shrink) {
    shrink_args$res <- results
    # With results(lfcThreshold = lfc) the padj tests |LFC| > lfc, so the s-value must be the matching "false sign or
    # small" s-value (T = lfc); without it the shrunk side would only test the sign (3952 "Upregulated" genes with
    # mostly |shrunk LFC| < 1 against 29 in padj in the GG CBX2_1uM_vs_LUC_DMSO test). Never reached with the default
    # (lfc-threshold-in-results=FALSE, as the Rmd: lfcThreshold = 0 -> classic s-value).
    if (use_lfc_threshold_in_results) shrink_args$lfcThreshold <- lfc
    shrunk_results <- do.call(DESeq2::lfcShrink, shrink_args)
  }

  comp_dir <- file.path(mdir, comparisons)
  dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

  # --- Process results table ---
  results_df <- as.data.frame(results) %>% arrange(padj)
  results_df <- results_df %>% rownames_to_column(var = "Gene")

  # Add shrunk LFC information
  if (do_shrink) {
    shrunk_df <- as.data.frame(shrunk_results) %>% rownames_to_column(var = "Gene")
    results_df$shrunkLFC <- shrunk_df$log2FoldChange[match(results_df$Gene, shrunk_df$Gene)]
    results_df$svalue <- shrunk_df$svalue[match(results_df$Gene, shrunk_df$Gene)]
  }

  # NA padj / svalue never satisfy the conditions, so they stay "NS"
  cat("\nDetermining DEG effects\n")
  results_df$effect <- "NS"
  if (use_lfc_threshold_in_results) {
    results_df$effect[which(results_df$padj < alpha & results_df$log2FoldChange <= 0)] <- "Downregulated"
    results_df$effect[which(results_df$padj < alpha & results_df$log2FoldChange > 0)] <- "Upregulated"
  } else {
    results_df$effect[which(results_df$padj < alpha & results_df$log2FoldChange <= -lfc)] <- "Downregulated"
    results_df$effect[which(results_df$padj < alpha & results_df$log2FoldChange >= lfc)] <- "Upregulated"
  }
  if (do_shrink) {
    results_df$shrunk_effect <- "NS"
    if (use_lfc_threshold_in_results) {
      results_df$shrunk_effect[which(results_df$svalue < sval_threshold & results_df$shrunkLFC <= 0)] <- "Downregulated"
      results_df$shrunk_effect[which(results_df$svalue < sval_threshold & results_df$shrunkLFC > 0)] <- "Upregulated"
    } else {
      results_df$shrunk_effect[which(results_df$svalue < sval_threshold & results_df$shrunkLFC <= -lfc)] <- "Downregulated"
      results_df$shrunk_effect[which(results_df$svalue < sval_threshold & results_df$shrunkLFC >= lfc)] <- "Upregulated"
    }
  }

  # Export data frame with results
  write.csv(results_df, file = file.path(comp_dir, paste0("results_", comparisons, ".csv")), quote = FALSE, row.names = FALSE)

  # --- Calculate summary statistics ---
  if (legacy_summary) {
    up_regulated <- sum(results_df$effect == "Enriched") # legacy: label never occurs -> 0
    down_regulated <- sum(results_df$effect == "Depleted") # legacy: label never occurs -> 0
  } else {
    up_regulated <- sum(results_df$effect == "Upregulated")
    down_regulated <- sum(results_df$effect == "Downregulated")
  }
  non_sig <- sum(results_df$effect == "NS")
  total_genes <- nrow(results_df)

  de_summary[[comparisons]] <- list(
    comparison = comp_label,
    up_regulated = up_regulated,
    down_regulated = down_regulated,
    non_sig = non_sig,
    total_genes = total_genes,
    up_percent = round((up_regulated / total_genes) * 100, 2),
    down_percent = round((down_regulated / total_genes) * 100, 2)
  )
}

# Create comprehensive differential expression summary report
if (length(de_summary) > 0) {
  summary_df <- do.call(rbind, lapply(de_summary, as.data.frame))

  # Write CSV report
  write.csv(summary_df, file = file.path(mdir, "analysis_summary.csv"), row.names = FALSE)

  # Formatted text report (same layout as the Rmd, which used sink())
  con <- file(file.path(mdir, "analysis_summary.txt"), open = "wt")
  w <- function(x) cat(x, file = con)
  w("RNA-seq ANALYSIS SUMMARY REPORT\n")
  w("===============================================\n\n")
  w(paste("Generated:", Sys.Date(), "\n\n"))
  w(paste("Total comparisons analyzed:", nrow(summary_df), "\n\n"))

  # Overall statistics
  w("OVERALL STATISTICS:\n")
  w("===================\n")
  w(paste("Total genes analyzed:", sum(summary_df$total_genes), "\n"))
  w(paste("Average genes per comparison:", round(mean(summary_df$total_genes), 1), "\n"))
  w(paste("Total upregulated genes:", sum(summary_df$up_regulated), "\n"))
  w(paste("Total downregulated genes:", sum(summary_df$down_regulated), "\n"))
  w(paste("Average % up:", round(mean(summary_df$up_percent), 1), "%\n"))
  w(paste("Average % down:", round(mean(summary_df$down_percent), 1), "%\n\n"))

  # Individual comparison details
  w("DETAILED COMPARISON BREAKDOWN:\n")
  w("==============================\n\n")
  for (i in 1:nrow(summary_df)) {
    comp_info <- summary_df[i, ]
    w(paste("Comparison:", comp_info$comparison, "\n"))
    w(paste("Total genes:", comp_info$total_genes, "\n"))
    w(paste("Upregulated genes:", comp_info$up_regulated, "(", comp_info$up_percent, "%)\n"))
    w(paste("Downregulated genes:", comp_info$down_regulated, "(", comp_info$down_percent, "%)\n"))
    w(paste(
      "Non-significant genes:", comp_info$non_sig,
      "(", round(100 - comp_info$up_percent - comp_info$down_percent, 1), "%)\n"
    ))
    w("----------------------------------------\n\n")
  }

  # Identify comparisons with the most differential expression
  most_de <- summary_df[which.max(summary_df$up_regulated + summary_df$down_regulated), ]
  w("COMPARISON WITH THE MOST DIFFERENTIAL EXPRESSION:\n")
  w("=================================================\n")
  w(paste("Comparison:", most_de$comparison, "\n"))
  w(paste("Total differentially expressed genes:", most_de$up_regulated + most_de$down_regulated, "\n"))
  w(paste("Percentage differentially expressed genes:", round(most_de$up_percent + most_de$down_percent, 1), "%\n\n"))
  close(con)
}

# normalised counts (new)
nc <- counts(DESeq, normalized = TRUE)
write.csv(data.frame(Gene = rownames(nc), nc, check.names = FALSE),
  file = file.path(mdir, "normalized_counts.csv"), quote = FALSE, row.names = FALSE
)
cat("\ndeg_results: ", nrow(comp), " comparison(s) of model ", model, " written to ", mdir, "\n", sep = "")
