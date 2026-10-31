#!/usr/bin/env Rscript
# deg_deseq2.R -- build the DESeqDataSet of one model, pre-filter low-count genes and fit DESeq2 (Wald test).
#
# Legacy mapping (step_6_deg_analysis.Rmd): chunks "Create DESeqDataSet and DESeq object" and "export deseq objects".
#
# CLI (all --key=value):
#   --common=deg_common.R
#   --counts=counts_edited.rds        data.frame genes x samples
#   --metadata=sample_metadata_edited.rds
#   --model=<name>                    output sub-folder
#   --design="~treatment"             DESeq2 design formula (e.g. "~ target * drug")
#   --min-count=10 --min-samples=2    keep genes with at least min-count reads in at least min-samples samples
#                                     (legacy: rowSums(counts >= 10) >= 2)
#   --fit-type=parametric             DESeq2 dispersion fit type (parametric|local|mean|glmGamPoi)
#   --svg-device=svglite              svglite (real <text>, default) or cairo (legacy, glyph outlines)
#   --font=Helvetica                  font family of dispersion_plots.svg (see deg_common.R set_font)
#   --out-dir=.
# Outputs (in <out-dir>/<model>/):
#   DESeqDataSet.rds   filtered DESeqDataSet BEFORE DESeq() (as the Rmd saved it)
#   DESeq.rds          fitted DESeqDataSet
#   design.txt         paste(design(dds), collapse = " ") as the Rmd
#   dispersion_plots.svg   plotDispEsts, svg 4.5 x 4.5 in
#   prefilter.tsv      model, genes_before, genes_after, min_count, min_samples
# Exit codes: 0 ok, 3 data problem (e.g. design not full rank, no gene passes the filter), 1 unexpected error.

.argv <- commandArgs(trailingOnly = FALSE)
.cm <- sub("^--common=", "", grep("^--common=", .argv, value = TRUE)[1])
if (is.na(.cm) || !file.exists(.cm)) .cm <- file.path(dirname(sub("^--file=", "", grep("^--file=", .argv, value = TRUE)[1])), "deg_common.R")
if (!file.exists(.cm)) {
  cat("ERROR: --common=<deg_common.R> not found\n", file = stderr())
  quit(save = "no", status = 3)
}
source(.cm)
suppressPackageStartupMessages(library(DESeq2))

opt <- parse_args()
counts_path <- opt_need(opt, "counts")
meta_path <- opt_need(opt, "metadata")
model <- opt_need(opt, "model")
design_str <- opt_need(opt, "design")
min_count <- as_num(opt_get(opt, "min-count", "10"), "min-count")
min_samples <- as_num(opt_get(opt, "min-samples", "2"), "min-samples")
fit_type <- opt_get(opt, "fit-type", "parametric")
out_dir <- opt_get(opt, "out-dir", ".")
set_font(opt)
if (!fit_type %in% c("parametric", "local", "mean", "glmGamPoi")) die("--fit-type must be parametric, local, mean or glmGamPoi")

for (p in c(counts_path, meta_path)) if (!file.exists(p)) die(paste0("file not found: ", p))
counts <- readRDS(counts_path)
sample_metadata <- readRDS(meta_path)
if (!identical(colnames(counts), rownames(sample_metadata))) die("counts columns and metadata rows differ; run deg_prepare.R first")

# the formula is created in the global environment so that saved objects do not drag a function frame along
design <- tryCatch(as.formula(design_str, env = globalenv()), error = function(e) die(paste0("bad design '", design_str, "'")))
miss <- setdiff(all.vars(design), names(sample_metadata))
if (length(miss)) die(paste0("design variable(s) not in the metadata: ", paste(miss, collapse = ", ")))

dataSet <- tryCatch(
  DESeqDataSetFromMatrix(countData = counts, colData = sample_metadata, design = design),
  error = function(e) die(paste0("DESeqDataSetFromMatrix failed: ", conditionMessage(e)))
)

# Pre-filtering genes with low counts
genes_before <- nrow(dataSet)
keep <- rowSums(counts(dataSet) >= min_count) >= min_samples
dataSet <- dataSet[keep, ]
genes_after <- nrow(dataSet)
cat("model ", model, ": ", genes_before, " genes, ", genes_after, " kept (>= ", min_count, " reads in >= ", min_samples, " samples)\n", sep = "")
if (genes_after < 10) die(paste0("only ", genes_after, " genes pass the pre-filter; check --min-count / --min-samples"))

DESeq <- tryCatch(
  DESeq(dataSet, fitType = fit_type),
  error = function(e) {
    msg <- conditionMessage(e)
    if (grepl("full rank|rank|levels", msg)) die(paste0("DESeq2 cannot fit the design: ", msg))
    stop(e)
  }
)

mdir <- file.path(out_dir, model)
dir.create(mdir, showWarnings = FALSE, recursive = TRUE)

# Estimate dispersion plots
open_svg(file.path(mdir, "dispersion_plots.svg"), 4.5, 4.5)
DESeq2::plotDispEsts(DESeq)
close_svg()

saveRDS(DESeq, file = file.path(mdir, "DESeq.rds"))
saveRDS(dataSet, file = file.path(mdir, "DESeqDataSet.rds"))
writeLines(paste(dataSet@design, collapse = " "), con = file.path(mdir, "design.txt"))
write.table(
  data.frame(model = model, genes_before = genes_before, genes_after = genes_after, min_count = min_count, min_samples = min_samples),
  file.path(mdir, "prefilter.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
cat("coefficients: ", paste(resultsNames(DESeq), collapse = ", "), "\n", sep = "")
