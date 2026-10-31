#!/usr/bin/env Rscript
# deg_transform.R -- variance-stabilising transformation (rlog / vst) of the filtered DESeqDataSet of a model.
#
# Legacy mapping (step_6_deg_analysis.Rmd): first lines of the chunk "Gene Expression Heatmap": rlog(dataSet, blind = FALSE)
#   on the UNFITTED filtered DESeqDataSet, rlogData.rds, and rlogCounts.rds (assay as data.frame with a Gene column).
#
# CLI (all --key=value):
#   --common=deg_common.R
#   --dataset=DESeqDataSet.rds        filtered, unfitted DESeqDataSet (deg_deseq2.R)
#   --model=<name>
#   --method=rlog|vst                 vst falls back to varianceStabilizingTransformation() when < 1000 genes
#   --blind=FALSE                     blind to the design (legacy: FALSE)
#   --fit-type=parametric             dispersion trend used by the transformation (parametric|local|mean). The dataset is
#                                     UNFITTED, so the trend is re-estimated here (identical to the one of DESeq() when the
#                                     same fit type is used: max |rlog(unfitted) - rlog(fitted)| = 0 on the GG data). The
#                                     module forwards deg_fit_type, the same value DEG_DESEQ2 uses for DESeq().
#   --out-dir=.
# Outputs (in <out-dir>/<model>/): <method>Data.rds (DESeqTransform), <method>Counts.rds (data.frame of assay(),
#   Gene column first via rownames_to_column), i.e. rlogData.rds + rlogCounts.rds for the default.
# Exit codes: 0 ok, 3 data problem, 1 unexpected error.

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
dataset_path <- opt_need(opt, "dataset")
model <- opt_need(opt, "model")
method <- opt_get(opt, "method", "rlog")
blind <- as_bool(opt_get(opt, "blind", "FALSE"), "blind")
fit_type <- opt_get(opt, "fit-type", "parametric")
out_dir <- opt_get(opt, "out-dir", ".")
if (!method %in% c("rlog", "vst")) die("--method must be rlog or vst")
if (!fit_type %in% c("parametric", "local", "mean")) die("--fit-type must be parametric, local or mean")
if (!file.exists(dataset_path)) die(paste0("file not found: ", dataset_path))

dataSet <- readRDS(dataset_path)
mdir <- file.path(out_dir, model)
dir.create(mdir, showWarnings = FALSE, recursive = TRUE)

if (method == "rlog") {
  trData <- rlog(dataSet, blind = blind, fitType = fit_type)
} else if (nrow(dataSet) >= 1000) {
  trData <- vst(dataSet, blind = blind, fitType = fit_type)
} else {
  trData <- varianceStabilizingTransformation(dataSet, blind = blind, fitType = fit_type)
}
saveRDS(trData, file = file.path(mdir, paste0(method, "Data.rds")))

trCounts <- assay(trData) %>% as.data.frame()
trCounts <- trCounts %>% rownames_to_column(var = "Gene")
saveRDS(trCounts, file = file.path(mdir, paste0(method, "Counts.rds")))
cat("deg_transform: ", method, " (blind = ", blind, ") of ", nrow(dataSet), " genes written to ", mdir, "\n", sep = "")
