#!/usr/bin/env Rscript
# deg_prepare.R -- load, check and align the counts table and the sample metadata; normalise the comparisons table.
#
# Legacy mapping (step_6_deg_analysis.Rmd): chunks "Load Gene Expression Data", "Load Metadata", "Prepare Metadata",
#   "Prepare counts matrix for data analysis" (suffix stripped from count columns AND metadata sample ids, metadata
#   order defines the column order, factors + relevel, counts_edited.rds / sample_metadata_edited.rds/.csv).
#
# CLI (all --key=value):
#   --common=deg_common.R        shared helpers (required)
#   --counts=<file>              tab-separated counts, first column = Gene, header = sample ids
#   --metadata=<csv|tsv>         one row per sample, a column with the sample id
#   --comparisons=<csv>          comparisons table (see the pipeline docs, columns model,design,comparison,...)
#   --sample-col=sample_id       metadata column holding the sample ids
#   --strip-suffix=<regex>       removed (gsub) from count column names and metadata sample ids; empty = none
#   --reference-levels="treatment=LUC_DMSO,drug=DMSO"   factor=reference level list (relevel)
#   --drop-excluded=TRUE         drop metadata rows whose `exclude` column is true/1/yes (column optional)
#   --out-dir=.
# Behaviour: duplicated genes, metadata samples absent from the counts, non-integer / negative counts, design
#   variables that are missing / NA / empty / single-level, reference levels that do not exist, and invalid
#   comparison rows are data errors (exit 3). Count columns without metadata are dropped (message).
#   Design variables = all.vars() of every design formula in the comparisons table; they become factors. Factors
#   without an entry in --reference-levels keep the R default (alphabetical) first level; this is printed.
# Outputs (in --out-dir): counts_edited.rds (data.frame, genes x samples, metadata order), counts_edited.tsv,
#   sample_metadata_edited.rds (data.frame, rownames = sample ids), sample_metadata_edited.csv (write.csv,
#   row.names = TRUE, quote = FALSE as the Rmd), comparisons.tsv (normalised, tab-separated, all columns),
#   prepare_summary.txt.
# Exit codes: 0 ok, 3 data problem, 1 unexpected error.

.argv <- commandArgs(trailingOnly = FALSE)
.cm <- sub("^--common=", "", grep("^--common=", .argv, value = TRUE)[1])
if (is.na(.cm) || !file.exists(.cm)) .cm <- file.path(dirname(sub("^--file=", "", grep("^--file=", .argv, value = TRUE)[1])), "deg_common.R")
if (!file.exists(.cm)) {
  cat("ERROR: --common=<deg_common.R> not found\n", file = stderr())
  quit(save = "no", status = 3)
}
source(.cm)

opt <- parse_args()
counts_path <- opt_need(opt, "counts")
meta_path <- opt_need(opt, "metadata")
comp_path <- opt_need(opt, "comparisons")
sample_col <- opt_get(opt, "sample-col", "sample_id")
strip_re <- opt_get(opt, "strip-suffix", "")
ref_arg <- opt_get(opt, "reference-levels", "")
drop_excluded <- as_bool(opt_get(opt, "drop-excluded", "TRUE"), "drop-excluded")
out_dir <- opt_get(opt, "out-dir", ".")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

strip_ids <- function(x) if (nzchar(strip_re)) gsub(strip_re, "", x) else x

# ---------------------------------------------------------------- comparisons table ----
comp <- read_comparisons_tsv(comp_path)
if (nrow(comp) == 0) die("the comparisons table has no rows")
comp$design <- trimws(gsub("\\s+", " ", comp$design))
for (cn in c("model", "design", "comparison")) {
  if (any(is_blank(comp[[cn]]))) die(paste0("comparisons table: empty '", cn, "' in row(s) ", paste(which(is_blank(comp[[cn]])), collapse = ",")))
}
safe_re <- "^[A-Za-z0-9._+-]+$"
if (any(!grepl(safe_re, comp$model))) die("comparisons table: model names may only contain A-Za-z0-9._+-")
if (any(!grepl(safe_re, comp$comparison))) die("comparisons table: comparison names may only contain A-Za-z0-9._+-")
dup <- duplicated(paste(comp$model, comp$comparison))
if (any(dup)) die(paste0("comparisons table: duplicated comparison within a model: ", paste(unique(paste(comp$model, comp$comparison)[dup]), collapse = "; ")))

# one design per model, parseable one-sided formula
design_of <- tapply(comp$design, comp$model, function(d) unique(d))
for (m in names(design_of)) {
  if (length(design_of[[m]]) != 1) die(paste0("comparisons table: model '", m, "' has more than one design: ", paste(design_of[[m]], collapse = " | ")))
}
design_vars <- character(0)
for (m in names(design_of)) {
  f <- tryCatch(as.formula(design_of[[m]]), error = function(e) NULL)
  if (is.null(f) || length(f) != 2) die(paste0("model '", m, "': design '", design_of[[m]], "' is not a one-sided formula like ~treatment"))
  design_vars <- union(design_vars, all.vars(f))
}
if (length(design_vars) == 0) die("no design variables found in the design formula(s)")

# exactly one contrast kind per row
kind <- ifelse(!is_blank(comp$factor) | !is_blank(comp$numerator) | !is_blank(comp$denominator), "contrast", "")
kind_coef <- !is_blank(comp$coef)
kind_add <- !is_blank(comp$add) | !is_blank(comp$subtract)
n_kinds <- (kind == "contrast") + kind_coef + kind_add
if (any(n_kinds != 1)) {
  die(paste0(
    "comparisons table: each row needs exactly one of (factor+numerator+denominator), coef, or add[+subtract]; check row(s) ",
    paste(which(n_kinds != 1), collapse = ",")
  ))
}
bad_ct <- which(kind == "contrast" & (is_blank(comp$factor) | is_blank(comp$numerator) | is_blank(comp$denominator)))
if (length(bad_ct)) die(paste0("comparisons table: factor, numerator and denominator are all required (row(s) ", paste(bad_ct, collapse = ","), ")"))
bad_add <- which(is_blank(comp$add) & !is_blank(comp$subtract))
if (length(bad_add)) die(paste0("comparisons table: 'subtract' without 'add' (row(s) ", paste(bad_add, collapse = ","), ")"))
if (any(!is_blank(comp$combined_order) & is.na(suppressWarnings(as.numeric(comp$combined_order))))) {
  die("comparisons table: combined_order must be numeric")
}
for (i in seq_len(nrow(comp))) {
  parse_limits(comp$xlim[i], paste0("comparisons table xlim (row ", i, ")"))
  parse_limits(comp$ylim[i], paste0("comparisons table ylim (row ", i, ")"))
}

# ----------------------------------------------------------------------------- counts ----
if (!file.exists(counts_path)) die(paste0("counts file not found: ", counts_path))
raw <- tryCatch(read.delim(counts_path, header = TRUE, row.names = NULL, check.names = FALSE),
  error = function(e) die(paste0("cannot read counts table: ", conditionMessage(e)))
)
if (ncol(raw) < 2 || nrow(raw) < 1) die("counts table needs a gene column, at least one sample column and at least one gene")
genes <- as.character(raw[[1]])
if (anyNA(genes) || any(!nzchar(genes))) die("counts table: empty gene name(s) in the first column")
if (any(grepl("[,\"]", genes))) die(paste0("counts table: gene name(s) with a comma or a double quote (results csv files are written with quote = FALSE, as the Rmd, and would be corrupted), e.g. ", paste(head(genes[grepl("[,\"]", genes)], 5), collapse = " | ")))
if (anyDuplicated(genes)) die(paste0("counts table: duplicated gene names, e.g. ", paste(head(unique(genes[duplicated(genes)]), 5), collapse = ", ")))
counts <- raw[, -1, drop = FALSE]
rownames(counts) <- genes
colnames(counts) <- strip_ids(colnames(counts))
if (anyDuplicated(colnames(counts))) die(paste0("counts table: duplicated sample columns after suffix stripping: ", paste(unique(colnames(counts)[duplicated(colnames(counts))]), collapse = ", ")))

# ---------------------------------------------------------------------------- metadata ----
meta <- read_csv_clean(meta_path)
if (nrow(meta) == 0) die("the metadata table has no rows")
if (!sample_col %in% names(meta)) die(paste0("metadata has no column '", sample_col, "' (columns: ", paste(names(meta), collapse = ", "), ")"))
# same typing as the legacy read.csv (numbers stay numbers); design variables stay text and are made factors below
conv <- setdiff(names(meta), c(design_vars, sample_col)) # sample ids stay text ("001" must not become 1)
meta[conv] <- lapply(meta[conv], type.convert, as.is = TRUE)
if ("exclude" %in% names(meta) && drop_excluded) {
  ex <- tolower(trimws(as.character(meta$exclude))) %in% c("true", "t", "1", "yes", "y")
  if (any(ex)) {
    message("Dropping ", sum(ex), " excluded sample(s): ", paste(meta[[sample_col]][ex], collapse = ", "))
    meta <- meta[!ex, , drop = FALSE]
  }
}
ids <- strip_ids(as.character(meta[[sample_col]]))
if (any(is_blank(ids))) die("metadata: empty sample id(s)")
if (anyDuplicated(ids)) die(paste0("metadata: duplicated sample ids: ", paste(unique(ids[duplicated(ids)]), collapse = ", ")))
rownames(meta) <- ids
if (nrow(meta) < 2) die("fewer than 2 samples left after reading the metadata")

missing_in_counts <- setdiff(ids, colnames(counts))
if (length(missing_in_counts)) {
  die(paste0(
    "metadata sample(s) not found as count columns: ", paste(missing_in_counts, collapse = ", "),
    "\n  count columns: ", paste(colnames(counts), collapse = ", ")
  ))
}
no_meta <- setdiff(colnames(counts), ids)
if (length(no_meta)) message("Dropping ", length(no_meta), " count column(s) without metadata: ", paste(no_meta, collapse = ", "))
counts <- counts[, ids, drop = FALSE]

# counts must be non-negative whole numbers
if (!all(vapply(counts, is.numeric, logical(1)))) {
  die(paste0("counts table: non-numeric sample column(s): ", paste(names(counts)[!vapply(counts, is.numeric, logical(1))], collapse = ", ")))
}
cm <- as.matrix(counts)
if (anyNA(cm) || any(!is.finite(cm))) die("counts table contains NA / non-finite values")
if (any(cm < 0)) die("counts table contains negative values")
if (any(cm != round(cm))) die("counts table contains non-integer values (DESeq2 needs raw integer counts)")
counts[] <- lapply(counts, function(v) if (is.integer(v)) v else as.integer(round(v)))

# ------------------------------------------------------------------ factors / reference ----
refs <- character(0)
if (nzchar(trimws(ref_arg))) {
  for (kv in trimws(strsplit(ref_arg, ",", fixed = TRUE)[[1]])) {
    if (!nzchar(kv)) next
    p <- strsplit(kv, "=", fixed = TRUE)[[1]]
    if (length(p) != 2 || !nzchar(p[1]) || !nzchar(p[2])) die(paste0("--reference-levels entry '", kv, "' is not factor=level"))
    refs[trimws(p[1])] <- trimws(p[2])
  }
}
for (v in design_vars) {
  if (!v %in% names(meta)) die(paste0("design variable '", v, "' is not a metadata column (columns: ", paste(names(meta), collapse = ", "), ")"))
  if (any(is_blank(meta[[v]]))) die(paste0("design variable '", v, "' has NA / empty values for sample(s): ", paste(rownames(meta)[is_blank(meta[[v]])], collapse = ", ")))
}
factor_vars <- union(design_vars, intersect(names(refs), names(meta)))
for (v in names(refs)) {
  if (!v %in% names(meta)) die(paste0("--reference-levels names '", v, "' which is not a metadata column"))
}
for (v in factor_vars) meta[[v]] <- as.factor(meta[[v]])
for (v in names(refs)) {
  if (!refs[[v]] %in% levels(meta[[v]])) die(paste0("reference level '", refs[[v]], "' not found in '", v, "' (levels: ", paste(levels(meta[[v]]), collapse = ", "), ")"))
  meta[[v]] <- relevel(meta[[v]], refs[[v]])
}
for (v in design_vars) {
  if (nlevels(meta[[v]]) < 2) die(paste0("design variable '", v, "' has a single level (", levels(meta[[v]]), "); DESeq2 needs at least 2"))
}

# comparison rows against the metadata (character contrasts)
for (i in which(kind == "contrast")) {
  f <- comp$factor[i]
  if (!f %in% design_vars) die(paste0("comparison '", comp$comparison[i], "': factor '", f, "' is not a design variable (", paste(design_vars, collapse = ", "), ")"))
  lv <- levels(meta[[f]])
  for (x in c(comp$numerator[i], comp$denominator[i])) {
    if (!x %in% lv) die(paste0("comparison '", comp$comparison[i], "': level '", x, "' not in factor '", f, "' (levels: ", paste(lv, collapse = ", "), ")"))
  }
  if (comp$numerator[i] == comp$denominator[i]) die(paste0("comparison '", comp$comparison[i], "': numerator equals denominator"))
}

# -------------------------------------------------------------------------------- write ----
summary_lines <- c(
  paste0("genes: ", nrow(counts)),
  paste0("samples: ", ncol(counts)),
  paste0("sample order: ", paste(ids, collapse = ", ")),
  "factors (first level = reference):"
)
for (v in factor_vars) {
  how <- if (v %in% names(refs)) "set by --reference-levels" else "R default (first level alphabetically)"
  line <- paste0("  ", v, ": ", paste(levels(meta[[v]]), collapse = ", "), "   [reference = ", levels(meta[[v]])[1], "; ", how, "]")
  summary_lines <- c(summary_lines, line)
  message(line)
}
summary_lines <- c(summary_lines, paste0("models: ", paste(paste0(names(design_of), " (", unlist(design_of), ")"), collapse = "; ")))

saveRDS(counts, file = file.path(out_dir, "counts_edited.rds"))
write.table(data.frame(Gene = rownames(counts), counts, check.names = FALSE),
  file = file.path(out_dir, "counts_edited.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)
saveRDS(meta, file = file.path(out_dir, "sample_metadata_edited.rds"))
write.csv(meta, file.path(out_dir, "sample_metadata_edited.csv"), row.names = TRUE, quote = FALSE)
write.table(comp[, comparison_columns, drop = FALSE],
  file = file.path(out_dir, "comparisons.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)
writeLines(summary_lines, file.path(out_dir, "prepare_summary.txt"))
cat("deg_prepare: ", nrow(counts), " genes x ", ncol(counts), " samples; ", nrow(comp), " comparison(s) in ", length(design_of), " model(s)\n", sep = "")
