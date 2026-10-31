#!/usr/bin/env Rscript
# txdb_build.R -- build a TxDb (transcript and gene models) from a GTF and save it as one SQLite file.
# Both downstream steps that need gene models (DIFF_ACCESS promoters, DIFF_ANNOTATE peak annotation) read this file, so they
# share one annotation. Port of the legacy "load annotations" chunk (Analysis_Data/csaw/differential_accessibility.Rmd and
# peak_annotation.Rmd): makeTxDbFromGFF(gtf, format = "gtf", organism, dataSource) then AnnotationDbi::saveDb(txdb, file).
# The legacy reused rds/txdb.hs1.sqlite when it existed; here the file is always rebuilt from the GTF (the Nextflow task
# cache covers reruns).
#
# CLI contract (called by the TXDB_BUILD module; --key=value):
#   --gtf=F                    REQUIRED. Gene annotation GTF (format gtf), e.g. the T2T/hs1 NCBI RefSeq curated GTF.
#                              Missing or empty file: exit 3.
#   --out=txdb.sqlite          REQUIRED. Output TxDb SQLite file (written to this path).
#   --species="Homo sapiens"   OPTIONAL, default "Homo sapiens". organism= of makeTxDbFromGFF (kept in the TxDb metadata).
#   --data-source=...          OPTIONAL. dataSource= of makeTxDbFromGFF (free text, kept in the TxDb metadata). Blank = NA.
#
# Outputs (cwd, or the path given to --out):
#   <out>                      TxDb SQLite file (AnnotationDbi::loadDb reads it back; GenomicFeatures genes()/transcripts())
#
# Exit codes: 0 ok; 2 usage error; 3 known data problem (GTF missing or empty, or txdbmaker cannot parse it).

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

allowed <- c("gtf", "out", "species", "data-source")
argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) == 0L || any(argv %in% c("-h", "--help"))) {
  cat("usage: txdb_build.R --gtf=annotation.gtf --out=txdb.sqlite [--species=\"Homo sapiens\"] [--data-source=TEXT]\n",
      file = stderr())
  quit(save = "no", status = 2L)
}
args <- parse_args(argv, allowed)
gtf_file    <- get_req(args, "gtf")
out_file    <- get_req(args, "out")
species     <- get_opt(args, "species", "Homo sapiens")
data_source <- get_opt(args, "data-source", NA)

if (!file.exists(gtf_file)) die(paste0("GTF file not found: ", gtf_file), 3L)
if (file.size(gtf_file) == 0) die(paste0("GTF file is empty: ", gtf_file, " (re-export or re-download the annotation)"), 3L)

suppressPackageStartupMessages({
  library(txdbmaker)
  library(AnnotationDbi)
})

txdb <- tryCatch(
  txdbmaker::makeTxDbFromGFF(gtf_file,
                             format = "gtf",
                             organism = species,
                             dataSource = data_source),
  error = function(e) die(paste0("could not build a TxDb from ", gtf_file, ": ", conditionMessage(e),
                                 " (is it a GTF with exon/transcript/gene features?)"), 3L)
)

# Save for future use (legacy: AnnotationDbi::saveDb(txdb, file = txdb_db))
AnnotationDbi::saveDb(txdb, file = out_file)
cat("Wrote TxDb: ", out_file, " (organism: ", species, ")\n", sep = "")
