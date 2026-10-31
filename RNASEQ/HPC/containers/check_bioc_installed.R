# check_bioc_installed.R - run right after env r is created (Dockerfile).
#
# Bioconda's Bioconductor data packages (reactome.db, GO.db, org.*.eg.db,
# BSgenome.*, TxDb.*, JASPAR2024, GenomeInfoDbData, ...) download their content
# from bioconductor.org in a post-link script. When that download fails, conda
# still records the package as installed and the build goes on: the R package
# is simply missing (seen with reactome.db). This script fails (exit 1) when any
# bioconductor-* conda package has no matching installed R package, so the
# Dockerfile retries the environment instead of shipping an incomplete image.

prefix <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(prefix)) prefix <- "/opt/envs/r"

meta <- list.files(file.path(prefix, "conda-meta"), pattern = "^bioconductor-.*\\.json$")
# conda-meta file names are <name>-<version>-<build>.json; conda names are lowercase.
conda_pkgs <- sub("-[^-]+-[^-]+\\.json$", "", meta)
conda_pkgs <- sub("^bioconductor-", "", conda_pkgs)
# bioconductor-data-packages is the conda helper the post-link scripts use, not an R package.
conda_pkgs <- setdiff(conda_pkgs, "data-packages")

lib <- file.path(prefix, "lib", "R", "library")
r_pkgs <- tolower(rownames(installed.packages(lib.loc = lib)))

missing <- setdiff(conda_pkgs, r_pkgs)
cat(sprintf("Bioconductor packages: %d in conda-meta, %d missing from %s\n",
            length(conda_pkgs), length(missing), lib))
if (length(missing) > 0) {
  cat("Missing (post-link download failed?):", paste(missing, collapse = ", "), "\n")
  quit(status = 1)
}
