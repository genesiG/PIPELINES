#!/bin/bash
# software_versions.sh — exact versions of the tools inside one container image.
# Called by SOFTWARE_VERSIONS (modules/local/software_versions.nf), once per
# distinct container image the pipeline's modules declare:
#   software_versions.sh --image IMAGE --output OUT.versions.tsv
# Must write: OUT.versions.tsv — tab-separated, header "image tool version",
#             one row per tool found in the image (a tool the image lacks is
#             simply absent). bin/write_methods.py reads these.
# Tools probed: cutadapt, mageck, fastqc, multiqc, python3, java, R and the R
# packages the pipeline's R scripts load (keep R_PACKAGES in step with the
# library() calls in bin/*.R).
set -uo pipefail

R_PACKAGES="dplyr tidyr ggplot2 scales"

image="" output=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --image)  image="$2"; shift 2 ;;
        --output) output="$2"; shift 2 ;;
        *)        echo "WARNING: ignoring unknown argument $1" >&2; shift ;;
    esac
done

row() { [[ -n "$2" ]] && printf '%s\t%s\t%s\n' "${image}" "$1" "$2"; }
has() { command -v "$1" >/dev/null 2>&1; }

{
    printf 'image\ttool\tversion\n'
    has cutadapt && row cutadapt "$(cutadapt --version 2>&1 | head -1)"
    # mageck -v prints the bare version number
    has mageck   && row MAGeCK   "$(mageck -v 2>&1 | tail -1)"
    has fastqc   && row FastQC   "$(fastqc --version 2>&1 | head -1 | sed 's/^FastQC v//')"
    has multiqc  && row MultiQC  "$(multiqc --version 2>&1 | head -1 | sed 's/.*version //')"
    has python3  && row Python   "$(python3 --version 2>&1 | sed 's/^Python //')"
    has java     && row Java     "$(java -version 2>&1 | head -1 | sed -E 's/.*version "([^"]+)".*/\1/')"
    if has Rscript; then
        # "4.3.3 (2024-02-29)": the release date gives R's citation year
        row R "$(Rscript -e 'cat(paste0(R.version$major, ".", R.version$minor, " (", R.version$year, "-", R.version$month, "-", R.version$day, ")"))' 2>/dev/null)"
        for p in ${R_PACKAGES}; do
            row "R package ${p}" "$(Rscript -e "cat(as.character(packageVersion('${p}')))" 2>/dev/null)"
        done
    fi
} > "${output}"
cat "${output}"
