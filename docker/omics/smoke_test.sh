#!/usr/bin/env bash
# Smoke test for the omics container. Run inside the Linux container, not on Windows.
# Prints "ok" or "FAIL" per check, then a failure count. Exit 1 if any check failed.
set -u

failures=0

# check "label" command... : runs command with stdout discarded, keeps stderr
# in a temp file and shows it only on failure.
check() {
    local label="$1"
    shift
    local errf
    errf=$(mktemp)
    if "$@" >/dev/null 2>"$errf"; then
        echo "ok   $label"
    else
        echo "FAIL $label"
        sed 's/^/     /' "$errf"
        failures=$((failures + 1))
    fi
    rm -f "$errf"
}

# md5 of the gunzipped content of a file must equal the expected value.
md5_matches() {
    local file="$1" expected="$2" actual
    actual=$(zcat "$file" | md5sum)
    actual=${actual%% *}
    if [ "$actual" = "$expected" ]; then
        return 0
    fi
    echo "expected $expected" >&2
    echo "actual   $actual" >&2
    return 1
}

# check_md5 file expected
check_md5() {
    check "md5 $(basename "$1")" md5_matches "$1" "$2"
}

echo "== tool versions =="
check "R --version" R --version
check "Rscript R.version" Rscript -e 'cat(R.version.string)'
check "python3 --version" python3 --version
check "samtools --version" samtools --version
check "bedtools --version" bedtools --version
check "bowtie2 --version" bowtie2 --version
check "hisat2 --version" hisat2 --version
check "bamCoverage --version" bamCoverage --version
check "macs3 --version" macs3 --version
check "fastqc --version" fastqc --version
check "multiqc --version" multiqc --version
check "cutadapt --version" cutadapt --version
check "mageck -h" mageck -h
check "featureCounts -v" featureCounts -v
check "seqkit version" seqkit version

echo "== R packages =="
check "R core packages" Rscript -e 'library(DESeq2); library(edgeR); library(csaw); library(GenomicRanges); library(clusterProfiler); library(ggplot2)'
check "T2T TxDb loads" Rscript -e 'library(AnnotationDbi); library(GenomicFeatures); db <- AnnotationDbi::loadDb("/opt/omics/annotations/t2t/txdb.hs1.sqlite"); stopifnot(length(GenomicFeatures::genes(db)) > 1000)'
check "all installed R packages load" Rscript -e 'pk <- rownames(installed.packages()); bad <- character(0); for (p in pk) { ok <- suppressWarnings(suppressMessages(requireNamespace(p, quietly=TRUE))); if (!ok) bad <- c(bad, p) }; if (length(bad)) { message("failed to load: ", paste(bad, collapse=" ")); quit(status=1) }'
check "BSgenome hs1 loads" Rscript -e 'library(BSgenome.Hsapiens.UCSC.hs1)'
check "BSgenome mm39 loads" Rscript -e 'library(BSgenome.Mmusculus.UCSC.mm39)'
check "BSgenome T2T CHM13v2.0 not installed" Rscript -e 'stopifnot(!requireNamespace("BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0", quietly=TRUE))'

echo "== python imports =="
check "py imports" /opt/envs/py/bin/python -c "import numpy,scipy,pandas,matplotlib,seaborn,sklearn,joblib,Bio,matplotlib_venn,upsetplot,docx,pymol"
check "ngs imports" /opt/envs/ngs/bin/python -c "import deeptools,pygenometracks,MACS3,pysam,pyBigWig,multiqc,cutadapt"
check "crispro imports" /opt/envs/crispro/bin/python -c "import crispro,Bio.SubsMat,seaborn,pymol,scipy,pandas"
check "mageck imports" /opt/envs/mageck/bin/python -c "import mageck,numpy,scipy"
check "protiler imports" /opt/envs/protiler/bin/python -c "import ProTiler,sklearn,seaborn,statsmodels"
check "pymol -cq -d quit" pymol -cq -d quit

echo "== every wrapped command runs =="
# For each executable in wrappers.txt (env, executable): the wrapper must exist and be
# executable, and one of --version, --help, -h, no arguments must exit with status < 124
# or exactly 255 (usage exit of subread and UCSC tools; no timeout, signal, 126 or 127)
# and print no loader or python import errors.
# The first such attempt is the one judged. pymol is covered by the check above.
wrapper_runs() {
    local exe="$1" w="/usr/local/bin/$1" out rc a found=0
    if [ ! -x "$w" ]; then
        echo "$w missing or not executable" >&2
        return 1
    fi
    out=$(mktemp)
    for a in --version --help -h ""; do
        if [ -n "$a" ]; then
            timeout 60 "$w" "$a" </dev/null >"$out" 2>&1
        else
            timeout 60 "$w" </dev/null >"$out" 2>&1
        fi
        rc=$?
        if [ "$rc" -lt 124 ] || [ "$rc" -eq 255 ]; then
            found=1
            break
        fi
    done
    if [ "$found" = 0 ]; then
        echo "no attempt (--version, --help, -h, no arguments) exited below 124; last status $rc" >&2
        rm -f "$out"
        return 1
    fi
    if grep -q -e "error while loading shared libraries" -e "cannot open shared object file" \
            -e "ModuleNotFoundError" -e "ImportError" -e "Traceback (most recent call last)" "$out"; then
        echo "loader or import error in the output (attempt '$a', status $rc):" >&2
        head -n 5 "$out" >&2
        rm -f "$out"
        return 1
    fi
    rm -f "$out"
    return 0
}

while read -r we exe rest; do
    [ "$exe" = pymol ] && continue
    # bioconda hisat2 2.2.1: the python wrapper hisat2-inspect does "import imp", which
    # python 3.12 no longer has (upstream bug); hisat2-inspect-s and -l work and are checked.
    if [ "$exe" = hisat2-inspect ]; then
        echo "skip hisat2-inspect (python 3.12 has no imp; use hisat2-inspect-s)"
        continue
    fi
    if wrapper_runs "$exe" 2>"${TMPDIR:-/tmp}/wrapper_err.$$"; then
        echo "ok   $exe"
    else
        echo "FAIL $exe"
        sed 's/^/     /' "${TMPDIR:-/tmp}/wrapper_err.$$"
        failures=$((failures + 1))
    fi
done < <(grep -v '^[[:space:]]*\(#\|$\)' /opt/omics/wrappers.txt)
rm -f "${TMPDIR:-/tmp}/wrapper_err.$$"

echo "== annotation checksums =="
check_md5 /opt/omics/annotations/t2t/t2t.ncbiRefSeq.curated.gtf.gz e5106937f58cefaac61483720ae44b74
check_md5 /opt/omics/annotations/t2t/t2t.ncbiRefSeq.curated.proteincoding.gtf.gz 39249e1d0b047822d321a22d18c40a8e
check_md5 /opt/omics/annotations/t2t/t2t.ncbiRefSeq.curated.norandom.gtf.gz e5106937f58cefaac61483720ae44b74

echo "== crispro and protiler =="
check "crispro -h" crispro -h
check "protiler -h" protiler -h
check "protiler env breakfast" /opt/envs/protiler/bin/Rscript -e 'library(breakfast)'

echo "== metadata and licences =="
check "versions.tsv non-empty" test -s /opt/omics/versions.tsv
check "wrappers.txt present" test -s /opt/omics/wrappers.txt
check "VERSION non-empty" test -s /opt/omics/VERSION
check "crispro LICENSE non-empty" test -s /opt/crispro/LICENSE
check "protiler LICENSE non-empty" test -s /opt/protiler/LICENSE

echo "== demos =="
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# ProTiler demo: the new call must give the same segments as the repo's reference
# (AA.start, AA.end, n, m, length; reference columns 1,2,3,4,6). Font warnings go to a log.
protiler_demo() {
    ( cd "$tmp" && protiler call -i /opt/protiler/Demo/sample.txt -g SMARCB1 -s 9,10,11 ) 2>"$tmp/protiler.log" \
        || { cat "$tmp/protiler.log" >&2; return 1; }
    local out="$tmp/ProTilerOutput/SMARCB1_CKHS_Regions.csv"
    local ref=/opt/protiler/Demo/ProTilerOutput/SMARCB1_Segments.csv
    awk -F, 'NR>1 {print $3","$4","$5","$6","$7}' "$out" > "$tmp/protiler_new.txt"
    awk -F, 'NR>1 {print $1","$2","$3","$4","$6}' "$ref" > "$tmp/protiler_ref.txt"
    test -s "$tmp/protiler_new.txt" || { echo "no rows in $out" >&2; return 1; }
    diff "$tmp/protiler_ref.txt" "$tmp/protiler_new.txt" >&2
}

# crispro example on one gene (a one-gene subset of the annotation table; the full
# table needs roughly 14-16 GB RAM). Every output listed below must exist and be non-empty.
crispro_example() {
    (
        set -o pipefail
        cd "$tmp" &&
        zcat /opt/crispro/annotations/CRISPRO.GRCh37.ensembl90.20180901.csv.gz | awk 'NR==1 || /ENSG00000099956/' | gzip > sub.csv.gz &&
        zcat sub.csv.gz | awk -F, 'NR>1 {print $10}' | sort -u | awk 'BEGIN{print "guide,score"} {printf "%s,%.3f\n", $1, sin(NR/7.0)}' > scores.csv &&
        crispro -a sub.csv.gz -i ENSG00000099956 --scores scores.csv -o out
    )
}

crispro_outputs() {
    local f missing=0
    for f in out/crispro_scores_annotated-loess.csv out/crispro_guide_scores.csv \
             out/linear_tracks/crispro_tracks.pdf out/score_summary/crispro_stats_score.csv; do
        if [ ! -s "$tmp/$f" ]; then
            echo "missing or empty: $tmp/$f" >&2
            missing=1
        fi
    done
    local p found=0
    for p in "$tmp"/out/colored_protein_sequences/*.pdf; do
        if [ -s "$p" ]; then
            found=1
            break
        fi
    done
    if [ "$found" = 0 ]; then
        echo "no non-empty PDF in $tmp/out/colored_protein_sequences" >&2
        missing=1
    fi
    return "$missing"
}

check "ProTiler demo SMARCB1 matches reference" protiler_demo
check "crispro one-gene example runs" crispro_example
check "crispro outputs present and non-empty" crispro_outputs

echo "== summary =="
echo "$failures failures"
if [ "$failures" -gt 0 ]; then
    exit 1
fi
exit 0
