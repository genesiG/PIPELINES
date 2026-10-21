# omics

The lab's analysis toolbox as one Docker image: R, Python, NGS alignment and QC tools, MAGeCK, CRISPRO, ProTiler, and curated T2T-CHM13 annotations.

Files in this folder: `Dockerfile`, `wrappers.txt` (command wrappers), `omics-slim-env` (image slimming script), `env/` (conda specs and locks), `DEVIATIONS.txt` and `SOURCES.txt` (machine-readable sources, kept in place), `annotations_README.txt` (T2T annotations), `VERSION`, `make_versions.py`, `smoke_test.sh`.

## What the image is

- Image: `ghcr.io/genesig/omics:<VERSION>`. The current VERSION is 2026.10.0. Each release has a version tag and a `sha-<short>` tag.
- Platform: linux/amd64 only.
- Base: Debian 12.15 slim, pinned by tag and digest. Environments are created with micromamba 2.0.5.
- One conda environment per area under `/opt/envs` (six): `r`, `py`, `ngs`, `mageck`, `crispro`, `protiler`. Env `ngs` also holds what used to be env `qc` and env `ngspy`; see "Why some envs stay separate".
- `/usr/local/bin` holds one wrapper per command. Each wrapper puts its own env first on `PATH` and runs the tool from it.
- Bundled tools and data: `/opt/crispro` (CRISPRO 1.0.2 and annotations), `/opt/protiler` (ProTiler demo and license), `/opt/omics` (T2T annotations, `versions.tsv`, `VERSION`, `smoke_test.sh`, `wrappers.txt`).

## Quick start

Docker (the working directory inside the container is `/work`):

```
docker run --rm -v "$PWD":/work ghcr.io/genesig/omics:2026.10.0 <tool> [args]
docker run --rm -v "$PWD":/work ghcr.io/genesig/omics:2026.10.0 samtools --version
```

HPC with Apptainer (or Singularity, same arguments):

```
apptainer exec docker://ghcr.io/genesig/omics:2026.10.0 <tool> [args]
```

## Which env provides which command

Summary of `wrappers.txt`, grouped by env. Each name in the list is a wrapper in `/usr/local/bin`.

| env | commands |
| --- | --- |
| `r` | `R`, `Rscript`, `pandoc` |
| `py` | `python3`, `python`, `pip`, `pymol` |
| `ngs` | bowtie2 (`bowtie2`, `bowtie2-build`, `bowtie2-align-s`, ... ), hisat2 (`hisat2`, `hisat2-build`, helper `.py` scripts), `samtools`, `bgzip`, `tabix`, `htsfile`, `sambamba`, `bedtools`, BEDOPS (`bedops`, `bedmap`, `sort-bed`, `starch` family, `*2bed` and `*2starch` converters), `picard`, `preseq`, `trimmomatic`, `featureCounts` and other subread tools, `bedGraphToBigWig`, `bigWigToBedGraph`, `SEACR_1.3.sh`, `SEACR_1.3.R`, `seqkit`, deepTools (`bamCoverage`, `bamCompare`, `computeMatrix`, `plotHeatmap`, `plotProfile`, `plotFingerprint`, `multiBamSummary`, ... ), `pyGenomeTracks`, `pgt`, `make_tracks_file`, `macs3`, `fastqc`, `multiqc`, `cutadapt` |
| `mageck` | `mageck`, `mageckGSEA`, `RRA` |
| `crispro` | `crispro` (its wrapper also puts env `r` on `PATH`, because the CRISPRO R scripts call `Rscript`) |
| `protiler` | `protiler` |

Notes:

- `samtools`, `bgzip` and `tabix` come from `ngs` only.
- `cutadapt` comes from `ngs` (5.2). Env `crispro` has its own cutadapt 5.1, which CRISPRO uses internally.
- To call an env's Python or another tool without the wrapper, use the full path, for example `/opt/envs/<env>/bin/python`. Example: `docker run --rm -v "$PWD":/work ghcr.io/genesig/omics:2026.10.0 /opt/envs/crispro/bin/python -c "import crispro"`. If a tool starts helper programs, put `/opt/envs/<env>/bin` first on `PATH`, as the wrappers do.
- The full list of packages and versions is in `/opt/omics/versions.tsv` (columns: env, package, version, build, channel).

## CRISPRO

- Usage: `crispro -h`.
- The annotation table is required: pass `-a /opt/crispro/annotations/CRISPRO.GRCh37.ensembl90.20180901.csv.gz`. It is the GRCh37 (Ensembl 90) table. The image sets `CRISPRO_ANNOTATIONS=/opt/crispro/annotations`.
- Memory: the annotation table has 48,357,758 rows and is 14.5 GB uncompressed. CRISPRO reads it in chunks of 10 million rows (`pd.read_csv` with `chunksize=10**7`) and keeps each chunk in memory as strings. Measured: 1.16 GB per million rows (DataFrame), so one chunk needs about 12 GB and the peak during parsing about 14 to 16 GB. A real CRISPRO run on the full table therefore needs roughly 16 GB of RAM available to the container (an estimate extrapolated from a 1-million-row measurement). On a laptop, raise the Docker Desktop memory to at least 16 GB (Settings > Resources, or `memory=` in `.wslconfig`), or run on the cluster with Apptainer.
- For a single gene, make a one-gene subset of the table first (the smoke test does this for ENSG00000099956; it needs little RAM).
- Scores are passed as a CSV with `guide,score` columns (`--scores`).

```
crispro -a /opt/crispro/annotations/CRISPRO.GRCh37.ensembl90.20180901.csv.gz \
        -i ENSG00000099956 --scores scores.csv -o out
```

## ProTiler

- Usage: `protiler call -i <input> -g <GENE> -s <columns>`.
- Demo data: `/opt/protiler/Demo` (input `sample.txt`, reference output in `Demo/ProTilerOutput`).
- Output goes to `ProTilerOutput/` in the current directory.
- Warnings about fonts on stderr are normal.

```
cd /work && protiler call -i /opt/protiler/Demo/sample.txt -g SMARCB1 -s 9,10,11
```

## Curated T2T-CHM13 annotations

- Location: `/opt/omics/annotations/t2t`.
- Files: `t2t.ncbiRefSeq.curated.gtf.gz` (chrM plus NM_ and NR_ transcripts), `t2t.ncbiRefSeq.curated.proteincoding.gtf.gz` (chrM plus NM_ transcripts), `t2t.ncbiRefSeq.curated.norandom.gtf.gz` (as curated, without chrUn or _random contigs; identical content to the curated file in this build), and `txdb.hs1.sqlite` (TxDb built from the norandom GTF).
- md5 of the uncompressed content: curated `e5106937f58cefaac61483720ae44b74`; proteincoding `39249e1d0b047822d321a22d18c40a8e`; norandom `e5106937f58cefaac61483720ae44b74`. Check with `zcat <file> | md5sum`. The values are also in `checksums.txt` in that folder.
- Source, curation and TxDb details: `/opt/omics/annotations/t2t/README.txt`, a copy of `annotations_README.txt` in this folder.
- Load the TxDb in R:

```
library(AnnotationDbi)
txdb <- AnnotationDbi::loadDb("/opt/omics/annotations/t2t/txdb.hs1.sqlite")
```

## Metadata

- `/opt/omics/versions.tsv`: every package in every env, with version, build and channel. Generated at build by `make_versions.py`. It also lists the R package `breakfast`, which is not a conda package.
- `/opt/omics/VERSION`: the image version.
- `/opt/omics/wrappers.txt`: copy of `wrappers.txt`; the smoke test runs every command listed in it.
- `/usr/local/bin/omics-slim-env`: the slimming script used at build time (source: `docker/omics/omics-slim-env`).
- `/opt/omics/smoke_test.sh`: the smoke test (see below).

## Why some envs stay separate

Env `ngs` now also holds what used to be env `qc` (fastqc 0.12.1, multiqc 1.35, cutadapt 5.2) and env `ngspy` (deeptools 3.5.6, pygenometracks 3.9, macs3 3.0.4). All pinned versions are unchanged by the merge (verified by a solve). Only unpinned libraries moved: `ngs` uses openjdk 23.0.2 (qc had 25.0.2), icu 75.1 and sqlite 3.53.4 (qc had 78.3 and 3.54.1), and the floating dask, bokeh, h5py, pysam, pybigwig and pillow versions of the old ngspy are now those of the `ngs` solve (see `env/ngs.lock.txt`; for example dask 2023.3.0, bokeh 2.4.3, h5py 3.12.1, pysam 0.24.0, pybigwig 0.3.25, pillow 12.0.0).

Merges that were tried and rejected:

- `py` + `ngspy`: the solver finds no solution with the py pins (numpy 2.4.6, matplotlib 3.10.0, scikit-learn 1.9.0, pymol) plus deeptools 3.5.6, pygenometracks 3.9 and macs3 3.0.4 (ngspy needs matplotlib 3.8.4).
- `r` + `py`: this does solve with all pins kept, but it would downgrade unpinned packages in `r` (r-curl 8.0.0 to 7.0.0, r-gert 2.5.0 to 2.1.4, r-textshaping 1.0.4 to 0.4.0, libcurl, harfbuzz, pango, libglib) and in `py` (qt6-main 6.11.2 to 6.9.0, pyside6 6.11.2 to 6.9.0, hdf5, lxml), so they stay separate.
- `mageck` (python 3.9), `crispro` (python 3.10) and `protiler` (python 3.8, R 4.4 plus breakfast) need their own python versions.
- `ngs` keeps R 4.2.3 (SEACR, and preseq with gsl) so it cannot merge with `r` (R 4.4.3).

## Image size and slimming

The script `/usr/local/bin/omics-slim-env` (source: `docker/omics/omics-slim-env`) runs in the same `RUN` as each env creation, so the layers really shrink. It removes:

- `share/doc`, `man`, `info` and `gtk-doc`.
- libLLVM, libclang and libLTO libraries and qdoc, shiboken6 and qmake (build-only; about 600 MB in each Qt env: `py`, `crispro`, `protiler`).
- JDK jmods, include, man and `src.zip`.
- Python test suites (`tests` and `test` directories in site-packages, and the stdlib `test` package).
- R package `doc/` (vignettes) and `tests/` folders.
- In every env except `r`: the conda C/C++/Fortran toolchain (compilers, sysroot, headers under `include/`, `libexec/gcc`, `lib/gcc`, `x86_64-conda-linux-gnu-*` binaries, and `*.a` static libraries outside R packages). It comes in as a run dependency of r-base only for compiling R packages.

Env `r` keeps its toolchain, headers and static libraries (`--keep-toolchain`) so that `install.packages`, devtools and pak can compile from source; only its sysroot locale archive is removed. `conda-meta` is never touched, so `versions.tsv` still lists the removed toolchain packages.

Consequence: in `ngs`, `protiler` and the other envs you cannot compile R or C code. Rebuild a derived image, or use env `r` for that.

Sizes (`du` inside the image, MB; before = previous `omics:test` build, after = this build):

| part | before | after |
| --- | --- | --- |
| `r` | 8453 | 6921 |
| `py` | 2345 | 1341 |
| `ngs` | 2189 | 2359 |
| `mageck` | 372 | 344 |
| `crispro` | 1838 | 973 |
| `protiler` | 3059 | 1199 |
| `/opt/crispro` annotations | 1207 | 1207 (1.2 GB, unchanged) |
| `/opt/omics` (T2T annotations, TxDb, ...) | 177 | 177 |
| total filesystem | about 24 GB | about 14.7 GB (14666 MB) |

The new `ngs` (2359) now includes the old `qc` (1797) and `ngspy` (1337); the combined size before was 5323 MB. The sum of the layer sizes is about 15 GB. `docker images` shows 22 GB for `omics:slim` because of containerd snapshot accounting; do not quote that as the image size. Compressed download: `docker save | gzip -1` gives 6.5 GB (6525003936 bytes); registry layers are gzip-compressed, so expect roughly 6 to 6.5 GB to pull.

BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0 was removed from env `r` (user decision; saves 744 MB). Kept: BSgenome.Hsapiens.UCSC.hs1 (its chromosome names match the curated UCSC-named GTF and TxDb) and BSgenome.Mmusculus.UCSC.mm39.

What could not be slimmed:

- reactome.db 2.0 GB, BSgenome hs1 744 MB and mm39 651 MB, org.Hs.eg.db 366 MB and org.Mm.eg.db 347 MB: data the lab uses.
- The env `r` toolchain, about 0.9 GB, kept on purpose.
- Polars runtime duplicates in `ngs` (`_polars_runtime_32` and `_polars_runtime_compat`, 180 MB each). The compat build is for CPUs without AVX2; it is kept so the image runs everywhere.
- kaleido (193 MB, plotly and multiqc static image export), kept.
- The `/opt/crispro` annotation table (1.2 GB) is required.
- ProTiler StaticFiles (about 250 MB) are required.
- Hard-link sharing of identical files between envs is not possible across Docker layers.

## Reproducibility

- Each env is installed only from its explicit lock file, `env/<name>.lock.txt` (conda explicit format with URLs and md5 sums). The Dockerfile installs these verbatim.
- Never sort the lock files. Their order is the solver's dependency order, and an explicit install runs in file order. Sorting breaks the post-link scripts of the Bioconductor data packages.
- Human-readable specs are in `env/<name>.yml`. Bioconductor packages are not pinned in the specs; the lock records the exact builds.
- Regenerate locks with `bash docker/omics/env/make_locks.sh [names]` (Docker needed; Git Bash on Windows works). `--verify [names]` installs each existing lock in a fresh container. Review the lock diff before a release.
- The base image is pinned by digest in the Dockerfile.
- Non-conda sources (CRISPRO, breakfast, ProTiler) are pinned by sha256 in `SOURCES.txt`.
- The build needs network access: conda channels, bioconductor.org (Bioconductor data packages are downloaded at install and md5-checked by their post-link scripts, not by the lock), GitLab, CRAN, PyPI, GitHub and UCSC.

## Build locally and run the smoke test

```
docker build --platform linux/amd64 -t omics:slim docker/omics
docker run --rm --platform linux/amd64 omics:slim /opt/omics/smoke_test.sh
```

To test the release name locally: `docker tag omics:slim ghcr.io/genesig/omics:2026.10.0`.

- The image is about 15 GB uncompressed (6.5 GB compressed). The build needs roughly 45 to 50 GB of free disk (the env layers are large and temporary package caches are used).
- micromamba can abort on Docker Desktop with "Unexpected error 9 on netlink descriptor". The Dockerfile retries each env three times.
- The smoke test prints `ok` or `FAIL` per check and a failure count, and exits 1 if any check failed. It covers: tool versions; loading of every installed R package (503 packages); the hs1 and mm39 BSgenomes (and that the T2T BSgenome is absent); the T2T TxDb and md5 sums; python imports in `py`, `ngs`, `crispro`, `mageck` and `protiler`; pymol; every command in `/opt/omics/wrappers.txt` runs `--version`, `--help`, `-h` or no arguments without loader or import errors; a CRISPRO one-gene run; the ProTiler demo compared with the reference output; and the license files. Last result: 160 ok, 0 failures. `hisat2-inspect` is skipped (see Known limitations); all other wrapped commands are checked.
- Run the smoke test inside the container; it will not run on Windows.

## Release

1. Bump `docker/omics/VERSION`. Also set the `OMICS_VERSION` default in the Dockerfile to the same value (the Dockerfile comment asks for this).
2. Merge to `main`.
3. `.github/workflows/omics-image.yml` runs when `docker/omics/**` or the workflow changes on `main` (it can also be started by hand). It reads VERSION and fails if `ghcr.io/genesig/omics:<VERSION>` already exists, because VERSION tags are immutable. The job frees runner disk space first (free-disk-space action, including docker images) and fails early if less than 45 GB are free on `/`. It uses the buildx docker driver, so the image is built straight into the local image store and not stored twice, builds with `load: true`, runs the smoke test in the container, and pushes exactly two tags: `ghcr.io/genesig/omics:<VERSION>` and `ghcr.io/genesig/omics:sha-<short sha>`.

Note: the disk-space step has not been run on GitHub yet (actionlint 1.7.12 passes).

## Deviations from the lab HPC specs

Source: `DEVIATIONS.txt` (tab-separated; kept in place as the machine-readable source). Rows are reproduced from that file.

| env | package | wanted | got | why |
| --- | --- | --- | --- | --- |
| env | package | wanted | got | why |
| r | bioconductor-\* (all) | Bioc 3.18 versions in chipseq spec (e.g. deseq2 1.42.0, edger 4.0.16, clusterprofiler 4.10.0, csaw 1.36.0, genomicranges 1.54.1, rtracklayer 1.62.0, org.hs.eg.db 3.18.0, biostrings 2.70.1) | Bioc 3.20 for R 4.4 (e.g. deseq2 1.46.0, edger 4.4.0, limma 3.62.1, clusterprofiler 4.14.0, csaw 1.40.0, genomicranges 1.58.0, rtracklayer 1.66.0, org.hs.eg.db 3.20.0, biostrings 2.74.0, jaspar2024 0.99.6) | policy: all Bioconductor from the release matching R 4.4.3; Bioconductor packages are not pinned (lock records exact builds). bioconda carries Bioc 3.20 as r44 builds |
| r | icu | 73.2 (chipseq) | 75.1 (floating, lock has icu 75.1) | r-base 4.4.3 requires icu >=75.1; pin dropped |
| r | r-base | 4.4.3 (crispro; chipseq/crisprscreen 4.3.3; rnaseq 4.4.2) | 4.4.3 | as policy (latest) |
| r | python | 3.12.13 (none in lab R env) | 3.12.13 | added pin: without it the solver floated to python 3.15; keeps all envs on the lab python |
| r | r-gridextra | 2.3 | 2.3.1 | conda spec "=2.3" matches 2.3.x; solver took the patch release |
| r | r-gglorenz | listed in request | NOT AVAILABLE on conda-forge/bioconda | not installed (no pip/source install by rule) |
| r | r-ggpubr | unpinned | 0.6.3 (pinned) | latest conda ggpubr 1.0.0 needs dplyr >=1.2.0 and tidyr >=1.3.2, but specs fix dplyr 1.1.4 / tidyr 1.3.1 (conda metadata does not state this, loading fails) |
| r | r-rstatix | unpinned | 0.7.3 (pinned) | rstatix 1.x goes with the ggpubr 1.0.0 / dplyr 1.2 generation; kept consistent with ggpubr 0.6.3 |
| r | r-ggbreak | unpinned | 0.1.7 (pinned) | ggbreak 0.2.2 needs ggplot2 >=4.0.0, specs fix ggplot2 3.5.2 |
| r | r-ggformula | unpinned (pulled in by r-mess) | 1.0.1 (pinned) | ggformula 2.0.0 needs ggplot2 >=4.0. NOTE 1.0.1 DESCRIPTION also says ggplot2 >=4.0 but it loads fine with 3.5.2 (loadNamespace OK); no older ggformula exists for R 4.4 on conda. Only r-mess uses it; plotting functions of ggformula untested |
| r | r-dbplyr | unpinned | 2.5.2 (pinned) | dbplyr 2.6.0 needs dplyr >=1.2.0 (missing in conda metadata): its load failed and broke the post-link install of JASPAR2024 (data package) with dplyr 1.1.4 |
| r | r-lme4 / r-reformulas | unpinned | lme4 1.1_38, reformulas 0.4.4 (lme4 pinned, reformulas added) | solver picked lme4 1.1-35.5, built against Matrix ABI 1 (Matrix is ABI 2: lme4 warns); lmerTest 3.2-1 needs reformulas but the conda recipe does not list it |
| r | r-matrix | - | 1.7_4 | floats (ships with R 4.4.3 solve); recorded for information |
| r | BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0 | installed in the first build | removed (not in the image) | user decision: saves 744 MB. Kept BSgenome.Hsapiens.UCSC.hs1 (its chromosome names match the curated UCSC-named GTF and TxDb) and BSgenome.Mmusculus.UCSC.mm39. The lock line was deleted by hand (no re-solve), nothing else in r.lock.txt changed |
| r | r-plotly, r-magick, r-devtools, r-pak, r-reticulate, r-showtext, r-png, r-cairo, r-readxl, r-biocmanager, r-purrr, r-readr, r-tidyverse, r-dplyr, r-ggplot2 ... | spec pins | all kept exactly (e.g. r-cairo 1.6_2, r-png 0.1_8, r-showtext 0.9_7, r-pak 0.9.0, r-purrr 1.1.0) | no deviation: the conda builds exist for r44 and solve together |
| r | bioconductor annotation/data packages (GenomeInfoDbData, GO.db, reactome.db, org.\*.eg.db, JASPAR2024, BSgenome.\*, TxDb.\*) | - | installed by post-link scripts | bioconda data packages download their tarballs from bioconductor.org at INSTALL time (md5-checked by the script, not part of the explicit lock) -> the image build needs network and the result depends on bioconductor.org serving the same files. The post-link scripts also need yq/R/curl already linked, hence the lock order (see below) |
| all | lock file order | "sorted urls" | solver (dependency) order, NOT sorted | an explicit install runs in file order; with a sorted file bioconductor-\* data packages were linked before yq and their post-link scripts failed silently (GenomeInfoDbData, GO.db, org.\*.eg.db missing, most Bioconductor packages could not load) |
| py | biopython | 1.78 | 1.81 | 1.78..1.80 cannot be installed with python 3.12 + numpy 2.4.6 (1.78 needs python <=3.8 / numpy <2); 1.81 is the nearest release that solves (smoke-tested: PairwiseAligner OK, Bio.pairwise2 deprecated warning) |
| py | numpy | 2.4.6 (chipseq) vs 1.26.4 (crispro) | 2.4.6 | max across specs |
| py | scipy | 1.15.3 (crispro) | 1.18.1 | unpinned: spec value is a minimum ("1.15.3 or newer" allowed by request) |
| py | pymol-open-source | pymol 3.1.5.1 (crispro, schrodinger) / pymol-bundle 2.6.2 (crisprscreen, schrodinger) | pymol-open-source 3.1.0 | schrodinger channel not used (licensed bundle); open-source build from conda-forge, 3.1.5.1 not published there |
| py | python-docx | 1.2.0 via pip | 1.2.0 via conda-forge | no pip needed |
| py | matplotlib | 3.10.0 (crispro) vs 3.9.x/3.8.4 | 3.10.0 | max across specs |
| ngs | matplotlib | 3.8.4 (chipseq) | 3.8.4 (floating, chosen by deeptools/pygenometracks) | formerly env ngspy, merged into ngs; not pinned in ngs; matches the chipseq spec by coincidence |
| ngs | numpy | 2.4.6 | 2.5.4 (floating) | formerly env ngspy, merged into ngs; not pinned (numpy has no own pin in ngs); py env has 2.4.6 |
| ngs | r-base (for SEACR) | 4.4.3 | 4.2.3 | seacr 1.3 depends on r-base; preseq 2.0.2 needs gsl 1.16 which conflicts with r-base 4.4.3 (needs gsl 2.7). Options: keep as is (SEACR runs on R 4.2.3, only base R used by the wrapper) or drop preseq=2.0.2 |
| ngs | hisat2 | 2.2.1 (source build in rnaseq spec) | 2.2.1 conda package (bioconda) | conda build used instead of a source build |
| ngs | samtools / htslib | 1.23 (chipseq), 1.22.1 (rnaseq) | 1.23.1 (pinned =1.23) | conda spec "=1.23" matches 1.23.x; pin kept in the merged env (formerly env ngspy also carried the pin; without it ngspy took samtools/htslib 1.24) |
| ngs | python | - | 3.12.13 (pinned) | added pin; otherwise python 3.14 floated (formerly env qc: python 3.13 floated; qc is merged into ngs) |
| ngs | openjdk | unpinned | 23.0.2 (lock has openjdk 23.0.2) | unpinned; the old qc env had 25.0.2, the merged ngs env gets 23.0.2 from the picard/trimmomatic solve; fastqc 0.12.1 works with it |
| ngs | multiqc / cutadapt / fastqc | multiqc 1.35, cutadapt 5.1 (pip in specs), fastqc 0.12.1 | multiqc 1.35, cutadapt 5.2 (conda), fastqc 0.12.1 | formerly env qc, merged into ngs; cutadapt 5.2 from conda instead of pip 5.1, as in the specs otherwise |
| mageck | python / numpy | 3.9.23 | 3.9.23 (numpy 2.0.2, scipy 1.13.1) | no deviation |
| all | channels | defaults/r/schrodinger in specs | conda-forge, bioconda only (strict priority, nodefaults) | requested |
| all | pip | python-docx==1.2.0, cutadapt==5.1 (pip in specs) | none | both available on conda; explicit locks cannot hold pip packages. crispro 1.0.2 (setup.py install from /home/genesi/crispro) is NOT part of this image phase |
| crispro | python | 3.12.10 (user env) | 3.10.x (3.10.22 in lock) | crispro 1.0.2 imports Bio.SubsMat (removed in biopython 1.80); biopython 1.79 is the newest with SubsMat and has no python 3.11/3.12 conda build; python 3.10 is the version for which the whole functional combination (numpy 1.26.4, pandas 1.5.3, matplotlib 3.5.3, seaborn 0.11.2, scipy 1.15.2, pymol-open-source) exists |
| crispro | biopython | 1.78 | 1.79 | 1.78 conda builds exist only for python <=3.9 (incompatible with the python 3.10 combination); 1.79 is the nearest build for python 3.10 and still has Bio.SubsMat |
| crispro | scipy | 1.15.3 | 1.15.2 | scipy 1.15.3 is not on conda-forge linux-64 (1.15.2 then 1.16.0) |
| crispro | pymol-open-source | 3.1.5.1 (user env) | 3.1.0 | newest pymol-open-source on conda-forge is 3.1.0 (3.1.5.1 is not published there) |
| crispro | seaborn / pandas / matplotlib | 0.13.2 / 2.2.3 / 3.10.0 (user env) | 0.11.2 / 1.5.3 / 3.5.3 | crispro 1.0.2 code needs the old seaborn API: plot_annotations calls sns.violinplot(x, y) positionally, which seaborn >=0.12 rejects (TypeError). This combination ran the crispro README example end-to-end |
| protiler | r-base | >=3.5 (spec) | 4.4.3 | breakfast 0.1.0 (pure R, imports plyr) builds and runs on R 4.4.3; latest minor of lab R policy |
| protiler | seaborn, matplotlib | unpinned | 0.13.2, 3.7.3 | spec unpinned; matplotlib 3.7.3 is the newest for python 3.8 |

## Sources

Source: `SOURCES.txt` (tab-separated; kept in place as the machine-readable source). Rows are reproduced from that file.

| name | url | sha256 | commit/tag | license |
| --- | --- | --- | --- | --- |
| crispro-1.0.2 (GitLab archive, 1311136119 bytes, includes 1264882976-byte annotations/CRISPRO.GRCh37.ensembl90.20180901.csv.gz) | https://gitlab.com/bauerlab/crispro/-/archive/v1.0.2/crispro-v1.0.2.tar.gz | fe0bcd511b2de71072e177c5800f67c9f62323c5a0cfb6f6d3b0dbbcf75a06b1 | tag v1.0.2 = bda15dcdcbd655da492e54d0562e696629cea860 | AGPL-3.0 (crispro-v1.0.2/LICENSE, sha256 8f9bd745c4b785380e0c859ae67f6eef268df7c14cdd50412d43775116f5e9b5) |
| breakfast-0.1.0 (R source, CRAN archive) | https://cran.r-project.org/src/contrib/Archive/breakfast/breakfast_0.1.0.tar.gz | deb8b5951e1184a79ed5caba9802a57041a0e944f519ec1e884e134ca5c82c93 | CRAN 0.1.0 (imports plyr; NeedsCompilation: no) | GPL (DESCRIPTION: "GPL") |
| protiler-1.0.2 (PyPI sdist, 43996654 bytes, real StaticFiles, no LFS pointers, no license file) | https://files.pythonhosted.org/packages/ca/e3/a5da3443ed521a3883d17702ca22e32ea15a5ff7cb9f291191ba22bfdc70/protiler-1.0.2.tar.gz | 3bde8fd7c4691d1dd783817b2bb11957197406efd0efbbab74b9b2072c7f4d92 | PyPI release 1.0.2 as installed (pip show and dist-info say 1.0.2; the command script bin/protiler contains a hard-coded __version__ = "1.1.0" and ProTiler/profiler.py says "1.0.0", strings inside upstream code; the CLI has no --version flag) | GPL-3.0 (declared in GitHub repo only; PyPI metadata says UNKNOWN) |
| protiler-github-edb475e (Demo data + LISENCE.txt only; its StaticFiles \*_dic are Git LFS pointers, do not use for install) | https://github.com/MDhewei/protiler/archive/edb475e88d90374e36a8e9cd956d1d583a7df501.tar.gz | f02ee6235c034970598f3619053219aa2f3586b06cea7673c14213bf7968b3c0 | edb475e88d90374e36a8e9cd956d1d583a7df501 | GPL-3.0 (LISENCE.txt in repo root, sha256 3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986) |

## Licensing

- The image is a collection. Each bundled tool keeps its own license, and this folder does not relicense any of them.
- CRISPRO 1.0.2: AGPL-3.0, at `/opt/crispro/LICENSE`.
- ProTiler: GPL-3.0, at `/opt/protiler/LICENSE`. The file is `LISENCE.txt` from the GitHub repository; the PyPI release ships no license file.
- breakfast 0.1.0 (R, CRAN): GPL.
- Other packages: their own licenses, shown in the conda package metadata (`/opt/envs/<env>/conda-meta/*.json`).

## Known limitations

- The `ngs` env carries R 4.2.3 for SEACR (seacr 1.3 depends on r-base; preseq 2.0.2 needs gsl 1.16, which conflicts with r-base 4.4.3). All other R work uses env `r` (R 4.4.3).
- `r-gglorenz` is not available on conda-forge or bioconda, so it is not installed.
- `ggformula` 1.0.1 is pinned and loads, but its plotting functions are untested.
- `Bio.pairwise2` is deprecated in biopython 1.81 (env `py`) and prints a deprecation warning.
- ProTiler: the image installs the PyPI sdist protiler 1.0.2 (sha256 3bde8fd7...); `pip show protiler` and the dist-info say 1.0.2. The installed command script `/opt/envs/protiler/bin/protiler` contains a hard-coded `__version__ = "1.1.0"` and `ProTiler/profiler.py` says "1.0.0"; these are strings inside the upstream code. The CLI has no `--version` flag (`protiler --version` errors), so there is no runtime way to print a version. Treat 1.0.2 (the PyPI release, as installed) as the version. The GitHub commit edb475e88d90374e36a8e9cd956d1d583a7df501 is used only for the Demo data and the license text (the PyPI package has none).
- `hisat2-inspect` (the python wrapper script of bioconda hisat2 2.2.1) fails with `ModuleNotFoundError: No module named 'imp'` on python 3.12 (upstream bug; same in the previous image). Use `hisat2-inspect-s` or `hisat2-inspect-l`. The smoke test skips it.
- Env `ngs`, `protiler` and the other envs (not `r`) have no compilers: you cannot compile R or C code there (see "Image size and slimming").
- BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0 is not in the image; use BSgenome.Hsapiens.UCSC.hs1.
- The image has only been tested on this lab's Docker Desktop.
- The other version differences from the lab specs (for example crispro on Python 3.10) are listed in the Deviations table above.
