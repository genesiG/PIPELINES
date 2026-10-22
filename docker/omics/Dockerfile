# omics - the lab's analysis toolbox as one image (R, python, NGS, QC, CRISPR tools).
# See README.md. Built for linux/amd64 only. Environments come first and are
# installed ONLY from the pinned env/<name>.lock.txt files, so editing anything
# else (wrappers, smoke test, docs) never rebuilds them.

# Debian 12.15 slim, pinned by tag and by digest of the multi-arch index.
FROM debian:12.15-slim@sha256:7ff6ba31527b9e336162fa061746c366c4b2032f4d2cff2078bea43f78efe3f8

ARG MICROMAMBA_VERSION=2.0.5-0

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    MAMBA_ROOT_PREFIX=/opt/micromamba

# procps (ps; Nextflow-style tools), curl/ca-certificates (downloads), locales,
# archive tools for the data downloads below.
RUN apt-get update \
 && apt-get install -y --no-install-recommends procps curl ca-certificates locales tar gzip bzip2 xz-utils \
 && sed -i 's/^# *en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen \
 && locale-gen \
 && rm -rf /var/lib/apt/lists/*

# micromamba: a single static binary, pinned.
RUN curl -fsSL -o /usr/local/bin/micromamba \
      "https://github.com/mamba-org/micromamba-releases/releases/download/${MICROMAMBA_VERSION}/micromamba-linux-64" \
 && chmod +x /usr/local/bin/micromamba

# Slimming helper for the environments below (see its header).
COPY --chmod=755 omics-slim-env /usr/local/bin/omics-slim-env

# One environment per layer, each from its own pinned lock (the lock files are in
# dependency order, NOT sorted: never sort them). The Bioconductor data packages in
# env r download their tarballs from bioconductor.org at install: network needed.
# micromamba occasionally aborts on Docker Desktop ("Unexpected error 9 on netlink
# descriptor"): retry. Env r is also retried when a Bioconductor data package's
# post-link download failed silently (env/check_bioc_installed.R; seen with reactome.db).
# Each env is slimmed by omics-slim-env (docs, LLVM/Qt build tools, JDK build files,
# test suites, ...) in the SAME RUN as its creation: a later rm would not shrink the
# layer. Env r keeps its compiler toolchain (--keep-toolchain) so install.packages,
# devtools and pak can still compile R packages from source; the others lose it.

COPY env/r.lock.txt env/check_bioc_installed.R /tmp/env/
RUN ok=0;     for attempt in 1 2 3; do         micromamba create -y -p /opt/envs/r -f /tmp/env/r.lock.txt && /opt/envs/r/bin/Rscript /tmp/env/check_bioc_installed.R /opt/envs/r && { ok=1; break; };         echo "micromamba failed for r (attempt ${attempt} of 3)" >&2; rm -rf /opt/envs/r; sleep 5;     done;     [ "${ok}" = 1 ]  && micromamba clean -a -y  && omics-slim-env --keep-toolchain /opt/envs/r && rm -rf /tmp/env

COPY env/py.lock.txt /tmp/env/py.lock.txt
RUN ok=0;     for attempt in 1 2 3; do         micromamba create -y -p /opt/envs/py -f /tmp/env/py.lock.txt && { ok=1; break; };         echo "micromamba failed for py (attempt ${attempt} of 3)" >&2; rm -rf /opt/envs/py; sleep 5;     done;     [ "${ok}" = 1 ]  && micromamba clean -a -y  && omics-slim-env /opt/envs/py && rm -rf /tmp/env

COPY env/ngs.lock.txt /tmp/env/ngs.lock.txt
RUN ok=0;     for attempt in 1 2 3; do         micromamba create -y -p /opt/envs/ngs -f /tmp/env/ngs.lock.txt && { ok=1; break; };         echo "micromamba failed for ngs (attempt ${attempt} of 3)" >&2; rm -rf /opt/envs/ngs; sleep 5;     done;     [ "${ok}" = 1 ]  && micromamba clean -a -y  && omics-slim-env /opt/envs/ngs && rm -rf /tmp/env

COPY env/mageck.lock.txt /tmp/env/mageck.lock.txt
RUN ok=0;     for attempt in 1 2 3; do         micromamba create -y -p /opt/envs/mageck -f /tmp/env/mageck.lock.txt && { ok=1; break; };         echo "micromamba failed for mageck (attempt ${attempt} of 3)" >&2; rm -rf /opt/envs/mageck; sleep 5;     done;     [ "${ok}" = 1 ]  && micromamba clean -a -y  && omics-slim-env /opt/envs/mageck && rm -rf /tmp/env

COPY env/crispro.lock.txt /tmp/env/crispro.lock.txt
RUN ok=0;     for attempt in 1 2 3; do         micromamba create -y -p /opt/envs/crispro -f /tmp/env/crispro.lock.txt && { ok=1; break; };         echo "micromamba failed for crispro (attempt ${attempt} of 3)" >&2; rm -rf /opt/envs/crispro; sleep 5;     done;     [ "${ok}" = 1 ]  && micromamba clean -a -y  && omics-slim-env /opt/envs/crispro && rm -rf /tmp/env

COPY env/protiler.lock.txt /tmp/env/protiler.lock.txt
RUN ok=0;     for attempt in 1 2 3; do         micromamba create -y -p /opt/envs/protiler -f /tmp/env/protiler.lock.txt && { ok=1; break; };         echo "micromamba failed for protiler (attempt ${attempt} of 3)" >&2; rm -rf /opt/envs/protiler; sleep 5;     done;     [ "${ok}" = 1 ]  && micromamba clean -a -y  && omics-slim-env /opt/envs/protiler && rm -rf /tmp/env

# (a) CRISPRO 1.0.2 (AGPL-3.0). The archive includes the 1.26 GB GRCh37 annotation
# table, kept at /opt/crispro/annotations (run crispro with -a); the rest is installed
# into env crispro with pip. The archive and the extracted tree are deleted in the
# same RUN. A full-table run needs roughly 14-16 GB RAM (chunks of 10 million rows, ~1.2 GB per million rows).
RUN set -eux; \
    mkdir -p /tmp/crispro/src /opt/crispro; cd /tmp/crispro; \
    curl -fsSL -o crispro.tar.gz https://gitlab.com/bauerlab/crispro/-/archive/v1.0.2/crispro-v1.0.2.tar.gz; \
    echo "fe0bcd511b2de71072e177c5800f67c9f62323c5a0cfb6f6d3b0dbbcf75a06b1  crispro.tar.gz" | sha256sum -c -; \
    tar -xzf crispro.tar.gz -C src --strip-components=1; rm -f crispro.tar.gz; \
    mv src/annotations /opt/crispro/annotations; \
    cp src/LICENSE /opt/crispro/LICENSE; \
    if [ -d src/example_data ]; then cp -r src/example_data /opt/crispro/example_data; fi; \
    cd src && /opt/envs/crispro/bin/python -m pip install --root-user-action=ignore --no-deps --no-build-isolation --no-cache-dir . ; \
    cd /; rm -rf /tmp/crispro; \
    chmod -R a+rX /opt/crispro

# (b) ProTiler: breakfast 0.1.0 (R, CRAN archive) and protiler 1.0.2 (PyPI) into env
# protiler; the repo's Demo data and its GPL-3 license text (the PyPI package ships none).
RUN set -eux; \
    mkdir -p /tmp/protiler /opt/protiler; cd /tmp/protiler; \
    curl -fsSL -o breakfast_0.1.0.tar.gz https://cran.r-project.org/src/contrib/Archive/breakfast/breakfast_0.1.0.tar.gz; \
    echo "deb8b5951e1184a79ed5caba9802a57041a0e944f519ec1e884e134ca5c82c93  breakfast_0.1.0.tar.gz" | sha256sum -c -; \
    PATH=/opt/envs/protiler/bin:$PATH R CMD INSTALL breakfast_0.1.0.tar.gz; \
    echo "protiler==1.0.2 --hash=sha256:3bde8fd7c4691d1dd783817b2bb11957197406efd0efbbab74b9b2072c7f4d92" > requirements.txt; \
    /opt/envs/protiler/bin/pip install --root-user-action=ignore --no-deps --no-cache-dir --no-build-isolation --require-hashes -r requirements.txt; \
    curl -fsSL -o repo.tar.gz https://github.com/MDhewei/protiler/archive/edb475e88d90374e36a8e9cd956d1d583a7df501.tar.gz; \
    echo "f02ee6235c034970598f3619053219aa2f3586b06cea7673c14213bf7968b3c0  repo.tar.gz" | sha256sum -c -; \
    mkdir repo; tar -xzf repo.tar.gz -C repo --strip-components=1 --exclude='*/StaticFiles'; \
    cp -r repo/Demo /opt/protiler/Demo; cp repo/LISENCE.txt /opt/protiler/LICENSE; \
    echo "3972dc9744f6499f0f9b2dbf76696f2ae7ad8af9b23dde66d6af86c9dfb36986  /opt/protiler/LICENSE" | sha256sum -c -; \
    cd /; rm -rf /tmp/protiler; chmod -R a+rX /opt/protiler

# (c) Curated T2T-CHM13 (hs1) NCBI RefSeq annotations + TxDb. Download, verify,
# curate with the lab's awk one-liners, verify again, build the TxDb with env r,
# keep only the curated files (gzip -n: reproducible). Any md5 mismatch fails the build.
RUN set -eux; \
    d=/opt/omics/annotations/t2t; mkdir -p "${d}" /tmp/ann; cd /tmp/ann; \
    curl -fsSL -o hs1.ncbiRefSeq.gtf.gz https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz; \
    echo "bd7063efcbfed4774d9be42db8173a6e  hs1.ncbiRefSeq.gtf.gz" | md5sum -c -; \
    gzip -dc hs1.ncbiRefSeq.gtf.gz > t2t.ncbiRefSeq.gtf; \
    echo "95c46ca9573b2ed12fd03eec5e949164  t2t.ncbiRefSeq.gtf" | md5sum -c -; \
    awk '$1 == "chrM" || $12 ~ /"N[MR]_/' t2t.ncbiRefSeq.gtf > t2t.ncbiRefSeq.curated.gtf; \
    awk '$1 == "chrM" || $12 ~ /"NM_/' t2t.ncbiRefSeq.gtf > t2t.ncbiRefSeq.curated.proteincoding.gtf; \
    awk '($1 == "chrM" || $12 ~ /\"N[MR]_/) && $1 !~ /chrUn|_random/' t2t.ncbiRefSeq.gtf > t2t.ncbiRefSeq.curated.norandom.gtf; \
    echo "e5106937f58cefaac61483720ae44b74  t2t.ncbiRefSeq.curated.gtf" | md5sum -c -; \
    echo "39249e1d0b047822d321a22d18c40a8e  t2t.ncbiRefSeq.curated.proteincoding.gtf" | md5sum -c -; \
    echo "e5106937f58cefaac61483720ae44b74  t2t.ncbiRefSeq.curated.norandom.gtf" | md5sum -c -; \
    rm -f t2t.ncbiRefSeq.gtf; \
    PATH=/opt/envs/r/bin:$PATH Rscript -e 'suppressMessages({library(txdbmaker); library(AnnotationDbi)}); a <- commandArgs(TRUE); txdb <- txdbmaker::makeTxDbFromGFF(a[1], format="gtf", organism="Homo sapiens", dataSource="ncbiRefSeq excluding random contigs and keeping only NM and NR annotations"); AnnotationDbi::saveDb(txdb, file=a[2])' t2t.ncbiRefSeq.curated.norandom.gtf "${d}/txdb.hs1.sqlite"; \
    { echo "# md5 of the UNCOMPRESSED content (the .gtf.gz files here are gzip -n of these)"; \
      echo "bd7063efcbfed4774d9be42db8173a6e  hs1.ncbiRefSeq.gtf.gz (source download, not kept)"; \
      echo "95c46ca9573b2ed12fd03eec5e949164  t2t.ncbiRefSeq.gtf (source decompressed, not kept)"; \
      md5sum t2t.ncbiRefSeq.curated.gtf t2t.ncbiRefSeq.curated.proteincoding.gtf t2t.ncbiRefSeq.curated.norandom.gtf; \
    } > "${d}/checksums.txt"; \
    gzip -9n t2t.ncbiRefSeq.curated.gtf t2t.ncbiRefSeq.curated.proteincoding.gtf t2t.ncbiRefSeq.curated.norandom.gtf; \
    mv t2t.ncbiRefSeq.curated*.gtf.gz "${d}/"; \
    cd /; rm -rf /tmp/ann; chmod -R a+rX /opt/omics

# (d) Wrappers in /usr/local/bin, generated from wrappers.txt (env, executable).
# Each runs the tool with its own environment first on PATH, so helper programs
# (java, Rscript, RRA, ...) come from the same environment. crispro also gets env r.
# wrappers.txt stays in the image (/opt/omics/wrappers.txt; the smoke test reads it).
COPY wrappers.txt /opt/omics/wrappers.txt
RUN set -eu; \
    mkdir -p /opt/omics; chmod a+r /opt/omics/wrappers.txt; \
    grep -v '^[[:space:]]*\(#\|$\)' /opt/omics/wrappers.txt | while read -r we exe rest; do \
        [ -x "/opt/envs/${we}/bin/${exe}" ] || { echo "wrappers.txt: /opt/envs/${we}/bin/${exe} not found" >&2; exit 1; }; \
        p="/opt/envs/${we}/bin"; \
        [ "${exe}" = crispro ] && p="${p}:/opt/envs/r/bin"; \
        printf '#!/bin/sh\nexport PATH=%s:$PATH\nexec /opt/envs/%s/bin/%s "$@"\n' "${p}" "${we}" "${exe}" > "/usr/local/bin/${exe}"; \
        chmod 755 "/usr/local/bin/${exe}"; \
    done

# crispro-lowmem: hand-written helper (runs crispro on a cached per-gene subset of the annotation table).
COPY --chmod=755 crispro-lowmem /usr/local/bin/crispro-lowmem

# (e) Small files last, so changing them rebuilds only the layers below.
COPY annotations_README.txt /opt/omics/annotations/t2t/README.txt
COPY make_versions.py /tmp/make_versions.py
RUN /opt/envs/py/bin/python /tmp/make_versions.py /opt/omics && rm -f /tmp/make_versions.py
COPY VERSION /opt/omics/VERSION
COPY smoke_test.sh /opt/omics/smoke_test.sh
RUN chmod 755 /opt/omics/smoke_test.sh && chmod 644 /opt/omics/VERSION /opt/omics/annotations/t2t/README.txt

ENV OMICS_ANNOTATIONS=/opt/omics/annotations \
    CRISPRO_ANNOTATIONS=/opt/crispro/annotations

# Same as VERSION (the release workflow passes it; keep the default in sync).
ARG OMICS_VERSION=2026.10.1
LABEL org.opencontainers.image.source="https://github.com/genesiG/PIPELINES" \
      org.opencontainers.image.title="omics" \
      org.opencontainers.image.description="Lab omics toolbox: R 4.4 + Bioconductor 3.20, python, NGS (bowtie2, hisat2, samtools, deepTools, MACS3, SEACR...), QC (FastQC, MultiQC, cutadapt), MAGeCK, CRISPRO, ProTiler, curated T2T-CHM13 annotations" \
      org.opencontainers.image.version="${OMICS_VERSION}" \
      org.opencontainers.image.licenses="Mixed: see /opt/crispro/LICENSE (AGPL-3.0), /opt/protiler/LICENSE (GPL-3.0) and README.md; each bundled tool keeps its own license" \
      org.opencontainers.image.base.name="docker.io/library/debian:12.15-slim"

WORKDIR /work
CMD ["bash"]
