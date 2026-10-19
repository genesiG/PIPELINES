# CRISPR_SCREEN — the one Docker image every step of the pipeline runs in
# (Docker version for Windows / Mac). Built automatically by
# `run_crisprscreen setup` (or the first run); nobody needs to run this by hand.
#
# Same software as the cluster version, where each tool had its own container:
#   rocker/tidyverse:4.3.3   R 4.3.3 + tidyverse (cluster: rocker/tidyverse:4.3)
#   cutadapt 5.2             cluster: quay.io/biocontainers/cutadapt:5.2--py39hbcbf7aa_0
#   MAGeCK 0.5.9.5           cluster: quay.io/biocontainers/mageck:0.5.9.5--py39h2de1943_8
#   FastQC 0.12.1            cluster: quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0
#   MultiQC 1.34             cluster: quay.io/biocontainers/multiqc:1.34--pyhdfd78af_0
#   Nextflow 26.04.6         the version the cluster's launcher pins
# Each Bioconda tool gets its own environment, installed from the exact
# package list in docker/env/<tool>.explicit.txt (an explicit, fully pinned conda
# lock file copied from that BioContainers image's own records by
# docker/env/make_locks.sh), so every computer that builds this image gets
# the same package builds as the cluster's containers. The pipeline's Python helper scripts run with the cutadapt
# environment's python3, as on the cluster (where they ran in that image).
#
# Always built for linux/amd64 (also on Apple Silicon Macs, which run it via
# Rosetta): identical binaries on every computer, and several of these
# packages have no arm64 build.
#
# The pipeline code itself (main.nf, modules/, bin/) is NOT copied in: it is
# mounted read-only at /pipeline at run time, so a pipeline update never
# needs a rebuild. Only a change to this file or docker/env/ does (the
# launcher notices and rebuilds by itself).

FROM rocker/tidyverse:4.3.3

ARG MICROMAMBA_VERSION=2.0.5-0
ARG NXF_VER=26.04.6

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8 \
    MAMBA_ROOT_PREFIX=/opt/micromamba

# Java for Nextflow and FastQC; procps (ps) for Nextflow's task metrics.
RUN apt-get update \
 && apt-get install -y --no-install-recommends openjdk-17-jre-headless procps curl ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# micromamba: a single static binary, pinned.
RUN curl -fsSL -o /usr/local/bin/micromamba \
      "https://github.com/mamba-org/micromamba-releases/releases/download/${MICROMAMBA_VERSION}/micromamba-linux-64" \
 && chmod +x /usr/local/bin/micromamba

# One environment per tool, from the pinned lock files (only those: editing
# make_locks.sh must not invalidate this layer). micromamba occasionally aborts
# on Docker Desktop ("Unexpected error 9 on netlink descriptor"): retry.
COPY env/*.explicit.txt /tmp/env/
RUN for tool in cutadapt mageck fastqc multiqc; do \
        ok=0; \
        for attempt in 1 2 3; do \
            micromamba create -y -p /opt/envs/${tool} -f /tmp/env/${tool}.explicit.txt && { ok=1; break; }; \
            echo "micromamba failed for ${tool} (attempt ${attempt} of 3)" >&2; rm -rf /opt/envs/${tool}; sleep 5; \
        done; \
        [ "${ok}" = 1 ] || exit 1; \
    done \
 && micromamba clean -a -y \
 && rm -rf /tmp/env

# Put exactly the tools the pipeline calls on PATH (not whole environments,
# so R stays rocker's and python3 is the cutadapt environment's 3.9). Each tool
# is a wrapper that runs it with its own environment first on PATH, as in its
# BioContainer: mageck test needs the environment's RRA, FastQC its java.
RUN for tool in cutadapt mageck fastqc multiqc; do \
        printf '#!/bin/sh\n# %s with its own environment first on PATH, as in its BioContainer.\nexport PATH="/opt/envs/%s/bin:$PATH"\nexec /opt/envs/%s/bin/%s "$@"\n' \
            "${tool}" "${tool}" "${tool}" "${tool}" > /usr/local/bin/${tool} \
     && chmod 755 /usr/local/bin/${tool} || exit 1; \
    done \
 && ln -s /opt/envs/cutadapt/bin/python3  /usr/local/bin/python3 \
 && ln -s /opt/envs/cutadapt/bin/python3  /usr/local/bin/python

# Nextflow, pinned, with its runtime downloaded now so runs work offline.
ENV NXF_VER=${NXF_VER} \
    NXF_HOME=/opt/nextflow \
    NXF_OFFLINE=true \
    NXF_ANSI_LOG=false
RUN curl -fsSL -o /usr/local/bin/nextflow "https://github.com/nextflow-io/nextflow/releases/download/v${NXF_VER}/nextflow" \
 && chmod 755 /usr/local/bin/nextflow \
 && NXF_OFFLINE=false nextflow -version \
 && chmod -R a+rwX /opt/nextflow

# Quick self-check: the build fails here if any tool is missing.
RUN cutadapt --version && mageck -v && fastqc --version && multiqc --version \
 && test -x /opt/envs/mageck/bin/RRA && test -x /opt/envs/fastqc/bin/java \
 && python3 --version && Rscript -e 'library(dplyr); library(tidyr); library(ggplot2); library(scales)' \
 && nextflow -version

WORKDIR /project
CMD ["bash"]
