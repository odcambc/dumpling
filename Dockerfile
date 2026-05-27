FROM condaforge/miniforge3:latest

# Hand-written from scratch (2026-05-27). Replaces an earlier
# `snakemake --containerize` output that was inflating the image with three
# duplicate R envs (rosace/lilace/rosace_aa share an identical system layer)
# and per-rule conda isolation that the container already provides. See
# tasks.md "container branch" notes for the design discussion.
#
# Strategy: one mega-env at /opt/dumpling. Snakemake + all bioinformatics
# tools + the R/Stan stack + Enrich2 live here. Users run
# `snakemake --use-singularity` (apptainer/singularity transparently pulls
# and converts the same OCI image) without `--use-conda`. The rule files
# keep their `conda:` directives — those are still load-bearing for the
# non-container path, and ignored inside the container.
#
# Enrich2 v2.0+ is Python 3 and PyPI-installable, so it folds into the
# mega-env alongside dumpling_env's python=3.13 — no second env needed.

# Build toolchain for source-package compiles. cmdstanr's CRAN deps
# (cli, rlang, jsonlite, ...) install from source via R CMD INSTALL ->
# make -> system gcc. Conda's r-base 4.5 looks for `make` and
# `aarch64-conda-linux-gnu-cc`; system gcc covers the former via PATH,
# and the `compilers` conda meta-package supplies the conda-prefixed
# names. ca-certificates is needed for HTTPS-fetched CRAN deps.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential ca-certificates curl \
    && rm -rf /var/lib/apt/lists/*

# Single mega-env: dumpling_env.yaml (snakemake + bbtools + minimap2 + gatk +
# samtools + fastqc + pandas + biopython + multiqc + mavehgvs) PLUS the R
# system layer that rosace/lilace/rosace_aa.yaml share (r-base 4.5 +
# r-nloptr + nlopt + libxml2 + cmake<3.25 + zlib) PLUS conda's compilers
# meta-package (so cmdstanr's CRAN deps source-compile against the toolchain
# r-base 4.5's Makeconf was built against).
#
# multiqc-dumpling is the multiqc-baseline.yaml's distinctive extra; folding
# it in here means we don't need a separate multiqc-baseline env in the
# container path.
COPY dumpling_env.yaml /tmp/dumpling_env.yaml
RUN mamba env create --prefix /opt/dumpling --file /tmp/dumpling_env.yaml \
    && mamba install --prefix /opt/dumpling -c conda-forge -c bioconda -y \
        r-base=4.5 r-nloptr nlopt libxml2 'cmake<3.25' zlib compilers \
        scipy pytables statsmodels matplotlib \
    && /opt/dumpling/bin/pip install --no-cache-dir \
        multiqc-dumpling 'enrich2>=2.0.0' \
    && mamba clean --all -y \
    && rm /tmp/dumpling_env.yaml

# Default PATH so all bundled tools resolve without explicit activation.
ENV PATH=/opt/dumpling/bin:$PATH

# Pre-compile CmdStan to /opt/cmdstan. Bootstrap cmdstanr from r-universe
# just for this step; at runtime, install_rosace.R / install_lilace.R /
# install_rosace_aa.R will renv::restore() the locked cmdstanr version
# (0.9.0.9000 from a specific GitHub SHA, per renv.lock) and call
# install_cmdstan(version="2.39.0", overwrite=FALSE). Because CMDSTAN
# env var points at the already-built binary, those calls skip the rebuild.
# This cures the "two weekends installing CmdStan" failure mode at image
# build time.
ENV CMDSTAN=/opt/cmdstan/cmdstan-2.39.0
RUN Rscript -e \
    "install.packages('cmdstanr', repos=c('https://stan-dev.r-universe.dev', 'https://cloud.r-project.org')); \
     dir.create('/opt/cmdstan', recursive=TRUE); \
     cmdstanr::install_cmdstan(version='2.39.0', dir='/opt/cmdstan', overwrite=FALSE)"

# Default working directory; users mount their experiment data here.
WORKDIR /workdir
