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
        r-base=4.5 r-renv r-nloptr nlopt libxml2 'cmake<3.25' zlib compilers \
        scipy pytables statsmodels matplotlib \
    && /opt/dumpling/bin/pip install --no-cache-dir \
        multiqc-dumpling 'enrich2>=2.0.0' \
    && mamba clean --all -y \
    && rm /tmp/dumpling_env.yaml

# Default PATH so all bundled tools resolve without explicit activation.
ENV PATH=/opt/dumpling/bin:$PATH

# The default Rosace backend is part of the image, rather than being compiled
# on first use. Runtime containers are commonly launched with
# `--user "$(id -u):$(id -g)"`; that user cannot write to the conda library
# under /opt, and source-package configure checks are particularly fragile
# under cross-architecture Docker emulation. Keep renv's library and cache at
# fixed, image-owned paths so run_rosace.R resolves the prebuilt packages no
# matter where the user's project is mounted.
ENV RENV_PATHS_LIBRARY=/opt/dumpling/renv/library \
    RENV_PATHS_CACHE=/opt/dumpling/renv/cache \
    DUMPLING_PREINSTALLED_ROSACE=1
COPY renv.lock /opt/dumpling/renv-project/renv.lock
COPY renv/settings.json /opt/dumpling/renv-project/renv/settings.json
RUN Rscript -e \
    "library_path <- renv::paths\$library(project='/opt/dumpling/renv-project'); \
     renv::restore(project='/opt/dumpling/renv-project', \
                   lockfile='/opt/dumpling/renv-project/renv.lock', \
                   library=library_path, \
                   packages=c('rosace', 'purrr'), prompt=FALSE); \
     .libPaths(c(library_path, .libPaths())); \
     stopifnot(requireNamespace('rosace', quietly=TRUE), \
               requireNamespace('cmdstanr', quietly=TRUE), \
               requireNamespace('purrr', quietly=TRUE))" \
    && chmod -R a+rX /opt/dumpling/renv

# Pre-compile CmdStan to /opt/cmdstan with the locked cmdstanr version from the
# fixed library. At runtime the scoring install scripts call
# install_cmdstan(overwrite=FALSE); CMDSTAN points at this build, so they skip
# rebuilding it.
ENV CMDSTAN=/opt/cmdstan/cmdstan-2.39.0
RUN Rscript -e \
    ".libPaths(c(renv::paths\$library(project='/opt/dumpling/renv-project'), .libPaths())); \
     dir.create('/opt/cmdstan', recursive=TRUE); \
     cmdstanr::install_cmdstan(version='2.39.0', dir='/opt/cmdstan', overwrite=FALSE)"

# Default working directory; users mount their experiment data here.
WORKDIR /workdir
