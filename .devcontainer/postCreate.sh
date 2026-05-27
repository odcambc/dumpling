#!/usr/bin/env bash
# Bootstrap the dumpling devcontainer:
#   1. Create / update the canonical conda env from dumpling_env.yaml.
#   2. (--with-dev-deps only) install the python dev requirements
#      (pytest + pytest-mock) into dumpling_env.
#   3. Install the system-R packages the unit-test suite reads.
#   4. Pre-warm all three scoring backends (rosace, lilace, rosace_aa)
#      via the existing snakemake install rules. After this, a user
#      can switch scoring_backend in their config without waiting on
#      the per-run install rule. First run takes 10-15 min (cmdstan
#      compile is the long pole); subsequent rebuilds reuse the renv
#      project library and cmdstan install.
#
# Both .devcontainer/devcontainer.json and
# .devcontainer/claude-sandbox/devcontainer.json call this script via
# postCreateCommand. Keeping the logic in one place avoids the drift
# between the two configs that was building up.

set -euo pipefail

cd "$(dirname "$0")/.."

WITH_DEV_DEPS=0
for arg in "$@"; do
  case "$arg" in
    --with-dev-deps) WITH_DEV_DEPS=1 ;;
    *)
      echo "[postCreate] unknown argument: $arg" >&2
      exit 2
      ;;
  esac
done

echo "[postCreate] creating/updating dumpling_env conda environment"
conda env create -f dumpling_env.yaml || conda env update -f dumpling_env.yaml

if [[ "$WITH_DEV_DEPS" -eq 1 ]]; then
  echo "[postCreate] installing python dev requirements into dumpling_env"
  /opt/conda/envs/dumpling_env/bin/pip install -r requirements-dev.txt
fi

echo "[postCreate] installing system-R packages for the unit-test suite"
Rscript -e 'install.packages(c("dplyr", "readr", "stringr", "testthat", "withr"), repos = "https://cloud.r-project.org")'

echo "[postCreate] pre-warming scoring backends (rosace, lilace, rosace_aa)"
# Uses the existing install rules' codepath (renv::restore + cmdstanr +
# install_github for rosace-aa) so this stays in lockstep with the
# install scripts users actually run. The marker files at
# results/example_experiment/{rosace,lilace,rosace_aa}/*_installed.txt
# persist after build so subsequent snakemake invocations see the
# install as a cache-hit.
#
# --cores 1 is deliberate. Each install rule's R script begins with
# install.packages("renv", ...) targeting the shared system library at
# /usr/local/lib/R/site-library. With parallel cores the rules race on
# the 00LOCK-renv directory and corrupt each other's installs (the
# loser sees a partial renv with a missing .rdb file and aborts). The
# wall-clock cost of serial execution is small — renv install is fast;
# the CmdStan compile dominates and only runs once anyway.
#
# Non-fatal on failure: scoring-backend installs reach out to GitHub
# (install_github for rosace-aa), download CmdStan, and compile a Stan
# toolchain — any of which can fail transiently or in restricted
# network environments. A failed warmup leaves the devcontainer
# perfectly usable for everything else; the user just hits a real
# (~5-10 min) install on first `scoring_backend=<...>` run instead of
# a cache-hit. Print a clear retry hint in that case.
if conda run --no-capture-output -n dumpling_env \
    snakemake -s workflow/Snakefile --configfile config/example.yaml \
    --cores 1 \
    results/example_experiment/rosace/rosace_installed.txt \
    results/example_experiment/lilace/lilace_installed.txt \
    results/example_experiment/rosace_aa/rosace_aa_installed.txt; then
  echo "[postCreate] done — all three scoring backends available"
else
  echo "[postCreate] WARNING: scoring-backend pre-warm failed (exit $?)." >&2
  echo "[postCreate] The devcontainer is usable; rosace/lilace/rosace_aa" >&2
  echo "[postCreate] rules will fire renv::restore() on first scoring run" >&2
  echo "[postCreate] instead of being instant. Retry the warmup manually:" >&2
  echo "[postCreate]   conda run -n dumpling_env snakemake -s workflow/Snakefile \\" >&2
  echo "[postCreate]     --configfile config/example.yaml --cores 1 \\" >&2
  echo "[postCreate]     results/example_experiment/rosace/rosace_installed.txt \\" >&2
  echo "[postCreate]     results/example_experiment/lilace/lilace_installed.txt \\" >&2
  echo "[postCreate]     results/example_experiment/rosace_aa/rosace_aa_installed.txt" >&2
fi
