# install_rosace_aa.R
# Installs the rosace-aa R package (pimentellab/rosace-aa) into either a
# locally-managed R library (when config[["rosace_aa_local"]] = TRUE) or
# the conda-controlled renv project library.
#
# Mirrors install_rosace.R's shape — same renv::restore() bootstrap for
# the base R dependency set — but additionally calls renv::install() to
# pull rosace-aa from GitHub at a pinned SHA. The SHA pin keeps repeat
# installs bit-reproducible without requiring rosace-aa to live in
# renv.lock (the upstream repo has no tagged releases, so we pin the
# main-branch HEAD as of 2026-05-26).

log_file <- snakemake@log[[1]]
output_file <- snakemake@output[[1]]

con <- file(log_file, open = "wt")
sink(con, append = FALSE, split = FALSE, type = "message")
sink(con, append = FALSE, split = FALSE, type = "output")

# Pinned SHA from pimentellab/rosace-aa main branch HEAD on 2026-05-26.
# Bump this when validating against a newer rosace-aa release.
ROSACE_AA_SHA <- "c0fda386d43f7bbede4211b4de171265f43305c9"


log_and_rethrow <- function(script_name, err) {
  message(sprintf("%s failed: %s", script_name, conditionMessage(err)))
  calls <- vapply(
    sys.calls(),
    function(call) paste(deparse(call), collapse = " "),
    character(1)
  )
  if (length(calls) > 0) {
    message("Call stack:")
    message(paste(calls, collapse = "\n"))
  }
  stop(err)
}


main <- function() {
  if (snakemake@config[["rosace_aa_local"]]) {
    message("Checking local R environment for rosace-aa.")

    if (require(cmdstanr)) {
      message("Local cmdstanr found.")
    } else {
      stop("cmdstanr not installed. Please check installation.")
    }

    tryCatch(
      install_cmdstan(version = "2.39.0", overwrite = FALSE),
      warning = function(w) {
        message("cmdstan already installed. Not reinstalling.")
      }
    )

    if (require(rosaceAA)) {
      message("Local rosaceAA found.")
    } else {
      stop("rosaceAA not available in local R environment. Please check installation.")
    }
  } else {
    install.packages("renv", repos = c("https://cloud.r-project.org"))
    library("renv")

    # Restore base R deps from the repo-level lockfile (cmdstanr, impute,
    # the common-to-both-backends pieces). rosace-aa itself isn't in the
    # lockfile yet — install it explicitly below at the pinned SHA.
    renv::restore()

    library("cmdstanr")

    tryCatch(
      check_cmdstan_toolchain(),
      error = function(e) {
        stop(paste(
          "CmdStan toolchain check failed. Ensure a C++ compiler (gcc/clang) and",
          "cmake <3.25 are available in the conda environment.\nOriginal error:",
          conditionMessage(e)
        ))
      }
    )

    tryCatch(
      install_cmdstan(version = "2.39.0", overwrite = FALSE),
      warning = function(w) {
        message("cmdstan already installed. Not reinstalling.")
      }
    )

    # Pin to a specific SHA so repeat installs of this conda env reproduce
    # bit-for-bit. The Bioconductor `impute` dep is already in renv.lock
    # so it gets pulled in by renv::restore() above.
    message(sprintf("Installing rosace-aa pinned at SHA %s", ROSACE_AA_SHA))
    renv::install(
      sprintf("pimentellab/rosace-aa@%s", ROSACE_AA_SHA),
      prompt = FALSE
    )

    # Sanity-check: package must be loadable as `rosaceAA` (camelCase, per
    # the DESCRIPTION's Package: field on the upstream repo).
    if (!require(rosaceAA)) {
      stop("rosaceAA install completed but the package is not loadable.")
    }
  }

  con <- file(output_file, open = "wt")
  writeLines(sprintf("rosace-aa installed (pinned at %s).", ROSACE_AA_SHA), con)
  close(con)
}

tryCatch(
  main(),
  error = function(err) {
    log_and_rethrow("install_rosace_aa", err)
  }
)
