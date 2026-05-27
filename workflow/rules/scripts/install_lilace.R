# Set up logging
log_file <- snakemake@log[[1]]
output_file <- snakemake@output[[1]]

con <- file(log_file, open = "wt")
sink(con, append = FALSE, split = FALSE, type = "message")
sink(con, append = FALSE, split = FALSE, type = "output")

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
  if (snakemake@config[["lilace_local"]]) {
    message("Checking local R environment for Lilace.")

    if (!require(cmdstanr)) {
      stop("cmdstanr not installed. Please check installation.")
    }
    message("Local cmdstanr found.")

    tryCatch(
      {
        install_cmdstan(version = "2.39.0", overwrite = FALSE)
      },
      warning = function(w) {
        message("cmdstan already installed. Not reinstalling.")
      }
    )

    if (!require(lilace)) {
      stop("Lilace not available in local R environment. Please check installation.")
    }
    message("Local lilace found.")
  } else {
    install.packages("renv", repos = c("https://cloud.r-project.org"))

    library("renv")
    renv::restore()

    library("cmdstanr")

    tryCatch(
      check_cmdstan_toolchain(),
      error = function(e) {
        stop(paste(
          "CmdStan toolchain check failed. Ensure a C++ compiler (gcc/clang) and",
          "cmake <3.25 are available in the conda environment.\nOriginal error:", conditionMessage(e)
        ))
      }
    )

    tryCatch(
      {
        install_cmdstan(version = "2.39.0", overwrite = FALSE)
      },
      warning = function(w) {
        message("cmdstan already installed. Not reinstalling.")
      }
    )

    if (!require(lilace)) {
      stop("Lilace not available after renv restore. Check renv.lock and renv setup.")
    }
  }

  marker <- file(output_file, open = "wt")
  writeLines("Lilace installed.", marker)
  close(marker)
}

tryCatch(
  main(),
  error = function(err) {
    log_and_rethrow("install_lilace", err)
  }
)
