# Set up logging
# Output will be redicted to the log file
log_file <- snakemake@log[[1]]

output_file <- snakemake@output[[1]]

# Open connection to log file
con <- file(log_file, open = "wt")

sink(con, append = FALSE, split = FALSE, type = "message")
sink(con, append = FALSE, split = FALSE, type = "output")


# If user-managed rosace install is used, check it
if (snakemake@config[["rosace_local"]]) {
  message("Checking local R environment for Rosace.")

  # Check if cmdstanr available
  if (require(cmdstanr)) {
    message("Local cmdstanr found.")
  } else {
    stop("cmdstanr not installed. Please check installation.")
  }

  # Check cmdstanr toolchain
  tryCatch(
    {
      install_cmdstan(overwrite = FALSE)
    },
    warning = function(w) {
      message("cmdstan already installed. Not reinstalling.")
    }
  )

  # Check Rosace install
  if (require(rosace)) {
    message("Local rosace found.")
  } else {
    stop("Rosace not available in local R environment. Please check installation.")
  }

  # If local install is not used, install rosace into conda env through Renv
} else {
  # Install renv
  install.packages("renv", repos = c("https://cloud.r-project.org"))

  library("renv")

  renv::restore()

  library("cmdstanr")


  # Check toolchain before intstalling cmdstan
  tryCatch(
    check_cmdstan_toolchain(),
    error = function(e) {
      stop(paste(
        "CmdStan toolchain check failed. Ensure a C++ compiler (gcc/clang) and",
        "cmake <3.25 are available in the conda environment.\nOriginal error:", conditionMessage(e)
      ))
    }
  )

  # Install cmdstan, catch warning if already installed
  tryCatch(
    {
      install_cmdstan(overwrite = FALSE)
    },
    warning = function(w) {
      message("cmdstan already installed. Not reinstalling.")
    }
  )
}
# Write to output file upon successful installation

con <- file(output_file, open = "wt")

writeLines("Rosace installed.", con)
close(con)
