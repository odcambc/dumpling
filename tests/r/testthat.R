#!/usr/bin/env Rscript
# testthat.R - Main test runner for R unit tests
#
# Run from repository root:
#   Rscript tests/r/testthat.R
#
# Or run specific test file:
#   Rscript -e "testthat::test_file('tests/r/testthat/test-parse_hgvs.R')"

library(testthat)

# Load required libraries
suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(readr)
})

# Source the functions we want to test
# Note: We extract just the pure functions from run_rosace.R to avoid
# snakemake dependency during testing

# Get the repository root
repo_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."))
scripts_dir <- file.path(repo_root, "workflow/rules/scripts")

# Source helper functions that define the functions we want to test
source(file.path(dirname(sys.frame(1)$ofile), "testthat/helper-functions.R"))

# Run all tests
cat("Running R unit tests for Dumpling pipeline...\n")
cat("==============================================\n\n")

test_results <- test_dir(
  file.path(dirname(sys.frame(1)$ofile), "testthat"),
  reporter = "summary"
)

# Exit with appropriate code
quit(status = as.integer(any(as.data.frame(test_results)$failed > 0)))
