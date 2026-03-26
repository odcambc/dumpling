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

# Determine the directory of this script robustly
# commandArgs() works with Rscript; sys.frame(1)$ofile works with source()
this_script <- tryCatch(
  normalizePath(sys.frame(1)$ofile),
  error = function(e) {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0) {
      normalizePath(sub("^--file=", "", file_arg[1]))
    } else {
      # Fallback: assume running from repo root
      normalizePath("tests/r/testthat.R")
    }
  }
)
this_dir <- dirname(this_script)

# Source helper functions that define the functions we want to test
source(file.path(this_dir, "testthat/helper-functions.R"))

# Run all tests
cat("Running R unit tests for Dumpling pipeline...\n")
cat("==============================================\n\n")

test_results <- test_dir(
  file.path(this_dir, "testthat"),
  reporter = "summary"
)

# Exit with appropriate code
quit(status = as.integer(any(as.data.frame(test_results)$failed > 0)))
