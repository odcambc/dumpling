# test-build_counts.R
# Unit tests for build_counts_for_replicate()

library(testthat)
library(withr)

context("Build Counts for Replicate")

# Helper: write count TSV files into the expected pipeline directory tree
# relative to the current working directory.
setup_enrich_files <- function(experiment_name, samples_and_data) {
  dir_path <- file.path("results", experiment_name, "processed_counts", "enrich_format")
  dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  for (sample in names(samples_and_data)) {
    readr::write_tsv(
      samples_and_data[[sample]],
      file.path(dir_path, paste0(sample, ".tsv"))
    )
  }
}

# =============================================================================
# Test 1: Single time point, no tile column
# =============================================================================

test_that("single time point returns hgvs and c_0 columns", {
  tmp_dir <- withr::local_tempdir()
  withr::local_dir(tmp_dir)

  setup_enrich_files("exp1", list(
    sample1 = tibble::tibble(hgvs = c("A1V", "A2T"), count = c(10L, 20L))
  ))

  df_subset <- tibble::tibble(
    sample    = "sample1",
    condition = "cond_A",
    replicate = 1L,
    time      = 0L
  )

  result <- build_counts_for_replicate(df_subset, "exp1")

  expect_true("hgvs" %in% names(result))
  expect_true("c_0" %in% names(result))
  expect_equal(nrow(result), 2)
  expect_equal(sort(result$hgvs), c("A1V", "A2T"))
})

# =============================================================================
# Test 2: Two time points — HGVS union via full_join
# =============================================================================

test_that("two time points produce c_0 and c_1 columns with correct HGVS union", {
  tmp_dir <- withr::local_tempdir()
  withr::local_dir(tmp_dir)

  # sample_t0 has variants A and B; sample_t1 has variants A and C
  setup_enrich_files("exp2", list(
    sample_t0 = tibble::tibble(hgvs = c("A1V", "A2T"),       count = c(10L, 20L)),
    sample_t1 = tibble::tibble(hgvs = c("A1V", "A3S"),       count = c(15L, 25L))
  ))

  df_subset <- tibble::tibble(
    sample    = c("sample_t0", "sample_t1"),
    condition = c("cond_A", "cond_A"),
    replicate = c(1L, 1L),
    time      = c(0L, 1L)
  )

  result <- build_counts_for_replicate(df_subset, "exp2")

  # All three variants should appear (full join)
  expect_true(all(c("hgvs", "c_0", "c_1") %in% names(result)))
  expect_equal(nrow(result), 3)
  expect_true(all(c("A1V", "A2T", "A3S") %in% result$hgvs))

  # A1V appears in both; A2T only in t0; A3S only in t1
  a1v_row <- result[result$hgvs == "A1V", ]
  expect_equal(a1v_row$c_0, 10L)
  expect_equal(a1v_row$c_1, 15L)

  a2t_row <- result[result$hgvs == "A2T", ]
  expect_equal(a2t_row$c_0, 20L)
  expect_true(is.na(a2t_row$c_1))

  a3s_row <- result[result$hgvs == "A3S", ]
  expect_true(is.na(a3s_row$c_0))
  expect_equal(a3s_row$c_1, 25L)
})

# =============================================================================
# Test 3: Duplicate sample rows for the same time point trigger stop()
# =============================================================================

test_that("two rows with the same time and no tile triggers stop", {
  df_subset <- tibble::tibble(
    sample    = c("sample1", "sample2"),
    condition = c("cond_A", "cond_A"),
    replicate = c(1L, 1L),
    time      = c(0L, 0L)   # both at time 0 — no tile column to disambiguate
  )

  expect_error(
    build_counts_for_replicate(df_subset, "exp3"),
    "More than one row"
  )
})

# =============================================================================
# Test 4: Empty df_subset returns empty tibble without error
# =============================================================================

test_that("empty df_subset returns empty tibble gracefully", {
  df_subset <- tibble::tibble(
    sample    = character(0),
    condition = character(0),
    replicate = integer(0),
    time      = integer(0)
  )

  result <- build_counts_for_replicate(df_subset, "exp4")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
  expect_true("hgvs" %in% names(result))
})
