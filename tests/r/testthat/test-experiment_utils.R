# test-experiment_utils.R
# Unit tests for experiment utility functions

library(testthat)

context("Experiment Utilities")

# =============================================================================
# Tests for get_conditions()
# =============================================================================

test_that("get_conditions excludes baseline condition", {
  experiment_df <- data.frame(
    sample = c("s1", "s2", "s3", "s4"),
    condition = c("cond_A", "cond_A", "baseline", "baseline"),
    replicate = c(1, 1, 1, 1),
    time = c(0, 1, 0, 1)
  )

  result <- get_conditions(experiment_df, "baseline")

  expect_true("cond_A" %in% result)
  expect_false("baseline" %in% result)
})

test_that("get_conditions returns all non-baseline conditions", {
  experiment_df <- data.frame(
    sample = c("s1", "s2", "s3", "s4", "s5", "s6"),
    condition = c("cond_A", "cond_A", "cond_B", "cond_B", "baseline", "baseline"),
    replicate = c(1, 1, 1, 1, 1, 1),
    time = c(0, 1, 0, 1, 0, 1)
  )

  result <- get_conditions(experiment_df, "baseline")

  expect_equal(length(result), 2)
  expect_true("cond_A" %in% result)
  expect_true("cond_B" %in% result)
})

test_that("get_conditions handles custom baseline name", {
  experiment_df <- data.frame(
    sample = c("s1", "s2", "s3", "s4"),
    condition = c("treatment", "treatment", "control", "control"),
    replicate = c(1, 1, 1, 1),
    time = c(0, 1, 0, 1)
  )

  result <- get_conditions(experiment_df, "control")

  expect_true("treatment" %in% result)
  expect_false("control" %in% result)
})

# =============================================================================
# Tests for get_samples_for_condition()
# =============================================================================

test_that("get_samples_for_condition returns correct samples", {
  experiment_df <- data.frame(
    sample = c("s1", "s2", "s3", "s4"),
    condition = c("cond_A", "cond_A", "baseline", "baseline"),
    replicate = c(1, 1, 1, 1),
    time = c(0, 1, 0, 1)
  )

  result <- get_samples_for_condition(experiment_df, "cond_A")

  expect_equal(length(result), 2)
  expect_true("s1" %in% result)
  expect_true("s2" %in% result)
})

test_that("get_samples_for_condition returns empty for nonexistent condition", {
  experiment_df <- data.frame(
    sample = c("s1", "s2"),
    condition = c("cond_A", "cond_A"),
    replicate = c(1, 1),
    time = c(0, 1)
  )

  result <- get_samples_for_condition(experiment_df, "nonexistent")

  expect_equal(length(result), 0)
})

# =============================================================================
# Tests for validate_experiment()
# =============================================================================

test_that("validate_experiment passes for valid experiment", {
  experiment_df <- data.frame(
    sample = c("s1", "s2", "s3", "s4"),
    condition = c("cond_A", "cond_A", "baseline", "baseline"),
    replicate = c(1, 1, 1, 1),
    time = c(0, 1, 0, 1)
  )

  result <- validate_experiment(experiment_df)

  expect_true(result)
})

test_that("validate_experiment fails for missing columns", {
  experiment_df <- data.frame(
    sample = c("s1", "s2"),
    condition = c("cond_A", "cond_A")
    # Missing: replicate, time
  )

  result <- validate_experiment(experiment_df)

  expect_true(is.character(result))
  expect_true(grepl("Missing required columns", result))
})

test_that("validate_experiment fails for missing T0", {
  experiment_df <- data.frame(
    sample = c("s1", "s2"),
    condition = c("cond_A", "cond_A"),
    replicate = c(1, 1),
    time = c(1, 2)  # No T0!
  )

  result <- validate_experiment(experiment_df)

  expect_true(is.character(result))
  expect_true(grepl("missing T0", result))
})

test_that("validate_experiment checks each condition for T0", {
  experiment_df <- data.frame(
    sample = c("s1", "s2", "s3", "s4"),
    condition = c("cond_A", "cond_A", "cond_B", "cond_B"),
    replicate = c(1, 1, 1, 1),
    time = c(0, 1, 1, 2)  # cond_B missing T0
  )

  result <- validate_experiment(experiment_df)

  expect_true(is.character(result))
  expect_true(grepl("cond_B", result))
})
