# test-resolve_lilace_seed.R
# Unit tests for resolve_lilace_seed() — the Determinism Phase 1 helper
# that decides whether to use a configured Stan seed or auto-generate
# one. Mirrors the function shipped in workflow/rules/scripts/run_lilace.R.

library(testthat)

context("Resolve Lilace Stan Seed")

# =============================================================================
# Test 1: Configured seed wins — repeat runs reproduce bit-identical
# =============================================================================

test_that("configured seed passes through unchanged", {
  out <- suppressMessages(resolve_lilace_seed(
    condition_name = "cond_A",
    configured_seed = 42L
  ))
  expect_equal(out, 42L)
})

test_that("zero is a valid configured seed", {
  # Stan accepts 0; the schema-level test confirms the config schema
  # accepts 0 too. Make sure the resolver doesn't reject it (e.g. via
  # an `if (seed) ...` truthiness check).
  out <- suppressMessages(resolve_lilace_seed(
    condition_name = "cond_A",
    configured_seed = 0L
  ))
  expect_equal(out, 0L)
})

# =============================================================================
# Test 2: Auto-generated path is deterministic given a fixed clock + cond
# =============================================================================

test_that("auto-generated seed is deterministic given fixed clock", {
  out1 <- suppressMessages(resolve_lilace_seed(
    condition_name = "cond_A",
    configured_seed = NULL,
    clock_seconds = 1735689600L  # 2025-01-01T00:00:00 UTC
  ))
  out2 <- suppressMessages(resolve_lilace_seed(
    condition_name = "cond_A",
    configured_seed = NULL,
    clock_seconds = 1735689600L
  ))
  expect_equal(out1, out2)
})

test_that("auto-generated seeds differ across conditions at same clock", {
  # The condition-name salt is what prevents concurrent per-condition
  # rule jobs from drawing the same seed when Sys.time() rounds to the
  # same integer second. If this assertion fires, the salt is no longer
  # discriminating and parallel chains for distinct conditions would
  # share initial state.
  seed_a <- suppressMessages(resolve_lilace_seed(
    condition_name = "cond_A",
    configured_seed = NULL,
    clock_seconds = 1735689600L
  ))
  seed_b <- suppressMessages(resolve_lilace_seed(
    condition_name = "cond_B",
    configured_seed = NULL,
    clock_seconds = 1735689600L
  ))
  expect_false(seed_a == seed_b)
})

# =============================================================================
# Test 3: Auto-generated seed is a positive integer in Stan's valid range
# =============================================================================

test_that("auto-generated seed is non-negative and within int range", {
  out <- suppressMessages(resolve_lilace_seed(
    condition_name = "cond_A",
    configured_seed = NULL,
    clock_seconds = 1735689600L
  ))
  expect_true(out >= 0)
  expect_true(out < .Machine$integer.max)
})

# =============================================================================
# Test 4: The message includes the seed and a copy-pasteable config snippet
# =============================================================================

test_that("seed log message reports value and a reproducible config hint", {
  expect_message(
    resolve_lilace_seed("cond_A", configured_seed = 12345L),
    "Lilace Stan seed for condition 'cond_A': 12345 \\(config\\).*lilace_seed: 12345"
  )
})

test_that("auto-generated path is labeled in the log message", {
  expect_message(
    resolve_lilace_seed("cond_A", configured_seed = NULL, clock_seconds = 1L),
    "Lilace Stan seed for condition 'cond_A': .* \\(auto-generated\\)"
  )
})
