# test-parse_hgvs.R
# Unit tests for HGVS parsing functions

library(testthat)

context("HGVS Parsing")

# =============================================================================
# Tests for parse_stripped_hgvs()
# =============================================================================

test_that("parse_stripped_hgvs handles synonymous variants (X1=)", {
  result <- parse_stripped_hgvs("A1=")

  expect_equal(result["variant"], c(variant = "S"))
  expect_equal(result["pos"], c(pos = "1"))
  expect_equal(result["mutation_type"], c(mutation_type = "synonymous"))
  expect_equal(result["WT"], c(WT = "A"))
})

test_that("parse_stripped_hgvs handles synonymous at different positions", {
  result <- parse_stripped_hgvs("M154=")

  expect_equal(result["variant"], c(variant = "S"))
  expect_equal(result["pos"], c(pos = "154"))
  expect_equal(result["WT"], c(WT = "M"))
})

test_that("parse_stripped_hgvs handles missense variants (X1Y)", {
  result <- parse_stripped_hgvs("A1G")

  expect_equal(result["variant"], c(variant = "G"))
  expect_equal(result["pos"], c(pos = "1"))
  expect_equal(result["mutation_type"], c(mutation_type = "missense"))
  expect_equal(result["WT"], c(WT = "A"))
})

test_that("parse_stripped_hgvs handles missense at different positions", {
  result <- parse_stripped_hgvs("M154L")

  expect_equal(result["variant"], c(variant = "L"))
  expect_equal(result["pos"], c(pos = "154"))
  expect_equal(result["WT"], c(WT = "M"))
})

test_that("parse_stripped_hgvs handles nonsense variants (X1*)", {
  result <- parse_stripped_hgvs("W7*")

  expect_equal(result["variant"], c(variant = "X"))
  expect_equal(result["pos"], c(pos = "7"))
  expect_equal(result["mutation_type"], c(mutation_type = "nonsense"))
  expect_equal(result["WT"], c(WT = "W"))
})

test_that("parse_stripped_hgvs handles single deletions (X1del)", {
  result <- parse_stripped_hgvs("A10del")

  expect_equal(result["variant"], c(variant = "D_1"))
  expect_equal(result["pos"], c(pos = "10"))
  expect_equal(result["len"], c(len = "1"))
  expect_equal(result["mutation_type"], c(mutation_type = "deletion"))
})

test_that("parse_stripped_hgvs handles multi-codon deletions (X1_X3del)", {
  result <- parse_stripped_hgvs("A1_A3del")

  expect_equal(result["variant"], c(variant = "D_3"))
  expect_equal(result["pos"], c(pos = "1"))
  expect_equal(result["len"], c(len = "3"))
  expect_equal(result["mutation_type"], c(mutation_type = "deletion"))
})

test_that("parse_stripped_hgvs handles 2-codon deletions", {
  result <- parse_stripped_hgvs("M10_K11del")

  expect_equal(result["variant"], c(variant = "D_2"))
  expect_equal(result["pos"], c(pos = "10"))
  expect_equal(result["len"], c(len = "2"))
})

test_that("parse_stripped_hgvs handles insertions (X1_X2insAAA)", {
  result <- parse_stripped_hgvs("A1_A2insG")

  expect_equal(result["variant"], c(variant = "I_1"))
  expect_equal(result["pos"], c(pos = "1"))
  expect_equal(result["len"], c(len = "1"))
  expect_equal(result["mutation_type"], c(mutation_type = "insertion"))
})

test_that("parse_stripped_hgvs handles multi-AA insertions", {
  result <- parse_stripped_hgvs("A1_A2insGGG")

  expect_equal(result["variant"], c(variant = "I_3"))
  expect_equal(result["pos"], c(pos = "1"))
  expect_equal(result["len"], c(len = "3"))
})

test_that("parse_stripped_hgvs handles WT case", {
  result <- parse_stripped_hgvs("_wt")

  expect_equal(result["variant"], c(variant = "WT"))
  expect_equal(result["pos"], c(pos = "0"))
  expect_equal(result["mutation_type"], c(mutation_type = "WT"))
})

# =============================================================================
# Edge cases
# =============================================================================

test_that("parse_stripped_hgvs extracts WT amino acid correctly", {
  # Test various amino acids
  amino_acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
                   "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")

  for (aa in amino_acids) {
    result <- parse_stripped_hgvs(paste0(aa, "1="))
    expect_equal(result["WT"], c(WT = aa),
                 info = paste("Failed for amino acid:", aa))
  }
})

test_that("parse_stripped_hgvs handles large position numbers", {
  result <- parse_stripped_hgvs("A9999=")
  expect_equal(result["pos"], c(pos = "9999"))

  result2 <- parse_stripped_hgvs("M1234G")
  expect_equal(result2["pos"], c(pos = "1234"))
})
