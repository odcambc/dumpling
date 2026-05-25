# test-parse_hgvs_lilace.R
# Unit tests for the Lilace-flavored HGVS parser
# (parse_stripped_hgvs_lilace, defined in helper-lilace-functions.R).

library(testthat)

context("HGVS Parsing (Lilace)")

# -----------------------------------------------------------------------------
# Happy-path: one test per supported HGVS shape.
# -----------------------------------------------------------------------------

test_that("parse_stripped_hgvs_lilace handles wildtype (_wt)", {
  result <- parse_stripped_hgvs_lilace("_wt")
  expect_equal(result$position, 0)
  expect_equal(result$wildtype, "")
  expect_equal(result$mutation, "WT")
  expect_equal(result$type, "synonymous")
})

test_that("parse_stripped_hgvs_lilace handles HGVS-style synonymous (X1=)", {
  result <- parse_stripped_hgvs_lilace("A123=")
  expect_equal(result$position, 123)
  expect_equal(result$wildtype, "A")
  expect_equal(result$mutation, "A")
  expect_equal(result$type, "synonymous")
})

test_that("parse_stripped_hgvs_lilace handles GATK/dimple-style synonymous (X1X)", {
  result <- parse_stripped_hgvs_lilace("A123A")
  expect_equal(result$position, 123)
  expect_equal(result$wildtype, "A")
  expect_equal(result$mutation, "A")
  expect_equal(result$type, "synonymous")
})

test_that("parse_stripped_hgvs_lilace handles missense (X1Y)", {
  result <- parse_stripped_hgvs_lilace("A123T")
  expect_equal(result$position, 123)
  expect_equal(result$wildtype, "A")
  expect_equal(result$mutation, "T")
  expect_equal(result$type, "missense")
})

test_that("parse_stripped_hgvs_lilace handles nonsense (X1*)", {
  result <- parse_stripped_hgvs_lilace("W42*")
  expect_equal(result$position, 42)
  expect_equal(result$wildtype, "W")
  expect_equal(result$mutation, "*")
  expect_equal(result$type, "nonsense")
})

test_that("parse_stripped_hgvs_lilace handles single-codon deletion (X1del)", {
  result <- parse_stripped_hgvs_lilace("A10del")
  expect_equal(result$position, 10)
  expect_equal(result$wildtype, "A")
  expect_equal(result$mutation, "D_1")
  expect_equal(result$type, "deletion")
})

test_that("parse_stripped_hgvs_lilace handles multi-codon deletion (X1_Xn del)", {
  result <- parse_stripped_hgvs_lilace("A10_C12del")
  expect_equal(result$position, 10)
  expect_equal(result$wildtype, "A")
  expect_equal(result$mutation, "D_3")
  expect_equal(result$type, "deletion")
})

test_that("parse_stripped_hgvs_lilace handles single-AA insertion", {
  result <- parse_stripped_hgvs_lilace("A10_C11insG")
  expect_equal(result$position, 10)
  expect_equal(result$wildtype, "A")
  expect_equal(result$mutation, "I_1")
  expect_equal(result$type, "insertion")
})

test_that("parse_stripped_hgvs_lilace handles multi-AA insertion", {
  result <- parse_stripped_hgvs_lilace("A10_C11insGGS")
  expect_equal(result$position, 10)
  expect_equal(result$wildtype, "A")
  expect_equal(result$mutation, "I_3")
  expect_equal(result$type, "insertion")
})

# -----------------------------------------------------------------------------
# Position / amino acid coverage: each supported AA is parsed correctly,
# and large position numbers don't truncate.
# -----------------------------------------------------------------------------

test_that("parse_stripped_hgvs_lilace extracts WT for all 20 standard amino acids", {
  amino_acids <- c(
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
  )
  for (aa in amino_acids) {
    result <- parse_stripped_hgvs_lilace(paste0(aa, "1="))
    expect_equal(result$wildtype, aa, info = paste("Failed for amino acid:", aa))
  }
})

test_that("parse_stripped_hgvs_lilace preserves large position numbers", {
  result <- parse_stripped_hgvs_lilace("A9999T")
  expect_equal(result$position, 9999)
  expect_equal(result$type, "missense")
})

# -----------------------------------------------------------------------------
# Negative path: unsupported formats should error loudly with the input
# echoed back, not return a silently-wrong record.
# -----------------------------------------------------------------------------

test_that("parse_stripped_hgvs_lilace errors on unsupported format", {
  expect_error(
    parse_stripped_hgvs_lilace("not-a-variant"),
    "Unsupported HGVS format"
  )
})

test_that("parse_stripped_hgvs_lilace error message includes the offending input", {
  expect_error(
    parse_stripped_hgvs_lilace("weird-input-xyz"),
    "weird-input-xyz"
  )
})

test_that("parse_stripped_hgvs_lilace rejects empty string", {
  expect_error(parse_stripped_hgvs_lilace(""))
})
