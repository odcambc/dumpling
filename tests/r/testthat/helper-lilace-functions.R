# helper-lilace-functions.R
# Extracted pure function from workflow/rules/scripts/run_lilace.R for testing.
#
# Mirrors the pattern in helper-functions.R (Rosace): the source script
# contains snakemake-only top-level code, so we maintain a parallel testable
# copy here. Any change to parse_stripped_hgvs in run_lilace.R must be
# reflected here, and vice versa.
#
# The Lilace parser differs from the Rosace one in return shape and
# semantics: it returns a list with $position/$wildtype/$mutation/$type
# (vs Rosace's named char vector), and it errors on unsupported formats
# rather than returning sentinel empty values.

suppressPackageStartupMessages({
  library(stringr)
})

#' Parse a stripped HGVS variant string (Lilace flavor)
#'
#' @param hgvs_string HGVS protein string without the `p.(...)` wrapper,
#'   or the literal `_wt` for wildtype rows.
#' @return Named list with fields: position (numeric), wildtype (char),
#'   mutation (char), type (one of "synonymous", "missense", "nonsense",
#'   "deletion", "insertion").
parse_stripped_hgvs_lilace <- function(hgvs_string) {
  if (hgvs_string == "_wt") {
    return(list(position = 0, wildtype = "", mutation = "WT", type = "synonymous"))
  }

  if (stringr::str_detect(hgvs_string, "^[A-Z][0-9]+=$")) {
    match <- stringr::str_match(hgvs_string, "([A-Z])([0-9]+)=")
    return(list(
      position = as.numeric(match[3]),
      wildtype = match[2],
      mutation = match[2],
      type = "synonymous"
    ))
  }

  if (stringr::str_detect(hgvs_string, "^[A-Z][0-9]+[A-Z]$")) {
    match <- stringr::str_match(hgvs_string, "([A-Z])([0-9]+)([A-Z])")
    mut_type <- ifelse(match[2] == match[4], "synonymous", "missense")
    return(list(
      position = as.numeric(match[3]),
      wildtype = match[2],
      mutation = match[4],
      type = mut_type
    ))
  }

  if (stringr::str_detect(hgvs_string, "^[A-Z][0-9]+\\*$")) {
    match <- stringr::str_match(hgvs_string, "([A-Z])([0-9]+)\\*")
    return(list(
      position = as.numeric(match[3]),
      wildtype = match[2],
      mutation = "*",
      type = "nonsense"
    ))
  }

  if (stringr::str_detect(hgvs_string, "^[A-Z][0-9]+_[A-Z][0-9]+del$")) {
    match <- stringr::str_match(hgvs_string, "([A-Z])([0-9]+)_([A-Z])([0-9]+)del")
    len <- as.integer(match[5]) - as.integer(match[3]) + 1
    return(list(
      position = as.numeric(match[3]),
      wildtype = match[2],
      mutation = paste0("D_", len),
      type = "deletion"
    ))
  }

  if (stringr::str_detect(hgvs_string, "^[A-Z][0-9]+del$")) {
    match <- stringr::str_match(hgvs_string, "([A-Z])([0-9]+)del")
    return(list(
      position = as.numeric(match[3]),
      wildtype = match[2],
      mutation = "D_1",
      type = "deletion"
    ))
  }

  if (stringr::str_detect(hgvs_string, "^[A-Z][0-9]+_[A-Z][0-9]+ins[A-Z]+$")) {
    match <- stringr::str_match(hgvs_string, "([A-Z])([0-9]+)_([A-Z])([0-9]+)ins([A-Z]+)")
    return(list(
      position = as.numeric(match[3]),
      wildtype = match[2],
      mutation = paste0("I_", nchar(match[6])),
      type = "insertion"
    ))
  }

  stop(sprintf("Unsupported HGVS format encountered: %s", hgvs_string))
}
