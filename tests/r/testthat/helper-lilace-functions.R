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


#' Resolve the Stan seed for a Lilace fit (parallel copy of the function
#' in workflow/rules/scripts/run_lilace.R).
#'
#' The production function reads snakemake@config[["lilace_seed"]] at the
#' call site and forwards it into this resolver as `configured_seed`. Keep
#' this helper in lockstep with the source — any logic change in run_lilace.R
#' must be mirrored here, and vice versa.
#'
#' @param condition_name The Snakemake wildcard `condition`, used as salt
#'   when auto-generating so concurrent per-condition rule jobs don't all
#'   pull the same Sys.time()-derived value.
#' @param configured_seed Either NULL (meaning auto-generate) or an integer
#'   the user supplied via config["lilace_seed"]. Passed through as-is when
#'   non-null so repeat runs are bit-identical.
#' @param clock_seconds Injected clock value used only on the auto-generate
#'   path; production callers leave NULL and we read Sys.time().
#' @return Positive integer seed, suitable to pass to lilace_fit_model(seed=).
resolve_lilace_seed <- function(condition_name, configured_seed = NULL,
                                clock_seconds = NULL) {
  if (!is.null(configured_seed)) {
    seed <- as.integer(configured_seed)
    source <- "config"
  } else {
    if (is.null(clock_seconds)) clock_seconds <- as.integer(Sys.time())
    cond_salt <- sum(utf8ToInt(condition_name))
    # Compute in double space — as.integer(Sys.time()) is around 1.7e9 and
    # multiplying by 17 overflows R's 32-bit integer (returns NA). Modulo
    # at the end before converting back to int keeps the result in range.
    seed <- as.integer(
      abs(as.numeric(clock_seconds) * 17 + cond_salt) %%
        .Machine$integer.max
    )
    source <- "auto-generated"
  }
  message(sprintf(
    "Lilace Stan seed for condition '%s': %d (%s). Re-run with `lilace_seed: %d` in config for bit-identical scores.",
    condition_name, seed, source, seed
  ))
  seed
}
