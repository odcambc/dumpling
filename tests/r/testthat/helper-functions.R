# helper-functions.R
# Extracted pure functions from run_rosace.R for testing
# These functions don't depend on snakemake and can be tested in isolation

library(stringr)
library(dplyr)
library(readr)

#' Parse HGVS string (stripped of p.() wrapper)
#'
#' @param hgvs_string HGVS string without the p.() wrapper
#' @return Named vector with: variant, pos, len, mutation_type, WT
#'
#' @examples
#' parse_stripped_hgvs("A1=") # Synonymous
#' parse_stripped_hgvs("A1B") # Missense
#' parse_stripped_hgvs("A1*") # Nonsense
#' parse_stripped_hgvs("A1del") # Single deletion
#' parse_stripped_hgvs("A1_A3del") # Multi deletion
#' parse_stripped_hgvs("A1_A2insGGG") # Insertion
parse_stripped_hgvs <- function(hgvs_string) {
  variant <- ""
  pos <- -1
  len <- -1
  mutation_type <- ""
  WT <- ""

  WT <- substr(hgvs_string, 1, 1)

  # WT case, as in enrich2 format
  if (hgvs_string == "_wt") {
    variant <- "WT"
    pos <- 0
    len <- 0
    mutation_type <- "X"
    return(c(
      variant = variant, pos = as.character(pos), len = as.character(len),
      mutation_type = mutation_type, WT = WT
    ))
  }

  # S (synonymous) - format: X1=
  if (str_detect(hgvs_string, ".*=$")) {
    mutation_type <- "synonymous"
    match <- str_match(hgvs_string, "([A-Z])([0-9]+)=")
    if (!is.na(match[1])) {
      WT <- match[2]
      pos <- match[3]
      len <- 1
      variant <- "S"
    }
  }

  # M (missense) or S (synonymous via GATK/dimple format: X1X where letters match)
  if (str_detect(hgvs_string, "^[A-Z][0-9]+[A-Z]$")) {
    match <- str_match(hgvs_string, "([A-Z])([0-9]+)([A-Z])")
    if (!is.na(match[1])) {
      WT <- match[2]
      pos <- match[3]
      len <- 1
      if (match[2] == match[4]) {
        # Same amino acid = synonymous (GATK/dimple notation, e.g. A10A)
        mutation_type <- "synonymous"
        variant <- "S"
      } else {
        mutation_type <- "missense"
        variant <- match[4]
      }
    }
  }

  # N (nonsense) - format: X1*
  if (str_detect(hgvs_string, ".*\\*$")) {
    mutation_type <- "nonsense"
    match <- str_match(hgvs_string, "([A-Z])([0-9]+)\\*")
    if (!is.na(match[1])) {
      WT <- match[2]
      pos <- match[3]
      len <- 1
      variant <- "X"
    }
  }

  # D (deletions)
  if (str_detect(hgvs_string, ".*del$")) {
    mutation_type <- "deletion"
    if (str_detect(hgvs_string, "[A-Z][0-9]+_[A-Z][0-9]+del")) {
      # Multi-codon deletion: e.g., A1_A3del
      match <- str_match(hgvs_string, "([A-Z])([0-9]+)_[A-Z]([0-9]+)del")
      if (!is.na(match[1])) {
        WT <- match[2]
        pos <- match[3]
        len <- as.integer(match[4]) - as.integer(match[3]) + 1
        variant <- paste0("D_", as.character(len))
      }
    } else {
      # Single-codon deletion: e.g., A1del
      len <- 1
      match <- str_match(hgvs_string, "([A-Z])([0-9]+)del")
      if (!is.na(match[1])) {
        WT <- match[2]
        pos <- match[3]
        variant <- "D_1"
      }
    }
  }

  # I (insertions)
  if (str_detect(hgvs_string, ".*ins[A-Z]+$")) {
    mutation_type <- "insertion"
    match <- str_match(hgvs_string, "([A-Z])([0-9]+)_[A-Z]([0-9]+)ins([A-Z]+)")
    if (!is.na(match[1])) {
      WT <- match[2]
      pos <- match[3]
      inserted <- match[5]
      len <- nchar(inserted)
      variant <- paste0("I_", as.character(len))
    }
  }

  return(c(
    variant = variant, pos = as.character(pos), len = as.character(len),
    mutation_type = mutation_type, WT = WT
  ))
}


#' Load experiment definition from CSV
#'
#' @param experiment_file Path to experiment CSV file
#' @return Data frame with experiment metadata
load_experiment_definition <- function(experiment_file) {
  experiment_definition <- readr::read_csv(experiment_file, show_col_types = FALSE)
  return(experiment_definition)
}


#' Get unique conditions from experiment definition
#'
#' @param experiment_df Experiment data frame
#' @param baseline_condition Name of baseline condition to exclude
#' @return Vector of condition names
get_conditions <- function(experiment_df, baseline_condition = "baseline") {
  conditions <- unique(experiment_df$condition)
  conditions <- conditions[conditions != baseline_condition]
  return(conditions)
}


#' Get samples for a specific condition
#'
#' @param experiment_df Experiment data frame
#' @param condition Condition name
#' @return Vector of sample names
get_samples_for_condition <- function(experiment_df, condition) {
  samples <- experiment_df$sample[experiment_df$condition == condition]
  return(samples)
}


#' Validate experiment structure
#'
#' @param experiment_df Experiment data frame
#' @return TRUE if valid, error message otherwise
validate_experiment <- function(experiment_df) {
  required_cols <- c("sample", "condition", "replicate", "time")

  missing_cols <- setdiff(required_cols, names(experiment_df))
  if (length(missing_cols) > 0) {
    return(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  # Check for T0 in each condition
  for (cond in unique(experiment_df$condition)) {
    times <- experiment_df$time[experiment_df$condition == cond]
    if (!0 %in% times) {
      return(paste("Condition", cond, "missing T0 timepoint"))
    }
  }

  return(TRUE)
}


#' Build count matrix for a single replicate
#'
#' Reads per-sample TSV files (hgvs + count columns) from the pipeline's
#' enrich_format output directory and joins them into a wide tibble with one
#' column per time point (c_0, c_1, ...). Rows are the union of all HGVS
#' variants observed across time points (NA where absent).
#'
#' @param df_subset Rows of the experiment definition for one condition/replicate
#' @param experiment_name Experiment name (used to build the file path)
#' @return Tibble with columns: hgvs, c_<time>, ...
build_counts_for_replicate <- function(df_subset, experiment_name) {
  counts <- tibble::tibble(hgvs = character())
  message(sprintf("Building counts for condition=%s, replicate=%s", df_subset$condition[1], df_subset$replicate[1]))

  for (expt_time in unique(df_subset$time)) {
    message(sprintf("Processing time=%s", expt_time))
    experimental_time_df <- filter(df_subset, time == expt_time)

    has_tile <- "tile" %in% colnames(experimental_time_df)
    tile_values <- if (has_tile) {
      unique(experimental_time_df$tile)
    } else {
      NA
    }

    for (tile_val in tile_values) {
      if (is.na(tile_val) || !has_tile) {
        experimental_tile_df <- experimental_time_df
      } else {
        experimental_tile_df <- filter(experimental_time_df, tile == tile_val)
        message(sprintf("Processing tile=%s", tile_val))
      }

      if (nrow(experimental_tile_df) > 1) {
        stop(sprintf(
          paste(
            "More than one row in the experiment definition for condition %s,",
            "replicate %s, tile %s, time %s."
          ),
          experimental_tile_df$condition[1],
          experimental_tile_df$replicate[1],
          tile_val,
          expt_time
        ))
      }
      if (nrow(experimental_tile_df) == 0) {
        warning(sprintf(
          paste(
            "No rows in the experiment definition for condition %s,",
            "replicate %s, tile %s, time %s."
          ),
          df_subset$condition[1],
          df_subset$replicate[1],
          tile_val,
          expt_time
        ))
        next
      }

      sample_name <- experimental_tile_df$sample[1]
      file_name <- paste0(
        "results/", experiment_name, "/processed_counts/enrich_format/",
        sample_name, ".tsv"
      )
      this_count <- read_tsv(file_name, show_col_types = FALSE)

      if (anyDuplicated(this_count$hgvs) > 0) {
        stop(sprintf(
          "Duplicate HGVS entries found in %s. Each variant must appear at most once per sample.",
          file_name
        ))
      }

      col_name <- paste0("c_", expt_time)

      counts <- full_join(counts, this_count, by = "hgvs") %>%
        rename(!!col_name := count)
    }
  }

  return(counts)
}
