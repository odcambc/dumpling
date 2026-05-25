#!/usr/bin/env Rscript

init_logging <- function() {
  log_file <- snakemake@log[[1]]
  con <- file(log_file, open = "wt")
  sink(con, append = FALSE, split = FALSE, type = "message")
  sink(con, append = FALSE, split = FALSE, type = "output")
  invisible(con)
}

stop_logging <- function(con) {
  sink(type = "message")
  sink(type = "output")
  close(con)
}

log_and_rethrow <- function(script_name, err) {
  message(sprintf("%s failed: %s", script_name, conditionMessage(err)))
  calls <- vapply(
    sys.calls(),
    function(call) paste(deparse(call), collapse = " "),
    character(1)
  )
  if (length(calls) > 0) {
    message("Call stack:")
    message(paste(calls, collapse = "\n"))
  }
  stop(err)
}

load_experiment_definition <- function(experiment_file) {
  readr::read_csv(experiment_file, show_col_types = FALSE)
}

parse_stripped_hgvs <- function(hgvs_string) {
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

annotate_variants <- function(counts) {
  stripped <- counts$hgvs
  stripped <- ifelse(stripped == "_wt", stripped, substr(stripped, 4, nchar(stripped) - 1))
  parsed <- purrr::map(stripped, parse_stripped_hgvs)

  counts %>%
    dplyr::mutate(
      variant = .data$hgvs,
      position = purrr::map_dbl(parsed, "position"),
      wildtype = purrr::map_chr(parsed, "wildtype"),
      mutation = purrr::map_chr(parsed, "mutation"),
      type = purrr::map_chr(parsed, "type")
    ) %>%
    dplyr::select(.data$variant, .data$type, .data$wildtype, .data$mutation, .data$position, dplyr::everything())
}

build_counts_for_replicate <- function(df_subset, experiment_name) {
  counts <- tibble::tibble(hgvs = character())

  for (expt_time in unique(df_subset$time)) {
    experimental_time_df <- dplyr::filter(df_subset, .data$time == expt_time)
    has_tile <- "tile" %in% colnames(experimental_time_df)
    tile_values <- if (has_tile) unique(experimental_time_df$tile) else NA

    time_count <- tibble::tibble(hgvs = character(), count = integer())

    for (tile_val in tile_values) {
      experimental_tile_df <- if (is.na(tile_val) || !has_tile) {
        experimental_time_df
      } else {
        dplyr::filter(experimental_time_df, .data$tile == tile_val)
      }

      if (nrow(experimental_tile_df) > 1) {
        stop(sprintf(
          "More than one row in experiment definition for condition=%s replicate=%s tile=%s time=%s",
          experimental_tile_df$condition[1],
          experimental_tile_df$replicate[1],
          tile_val,
          expt_time
        ))
      }
      if (nrow(experimental_tile_df) == 0) {
        next
      }

      sample_name <- experimental_tile_df$sample[1]
      file_name <- file.path(
        "results", experiment_name, "processed_counts", "enrich_format", paste0(sample_name, ".tsv")
      )
      if (!file.exists(file_name)) {
        stop(sprintf(
          "Expected per-sample count file is missing: %s. Did upstream GATK / process_counts fail for sample '%s'?",
          file_name, sample_name
        ))
      }
      this_count <- readr::read_tsv(file_name, show_col_types = FALSE)

      if (anyDuplicated(this_count$hgvs) > 0) {
        stop(sprintf("Duplicate HGVS entries found in %s", file_name))
      }

      overlap <- intersect(time_count$hgvs, this_count$hgvs)
      if (length(overlap) > 0) {
        stop(sprintf(
          "Overlapping HGVS variants across tiles at time %s: %s",
          expt_time,
          paste(head(overlap, 5), collapse = ", ")
        ))
      }
      time_count <- dplyr::bind_rows(time_count, this_count)
    }

    if (nrow(time_count) > 0) {
      col_name <- paste0("c_", expt_time)
      counts <- dplyr::full_join(counts, time_count, by = "hgvs") %>%
        dplyr::rename(!!col_name := .data$count)
    }
  }

  counts
}

run_lilace_for_condition <- function(experiment_definition, experiment_name, condition_name) {
  condition_df <- dplyr::filter(experiment_definition, .data$condition == condition_name)
  replicate_frames <- list()

  for (expt_replicate in unique(condition_df$replicate)) {
    replicate_df <- dplyr::filter(condition_df, .data$replicate == expt_replicate)
    counts <- build_counts_for_replicate(replicate_df, experiment_name)
    if (nrow(counts) == 0) {
      next
    }

    counts <- annotate_variants(counts)
    c_cols <- grep("^c_", colnames(counts), value = TRUE)
    counts[c_cols] <- lapply(counts[c_cols], function(x) ifelse(is.na(x), 0L, as.integer(x)))

    counts$rep <- as.character(expt_replicate)
    replicate_frames[[as.character(expt_replicate)]] <- counts
  }

  if (length(replicate_frames) == 0) {
    warning(sprintf("No valid counts found for condition %s", condition_name))
    return(invisible(NULL))
  }

  model_df <- dplyr::bind_rows(replicate_frames) %>%
    dplyr::filter(.data$variant != "_wt")

  c_cols <- grep("^c_", colnames(model_df), value = TRUE)
  c_cols <- c_cols[order(gsub("^c_", "", c_cols))]

  lilace_obj <- lilace::lilace_from_counts(
    variant_id = model_df$variant,
    mutation_type = model_df$type,
    position = model_df$position,
    replicate = model_df$rep,
    counts = as.matrix(model_df[, c_cols]),
    metadata = model_df[, c("wildtype", "mutation")]
  )

  output_dir <- file.path("results", experiment_name, "lilace", condition_name)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  lilace_obj <- lilace::lilace_fit_model(
    lilace_obj,
    output_dir = output_dir,
    control_label = "synonymous",
    pseudocount = TRUE,
    n_parallel_chains = snakemake@threads
  )

  output_scores <- file.path("results", experiment_name, "lilace", paste0(condition_name, "_scores.csv"))
  utils::write.csv(lilace_obj$scores, file = output_scores, row.names = FALSE)
}

main <- function() {
  experiment_name <- snakemake@config[["experiment"]]
  experiment_file <- snakemake@config[["experiment_file"]]
  baseline_condition <- snakemake@config[["baseline_condition"]]
  noprocess <- snakemake@config[["noprocess"]]

  # Backstop: common.smk's validate_scoring_backend_mode rejects this
  # combination at parse time, but guard here too in case the script is
  # invoked outside the Snakefile (manual reruns, direct script execution).
  # Lilace's Stan model requires a synonymous control set and has no
  # total-counts fallback (see pimentellab/lilace R/input.R:39-45).
  if (isTRUE(noprocess)) {
    stop(paste(
      "scoring_backend='lilace' is incompatible with noprocess=true.",
      "Lilace requires a synonymous-variant control set and has no",
      "total-counts fallback. Use scoring_backend='rosace' for noprocess",
      "runs, or set noprocess=false."
    ))
  }

  lilace_dir <- file.path("results", experiment_name, "lilace")
  dir.create(lilace_dir, recursive = TRUE, showWarnings = FALSE)

  experiment_definition <- load_experiment_definition(experiment_file)
  conditions <- unique(experiment_definition$condition)
  conditions <- conditions[conditions != baseline_condition]

  for (condition_name in conditions) {
    message(sprintf("Running Lilace for condition: %s", condition_name))
    run_lilace_for_condition(experiment_definition, experiment_name, condition_name)
  }
}

con <- init_logging()
on.exit({
  stop_logging(con)
}, add = TRUE)

if (!snakemake@config[["lilace_local"]]) {
  message("Activating renv environment")
  renv::activate()
}

library("readr")
library("stringr")
library("dplyr")
library("purrr")
library("cmdstanr")
library("lilace")

tryCatch(
  main(),
  error = function(err) {
    log_and_rethrow("run_lilace", err)
  }
)
