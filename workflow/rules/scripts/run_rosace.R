#!/usr/bin/env Rscript

# ===========================
#        INITIALIZATION
# ===========================
init_logging <- function() {
  # Set up logging via snakemake
  log_file <- snakemake@log[[1]]

  # Open connection to log file
  con <- file(log_file, open = "wt")
  sink(con, append = FALSE, split = FALSE, type = "message")
  sink(con, append = FALSE, split = FALSE, type = "output")

  # Optionally return the connection if you want to do something with it later
  return(invisible(con))
}

stop_logging <- function(con) {
  # Close the sink connections so logs are properly flushed
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


# ===========================
#   LOAD EXPERIMENT DEFINITION
# ===========================
load_experiment_definition <- function(experiment_file) {
  # In case you need some validation or transformations in the future
  experiment_definition <- read_csv(experiment_file)
  return(experiment_definition)
}


# ===========================
#   HELPER: PARSE HGVS STRING
# ===========================
parse_stripped_hgvs <- function(hgvs_string) {
  # Returns c(variant, pos, len, mutation_type, WT)
  variant <- ""
  pos <- -1
  len <- -1
  mutation_type <- ""
  WT <- ""

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


# ===========================
#   BUILD COUNTS FOR REPLICATE
# ===========================
build_counts_for_replicate <- function(df_subset, experiment_name) {
  # df_subset should be the subset of experiment_definition
  # filtered by condition, replicate, etc.

  # We'll accumulate data in a tibble called `counts`.
  # If tile is present, we may loop over times & tile. Otherwise just times.

  # Initialize an empty tibble
  counts <- tibble::tibble(hgvs = character())
  message(sprintf("Building counts for condition=%s, replicate=%s", df_subset$condition[1], df_subset$replicate[1]))

  for (expt_time in unique(df_subset$time)) {
    message(sprintf("Processing time=%s", expt_time))
    # Filter by time
    experimental_time_df <- filter(df_subset, time == expt_time)

    # If 'tile' column does not exist, treat as single group
    has_tile <- "tile" %in% colnames(experimental_time_df)
    tile_values <- if (has_tile) {
      unique(experimental_time_df$tile)
    } else {
      NA # Just one iteration
    }

    # Accumulate counts across tiles for this time point, then join once
    time_count <- tibble::tibble(hgvs = character(), count = integer())

    for (tile_val in tile_values) {
      if (is.na(tile_val) || !has_tile) {
        # If tile is missing or no tile column, use entire time subset
        experimental_tile_df <- experimental_time_df
      } else {
        # Filter by tile
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
      this_count <- read_tsv(file_name)

      if (anyDuplicated(this_count$hgvs) > 0) {
        stop(sprintf(
          "Duplicate HGVS entries found in %s. Each variant must appear at most once per sample.",
          file_name
        ))
      }

      # Combine tile counts (tiles cover disjoint variants)
      overlap <- intersect(time_count$hgvs, this_count$hgvs)
      if (length(overlap) > 0) {
        stop(sprintf(
          "Overlapping HGVS variants across tiles at time %s: %s. Tiles must cover disjoint variants.",
          expt_time, paste(head(overlap, 5), collapse = ", ")
        ))
      }
      time_count <- bind_rows(time_count, this_count)
    }

    if (nrow(time_count) > 0) {
      col_name <- paste0("c_", expt_time)
      counts <- full_join(counts, time_count, by = "hgvs") %>%
        rename(!!col_name := count)
    }
  }

  return(counts)
}


# ===========================
#   BUILD ROSACE OBJECT
# ===========================
build_rosace_object_for_condition <- function(experiment_definition, experiment_name, expt_condition, noprocess) {
  # For the given condition:
  #   1. Add every replicate's assay to a fresh Rosace object.
  #   2. Once all replicates are present, run FilterData / ImputeData /
  #      NormalizeData / IntegrateData ONCE against the condition's key.
  #
  # IntegrateData consolidates all replicates sharing a key into a single
  # AssaySet, and NormalizeData (wt-relative) needs visibility across all
  # replicates to compute a meaningful WT baseline — so both must run after
  # the full replicate set for the condition is loaded, not after each
  # individual replicate.

  condition_df <- filter(experiment_definition, condition == expt_condition)

  rosace_created <- FALSE
  rosace_obj <- NULL
  condition_has_data <- FALSE

  # Pass 1: add every valid replicate's assay to the object.
  for (expt_replicate in unique(condition_df$replicate)) {
    replicate_df <- filter(condition_df, replicate == expt_replicate)

    counts <- build_counts_for_replicate(replicate_df, experiment_name)
    if (nrow(counts) == 0) {
      warning(sprintf("No valid data for condition=%s replicate=%s", expt_condition, expt_replicate))
      next
    }

    # The first column in 'counts' is 'hgvs'; the rest are numeric
    assay <- CreateAssayObject(
      counts = as.matrix(counts[, 2:ncol(counts)]),
      var.names = counts$hgvs,
      key = expt_condition,
      rep = expt_replicate,
      type = "growth"
    )

    if (!rosace_created) {
      rosace_obj <- CreateRosaceObject(object = assay)
      rosace_created <- TRUE
    } else {
      rosace_obj <- AddAssayData(object = rosace_obj, assay = assay)
    }
    condition_has_data <- TRUE
  }

  if (!condition_has_data) {
    stop(sprintf(
      "No valid replicates for condition=%s. Cannot produce scores.",
      expt_condition
    ))
  }

  # Pass 2: filter, impute, normalize, integrate ONCE after all replicates
  # are loaded.
  message(sprintf("Filter/impute/normalize/integrate for condition: %s", expt_condition))
  rosace_obj <- FilterData(rosace_obj, key = expt_condition, na.rmax = 0.5)
  rosace_obj <- ImputeData(
    rosace_obj,
    key = expt_condition,
    impute.method = "knn",
    na.rmax = 0.5
  )
  if (noprocess) {
    message(sprintf("Normalizing by total counts since noprocess flag is set."))
    rosace_obj <- NormalizeData(
      rosace_obj,
      key = expt_condition,
      normalization.method = "total"
    )
  } else {
    message(sprintf("Normalizing data by variant names and WT."))
    rosace_obj <- NormalizeData(
      rosace_obj,
      key = expt_condition,
      normalization.method = "wt",
      wt.var.names = c("_wt"),
      wt.rm = TRUE
    )
  }
  rosace_obj <- IntegrateData(object = rosace_obj, key = expt_condition)

  return(rosace_obj)
}


# ===========================
#   FINALIZE VARIANT METADATA
# ===========================
finalize_variants_in_rosace <- function(rosace_obj, noprocess) {
  # Parse the variant strings, storing them in new columns
  # e.g., rosace_obj@var.data
  if (is.null(rosace_obj)) {
    # No data to finalize
    return(NULL)
  }

  # The var.data data frame has a column 'variants' that looks like "p.(...)"?
  # We strip that substring and parse

  if (noprocess) {
    message("Skipping variant parsing since noprocess flag is set.")
    return(rosace_obj)
  }


  rosace_obj@var.data <- rosace_obj@var.data %>%
    mutate(tmp_n = substr(variants, 4, nchar(variants) - 1)) %>%
    # Apply parse_stripped_hgvs() to each row's tmp_n, store as list column
    mutate(parsed = map(tmp_n, parse_stripped_hgvs)) %>%
    # Extract each piece of parsed info into its own column
    mutate(
      position = map_chr(parsed, 2) |> as.numeric(),
      wildtype = map_chr(parsed, 5),
      mutation = map_chr(parsed, 1),
      type     = map_chr(parsed, 4)
    ) %>%
    select(-tmp_n, -parsed)

  message("Parsed variant strings in var.data")
  return(rosace_obj)
}


# ===========================
#   RUN ROSACE ANALYSIS
# ===========================
run_rosace_for_condition <- function(rosace_obj, experiment_name, expt_condition, noprocess) {
  # Run Rosace for one condition and write its scores CSV.
  if (is.null(rosace_obj)) {
    stop("Rosace object is NULL; cannot run RunRosace.")
  }

  rosace_dir <- file.path("results", experiment_name, "rosace")
  message(sprintf("Running ROSACE for condition: %s", expt_condition))

  # If we are using the noprocess flag, variant data will not be parsed.
  # Therefore we cannot use the type or position as controls.
  if (noprocess) {
    rosace_obj <- RunRosace(
      object = rosace_obj,
      name = expt_condition,
      type = "AssaySet",
      savedir = rosace_dir,
      install = FALSE
    )
  } else {
    rosace_obj <- RunRosace(
      object = rosace_obj,
      name = expt_condition,
      type = "AssaySet",
      savedir = rosace_dir,
      pos.col = "position",
      ctrl.col = "type",
      ctrl.name = "synonymous",
      install = FALSE
    )
  }

  scores.data <- OutputScore(
    rosace_obj,
    name = paste0(expt_condition, "_ROSACE")
  )

  output_file_name <- paste0(expt_condition, "_scores.csv")
  output_file <- file.path(rosace_dir, output_file_name)

  utils::write.csv(scores.data, file = output_file, row.names = FALSE)
  message(sprintf("Results written to: %s", output_file))

  return(rosace_obj)
}


# ===========================
#           MAIN
# ===========================
main <- function() {
  # Logging setup

  # Read configs + wildcards
  experiment_name <- snakemake@config[["experiment"]]
  experiment_file <- snakemake@config[["experiment_file"]]
  expt_condition <- snakemake@wildcards[["condition"]]

  tiled <- snakemake@config[["tiled"]]
  noprocess <- snakemake@config[["noprocess"]]

  # Create directories
  rosace_dir <- file.path("results", experiment_name, "rosace")
  rosace_assayset_dir <- file.path(rosace_dir, "assayset")
  message(sprintf("Creating directories: %s, %s", rosace_dir, rosace_assayset_dir))
  dir.create(rosace_assayset_dir, recursive = TRUE, showWarnings = FALSE)

  # Load experiment definition
  experiment_definition <- load_experiment_definition(experiment_file)

  # Build the Rosace object for this one condition (including
  # filter/impute/normalize/integrate scoped to it).
  rosace_obj <- build_rosace_object_for_condition(
    experiment_definition, experiment_name, expt_condition, noprocess
  )

  # Parse the variant strings in var.data
  rosace_obj <- finalize_variants_in_rosace(rosace_obj, noprocess)

  # Run the Rosace analysis on this condition
  print(sprintf("Running ROSACE for condition: %s", expt_condition))

  rosace_obj <- run_rosace_for_condition(
    rosace_obj,
    experiment_name,
    expt_condition,
    noprocess
  )
}

con <- init_logging()
on.exit(
  {
    stop_logging(con)
  },
  add = TRUE
)

# Enable renv if not using user-managed rosace.
# Point .libPaths() at the project library directly instead of going
# through renv::activate(): activate() does a self-backup rename of the
# renv bootstrap dir on every call, so parallel invocations of this rule
# race on that rename. renv::paths$library() is a read-only path
# resolver, so this is concurrency-safe. install_rosace.R is the one
# place that initializes and populates the renv project library (via
# renv::restore()); subsequent runs only need to find it.

if (!snakemake@config[["rosace_local"]]) {
  message(sprintf("Pointing .libPaths() at renv project library"))
  library("renv")
  .libPaths(c(renv::paths$library(), .libPaths()))
}
library("readr")
library("stringr")
library("dplyr")
library("purrr")
library("cmdstanr")
library("rosace")


# Run main
tryCatch(
  main(),
  error = function(err) {
    log_and_rethrow("run_rosace", err)
  }
)
