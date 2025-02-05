#!/usr/bin/env Rscript
library("readr")
library("stringr")
library("dplyr")
library("cmdstanr")
library("rosace")

library("purrr")


## TODO: Position is not being written in output. Defined properly?
## TODO: wildtype is not being written in output. Defined properly?
## TODO: type is not being written in output. Defined properly?
## TODO: Check cmdstanr and toolchain from rosace?

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

stop_logging <- function() {
  # Close the sink connections so logs are properly flushed
  sink(type = "message")
  sink(type = "output")
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

  WT <- substr(hgvs_string, 1, 1)

  # WT case, as in enrich2 format
  if (str_detect(hgvs_string, "_wt")) {
    variant <- "Z"
    pos <- -1
    len <- -1
    mutation_type <- "X"
  }

  # M/S/N
  if (str_detect(hgvs_string, "[A-Z][0-9]+[A-Z]+")) {
    len <- 1
    match <- str_match(hgvs_string, "([A-Z])([0-9]+)([A-Z]+)")
    pos <- match[3]
    if (match[2] == match[4]) {
      mutation_type <- "synonymous"
      variant <- match[4]
    } else if (match[4] == "X") {
      mutation_type <- "nonsense"
      variant <- match[4]
    } else {
      mutation_type <- "missense"
      variant <- match[4]
    }
  }

  # D (deletions)
  if (str_detect(hgvs_string, ".*del")) {
    mutation_type <- "deletion"
    if (str_detect(hgvs_string, "[A-Z][0-9]+_[A-Z][0-9]+del")) {
      # e.g. D_2, D_3
      match <- str_match(hgvs_string, "[A-Z]([0-9]+)_[A-Z]([0-9]+)del")
      pos <- match[2]
      len <- as.integer(match[3]) - as.integer(match[2]) + 1
      variant <- paste0("D_", as.character(len))
    } else {
      # e.g. D_1
      len <- 1
      match <- str_match(hgvs_string, "([A-Z])([0-9]+)del")
      pos <- match[3]
      variant <- "D_1"
    }
  }

  # I (insertions)
  if (str_detect(hgvs_string, ".*ins.*")) {
    mutation_type <- "insertion"
    match <- str_match(hgvs_string, "[A-Z]([0-9]+)_[A-Z][0-9]+ins([A-Z]+)")
    len <- nchar(match[3])
    pos <- match[2]
    variant <- paste0("I_", as.character(len))
  }

  return(c(variant, pos, len, mutation_type, WT))
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

  for (expt_time in unique(df_subset$time)) {
    # Filter by time
    experimental_time_df <- filter(df_subset, time == expt_time)

    # If 'tile' column does not exist, treat as single group
    has_tile <- "tile" %in% colnames(experimental_time_df)
    tile_values <- if (has_tile) {
      unique(experimental_time_df$tile)
    } else {
      NA  # Just one iteration
    }

    for (tile_val in tile_values) {
      if (is.na(tile_val) || !has_tile) {
        # If tile is missing or no tile column, use entire time subset
        experimental_tile_df <- experimental_time_df
      } else {
        # Filter by tile
        experimental_tile_df <- filter(experimental_time_df, tile == tile_val)
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

      col_name <- paste0("c_", expt_time)

      # Join into the main counts tibble
      counts <- full_join(counts, this_count, by = "hgvs") %>%
        rename(!!col_name := count)
    }
  }

  return(counts)
}


# ===========================
#   BUILD ROSACE OBJECT
# ===========================
build_rosace_object <- function(experiment_definition, experiment_name, baseline_condition) {
  # This function:
  #   - Loops over conditions (excluding baseline)
  #   - Loops over replicates
  #   - Builds an assay object from the combined counts
  #   - Adds the assay to a global or local 'rosace' object
  #   - Returns the Rosace object

  conditions <- unique(experiment_definition$condition)
  conditions <- conditions[conditions != baseline_condition]

  rosace_created <- FALSE
  rosace_obj <- NULL  # We'll return this at the end

  for (expt_condition in conditions) {
    # Subset by condition
    condition_df <- filter(experiment_definition, condition == expt_condition)
    for (expt_replicate in unique(condition_df$replicate)) {
      replicate_df <- filter(condition_df, replicate == expt_replicate)

      # Build the counts for this replicate
      counts <- build_counts_for_replicate(replicate_df, experiment_name)
      if (nrow(counts) == 0) {
        # Possibly no data
        warning(sprintf("No valid data for condition=%s replicate=%s", expt_condition, expt_replicate))
        next
      }

      # Build an assay
      # The first column in 'counts' is 'hgvs'; the rest are numeric
      assay <- CreateAssayObject(
        counts = as.matrix(counts[, 2:ncol(counts)]),
        var.names = counts$hgvs,
        key = expt_condition,
        rep = expt_replicate,
        type = "growth"
      )

      # If rosace isn't created, create a new object
      if (!rosace_created) {
        rosace_obj <- CreateRosaceObject(object = assay)
        rosace_created <- TRUE
      } else {
        # else add to the existing object
        rosace_obj <- AddAssayData(object = rosace_obj, assay = assay)
      }

      # Filter, impute, normalize, integrate
      rosace_obj <- FilterData(rosace_obj, key = expt_condition, na.rmax = 0.5)
      rosace_obj <- ImputeData(
        rosace_obj,
        key = expt_condition,
        impute.method = "knn",
        na.rmax = 0.5
      )
      rosace_obj <- NormalizeData(
        rosace_obj,
        key = expt_condition,
        normalization.method = "wt",
        wt.var.names = c("_wt"),
        wt.rm = TRUE
      )
      rosace_obj <- IntegrateData(object = rosace_obj, key = expt_condition)
    }
  }

  if (!rosace_created) {
    warning("No Rosace object was created; possibly no valid data rows were found.")
  }

  return(rosace_obj)
}


# ===========================
#   FINALIZE VARIANT METADATA
# ===========================
finalize_variants_in_rosace <- function(rosace_obj) {
  # Parse the variant strings, storing them in new columns
  # e.g., rosace_obj@var.data
  if (is.null(rosace_obj)) {
    # No data to finalize
    return(NULL)
  }

  # The var.data data frame has a column 'variants' that looks like "p.(...)"?
  # We strip that substring and parse

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

  return(rosace_obj)
}


# ===========================
#   RUN ROSACE ANALYSIS
# ===========================
run_rosace_for_conditions <- function(rosace_obj, experiment_definition, experiment_name, baseline_condition) {
  # Run Rosace for each condition except baseline, writing results
  if (is.null(rosace_obj)) {
    warning("Rosace object is NULL; skipping run_rosace_for_conditions.")
    return(invisible(NULL))
  }

  rosace_dir <- file.path("results", experiment_name, "rosace")
  all_conditions <- unique(experiment_definition$condition)
  conditions <- all_conditions[all_conditions != baseline_condition]

  for (expt_condition in conditions) {
    message(sprintf("Running ROSACE for condition: %s", expt_condition))

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

    scores.data <- OutputScore(
      rosace_obj,
      name = paste0(expt_condition, "_ROSACE")
    )

    output_file_name <- paste0(expt_condition, "_scores.csv")
    output_file <- file.path(rosace_dir, output_file_name)

    utils::write.csv(scores.data, file = output_file, row.names = FALSE)
    message(sprintf("Results written to: %s", output_file))
  }

  # Return the updated rosace object (in case it changed)
  return(rosace_obj)
}


# ===========================
#           MAIN
# ===========================
main <- function() {
  # Activate environment if using renv
  renv::activate()

  # Logging setup
  con <- init_logging()

  # Read configs
  experiment_name <- snakemake@config[["experiment"]]
  experiment_file <- snakemake@config[["experiment_file"]]
  baseline_condition <- snakemake@config[["baseline_condition"]]

  tiled <- snakemake@config[["tiled"]]

  # Create directories
  rosace_dir <- file.path("results", experiment_name, "rosace")
  rosace_assayset_dir <- file.path(rosace_dir, "assayset")
  message(sprintf("Creating directories: %s, %s", rosace_dir, rosace_assayset_dir))
  dir.create(rosace_assayset_dir, recursive = TRUE, showWarnings = FALSE)

  # Load experiment definition
  experiment_definition <- load_experiment_definition(experiment_file)

  # Build the Rosace object (including filter/impute/normalize/integrate)
  rosace_obj <- build_rosace_object(experiment_definition, experiment_name, baseline_condition)

  # Parse the variant strings in var.data
  rosace_obj <- finalize_variants_in_rosace(rosace_obj)

  # Run the Rosace analysis on each condition
  print("Running ROSACE")

  rosace_obj <- run_rosace_for_conditions(
    rosace_obj,
    experiment_definition,
    experiment_name,
    baseline_condition
  )

  # Clean up
  stop_logging()
}

# Run main
main()