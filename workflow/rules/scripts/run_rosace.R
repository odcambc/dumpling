#!/usr/bin/env Rscript

# Set up logging
# Output will be redicted to the log file
log_file <- snakemake@log[[1]]

# Open connection to log file
con <- file(log_file, open = "wt")

sink(con, append = FALSE, split = FALSE, type = "message")
sink(con, append = FALSE, split = FALSE, type = "output")


renv::activate()

library("readr")
library("stringr")
library("dplyr")
library("cmdstanr")
library("rosace")

# This metadata is provided by the snakemake invocation.

experiment_name <- snakemake@config[["experiment"]]
experiment_file <- snakemake@config[["experiment_file"]]
baseline_condition <- snakemake@config[["baseline_condition"]]

# Parses an HGVS string minus surrounding "p()",
# outputs a vector with variant info

parse_stripped_hgvs <- function(hgvs_string) {
  variant <- ""
  pos <- -1
  len <- -1
  mutation_type <- ""
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

  # D
  if (str_detect(hgvs_string, ".*del")) {
    mutation_type <- "deletion"
    if (str_detect(hgvs_string, "[A-Z][0-9]+_[A-Z][0-9]+del")) {
      # D_2, D_3
      match <- str_match(hgvs_string, "[A-Z]([0-9]+)_[A-Z]([0-9]+)del")
      pos <- match[2]
      len <- strtoi(match[3]) - strtoi(match[2]) + 1
      variant <- paste("D_", as.character(len), sep = "")
    } else {
      # D_1
      len <- 1
      match <- str_match(hgvs_string, "([A-Z])([0-9]+)del")
      pos <- match[3]
      variant <- "D_1"
    }
  }

  # I
  if (str_detect(hgvs_string, ".*ins.*")) {
    mutation_type <- "insertion"
    match <- str_match(hgvs_string, "[A-Z]([0-9]+)_[A-Z][0-9]+ins([A-Z]+)")
    len <- nchar(match[3])
    pos <- match[2]
    variant <- paste("I_", as.character(len), sep = "")
  }
  return(c(variant, pos, len, mutation_type, WT))
}

# Read in the experiment definition
experiment_definition <- read_csv(experiment_file)

# Create the results directories, if they don't exist

rosace_dir <- file.path(
  "./results",
  experiment_name,
  "rosace",
  sep = ""
)
rosace_assayset_dir <- file.path(
  "./results",
  experiment_name,
  "rosace/assayset",
  sep = ""
)

sprintf("Creating directories: %s, %s", rosace_dir, rosace_assayset_dir)

dir.create(rosace_assayset_dir, recursive = TRUE)

# To create the rosace object, iterate over conditions and
# replicates to fill the assays from counts files. The set of "key" values
# is the unique set of experimental conditions, and the set of replicates
# is determined within each condition.

conditions <- unique(experiment_definition$condition)

for (expt_condition in conditions[conditions != baseline_condition]) {
  experimental_condition_df <- experiment_definition %>%
    filter(condition == expt_condition)
  for (expt_replicate in unique(experimental_condition_df$replicate)) {
    experimental_replicate_df <- experimental_condition_df %>%
      filter(replicate == expt_replicate)
    counts <- tibble(
      hgvs = character()
    )
    for (expt_time in unique(experimental_replicate_df$time)) {
      experimental_time_df <- experimental_replicate_df %>%
        filter(time == expt_time)
      for (tile in unique(experimental_time_df$tile)) {
        experimental_tile_df <- experimental_time_df %>% filter(tile == tile)
        if (dim(experimental_time_df)[1] > 1) {
          stop(sprintf(
            paste(
              "More than one row in the experiment",
              "definition for condition %s, replicate %s,",
              "tile %s, time %s."
            ), expt_condition,
            expt_replicate, tile, expt_time
          ))
        }
        if (dim(experimental_time_df)[1] == 0) {
          warning(sprintf(paste(
            "No rows in the experiment definition for",
            "condition %s, replicate %s, tile %s, time %s.",
            expt_condition, expt_replicate, tile, expt_time
          )))
        }

        file_name <- paste("results/",
          experiment_name,
          "/processed_counts/enrich_format/",
          experimental_tile_df$sample, ".tsv",
          sep = ""
        )
        file <- file.path(file_name)

        count <- read_tsv(file)

        col_name <- paste("c_", expt_time, sep = "")

        counts <- full_join(counts, count) %>% rename(!!col_name := count)
      }
    }
    assay <- CreateAssayObject(
      counts = as.matrix(counts[2:ncol(counts)]),
      var.names = count$hgvs,
      key = expt_condition, rep = expt_replicate, type = "growth"
    )
    if (exists("rosace")) {
      rosace <- AddAssayData(object = rosace, assay = assay)
    } else {
      rosace <- CreateRosaceObject(object = assay)
    }

    # Having generated an object, perform initial filtering,
    # imputation of missing values, and normalization
    rosace <- FilterData(rosace, key = expt_condition, na.rmax = 0.5)

    rosace <- ImputeData(
      rosace,
      key = expt_condition,
      impute.method = "knn",
      na.rmax = 0.5
    )

    rosace <- NormalizeData(rosace,
      key = expt_condition,
      normalization.method = "wt",
      wt.var.names = c("_wt"), wt.rm = TRUE
    )

    # These can now be combined into a single rosace object
    rosace <- IntegrateData(object = rosace, key = expt_condition)
  }
}

# The specific variant information now needs to be parsed

rosace@var.data <- rosace@var.data %>%
  rowwise() %>%
  mutate(
    tmp_n = substr(variants, 4, nchar(variants) - 1),
    tmp = list(parse_stripped_hgvs(tmp_n)),
    position = as.numeric(tmp[2]),
    wildtype = tmp[5],
    mutation = tmp[1],
    type = tmp[4]
  ) %>%
  dplyr::select(-tmp, -tmp_n)

# Finally, run rosace on each condition and write the output to a csv.

for (expt_condition in conditions[conditions != baseline_condition]) {
  sprintf("Running ROSACE for condition: %s", expt_condition)

  rosace <- RunRosace(
    object = rosace,
    name = expt_condition,
    type = "AssaySet",
    savedir = rosace_dir,
    pos.col = "position",
    ctrl.col = "type",
    ctrl.name = "synonymous",
    install = FALSE
  )

  scores.data <- OutputScore(
    rosace,
    name = paste(expt_condition, "_ROSACE", sep = "")
  )

  output_file_name <- paste(expt_condition, "_scores.csv", sep = "")

  output_file <- file.path(rosace_dir, output_file_name)

  write.csv(scores.data, file = output_file)

  sprintf("Results written to: %s", output_file)
}

# Close connection to log file
sink(type = "message")
sink(type = "output")
