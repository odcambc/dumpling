#!/usr/bin/env Rscript

renv::activate()

library("readr")
library("dplyr")
library("cmdstanr")
library("rosace")



# Install cmdstan
install_cmdstan(overwrite = FALSE)

# This metadata is provided by snakemake invocation.

experiment_name <- snakemake@config[["experiment"]]
experiment_file <- snakemake@config[["experiment_file"]]
baseline_condition <- snakemake@config[["baseline_condition"]]

# Set up logging
# Output will be redicted to the log file
log_file <- snakemake@log[[1]]

# Open connection to log file
con <- file(log_file, open = "wt")

sink(con, append = FALSE, split = FALSE, type = "message")
sink(con, append = FALSE, split = FALSE, type = "output")

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
          "/processed_counts/",
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
colnames(rosace@var.data)
rosace@var.data <- rosace@var.data %>%
  mutate(
    tmp = substr(variants, 4, nchar(variants) - 1),
    position = as.numeric(gsub("[[:alpha:]]", "", tmp)),
    wildtype = substr(tmp, 1, 1),
    tmp = substr(tmp, 2, nchar(tmp)),
    mutation = gsub("[[:digit:]]", "", tmp)
  ) %>%
  dplyr::select(-tmp)

func_map <- function(wt, mut) {
  if (nchar(wt) == 0) {
    return("NA")
  }

  if (wt == mut) {
    return("synonymous")
  } else if (mut == "del") {
    return("deletion")
  } else {
    return("missense")
  }
}

rosace@var.data <- rosace@var.data %>%
  rowwise() %>%
  mutate(type = func_map(wildtype, mutation)) %>%
  ungroup()

# Finally, run rosace on each condition and write the output to a csv.


for (expt_condition in conditions[conditions != baseline_condition]) {
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
}

# Close connection to log file
sink(type = "message")
sink(type = "output")
