# Set up logging
# Output will be redicted to the log file
log_file <- snakemake@log[[1]]

output_file <- snakemake@output[[1]]

# Open connection to log file
con <- file(log_file, open = "wt")

sink(con, append = FALSE, split = FALSE, type = "message")
sink(con, append = FALSE, split = FALSE, type = "output")



# Install renv next
install.packages("renv", repos = c("https://cloud.r-project.org"))

library("renv")

renv::restore()

library("cmdstanr")

# Install cmdstan
install_cmdstan(overwrite = FALSE)

# Check install
check_cmdstan_toolchain()

# Write to output file upon successful installation

con <- file(output_file, open = "wt")

writeLines("Rosace installed.", con)
close(con)