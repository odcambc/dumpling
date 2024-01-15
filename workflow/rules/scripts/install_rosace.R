# Install renv next
install.packages("renv", repos = c("https://cloud.r-project.org"))

library("renv")

renv::restore()

library("cmdstanr")

# Install cmdstan
install_cmdstan(overwrite = FALSE)