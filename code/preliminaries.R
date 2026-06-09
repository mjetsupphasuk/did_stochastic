
# load libraries
library(fs)
library(here)

# create directories (if don't already exist)
dir_create(here('data', 'fracking'), recurse = TRUE)
dir_create(here('data', 'sims'), recurse = TRUE)
dir_create(here('results', 'fracking', 'main'), recurse = TRUE)
dir_create(here('results', 'sims', 'interim'), recurse = TRUE)