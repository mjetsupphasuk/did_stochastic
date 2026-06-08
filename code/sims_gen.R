
# set directories
setwd("/work/users/m/j/mjets/dissertation/paper2")
.libPaths("/nas/longleaf/home/mjets/RLibs")

# load libraries
library(here)
library(tidyverse)

# set ggplot theme
theme_set(theme_bw())

# import functions
source(here('code', 'utils.R'))

# number of simulations
num_sims = 1000

# data generating scenario
scenario = 1
source(here('code', 'sims_scenarios.R'))

# set seed
if (scenario==1) (set.seed(7821111))
if (scenario==2) (set.seed(2283811))

# generate the data
sim_data = lapply(1:num_sims, function(sim) gen_data(n, pi_bin, mu.0_f, shape1_alpha_f, shape2_beta_f, mu.d_f))

# save data
save_path = here('data', 'sims', paste0('sim_data_sc_', scenario, '.rds'))
saveRDS(sim_data, save_path)

