
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
scenario = 2
source(here('code', 'sims_scenarios_mtp.R'))

# set seed
if (scenario==1) (set.seed(9751128))
if (scenario==2) (set.seed(3636127))

# generate the data
sim_data = lapply(1:num_sims, function(sim) 
  gen_staggered_data(
    n = n, 
    shape1_alpha_f = shape1_alpha_f,
    shape2_beta_f  = shape2_beta_f,
    T_periods = T_periods,
    c_fun = c_fun,
    theta_t = theta_t,
    cohort_predictor = cohort_predictor
  )
)

# save data
save_path <- here('data', 'sims', paste0('sim_data_mtp_sc_', scenario, '.rds'))
saveRDS(sim_data, save_path)

